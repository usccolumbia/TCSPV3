import csv
import hashlib
import itertools
import json
import math
from .mmap_cifs import MMapCIFStore
import pickle
import re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
from pymatgen.core import Composition, Element, Species, Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from .oxidation import legacy_evidence

MATCH_SETTINGS = dict(ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True,
                      scale=True, attempt_supercell=False)
MATCHER = StructureMatcher(**MATCH_SETTINGS)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False))
    tmp.replace(path)


@lru_cache(maxsize=200000)
def props(formula):
    comp = Composition(formula).reduced_composition
    d = comp.get_el_amt_dict()
    if any(n <= 0 or abs(n-round(n)) > 1e-7 for n in d.values()):
        raise ValueError('Requires positive integer stoichiometry')
    return d, tuple(sorted(d.values())), comp.reduced_formula


def symmetry(structure):
    try:
        return int(SpacegroupAnalyzer(structure, symprec=0.01, angle_tolerance=5).get_space_group_number())
    except Exception:
        return None


def simultaneous_substitute(source, mapping, query):
    result = source.copy()
    result.remove_oxidation_states()
    if not result.is_ordered:
        raise ValueError('Disordered source')
    if set(mapping) != set(result.composition.get_el_amt_dict()):
        raise ValueError('Mapping does not cover source species')
    result.replace_species(mapping)
    if props(result.composition.formula)[2] != props(query)[2]:
        raise ValueError('Composition changed incorrectly')
    return result


def merge_routes(a, b):
    seen = set()
    out = []
    for pair in itertools.zip_longest(a, b):
        for row in pair:
            if row is None:
                continue
            key = (row['formula'], tuple(sorted(row['mapping'].items())))
            if key not in seen:
                seen.add(key)
                out.append(row)
    return out


def unique(candidates, limit=20, same_sg_only=False):
    kept, structures = [], []
    duplicates = []
    for c in candidates:
        s = Structure.from_dict(c['structure'])
        duplicate = None
        for k, previous in zip(kept, structures):
            if same_sg_only and c['sg'] != k['sg']:
                continue
            if MATCHER.fit(s, previous):
                duplicate = k['candidate_id']
                break
        if duplicate:
            duplicates.append({'candidate_id': c['candidate_id'], 'representative': duplicate})
        else:
            kept.append(c)
            structures.append(s)
        if len(kept) == limit:
            break
    return kept, duplicates


def radius(element, oxidation=None, coordination=None):
    roman = {2: 'II', 3: 'III', 4: 'IV', 5: 'V', 6: 'VI', 7: 'VII', 8: 'VIII', 9: 'IX', 10: 'X', 12: 'XII'}
    if oxidation is not None and oxidation != 0 and coordination in roman:
        for spin in ('', 'High Spin', 'Low Spin'):
            try:
                val = float(Species(element, oxidation).get_shannon_radius(roman[coordination], spin=spin))
                if val > 0:
                    return val, 'shannon'
            except (KeyError, ValueError):
                pass
    e = Element(element)
    for attr in ('atomic_radius', 'atomic_radius_calculated'):
        val = getattr(e, attr, None)
        if val is not None and np.isfinite(float(val)) and val > 0:
            return float(val), attr
    return None, 'missing'


def geometry(source, mapping, query_evidence, template_evidence):
    """Nearest-shell CN is a geometric proxy, not a bond-valence assignment."""
    qa = query_evidence.get('assignments', [])
    ta = template_evidence.get('assignments', [])
    qo = qa[0]['states'] if qa else {}
    to = ta[0]['states'] if ta else {}
    ratios, coverage = [], defaultdict(int)
    nearest = source.get_all_neighbors(5.0)
    min_contact = float('inf')
    for i, site in enumerate(source):
        ds = [float(n.nn_distance) for n in nearest[i] if n.nn_distance > 1e-8]
        mind = min(ds) if ds else None
        if mind is not None:
            min_contact = min(min_contact, mind)
        cn = sum(d <= 1.25*mind for d in ds) if mind else None
        a, b = site.specie.symbol, mapping[site.specie.symbol]
        rt, st = radius(a, to.get(a), cn)
        rq, sq = radius(b, qo.get(b), cn)
        # Do not compare an ionic radius to an atomic radius.
        if (st == 'shannon') != (sq == 'shannon'):
            rt, st = radius(a)
            rq, sq = radius(b)
        coverage['ionic' if st == sq == 'shannon' else 'atomic' if rt and rq else 'missing'] += 1
        if rt and rq:
            ratios.append(rq/rt)
    # Normal separation between successive planes along each lattice direction.
    gaps = []
    for axis in range(3):
        coords = np.sort(source.frac_coords[:, axis] % 1.)
        gap = float(np.max(np.diff(np.r_[coords, coords[0]+1.])))
        height = 1. / np.linalg.norm(source.lattice.reciprocal_lattice_crystallographic.matrix[axis])
        gaps.append(gap*height)
    return {'radius_log_mismatch': float(np.mean(np.abs(np.log(ratios)))) if ratios else 0.,
            'radius_coverage': dict(coverage),
            'radius_scale': float(np.clip(np.median(ratios), .85, 1.18)) if ratios else 1.,
            'shortest_contact_A': min_contact if math.isfinite(min_contact) else None,
            'max_axis_gap_A': max(gaps), 'possible_vacuum': bool(max(gaps) > 8.),
            'coordination_method': 'neighbors within 1.25 times nearest distance, cutoff 5 A'}


class Predictor:
    def __init__(self, database, index, exclude=(), cif_store=None):
        self.database = Path(database)
        self.index = Path(index)
        self.cif_store = None
        default_data = self.database/'data'
        store_path = Path(cif_store) if cif_store else default_data/'source_cifs.bin'
        if not store_path.exists() and not cif_store:
            store_path = default_data/'source_cifs.pkl'
        if store_path.exists():
            if store_path.suffix == '.bin':
                self.cif_store = MMapCIFStore(store_path)
            elif store_path.suffix == '.pkl':
                with store_path.open('rb') as handle:
                    packed = pickle.load(handle)  # Trusted locally generated pickle only.
                if packed.get('format') != 'tcspv3-cif-store-v1':
                    raise ValueError('Unsupported CIF store format')
                self.cif_store = packed['cifs']
            else:
                raise ValueError('CIF store must be .bin or .pkl')
        elif cif_store:
            raise FileNotFoundError(store_path)
        self.blocked = {props(f)[2] for f in exclude}
        self.embedding = {k: np.asarray(v) for k, v in json.loads((self.database/'data/matscholar-embedding.json').read_text()).items()}
        with (Path(__file__).parent/'data/villars_pn.csv').open(encoding='utf-8-sig') as handle:
            self.pn = {r['Symbol']: int(r['PN']) for r in csv.DictReader(handle)}
        self.buckets = defaultdict(list)
        self.form_ids = defaultdict(list)
        idx = pickle.loads(self.index.read_bytes())  # Only a trusted local index.
        self.index_rows = sum(map(len, idx.values()))
        seen = set()
        for key, rows in idx.items():
            signature = tuple(sorted(float(x or 1) for x in re.findall(r'[A-Z](\d*)', key)))
            for mid, formula in rows:
                if (mid, formula) not in seen:
                    self.buckets[signature].append((mid, formula))
                    self.form_ids[formula].append(mid)
                    seen.add((mid, formula))

    @lru_cache(maxsize=None)
    def distance(self, a, b):
        if a == b:
            return 0.
        if a not in self.embedding or b not in self.embedding:
            # Explicit fallback rather than silently zeroing missing embeddings.
            return abs(self.pn[a]-self.pn[b])/6. + float(Element(a).group != Element(b).group)
        return float(np.linalg.norm(self.embedding[a]-self.embedding[b])) + float(Element(a).group != Element(b).group)

    def xscore(self, q, t):
        key = lambda a: (float(Element(a).X or 0), Element(a).Z)
        pairs = list(zip(sorted(props(q)[0], key=key), sorted(props(t)[0], key=key)))
        changed = [(a, b) for a, b in pairs if a != b]
        return sum(self.distance(a, b) for a, b in changed)/max(1, len(changed))

    def assignments(self, q, t, qe=None, te=None, pn_radius=None, alternatives=1, worst_pn=False):
        qd, td = props(q)[0], props(t)[0]
        qs, ts = sorted(qd), sorted(td)
        qe, te = qe or {}, te or {}
        use_os = bool(qe.get('assignments')) and bool(te.get('assignments'))
        cost = np.full((len(qs), len(ts)), np.inf)
        for i, a in enumerate(qs):
            for j, b in enumerate(ts):
                if qd[a] != td[b] or (pn_radius is not None and abs(self.pn[a]-self.pn[b]) > pn_radius):
                    continue
                value = self.distance(a, b)
                if use_os:
                    penalty = sum(x['weight']*y['weight']*min(abs(x['states'][a]-y['states'][b])/2., 2.)
                                  for x in qe['assignments'] for y in te['assignments'])
                    value += .5 * min(qe['neutral_mass'], te['neutral_mass']) * penalty
                cost[i, j] = qd[a]/sum(qd.values()) * value
        try:
            ii, jj = linear_sum_assignment(cost)
        except ValueError:
            return []
        if not np.isfinite(cost[ii, jj]).all():
            return []
        choices = [tuple(jj)]
        if alternatives > 1:
            # Equal-multiplicity groups reduce factorial enumeration drastically.
            groups = sorted(set(qd.values()))
            count = math.prod(math.factorial(sum(td[t] == mult for t in ts)) for mult in groups)
            if count > 40320:
                raise ValueError('Too many mappings; request one assignment')
            blocks = [list(itertools.permutations([j for j, t in enumerate(ts) if td[t] == mult])) for mult in groups]
            choices = []
            for selected in itertools.product(*blocks):
                js = [None]*len(qs)
                for mult, block in zip(groups, selected):
                    for i, j in zip([i for i, e in enumerate(qs) if qd[e] == mult], block):
                        js[i] = j
                if np.isfinite(cost[np.arange(len(qs)), js]).all():
                    choices.append(tuple(js))
        result = []
        for js in choices:
            mp = {ts[j]: qs[i] for i, j in enumerate(js)}
            maxpn = max(abs(self.pn[a]-self.pn[b]) for a, b in mp.items())
            pn_penalty = .1 * min(maxpn/6., 2.) if worst_pn else 0.
            result.append({'formula': t, 'mapping': mp, 'score': float(cost[np.arange(len(qs)), js].sum())+pn_penalty,
                           'pn_max': maxpn, 'pn_penalty': pn_penalty})
        result.sort(key=lambda r: (r['score'], tuple(sorted(r['mapping'].items()))))
        return [r for r in result if r['score'] <= result[0]['score']+.25][:alternatives]

    def retrieve(self, query, limit=200):
        good = sorted({f for _, f in self.buckets[props(query)[1]] if props(f)[2] not in self.blocked and props(f)[2] != props(query)[2]})
        xs = sorted(((self.xscore(query, f), f) for f in good))[:limit]
        aa = sorted(((self.assignments(query, f)[0]['score'], f) for f in good))[:limit]
        pn = []
        radius_used = None
        for rad in (4, 6):
            for f in good:
                mappings = self.assignments(query, f, pn_radius=rad)
                if mappings:
                    pn.append(mappings[0])
            if pn:
                radius_used = rad
                break
        pn.sort(key=lambda r: (r['score'], r['formula']))
        return {'query': query, 'formula_count': len(good), 'embedding': sorted({f for _, f in xs+aa}),
                'pn': [r['formula'] for r in pn[:limit]], 'pn_radius': radius_used,
                'pn_eligible_formulas': len(pn)}

    def rank(self, plan, mode, evidence, legacy):
        q = plan['query']
        v3 = mode == 'v3'
        def route(formulas, pn=False):
            result = []
            for f in formulas:
                qe = evidence.get(q, {}) if v3 else legacy_evidence(q, legacy.get(q, {}))
                te = evidence.get(f, {}) if v3 else legacy_evidence(f, legacy.get(f, {}))
                result.extend(self.assignments(q, f, qe, te, plan['pn_radius'] if pn else None,
                                                alternatives=3 if v3 else 1, worst_pn=v3))
            return sorted(result, key=lambda r: (r['score'], r['formula'], tuple(sorted(r['mapping'].items()))))
        if mode == 'baseline':
            return route(plan['embedding'])
        if mode == 'pn':
            return route(plan['pn'], True)
        return merge_routes(route(plan['embedding']), route(plan['pn'], True))

    @lru_cache(maxsize=4096)
    def source(self, mid):
        if Path(mid).name != mid:
            raise ValueError('Invalid template ID')
        path = self.database/'data/mixed_data4'/(mid+'.cif')
        if self.cif_store is not None:
            raw = self.cif_store[mid]
            s = Structure.from_str(raw.decode('utf-8'), fmt='cif')
            digest = hashlib.sha256(raw).hexdigest()
        else:
            s = Structure.from_file(path)
            digest = sha256(path)
        s.remove_oxidation_states()
        return s, digest

    def generate(self, plan, mode, evidence, legacy, pool_size=100, topk=20, target_atoms=None):
        ranked = self.rank(plan, mode, evidence, legacy)
        q = plan['query']
        pool, failures, seen = [], [], set()
        for rank in ranked:
            for mid in sorted(set(self.form_ids[rank['formula']])):
                cid = hashlib.sha256(json.dumps([mid, sorted(rank['mapping'].items())]).encode()).hexdigest()[:20]
                if cid in seen:
                    continue
                seen.add(cid)
                try:
                    source, digest = self.source(mid)
                    if target_atoms is not None and len(source) != target_atoms:
                        continue
                    actual = props(source.composition.formula)[2]
                    if actual in self.blocked or actual == props(q)[2]:
                        raise ValueError('Actual source composition excluded')
                    if actual != props(rank['formula'])[2]:
                        raise ValueError('Index/source composition mismatch')
                    s = simultaneous_substitute(source, rank['mapping'], q)
                    features = geometry(source, rank['mapping'], evidence.get(q, {}), evidence.get(rank['formula'], {}))
                    if features['shortest_contact_A'] is not None and features['shortest_contact_A'] < .5:
                        raise ValueError('Contact below 0.5 A')
                    penalty = .25*features['radius_log_mismatch'] if mode == 'v3' else 0.
                    atoms_per_formula = int(round(len(s) / sum(props(q)[0].values())))
                    pool.append(dict(rank, candidate_id=cid, template_id=mid, source_sha256=digest,
                                     nsites=len(s), formula_units=atoms_per_formula,
                                     source_formula=actual, chemical_score=rank['score'], score=rank['score']+penalty,
                                     geometry_penalty=penalty, geometry=features, sg=symmetry(s), structure=s.as_dict()))
                except Exception as exc:
                    failures.append({'template_id': mid, 'formula': rank['formula'], 'error': str(exc)})
                if len(pool) >= pool_size:
                    break
            if len(pool) >= pool_size:
                break
        if mode == 'v3':
            pool.sort(key=lambda r: (r['score'], r['formula'], r['template_id'], r['candidate_id']))
        selected, duplicates = unique(pool, topk, same_sg_only=mode == 'baseline')
        return {'query': q, 'mode': mode, 'plan': plan, 'target_atoms': target_atoms, 'pool': pool, 'selected': selected,
                'duplicates': duplicates, 'failures': failures}


def symmetry_vote(candidates, strength=.1, minority_slots=2):
    """Optional family-weighted vote; candidates must already be deduplicated."""
    weights = defaultdict(float)
    for c in candidates:
        if c['sg'] is not None:
            weights[c['sg']] += math.exp(-min(c['score'], 100))
    total = sum(weights.values())
    if not total:
        return candidates
    ranked = sorted(candidates, key=lambda c: c['score']-strength*weights[c['sg']]/total)
    majority = max(weights, key=weights.get)
    minority = [c for c in ranked if c['sg'] != majority][:minority_slots]
    ids = {c['candidate_id'] for c in minority}
    rest = [c for c in ranked if c['candidate_id'] not in ids]
    split = max(0, min(5-len(minority), len(rest)))
    return rest[:split]+minority+rest[split:]

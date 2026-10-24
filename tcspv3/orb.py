"""ORB-v3 relaxation and energy ranking: the default final stage of TCSPV3.

The template stage returns up to 20 candidates in its native chemical-score order. This stage relaxes every candidate with ORB-v3
(``orb-v3-conservative-inf-omat``; FIRE with ASE ``ExpCellFilter``, at most 500 steps, force threshold 0.01 eV/A) and ranks the relaxed
structures by ORB-v3 energy per atom, the same final stage used for the deep generative models in the benchmark. Nonconverged structures
are kept. Candidates whose relaxation fails or changes the composition are dropped from the ranking and recorded.
Requires the optional dependency ``orb-models`` (listed in requirements.txt); use ``--ranker native`` to keep the chemical-score order.
Reference structures are never read.
"""
import importlib.util
import subprocess
import time
import warnings

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from .core import props

MODEL = 'orb_v3_conservative_inf_omat'
STEPS = 500
FMAX = 0.01


def available():
    """True when the optional ORB-v3 dependencies can be imported."""
    return all(importlib.util.find_spec(m) is not None for m in ('orb_models', 'ase', 'torch'))


def require():
    if not available():
        raise RuntimeError('The default ORB-v3 final stage needs the optional dependency orb-models '
                           '(pip install "orb-models>=0.5.5" ase). Use --ranker native to keep the chemical-score order.')


def smi_load(torch):
    """GPU utilization from nvidia-smi, matched to torch device order by UUID."""
    out = subprocess.run(['nvidia-smi', '--query-gpu=uuid,utilization.gpu', '--format=csv,noheader,nounits'],
                         capture_output=True, text=True, timeout=10, check=True).stdout
    by_uuid = {u.strip().removeprefix('GPU-'): float(v) for u, v in (line.split(',') for line in out.strip().splitlines())}
    return [by_uuid[str(torch.cuda.get_device_properties(i).uuid)] for i in range(torch.cuda.device_count())]


def pick_device(torch):
    """Least-busy visible GPU as 'cuda:N', or 'cpu' without CUDA."""
    if not torch.cuda.is_available():
        return 'cpu'
    for read in (lambda: [torch.cuda.utilization(i) for i in range(torch.cuda.device_count())], lambda: smi_load(torch)):
        try:
            load = read()
        except Exception:  # NVML or nvidia-smi unavailable
            continue
        return f'cuda:{load.index(min(load))}'
    return 'cuda'


class OrbStage:
    def __init__(self, device='auto', model=MODEL, weights_path=None, precision='float32-high', steps=STEPS, fmax=FMAX):
        require()
        import torch
        from orb_models.forcefield import pretrained
        if device == 'auto':
            device = pick_device(torch)
        if device.startswith('cuda') and not torch.cuda.is_available():
            raise RuntimeError('CUDA was requested but is unavailable')
        factory = getattr(pretrained, model, None)
        if factory is None:
            raise ValueError('Unknown ORB model: ' + model)
        kwargs = {'device': device, 'precision': precision}
        if weights_path:
            kwargs['weights_path'] = str(weights_path)
        self.model_name, self.device, self.precision, self.steps, self.fmax = model, device, precision, steps, fmax
        self.weights_path = str(weights_path) if weights_path else None
        loaded = factory(**kwargs)
        # orb-models < 0.5 returns (model, atoms_adapter); 0.5 and later return the model alone
        self.orbff, self.adapter = loaded if isinstance(loaded, tuple) else (loaded, None)
        self.orbff.eval()

    def provenance(self):
        return {'stage': 'orb-v3 relaxation and energy ranking', 'model': self.model_name, 'weights_path': self.weights_path,
                'device': self.device, 'precision': self.precision, 'optimizer': 'FIRE', 'filter': 'ExpCellFilter',
                'steps': self.steps, 'fmax_eV_A': self.fmax, 'nonconverged_structures_kept': True}

    def relax(self, structure):
        """Relax a pymatgen Structure; returns (relaxed Structure, converged, energy per atom in eV)."""
        from ase.optimize import FIRE
        try:
            from ase.filters import ExpCellFilter
        except ImportError:                                   # older ASE
            from ase.constraints import ExpCellFilter
        atoms = AseAtomsAdaptor.get_atoms(structure)
        try:                                                  # orb-models >= 0.5
            from orb_models.forcefield.calculator import ORBCalculator
            atoms.calc = ORBCalculator(self.orbff, device=self.device)
        except ImportError:                                   # older releases
            from orb_models.forcefield.inference.calculator import ORBCalculator
            atoms.calc = ORBCalculator(model=self.orbff, atoms_adapter=self.adapter, device=self.device)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            converged = bool(FIRE(ExpCellFilter(atoms), logfile=None).run(fmax=self.fmax, steps=self.steps))
        energy = float(atoms.get_potential_energy())
        return AseAtomsAdaptor.get_structure(atoms), converged, energy / len(atoms)


def rank(candidates, relax, topk=None, progress=None):
    """Relax and rank candidate dicts (core.Predictor output). ``relax(structure) -> (structure, converged, energy_per_atom)``.

    Returns (ranked candidates, failures). Each ranked candidate keeps its keys; ``structure`` becomes the relaxed structure and
    ``orb_v3`` records the rank, energy per atom, convergence flag and the candidate's native rank.
    ``progress()`` is called after each candidate.
    """
    scored, failures = [], []
    for native_rank, c in enumerate(candidates, 1):
        started = time.time()
        try:
            source = Structure.from_dict(c['structure'])
            relaxed, converged, energy = relax(source)
            if props(relaxed.composition.formula)[2] != props(source.composition.formula)[2]:
                raise ValueError('ORB-v3 relaxation changed the composition')
            if not energy == energy or energy in (float('inf'), float('-inf')):
                raise ValueError('Nonfinite ORB-v3 energy')
            out = dict(c)
            out['structure'] = relaxed.as_dict()
            out['orb_v3'] = {'energy_per_atom_eV': energy, 'converged': converged, 'native_rank': native_rank,
                             'elapsed_s': round(time.time() - started, 2)}
            scored.append(out)
        except Exception as exc:
            failures.append({'candidate_id': c.get('candidate_id'), 'native_rank': native_rank, 'error': type(exc).__name__ + ': ' + str(exc)})
        if progress:
            progress()
    scored.sort(key=lambda c: (c['orb_v3']['energy_per_atom_eV'], c['orb_v3']['native_rank']))
    for i, c in enumerate(scored, 1):
        c['orb_v3']['rank'] = i
    return (scored[:topk] if topk else scored), failures


def finalize_directory(destination, stage, topk):
    """Apply the ORB-v3 stage to an existing native prediction directory (prediction.json + tcspv3_###.cif).

    Used by the parallel CSV tool, where the native stage runs in CPU workers and the ORB-v3 stage runs afterwards in one process.
    The unrelaxed native candidates are kept in ``native/``. Returns (candidates kept, ORB-v3 failures).
    """
    import json
    from pathlib import Path
    from .core import dump
    folder = Path(destination)
    path = folder / 'prediction.json'
    result = json.loads(path.read_text())
    native = result['selected']
    ranked, failures = rank(native, stage.relax, topk)
    if not ranked:
        raise RuntimeError('ORB-v3 stage produced no valid candidate: ' + json.dumps(failures)[:500])
    (folder / 'native').mkdir(exist_ok=True)
    for i, c in enumerate(native, 1):
        Structure.from_dict(c['structure']).to(filename=str(folder / 'native' / ('tcspv3_%03d.cif' % i)))
    for old in folder.glob('tcspv3_*.cif'):
        old.unlink()
    for i, c in enumerate(ranked, 1):
        Structure.from_dict(c['structure']).to(filename=str(folder / ('tcspv3_%03d.cif' % i)))
    result['native_selected_ids'] = [c['candidate_id'] for c in native]
    result['selected'] = ranked
    result['orb_v3_failures'] = failures
    result.setdefault('provenance', {}).update({'ranker': 'orb-v3', 'final_stage': stage.provenance()})
    dump(path, result)
    return len(ranked), len(failures)

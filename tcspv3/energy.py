"""ASE energy ranking/limited relaxation; consumes predictions, never references."""
import argparse
import hashlib
import importlib.metadata
import json
import signal
import time
import warnings
from pathlib import Path
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from .core import dump, sha256, props


def timeout_handler(signum, frame):
    raise TimeoutError('Per-candidate 180 second limit')


def rank_energy(candidates):
    return sorted(candidates, key=lambda c: (c.get('energy_per_atom') is None,
                                             c.get('energy_per_atom') if c.get('energy_per_atom') is not None else float('inf'),
                                             c.get('score', 0.), c['candidate_id']))


def consensus_rank(arms):
    """Optional rank fusion for complete model outputs; no energy-offset mixing."""
    if not arms:
        return []
    ids = [{c['candidate_id'] for c in candidates} for candidates in arms]
    if any(s != ids[0] for s in ids):
        raise ValueError('Consensus requires identical candidate IDs')
    ranks = [{c['candidate_id']: i for i, c in enumerate(rank_energy(cs))} for cs in arms]
    return sorted(arms[0], key=lambda c: (sum(r[c['candidate_id']] for r in ranks), c['candidate_id']))


class EnergyStage:
    def __init__(self, device='cuda', checkpoint=None, threads=4):
        import torch
        from mattersim.forcefield import MatterSimCalculator
        from mattersim.applications.relax import Relaxer
        if threads < 1:
            raise ValueError('threads must be positive')
        torch.set_num_threads(threads)
        self.calculator = MatterSimCalculator(device=device, **({'load_path': checkpoint} if checkpoint else {}))
        self.relaxer = Relaxer(optimizer='FIRE', filter='ExpCellFilter', constrain_symmetry=False)
        signal.signal(signal.SIGALRM, timeout_handler)
        # The loaded parameter bytes identify the model even if no checkpoint path is given.
        digest = hashlib.sha256()
        for name, tensor in sorted(self.calculator.potential.model.state_dict().items()):
            digest.update(name.encode())
            digest.update(tensor.detach().cpu().contiguous().numpy().tobytes())
        self.provenance = {'calculator': 'MatterSimCalculator', 'device': device,
                           'torch_threads': torch.get_num_threads(),
                           'stage_version': 2,
                           'state_dict_sha256': digest.hexdigest(),
                           'checkpoint': checkpoint, 'checkpoint_sha256': sha256(checkpoint) if checkpoint else None,
                           'packages': {p: importlib.metadata.version(p) for p in ['mattersim','torch','ase','pymatgen']},
                           'relaxation': {'optimizer':'FIRE','filter':'ExpCellFilter','steps':100,'fmax':.05,'timeout_s':180,'symmetry_constrained':False}}

    def compute(self, candidate, scale=False, relax=False):
        started = time.time()
        result = {k: v for k, v in candidate.items() if k != 'structure'}
        result.update(energy_per_atom=None, max_force_eV_A=None, converged=False, optimizer_converged=False, error=None)
        try:
            signal.alarm(180)
            s = Structure.from_dict(candidate['structure'])
            if scale:
                factor = candidate['geometry']['radius_scale']
                s.scale_lattice(s.volume*factor**3)
                result['applied_linear_scale'] = factor
            atoms = AseAtomsAdaptor.get_atoms(s)
            atoms.calc = self.calculator
            if relax:
                converged, atoms = self.relaxer.relax(atoms, steps=100, fmax=.05, verbose=False)
                result['optimizer_converged'] = bool(converged)
            energy = float(atoms.get_potential_energy()) / len(atoms)
            forces = atoms.get_forces()
            if not np.isfinite(energy) or not np.isfinite(forces).all():
                raise ValueError('Nonfinite model output')
            output = AseAtomsAdaptor.get_structure(atoms)
            if props(output.composition.formula)[2] != props(s.composition.formula)[2]:
                raise ValueError('Physical stage composition mismatch')
            maximum_force = float(np.linalg.norm(forces,axis=1).max())
            if not np.isfinite(output.volume) or output.volume <= 0 or not np.isfinite(output.cart_coords).all():
                raise ValueError('Invalid output geometry')
            # Cell filters can make their force criterion differ from atomic forces.
            result['converged'] = result['optimizer_converged'] and maximum_force <= .05
            result.update(energy_per_atom=energy, max_force_eV_A=maximum_force,
                          structure=output.as_dict())
        except Exception as exc:
            result.update(error=type(exc).__name__+': '+str(exc), structure=None if relax else candidate['structure'])
        finally:
            signal.alarm(0)
        result['elapsed_s'] = time.time()-started
        return result


def main():
    warnings.filterwarnings('ignore')
    p = argparse.ArgumentParser()
    p.add_argument('--run', required=True, help='Benchmark run folder, or standalone prediction directory')
    p.add_argument('--shard', type=int, default=0)
    p.add_argument('--shards', type=int, default=1)
    p.add_argument('--device', default='cpu')
    p.add_argument('--checkpoint', default=str(Path(__file__).resolve().parent.parent/'checkpoints/mattersim-v1.0.0-1M.pth'))
    args = p.parse_args()
    if not 0 <= args.shard < args.shards:
        p.error('Invalid shard')
    root = Path(args.run)
    if (root/'prediction.json').exists():
        tasks = [('prediction', root/'prediction.json')]
    else:
        queries = json.loads((root/'queries.json').read_text())
        tasks = [(q['id'], root/'predictions/v3'/(q['id']+'.json')) for i,q in enumerate(queries) if i % args.shards == args.shard]
    stage = EnergyStage(args.device, args.checkpoint)
    dump(root/'physical'/('model_%s.json' % args.shard), stage.provenance)
    for case, path in tasks:
        out = root/'physical'/(case+'.json')
        identity = {'prediction_sha256':sha256(path), 'model_state_sha256':stage.provenance['state_dict_sha256'],
                    'stage_version':stage.provenance['stage_version']}
        if out.exists():
            previous = json.loads(out.read_text())
            if previous.get('identity') != identity:
                raise ValueError('Existing physical output has different inputs/model/code; use a fresh run directory: '+str(out))
            continue
        source = json.loads(path.read_text())
        candidates = source['selected']
        checkpoint = root/'physical'/(case+'_progress.json')
        progress = json.loads(checkpoint.read_text()) if checkpoint.exists() else {'raw':{},'scaled':{},'relaxed':{},'identity':identity}
        if progress.get('identity') != identity:
            raise ValueError('Progress checkpoint identity mismatch: '+str(checkpoint))
        for c in candidates:
            cid = c['candidate_id']
            for key, scale in [('raw',False),('scaled',True)]:
                if cid not in progress[key]:
                    progress[key][cid] = stage.compute(c, scale=scale)
                    dump(checkpoint, progress)
        raw = rank_energy([progress['raw'][c['candidate_id']] for c in candidates])
        scaled = rank_energy([progress['scaled'][c['candidate_id']] for c in candidates])
        heads = {'v3_relaxed':candidates[:5], 'singlepoint_relaxed':raw[:5]}
        for cs in heads.values():
            for c in cs:
                cid = c['candidate_id']
                if cid not in progress['relaxed']:
                    # Always relax the same original geometry, including shared IDs.
                    original = next(x for x in candidates if x['candidate_id'] == cid)
                    progress['relaxed'][cid] = stage.compute(original, relax=True)
                    dump(checkpoint, progress)
        arms = {'singlepoint':raw,'scaled_singlepoint':scaled}
        for key, cs in heads.items():
            arms[key] = rank_energy([progress['relaxed'][c['candidate_id']] for c in cs])
        dump(out, {'id':case,'identity':identity,'model':stage.provenance,'arms':arms,'unique_relaxed':len(progress['relaxed']),
                   'source_prediction_sha256':sha256(path)})
        for name, cs in arms.items():
            folder = root/'cifs'/name/case
            folder.mkdir(parents=True, exist_ok=True)
            for i, c in enumerate(cs, 1):
                if c.get('structure') is not None:
                    Structure.from_dict(c['structure']).to(filename=str(folder/('%03d.cif' % i)))
        print('PHYSICAL_DONE',case,'candidates',len(candidates),'relaxed',len(progress['relaxed']),
              'converged',sum(c['converged'] for c in progress['relaxed'].values()),flush=True)
    (root/('PHYSICAL_DONE_%s' % args.shard)).write_text('All assigned cases processed; inspect failures/convergence\n')


if __name__ == '__main__':
    main()

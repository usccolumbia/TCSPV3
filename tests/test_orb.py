import sys
import unittest
from unittest import mock

from pymatgen.core import Lattice, Structure

from tcspv3 import orb
from tcspv3.__main__ import build_parser


def candidate(cid, a, species=('Na', 'Cl')):
    s = Structure(Lattice.cubic(a), list(species), [[0, 0, 0], [.5, .5, .5]])
    return {'candidate_id': cid, 'score': 0.0, 'structure': s.as_dict()}


class RankTest(unittest.TestCase):
    def test_ranks_by_relaxed_energy_and_keeps_native_rank(self):
        cands = [candidate('a', 4.0), candidate('b', 4.5), candidate('c', 5.0)]
        energies = {4.0: -3.0, 4.5: -5.0, 5.0: -4.0}

        def fake_relax(structure):
            return structure, True, energies[round(structure.lattice.a, 1)]
        ranked, failures = orb.rank(cands, fake_relax)
        self.assertEqual([c['candidate_id'] for c in ranked], ['b', 'c', 'a'])
        self.assertEqual([c['orb_v3']['rank'] for c in ranked], [1, 2, 3])
        self.assertEqual([c['orb_v3']['native_rank'] for c in ranked], [2, 3, 1])
        self.assertEqual(failures, [])

    def test_topk_truncates_and_failures_are_recorded(self):
        cands = [candidate('a', 4.0), candidate('b', 4.5), candidate('c', 5.0)]

        def fake_relax(structure):
            if round(structure.lattice.a, 1) == 4.5:
                raise RuntimeError('boom')
            return structure, False, -1.0 * structure.lattice.a
        ranked, failures = orb.rank(cands, fake_relax, topk=1)
        self.assertEqual([c['candidate_id'] for c in ranked], ['c'])
        self.assertFalse(ranked[0]['orb_v3']['converged'])      # nonconverged structures are kept
        self.assertEqual([f['candidate_id'] for f in failures], ['b'])

    def test_composition_change_is_rejected(self):
        cands = [candidate('a', 4.0)]
        bad = Structure(Lattice.cubic(4.0), ['Na', 'Br'], [[0, 0, 0], [.5, .5, .5]])
        ranked, failures = orb.rank(cands, lambda s: (bad, True, -1.0))
        self.assertEqual(ranked, [])
        self.assertIn('composition', failures[0]['error'])


class CliTest(unittest.TestCase):
    def test_orb_is_the_default_ranker(self):
        args = build_parser().parse_args(['--formula', 'NaCl'])
        self.assertEqual(args.ranker, 'orb-v3')
        self.assertEqual((args.orb_steps, args.orb_fmax), (500, 0.01))

    def test_native_ranker_is_selectable(self):
        args = build_parser().parse_args(['--formula', 'NaCl', '--ranker', 'native'])
        self.assertEqual(args.ranker, 'native')

    def test_missing_dependency_gives_a_clear_message(self):
        with mock.patch.object(orb, 'available', return_value=False):
            with self.assertRaises(RuntimeError) as ctx:
                orb.OrbStage()
        self.assertIn('--ranker native', str(ctx.exception))
        self.assertIn('orb-models', str(ctx.exception))


if __name__ == '__main__':
    unittest.main()


class DirectoryTest(unittest.TestCase):
    def test_finalize_directory_rewrites_outputs(self):
        import json
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            cands = [candidate('a', 4.0), candidate('b', 4.5)]
            (folder / 'prediction.json').write_text(json.dumps({'selected': cands, 'provenance': {}}))
            for i, c in enumerate(cands, 1):
                Structure.from_dict(c['structure']).to(filename=str(folder / ('tcspv3_%03d.cif' % i)))

            class FakeStage:
                def relax(self, structure):
                    return structure, True, -structure.lattice.a   # larger cell = lower energy
                def provenance(self):
                    return {'stage': 'fake'}
            kept, failed = orb.finalize_directory(folder, FakeStage(), 1)
            self.assertEqual((kept, failed), (1, 0))
            result = json.loads((folder / 'prediction.json').read_text())
            self.assertEqual([c['candidate_id'] for c in result['selected']], ['b'])
            self.assertEqual(result['native_selected_ids'], ['a', 'b'])
            self.assertEqual(result['provenance']['ranker'], 'orb-v3')
            self.assertEqual(len(list(folder.glob('tcspv3_*.cif'))), 1)
            self.assertEqual(len(list((folder / 'native').glob('*.cif'))), 2)


class PickDeviceTest(unittest.TestCase):
    def fake_torch(self, available=True, load=(99, 0)):
        cuda = mock.Mock(is_available=lambda: available, device_count=lambda: len(load),
                         utilization=lambda i: load[i])
        return mock.Mock(cuda=cuda)

    def test_picks_least_busy_gpu(self):
        self.assertEqual(orb.pick_device(self.fake_torch(load=(99, 0))), 'cuda:1')
        self.assertEqual(orb.pick_device(self.fake_torch(load=(5, 40, 5))), 'cuda:0')

    def test_cpu_without_cuda_and_plain_cuda_without_nvml(self):
        self.assertEqual(orb.pick_device(self.fake_torch(available=False)), 'cpu')
        t = self.fake_torch()
        t.cuda.utilization = mock.Mock(side_effect=ModuleNotFoundError('pynvml'))
        with mock.patch.object(orb, 'smi_load', side_effect=FileNotFoundError('nvidia-smi')):
            self.assertEqual(orb.pick_device(t), 'cuda')

    def test_nvidia_smi_fallback_matches_uuids(self):
        t = self.fake_torch(load=(0, 0))
        t.cuda.utilization = mock.Mock(side_effect=ModuleNotFoundError('pynvml'))
        t.cuda.get_device_properties = lambda i: mock.Mock(uuid=['aaa', 'bbb'][i])
        smi = mock.Mock(stdout='GPU-bbb, 0\nGPU-aaa, 100\n')   # nvidia-smi order differs from torch order
        with mock.patch.object(orb.subprocess, 'run', return_value=smi):
            self.assertEqual(orb.pick_device(t), 'cuda:1')

    def test_orb_device_values(self):
        from tcspv3.__main__ import build_parser
        parse = lambda *a: build_parser().parse_args(['--formula', 'NaCl', *a]).device
        self.assertEqual(parse(), 'auto')
        self.assertEqual(parse('--device', 'cuda:1'), 'cuda:1')
        self.assertEqual(parse('--orb-device', 'cpu'), 'cpu')   # old name still works

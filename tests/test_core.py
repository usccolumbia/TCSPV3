import unittest
from unittest.mock import patch
from pymatgen.core import Lattice, Structure
from tcspv3.core import Predictor, props, simultaneous_substitute, unique
from tcspv3.oxidation import evidence_from_probabilities


class Invariants(unittest.TestCase):
    def engine(self):
        e = Predictor.__new__(Predictor)
        e.pn = {'Li': 1, 'Na': 2, 'K': 3, 'Rb': 4, 'O': 8, 'S': 9}
        e.distance = lambda a, b: 0. if a == b else abs(e.pn[a]-e.pn[b])
        return e

    def test_formula_order_and_scale(self):
        self.assertEqual(props('O6Fe4')[1:], props('Fe2O3')[1:])

    def test_simultaneous_cyclic_substitution(self):
        s = Structure(Lattice.cubic(5), ['Na','K','O'], [[0,0,0],[.5,.5,0],[.5,0,.5]])
        x = simultaneous_substitute(s, {'Na':'K','K':'Na','O':'O'}, 'NaKO')
        self.assertEqual([site.specie.symbol for site in x], ['K','Na','O'])
        self.assertEqual(s[0].specie.symbol, 'Na')

    def test_stoichiometric_assignment(self):
        e = self.engine()
        mappings = e.assignments('Li2O', 'NaS2', alternatives=3)
        self.assertEqual(mappings[0]['mapping'], {'Na':'O','S':'Li'})

    def test_radius_is_constraint_before_optimization(self):
        e = self.engine()
        e.pn.update({'Li':0, 'Na':4, 'K':3, 'Rb':7})
        # Crossed mapping has cheaper mean cost but violates max distance 4.
        e.distance = lambda a,b: 0. if (a,b) in [('Li','Rb'),('Na','K')] else 2.
        result = e.assignments('LiNa', 'KRb', pn_radius=4)
        self.assertEqual(result[0]['mapping'], {'K':'Li','Rb':'Na'})
        self.assertEqual(e.assignments('LiNa','KRb',pn_radius=2), [])

    def test_neutral_evidence_and_missing(self):
        ev = evidence_from_probabilities('NaCl', {'Na':{1:.9,0:.1}, 'Cl':{-1:.8,0:.2}})
        self.assertAlmostEqual(ev['neutral_mass'], .74)
        self.assertEqual(ev['assignments'][0]['states'], {'Cl':-1,'Na':1})
        self.assertEqual(evidence_from_probabilities('NaCl', {'Na':{1:1.}})['neutral_mass'],0.)
        self.assertEqual(evidence_from_probabilities('LiNa', {'Li':{1:1.},'Na':{-1:1.}})['reason'],'all_metal')

    def test_large_mapping_budget_checked_before_allocation(self):
        e = self.engine()
        e.distance = lambda a, b: 0.
        with patch('tcspv3.core.itertools.permutations', side_effect=AssertionError('Allocated before budget check')):
            with self.assertRaisesRegex(ValueError, 'Too many mappings'):
                e.assignments('LiNaKBeMgCaScTiVO9', 'LiNaKBeMgCaScTiVO9', alternatives=3)

    def test_geometry_dedup_not_spacegroup_only(self):
        a = Structure(Lattice.cubic(5), ['Na','Cl'], [[0,0,0],[.5,.5,.5]])
        b = a.copy(); b.translate_sites([0,1], [.2,.2,.2])
        cs = [{'candidate_id':'a','sg':1,'structure':a.as_dict()}, {'candidate_id':'b','sg':2,'structure':b.as_dict()}]
        selected, dup = unique(cs)
        self.assertEqual(len(selected),1)
        self.assertEqual(dup[0]['representative'],'a')


if __name__ == '__main__':
    unittest.main()

"""Physics checks for measured visible tau observables."""
import math
import unittest
import numpy as np
from analyze_visible_taus import observables, choose_matches, reco_four


class VisibleTests(unittest.TestCase):
    def test_back_to_back_photons(self):
        r=observables([[1,1,0,0],[1,-1,0,0]],[0,0,1])
        self.assertAlmostEqual(r['mass_GeV'],2)
        self.assertEqual(r['scalar_pt_GeV'],2)
        self.assertEqual(r['vector_pt_GeV'],0)
        self.assertTrue(math.isnan(r['angle_deg']))

    def test_collinear_and_perpendicular(self):
        r=observables([[5,3,4,0]],[3,4,0])
        self.assertAlmostEqual(r['mass_GeV'],0)
        self.assertAlmostEqual(r['angle_deg'],0)
        self.assertAlmostEqual(observables([[5,3,4,0]],[0,0,1])['angle_deg'],90)

    def test_neutral_energy_direction(self):
        np.testing.assert_allclose(reco_four([0,0,2],10,0,0),[10,0,0,10])
        self.assertIsNone(reco_four([0,0,0],10,0,0))

    def test_charged_track_mass(self):
        np.testing.assert_allclose(reco_four([3,4,0],99,2,1),[math.sqrt(29),3,4,0])

    def test_duplicate_track_preferred_over_neutral(self):
        selected,rejected=choose_matches([0,1,2],[4,4,5],[1,.8,1],[0,1,0])
        self.assertEqual(selected[4],(1,.8))
        self.assertEqual(selected[5],(2,1))
        self.assertEqual(rejected,[(4,0,1)])

    def test_one_truth_per_reco_and_nonpositive_weights(self):
        selected,_=choose_matches([0,0,1],[3,4,5],[.1,.9,0],[1,0])
        self.assertEqual(selected,{4:(0,.9)})

    def test_empty_and_spacelike(self):
        self.assertTrue(math.isnan(observables([],[1,0,0])['mass_GeV']))
        self.assertTrue(math.isnan(observables([[1,2,0,0]],[1,0,0])['mass_GeV']))


if __name__=='__main__': unittest.main()

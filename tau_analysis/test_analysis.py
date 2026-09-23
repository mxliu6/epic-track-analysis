"""Focused checks for ancestry and beam-axis kinematics."""
import math
import unittest
from analyze_taus import kinematics, relations, stable_descendants

class PhysicsLogicTests(unittest.TestCase):
    def test_stops_at_generator_final_state(self):
        # A stable generator pion may have detector-generated children.
        stable, unresolved = stable_descendants(0, [[1,2],[3],[],[]], [2,1,1,0])
        self.assertEqual(stable,[1,2]); self.assertEqual(unresolved,[])

    def test_cascade_and_shared_descendant(self):
        stable, unresolved = stable_descendants(0,[[1,2],[3],[3,4],[],[]],[2,2,2,1,1])
        self.assertEqual(stable,[3,4]); self.assertEqual(unresolved,[])

    def test_cycle_rejected(self):
        with self.assertRaises(ValueError): stable_descendants(0,[[1],[0]],[2,2])

    def test_wrong_collection_rejected(self):
        event={'PDG':[15,16],'daughters.index':[1],'daughters.collectionID':[8],
               'daughters_begin':[0,1],'daughters_end':[1,1]}
        with self.assertRaises(ValueError): relations(event,'daughters',7)
        event['daughters.collectionID']=[7]
        self.assertEqual(relations(event,'daughters',7),[[1],[]])

    def test_kinematics(self):
        pt,eta=kinematics(3,4,5)
        self.assertEqual(pt,5); self.assertAlmostEqual(eta,math.asinh(1))
        self.assertEqual(kinematics(0,0,-1)[1],-math.inf)
        self.assertTrue(math.isnan(kinematics(0,0,0)[1]))

if __name__=='__main__': unittest.main()

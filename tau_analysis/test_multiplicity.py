import math
import unittest
from daughter_multiplicity import accepted_eta

class AcceptanceTests(unittest.TestCase):
    def test_strict_symmetric_eta_window(self):
        for eta in [-4,-3.5,-3.49,0,3.49,3.5,4]:
            value,keep=accepted_eta(1,0,math.sinh(eta))
            self.assertAlmostEqual(value,eta)
            self.assertEqual(keep,-3.5<eta<3.5)

    def test_undefined_direction_rejected(self):
        self.assertFalse(accepted_eta(0,0,1)[1])
        self.assertFalse(accepted_eta(0,0,0)[1])

if __name__=='__main__':unittest.main()

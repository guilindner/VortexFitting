import unittest
import numpy as np
import sys

sys.path.insert(1,'../vortexfitting')

import schemes

"""
test using "nosetests /tests/test_schemes.py"
"""
class SchemeTest(unittest.TestCase):
    

    def test_my(self):
        x = np.arange(0.0,1.0,0.01)
        func = np.sin(np.pi*x) 
        res = np.pi*(np.cos(np.pi*x))
        res_my = schemes.second_order_diff(func)
        self.assertAlmostEqual(res,res_my,6)

#    def test_numpy(self):
#        res_np = np.gradient(func,0.01)
        
if __name__ == '__main__':
    unittest.main()

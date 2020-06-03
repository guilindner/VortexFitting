import unittest
from nose import with_setup # optional

import numpy as np
import sys
sys.path.insert(1,'../src')

import fitting

"""
test using "nosetests -v tests/test_tools.py"
"""


# Dummy field to use
sample_field = np.array([[1.1, 0.9, 1.3, 0.7],
                         [2.1, 1.9, 2.3, 1.7],
                         [3.1, 2.9, 3.3, 2.7],
                         [4.1, 3.9, 4.3, 3.7]])

class FittingTest(unittest.TestCase):
    

#    def my_teardown_function():
#	print ("my_teardown_function")
    
    def test_correlation_coef_1_1_1_1(self):
        # Exact same u and v compared with exact same u and v
        result_calc = fitting.correlation_coef(sample_field,
                                               sample_field,
                                               sample_field,
                                               sample_field)

        np.testing.assert_almost_equal(result_calc, 1.0)
        
    def test_correlation_coef_2_05_2_05(self):
        # Different u and v compared with the same u and v
        result_calc = fitting.correlation_coef(sample_field*2.0,
                                               sample_field*0.5,
                                               sample_field*2.0,
                                               sample_field*0.5)

        np.testing.assert_almost_equal(result_calc, 1.0)
        
    def test_correlation_coef_2_2_0_0(self):
        # Different u and v compared with u=0 and v=v
        result_calc = fitting.correlation_coef(sample_field*2.0,
                                               sample_field*2.0,
                                               sample_field,
                                               sample_field)

        np.testing.assert_almost_equal(result_calc, 0.5)
        
    def test_correlation_coef_05_05_1_1(self):
        # Different u and v compared with u=0 and v=v
        result_calc = fitting.correlation_coef(sample_field*0.5,
                                               sample_field*0.5,
                                               sample_field,
                                               sample_field)

        np.testing.assert_almost_equal(result_calc, 0.5)

    def test_correlation_coef_025_025_1_1(self):
        # Different u and v compared with u=0 and v=v
        result_calc = fitting.correlation_coef(sample_field*0.25,
                                               sample_field*0.25,
                                               sample_field,
                                               sample_field)

        np.testing.assert_almost_equal(result_calc, 0.25)

    def test_correlation_coef_1_0_1_1(self):
        # Different u and v compared with u=0 and v=v
        result_calc = fitting.correlation_coef(sample_field,
                                               sample_field*0.,
                                               sample_field,
                                               sample_field)

        np.testing.assert_almost_equal(result_calc, 0.5)

    def test_correlation_coef_1_1_0_1(self):
        # Different u and v compared with u=0 and v=v
        result_calc = fitting.correlation_coef(sample_field,
                                               sample_field,
                                               sample_field*0.0,
                                               sample_field)

        np.testing.assert_almost_equal(result_calc, 0.5)

if __name__ == '__main__':
    unittest.main()

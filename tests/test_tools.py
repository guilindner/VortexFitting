import unittest
import numpy as np

import tools

"""
test using "nosetests tests/test_tools.py"
"""


class ToolsTest(unittest.TestCase):

    def test_get_fluc(self):
        sample_field = np.array([[1.1, 0.9, 1.3, 0.7],
                                 [2.1, 1.9, 2.3, 1.7],
                                 [3.1, 2.9, 3.3, 2.7],
                                 [4.1, 3.9, 4.3, 3.7]])

        # Test in hom. in x axis
        mean = np.array([1.0, 2.0, 3.0, 4.0])

        result_calc = tools.get_fluc(sample_field, mean, 'x')

        result = np.array([[0.1, -0.1, 0.3, -0.3],
                           [0.1, -0.1, 0.3, -0.3],
                           [0.1, -0.1, 0.3, -0.3],
                           [0.1, -0.1, 0.3, -0.3]])

        np.testing.assert_array_almost_equal(result_calc, result)

        # Test in hom. in y axis
        mean = np.array([2.6, 2.4, 2.8, 2.2])

        result_calc = tools.get_fluc(sample_field, mean, 'y')

        result = np.array([[-1.5, -1.5, -1.5, -1.5],
                           [-0.5, -0.5, -0.5, -0.5],
                           [0.5, 0.5, 0.5, 0.5],
                           [1.5, 1.5, 1.5, 1.5]])

        np.testing.assert_array_almost_equal(result_calc, result)


if __name__ == '__main__':
    unittest.main()

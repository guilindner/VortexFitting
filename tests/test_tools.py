import unittest
import numpy as np

import tools
"""
test using "nosetests /tests/test_tools.py"
"""
class ToolsTest(unittest.TestCase):
    

    def test_sub_mean_x(self):
        field2D = np.array([[1.1, 0.9, 1.3, 0.7],
                            [2.1, 1.9, 2.3, 1.7],
                            [3.1, 2.9, 3.3, 2.7],
                            [4.1, 3.9, 4.3, 3.7]])
        sub_field_x = np.array([[0.1, -0.1, 0.3, -0.3],
                                [0.1, -0.1, 0.3, -0.3],
                                [0.1, -0.1, 0.3, -0.3],
                                [0.1, -0.1, 0.3, -0.3]])
        
        my_field_x = tools.sub_mean(field2D,1)
   
        np.testing.assert_array_almost_equal(sub_field_x, my_field_x)

    def test_sub_mean_y(self):
        field2D = np.array([[1.1, 0.9, 1.3, 0.7],
                            [2.1, 1.9, 2.3, 1.7],
                            [3.1, 2.9, 3.3, 2.7],
                            [4.1, 3.9, 4.3, 3.7]])
        sub_field_y = np.array([[-1.5, -1.5, -1.5, -1.5],
                                [-0.5, -0.5, -0.5, -0.5],
                                [0.5, 0.5, 0.5, 0.5],
                                [1.5, 1.5, 1.5, 1.5]])
        
        my_field_y = tools.sub_mean(field2D,0)

        np.testing.assert_array_almost_equal(sub_field_y, my_field_y)
        
if __name__ == '__main__':
    unittest.main()

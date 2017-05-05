import numpy as np

def sub_mean(x, hom_axis):
    mean = np.mean(x, axis=hom_axis)
    x = x - mean[:,None]
    return x

def normalize(x, hom_axis):
    mean = np.mean(x, axis=hom_axis)
    x = x/mean[:,None]
    where_are_NaNs = np.isnan(x)
    x[where_are_NaNs] = 0
    return x

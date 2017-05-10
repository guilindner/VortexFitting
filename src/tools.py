import numpy as np
from classes import VelocityField

np.seterr(divide='ignore', invalid='ignore')

def sub_mean(x, hom_axis):
    mean = np.mean(x, axis=hom_axis)
    if (hom_axis == 0):
        x = x - mean
    else:
        x = x - mean[:,None]
    return x
      
def normalize(x, hom_axis):
    mean = np.mean(x, axis=hom_axis)
    x = x/mean
    where_are_NaNs = np.isnan(x)
    x[where_are_NaNs] = 0
    return x

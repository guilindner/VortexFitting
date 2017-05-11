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
    mean = np.mean(x**2, axis=hom_axis)
    x = x/np.sqrt(mean)
    return x

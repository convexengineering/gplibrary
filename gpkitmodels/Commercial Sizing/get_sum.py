"Function to simply integrate the dh values to get height values"
import numpy as np

def get_hft(dhft):
    """
    Returns an altitude vector with alt in feet,
    corresponds to altitudes at each segemnt of discretized climb
    """
    hft = []
    for i in range(0,len(dhft)):
        if i == 0:
            hft[i] = dhft[i]
        else:
            hft[i] = dhft[i]+hft[i-1]
    return hft

def get_cruiserng(ReqRng, RngClimb):
    """
    returns the required cruise range
    """
    return ReqRng - RngClimb
    

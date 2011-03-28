"""
This module wraps the dynamic library clustersearch.dylib.

It defines three functions for external use:
  calculate_measures_as_tuple
  calculate_measures
  reseed
"""

from ctypes.util import find_library
from ctypes import cdll
from ctypes import c_uint
from ctypes import c_double

libclustersearch = cdll.LoadLibrary(find_library('clustersearch'))

# load function reseed
_reseed = libclustersearch.reseed
_reseed.argtypes = [ c_uint ]

# load function calculate_measures
_calculate_measures = libclustersearch.calculate_measures
_calculate_measures.argtypes = [ c_uint, c_uint, c_uint, c_double ]

# define result struct for calculate_measure
from ctypes import Structure
class cluster_measure(Structure):
    _fields_ = [("cluster_size", c_uint),
                ("perimeter_size", c_uint),
                ("colors", c_uint),
                ("exits_size", c_uint),
                ("robustness", c_double)]
    def __repr__(self):
        s= """cluster_size   = %s
perimeter_size = %s
colors         = %s
exits_size     = %s
robustness     = %s""" % (self.cluster_size, self.perimeter_size, self.colors, self.exits_size, self.robustness)
        return s
    pass

_calculate_measures.restype = cluster_measure

def reseed(seed):
    """Seeds the random number generator.

    Reseeding to the same start value (e.g., 0), should cause
    subsequent searches to produce the same result. But it does not. I
    don't know why.
    """
    _reseed(seed)

def calculate_measures_as_tuple(alphabetsize,length,colors,gray=0.0):
    "Returns tuple of basic measures from a single cluster search"
    x = _calculate_measures(length,alphabetsize,colors,gray)
    return (x.cluster_size ,x.perimeter_size ,x.colors ,x.exits_size, x.robustness)

def calculate_measures(alphabetsize,length,colors,gray=0.0):
    "Returns dict of basic measures from a single cluster search"
    measurenames = ["cluster_size", "perimeter_size", "colors", "exits_size", "robustness"]
    t = calculate_measures_as_tuple(alphabetsize,length,colors,gray)
    d = dict(zip(measurenames,t))
    return d

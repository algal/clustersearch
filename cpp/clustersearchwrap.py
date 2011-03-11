"""
This module wraps the dynamic library clustersearch.dylib.

It defines three functions for external use:
  calculate_measures_as_tuple
  calculate_measures
  srand
"""

from ctypes.util import find_library
from ctypes import cdll
from ctypes import c_uint

libclustersearch = cdll.LoadLibrary(find_library('clustersearch'))

# load function srand
_srand = libclustersearch.srand
_srand.argtypes = [ c_uint ]

# load function calculate_measures
_calculate_measures = libclustersearch.calculate_measures
_calculate_measures.argtypes = [ c_uint,c_uint,c_uint]

#from ctypes import *
from ctypes import Structure
class cluster_measure(Structure):
    _fields_ = [("cluster_size", c_uint),
                ("perimeter_size", c_uint),
                ("colors", c_uint),
                ("exits_size", c_uint)]
    def __repr__(self):
        s= """cluster_size   = %s
perimeter_size = %s
colors         = %s
exits_size     = %s""" % (self.cluster_size ,self.perimeter_size ,self.colors ,self.exits_size)
        return s
    pass

_calculate_measures.restype = cluster_measure

def calculate_measures_as_tuple(length,alphabetsize,colors):
    "Returns tuple of basic measures from a single cluster search"
    x = _calculate_measures(length,alphabetsize,colors)
    return (x.cluster_size ,x.perimeter_size ,x.colors ,x.exits_size)

def srand(seed):
    """Seeds the random number generator.

    Reseeding to the same start value (e.g., 0), should cause
    subsequent searches to produce the same result. But it does not. I
    don't know why.
    """
    _srand(seed)

def calculate_measures(length,alphabetsize,colors):
    "Returns dict of basic measures from a single cluster search"
    t = calculate_measures_as_tuple(length,alphabetsize,colors)
    d = dict()
    d["cluster_size"] = t[0]
    d["perimeter_size"] = t[1]
    d["colors"] = t[2]
    d["exits_size"] = t[3]
    return d


    

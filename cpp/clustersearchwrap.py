from ctypes.util import find_library
from ctypes import cdll
from ctypes import c_uint

libclustersearch = cdll.LoadLibrary(find_library('clustersearch'))

# define wrapper for C func cluster_size
foo = libclustersearch.cluster_size
foo.argtypes = [ c_uint,c_uint,c_uint]

calculate_measures = libclustersearch.calculate_measures
calculate_measures.argtypes = [ c_uint,c_uint,c_uint]

srand = libclustersearch.srand
srand.argtypes = [ c_uint ]

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
exits_size     = %s""" %        (self.cluster_size   ,self.perimeter_size ,self.colors         ,self.exits_size     )
        return s
    pass



calculate_measures.restype = cluster_measure

def display(x):
    print x.cluster_size
    print x.perimeter_size
    print x.colors
    print x.exits_size
    





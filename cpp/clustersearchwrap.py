from ctypes.util import find_library
from ctypes import cdll
from ctypes import c_uint

libclustersearch = cdll.LoadLibrary(find_library('clustersearch'))

# define wrapper for C func cluster_size
foo = libclustersearch.cluster_size
foo.argtypes = [ c_uint,c_uint,c_uint]

calculate_measures = libclustersearch.calculate_measures
calculate_measures.argtypes = [ c_uint,c_uint,c_uint]







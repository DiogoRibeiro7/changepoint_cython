import ctypes

# Load the shared library
stat_methods = ctypes.CDLL('./cost_functions.so')

# Define the argument and return types for the C functions
stat_methods.mll_mean.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_longlong]
stat_methods.mll_mean.restype = ctypes.c_double
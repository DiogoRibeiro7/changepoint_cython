import ctypes

# Load the shared library
stat_methods = ctypes.CDLL('./cost_functions.so')

# Define the argument and return types for the C functions
stat_methods.mll_mean.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_longlong]
stat_methods.mll_mean.restype = ctypes.c_double

# Now you can call it like a regular Python function
result = stat_methods.mll_mean(1.0, 2.0, 3.0, 4)
print("Result:", result)

import os, ctypes
import numpy as np
from scipy import integrate, LowLevelCallable, special

class ARGS(ctypes.Structure):
    _fields_ = [("q", ctypes.c_double),
                ("k", ctypes.c_double),
                ("v", ctypes.c_double)]

comppath = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + 'studentizedrange/studentized_range.so'
lib = ctypes.CDLL(comppath)

lib._Z19studentized_range_piPdPv.restype = ctypes.c_double
lib._Z19studentized_range_piPdPv.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_void_p)

lib._Z30studentized_range_p_asymptoticiPdPv.restype = ctypes.c_double
lib._Z30studentized_range_p_asymptoticiPdPv.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_void_p)

def studentized_cdf(q, k, v):
    if v > 40:
        user_data = ctypes.cast(ctypes.pointer(ARGS(q, k, v)), ctypes.c_void_p)
        res = integrate.quad(LowLevelCallable(lib._Z30studentized_range_p_asymptoticiPdPv, user_data), -np.inf, np.inf)[0]
        return k * res
    else:

        user_data = ctypes.cast(ctypes.pointer(ARGS(q, k, v)), ctypes.c_void_p)
        res = integrate.dblquad(LowLevelCallable(lib._Z19studentized_range_piPdPv, user_data), 0, np.inf,
                          gfun=-np.inf, hfun=np.inf)[0]
        return np.sqrt(2 * np.pi) * k * v ** (v / 2) / (
                special.gamma(v / 2) * 2 ** (v / 2 - 1)) * res





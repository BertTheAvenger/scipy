# distutils: language = c++

from .cy_studentized_range cimport studentized_range_p

def cy_studentized_range_p(double q, double k, double v):
    return studentized_range_p(q, k, v)

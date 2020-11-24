//Run "gcc -shared -fPIC -o studentized_range.so studentized_range.cpp" to compile
#include <cmath>

struct RangeArg {
    double q;
    double k;
    double v;
};

//PDF approximation
double Phi(double value) {
   return 0.5 * erfc(-value * M_SQRT1_2);
}

//PDF approximation
double phi(double z)
{
    static const double inv_sqrt_2pi = 0.3989422804014327;
    return inv_sqrt_2pi * std::exp(-0.5f * z * z);
}

double studentized_range_p(int n, double *x, void *user_data) {
    RangeArg data = *(RangeArg *)user_data;

    double s = x[1];
    double z = x[0];

    double q = data.q;
    double k = data.k;
    double v = data.v;

    return pow(s, (v - 1)) * phi(sqrt(v) * s) *
        phi(z) * pow(Phi(z + q * s) - Phi(z), k - 1);
}

double studentized_range_p_asymptotic(int n, double *x, void *user_data) {
    RangeArg data = *(RangeArg *)user_data;

    double z = x[0];

    double q = data.q;
    double k = data.k;
    double v = data.v;

    return phi(z) * pow(Phi(z + q) - Phi(z), k - 1);
}
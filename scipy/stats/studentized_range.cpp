#include <math.h>
#include <random>
#include <tgmath.h>
#include <iostream>
#include <iomanip>
#include "studentized_range.h"

using namespace std;

//Approximate using trapazoidal rule.
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


//Inner integral. Confirmed
double inner(double z, double s, double q, double k) {
    return phi(z) * pow( Phi(z + q * s) - Phi(z), k - 1 );
}

//Confirmed using q = 1, k = 3, s = 1
//Uses Simpson's rule - integral ~ dx / 3 * [ f(x0) + 4f(x1) + 2f(x2) ... + 4f(xn-1) + f(xn) ]
double integrate_inner(double s, double q, double k, double lower, double upper, int num) {

    double dX = (upper - lower) / (double) num;
    double sFactor = dX / 3.0; //Simpson's rule factor

    double inner_int = inner(lower, s, q, k) + inner(upper, s, q, k);

    for(int i = 1; i < num; i++) {
        double z = i * dX + lower;

        inner_int += inner(z, s, q, k) * (i % 2 == 0 ? 2 : 4); //Follow simpson's rule when multiplying. Alternate 2 and 4
    }
    //cout << inner_int << endl;
    return sFactor * inner_int;
}


// Confirmed
double outer(double s, double q, double k, double v, double lower, double upper, int num) {
    return pow(s, v - 1) * phi(sqrt(v) * s) * integrate_inner(s, q, k, lower, upper, num);
}


//TODO: Possibly break out integration into own function? integrate_outer and integrate_inner do same thing.
//Uses Simpson's rule - integral ~ dx / 3 * [ f(x1) + 4f(x2) + 2f(x3) ... + 4f(xn-1) + f(xn) ]
double integrate_outer(double q, double k, double v, double lower, double upper, int num) {
    double dX = upper / (double) num;
    double sFactor = dX / 3.0; //Simpson's rule factor

    double outer_int = outer(0, q, k, v, lower, upper, num) + outer(upper, q, k, v, lower, upper, num);
    for(int i = 1; i < num; i++) {
        double s = (double) i * dX;
        outer_int += outer(s, q, k, v, lower, upper, num) * (i % 2 == 0 ? 2.0 : 4.0);
    }

    return sFactor * outer_int;
}


double studentized_range_p(double q, double k, double v) {
    //TODO: Dynamically determine lower/upper bounds + num of integrations based on inputs. Should probably be able to be calculated using simple-ish formula.
    double lower = -30;
    double upper = 30;
    double num = 200;

    //numerator / denom is confirmed.
    double numerator = sqrt(2 * M_PI) * k * pow(v, v / 2);
    double denom = tgamma(v / 2) * pow(2, v / 2 - 1);

    return numerator / denom * integrate_outer(q, k, v, lower, upper, num);

}
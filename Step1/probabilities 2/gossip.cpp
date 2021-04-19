#include <iostream>
#include <limits>
#include <iomanip>
#include <fstream>

#include "math/factorial.h"

long double disconnect_probability(const int & N, const long double & f, const int & G)
{
    int H = (1.l - f) * N;
    long double p = 1.l - powl(1.l - ((long double) G) / ((long double) N), 2);

    long double dp = 0.l;

    for(int k = 1; k <= H / 2; k++)
    {
        long double termlog = factorial_log(H) - factorial_log(k) - factorial_log(H - k);
        termlog += logl(1.l - p) * ((long double) k) * ((long double) (H - k));

        dp += expl(termlog);
    }

    return dp;
}

int main()
{
    int N = 524288;
    long double f = 0.20l;
    for(int G = 0; G < 128; G++)
        std :: cout << "(" << G << ", " << disconnect_probability(N, f, G) << ")" << std :: endl;
}

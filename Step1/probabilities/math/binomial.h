#ifndef __probabilities__math__binomial__h
#define __probabilities__math__binomial__h

// Libraries

#include <math.h>

// Includes

#include "factorial.h"

long double binomial(const long double & p, const int & n, const int & k)
{
    if(k >= 0 && k <= n)
    {
        if(p <= 0)
            return (k == 0 ? 1.l : 0.l);

        if(p >= 1)
            return (k == n ? 1.l : 0.l);

        return expl(factorial_log(n) - factorial_log(k) - factorial_log(n - k) + logl(p) * k + logl(1.l - p) * (n - k));
    }
    else
        return 0;
}

#endif

#ifndef __probabilities__math__any__h
#define __probabilities__math__any__h

// Libraries

#include <math.h>

long double any(const long double & p, const int & N)
{
    long double log_one_minus_p = 0;

    for(int i = 1;; i++)
    {
        long double term = powl(p, i) / i;
        long double next = log_one_minus_p - term;

        if(log_one_minus_p == next)
            break;

        log_one_minus_p = next;

        if(log_one_minus_p * N < -100.l)
            return 1;
    }

    return 1.l - expl(log_one_minus_p * N);
}

#endif

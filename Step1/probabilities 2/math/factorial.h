#ifndef __probabilities__math__factorial__h
#define __probabilities__math__factorial__h

// Libraries

#include <math.h>

class factorial_log_lookup
{
    // Settings

    static constexpr size_t size = 1048576;

    // Members

    long double _values[size];

public:

    // Constructors

    factorial_log_lookup()
    {
        this->_values[0] = 0;
        for(size_t i = 1; i < size; i++)
            this->_values[i] = this->_values[i - 1] + logl((long double) i);
    }

    // Operators

    const long double & operator () (const int & value)
    {
        if(value < size)
            return this->_values[value];
        else
            throw "Out of range";
    }
};

factorial_log_lookup factorial_log;

#endif

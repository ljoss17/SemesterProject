#ifndef __probabilities__utils__distribution__h
#define __probabilities__utils__distribution__h

#include <vector>

class distribution
{
    struct pool
    {
        size_t size = 0;
        std :: vector <long double *> blocks;
    };

    // Static members

    static pool * pools;

    // Members

    pool * _pool;
    long double * _block;

public:

    // Constructors

    distribution(const size_t & size) : _pool(nullptr)
    {
        for(this->_pool = pools; this->_pool->size; this->_pool++)
            if(this->_pool->size == size)
                break;

        if(!(this->_pool->size))
            this->_pool->size = size;

        if(this->_pool->blocks.size())
        {
            this->_block = this->_pool->blocks.back();
            this->_pool->blocks.pop_back();
        }
        else
            this->_block = new long double[size + 1];
    }

    distribution(const distribution & that) : distribution(that._pool->size)
    {
        for(size_t i = 0; i < this->_pool->size; i++)
            this->_block[i] = that._block[i];
    }

    distribution(distribution && that) : _pool(that._pool), _block(that._block)
    {
        that._block = nullptr;
    }

    // Destructor

    ~distribution()
    {
        if(this->_block)
            this->_pool->blocks.push_back(this->_block);
    }

    // Operators

    long double & operator [] (const size_t & index)
    {
        return this->_block[index];
    }

    const long double & operator [] (const size_t & index) const
    {
        return this->_block[index];
    }
};

// Static member declarations

distribution :: pool * distribution :: pools = new distribution :: pool[16]{};

#endif

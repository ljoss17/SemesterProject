#ifndef __probabilities__pebbling__pebbling__h
#define __probabilities__pebbling__pebbling__h

#include "../math/binomial.h"
#include "../utils/distribution.h"

class pebbling
{
    // Members

    int N;
    long double f;
    int D;
    int R;
    int B;
    int P;

public:

    // Constructors

    pebbling(const int & N, const long double & f, const int & D, const int & R, const int & B, const int & P) : N(N), f(f), D(D), R(R), B(B), P(P)
    {
    }

    // Methods

    // epsilon 0 Contagion validity with f
    distribution F()
    {
        distribution response(D);

        for(int F = 0; F <= D; F++)
            response[F] = binomial(f, D, F);

        return response;
    }

    // epsilon 0 Contagion validity with p = (1-f)*Nim/N
    distribution Vim1$Ni_Ui(const int & Ni, const int & Ui)
    {
        int Nim1 = Ni - Ui;
        long double p = (1.l - f) * ((long double) Nim1) / ((long double) N);

        distribution response(D);

        for(int V = 0; V <= D; V++)
            response[V] = binomial(p, D, V);

        return response;
    }

    distribution Vim1$nWi_Ni_Ui(const int & Ni, const int & Ui)
    {
        int Nim1 = Ni - Ui;
        long double p = (1.l - f) * ((long double) Nim1) / ((long double) N);

        auto term = [&](const int & V)
        {
            return binomial(p, D, V);
        };

        distribution response(D);
        long double total = 0;

        for(int V = 0; V < P; V++)
        {
            response[V] = term(V);
            total += response[V];
        }

        for(int V = 0; V < P; V++)
            response[V] /= total;

        for(int V = P; V <= D; V++)
            response[V] = 0;

        return response;
    }

    distribution Vim1$F_Ni_Ui(const int & F, const int & Ni, const int & Ui)
    {
        int Nim1 = Ni - Ui;
        long double p = ((long double) Nim1) / ((long double) N);

        distribution response(D);
        for(int V = 0; V <= D; V++)
            response[V] = binomial(p, D - F, V);

        return response;
    }

    distribution Vim1$F_nWi_Ni_Ui(const int & F, const int & Ni, const int & Ui)
    {
        int Nim1 = Ni - Ui;
        long double p = ((long double) Nim1) / ((long double) N);

        auto term = [&](const int & V)
        {
            return binomial(p, D - F, V);
        };

        distribution response(D);
        long double total = 0;

        for(int V = 0; V < P; V++)
        {
            response[V] = term(V);
            total += response[V];
        }

        for(int V = 0; V < P; V++)
            response[V] /= total;

        for(int V = P; V <= D; V++)
            response[V] = 0;

        return response;
    }

    /*distribution F$nWi_Ni_Ui(const int & Ni, const int & Ui)
    {
        distribution Fdist = this->F();
        distribution Vim1$Ni_Ui = this->Vim1$Ni_Ui(Ni, Ui);

        long double unconditioned = 0.l;
        for(int V = 0; V < P; V++)
            unconditioned += Vim1$Ni_Ui[V];

        distribution response(D);

        for(int F = 0; F <= D; F++)
        {
            long double conditioned = 0;
            distribution Vim1$F_Ni_Ui = this->Vim1$F_Ni_Ui(F, Ni, Ui);

            for(int V = 0; V < P; V++)
                conditioned += Vim1$F_Ni_Ui[V];

            response[F] = conditioned * Fdist[F] / unconditioned;
        }

        return response;
    }*/

    distribution Vi$Vim1_F_nWi_Ni_Ui(const int & Vim1, const int & F, const int & Ni, const int & Ui)
    {
        long double p = ((long double)(Ui) / (long double)(N - (Ni - Ui)));

        distribution response(D);
        for(int V = 0; V <= D; V++) {
            response[V] = binomial(p, D - F - Vim1, V - Vim1);
        }

        return response;
    }

    long double Wip1_nWi_Ni_Ui(const int & Ni, const int & Ui)
    {
        long double total = 0.l;
        distribution Fdist = this->F();
        for(int F = 0; F <= D; F++)
        {
            long double sumF = 0.l;
            distribution Vim1$F_nWi_Ni_Ui = this->Vim1$F_nWi_Ni_Ui(F, Ni, Ui);

            for(int Vim1 = 0; Vim1 < P; Vim1++)
            {
                long double sumV = 0.l;
                distribution Vi$Vim1_F_nWi_Ni_Ui = this->Vi$Vim1_F_nWi_Ni_Ui(Vim1, F, Ni, Ui);

                for(int Vi = P; Vi <= D; Vi++) {
                    sumV += Vi$Vim1_F_nWi_Ni_Ui[Vi];
                }
                sumV *= Vim1$F_nWi_Ni_Ui[Vim1];
                sumF += sumV;
            }

            sumF *= Fdist[F];
            total += sumF;
        }

        return total;
    }
};

#endif

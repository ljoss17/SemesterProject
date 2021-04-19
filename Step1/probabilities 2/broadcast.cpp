#include <iostream>
#include <iomanip>

#include "pebbling/pebbling.h"

int main()
{
    int N = 8192;

    int D = 128;
    int P = 48;

    long double f = 0.15l;

    int R = 1;
    int B = 2048;

    pebbling my_pebbling(N, f, D, R, B, P);

    for(int Ni = 2000; Ni <= 2000; Ni++)
        for(int Ui = 0; Ui <= Ni; Ui++)
            std :: cout << "(" << Ni << "/" << Ui << "): " << my_pebbling.Wip1_nWi_Ni_Ui(Ni, Ui) << std :: endl;
}

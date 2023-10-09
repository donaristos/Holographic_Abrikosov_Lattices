//
//  main.cpp
//  TestData
//
//  Created by Aristomenis Donos on 22/11/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#define _USE_MATH_DEFINES
//#define EIGEN_NO_DEBUG

#include <iostream>
#include "function.h"
#include "omp.h"
#include "MyTypes.h"

#include <paralution/paralution.hpp>

#ifdef F128
typedef float128 MyType;
#elif MP
typedef mpreal MyType;
#elif LONGD
typedef long double MyType;
#else
typedef double MyType;
#endif


typedef std::complex<MyType> MyComplex;

MyType **cheby_d1;
MyType **cheby_d2;
MyType **fourierX_d1;
MyType **fourierX_d2;
MyType **fourierY_d1;
MyType **fourierY_d2;
MyType **idr=0;
MyType **idx=0;
MyType **idy=0;

clock_t start1,start2;
double omp_start, omp_end;
MyType diff;


int main(int argc, const char * argv[])
{
    
//    std::cout << std::numeric_limits<int>::max() << std::endl;
//    std::cout << std::numeric_limits<long>::max() << std::endl;
//    std::cout << std::numeric_limits<long long>::max() << std::endl;
//    
//    std::cout << sizeof(int) << std::endl;
//    std::cout << sizeof(long) << std::endl;
//    std::cout << sizeof(long long) << std::endl;
    
    for(int ii=0; (ii<10); ii++){
        if (ii==3) {
            continue;
        }
        
        std::cout << ii << std::endl;
    }
    
    return 0;
}
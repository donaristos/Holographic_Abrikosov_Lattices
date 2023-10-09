//
//  MyTypes.h
//  OpticalConductivity
//
//  Created by Aristomenis Donos on 06/12/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#ifndef OpticalConductivity_MyTypes_h
#define OpticalConductivity_MyTypes_h

#include <boost/multiprecision/float128.hpp>
#include <unsupported/Eigen/MPRealSupport>
//#include <mpfr.h>
#define MKL_Complex16 std::complex<double>

using namespace boost::multiprecision;


typedef int64_t My_Int;
// typedef int My_Int;

#define UF_long My_Int

typedef std::complex<float128> ComplexQ;
typedef std::complex<double> Complex;
typedef std::complex<mpfr::mpreal> ComplexMP;
typedef mpfr::mpreal mpreal;
typedef std::complex<long double> ComplexLD;

//#include <mkl_types.h>

//typedef std::complex<long double> ComplexMP;
//typedef long double mpreal;

#endif

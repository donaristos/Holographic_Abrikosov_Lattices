//
//  CustomFunctions.h
//  2DEM
//
//  Created by Aristomenis Donos on 11/12/2014.
//  Copyright (c) 2014 Aristomenis Donos. All rights reserved.
//

#ifndef ___DEM__CustomFunctions__
#define ___DEM__CustomFunctions__

#include <iostream>
#include <complex>
#include "math.h"
#include "MyTypes.h"

//Oveload complex

std::complex<double> operator*(const int i, const std::complex<double>& z);

std::complex<double> operator+(const int i, const std::complex<double>& z);

std::complex<double> operator/(const int i, const std::complex<double>& z);

std::complex<double> operator/(const std::complex<double>& z,const int i);

std::complex<double> operator-(const int i, const std::complex<double>& z);

std::complex<float128> operator*(const double& i,const std::complex<float128>& z);

std::complex<float128> operator*(const std::complex<float128>& z,const double& i);

std::complex<float128> operator+(const double& i, const std::complex<float128>& z);

std::complex<float128> operator+(const std::complex<float128>& z,const double& i);

std::complex<float128> operator/(const double& i, const std::complex<float128>& z);

std::complex<float128> operator/(const std::complex<float128>& z,const double& i);

std::complex<float128> operator-(const double& i, const std::complex<float128>& z);

std::complex<float128> operator-(const std::complex<float128>& z,const double& i);

std::complex<mpreal> operator*(const double& i, const std::complex<mpreal>& z);

std::complex<mpreal> operator*(const std::complex<mpreal>& z,const double& i);

std::complex<mpreal> operator+(const double& i, const std::complex<mpreal>& z);

std::complex<mpreal> operator+(const std::complex<mpreal>& z,const double& i);

std::complex<mpreal> operator/(const double& i, const std::complex<mpreal>& z);

std::complex<mpreal> operator/(const std::complex<mpreal> &z,const double& i);

std::complex<mpreal> operator-(const double& i, const std::complex<mpreal>& z);

std::complex<mpreal> operator-(const std::complex<mpreal>& z,const double& i);

std::complex<long double> operator*(const double& i, const std::complex<long double>& z);

std::complex<long double> operator*(const std::complex<long double>& z,const double& i);

std::complex<long double> operator+(const double& i, const std::complex<long double>& z);

std::complex<long double> operator+(const std::complex<long double>& z,const double& i);

std::complex<long double> operator/(const double& i, const std::complex<long double>& z);

std::complex<long double> operator/(const std::complex<long double>& z,const double& i);

std::complex<long double> operator-(const double& i, const std::complex<long double>& z);

std::complex<long double> operator-(const std::complex<long double>& z,const double& i);

std::complex<long double> operator*(const int i, const std::complex<long double>& z);

std::complex<long double> operator+(const int i, const std::complex<long double>& z);

std::complex<long double> operator/(const int i, const std::complex<long double>& z);

std::complex<long double> operator/(const std::complex<long double>& z, const int i);

std::complex<long double> operator-(const int i, const std::complex<long double>& z);

//Overload abs function

template<typename dType, typename Type> dType modulus (const Type& z);

template<> double modulus<double, double> (const double& z);

template<> double modulus<double, std::complex<double> > (const std::complex<double>& z);


template<> float128 modulus<float128, float128> (const float128& z);

template<> float128 modulus<float128, std::complex<float128> > (const std::complex<float128>& z);

template<> mpfr::mpreal modulus<mpfr::mpreal, mpfr::mpreal> (const mpfr::mpreal& z);

template<> mpfr::mpreal modulus<mpfr::mpreal, std::complex<mpfr::mpreal> > (const std::complex<mpfr::mpreal>& z);

template<> long double modulus<long double, long double> (const long double& z);

template<> long double modulus<long double, std::complex<long double> > (const std::complex<long double>& z);

//Overload long double functions

long double sqrt(const long double&);

long double exp(const long double&);

long double sinh(const long double&);

long double cosh(const long double&);

long double fabs(const long double&);

//Implement power laws

template<typename Type1, typename Type2>
Type1 Power(const Type1& a, const Type2 b){
    return pow(a,b);
}

template<>
double Power<double, int>(const double& a, const int b);

template<>
double Power<double, double>(const double& a, const double b);

template<>
long double Power<long double, int>(const long double& a,const int b);

template<>
long double Power<long double, double>(const long double& a, const double b);

template<>
std::complex<float128> Power<std::complex<float128> , int >(const std::complex<float128>& z, const int i);

template<>
std::complex<mpreal> Power< std::complex<mpreal> , int >(const std::complex<mpreal>& z, const int i);

//Implement hyperbolic functions

//template<typename Type> Type Sech(Type z){return 1/cosh(z);}
//
//template<typename Type> Type Csch(Type z){return 1/sinh(z);}
//
//template<typename Type> Type Sinh(Type z){return sinh(z);}
//
//template<typename Type> Type Cosh(Type z){return cosh(z);}
//
//template<typename Type> Type Tanh(Type z){return tanh(z);}
//
//template<typename Type> Type Coth(Type z){return 1/tanh(z);}


//double Sech(double z);
//
//double Csch(double z);
//
//double Sinh(double z);
//
//double Cosh(double z);
//
//double Tanh(double z);
//
//double Coth(double z);

#endif /* defined(___DEM__CustomFunctions__) */

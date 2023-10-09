//
//  CustomFunctions.cpp
//  2DEM
//
//  Created by Aristomenis Donos on 11/12/2014.
//  Copyright (c) 2014 Aristomenis Donos. All rights reserved.
//

#include "CustomFunctions.h"

std::complex<double> operator*(const int i, const std::complex<double>& z){
    return (double)i*z;
}

std::complex<double> operator+(const int i, const std::complex<double>& z){
    return (double)i+z;
}

std::complex<double> operator/(const int i, const std::complex<double>& z){
    return (double)i/z;
}

std::complex<double> operator/(const std::complex<double>& z,const int i){
    return z/(double)i;
}

std::complex<double> operator-(const int i, const std::complex<double>& z){
    return (double)i-z;
}

std::complex<float128> operator*(const double& i,const std::complex<float128>& z){
    return (float128)i*z;
}

std::complex<float128> operator*(const std::complex<float128>& z,const double& i){
    return (float128)i*z;
}

std::complex<float128> operator+(const double& i, const std::complex<float128>& z){
    return (float128)i+z;
}

std::complex<float128> operator+(const std::complex<float128>& z,const double& i){
    return (float128)i+z;
}

std::complex<float128> operator/(const double& i, const std::complex<float128>& z){
    return (float128)i/z;
}

std::complex<float128> operator/(const std::complex<float128>& z,const double& i){
    return z/(float128)i;
}

std::complex<float128> operator-(const double& i, const std::complex<float128>& z){
    return (float128)i-z;
}

std::complex<float128> operator-(const std::complex<float128>& z,const double& i){
    return z-(float128)i;
}

std::complex<mpreal> operator*(const double& i, const std::complex<mpreal>& z){
    return (mpreal)i*z;
}

std::complex<mpreal> operator*(const std::complex<mpreal>& z,const double& i){
    return (mpreal)i*z;
}

std::complex<mpreal> operator+(const double& i, const std::complex<mpreal>& z){
    return (mpreal)i+z;
}

std::complex<mpreal> operator+(const std::complex<mpreal>& z,const double& i){
    return (mpreal)i+z;
}

std::complex<mpreal> operator/(const double& i, const std::complex<mpreal>& z){
    return (mpreal)i/z;
}

std::complex<mpreal> operator/(const std::complex<mpreal> &z,const double& i){
    return z/(mpreal)i;
}

std::complex<mpreal> operator-(const double& i, const std::complex<mpreal>& z){
    return (mpreal)i-z;
}

std::complex<mpreal> operator-(const std::complex<mpreal>& z,const double& i){
    return z-(mpreal)i;
}

std::complex<long double> operator*(const double& i, const std::complex<long double>& z){
    return (long double)i*z;
}

std::complex<long double> operator*(const std::complex<long double>& z,const double& i){
    return (long double)i*z;
}

std::complex<long double> operator+(const double& i, const std::complex<long double>& z){
    return (long double)i+z;
}

std::complex<long double> operator+(const std::complex<long double>& z,const double& i){
    return (long double)i+z;
}

std::complex<long double> operator/(const double& i, const std::complex<long double>& z){
    return (long double)i/z;
}

std::complex<long double> operator/(const std::complex<long double>& z,const double& i){
    return z/(long double)i;
}

std::complex<long double> operator-(const double& i, const std::complex<long double>& z){
    return (long double)i-z;
}

std::complex<long double> operator-(const std::complex<long double>& z,const double& i){
    return z-(long double)i;
}

std::complex<long double> operator*(const int i, const std::complex<long double>& z){
    return (long double)i*z;
}

std::complex<long double> operator+(const int i, const std::complex<long double>& z){
    return (long double)i+z;
}

std::complex<long double> operator/(const int i, const std::complex<long double>& z){
    return (long double)i/z;
}

std::complex<long double> operator/(const std::complex<long double>& z, const int i){
    return z/(long double)i;
}

std::complex<long double> operator-(const int i, const std::complex<long double>& z){
    return (long double)i-z;
}

// Overload abs function

template<> double modulus<double, double> (const double& z){
    return fabs(z);
}

template<> double modulus<double, std::complex<double> > (const std::complex<double>& z){
    return abs(z);
}


template<> float128 modulus<float128, float128> (const float128& z){
    return fabs(z);
}

template<> float128 modulus<float128, std::complex<float128> > (const std::complex<float128>& z){
    return abs(z);
}

template<> mpfr::mpreal modulus<mpfr::mpreal, mpfr::mpreal> (const mpfr::mpreal& z){
    return fabs(z);
}

template<> mpfr::mpreal modulus<mpfr::mpreal, std::complex<mpfr::mpreal> > (const std::complex<mpfr::mpreal>& z){
    return abs(z);
}

template<> long double modulus<long double, long double> (const long double& z){
    return fabsl(z);
}

template<> long double modulus<long double, std::complex<long double> > (const std::complex<long double>& z){
    return abs(z);
}

//Overload long double functions

long double sqrt(const long double& x){
    return sqrtl(x);
}

long double exp(const long double& x){
    return expl(x);
}

long double sinh(const long double& x){
    return sinhl(x);
}

long double cosh(const long double& x){
    return coshl(x);
}

long double fabs(const long double& x){
    return fabsl(x);
}

template<>
double Power(const double& a, const int b){
    return pow(a,b);
}

template<>
double Power(const double& a, const double b){
    return pow(a, b);
}

template<>
long double Power(const long double& a, const int b){
    return powl(a,b);
}

template<>
long double Power(const long double& a, const double b){
    return powl(a, (long double)b);
}

template<>
std::complex<float128> Power<std::complex<float128> , int >(const std::complex<float128>& __x, const int __m){
    
    int __n(__m);
    
    std::complex<float128> xr(__x);
    
    if(__n >= 0){
        std::complex<float128> __y(xr);
        __y = (__n%2)?xr:std::complex<float128>(1);
        
        while (__n >>= 1)
        {
            xr = xr * xr;
            if (__n % 2)
                __y = __y * xr;
        }
        
        return __y;}
    else return 1/Power(xr, -__n);
}

template<>
std::complex<mpreal> Power< std::complex<mpreal> , int >(const std::complex<mpreal>& __x, const int __m){
    
    
    int __n(__m);
    
    std::complex<mpreal> xr(__x);
    
    if(__n >= 0){
        std::complex<mpreal> __y(xr);
        __y = (__n%2)?xr:std::complex<mpreal>(1);
        
        while (__n >>= 1)
        {
            xr = xr * xr;
            if (__n % 2)
                __y = __y * xr;
        }
        
        return __y;}
    else return 1/Power(xr, -__n);
}

//Implement hyperbolic functions

//double Sech(double z){return 1/cosh(z);}
//
//double Csch(double z){return 1/sinh(z);}
//
//double Sinh(double z){return sinh(z);}
//
//double Cosh(double z){return cosh(z);}
//
//double Tanh(double z){return tanh(z);}
//
//double Coth(double z){return 1/tanh(z);}
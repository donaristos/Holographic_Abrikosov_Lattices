//
//  EOMS.h
//  PDEs
//
//  Created by Aristomenis Donos on 05/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#ifndef __PDEs__EOMS__
#define __PDEs__EOMS__

#include <iostream>
#include "function.h"
#include "CustomFunctions.h"
#include "MyTypes.h"

template<typename dType, typename Type>
void equations(Type* elem, Cfunction<dType, Type> *F1, grid<dType>& grid1, grid<dType>& grid2, grid<dType>& grid3, dType* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<typename dType, typename Type>
void bconditions(Type* elem, Cfunction<dType, Type> *F1, grid<dType>& grid1, grid<dType>& grid2, grid<dType>& grid3,dType* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<typename dType, typename Type>
void ksi2(Type* elem, Cfunction<dType, Type> *F1, grid<dType>& grid1, grid<dType>& grid2, grid<dType>& grid3, dType* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<typename dType, typename Type>
void psi(Type* elem, Cfunction<dType, Type> *F1, grid<dType>& grid1, grid<dType>& grid2, grid<dType>& grid3, dType* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void equations<double, double>(double* elem, Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3,  double* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void bconditions<double, double>(double* elem, Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3,double* params,  const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void ksi2<double, double>(double* elem, Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void psi<double, double>(double* elem, Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void equations<float128, float128>(float128* elem, Cfunction<float128, float128> *F1, grid<float128>& grid1, grid<float128>& grid2, grid<float128>& grid3,  float128* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void bconditions<float128, float128>(float128* elem, Cfunction<float128, float128> *F1, grid<float128>& grid1, grid<float128>& grid2, grid<float128>& grid3,float128* params,  const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void ksi2<float128, float128>(float128* elem, Cfunction<float128, float128> *F1, grid<float128>& grid1, grid<float128>& grid2, grid<float128>& grid3, float128* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void psi<float128, float128>(float128* elem, Cfunction<float128, float128> *F1, grid<float128>& grid1, grid<float128>& grid2, grid<float128>& grid3, float128* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void equations<mpreal, mpreal>(mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1, grid<mpreal>& grid2, grid<mpreal>& grid3,  mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void bconditions<mpreal, mpreal>(mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1, grid<mpreal>& grid2, grid<mpreal>& grid3,mpreal* params,  const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void ksi2<mpreal, mpreal>(mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1, grid<mpreal>& grid2, grid<mpreal>& grid3, mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void psi<mpreal, mpreal>(mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1, grid<mpreal>& grid2, grid<mpreal>& grid3, mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void equations<long double, long double>(long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid1, grid<long double>& grid2, grid<long double>& grid3,  long double* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void bconditions<long double, long double>(long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid1, grid<long double>& grid2, grid<long double>& grid3,long double* params,  const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void ksi2<long double, long double>(long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid1, grid<long double>& grid2, grid<long double>& grid3, long double* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void psi<long double, long double>(long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid1, grid<long double>& grid2, grid<long double>& grid3, long double* params, const My_Int& i, const My_Int& j, const My_Int& k);

#endif /* defined(__PDEs__EOMS__) */

//
//  InitLinearOp.h
//  TimeDependence
//
//  Created by Aristomenis Donos on 17/06/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#ifndef __TimeDependence__InitLinearOp__
#define __TimeDependence__InitLinearOp__

#include <iostream>
#include "function.h"
#include "CustomFunctions.h"
#include "MyTypes.h"


template<typename dType, typename Type>
void BulkLinearOp(Type *elem, derivativeCol<dType>& Der,Cfunction<dType, Type> *F1, grid<dType>& grid1, grid<dType>& grid2, grid<dType>& grid3,dType* params ,const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k);

template<typename dType, typename Type>
void BoundaryLinearOp(Type *elem, derivativeCol<dType>& Der,Cfunction<dType, Type> *F1, grid<dType>& grid1, grid<dType>& grid2, grid<dType>& grid3, Type* params,  const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void BulkLinearOp<double, double>(double *elem, derivativeCol<double>& D,Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3,  double* params, const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1);

template<>
void BoundaryLinearOp<double, double>(double *elem, derivativeCol<double>& D,Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params,  const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1);

template<>
void BulkLinearOp<float128, float128>(float128 *elem, derivativeCol<float128>& D,Cfunction<float128, float128> *F1, grid<float128>& grid1, grid<float128>& grid2, grid<float128>& grid3, float128* params, const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1);

template<>
void BoundaryLinearOp<float128, float128>(float128 *elem, derivativeCol<float128>& D,Cfunction<float128, float128> *F1, grid<float128>& grid1, grid<float128>& grid2, grid<float128>& grid3, float128* params,  const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1);

template<>
void BulkLinearOp<mpreal, mpreal>(mpreal *elem, derivativeCol<mpreal>& D,Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1, grid<mpreal>& grid2, grid<mpreal>& grid3, mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1);

template<>
void BoundaryLinearOp<mpreal, mpreal>(mpreal *elem, derivativeCol<mpreal>& D,Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1, grid<mpreal>& grid2, grid<mpreal>& grid3, mpreal* params,  const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1);

template<>
void BulkLinearOp<long double, long double>(long double *elem, derivativeCol<long double>& D,Cfunction<long double, long double> *F1, grid<long double>& grid1, grid<long double>& grid2, grid<long double>& grid3, long double* params, const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1);

template<>
void BoundaryLinearOp<long double, long double>(long double *elem, derivativeCol<long double>& D,Cfunction<long double, long double> *F1, grid<long double>& grid1, grid<long double>& grid2, grid<long double>& grid3, long double* params,  const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1);

#endif /* defined(__TimeDependence__InitLinearOp__) */

//
//  InitLinearOp.cpp
//  TimeDependence
//
//  Created by Aristomenis Donos on 17/06/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#include "LatticeLinearOp.h"

//using namespce mpfr;

template<>
void BulkLinearOp<float128, float128>(float128 *elem, derivativeCol<float128>& D,Cfunction<float128, float128> *F1, grid<float128>& grid1,  float128* params, const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1){
    Cfunction<float128,float128> &Qtt=F1[0];
    Cfunction<float128,float128> &Qrr=F1[1];
    Cfunction<float128,float128> &Q=F1[2];
    Cfunction<float128,float128> &B=F1[3];
    Cfunction<float128,float128> &R=F1[4];
    Cfunction<float128,float128> &Qr1=F1[5];
    Cfunction<float128,float128> &Qr2=F1[6];
    Cfunction<float128,float128> &Qt1=F1[7];
    Cfunction<float128,float128> &Qt2=F1[8];
    Cfunction<float128,float128> &Qtr=F1[9];
    Cfunction<float128,float128> &a0=F1[10];
    Cfunction<float128,float128> &a1=F1[11];
    Cfunction<float128,float128> &a2=F1[12];
    Cfunction<float128,float128> &ar=F1[13];
    Cfunction<float128,float128> &h=F1[14];
    grid<float128>& rgrid=grid1;
    float128 mu=params[0];
    float128 rh=params[1];
    float128 H=params[2];
    float128 R0=params[3];
    float128 L1=params[4];
    float128 L2=params[5];
    float128 S=params[6];
    float128 X=params[7];
    My_Int d100(1),d200(2),d110(3),d010(4),d020(5),d001(6),d002(7),d101(8),d011(9);
    
#include "SourceFiles/LO.c"
}

template<>
void BoundaryLinearOp<float128, float128>(float128 *elem, derivativeCol<float128>& D,Cfunction<float128, float128> *F1, float128* params,  const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1){
    Cfunction<float128,float128> &Qtt=F1[0];
    Cfunction<float128,float128> &Qrr=F1[1];
    Cfunction<float128,float128> &Q=F1[2];
    Cfunction<float128,float128> &B=F1[3];
    Cfunction<float128,float128> &R=F1[4];
    Cfunction<float128,float128> &Qr1=F1[5];
    Cfunction<float128,float128> &Qr2=F1[6];
    Cfunction<float128,float128> &Qt1=F1[7];
    Cfunction<float128,float128> &Qt2=F1[8];
    Cfunction<float128,float128> &Qtr=F1[9];
    Cfunction<float128,float128> &a0=F1[10];
    Cfunction<float128,float128> &a1=F1[11];
    Cfunction<float128,float128> &a2=F1[12];
    Cfunction<float128,float128> &ar=F1[13];
    Cfunction<float128,float128> &h=F1[14];
    float128 mu=params[0];
    float128 rh=params[1];
    float128 H=params[2];
    float128 R0=params[3];
    float128 L1=params[4];
    float128 L2=params[5];
    float128 S=params[6];
    float128 X=params[7];
    My_Int d100(1),d200(2),d110(3),d010(4),d020(5),d001(6),d002(7),d101(8),d011(9);
    
#include "SourceFiles/BCLO.c"
    
}
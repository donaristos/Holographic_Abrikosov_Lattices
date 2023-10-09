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
void BulkLinearOp<mpreal, mpreal>(mpreal *elem, derivativeCol<mpreal>& D,Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1,  mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1){
    Cfunction<mpreal,mpreal> &Qtt=F1[0];
    Cfunction<mpreal,mpreal> &Qrr=F1[1];
    Cfunction<mpreal,mpreal> &Q=F1[2];
    Cfunction<mpreal,mpreal> &B=F1[3];
    Cfunction<mpreal,mpreal> &R=F1[4];
    Cfunction<mpreal,mpreal> &Qr1=F1[5];
    Cfunction<mpreal,mpreal> &Qr2=F1[6];
    Cfunction<mpreal,mpreal> &Qt1=F1[7];
    Cfunction<mpreal,mpreal> &Qt2=F1[8];
    Cfunction<mpreal,mpreal> &Qtr=F1[9];
    Cfunction<mpreal,mpreal> &a0=F1[10];
    Cfunction<mpreal,mpreal> &a1=F1[11];
    Cfunction<mpreal,mpreal> &a2=F1[12];
    Cfunction<mpreal,mpreal> &ar=F1[13];
    Cfunction<mpreal,mpreal> &h=F1[14];
    grid<mpreal>& rgrid=grid1;
    mpreal mu=params[0];
    mpreal rh=params[1];
    mpreal H=params[2];
    mpreal R0=params[3];
    mpreal L1=params[4];
    mpreal L2=params[5];
    mpreal S=params[6];
    mpreal X=params[7];
    My_Int d100(1),d200(2),d110(3),d010(4),d020(5),d001(6),d002(7),d101(8),d011(9);
    
#include "SourceFiles/LO.c"
}

template<>
void BoundaryLinearOp<mpreal, mpreal>(mpreal *elem, derivativeCol<mpreal>& D,Cfunction<mpreal, mpreal> *F1, mpreal* params,  const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1){
    Cfunction<mpreal,mpreal> &Qtt=F1[0];
    Cfunction<mpreal,mpreal> &Qrr=F1[1];
    Cfunction<mpreal,mpreal> &Q=F1[2];
    Cfunction<mpreal,mpreal> &B=F1[3];
    Cfunction<mpreal,mpreal> &R=F1[4];
    Cfunction<mpreal,mpreal> &Qr1=F1[5];
    Cfunction<mpreal,mpreal> &Qr2=F1[6];
    Cfunction<mpreal,mpreal> &Qt1=F1[7];
    Cfunction<mpreal,mpreal> &Qt2=F1[8];
    Cfunction<mpreal,mpreal> &Qtr=F1[9];
    Cfunction<mpreal,mpreal> &a0=F1[10];
    Cfunction<mpreal,mpreal> &a1=F1[11];
    Cfunction<mpreal,mpreal> &a2=F1[12];
    Cfunction<mpreal,mpreal> &ar=F1[13];
    Cfunction<mpreal,mpreal> &h=F1[14];
    mpreal mu=params[0];
    mpreal rh=params[1];
    mpreal H=params[2];
    mpreal R0=params[3];
    mpreal L1=params[4];
    mpreal L2=params[5];
    mpreal S=params[6];
    mpreal X=params[7];
    My_Int d100(1),d200(2),d110(3),d010(4),d020(5),d001(6),d002(7),d101(8),d011(9);
    
#include "SourceFiles/BCLO.c"
    
}
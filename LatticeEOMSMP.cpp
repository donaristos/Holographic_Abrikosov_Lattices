//
//  EOMS.cpp
//  PDEs
//
//  Created by Aristomenis Donos on 05/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#include "LatticeEOMS.h"

//using namespace mpfr;


template<>
void equations<mpreal, mpreal>(mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1,  mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k){
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
    
#include "SourceFiles/EOM.c"
}

template<>
void bconditions<mpreal, mpreal>(mpreal* elem, Cfunction<mpreal, mpreal> *F1,mpreal* params,  const My_Int& i, const My_Int& j, const My_Int& k){
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
    
#include "SourceFiles/BC.c"
}

template<>
void ksi2<mpreal, mpreal>(mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1, mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k){
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
    
#include "SourceFiles/Ksi2.c"
}

template<>
void psi<mpreal, mpreal>(mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1, mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k){
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
    
#include "SourceFiles/Psi.c"
}
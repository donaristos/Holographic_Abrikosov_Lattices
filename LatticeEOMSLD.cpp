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
void equations<long double, long double>(long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid1,  long double* params, const My_Int& i, const My_Int& j, const My_Int& k){
    Cfunction<long double,long double> &Qtt=F1[0];
    Cfunction<long double,long double> &Qrr=F1[1];
    Cfunction<long double,long double> &Q=F1[2];
    Cfunction<long double,long double> &B=F1[3];
    Cfunction<long double,long double> &R=F1[4];
    Cfunction<long double,long double> &Qr1=F1[5];
    Cfunction<long double,long double> &Qr2=F1[6];
    Cfunction<long double,long double> &Qt1=F1[7];
    Cfunction<long double,long double> &Qt2=F1[8];
    Cfunction<long double,long double> &Qtr=F1[9];
    Cfunction<long double,long double> &a0=F1[10];
    Cfunction<long double,long double> &a1=F1[11];
    Cfunction<long double,long double> &a2=F1[12];
    Cfunction<long double,long double> &ar=F1[13];
    Cfunction<long double,long double> &h=F1[14];
    grid<long double>& rgrid=grid1;
    long double mu=params[0];
    long double rh=params[1];
    long double H=params[2];
    long double R0=params[3];
    long double L1=params[4];
    long double L2=params[5];
    long double S=params[6];
    long double X=params[7];
    My_Int d100(1),d200(2),d110(3),d010(4),d020(5),d001(6),d002(7),d101(8),d011(9);
    
#include "SourceFiles/EOM.c"
}

template<>
void bconditions<long double, long double>(long double* elem, Cfunction<long double, long double> *F1,long double* params,  const My_Int& i, const My_Int& j, const My_Int& k){
    Cfunction<long double,long double> &Qtt=F1[0];
    Cfunction<long double,long double> &Qrr=F1[1];
    Cfunction<long double,long double> &Q=F1[2];
    Cfunction<long double,long double> &B=F1[3];
    Cfunction<long double,long double> &R=F1[4];
    Cfunction<long double,long double> &Qr1=F1[5];
    Cfunction<long double,long double> &Qr2=F1[6];
    Cfunction<long double,long double> &Qt1=F1[7];
    Cfunction<long double,long double> &Qt2=F1[8];
    Cfunction<long double,long double> &Qtr=F1[9];
    Cfunction<long double,long double> &a0=F1[10];
    Cfunction<long double,long double> &a1=F1[11];
    Cfunction<long double,long double> &a2=F1[12];
    Cfunction<long double,long double> &ar=F1[13];
    Cfunction<long double,long double> &h=F1[14];
    long double mu=params[0];
    long double rh=params[1];
    long double H=params[2];
    long double R0=params[3];
    long double L1=params[4];
    long double L2=params[5];
    long double S=params[6];
    long double X=params[7];
    My_Int d100(1),d200(2),d110(3),d010(4),d020(5),d001(6),d002(7),d101(8),d011(9);
    
#include "SourceFiles/BC.c"
}

template<>
void ksi2<long double, long double>(long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid1, long double* params, const My_Int& i, const My_Int& j, const My_Int& k){
    Cfunction<long double,long double> &Qtt=F1[0];
    Cfunction<long double,long double> &Qrr=F1[1];
    Cfunction<long double,long double> &Q=F1[2];
    Cfunction<long double,long double> &B=F1[3];
    Cfunction<long double,long double> &R=F1[4];
    Cfunction<long double,long double> &Qr1=F1[5];
    Cfunction<long double,long double> &Qr2=F1[6];
    Cfunction<long double,long double> &Qt1=F1[7];
    Cfunction<long double,long double> &Qt2=F1[8];
    Cfunction<long double,long double> &Qtr=F1[9];
    Cfunction<long double,long double> &a0=F1[10];
    Cfunction<long double,long double> &a1=F1[11];
    Cfunction<long double,long double> &a2=F1[12];
    Cfunction<long double,long double> &ar=F1[13];
    Cfunction<long double,long double> &h=F1[14];
    grid<long double>& rgrid=grid1;
    long double mu=params[0];
    long double rh=params[1];
    long double H=params[2];
    long double R0=params[3];
    long double L1=params[4];
    long double L2=params[5];
    long double S=params[6];
    long double X=params[7];
    My_Int d100(1),d200(2),d110(3),d010(4),d020(5),d001(6),d002(7),d101(8),d011(9);
    
#include "SourceFiles/Ksi2.c"
}

template<>
void psi<long double, long double>(long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid1, long double* params, const My_Int& i, const My_Int& j, const My_Int& k){
    Cfunction<long double,long double> &Qtt=F1[0];
    Cfunction<long double,long double> &Qrr=F1[1];
    Cfunction<long double,long double> &Q=F1[2];
    Cfunction<long double,long double> &B=F1[3];
    Cfunction<long double,long double> &R=F1[4];
    Cfunction<long double,long double> &Qr1=F1[5];
    Cfunction<long double,long double> &Qr2=F1[6];
    Cfunction<long double,long double> &Qt1=F1[7];
    Cfunction<long double,long double> &Qt2=F1[8];
    Cfunction<long double,long double> &Qtr=F1[9];
    Cfunction<long double,long double> &a0=F1[10];
    Cfunction<long double,long double> &a1=F1[11];
    Cfunction<long double,long double> &a2=F1[12];
    Cfunction<long double,long double> &ar=F1[13];
    Cfunction<long double,long double> &h=F1[14];
    grid<long double>& rgrid=grid1;
    long double mu=params[0];
    long double rh=params[1];
    long double H=params[2];
    long double R0=params[3];
    long double L1=params[4];
    long double L2=params[5];
    long double S=params[6];
    long double X=params[7];
    My_Int d100(1),d200(2),d110(3),d010(4),d020(5),d001(6),d002(7),d101(8),d011(9);
    
#include "SourceFiles/Psi.c"
}
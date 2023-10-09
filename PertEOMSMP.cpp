//
//  TimeEvoEOMS.cpp
//  TimeDependence
//
//  Created by Aristomenis Donos on 13/05/2.013.
//  Copyright (c) 2.013 Aristomenis Donos. All rights reserved.
//

#include "PertEOMS.h"


template<>
void equations<mpreal, ComplexMP>(std::complex<mpreal>* elem, Cfunction<mpreal, std::complex<mpreal> > *F1, grid<mpreal>& grid1,  mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k){
	Cfunction<mpreal, std::complex<mpreal> > &Ax=F1[0];
	Cfunction<mpreal, std::complex<mpreal> > &h11=F1[1];
	Cfunction<mpreal, std::complex<mpreal> > &Dh=F1[2];
	Cfunction<mpreal, std::complex<mpreal> > &htt=F1[3];
	Cfunction<mpreal, std::complex<mpreal> > &h22=F1[4];
	Cfunction<mpreal, std::complex<mpreal> > &ht1=F1[5];
	Cfunction<mpreal, std::complex<mpreal> > &A0=F1[6];
	Cfunction<mpreal, std::complex<mpreal> > &Qtt=F1[7];
	Cfunction<mpreal, std::complex<mpreal> > &Qrr=F1[8];
	Cfunction<mpreal, std::complex<mpreal> > &Q11=F1[9];
	Cfunction<mpreal, std::complex<mpreal> > &Qr1=F1[10];
	Cfunction<mpreal, std::complex<mpreal> > &Q22=F1[11];
	Cfunction<mpreal, std::complex<mpreal> > &a0=F1[12];
	Cfunction<mpreal, std::complex<mpreal> > &h=F1[13];
	typedef std::complex<mpreal> complex;
	grid<mpreal>& rgrid=grid1;
	mpreal mu=params[0];
	mpreal L=params[1]/mu;
	mpreal w=params[7];
	mpreal muj=1;
	mpreal &B1=params[6];
	mpreal &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/PEOM.c"
	
}


template<>
void bconditions<mpreal, ComplexMP>(std::complex<mpreal>* elem, Cfunction<mpreal, std::complex<mpreal> > *F1,mpreal* params,  const My_Int& i, const My_Int& j, const My_Int& k){
	Cfunction<mpreal, std::complex<mpreal> > &Ax=F1[0];
	Cfunction<mpreal, std::complex<mpreal> > &h11=F1[1];
	Cfunction<mpreal, std::complex<mpreal> > &Dh=F1[2];
	Cfunction<mpreal, std::complex<mpreal> > &htt=F1[3];
	Cfunction<mpreal, std::complex<mpreal> > &h22=F1[4];
	Cfunction<mpreal, std::complex<mpreal> > &ht1=F1[5];
	Cfunction<mpreal, std::complex<mpreal> > &A0=F1[6];
	Cfunction<mpreal, std::complex<mpreal> > &Qtt=F1[7];
	Cfunction<mpreal, std::complex<mpreal> > &Qrr=F1[8];
	Cfunction<mpreal, std::complex<mpreal> > &Q11=F1[9];
	Cfunction<mpreal, std::complex<mpreal> > &Qr1=F1[10];
	Cfunction<mpreal, std::complex<mpreal> > &Q22=F1[11];
	Cfunction<mpreal, std::complex<mpreal> > &a0=F1[12];
	Cfunction<mpreal, std::complex<mpreal> > &h=F1[13];
	typedef std::complex<mpreal> complex;
	typedef mpreal lType;
	mpreal mu=params[0];
	mpreal L=params[1]/mu;
	mpreal w=params[7];
	mpreal muj=1;
	mpreal &B1=params[6];
	mpreal &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5), d12(6);
	
#include "SourceFiles/PBC.c"
#include "SourceFiles/PAsBC.c"
	
	
}


template<>
void SecOrderEqs<mpreal, ComplexMP>(std::complex<mpreal>* elem, Cfunction<mpreal, std::complex<mpreal> > *F1, grid<mpreal>& grid1,  mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k){
	Cfunction<mpreal, std::complex<mpreal> > &Ax=F1[0];
	Cfunction<mpreal, std::complex<mpreal> > &h11=F1[1];
	Cfunction<mpreal, std::complex<mpreal> > &Dh=F1[2];
	Cfunction<mpreal, std::complex<mpreal> > &htt=F1[3];
	Cfunction<mpreal, std::complex<mpreal> > &h22=F1[4];
	Cfunction<mpreal, std::complex<mpreal> > &ht1=F1[5];
	Cfunction<mpreal, std::complex<mpreal> > &A0=F1[6];
	Cfunction<mpreal, std::complex<mpreal> > &Qtt=F1[7];
	Cfunction<mpreal, std::complex<mpreal> > &Qrr=F1[8];
	Cfunction<mpreal, std::complex<mpreal> > &Q11=F1[9];
	Cfunction<mpreal, std::complex<mpreal> > &Qr1=F1[10];
	Cfunction<mpreal, std::complex<mpreal> > &Q22=F1[11];
	Cfunction<mpreal, std::complex<mpreal> > &a0=F1[12];
	Cfunction<mpreal, std::complex<mpreal> > &h=F1[13];
	typedef std::complex<mpreal> complex;
	grid<mpreal>& rgrid=grid1;
	mpreal mu=params[0];
	mpreal L=params[1]/mu;
	mpreal w=params[7];
	mpreal muj=1;
	mpreal &B1=params[6];
	mpreal &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/PEOMSO.c"
}
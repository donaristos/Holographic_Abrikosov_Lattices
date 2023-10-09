//
//  TimeEvoEOMS.cpp
//  TimeDependence
//
//  Created by Aristomenis Donos on 13/05/2.013.
//  Copyright (c) 2.013 Aristomenis Donos. All rights reserved.
//

#include "PertEOMS.h"

template<>
void equations<float128, ComplexQ>(std::complex<float128>* elem, Cfunction<float128, std::complex<float128> > *F1, grid<float128>& grid1,  float128* params, const My_Int& i, const My_Int& j, const My_Int& k){
	Cfunction<float128, std::complex<float128> > &Ax=F1[0];
	Cfunction<float128, std::complex<float128> > &h11=F1[1];
	Cfunction<float128, std::complex<float128> > &Dh=F1[2];
	Cfunction<float128, std::complex<float128> > &htt=F1[3];
	Cfunction<float128, std::complex<float128> > &h22=F1[4];
	Cfunction<float128, std::complex<float128> > &ht1=F1[5];
	Cfunction<float128, std::complex<float128> > &A0=F1[6];
	Cfunction<float128, std::complex<float128> > &Qtt=F1[7];
	Cfunction<float128, std::complex<float128> > &Qrr=F1[8];
	Cfunction<float128, std::complex<float128> > &Q11=F1[9];
	Cfunction<float128, std::complex<float128> > &Qr1=F1[10];
	Cfunction<float128, std::complex<float128> > &Q22=F1[11];
	Cfunction<float128, std::complex<float128> > &a0=F1[12];
	Cfunction<float128, std::complex<float128> > &h=F1[13];
	typedef std::complex<float128> complex;
	grid<float128>& rgrid=grid1;
	float128 mu=params[0];
	float128 L=params[1]/mu;
	float128 w=params[7];
	float128 muj=1;
	float128 &B1=params[6];
	float128 &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/PEOM.c"
	
}


template<>
void bconditions<float128, ComplexQ>(std::complex<float128>* elem, Cfunction<float128, std::complex<float128> > *F1,float128* params,  const My_Int& i, const My_Int& j, const My_Int& k){
	Cfunction<float128, std::complex<float128> > &Ax=F1[0];
	Cfunction<float128, std::complex<float128> > &h11=F1[1];
	Cfunction<float128, std::complex<float128> > &Dh=F1[2];
	Cfunction<float128, std::complex<float128> > &htt=F1[3];
	Cfunction<float128, std::complex<float128> > &h22=F1[4];
	Cfunction<float128, std::complex<float128> > &ht1=F1[5];
	Cfunction<float128, std::complex<float128> > &A0=F1[6];
	Cfunction<float128, std::complex<float128> > &Qtt=F1[7];
	Cfunction<float128, std::complex<float128> > &Qrr=F1[8];
	Cfunction<float128, std::complex<float128> > &Q11=F1[9];
	Cfunction<float128, std::complex<float128> > &Qr1=F1[10];
	Cfunction<float128, std::complex<float128> > &Q22=F1[11];
	Cfunction<float128, std::complex<float128> > &a0=F1[12];
	Cfunction<float128, std::complex<float128> > &h=F1[13];
	typedef std::complex<float128> complex;
	typedef float128 lType;
	float128 mu=params[0];
	float128 L=params[1]/mu;
	float128 w=params[7];
	float128 muj=1;
	float128 &B1=params[6];
	float128 &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5), d12(6);
	
#include "SourceFiles/PBC.c"
#include "SourceFiles/PAsBC.c"
	
	
}


template<>
void SecOrderEqs<float128, ComplexQ>(std::complex<float128>* elem, Cfunction<float128, std::complex<float128> > *F1, grid<float128>& grid1,  float128* params, const My_Int& i, const My_Int& j, const My_Int& k){
	Cfunction<float128, std::complex<float128> > &Ax=F1[0];
	Cfunction<float128, std::complex<float128> > &h11=F1[1];
	Cfunction<float128, std::complex<float128> > &Dh=F1[2];
	Cfunction<float128, std::complex<float128> > &htt=F1[3];
	Cfunction<float128, std::complex<float128> > &h22=F1[4];
	Cfunction<float128, std::complex<float128> > &ht1=F1[5];
	Cfunction<float128, std::complex<float128> > &A0=F1[6];
	Cfunction<float128, std::complex<float128> > &Qtt=F1[7];
	Cfunction<float128, std::complex<float128> > &Qrr=F1[8];
	Cfunction<float128, std::complex<float128> > &Q11=F1[9];
	Cfunction<float128, std::complex<float128> > &Qr1=F1[10];
	Cfunction<float128, std::complex<float128> > &Q22=F1[11];
	Cfunction<float128, std::complex<float128> > &a0=F1[12];
	Cfunction<float128, std::complex<float128> > &h=F1[13];
	typedef std::complex<float128> complex;
	grid<float128>& rgrid=grid1;
	float128 mu=params[0];
	float128 L=params[1]/mu;
	float128 w=params[7];
	float128 muj=1;
	float128 &B1=params[6];
	float128 &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/PEOMSO.c"
}
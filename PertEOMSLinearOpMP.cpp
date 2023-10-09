//
//  TimeEvoLinearOp.cpp
//  TimeDependence
//
//  Created by Aristomenis Donos on 17/06/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#include "PertEOMSLinearOp.h"

template<>
void BulkLinearOp<mpreal, ComplexMP>(ComplexMP **elem, derivativeCol<mpreal>& D,Cfunction<mpreal, std::complex<mpreal> > *F1, grid<mpreal>& grid1,  mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1){
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
	
#include "SourceFiles/PLO.c"
	
}

template<>
void BoundaryLinearOp<mpreal, ComplexMP>(ComplexMP **elem, derivativeCol<mpreal>& D,Cfunction<mpreal, std::complex<mpreal> > *F1, mpreal* params,  const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1){
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
	mpreal mu=params[0];
	mpreal L=params[1]/mu;
	mpreal w=params[7];
	mpreal muj=1;
	mpreal &B1=params[6];
	mpreal &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5), d12(6);
	
#include "SourceFiles/PBCLO.c"
#include "SourceFiles/PAsBCLO.c"
	
}
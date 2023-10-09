//
//  TimeEvoLinearOp.cpp
//  TimeDependence
//
//  Created by Aristomenis Donos on 17/06/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#include "PertEOMSLinearOp.h"

template<>
void BulkLinearOp<long double, ComplexLD>(ComplexLD **elem, derivativeCol<long double>& D,Cfunction<long double, std::complex<long double> > *F1, grid<long double>& grid1,  long double* params, const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1){
	Cfunction<long double, std::complex<long double> > &Ax=F1[0];
	Cfunction<long double, std::complex<long double> > &h11=F1[1];
	Cfunction<long double, std::complex<long double> > &Dh=F1[2];
	Cfunction<long double, std::complex<long double> > &htt=F1[3];
	Cfunction<long double, std::complex<long double> > &h22=F1[4];
	Cfunction<long double, std::complex<long double> > &ht1=F1[5];
	Cfunction<long double, std::complex<long double> > &A0=F1[6];
	Cfunction<long double, std::complex<long double> > &Qtt=F1[7];
	Cfunction<long double, std::complex<long double> > &Qrr=F1[8];
	Cfunction<long double, std::complex<long double> > &Q11=F1[9];
	Cfunction<long double, std::complex<long double> > &Qr1=F1[10];
	Cfunction<long double, std::complex<long double> > &Q22=F1[11];
	Cfunction<long double, std::complex<long double> > &a0=F1[12];
	Cfunction<long double, std::complex<long double> > &h=F1[13];
	typedef std::complex<long double> complex;
	grid<long double>& rgrid=grid1;
	long double mu=params[0];
	long double L=params[1]/mu;
	long double w=params[7];
	long double muj=1;
	long double &B1=params[6];
	long double &A1=params[5];
    My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/PLO.c"
	
}

template<>
void BoundaryLinearOp<long double, ComplexLD>(ComplexLD **elem, derivativeCol<long double>& D,Cfunction<long double, std::complex<long double> > *F1, long double* params,  const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1){
	Cfunction<long double, std::complex<long double> > &Ax=F1[0];
	Cfunction<long double, std::complex<long double> > &h11=F1[1];
	Cfunction<long double, std::complex<long double> > &Dh=F1[2];
	Cfunction<long double, std::complex<long double> > &htt=F1[3];
	Cfunction<long double, std::complex<long double> > &h22=F1[4];
	Cfunction<long double, std::complex<long double> > &ht1=F1[5];
	Cfunction<long double, std::complex<long double> > &A0=F1[6];
	Cfunction<long double, std::complex<long double> > &Qtt=F1[7];
	Cfunction<long double, std::complex<long double> > &Qrr=F1[8];
	Cfunction<long double, std::complex<long double> > &Q11=F1[9];
	Cfunction<long double, std::complex<long double> > &Qr1=F1[10];
	Cfunction<long double, std::complex<long double> > &Q22=F1[11];
	Cfunction<long double, std::complex<long double> > &a0=F1[12];
	Cfunction<long double, std::complex<long double> > &h=F1[13];
	typedef std::complex<long double> complex;
	long double mu=params[0];
	long double L=params[1]/mu;
	long double w=params[7];
	long double muj=1;
	long double &B1=params[6];
	long double &A1=params[5];
    My_Int d10(1),d20(2),d11(3),d01(4),d02(5), d12(6);
	
#include "SourceFiles/PBCLO.c"
#include "SourceFiles/PAsBCLO.c"
	
}
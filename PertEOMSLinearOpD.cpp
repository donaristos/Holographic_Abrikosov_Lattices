//
//  TimeEvoLinearOp.cpp
//  TimeDependence
//
//  Created by Aristomenis Donos on 17/06/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#include "PertEOMSLinearOp.h"

template<>
void BulkLinearOp<double, Complex >(Complex **elem, derivativeCol<double>& D,Cfunction<double, Complex > *F1, grid<double>& grid1, double* params, const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1)
{   Cfunction<double, std::complex<double> > &Ax=F1[0];
	Cfunction<double, std::complex<double> > &h11=F1[1];
	Cfunction<double, std::complex<double> > &Dh=F1[2];
	Cfunction<double, std::complex<double> > &htt=F1[3];
	Cfunction<double, std::complex<double> > &h22=F1[4];
	Cfunction<double, std::complex<double> > &ht1=F1[5];
	Cfunction<double, std::complex<double> > &A0=F1[6];
	Cfunction<double, std::complex<double> > &Qtt=F1[7];
	Cfunction<double, std::complex<double> > &Qrr=F1[8];
	Cfunction<double, std::complex<double> > &Q11=F1[9];
	Cfunction<double, std::complex<double> > &Qr1=F1[10];
	Cfunction<double, std::complex<double> > &Q22=F1[11];
	Cfunction<double, std::complex<double> > &a0=F1[12];
	Cfunction<double, std::complex<double> > &h=F1[13];
	typedef std::complex<double> complex;
	grid<double>& rgrid=grid1;
	double mu=params[0];
	double L=params[1]/mu;
	double w=params[7];
	double muj=1.;
	double &B1=params[6];
	double &A1=params[5];
    My_Int d10(1),d20(2),d11(3),d01(4),d02(5);

#include "SourceFiles/PLO.c"

}

template<>
void BoundaryLinearOp<double, Complex>(Complex **elem, derivativeCol<double>& D,Cfunction<double, Complex > *F1, double* params,  const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1)
{   Cfunction<double, std::complex<double> > &Ax=F1[0];
	Cfunction<double, std::complex<double> > &h11=F1[1];
	Cfunction<double, std::complex<double> > &Dh=F1[2];
	Cfunction<double, std::complex<double> > &htt=F1[3];
	Cfunction<double, std::complex<double> > &h22=F1[4];
	Cfunction<double, std::complex<double> > &ht1=F1[5];
	Cfunction<double, std::complex<double> > &A0=F1[6];
	Cfunction<double, std::complex<double> > &Qtt=F1[7];
	Cfunction<double, std::complex<double> > &Qrr=F1[8];
	Cfunction<double, std::complex<double> > &Q11=F1[9];
	Cfunction<double, std::complex<double> > &Qr1=F1[10];
	Cfunction<double, std::complex<double> > &Q22=F1[11];
	Cfunction<double, std::complex<double> > &a0=F1[12];
	Cfunction<double, std::complex<double> > &h=F1[13];
	typedef std::complex<double> complex;
	double mu=params[0];
	double L=params[1]/mu;
	double w=params[7];
	double muj=1.;
	double &B1=params[6];
	double &A1=params[5];
    My_Int d10(1),d20(2),d11(3),d01(4),d02(5), d12(6);

#include "SourceFiles/PBCLO.c"
#include "SourceFiles/PAsBCLO.c"

}
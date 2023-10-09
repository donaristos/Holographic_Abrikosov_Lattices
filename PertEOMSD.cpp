//
//  TimeEvoEOMS.cpp
//  TimeDependence
//
//  Created by Aristomenis Donos on 13/05/2.013.
//  Copyright (c) 2.013 Aristomenis Donos. All rights reserved.
//

#include "PertEOMS.h"

template<>
void equations<double, Complex>(std::complex<double>* elem, Cfunction<double, std::complex<double> > *F1, grid<double>& grid1,  double* params, const My_Int& i, const My_Int& j, const My_Int& k){
	Cfunction<double, std::complex<double> > &Ax=F1[0];
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
	double muj=1.L;
	double &B1=params[6];
	double &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);

#include "SourceFiles/PEOM.c"

}


template<>
void bconditions<double, Complex>(std::complex<double>* elem, Cfunction<double, std::complex<double> > *F1,double* params,  const My_Int& i, const My_Int& j, const My_Int& k){
	Cfunction<double, std::complex<double> > &Ax=F1[0];
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
	typedef double lType;
	double mu=params[0];
	double L=params[1]/mu;
	double w=params[7];
	double muj=1.L;
	double &B1=params[6];
	double &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5), d12(6);

#include "SourceFiles/PBC.c"
#include "SourceFiles/PAsBC.c"


}


template<>
void SecOrderEqs<double, Complex>(std::complex<double>* elem, Cfunction<double, std::complex<double> > *F1, grid<double>& grid1,  double* params, const My_Int& i, const My_Int& j, const My_Int& k){
	Cfunction<double, std::complex<double> > &Ax=F1[0];
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
	double muj=1.L;
	double &B1=params[6];
	double &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);

	#include "SourceFiles/PEOMSO.c"
}
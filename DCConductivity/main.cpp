#include <iostream>
#include <fstream>
#include <string>
//#define EIGEN_NO_DEBUG
#define _USE_MATH_DEFINES

#include "../function.h"
#include "../NewtonMethod.h"
#include "../ReadWrite.h"
#include <ctime>

#include "../PertEOMS.h"
#include "../PertEOMSLinearOp.h"

#include "../MyTypes.h"

#ifdef F128
typedef float128 MyType;
#elif MP
typedef mpreal MyType;
#elif LONGD
typedef long double MyType;
#else
typedef double MyType;
#endif

typedef std::complex<MyType> MyComplex;

MyType **cheby_d1;
MyType **cheby_d2;
MyType **cheby_d3;
MyType **fourier_d1;
MyType **fourier_d2;
MyType **finitediff1;
MyType **finitediff2;
MyType **idr=0;
MyType **idx=0;

clock_t start1,start2;
double omp_start, omp_end;
MyType diff;


int main(int argc, const char * argv[])
{
    mpfr::mpreal::set_default_prec(76);
	
	My_Int* Rdimensions;
	My_Int* ROrders;
	My_Int* Xdimensions;
	
	My_Int rlength=GetLength<My_Int>(1, "Rdimensions.dat");
	
	Rdimensions= new My_Int[rlength];
	ROrders= new My_Int[rlength];
	Xdimensions= new My_Int[1];
	
	ReadIntParameters(Rdimensions, rlength, "Rdimensions.dat");
	ReadIntParameters(ROrders, rlength, "ROrders.dat");
	ReadIntParameters(Xdimensions, 1, "Xdimensions.dat");
	
	My_Int tnp(0);
	My_Int tnx=Xdimensions[0];
	
	for(My_Int i=0;i<rlength;++i){
		tnp+=Rdimensions[i];
	}
	
	std::cout << "Np= " << tnp << "\t, Nx= " << tnx << std::endl;
	
	My_Int XOrders[]={Xdimensions[0]};
	
	MyType* trgrid;
	MyType* txgrid;
	
	trgrid= new MyType[tnp];
	txgrid= new MyType[tnx];
	
    //	setxgrid(txgrid, tnx);
	
	ReadDoubleParameters(trgrid, tnp, "rgrid.dat");
	ReadDoubleParameters(txgrid, tnx, "xgrid.dat");
	
	//Define the Chebyshev differentiation matrix and its square
	
	cheby_d1= new MyType*[tnp];
	cheby_d2= new MyType*[tnp];
	for (My_Int i1=0; i1<tnp;i1++){cheby_d1[i1]=new MyType[tnp];cheby_d2[i1]=new MyType[tnp];}
	
	//Define the Fourier differenations matrix and its square
	
	fourier_d1= new MyType*[tnx];
	fourier_d2= new MyType*[tnx];
	for (My_Int i1=0; i1<tnx ;i1++){fourier_d1[i1]=new MyType[tnx];fourier_d2[i1]=new MyType[tnx];}
	
	ReadDiffFromBin<MyType>(cheby_d1,  tnp, "cheby_d1.dat");
	ReadDiffFromBin<MyType>(cheby_d2,  tnp, "cheby_d2.dat");
	ReadDiffFromBin<MyType>(fourier_d1, tnx, "fourier_d1.dat");
	ReadDiffFromBin<MyType>(fourier_d2, tnx, "fourier_d2.dat");
	
	grid<MyType> rgrid(Rdimensions, trgrid, ROrders,rlength);
	grid<MyType> xgrid(Xdimensions, txgrid, XOrders,1);
	
	delete [] trgrid;
    delete [] txgrid;
	
	derivative<MyType> d10(cheby_d1,idx,&rgrid,&xgrid);
	derivative<MyType> d20(cheby_d2,idx,&rgrid,&xgrid);
	derivative<MyType> d01(idr,fourier_d1,&rgrid,&xgrid);
	derivative<MyType> d02(idr,fourier_d2,&rgrid,&xgrid);
	derivative<MyType> d11(cheby_d1,fourier_d1,&rgrid,&xgrid);
    derivative<MyType> d12(cheby_d1,fourier_d2,&rgrid,&xgrid);
	
	//Start playing
	
	MyType ps[8];
	
	ReadDoubleParameters(ps, 7, "ps.dat");
	
	MyComplex **test;
	MyComplex **test2;
	
	test= new MyComplex*[tnp];
	test2= new MyComplex*[tnp];
	
    for (My_Int i=0; i<tnp; i++) {
        test[i]=new MyComplex[tnx];
		test2[i]=new MyComplex[tnx];
    }
	
    for (My_Int i=0; i<tnp; i++) {
        for (My_Int j=0; j<tnx; j++) {
			test[i][j]=0.;
			test2[i][j]=-1.;
        }}
	
	Cfunction<MyType, MyComplex> Q[14];
	
    function<MyType, MyComplex> testf(test,&rgrid,&xgrid);
	
	derivativeCol<MyType> D(6);
	
	D[0]=d10; D[1]=d20; D[2]=d11; D[3]=d01; D[4]=d02, D[5]=d12;
	
	
	My_Int index[]={7,8,9,10,11,12,13};
	
	My_Int index2[]={0,1,3,4,5,6,2};
	
	for (My_Int ii=0; ii<14; ++ii) {
		Q[ii]=Cfunction<MyType, MyComplex>(testf,D);}
	
	long NN(GetLength<MyType>(1, "freqmesh.dat"));
	long MM(GetLength<MyType>(1, "tempstamp"));
	
	std::cout<< NN << std::endl;
	
	
	for(My_Int ii=0; ii<NN;++ii){
		ReadStamp(ps[7],ii, "freqmesh.dat");
		std::cout << ps[7] << "\t";
	}
	
	std::cout << std::endl << MM << std::endl;
	
	for(My_Int ii=0; ii<MM;++ii){
		ReadStamp(ps[0],ii, "tempstamp");
		std::cout << ps[0] << "\t";
	}
	
	std::cout << std::endl;
	
	My_Int count(0);
	
	for(My_Int jj=0; jj<MM; ++jj){
		
	ReadFromBin(Q, index, 7, tnp, tnx, "thermodata",jj);
	ReadStamp(ps[0], jj, "tempstamp");
		
	for (My_Int ii=0; ii<14; ++ii) {
			Q[ii].update(D);}
	
	count=0;
	
	{My_Int nn(0);MyType w1=0, w2=0,ww=0;
		while(nn<NN){
			ReadStamp(ps[7],nn, "freqmesh.dat");
			w2=MaxEOM<MyType, MyComplex>(equations, 7, Q, rgrid,xgrid, ps);
			w1=MaxBCond<MyType, MyComplex>(bconditions,14, Q, ps,  xgrid, rgrid.TotalLength()-1);
			w2=w1>w2?w1:w2;
			std::cout << w2 << std::endl;
			ww=0;
			ww=w2;
			start1 = clock();
            omp_start = omp_get_wtime();
			UpdateFunctions(equations, bconditions, BulkLinearOp, BoundaryLinearOp, Q, D, rgrid,xgrid, 7, 7, 14, ps );
            omp_end = omp_get_wtime();
			diff = ( std::clock() - start1) / (MyType)CLOCKS_PER_SEC;
			std::cout<<"CPU time: "<< diff <<'\n';
            std::cout<<"Wall time: "<< omp_end-omp_start <<'\n';
			w2=MaxEOM(equations, 7, Q, rgrid,xgrid, ps);
			std::cout << w2 << std::endl;
			w1=MaxBCond<MyType, MyComplex>(bconditions,14, Q, ps,  xgrid, rgrid.TotalLength()-1);
			std::cout << w1 << std::endl;
			w2=w1>w2?w1:w2;
			std::cout << "MaxError=" << w2 << ", mu=" << ps[0] << ", omega=" << ps[7] << std::endl;
			
			if(w2< 1.E-6){
				AppendToBin(Q, index2, 7, tnp, tnx, "DCCdata");
				WriteStamp(ps[7], "DCfreqstamp");
				++count;
			}
			nn++;
		}
	}
		WriteStamp<My_Int>(count, "countindex");
	}
	
    return 0;
}


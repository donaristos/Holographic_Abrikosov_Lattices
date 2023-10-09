//
//  main.cpp
//  TestData
//
//  Created by Aristomenis Donos on 22/11/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#define _USE_MATH_DEFINES
//#define EIGEN_NO_DEBUG

#include <mpi.h>
#include <iostream>
#include "omp.h"
#include "../ReadWrite.h"
#include "../function.h"
#include "../CustomFunctions.h"
#include "../MyTypes.h"
#include "../NewtonMethod.h"
#include "../ReadWrite.h"

#include "../LatticeEOMS.h"
#include "../LatticeLinearOp.h"

#ifdef F128
typedef float128 MyType;
#elif MP
typedef mpreal MyType;
#elif LONGD
typedef long double MyType;
#else
typedef double MyType;
#endif

#define MySize sizeof(MyType)

typedef std::complex<MyType> MyComplex;

MyType **cheby_d1;
MyType **cheby_d2;
MyType **fourierX_d1;
MyType **fourierX_d2;
MyType **fourierY_d1;
MyType **fourierY_d2;
MyType **idr=0;
MyType **idx=0;
MyType **idy=0;

clock_t start1,start2;
double omp_start, omp_end;
MyType diff;


int main(int argc, char * argv[])
{
    
    int numtasks, rank;
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
//    paralution::init_paralution();
    mpfr::mpreal::set_default_prec(76);
    
    My_Int* Rdimensions;
    My_Int* ROrders;
    My_Int* Xdimensions;
    My_Int* Ydimensions;
    
    My_Int rlength=GetLength<My_Int>(1, "Rdimensions.dat");
    
    Rdimensions= new My_Int[rlength];
    ROrders= new My_Int[rlength];
    Xdimensions= new My_Int[1];
    Ydimensions= new My_Int[1];
    
    ReadIntParameters(Rdimensions, rlength, "Rdimensions.dat");
    ReadIntParameters(ROrders, rlength, "ROrders.dat");
    ReadIntParameters(Xdimensions, 1, "Xdimensions.dat");
    ReadIntParameters(Ydimensions, 1, "Ydimensions.dat");
    
    My_Int tnp(0);
    My_Int tnx=Xdimensions[0];
    My_Int tny=Ydimensions[0];
    
    for(My_Int i=0;i<rlength;++i){
        tnp+=Rdimensions[i];
    }
    
    My_Int XOrders[]={Xdimensions[0]};
    My_Int YOrders[]={Ydimensions[0]};
    
    MyType* trgrid;
    MyType* txgrid;
    MyType* tygrid;
    
    trgrid= new MyType[tnp];
    txgrid= new MyType[tnx];
    tygrid= new MyType[tny];
    
    //	setxgrid(txgrid, tnx);
    
    ReadArray(trgrid, tnp, "rgrid.dat");
    ReadArray(txgrid, tnx, "xgrid.dat");
    ReadArray(tygrid, tny, "ygrid.dat");
    
    //Define the Chebyshev differentiation matrix and its square
    
    cheby_d1= new MyType*[tnp];
    cheby_d2= new MyType*[tnp];
    for (My_Int i1=0; i1<tnp;i1++){cheby_d1[i1]=new MyType[tnp];cheby_d2[i1]=new MyType[tnp];}
    
    //Define the Fourier differenations matrix and its square
    
    fourierX_d1= new MyType*[tnx];
    fourierX_d2= new MyType*[tnx];
    for (My_Int i1=0; i1<tnx ;i1++){fourierX_d1[i1]=new MyType[tnx];fourierX_d2[i1]=new MyType[tnx];}
    
    fourierY_d1= new MyType*[tny];
    fourierY_d2= new MyType*[tny];
    for (My_Int i1=0; i1<tny ;i1++){fourierY_d1[i1]=new MyType[tny];fourierY_d2[i1]=new MyType[tny];}
    
    ReadDiffFromBin<MyType>(cheby_d1,  tnp, "cheby_d1.dat");
    ReadDiffFromBin<MyType>(cheby_d2,  tnp, "cheby_d2.dat");
    ReadDiffFromBin<MyType>(fourierX_d1, tnx, "fourierX_d1.dat");
    ReadDiffFromBin<MyType>(fourierX_d2, tnx, "fourierX_d2.dat");
    ReadDiffFromBin<MyType>(fourierY_d1, tny, "fourierY_d1.dat");
    ReadDiffFromBin<MyType>(fourierY_d2, tny, "fourierY_d2.dat");
    
    
    grid<MyType> rgrid(Rdimensions, trgrid, ROrders,rlength);
    grid<MyType> xgrid(Xdimensions, txgrid, XOrders,1);
    grid<MyType> ygrid(Ydimensions, tygrid, YOrders,1);
    
    delete [] trgrid;
    delete [] txgrid;
    delete [] tygrid;
    
    derivative<MyType> d100(cheby_d1,idx,idy,&rgrid,&xgrid,&ygrid);
    derivative<MyType> d200(cheby_d2,idx,idy,&rgrid,&xgrid,&ygrid);
    derivative<MyType> d010(idr,fourierX_d1,idy,&rgrid,&xgrid,&ygrid);
    derivative<MyType> d020(idr,fourierX_d2,idy,&rgrid,&xgrid,&ygrid);
    derivative<MyType> d110(cheby_d1,fourierX_d1,idy,&rgrid,&xgrid,&ygrid);
    derivative<MyType> d001(idr,idx,fourierY_d1,&rgrid,&xgrid,&ygrid);
    derivative<MyType> d002(idr,idx,fourierY_d2,&rgrid,&xgrid,&ygrid);
    derivative<MyType> d101(cheby_d1,idx,fourierY_d1,&rgrid,&xgrid,&ygrid);
    derivative<MyType> d011(idr,fourierX_d1,fourierY_d1,&rgrid,&xgrid,&ygrid);
    
    derivativeCol<MyType> D(9);
    
    D[0]=d100; D[1]=d200; D[2]=d110; D[3]=d010; D[4]=d020; D[5]=d001; D[6]=d002; D[7]=d101; D[8]=d011;
    
    //Start playing
    
    My_Int index[15]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    
    MyType *ps;
    
    MyType PI=acos(MyType(-1));
    
    MyType ***test0;
    
    test0=new MyType**[tnp];
    
    
    for (My_Int i=0; i<tnp; ++i){
        test0[i]=new MyType*[tnx];
        for (My_Int j=0; j<tnx; ++j){
            test0[i][j]=new MyType[tny];
        }
    }
    
    
    
#pragma omp parallel for
    for (My_Int i=0; i<tnp; ++i){
        for (My_Int j=0; j<tnx; ++j){
            for (My_Int k=0; k<tny; ++k){
                test0[i][j][k]=0;
            }
        }
    }
    
    
    function<MyType, MyType> testf0(test0,&rgrid,&xgrid,&ygrid);
    
    Cfunction<MyType, MyType> testC0(testf0,D);
    
    Cfunction<MyType, MyType> Q[15];
    
    Q[0]=testC0;
    Q[1]=testC0;
    Q[2]=testC0;
    Q[3]=testC0;
    Q[4]=testC0;
    Q[5]=testC0;
    Q[6]=testC0;
    Q[7]=testC0;
    Q[8]=testC0;
    Q[9]=testC0;
    Q[10]=testC0;
    Q[11]=testC0;
    Q[12]=testC0;
    Q[13]=testC0;
    Q[14]=testC0;
    
    ReadFromBin(Q, index, 15, tnp, tnx, tny, "tempdata");
    
    for (My_Int i=0; i<15; ++i) {
        Q[i].update(D);
    }
	
   	int NN(GetLength<MyType>(9, "psmesh.dat"));
    
    int Njobs=NN/(numtasks-1);
    
    int Ntodo=Njobs*(numtasks-1);
    
    MyType *psArray;
    
    if (rank==0) {
    
   std::cout << "Np= " << tnp << " , Nx= " << tnx << " , Ny= " << tny << std::endl;
        
   psArray= new MyType[NN*9];
    MyType *psArrayCopy;
    psArrayCopy=new MyType[NN*9];
        
    MyType **ArrayOfSols;
    ArrayOfSols=new MyType*[numtasks-1];
        for (int ii=0; ii<numtasks-1; ii++) {
         ArrayOfSols[ii]=new MyType[15*tnp*tnx*tny];
        }
    
    MPI_Request* Reqs;
    Reqs=new MPI_Request[Ntodo];
    int MPI_Index;
        
    int *DoneIndex;
    DoneIndex=new int[NN];
        for (int ii=0; ii<NN; ii++) {
            DoneIndex[ii]=0;
        }
    
    ReadArray(psArray, 9*NN, "psmesh.dat");
    ReadArray(psArrayCopy, 9*NN, "psmesh.dat");
	std::cout<< Ntodo << std::endl;
	
    for (My_Int ii=0; ii<Ntodo; ii++) {
        std::cout << "mu= " << psArray[ii*9]<< ", rh= "<< psArray[ii*9+1] << ", H= "<< psArray[ii*9+2]
         << ", R0= "<< psArray[ii*9+3]<< ", L1= "<< psArray[ii*9+4]<< ", L2= "<< psArray[ii*9+5]
        << std::endl;}
        
        for (int ii=0; ii<numtasks-1 ; ii++) {
            MPI_Send ((psArray +ii*Njobs*9),9*Njobs*MySize,MPI_BYTE,1+ii,0 ,MPI_COMM_WORLD);
        }
        
        for (int ii=0; ii<Ntodo; ii++) {
            MPI_Irecv( ArrayOfSols[ii/Njobs], 15*tnp*tnx*tny*MySize, MPI_BYTE, MPI_ANY_SOURCE,
                      ii, MPI_COMM_WORLD, Reqs+ii );}
        
        int count(0);
        int DummyIndex(0);
        do {count++;
            MPI_Waitany(Ntodo, Reqs, &MPI_Index, &status);
            AppendArray(ArrayOfSols[status.MPI_TAG/Njobs], 15*tnp*tnx*tny, "thermodata");
            AppendArray(psArray+9*status.MPI_TAG,9,"psstamp.dat");
            DoneIndex[status.MPI_TAG]=1;
            for (int ii=0; ii<NN ; ii++) {
                if(!DoneIndex[ii]){
                    for(int jj=0; jj<9;jj++){
                        psArrayCopy[DummyIndex*9+jj]=psArray[9*ii+jj];
                    }
                    DummyIndex++;
                }
            }
            WriteArray(psArrayCopy,9*(NN-count),"psmesh.dat");
        } while (count<Ntodo);
        
//        WriteArray(psArray+(nn+1)*9,9*(NN-nn-1),"psmesh.dat");
        
        delete[] psArray;
        for (int ii=0; ii<numtasks-1; ii++) {
            delete[] ArrayOfSols[ii];
        }
        delete[] ArrayOfSols;
        delete[] Reqs;
    }
    else if(rank<Ntodo+1){
       paralution::init_paralution();
	std::cout << std::endl;
      psArray=new MyType[9*Njobs];
        MyType *ArrayOfSol;
        ArrayOfSol=new MyType[15*tnp*tnx*tny];
        
      MPI_Recv(psArray,9*Njobs*MySize,MPI_BYTE,0,0,MPI_COMM_WORLD,&status);
        
	{My_Int nn(0); MyType w1=0, w2=0, w3=0, w4=0, tw3=0,tw4=0;
		while (nn<Njobs) {
        ps=psArray+nn*9;
		std::cout << "==========" << std::endl;
            std::cout << "mu= " << ps[0]<< ", rh= "<< ps[1] << ", H= "<< ps[2]
            << ", R0= "<< ps[3]<< ", L1= "<< ps[4]<< ", L2= "<< ps[5]
            << std::endl;
        std::cout << "==========" << std::endl;
        w2=MaxEOM(equations, 15, Q, rgrid, xgrid, ygrid, ps);
        w1=MaxBCond(bconditions, 15, Q, ps, xgrid, ygrid, rgrid.TotalLength()-1);
		w2=w1>w2?w1:w2;
		std::cout << w2 << std::endl;
        w3=MaxEOM(ksi2, 1, Q, rgrid, xgrid, ygrid, ps);
        w4=MaxEOM(psi, 1, Q, rgrid, xgrid, ygrid, ps);
		do
		{	tw3=w3;
            tw4=w4;
			start1 = clock();
            omp_start = omp_get_wtime();
            UpdateFunctions<MyType,MyType>(equations, bconditions, BulkLinearOp, BoundaryLinearOp, Q, D, rgrid, xgrid, ygrid,15, ps );
            omp_end = omp_get_wtime();
			diff = ( std::clock() - start1) / (MyType)CLOCKS_PER_SEC;
			std::cout<<"CPU time: "<< diff <<'\n';
            std::cout<<"Wall time: "<< omp_end-omp_start <<'\n';
            w2=MaxEOM(equations, 15, Q, rgrid, xgrid, ygrid, ps);
            w1=MaxBCond(bconditions, 15, Q, ps, xgrid, ygrid, rgrid.TotalLength()-1);
            w2=w1>w2?w1:w2;
            w3=MaxEOM(ksi2, 1, Q, rgrid, xgrid, ygrid, ps);
            w4=MaxEOM(psi, 1, Q, rgrid, xgrid, ygrid, ps);
            std::cout << "MaxError= " << w2 << std::endl;
            std::cout << "ksi= " << w3 << std::endl;
            std::cout << "psi= " << w4 << std::endl;
            WriteToBin(Q, index, 15, tnp, tnx, tny, "tempdata");}
		while((w2>1.E-7) || fabs(1-(tw3/w3))>1. || fabs(1-(tw4/w4))>1. );
	
            for (My_Int i1=0; i1<15; i1++) {
#pragma omp parallel for
                for (My_Int i2=0; i2<tnp; i2++) {
                    for (My_Int i3=0; i3<tnx; i3++) {
                        for (My_Int i4=0; i4<tny; i4++) {
                        
                            ArrayOfSol[i1*tnp*tnx*tny+i2*tnx*tny+i3*tny+i4]=Q[i1](i2,i3,i4);
                            
                        }}}}
            
            MPI_Send(ArrayOfSol,15*tnp*tnx*tny*MySize,MPI_BYTE,0,nn+Njobs*(rank-1),MPI_COMM_WORLD);
            
			++nn;
	}
	}
        delete[] ArrayOfSol;
        paralution::stop_paralution();
    }
    
    MPI_Finalize();
//paralution::stop_paralution();
    return 0;
}


//
//  main.cpp
//  TestData
//
//  Created by Aristomenis Donos on 22/11/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#define _USE_MATH_DEFINES
//#define EIGEN_NO_DEBUG

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

typedef std::complex<MyType> MyComplex;

MyType *cheby_d1;
MyType *cheby_d2;
MyType *fourierX_d1;
MyType *fourierX_d2;
MyType *fourierY_d1;
MyType *fourierY_d2;
MyType *idr=0;
MyType *idx=0;
MyType *idy=0;

clock_t start1,start2;
double omp_start, omp_end, ltime, mtime;
MyType ldiff, mdiff;


int main(int argc, char * argv[])
{
    mpfr::mpreal::set_default_prec(76);
    
    int numtasks;
    int rank;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    PetscInitialize(&argc,&argv,(char*)0,NULL);
    
    My_Int* Rdimensions;
    My_Int* ROrders;
    My_Int* YOrders;
    My_Int* Xdimensions;
    My_Int* Ydimensions;
    
    My_Int rlength=GetLength<My_Int>(1, "Rdimensions.dat");
    
    Rdimensions= new My_Int[rlength];
    ROrders= new My_Int[rlength];
    YOrders= new My_Int[1];
    Xdimensions= new My_Int[1];
    Ydimensions= new My_Int[1];
    
    ReadIntParameters(Rdimensions, rlength, "Rdimensions.dat");
    ReadIntParameters(ROrders, rlength, "ROrders.dat");
    ReadIntParameters(Xdimensions, 1, "Xdimensions.dat");
    ReadIntParameters(YOrders, 1, "YOrders.dat");
    ReadIntParameters(Ydimensions, 1, "Ydimensions.dat");
    
    My_Int tnp(0);
    My_Int tnx=Xdimensions[0];
    My_Int tny=Ydimensions[0];
    
    for(My_Int i=0;i<rlength;++i){
        tnp+=Rdimensions[i];
    }
    
    if(rank==0){
        std::cout << "Np= " << tnp << " , Nx= " << tnx << " , Ny= " << tny << std::endl;}
    
    My_Int XOrders[]={Xdimensions[0]};
    
    MyType* trgrid;
    MyType* txgrid;
    MyType* tygrid;
    
    trgrid= new MyType[tnp];
    txgrid= new MyType[tnx];
    tygrid= new MyType[tny];
    
    ReadArray(trgrid, tnp, "rgrid.dat");
    ReadArray(txgrid, tnx, "xgrid.dat");
    ReadArray(tygrid, tny, "ygrid.dat");
    
    //Define the Chebyshev differentiation matrix and its square
    
    cheby_d1= new MyType[tnp*tnp];
    cheby_d2= new MyType[tnp*tnp];
    
    //Define the Fourier differenations matrix and its square
    
    fourierX_d1= new MyType[tnx*tnx];
    fourierX_d2= new MyType[tnx*tnx];
    
    fourierY_d1= new MyType[tny*tny];
    fourierY_d2= new MyType[tny*tny];
    
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
    
    My_Int index[16]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    
    MyType ps[9];
    
    ReadArray(ps, 9, "ps.dat");
    
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
   
    Cfunction<MyType, MyType> Q[16];
    
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
    Q[15]=testC0;
    
    
    ReadFromBin(Q, index, 16, tnp, tnx, tny, "data");
    
    for (My_Int i=0; i<16; ++i) {
        Q[i].update(D);
    }

    
    {MyType w1=0, w2=0;
        w2=MaxEOM(equations, 16, Q, rgrid, xgrid, ygrid, ps);
        w1=MaxBCond(bconditions, 16, Q, ps, rgrid, xgrid, ygrid, rgrid.TotalLength()-1);
        w2=w1>w2?w1:w2;
        
        if(rank==0){
        std::cout <<"MaxError= " << w2 << std::endl;}}
    
    {MyType w1=0, w2=0, w3=0, w4=0, tw2=0, tw3=0,tw4=0;
        w2=MaxEOM(equations, 16, Q, rgrid, xgrid, ygrid, ps);
        w1=MaxBCond(bconditions, 16, Q, ps, rgrid, xgrid, ygrid, rgrid.TotalLength()-1);
        w3=MaxEOM(ksi2, 1, Q, rgrid, xgrid, ygrid, ps);
        w4=MaxEOM(psi, 1, Q, rgrid, xgrid, ygrid, ps);
        do
        {   tw2=w2;
            tw3=w3;
            tw4=w4;
            start1 = clock();
            omp_start = omp_get_wtime();
            UpdateFunctions<MyType,MyType>(equations, bconditions, BulkLinearOp, BoundaryLinearOp, Q, D, rgrid, xgrid, ygrid,16, ps );
            omp_end = omp_get_wtime();
            ldiff = ( std::clock() - start1) / (double)CLOCKS_PER_SEC;
            ltime = omp_end - omp_start;
            MPI_Reduce(&ldiff,&mdiff,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);
            MPI_Reduce(&ltime,&mtime,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);
            if(rank==0){
            std::cout<<"CPU time: "<< ldiff <<'\n';
                std::cout<<"Wall time: "<< mtime <<'\n';}
            
            w2=MaxEOM(equations, 16, Q, rgrid, xgrid, ygrid, ps);
	        w1=MaxBCond(bconditions, 16, Q, ps,rgrid, xgrid, ygrid, rgrid.TotalLength()-1);
            w2=w1>w2?w1:w2;
            w3=MaxEOM(ksi2, 1, Q, rgrid, xgrid, ygrid, ps);
            w4=MaxEOM(psi, 1, Q, rgrid, xgrid, ygrid, ps);
            if(rank==0){
            std::cout << "MaxError= " << w2 << std::endl;
            std::cout << "ksi= " << w3 << std::endl;
            std::cout << "psi= " << w4 << std::endl;
            WriteToBin(Q, index, 16, tnp, tnx, tny, "data");}
        }while(fabs(1-(tw2/w2))> .2 || fabs(1-(tw3/w3))> .2 || fabs(1-((1+tw4)/(1+w4)))> .2);}
    
    
    PetscFinalize();
    MPI_Finalize();
    return 0;
}
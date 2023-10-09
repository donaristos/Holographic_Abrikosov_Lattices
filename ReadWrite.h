//
//  ReadWrite.h
//  TimeDependence
//
//  Created by Aristomenis Donos on 13/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#ifndef __TimeDependence__ReadWrite__
#define __TimeDependence__ReadWrite__

#include <iostream>
#include <fstream>
#include <vector>
#include "function.h"
#include "MyTypes.h"


void ReadIntParameters(My_Int *parameters, My_Int length,const std::string& file_path);

//General templates

template<typename Type>
int GetLength(My_Int CellLength,const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	long begin(0),end(0);
	if (is.is_open()) {
		is.seekg(0,is.beg);
		begin = is.tellg();
		is.seekg(0,is.end);
		end= is.tellg();
	}
	return (end-begin)/(CellLength*(std::streamsize(sizeof(Type))));
}

template<typename Type>
void ReadArray(Type *parameters, My_Int length, const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	if(is.is_open()){
//		std::cout << "It is open now!"<< std::endl;
		is.seekg(0,is.beg);
		is.read(reinterpret_cast<char *>(parameters),std::streamsize(length*sizeof(Type)));
	}
}

template<typename Type>
void WriteArray(Type *ArPointer, My_Int length, const std::string& file_path, My_Int prec=30){
    std::ofstream os(file_path.c_str(),std::ios::out|std::ios::binary);
    os.write(reinterpret_cast<char *>(ArPointer), std::streamsize(length*sizeof(Type)));
    os.close();
}

template<typename Type>
void AppendArray(Type *ArPointer, My_Int length, const std::string& file_path, My_Int prec=30){
    std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out|std::ios::binary);
    os.write(reinterpret_cast<char *>(ArPointer), std::streamsize(length*sizeof(Type)));
    os.close();
}

template<typename Type>
void ReadDiffFromBin(Type* diff_m,My_Int dim, const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	if(is.is_open()){
//		std::cout << "Diff file is open" << std::endl;
		is.seekg(0,is.beg);
			is.read(reinterpret_cast<char *>(diff_m), std::streamsize(dim*dim*sizeof(Type)));
	}
}



template<typename dType, typename Type>
void WriteToBin(Cfunction<dType, Type> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec=30){
	Type *tarray;
	
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
	tarray=new Type[array_size];
	
	
	for (My_Int i=0; i<Nf; ++i){
        temp_int=i*Np*Nx*Ny;
		for (My_Int j=0; j<Np; ++j) {
			for(My_Int k=0; k<Nx;++k){
                for(My_Int l=0; l<Ny;++l){
                 tarray[temp_int]=F[ind[i]](j,k,l);
                    temp_int+=1;
                }
			}
		}
	}
	
    std::ofstream os(file_path.c_str(),std::ios::out|std::ios::binary);
	

    os.write(reinterpret_cast<char *>(tarray), std::streamsize(array_size*sizeof(Type)));
    
    delete[] tarray;
    
    os.close();
}


template<typename dType, typename Type>
void AppendToBin(Cfunction<dType, Type> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec=30){
   
    Type *tarray;
    
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
    tarray=new Type[array_size];
    
    
    for (My_Int i=0; i<Nf; ++i){
        temp_int=i*Np*Nx*Ny;
        for (My_Int j=0; j<Np; ++j) {
            for(My_Int k=0; k<Nx;++k){
                for(My_Int l=0; l<Ny;++l){
                    tarray[temp_int]=F[ind[i]](j,k,l);
                    temp_int+=1;
                }
            }
        }
    }
	
    std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out|std::ios::binary);
	
   os.write(reinterpret_cast<char *>(tarray), std::streamsize(array_size*sizeof(Type)));
    
    
    delete [] tarray;
    
    os.close();
}

template<typename dType, typename Type>
void ReadFromBin(Cfunction<dType, Type> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int pos=-1){
    
	dType *tarray;
	
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
	tarray=new dType[array_size];
	
    std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	
	if (is.is_open()) {
		if (pos<0) {
			is.seekg(pos*std::streamsize(Nf*Np*Nx*Ny*sizeof(dType)),is.end);
		}
		else
		{
			is.seekg(pos*std::streamsize(Nf*Np*Nx*Ny*sizeof(dType)),is.beg);
		}
		
 
        
        is.read(reinterpret_cast<char *>(tarray), std::streamsize(array_size*sizeof(dType)));
		
        
   
        for (My_Int i=0; i<Nf; ++i){
            temp_int=i*Np*Nx*Ny;
            for (My_Int j=0; j<Np; ++j) {
                for(My_Int k=0; k<Nx;++k){
                    for(My_Int l=0; l<Ny;++l){
					F[ind[i]](j,k,l)=tarray[temp_int];
                    temp_int+=1;
                    }
				}
			}
		}
    }
    
    delete [] tarray;
    
    is.close();
}

//Specialization to mpreal type

template<>
int GetLength<mpreal>(My_Int CellLength,const std::string& file_path);

template<>
void ReadArray<mpreal>(mpreal *parameters, My_Int length, const std::string& file_path);

template<>
void WriteArray<mpreal>(mpreal *ArPointer, My_Int length, const std::string& file_path, My_Int prec);

template<>
void AppendArray<mpreal>(mpreal *ArPointer, My_Int length, const std::string& file_path, My_Int prec);

template<>
void ReadDiffFromBin<mpreal>(mpreal* diff_m,My_Int dim, const std::string& file_path);


template<>
void WriteToBin<mpreal, mpreal>(Cfunction<mpreal, mpreal> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec);

template<>
void AppendToBin<mpreal, mpreal>(Cfunction<mpreal, mpreal> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec);

template<>
void ReadFromBin<mpreal, mpreal>(Cfunction<mpreal, mpreal> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int pos);

template<>
void WriteToBin<mpreal, ComplexMP>(Cfunction<mpreal, ComplexMP> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec);

template<>
void AppendToBin<mpreal, ComplexMP>(Cfunction<mpreal, ComplexMP> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec);

template<>
void ReadFromBin<mpreal, ComplexMP>(Cfunction<mpreal, ComplexMP> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int pos);

//Specialization to long double type

template<>
int GetLength<long double>(My_Int CellLength,const std::string& file_path);

template<>
void ReadArray<long double>(long double *parameters, My_Int length, const std::string& file_path);

template<>
void WriteArray<long double>(long double *ArPointer, My_Int length, const std::string& file_path, My_Int prec);

template<>
void AppendArray<long double>(long double *ArPointer, My_Int length, const std::string& file_path, My_Int prec);

template<>
void ReadDiffFromBin<long double>(long double* diff_m,My_Int dim, const std::string& file_path);

template<>
void WriteToBin<long double, long double>(Cfunction<long double, long double> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec);

template<>
void AppendToBin<long double, long double>(Cfunction<long double, long double> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec);

template<>
void ReadFromBin<long double, long double>(Cfunction<long double, long double> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int pos);

template<>
void WriteToBin<long double, ComplexLD>(Cfunction<long double, ComplexLD> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec);

template<>
void AppendToBin<long double, ComplexLD>(Cfunction<long double, ComplexLD> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec);

template<>
void ReadFromBin<long double, ComplexLD>(Cfunction<long double, ComplexLD> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int pos);

#endif /* defined(__TimeDependence__ReadWrite__) */

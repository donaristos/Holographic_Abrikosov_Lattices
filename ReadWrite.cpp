//
//  ReadWrite.cpp
//  TimeDependence
//
//  Created by Aristomenis Donos on 13/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#include "ReadWrite.h"

void ReadIntParameters(My_Int *parameters, My_Int length, const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	if(is.is_open()){
//		std::cout << "It is open now!"<< std::endl;
		is.seekg(0,is.beg);
		is.read(reinterpret_cast<char *>(parameters),std::streamsize(length*sizeof(My_Int)));
	}
}

//Implementation of full mpreal specialization

template<>
int GetLength<mpreal>(My_Int CellLength,const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	My_Int i(0);
	if(is.is_open()){
		mpreal temp;
		while(is >> temp){
			++i;}
	}
	return i;
}

template<>
void ReadArray<mpreal>(mpreal *parameters, My_Int length, const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in);
	if(is.is_open()){
//		std::cout << "It is open now!"<< std::endl;
		is.seekg(0,is.beg);
		for (My_Int ii=0; ii<length; ++ii) {
			is >> parameters[ii];
		}
	}
}

template<>
void WriteArray<mpreal>(mpreal *ArPointer, My_Int length, const std::string& file_path, My_Int prec){
    std::ofstream os(file_path.c_str(),std::ios::out);
    os << std::fixed;
    os.precision(prec);
    for (My_Int ii=0; ii<length; ii++) {
	os <<  *(ArPointer+ii) << "\t" ;
    }
    os.close();
}

template<>
void AppendArray<mpreal>(mpreal *ArPointer, My_Int length, const std::string& file_path, My_Int prec){
    std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out);
    os << std::fixed;
    os.precision(prec);
    for (My_Int ii=0; ii<length; ii++) {
        os <<  *(ArPointer+ii) << "\t" ;
    }
    os.close();
}


template<>
void ReadDiffFromBin<mpreal>(mpreal* diff_m,My_Int dim, const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in);
	if(is.is_open()){
//		std::cout << "Diff file is open" << std::endl;
		is.seekg(0,is.beg);
		for (My_Int ii=0; ii<dim*dim; ++ii) {
				is >> diff_m[ii];
			}
	}
}


template<>
void WriteToBin<mpreal, mpreal>(Cfunction<mpreal, mpreal> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec){
    
    mpreal *tarray;
    
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
    tarray=new mpreal[array_size];
    
    
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
	
    std::ofstream os(file_path.c_str(),std::ios::out);
	
	os << std::fixed;
	os.precision(prec);
	
	for (My_Int i=0; i<array_size; ++i) {
				os << tarray[i] << "\t";
                }
    
    delete [] tarray;
    
    os.close();
}

template<>
void AppendToBin<mpreal, mpreal>(Cfunction<mpreal, mpreal> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec){
    
    mpreal *tarray;
    
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
    tarray=new mpreal[array_size];
    
    
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
	
    std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out);
	
	os << std::fixed;
	os.precision(prec);
	
	
    for (My_Int i=0; i<array_size; ++i) {
        os << tarray[i] << "\t";
    }
    
    delete [] tarray;
    
    os.close();
}

template<>
void ReadFromBin<mpreal, mpreal>(Cfunction<mpreal, mpreal> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int pos){

    
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
	
    std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	
	if(is.is_open()){
		
		is.seekg(0, is.beg);
		
		mpreal tempvar;
		std::vector<mpreal> tempvec;
		My_Int length;
		My_Int position;
		
		while(is >> tempvar){
			tempvec.push_back(tempvar);
		}
		
		length=tempvec.size();
		
		if (pos<0) {
			position=length+ pos*Ny*Nx*Np*Nf;
		}
		else
		{
			position=pos*Ny*Nx*Np*Nf;
		}
		
        for (My_Int i=0; i<Nf; ++i){
            temp_int=position+i*Np*Nx*Ny;
            for (My_Int j=0; j<Np; ++j) {
                for(My_Int k=0; k<Nx;++k){
                    for(My_Int l=0; l<Ny;++l){
                        F[ind[i]](j,k,l)=tempvec[temp_int];
                        temp_int+=1;
                    }
                }
            }
        }
    }
		
  
    is.close();
}

template<>
void WriteToBin<mpreal, ComplexMP>(Cfunction<mpreal, ComplexMP> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec){

   
    ComplexMP *tarray;
    
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
    tarray=new ComplexMP[array_size];
    
    
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
    
    std::ofstream os(file_path.c_str(),std::ios::out);
    
    os << std::fixed;
    os.precision(prec);
    
    for (My_Int i=0; i<array_size; ++i) {
        os << tarray[i] << "\t";
    }
    
    delete [] tarray;
    
    os.close();
}

template<>
void AppendToBin<mpreal, ComplexMP>(Cfunction<mpreal, ComplexMP> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec){

    ComplexMP *tarray;
    
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
    tarray=new ComplexMP[array_size];
    
    
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
    
    std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out);
    
    os << std::fixed;
    os.precision(prec);
    
    
    for (My_Int i=0; i<array_size; ++i) {
        os << tarray[i] << "\t";
    }
    
    delete [] tarray;
    
    os.close();
}

template<>
void ReadFromBin<mpreal, ComplexMP>(Cfunction<mpreal, ComplexMP> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int pos){

    
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
    
    std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
    
    if(is.is_open()){
        
        is.seekg(0, is.beg);
        
        ComplexMP tempvar;
        std::vector<ComplexMP> tempvec;
        My_Int length;
        My_Int position;
        
        while(is >> tempvar){
            tempvec.push_back(tempvar);
        }
        
        length=tempvec.size();
        
        if (pos<0) {
            position=length+ pos*Ny*Nx*Np*Nf;
        }
        else
        {
            position=pos*Ny*Nx*Np*Nf;
        }
        
        for (My_Int i=0; i<Nf; ++i){
            temp_int=position+i*Np*Nx*Ny;
            for (My_Int j=0; j<Np; ++j) {
                for(My_Int k=0; k<Nx;++k){
                    for(My_Int l=0; l<Ny;++l){
                        F[ind[i]](j,k,l)=tempvec[temp_int];
                        temp_int+=1;
                    }
                }
            }
        }
    }
    
    
    
    is.close();
}

//Implementation of full long double specialization


template<>
int GetLength<long double>(My_Int CellLength,const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	My_Int i(0);
	if(is.is_open()){
		long double temp;
		while(is >> temp){
			++i;}
	}
	return i;
}

template<>
void ReadArray<long double>(long double *parameters, My_Int length, const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in);
	if(is.is_open()){
//		std::cout << "It is open now!"<< std::endl;
		is.seekg(0,is.beg);
		for (My_Int ii=0; ii<length; ++ii) {
			is >> parameters[ii];
		}
	}
}

template<>
void WriteArray<long double>(long double *ArPointer, My_Int length, const std::string& file_path, My_Int prec){
    std::ofstream os(file_path.c_str(),std::ios::out);
    os << std::fixed;
    os.precision(prec);
    for (My_Int ii=0; ii<length; ii++) {
        os <<  *(ArPointer+ii) << "\t" ;
    }
    os.close();
}

template<>
void AppendArray<long double>(long double *ArPointer, My_Int length, const std::string& file_path, My_Int prec){
    std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out);
    os << std::fixed;
    os.precision(prec);
    for (My_Int ii=0; ii<length; ii++) {
        os <<  *(ArPointer+ii) << "\t" ;
    }
    os.close();
}

template<>
void ReadDiffFromBin<long double>(long double* diff_m,My_Int dim, const std::string& file_path){
    std::ifstream is(file_path.c_str(),std::ios::in);
    if(is.is_open()){
        //		std::cout << "Diff file is open" << std::endl;
        is.seekg(0,is.beg);
        for (My_Int ii=0; ii<dim*dim; ++ii) {
            is >> diff_m[ii];
        }
    }
}


template<>
void WriteToBin<long double, long double>(Cfunction<long double, long double> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec){
    
    long double *tarray;
    
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
    tarray=new long double[array_size];
    
    
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
    
    std::ofstream os(file_path.c_str(),std::ios::out);
    
    os << std::fixed;
    os.precision(prec);
    
    for (My_Int i=0; i<array_size; ++i) {
        os << tarray[i] << "\t";
    }
    
    delete [] tarray;
    
    os.close();
}

template<>
void AppendToBin<long double, long double>(Cfunction<long double, long double> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec){
    
    long double *tarray;
    
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
    tarray=new long double[array_size];
    
    
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
    
    std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out);
    
    os << std::fixed;
    os.precision(prec);
    
    
    for (My_Int i=0; i<array_size; ++i) {
        os << tarray[i] << "\t";
    }
    
    delete [] tarray;
    
    os.close();
}

template<>
void ReadFromBin<long double, long double>(Cfunction<long double, long double> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int pos){
    
    
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
    
    std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
    
    if(is.is_open()){
        
        is.seekg(0, is.beg);
        
        long double tempvar;
        std::vector<long double> tempvec;
        My_Int length;
        My_Int position;
        
        while(is >> tempvar){
            tempvec.push_back(tempvar);
        }
        
        length=tempvec.size();
        
        if (pos<0) {
            position=length+ pos*Ny*Nx*Np*Nf;
        }
        else
        {
            position=pos*Ny*Nx*Np*Nf;
        }
        
        for (My_Int i=0; i<Nf; ++i){
            temp_int=position+i*Np*Nx*Ny;
            for (My_Int j=0; j<Np; ++j) {
                for(My_Int k=0; k<Nx;++k){
                    for(My_Int l=0; l<Ny;++l){
                        F[ind[i]](j,k,l)=tempvec[temp_int];
                        temp_int+=1;
                    }
                }
            }
        }
    }
    
    
    is.close();
}

template<>
void WriteToBin<long double, ComplexLD>(Cfunction<long double, ComplexLD> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec){
    
    
    ComplexLD *tarray;
    
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
    tarray=new ComplexLD[array_size];
    
    
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
    
    std::ofstream os(file_path.c_str(),std::ios::out);
    
    os << std::fixed;
    os.precision(prec);
    
    for (My_Int i=0; i<array_size; ++i) {
        os << tarray[i] << "\t";
    }
    
    delete [] tarray;
    
    os.close();
}

template<>
void AppendToBin<long double, ComplexLD>(Cfunction<long double, ComplexLD> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int prec){
    
    ComplexLD *tarray;
    
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
    tarray=new ComplexLD[array_size];
    
    
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
    
    std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out);
    
    os << std::fixed;
    os.precision(prec);
    
    
    for (My_Int i=0; i<array_size; ++i) {
        os << tarray[i] << "\t";
    }
    
    delete [] tarray;
    
    os.close();
}

template<>
void ReadFromBin<long double, ComplexLD>(Cfunction<long double, ComplexLD> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx, My_Int Ny,const std::string& file_path, My_Int pos){
    
    
    My_Int array_size=Nf*Np*Nx*Ny;
    My_Int temp_int;
    
    
    std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
    
    if(is.is_open()){
        
        is.seekg(0, is.beg);
        
        ComplexLD tempvar;
        std::vector<ComplexLD> tempvec;
        My_Int length;
        My_Int position;
        
        while(is >> tempvar){
            tempvec.push_back(tempvar);
        }
        
        length=tempvec.size();
        
        if (pos<0) {
            position=length+ pos*Ny*Nx*Np*Nf;
        }
        else
        {
            position=pos*Ny*Nx*Np*Nf;
        }
        
        for (My_Int i=0; i<Nf; ++i){
            temp_int=position+i*Np*Nx*Ny;
            for (My_Int j=0; j<Np; ++j) {
                for(My_Int k=0; k<Nx;++k){
                    for(My_Int l=0; l<Ny;++l){
                        F[ind[i]](j,k,l)=tempvec[temp_int];
                        temp_int+=1;
                    }
                }
            }
        }
    }
    
    
    is.close();
}
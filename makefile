PETSC_DIR=$(HOME)/Downloads/C++/petsc
PETSC_ARCH=(your-PETSC_ARCH)


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

CC=icpc
CC=$(CLINKER)

CPPFLAGS= -ipo -c -DEIGEN_NO_DEBUG -DNDEBUG -DBOOST_MP_USE_QUAD -fPIC \
-xHost -m64 -qopenmp\
-Qoption,cpp,--extended_float_type -stdlib=libstdc++

HPATHS= -I $(HOME)/Downloads/C++/boost -I $(HOME)/Downloads/C++/eigen \
-I $(HOME)/Downloads/C++/include \
-I $(HOME)/Downloads/C++/eigen/unsupported/test/mpreal \
-I $(PETSC_DIR)/$(PETSC_ARCH)/include -I /sw/include -I $(PETSC_DIR)/include \
-I /opt/intel/mkl/include

LPATHS= -L $(HOME)/Downloads/C++/lib 

LFLAGS= -m64 -fPIC -DEIGEN_NO_DEBUG -DNDEBUG -Wall -multiple-processes \
-xHost -align -ansi-alias -restrict -qopenmp -shared-intel -Qoption,cpp,--extended_float_type \
-limf -lmpfr -stdlib=libstdc++ ${PETSC_LIB} -lpthread \
-Xlinker -rpath -Xlinker ../Libs

OBJS=function.o ReadWrite.o NewtonMethod.o CustomFunctions.o

LATTICEFILES=LatticeEOMSD.dylib #LatticeEOMSLD.dylib LatticeEOMSF128.dylib LatticeEOMSMP.dylib

LATTICELINFILES=LatticeLinearOpD.dylib #LatticeLinearOpLD.dylib LatticeLinearOpF128.dylib LatticeLinearOpMP.dylib

all:: Background Thermod

Background: BackgroundDI #BackgroundLDI BackgroundQI BackgroundMPI

Thermod: ThermodDI #ThermodDMPI ThermodLDI ThermodQI ThermodMPI

ThermodDI:: thermoD.o $(LATTICEFILES) $(LATTICELINFILES) ProjectLibs.dylib
	$(CC) -o ThermodDI thermoD.o ProjectLibs.dylib $(LATTICEFILES) $(LATTICELINFILES) $(LFLAGS) $(HPATHS) $(LPATHS)
	
ThermodDMPI:: thermoDMP.o $(LATTICEFILES) $(LATTICELINFILES) ProjectLibs.dylib
	$(CC) -o ThermodDMPI thermoDMP.o ProjectLibs.dylib $(LATTICEFILES) $(LATTICELINFILES) $(LFLAGS) $(HPATHS) $(LPATHS)
	
ThermodLDI:: thermoLD.o $(LATTICEFILES) $(LATTICELINFILES) ProjectLibs.dylib
	$(CC) -o ThermodLDI thermoLD.o ProjectLibs.dylib $(LATTICEFILES) $(LATTICELINFILES) $(LFLAGS) $(HPATHS) $(LPATHS)
	
ThermodQI:: thermoQ.o $(LATTICEFILES) $(LATTICELINFILES) ProjectLibs.dylib
	$(CC) -o ThermodQI thermoQ.o ProjectLibs.dylib $(LATTICEFILES) $(LATTICELINFILES) $(LFLAGS) $(HPATHS) $(LPATHS)
	
ThermodMPI:: thermoMP.o $(LATTICEFILES) $(LATTICELINFILES) ProjectLibs.dylib
	$(CC) -o ThermodMPI thermoMP.o ProjectLibs.dylib $(LATTICEFILES) $(LATTICELINFILES) $(LFLAGS) $(HPATHS) $(LPATHS)

BackgroundDI:: backgroundD.o $(LATTICEFILES) $(LATTICELINFILES) ProjectLibs.dylib
	$(CC) -o BackgroundDI backgroundD.o ProjectLibs.dylib $(LATTICEFILES) $(LATTICELINFILES) $(LFLAGS) $(HPATHS) $(LPATHS)
	
BackgroundLDI:: backgroundLD.o $(LATTICEFILES) $(LATTICELINFILES) ProjectLibs.dylib
	$(CC) -o BackgroundLDI backgroundLD.o ProjectLibs.dylib $(LATTICEFILES) $(LATTICELINFILES) $(LFLAGS) $(HPATHS) $(LPATHS)

BackgroundQI:: backgroundQ.o $(LATTICEFILES) $(LATTICELINFILES) ProjectLibs.dylib
	$(CC) -o BackgroundQI backgroundQ.o ProjectLibs.dylib $(LATTICEFILES) $(LATTICELINFILES) $(LFLAGS) $(HPATHS) $(LPATHS)
	
BackgroundMPI:: backgroundMP.o $(LATTICEFILES) $(LATTICELINFILES) ProjectLibs.dylib
	$(CC) -o BackgroundMPI backgroundMP.o ProjectLibs.dylib $(LATTICEFILES) $(LATTICELINFILES) $(LFLAGS) $(HPATHS) $(LPATHS)

backgroundD.o:: Background/main.cpp
	$(CC) $(CPPFLAGS) $(LPATHS) $(HPATHS) Background/main.cpp -o backgroundD.o
	
backgroundLD.o:: Background/main.cpp
	$(CC) $(CPPFLAGS) -DLONGD $(LPATHS) $(HPATHS) Background/main.cpp -o backgroundLD.o

backgroundQ.o:: Background/main.cpp
	$(CC) $(CPPFLAGS) -DF128 $(LPATHS) $(HPATHS) Background/main.cpp -o backgroundQ.o
	
backgroundMP.o:: Background/main.cpp
	$(CC) $(CPPFLAGS) -DMP $(LPATHS) $(HPATHS) Background/main.cpp -o backgroundMP.o

thermoD.o:: Thermo/main.cpp
	$(CC) $(CPPFLAGS) $(LPATHS) $(HPATHS) Thermo/main.cpp -o thermoD.o
	
thermoDMP.o:: Thermo/MPImain.cpp
	$(CC) $(CPPFLAGS) $(LPATHS) $(HPATHS) Thermo/MPImain.cpp -o thermoDMP.o
	
thermoLD.o:: Thermo/main.cpp
	$(CC) $(CPPFLAGS) -DLONGD $(LPATHS) $(HPATHS) Thermo/main.cpp -o thermoLD.o
	
thermoQ.o:: Thermo/main.cpp
	$(CC) $(CPPFLAGS) -DF128 $(LPATHS) $(HPATHS) Thermo/main.cpp -o thermoQ.o
	
thermoMP.o:: Thermo/main.cpp
	$(CC) $(CPPFLAGS) -DMP $(LPATHS) $(HPATHS) Thermo/main.cpp -o thermoMP.o

LatticeLinearOp%.o:: LatticeLinearOp%.cpp
	$(CC) $(CPPFLAGS) $(LPATHS) $(HPATHS) $< -o $@
	
LatticeLinearOp%.dylib: LatticeLinearOp%.o #$(OBJS)
	$(CC) -dynamiclib -o $@ $< $(OBJS) $(LFLAGS) $(LPATHS) $(HPATHS) -install_name ../Libs/$@
	
LatticeEOMS%.o: LatticeEOMS%.cpp
	$(CC) $(CPPFLAGS) $(LPATHS) $(HPATHS) $< -o $@
	
LatticeEOMS%.dylib: LatticeEOMS%.o #$(OBJS)
	$(CC) -dynamiclib -o $@ $< $(OBJS) $(LFLAGS) $(LPATHS) $(HPATHS) -install_name ../Libs/$@
	
LatticeLinearOpMP.o: LatticeLinearOpMP.cpp
	$(CC) $(CPPFLAGS0) $(LPATHS) $(HPATHS) $< -o $@
	
LatticeLinearOpMP.dylib: LatticeLinearOpMP.o $(OBJS)
	$(CC) -dynamiclib -o $@ $< $(OBJS) $(LFLAGS) $(LPATHS) $(HPATHS) -install_name ../Libs/$@
	
LatticeEOMSMP.o: LatticeEOMSMP.cpp
	$(CC) $(CPPFLAGS0) $(LPATHS) $(HPATHS) $< -o $@
	
LatticeEOMSMP.dylib: LatticeEOMSMP.o $(OBJS)
	$(CC) -dynamiclib -o $@ $< $(OBJS) $(LFLAGS) $(LPATHS) $(HPATHS) -install_name ../Libs/$@

#LEOMS.dylib: $(LATTICEFILES) $(LATTICELINFILES) ProjectLibs.dylib
#	$(CC) -dynamiclib -o LEOMS.dylib $(LATTICEFILES) $(LATTICELINFILES) ProjectLibs.dylib $(LFLAGS2) $(HPATHS) $(LPATHS)
	
ProjectLibs.dylib: $(OBJS)
	$(CC) -dynamiclib -o ProjectLibs.dylib $(OBJS) $(LFLAGS) $(HPATHS) $(LPATHS) -install_name ../Libs/$@

function.o: function.cpp function.h
	$(CC) $(CPPFLAGS) $(LPATHS) $(HPATHS) function.cpp
	
CustomFunctions.o: CustomFunctions.cpp CustomFunctions.h
	$(CC) $(CPPFLAGS) $(LPATHS) $(HPATHS) CustomFunctions.cpp

NewtonMethod.o: NewtonMethod.cpp NewtonMethod.h
	$(CC) $(CPPFLAGS) $(LPATHS) $(HPATHS) NewtonMethod.cpp

ReadWrite.o: ReadWrite.cpp ReadWrite.h
	$(CC) $(CPPFLAGS) $(LPATHS) $(HPATHS) ReadWrite.cpp

clean::
	rm -rf *.o *.dylib BackgroundDI ThermodDI


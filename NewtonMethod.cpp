//
//  NewtonMethod.cpp
//  PDEs
//
//  Created by Aristomenis Donos on 04/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#include "NewtonMethod.h"

PetscErrorCode PCMGSetupViaCoarsen(PC pc,DM da_fine)
{
    PetscInt       nlevels,k,PETSC_UNUSED finest;
    DM             *da_list,*daclist;
    Mat            R;
    PetscErrorCode ierr;
    
    PetscFunctionBeginUser;
    nlevels = 1;
    PetscOptionsGetInt(NULL,NULL,"-levels",&nlevels,0);
    
    ierr = PetscMalloc1(nlevels,&da_list);CHKERRQ(ierr);
    for (k=0; k<nlevels; k++) da_list[k] = NULL;
    ierr = PetscMalloc1(nlevels,&daclist);CHKERRQ(ierr);
    for (k=0; k<nlevels; k++) daclist[k] = NULL;
    
    /* finest grid is nlevels - 1 */
    finest     = nlevels - 1;
    daclist[0] = da_fine;
    PetscObjectReference((PetscObject)da_fine);
    ierr = DMCoarsenHierarchy(da_fine,nlevels-1,&daclist[1]);CHKERRQ(ierr);
    for (k=0; k<nlevels; k++) {
        da_list[k] = daclist[nlevels-1-k];
//        DMDASetUniformCoordinates(da_list[k],0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr);
    }
    
    ierr = PCMGSetLevels(pc,nlevels,NULL);CHKERRQ(ierr);
    ierr = PCMGSetType(pc,PC_MG_MULTIPLICATIVE);CHKERRQ(ierr);
    ierr = PCMGSetGalerkin(pc,PC_MG_GALERKIN_NONE);CHKERRQ(ierr);
    
    for (k=1; k<nlevels; k++) {
        ierr = DMCreateInterpolation(da_list[k-1],da_list[k],&R,NULL);CHKERRQ(ierr);
        ierr = PCMGSetInterpolation(pc,k,R);CHKERRQ(ierr);
        ierr = MatDestroy(&R);CHKERRQ(ierr);
    }
    
    /* tidy up */
    for (k=0; k<nlevels; k++) {
        ierr = DMDestroy(&da_list[k]);CHKERRQ(ierr);
    }
    ierr = PetscFree(da_list);CHKERRQ(ierr);
    ierr = PetscFree(daclist);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

  template<>
  void UpdateFunctions<double, double>( void (*eom) (double* elem, Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i, const My_Int& j, const My_Int& k)
                                       ,void (*b_cond) (double* elem, Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3,double* params, const My_Int& i, const My_Int& j, const My_Int& k)
                                       ,void (*l_eom)(double *elem, derivativeCol<double>& Der,Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k)
                                       ,void (*l_b_cond)(double *elem,  derivativeCol<double>& Der,Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k),Cfunction<double, double>* f, derivativeCol<double>& D, grid<double>& Rgrid,grid<double>& Xgrid, grid<double>& Ygrid,My_Int N_functions, double* params){
    

      
      int numtasks;
      int rank;
      
      double omp_start, omp_end, time, mtime;
      double omp_start_2, omp_end_2, time_2, mtime_2;
      
      MPI_Comm_size(PETSC_COMM_WORLD, &numtasks);
      MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
      
    My_Int Np=Rgrid.TotalLength();
    My_Int Nx=Xgrid.TotalLength();
    My_Int Ny=Ygrid.TotalLength();
    My_Int nbdof = (Np-1)*Nx*Ny*N_functions; // number of degrees of freedom.
    My_Int N_lattice_sites=(Np-1)*Nx*Ny;
 	My_Int nvars = Ny*Nx+ (2*Rgrid.MaxOrder()+1)*(Nx+Ny)+N_functions*(2*Rgrid.MaxOrder()+1+Nx+Ny);
      
      std::vector<Eigen::Triplet<double, My_Int> > TriplVec;
      Eigen::Matrix<double, Eigen::Dynamic, 1> B(nbdof);
      
      TriplVec.reserve((long long)(nvars)* (long long)(nbdof)/(long long)numtasks);
      
      double tole=1.E-9;
      
      Vec            x, b, u, X;          /* approx solution, RHS, exact solution */
      VecScatter    Scatter;
      Mat            A;                /* linear system matrix */
      KSP            ksp;              /* linear solver context */
      PC             pc;               /* preconditioner context */
      PetscReal      norm,tol=1.e-11;  /* norm of solution error */
      PetscErrorCode ierr;
      PetscInt       n(1), ti, tj,rstart,rend,nlocal, its, sten_size;
      PetscInt       const *rstarts;
      PetscScalar    tvalue;
      IS ix,iX;
      
      DM    da;
      PetscInt  *i_map;
      DMDALocalInfo info;
      PetscInt  rs, xs, ys, re, xe, ye, rm, xm, ym;
      AO    ao;
      
      DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_PERIODIC,DM_BOUNDARY_PERIODIC,DMDA_STENCIL_BOX, Np-1,
                   Nx,Ny,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,N_functions,1,NULL,NULL,NULL, &da);
      
      
      DMSetUp(da);
      
      DMDAGetLocalInfo(da,&info);
      
      DMDAGetAO(da,&ao);
      
      rs=info.xs; re=rs+ info.xm;
      xs=info.ys; xe=xs+ info.ym;
      ys=info.zs; ye=ys+ info.zm;

      DMDAGetCorners(da, &rs,&xs,&ys,&rm,&xm,&ym);
	
   
      
      re=rs+rm; xe=xs+xm; ye=ys+ym;
      
      N_lattice_sites= info.xm*info.ym*info.zm;
                 
      rstarts= new PetscInt[numtasks+1];
      
      ierr=PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);
      
      DMCreateGlobalVector(da, &x);
      
      DMCreateGlobalVector(da, &b);
      
      DMCreateGlobalVector(da, &u);
      
      i_map = new PetscInt[nbdof];
      
      for (PetscInt ii=0; ii< nbdof; ++ii) {
          i_map[ii]=ii;
      }
      
     
      AOApplicationToPetsc(ao, nbdof, i_map);
      
      VecGetOwnershipRange(x,&rstart,&rend);
      VecGetOwnershipRanges(x,&rstarts);
      
//      std::cout << "Rank\t" << rank << "\t" << rstart << "\t" << rend-1 << std::endl;
      
      nlocal=rend-rstart;
      
      std::vector<My_Int> lindices;
      lindices.reserve(N_lattice_sites);
      
      for (My_Int kk=ys; kk<ye; ++kk) {
          for (My_Int jj=xs; jj<xe; ++jj){
              for (My_Int ii=rs; ii<re; ++ii) {
                  lindices.push_back(ii+(Np-1)*jj+(Np-1)*Nx*kk);
              }}}
      
      
//      std::cout << "Rank\t" <<rank  << "\tapp start\t" << i_map[N_functions*lindices[0]] << "\tapp end\t" << i_map[N_functions*lindices[lindices.size()-1]+N_functions-1] << std::endl;
      
      My_Int temp_int;
      
      
      omp_start = omp_get_wtime();
      
      ConstructLinearOp(TriplVec, B, eom, b_cond, l_eom, l_b_cond,f,D,Rgrid, Xgrid, Ygrid, N_functions,params, lindices, i_map , tole );

      omp_end = omp_get_wtime();
      
      time = omp_end - omp_start;
      
      MPI_Reduce(&time,&mtime,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);
      
            if (rank==0) {
      std::cout<<"Got linear operator in: "<< mtime <<'\n';
            }
      
      omp_start = omp_get_wtime();
      
      long long nz=TriplVec.size();
    
//      std::cout << "Rank: " << rank << "\tReserved:" << (long long)(nvars)* (long long)(nbdof)/(long long)numtasks << "\tNeeded: " << nz << std::endl;

      My_Int *dnz;
      My_Int *onz;
      My_Int *totnz;
      My_Int *tdnz;
      My_Int *tonz;
      
      dnz=new My_Int[nbdof];
      onz=new My_Int[nbdof];
      totnz=new My_Int[nbdof];
      tdnz=new My_Int[nbdof];
      tonz=new My_Int[nbdof];
      
#pragma omp parallel for
for (My_Int ii=0; ii<nbdof; ++ii) {
          dnz[ii]=0; onz[ii]=0;
          tdnz[ii]=0; tonz[ii]=0;
          totnz[ii]=0;
      }

      
My_Int *ldnz, *ltotnz;
      
#pragma omp parallel
{
    const int nthreads = omp_get_num_threads();
    const int ithread = omp_get_thread_num();
    
#pragma omp single
{
    ldnz= new My_Int[nthreads*nbdof];
    ltotnz= new My_Int[nthreads*nbdof];
}
    

#pragma omp for
for (My_Int ii=0; ii<nthreads*nbdof; ++ii) {
        ldnz[ii]=0;
        ltotnz[ii]=0;
    }
    
#pragma omp for
for (long long ii=0; ii<nz; ++ii) {
//          if((rstart-1< TriplVec[ii].col() && TriplVec[ii].col()  < rend && rstart-1< TriplVec[ii].row() && TriplVec[ii].row()  <rend )){
//              ldnz[TriplVec[ii].row()+ithread*nbdof]+=1;
//          }
    if((rstart-1< TriplVec[ii].col() && TriplVec[ii].col()  < rend )){
        ldnz[TriplVec[ii].row()+ithread*nbdof]+=1;
    }
         ltotnz[TriplVec[ii].row()+ithread*nbdof]+=1;
    
      }
    
#pragma omp for
for (My_Int ii=rstart; ii<rend; ++ii) {
            for (My_Int nn=0; nn<nthreads; nn++) {
            dnz[ii]+=ldnz[ii+nn*nbdof];
            totnz[ii]+=ltotnz[ii+nn*nbdof];
            }
        }
#pragma omp single
{
    delete[] ldnz;
    delete[] ltotnz;
}
    
}
      
      My_Int N_non_empty_rows(0);
      
#pragma omp parallel for reduction(+:N_non_empty_rows)
      for (My_Int ii=rstart; ii<rend; ii++) {
          if (totnz[ii]) {
          onz[ii]=totnz[ii]-dnz[ii];
            N_non_empty_rows++;
          }
      }
      
      My_Int **Col_Indices;
      PetscScalar **Vals;
      
      Col_Indices= new My_Int*[nbdof];
      Vals= new PetscScalar*[nbdof];
      
      My_Int temp_int1(0);
      
#pragma omp parallel private(temp_int1)
{
#pragma omp for
for (My_Int ii=0; ii<nbdof; ii++) {
          temp_int1=totnz[ii];
          Col_Indices[ii]=0;
          Vals[ii]=0;
          if (temp_int1) {
              Col_Indices[ii]=new My_Int[temp_int1];
              Vals[ii]=new PetscScalar[temp_int1];
          }
      }
      }
      
      My_Int *temp_ints;
      
      temp_ints=new My_Int[nbdof];

#pragma omp parallel for
for (My_Int ii=0; ii<nbdof; ii++) {
          temp_ints[ii]=0;
      }
      
      omp_lock_t *lock;

      lock=new omp_lock_t[nlocal];
      
      
#pragma omp parallel for shared(lock)
      for (My_Int ii=0; ii<nlocal; ii++) {
          omp_init_lock(lock+ii);
      }
      

#pragma omp parallel for private(ti,tj, temp_int1) shared(lock, TriplVec, Col_Indices, Vals, temp_ints)
for(My_Int ii=0; ii<nz; ii++) {
          ti=TriplVec[ii].row();
          tj=TriplVec[ii].col();
        omp_set_lock(&(lock[ti-rstart]));
          temp_int1=temp_ints[ ti ]++;
        omp_unset_lock(&(lock[ti-rstart]));
        Col_Indices[ ti ][ temp_int1 ] = tj;
        Vals[ ti ][ temp_int1 ] = TriplVec[ii].value();
      }
      
      delete[] temp_ints;
      
      TriplVec.erase(TriplVec.begin(), TriplVec.end());
      
#pragma omp parallel for shared(lock)
      for (My_Int ii=0; ii<nlocal; ii++) {
          omp_destroy_lock(lock+ii);
      }
      
      delete []lock;
      
      omp_end = omp_get_wtime();
      
      time = omp_end - omp_start;
      
      MPI_Reduce(&time,&mtime,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);
      
      
      if (rank==0) {
          std::cout<<"Did loops in: "<< mtime  <<'\n';}
      
      MPI_Allreduce(dnz,tdnz,nbdof,MPI_LONG_LONG,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(onz,tonz,nbdof,MPI_LONG_LONG,MPI_SUM,PETSC_COMM_WORLD);
      
      My_Int av1(0),av2(0);
      
     MatCreateAIJ(PETSC_COMM_WORLD, nlocal,nlocal ,nbdof,nbdof, NULL, tdnz+rstart , NULL , tonz+rstart, &A);
      

      omp_start = omp_get_wtime();
      
      PetscInt i1,i2,i3;
      
      My_Int tii;
      
#pragma omp parallel private(tii, i3, i2, i1, ti, tvalue)
      {
#pragma omp for
          for(My_Int ii=0; ii<lindices.size(); ii++){
              tii=lindices[ii];
              i3= tii/((Np-1)*Nx);
              i2= (tii%((Np-1)*Nx))/(Np-1);
              i1=1+(tii%((Np-1)*Nx))%(Np-1);
              for(My_Int l=0; l<N_functions; l++){
                  ti=i_map[i3*N_functions*(Np-1)*Nx +i2*N_functions*(Np-1) +(i1-1)*N_functions +l];
              if (totnz[ti]) {
                  MatSetValues(A,1,&ti, totnz[ti] ,Col_Indices[ti],Vals[ti],INSERT_VALUES);
                  delete[] Col_Indices[ti];
                  delete[] Vals[ti];
              }
                  tvalue= B[ti];
                  VecSetValues(b,1, &ti,&tvalue,INSERT_VALUES);
              }
          }
      }

      
      delete[] Col_Indices;
      delete[] Vals;
      
      
      time = omp_end - omp_start;
      
      MPI_Reduce(&time,&mtime,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);
      
      
      MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
      VecAssemblyBegin(b);
      MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
      VecAssemblyEnd(b);
      
      omp_end = omp_get_wtime();
      
      time = omp_end - omp_start;
      
      MPI_Reduce(&time,&mtime,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);
      
      if (rank==0) {
          std::cout<<"Assembled linear operator in: "<< mtime <<'\n';}

      
      KSPCreate(PETSC_COMM_WORLD,&ksp);
      KSPGetPC(ksp,&pc);
      KSPSetOperators(ksp,A,A);
      KSPSetDM(ksp,da);
      KSPSetDMActive(ksp,PETSC_FALSE);
//      KSPSetType(ksp, KSPBCGS);
//      PCSetType(pc,PCBJACOBI);
      PCSetType(pc,PCMG);

      KSPSetTolerances(ksp,1.e-2,1.e-12,1.e+8,3000);
      
      PetscBool tf;
      
      PCSetFromOptions(pc);
      KSPSetFromOptions(ksp);
      
      
      PetscObjectTypeCompare((PetscObject)pc, PCMG,&tf);
      
      if (tf) {
          PCMGSetupViaCoarsen(pc,da);}
      
      KSPSolve(ksp,b,x);
      
      KSPConvergedReason conv;
      KSPGetConvergedReason(ksp,&conv);

      KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
      
      ierr = MatMult(A,x,u);
      ierr = VecAXPY(u,-1.,b);
      ierr = VecNorm(u,NORM_2,&norm);
      ierr = KSPGetIterationNumber(ksp,&its);
      
      PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);
      
      VecCreateSeq(MPI_COMM_SELF,nbdof,&X);
      VecSetType(X,VECSEQ);
      VecSetSizes(X,nbdof,nbdof);
      
      

      ISCreateStride(PETSC_COMM_WORLD ,nbdof,0,1, &iX);
      
      VecScatterCreate(x, NULL, X, iX, &Scatter);
      
      VecScatterBegin(Scatter,x,X,INSERT_VALUES,SCATTER_FORWARD);
      VecScatterEnd(Scatter,x,X,INSERT_VALUES,SCATTER_FORWARD);
      
    PetscScalar *xx;
    VecGetArray(X,&xx);
      
      if(conv>0){
      for (My_Int l3=0; l3<N_functions; l3++) {
          for (My_Int i1=1; i1 < Np; i1++) {
              for (My_Int i2=0; i2<Nx; i2++) {
                  for (My_Int i3=0; i3<Ny; i3++){
                      f[l3](i1,i2,i3)=f[l3](i1,i2,i3)-xx[i_map[i3*N_functions*(Np-1)*Nx +i2*N_functions*(Np-1) +(i1-1)*N_functions +l3]];
                  }}}}
      
      
      for (My_Int l3=0; l3<N_functions; ++l3) {
          f[l3].update(D);
      }
      }
      else{
          std::cout << "Solver did not converge!" << std::endl;
          exit(1);
      }
      
      ierr = VecDestroy(&x); ierr = VecDestroy(&u);
      ierr = VecDestroy(&b); ierr = MatDestroy(&A);
      ierr = VecDestroy(&X); ISDestroy(&iX);
      ierr = KSPDestroy(&ksp);
//      AODestroy(&ao);
      DMDestroy(&da);
      VecScatterDestroy(&Scatter);
//      delete[] rstarts;
      delete[] dnz;
      delete[] onz;
      delete[] totnz;
      delete[] tdnz;
      delete[] tonz;
      delete[] i_map;

  }
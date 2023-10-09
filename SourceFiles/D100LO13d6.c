
elem[197]+=
(-2*a2.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]) - 
    2*H*L1*pow(1 + pow(R0,2),0.25)*(1 - pow(1 - rgrid[i],2))*
     (1 - rgrid[i]) - (2*a0.diff(d001,i,j,k)*L1*
       pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],3)*Qt1(i,j,k))/L2 + 
    (2*a1.diff(d001,i,j,k)*L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
       sqrt(1 + pow(R0,2)))/L2)*D[-1 + d100](i,j,k,i1,j1,k1)
;

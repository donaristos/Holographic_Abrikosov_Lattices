
elem[182]+=
(-2*a1.diff(d001,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]) + 
    (2*a2.diff(d010,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
       (1 - rgrid[i]))/(L1*sqrt(1 + pow(R0,2))) + 
    (2*H*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]))/
     pow(1 + pow(R0,2),0.25) - 
    (2*a0.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],3)*Qt2(i,j,k))/(L1*pow(1 + pow(R0,2),0.25)))*
  D[-1 + d100](i,j,k,i1,j1,k1)
;

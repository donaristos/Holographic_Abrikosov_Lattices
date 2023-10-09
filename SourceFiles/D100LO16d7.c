
elem[246]+=
(-2*h2.diff(d001,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]) + 
    2*q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
     (-((L2*a2(i,j,k))/pow(1 + pow(R0,2),0.25)) + 
       L2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        Qt2(i,j,k)))*D[-1 + d100](i,j,k,i1,j1,k1)
;

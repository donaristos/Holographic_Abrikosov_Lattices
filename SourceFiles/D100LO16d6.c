
elem[245]+=
(-2*h2.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]) + 
    2*H*L1*L2*q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
     ygrid[k] + (q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
       (-4*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k) + 
         4*L1*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
          Qt1(i,j,k)))/2.)*D[-1 + d100](i,j,k,i1,j1,k1)
;

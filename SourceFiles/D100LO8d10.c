
elem[121]+=
(-2*Qrr.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
    D[-1 + d100](i,j,k,i1,j1,k1))/
  (L1 + L1*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))
;

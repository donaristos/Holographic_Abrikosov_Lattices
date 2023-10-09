
elem[38]+=
(-2*Qrr.diff(d001,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
    (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*D[-1 + d100](i,j,k,i1,j1,k1))/
  (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))
;

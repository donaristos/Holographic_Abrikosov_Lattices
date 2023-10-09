
elem[205]+=
(-2*Qrr.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*(1 - rgrid[i])*
    D[-1 + d100](i,j,k,i1,j1,k1))/
  (L2 + L2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))
;


elem[250]+=
2*q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],3)*
  ((1 - pow(1 - rgrid[i],2))*
     (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k))*
  D[-1 + d100](i,j,k,i1,j1,k1)
;

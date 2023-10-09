
elem[188]+=
(-8*Qr2.diff(d010,i,j,k)*L2*pow(1 - rgrid[i],2))/
   (L1*sqrt(1 + pow(R0,2))) + 
  (8*Qrr.diff(d010,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*Qr2(i,j,k))/
   (L1*sqrt(1 + pow(R0,2))*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))
;

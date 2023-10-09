
elem[203]+=
(-8*Qr1.diff(d001,i,j,k)*L1*pow(1 - rgrid[i],2)*sqrt(1 + pow(R0,2)))/
   L2 + (8*Qrr.diff(d001,i,j,k)*L1*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))/
   (L2*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))
;

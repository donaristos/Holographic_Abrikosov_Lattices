
elem[193]+=
(-((a2.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)))/
       (2 + 2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
    (2*a2.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr1(i,j,k))/
     (2 + 2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
    (2*H*L1*pow(1 + pow(R0,2),0.25)*
       pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*Qr1(i,j,k))/
     (2 + 2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
    (2*ar.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*(1 - rgrid[i]))/
     (2*L2 + 2*L2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
    (a0.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
       pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*
       ((1 - pow(1 - rgrid[i],2))*
          (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
         Qtr(i,j,k)))/(L2 + L2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
    (a1.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))/
     (L2 + L2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))*
  D[-1 + d100](i,j,k,i1,j1,k1)
;

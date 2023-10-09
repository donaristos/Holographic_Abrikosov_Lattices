
elem[254]+=
(4*H*L1*L2*q*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*ygrid[k]*
     Qr1(i,j,k) + (q*(1 - rgrid[i])*
       (4*ar(i,j,k) - 4*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
          pow(1 - pow(1 - rgrid[i],2),2)*Qr1(i,j,k) - 
         (4*L2*a2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*Qr2(i,j,k))/
          pow(1 + pow(R0,2),0.25) + 
         4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
          pow(1 - rgrid[i],2)*
          ((1 - pow(1 - rgrid[i],2))*
             (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
            Qtr(i,j,k))))/(1 - pow(1 - rgrid[i],2)))*
  D[-1 + d100](i,j,k,i1,j1,k1)
;

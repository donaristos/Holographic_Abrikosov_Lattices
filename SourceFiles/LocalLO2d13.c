
elem[28]+=
(-4*h2.diff(d100,i,j,k)*L2*q*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*Qr2(i,j,k)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   pow(1 + pow(R0,2),0.25) + 
  (4*h1.diff(d100,i,j,k)*L2*q*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*Qr2(i,j,k)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   pow(1 + pow(R0,2),0.25) + 
  (8*h2.diff(d010,i,j,k)*L2*q*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*Qr1(i,j,k)*
     Qr2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   pow(1 + pow(R0,2),0.25) - 
  (8*h1.diff(d010,i,j,k)*L2*q*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*Qr1(i,j,k)*
     Qr2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   pow(1 + pow(R0,2),0.25) - 
  (8*H*L1*pow(L2,2)*pow(q,2)*(pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*ygrid[k]*
     Qr1(i,j,k)*Qr2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   pow(1 + pow(R0,2),0.25) + 
  (8*h2.diff(d001,i,j,k)*L2*q*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
     pow(Qr2(i,j,k),2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   pow(1 + pow(R0,2),0.25) - 
  (8*h1.diff(d001,i,j,k)*L2*q*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
     pow(Qr2(i,j,k),2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   pow(1 + pow(R0,2),0.25) + 
  (2*L2*pow(q,2)*(pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
     (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qr2(i,j,k)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (-4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*
        (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
          4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*pow(1 - rgrid[i],2)*
           (-((1 - pow(1 - rgrid[i],2))*
                (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))) + 
             Qtr(i,j,k)) + 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))\
)*sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
   (sqrt(1 + pow(R0,2))*sqrt(1 + 
       pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))
;


elem[176]+=
((a1.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)))/
     (2 + 2*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
    (2*a1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr2(i,j,k))/
     (2 + 2*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
    (2*H*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*Qr2(i,j,k))/
     (pow(1 + pow(R0,2),0.25)*
       (2 + 2*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
    (4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
       Qt1(i,j,k))/
     (pow(1 + pow(R0,2),0.25)*
       (2 + 2*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
    (4*a0.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],2)*Qt1(i,j,k))/
     (pow(1 + pow(R0,2),0.25)*
       (4 + 4*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
    (8*a0.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
       pow(1 - rgrid[i],3)*Qr2(i,j,k)*Qt1(i,j,k))/
     (pow(1 + pow(R0,2),0.25)*
       (4 + 4*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
    (ar.diff(d010,i,j,k)*(1 - rgrid[i]))/
     (pow(1 + pow(R0,2),0.25)*
       (L1 + L1*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
    (a0.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],3)*((1 - pow(1 - rgrid[i],2))*
          (L1*Qr1(i,j,k)*Qt1(i,j,k) - L2*Qr2(i,j,k)*Qt2(i,j,k)) + 
         Qtr(i,j,k)))/
     (pow(1 + pow(R0,2),0.25)*
       (L1 + L1*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
    (2*a2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr2(i,j,k))/
     (sqrt(1 + pow(R0,2))*(2*L1 + 
         2*L1*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))))*
  D[-1 + d100](i,j,k,i1,j1,k1)
;

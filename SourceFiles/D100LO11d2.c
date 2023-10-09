
elem[161]+=
((a0(i,j,k)*(1 - pow(1 - rgrid[i],2)))/
     ((1 - rgrid[i])*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) - 
    (a0.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)))/
     (2 + 2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
    (2*a0.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr1(i,j,k))/
     (2 + 2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
    (2*a0.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr2(i,j,k))/
     (2 + 2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))*
  D[-1 + d100](i,j,k,i1,j1,k1)
;

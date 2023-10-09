
elem[166]+=
(4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2)) - 
    2*a0.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]) + 
    4*a0.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
     pow(1 - rgrid[i],2)*Qr1(i,j,k) + 
    8*a0.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
     pow(1 - rgrid[i],2)*Qr2(i,j,k))*D[-1 + d001](i,j,k,i1,j1,k1)
;

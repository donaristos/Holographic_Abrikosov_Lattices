
elem[229]+=
(4*h1.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
     pow(1 - rgrid[i],2)*Qr2(i,j,k) + 
    4*H*L1*L2*q*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
     pow(1 - rgrid[i],2)*ygrid[k]*Qr2(i,j,k) + 
    q*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
     Qr2(i,j,k)*(-4*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k) + 
       4*L1*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        Qt1(i,j,k)))*D[-1 + d001](i,j,k,i1,j1,k1)
;

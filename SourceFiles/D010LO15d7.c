
elem[230]+=
(4*h1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
     pow(1 - rgrid[i],2)*Qr1(i,j,k) + 
    4*q*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
     Qr1(i,j,k)*(-((L2*a2(i,j,k))/pow(1 + pow(R0,2),0.25)) + 
       L2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        Qt2(i,j,k)))*D[-1 + d010](i,j,k,i1,j1,k1)
;

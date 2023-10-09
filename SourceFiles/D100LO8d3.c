
elem[114]+=
(-((Qt1.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)))/
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
    (2*Qtr.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]))/
     (L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
    (2*Qt1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr2(i,j,k))/
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
    (2*Qt2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr2(i,j,k))/
     (L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
    (2*(1 - rgrid[i])*Qt1(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))\
)*D[-1 + d100](i,j,k,i1,j1,k1)
;

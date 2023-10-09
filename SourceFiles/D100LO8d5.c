
elem[116]+=
(-((Qt2.diff(d100,i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
         (1 - pow(1 - rgrid[i],2)))/
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (2*Qtr.diff(d001,i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]))/
     (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (2*Qt2.diff(d010,i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*Qr1(i,j,k))/
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) - 
    (2*Qt1.diff(d001,i,j,k)*L1*
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*Qr1(i,j,k))/
     (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
    (2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*(1 - rgrid[i])*
       Qt2(i,j,k))/
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))*
  D[-1 + d100](i,j,k,i1,j1,k1)
;

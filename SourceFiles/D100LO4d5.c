
elem[52]+=
((2*B.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2))*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
     (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
    (4*B(i,j,k)*(1 - rgrid[i])*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
     (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
    (4*Qr1.diff(d010,i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
     (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
    (4*Qr2.diff(d001,i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
     (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
    (4*B.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr1(i,j,k)*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
     (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
    (4*B.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr2(i,j,k)*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
     (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
    (4*Qr2.diff(d010,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
       (1 - rgrid[i]))/
     (L1*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (4*Qr1.diff(d001,i,j,k)*L1*
       pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
       (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]))/
     (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))*
  D[-1 + d100](i,j,k,i1,j1,k1)
;

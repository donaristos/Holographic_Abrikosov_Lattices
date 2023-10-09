
elem[178]+=
((-4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
       Qt1(i,j,k))/
     (pow(1 + pow(R0,2),0.25)*
       (2 + 2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
    (4*a0.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],2)*Qt1(i,j,k))/
     (pow(1 + pow(R0,2),0.25)*
       (4 + 4*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
    (8*a0.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
       pow(1 - rgrid[i],3)*Qr1(i,j,k)*Qt1(i,j,k))/
     (pow(1 + pow(R0,2),0.25)*
       (4 + 4*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
    (8*a0.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
       pow(1 - rgrid[i],3)*Qr2(i,j,k)*Qt1(i,j,k))/
     (pow(1 + pow(R0,2),0.25)*
       (4 + 4*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))))*
  D[-1 + d100](i,j,k,i1,j1,k1)
;

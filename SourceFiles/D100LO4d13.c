
elem[60]+=
(-((a2.diff(d100,i,j,k)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
          2)*(1 - pow(1 - rgrid[i],2)))/
       (sqrt(1 + pow(R0,2))*pow(rh,2)*
         (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
         sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
    (2*ar.diff(d001,i,j,k)*pow(1 + 
         B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*(1 - rgrid[i]))/
     (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
       (L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (2*a2.diff(d010,i,j,k)*pow(1 + 
         B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
       pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*Qr1(i,j,k))/
     (sqrt(1 + pow(R0,2))*pow(rh,2)*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (2*H*L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
       pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*Qr1(i,j,k))/
     (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
    (2*a1.diff(d001,i,j,k)*L1*
       pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
       pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*Qr1(i,j,k))/
     (pow(rh,2)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
    (2*a0(i,j,k)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
       pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*Qt2(i,j,k))/
     (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (a0.diff(d100,i,j,k)*pow(1 + 
         B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
       pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*Qt2(i,j,k)\
)/(pow(1 + pow(R0,2),0.25)*pow(rh,2)*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
    (2*a0.diff(d010,i,j,k)*pow(1 + 
         B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
       pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],3)*
       Qr1(i,j,k)*Qt2(i,j,k))/
     (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (2*a0.diff(d001,i,j,k)*pow(1 + 
         B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
       pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*
       (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - Qtr(i,j,k)))/
     (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
       (L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))*
  D[-1 + d100](i,j,k,i1,j1,k1)
;

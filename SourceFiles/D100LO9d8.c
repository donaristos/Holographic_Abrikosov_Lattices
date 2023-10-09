
elem[135]+=
((2*Qrr.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr1(i,j,k))/
     (L2 + L2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
    (2*Qr1.diff(d001,i,j,k)*L1*(1 - pow(1 - rgrid[i],2))*
       (1 - rgrid[i])*(1 + pow(R0,2) + 
         2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/L2 + 
    (2*Qr2.diff(d010,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
       (1 - rgrid[i])*(1 + pow(R0,2) + 
         2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
     (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) - 
    (R.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)))/
     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (2*R.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr1(i,j,k))/
     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (2*R.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr2(i,j,k))/
     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (2*Qr1.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
       (R0 + pow(R0,3) + (1 + 3*pow(R0,2))*
          (1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
         pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)))/
     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
    (2*Qr2.diff(d001,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
       (R0 + pow(R0,3) + (1 + 3*pow(R0,2))*
          (1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
         pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)))/
     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
    (2*(1 - rgrid[i])*(R(i,j,k) + 
         B(i,j,k)*(R0 + pow(R0,3) + 
            (2 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
            pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3))))/
     (pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
    (B.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2))*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
     pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) + 
    (2*B.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr1(i,j,k)*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
     pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) + 
    (2*B.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr2(i,j,k)*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
     pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2))*
  D[-1 + d100](i,j,k,i1,j1,k1)
;

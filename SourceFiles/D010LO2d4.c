
elem[19]+=
((-4*B.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr1(i,j,k)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
     pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) + 
    (8*Qr1.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],2)*Qr1(i,j,k)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) - 
    (8*Qr2.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],2)*Qr1(i,j,k)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
    (8*B.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
       pow(1 - rgrid[i],2)*pow(Qr1(i,j,k),2)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
     pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) + 
    (8*B.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
       pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qr2(i,j,k)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
     pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) - 
    (8*Qr1.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],2)*Qr1(i,j,k)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (R0 + pow(R0,3) + (1 + 3*pow(R0,2))*
          (1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
         pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)))/
     (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (8*Qr2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],2)*Qr1(i,j,k)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (R0 + pow(R0,3) + (1 + 3*pow(R0,2))*
          (1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
         pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)))/
     (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (2*(1 - rgrid[i])*Qr1(i,j,k)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (-4*B(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
          (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
         ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
            (pow(L1,2)*pow(1 - pow(1 - rgrid[i],2),3)*
               ((pow(-1 + rgrid[i],2)*
                    (3*pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],2) + 
                      3*pow(H,2)*pow(-2 + rgrid[i],2)*
                       pow(rgrid[i],3) + 
                      pow(rh,2)*
                       (3*pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],2) + 
                        3*pow(mu,2)*pow(-2 + rgrid[i],2)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (2 + 6*rgrid[i] - 12*pow(rgrid[i],2) + 
                        4*pow(rgrid[i],3)))))/(4.*pow(rh,4)) + 
                 ((-1 + rgrid[i])*
                    (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      pow(rh,2)*
                       (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
                  (2.*pow(rh,4)))*
               (-pow(Qt1(i,j,k),2) + 
                 pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                  pow(Qt2(i,j,k),2)) + 
              4*(1 - rgrid[i])*
               (pow(L1,2)*((pow(1 - pow(1 - rgrid[i],2),2)*
                       pow(-1 + rgrid[i],2)*
                       (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                       pow(Qt1(i,j,k),2))/(4.*pow(rh,4)) - 
                    sqrt(1 + pow(R0,2))) + 
                 pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                  (-(pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(-1 + rgrid[i],2)*
                        (pow(H,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                        pow(Qt2(i,j,k),2))/(4.*pow(rh,4)) + 
                    sqrt(1 + pow(R0,2)))))*
            sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
          (pow(L1,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/
     pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2))*
  D[-1 + d010](i,j,k,i1,j1,k1)
;

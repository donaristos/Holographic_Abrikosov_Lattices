
elem[21]+=
(-8*Qr2.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (-1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
    (8*Qr1.diff(d001,i,j,k)*pow(L1,2)*
       pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
       (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/pow(L2,2) \
- (4*R.diff(d100,i,j,k)*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
     (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (8*R.diff(d010,i,j,k)*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
       Qr1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
     (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (8*R.diff(d001,i,j,k)*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
       Qr2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
     (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (4*B.diff(d100,i,j,k)*L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (R0 + pow(R0,3) + (1 + 3*pow(R0,2))*
          (1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
         pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)))/
     (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
    (8*B.diff(d010,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],2)*Qr1(i,j,k)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (R0 + pow(R0,3) + (1 + 3*pow(R0,2))*
          (1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
         pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)))/
     (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
    (8*B.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],2)*Qr2(i,j,k)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (R0 + pow(R0,3) + (1 + 3*pow(R0,2))*
          (1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
         pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)))/
     (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (4*L1*(1 - rgrid[i])*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (2*(1 - rgrid[i])*(-R(i,j,k) + 
            B(i,j,k)*(R0 + pow(R0,3) + 
               3*pow(R0,2)*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               3*R0*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(R(i,j,k),2) + 
               pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3))) - 
         (L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
            pow(-1 + rgrid[i],4)*
            pow(pow(H,2)*pow(-2 + rgrid[i],3)*
               pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
            Qr2(i,j,k)*Qt1(i,j,k)*
            ((1 - pow(1 - rgrid[i],2))*
               (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
              Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
              2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
              pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
          (4.*pow(rh,6)) + (L2*(1 - pow(1 - rgrid[i],2))*
            pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
            Qr2(i,j,k)*(pow(1 - pow(1 - rgrid[i],2),3)*
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
                  (2.*pow(rh,4)))*Qt1(i,j,k)*
               (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
              pow(1 - pow(1 - rgrid[i],2),2)*
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
                  (2.*pow(rh,4)))*Qt1(i,j,k)*Qtr(i,j,k) + 
              4*(1 - rgrid[i])*
               (L2*R0*Qr2(i,j,k) + L1*Qr1(i,j,k)*sqrt(1 + pow(R0,2))))*
            sqrt(1 + pow(R0,2) + 
              2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
              pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
          (4.*pow(rh,2))))/
     (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
    (8*Qr1.diff(d010,i,j,k)*L1*
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/L2 + 
    (8*Qr2.diff(d001,i,j,k)*L1*
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/L2)*
  D[-1 + d001](i,j,k,i1,j1,k1)
;

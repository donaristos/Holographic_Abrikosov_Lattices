
elem[114]+=
((2*Qt1.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr1(i,j,k))/
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
    (4*Qtr.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],2)*Qr1(i,j,k))/
     (L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
    (4*Qt2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),3)*
       pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qr2(i,j,k))/
     (L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
    (16*Qt1.diff(d001,i,j,k)*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*
       pow(1 - rgrid[i],2)*((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
            pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
            (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
            Qr2(i,j,k))/(4.*pow(rh,2)) - 
         2*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))/
     (L1*L2*pow(-1 + rgrid[i],2)*
       (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
    (16*Qt1.diff(d010,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
       pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
     (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
          pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
    (16*Qt2.diff(d001,i,j,k)*pow(rh,4)*
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
     (L1*L2*pow(-1 + rgrid[i],2)*
       (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
    ((1 - rgrid[i])*(4*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*Qt1(i,j,k) \
+ ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qtr(i,j,k)*
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
               (pow(Qt1(i,j,k),2)*
                  (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),
                     2)) + pow(1 + 
                    B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                  pow(Qt2(i,j,k),2)*
                  (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),
                     2)) - 2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                  Qt1(i,j,k)*Qt2(i,j,k)*
                  (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                  sqrt(1 + pow(R0 + 
                      (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
              4*(1 - rgrid[i])*
               (pow(L1,2)*(1 + 
                    pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
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
                        pow(Qt1(i,j,k),2))/(4.*pow(rh,4)) + 
                    sqrt(1 + pow(R0,2))) + 
                 pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                  (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),
                     2))*(-(pow(1 - pow(1 - rgrid[i],2),2)*
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
                    sqrt(1 + pow(R0,2))) - 
                 2*pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                  (R0 - (pow(1 - pow(1 - rgrid[i],2),2)*
                       pow(-1 + rgrid[i],2)*
                       (pow(H,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                       Qt1(i,j,k)*Qt2(i,j,k))/(4.*pow(rh,4)))*
                  (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                  sqrt(1 + pow(R0 + 
                      (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))))/
          (pow(L1,3)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
            sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))))/
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))*
  D[-1 + d010](i,j,k,i1,j1,k1)
;

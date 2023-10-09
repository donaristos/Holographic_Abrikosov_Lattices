
elem[220]+=
((4*Qr2.diff(d100,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),3)*
       (1 - rgrid[i])*Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) - 
    (2*Qrr.diff(d100,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),4)*
       (1 - rgrid[i])*pow(Qr2(i,j,k),2))/
     (pow(1 + pow(R0,2),0.25)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
    (8*R.diff(d100,i,j,k)*pow(rh,2)*
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
     (L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
       (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (8*B.diff(d100,i,j,k)*pow(rh,2)*
       pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
     (L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
       (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
    (8*Q.diff(d100,i,j,k)*pow(rh,4)*
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
     (L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
       (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
    ((1 - pow(1 - rgrid[i],2))*
       (4*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
          pow(Qr2(i,j,k),2) + 
         2*pow(L2,2)*(1 - pow(1 - rgrid[i],2))*
          (-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
          pow(Qr2(i,j,k),2) + 
         8*pow(L2,2)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
          pow(Qr2(i,j,k),2) + 
         (8*pow(L2,2)*pow(rh,4)*pow(1 - pow(1 - rgrid[i],2),2)*
            (1 - rgrid[i])*((pow(-1 + rgrid[i],2)*
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
               (2.*pow(rh,4)))*pow(Qr2(i,j,k),2))/
          (pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))))) - 
         (4*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
            pow(1 - rgrid[i],2)*pow(Qr2(i,j,k),2)*Qrr(i,j,k))/
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
         (16*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
            pow(1 - rgrid[i],2)*
            (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*R(i,j,k)*
            (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
          (pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
            (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
            sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) \
+ (16*pow(rh,2)*B(i,j,k)*pow(1 - rgrid[i],2)*
            (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
            sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
          (pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
            (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
         (32*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
            pow(1 - rgrid[i],2)*
            (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
            sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
          ((1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
            (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
         (16*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
            pow(1 - rgrid[i],2)*Q(i,j,k)*
            (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
            sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
          (pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
            pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))))/
     (L2*pow(1 + pow(R0,2),0.25)))*D[-1 + d001](i,j,k,i1,j1,k1)
;

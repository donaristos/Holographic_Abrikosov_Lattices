
elem[57]+=
(-((Qtr.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*
         pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
         (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
           pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
               pow(rgrid[i],3) + 
              4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                 4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
         (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
       (pow(L1,2)*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
         sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
    (Qt1.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
       pow(-1 + rgrid[i],2)*
       (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
     (2.*pow(rh,4)*(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (Qt2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
       (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qr2(i,j,k)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
     (pow(L1,2)*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
    (Qt1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
       (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qr2(i,j,k)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
     (pow(rh,4)*(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    ((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
       (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       ((2*(1 - rgrid[i])*Qt1(i,j,k)*
            (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
         pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
          (-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                 (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                   pow(rh,2)*
                    (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                         4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
               pow(rh,4)) + 
            (1 - pow(1 - rgrid[i],2))*
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
                  (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                         4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
                (2.*pow(rh,4))))*
          ((1 - pow(1 - rgrid[i],2))*
             (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
            Qtr(i,j,k))*(L2*R0*Qr2(i,j,k) + 
            L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*R(i,j,k) + 
            L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))\
)/(2.*L1*pow(rh,4)*sqrt(1 + pow(R0 + 
           (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))*
  D[-1 + d010](i,j,k,i1,j1,k1)
;

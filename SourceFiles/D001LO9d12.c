
elem[139]+=
((8*L1*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
       pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*Qr1(i,j,k)\
)/((L2*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
        (2.*pow(rh,2)) + (L2*(1 - pow(1 - rgrid[i],2))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qtt(i,j,k))/
        (2.*pow(rh,2))) - (8*a0.diff(d100,i,j,k)*L1*
       pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],3)*Qr1(i,j,k))/
     ((L2*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/pow(rh,2) \
+ (L2*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qtt(i,j,k))/
        pow(rh,2)) - (64*a0.diff(d001,i,j,k)*
       pow(1 + pow(R0,2),0.25)*pow(rh,4)*
       (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],4)*
       (-(L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qr2(i,j,k))/(4.*pow(rh,2)) + 
         (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))/
     (pow(L2,2)*pow(-1 + rgrid[i],4)*
       pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
    (64*a0.diff(d010,i,j,k)*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
       (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],4)*
       ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
         (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
            pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
            (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
            pow(Qr1(i,j,k),2)*
            sqrt(1 + pow(R0,2) + 
              2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
              pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
          (4.*pow(rh,2))))/
     (L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       pow(-1 + rgrid[i],4)*pow(pow(H,2)*pow(-2 + rgrid[i],3)*
          pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))*
  D[-1 + d001](i,j,k,i1,j1,k1)
;

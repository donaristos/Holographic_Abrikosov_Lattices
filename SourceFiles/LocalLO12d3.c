
elem[178]+=
(4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
     Q(i,j,k)*Qt1(i,j,k))/
   (pow(1 + pow(R0,2),0.25)*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  Q.diff(d100,i,j,k)*((2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qt1(i,j,k))/
      (pow(1 + pow(R0,2),0.25)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (a0.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qt1(i,j,k))/
      (pow(1 + pow(R0,2),0.25)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (2*a0.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],3)*Qr1(i,j,k)*Qt1(i,j,k))/
      (pow(1 + pow(R0,2),0.25)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (2*a0.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],3)*Qr2(i,j,k)*Qt1(i,j,k))/
      (pow(1 + pow(R0,2),0.25)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  (32*a1.diff(d011,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   (L1*L2*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  a0.diff(d001,i,j,k)*R.diff(d001,i,j,k)*
   ((-4*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
        Qr2(i,j,k)*(L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
           Qt1(i,j,k) - Qtr(i,j,k)))/
      (pow(1 + pow(R0,2),0.25)*
        (L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (16*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
        ((L2*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
             Qtr(i,j,k))/(4.*pow(rh,2)) + 
          Qt1(i,j,k)*(-(L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                 pow(-1 + rgrid[i],2)*
                 (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                   pow(rh,2)*
                    (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                         4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                 (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                 Qr2(i,j,k))/(4.*pow(rh,2)) + 
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))))/
      (pow(L2,2)*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  (a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
     Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
     (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   (pow(1 + pow(R0,2),0.25)*pow(rh,4)*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  a1.diff(d001,i,j,k)*(-((pow(1 - pow(1 - rgrid[i],2),3)*
          pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
          Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
            L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (L2*pow(rh,4)*pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) \
+ (pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(Qt1(i,j,k),2)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(rh,4)*(L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  (4*pow(L1,2)*L2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*Qt1(i,j,k)*
     (-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
       (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (Qt1(i,j,k)*(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (4.*pow(rh,4))))/
   (pow(1 + pow(R0,2),0.25)*
     (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  R.diff(d010,i,j,k)*((16*a1.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (16*a2.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
      (pow(L1,2)*sqrt(1 + pow(R0,2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (16*a0.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k))/
      (pow(L1,2)*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (16*a1.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (16*a2.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*sqrt(1 + pow(R0,2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     a0.diff(d001,i,j,k)*((-4*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
           Qr1(i,j,k)*(L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
              Qt1(i,j,k) - Qtr(i,j,k)))/
         (pow(1 + pow(R0,2),0.25)*
           (L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (16*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
           ((L1*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
                  Qtr(i,j,k)))/(4.*pow(rh,2)) + 
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))/
         (L1*L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     (ar(i,j,k)*(1 - rgrid[i])*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
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
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
              (2.*pow(rh,4)))*
           (pow(Qt1(i,j,k),2)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              pow(Qt2(i,j,k),2)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt1(i,j,k)*
              Qt2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
          4*(1 - rgrid[i])*(pow(L1,2)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
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
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
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
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
      (pow(L1,3)*pow(1 + pow(R0,2),0.25)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  a0.diff(d001,i,j,k)*((4*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],4)*Q(i,j,k)*Qr2(i,j,k)*Qt1(i,j,k))/
      (pow(1 + pow(R0,2),0.25)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
             Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (4*pow(L1,2)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],4)*Qt1(i,j,k)*
        (L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          ((1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
               (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
                  Qtr(i,j,k))*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (4.*pow(rh,4))))/
      (pow(1 + pow(R0,2),0.25)*
        (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  Qrr.diff(d010,i,j,k)*((8*a1.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (8*a2.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(L1,2)*sqrt(1 + pow(R0,2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     a0.diff(d010,i,j,k)*((2*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],4)*Qr1(i,j,k)*
           ((1 - pow(1 - rgrid[i],2))*
              (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
             Qtr(i,j,k)))/
         (L1*pow(1 + pow(R0,2),0.25)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) - 
        (8*pow(rh,4)*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],4)*
           ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (-(L1*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                    (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      pow(rh,2)*
                       (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                    (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                    Qr1(i,j,k)*
                    ((1 - pow(1 - rgrid[i],2))*
                       (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                        L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)))/
                 (4.*pow(rh,2)) + 
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
                 (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L1,2)*pow(1 + pow(R0,2),0.25)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))) + 
  Qtt.diff(d010,i,j,k)*((8*a1.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
     (8*a2.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(L1,2)*sqrt(1 + pow(R0,2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
     (16*a0.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
     a0.diff(d010,i,j,k)*((-2*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],4)*Qr1(i,j,k)*
           ((1 - pow(1 - rgrid[i],2))*
              (L1*Qr1(i,j,k)*Qt1(i,j,k) - L2*Qr2(i,j,k)*Qt2(i,j,k)) + 
             Qtr(i,j,k)))/
         (L1*pow(1 + pow(R0,2),0.25)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
        (8*pow(rh,4)*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],4)*
           (-((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              ((L1*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                   (pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   Qr1(i,j,k)*
                   ((1 - pow(1 - rgrid[i],2))*
                      (L1*Qr1(i,j,k)*Qt1(i,j,k) - 
                        L2*Qr2(i,j,k)*Qt2(i,j,k)) + Qtr(i,j,k)))/
                 (4.*pow(rh,2)) + 
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
                 (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L1,2)*pow(1 + pow(R0,2),0.25)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))) - 
  (8*a1.diff(d001,i,j,k)*Qrr.diff(d001,i,j,k)*pow(rh,2)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (pow(L2,2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
  (16*a1.diff(d020,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
  (16*a1.diff(d002,i,j,k)*pow(rh,2)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (pow(L2,2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
  (16*a0.diff(d001,i,j,k)*B.diff(d001,i,j,k)*pow(rh,2)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (pow(L2,2)*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  Q.diff(d010,i,j,k)*((-4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qt1(i,j,k))/
      (pow(1 + pow(R0,2),0.25)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (64*a1.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)) + 
     a0.diff(d010,i,j,k)*((4*pow(1 - pow(1 - rgrid[i],2),5)*
           pow(1 - rgrid[i],4)*pow(Qr1(i,j,k),2)*Qt1(i,j,k))/
         (pow(1 + pow(R0,2),0.25)*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
        (32*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],4)*Qt1(i,j,k)*
           ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
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
         (pow(L1,2)*pow(1 + pow(R0,2),0.25)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     a0.diff(d001,i,j,k)*((4*pow(1 - pow(1 - rgrid[i],2),5)*
           pow(1 - rgrid[i],4)*Qr1(i,j,k)*Qr2(i,j,k)*Qt1(i,j,k))/
         (pow(1 + pow(R0,2),0.25)*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
        (32*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],4)*
           (Qt1(i,j,k)*((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                   pow(-1 + rgrid[i],2)*
                   (pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   Qr1(i,j,k)*Qr2(i,j,k))/(4.*pow(rh,2)) - 
                2*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                 (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3))) - 
     (2*ar(i,j,k)*(1 - rgrid[i])*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
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
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
              (2.*pow(rh,4)))*
           (pow(Qt1(i,j,k),2)*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              pow(Qt2(i,j,k),2)*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt1(i,j,k)*
              Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
          4*(1 - rgrid[i])*(pow(L1,2)*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
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
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
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
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
      (pow(1 + pow(R0,2),0.25)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (32*a1.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)) + 
     (32*a2.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*L2*sqrt(1 + pow(R0,2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3))) + 
  Qt1.diff(d010,i,j,k)*((2*pow(rh,4)*
        pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        (((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           pow(rh,4) - (1 - pow(1 - rgrid[i],2))*
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
              (2.*pow(rh,4))))*Qr1(i,j,k)*
        (4*a0(i,j,k)*pow(1 - rgrid[i],2) + 
          (ar(i,j,k)*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             ((1 - pow(1 - rgrid[i],2))*
                (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
               Qtr(i,j,k)))/(4.*pow(rh,6))))/
      (pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (32*a0.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (2*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        (-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
             pow(rh,4)) + (1 - pow(1 - rgrid[i],2))*
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
        (4*L1*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          (ar(i,j,k)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
- (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (-(L1*(1 - pow(1 - rgrid[i],2))*
                      pow(-1 + rgrid[i],2)*
                      (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                      Qr1(i,j,k)*
                      ((1 - pow(1 - rgrid[i],2))*
                        (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                        L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)))/
                   (4.*pow(rh,2)) + 
                  (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                   Qt2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           pow(rh,4)))/
      (L1*pow(1 + pow(R0,2),0.25)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (16*a0.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,2)*pow(1 + pow(R0,2),0.25)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  Qtt.diff(d001,i,j,k)*((-8*a1.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
     (16*a0.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))) + 
  Qt1.diff(d100,i,j,k)*((4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i]))/
      (pow(1 + pow(R0,2),0.25)*
        (2 + 2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (a2.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (4.*sqrt(1 + pow(R0,2))*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (ar.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (2.*L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (a1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*L2*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (a2.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*L1*sqrt(1 + pow(R0,2))*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (a0.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],3)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          Qtr(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*L1*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (2*pow(L1,2)*L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*(-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
          (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (Qt1(i,j,k)*(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
- (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (4.*pow(rh,4))))/
      (pow(1 + pow(R0,2),0.25)*
        (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     a0.diff(d100,i,j,k)*((-4*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2))/
         (pow(1 + pow(R0,2),0.25)*
           (4 + 4*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
        (pow(L1,2)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) \
- (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                (Qt1(i,j,k)*(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) - 
                  (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
)))/(4.*pow(rh,4))))/
         (pow(1 + pow(R0,2),0.25)*
           (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     a0.diff(d001,i,j,k)*((8*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],3)*Qr2(i,j,k))/
         (pow(1 + pow(R0,2),0.25)*
           (4 + 4*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
        (2*pow(L1,2)*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],3)*
           (L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             ((1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
                   (1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) + 
                  (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                   (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                      Qt1(i,j,k) - Qtr(i,j,k))*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
)))/(4.*pow(rh,4))))/
         (pow(1 + pow(R0,2),0.25)*
           (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) - 
     (a1.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (4.*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (ar.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (2.*L1*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  Qtr.diff(d010,i,j,k)*((4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2))/
      (pow(1 + pow(R0,2),0.25)*
        (L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (ar.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (a1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*L2*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (a2.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(L1,2)*sqrt(1 + pow(R0,2))*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (a0.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],4)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          Qtr(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(L1,2)*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (4*L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        (-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
          (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (Qt1(i,j,k)*(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
- (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (4.*pow(rh,4))))/
      (L1*pow(1 + pow(R0,2),0.25)*
        (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     a0.diff(d001,i,j,k)*((16*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],4)*Qr2(i,j,k))/
         (pow(1 + pow(R0,2),0.25)*
           (4*L1 + 4*L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
        (4*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
           (L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             ((1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
                   (1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) + 
                  (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                   (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                      Qt1(i,j,k) - Qtr(i,j,k))*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
)))/(4.*pow(rh,4))))/
         (L1*pow(1 + pow(R0,2),0.25)*
           (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) - 
     (ar.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,2)*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  a1.diff(d100,i,j,k)*(-(Qtr.diff(d001,i,j,k)*
         pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
         pow(-1 + rgrid[i],2)*
         (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
           pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
               pow(rgrid[i],3) + 
              4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                 4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
         (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
         (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (2.*L2*pow(rh,4)*pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) \
+ (pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (2.*pow(rh,4)*pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (Qt2.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (Qt1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*L2*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (Qtr.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (2.*L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(Qt1(i,j,k),2)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (2.*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  ar.diff(d010,i,j,k)*((Qtr.diff(d001,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     ((1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (Qt1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*pow(rh,4)*(1 - rgrid[i])*
        (((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (Qrr(i,j,k) - Qtt(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           pow(rh,4) + (1 - pow(1 - rgrid[i],2))*
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
              (2.*pow(rh,4)))*(-Qrr(i,j,k) + Qtt(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          ((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             pow(-1 + rgrid[i],6)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),3)*
             pow(-((1 - pow(1 - rgrid[i],2))*
                  (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                    L2*Qr2(i,j,k)*Qt2(i,j,k))) + Qtr(i,j,k),2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (16.*pow(rh,10)) + 
          ((1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             ((1 - pow(1 - rgrid[i],2))*
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
                pow(-((1 - pow(1 - rgrid[i],2))*
                     (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                       L2*Qr2(i,j,k)*Qt2(i,j,k))) + Qtr(i,j,k),2)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
+ 4*(1 - rgrid[i])*(2*L1*L2*R0*Qr1(i,j,k)*Qr2(i,j,k)*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) + 
                  pow(L1,2)*pow(Qr1(i,j,k),2)*
                   sqrt((1 + pow(R0,2))*
                     (1 + pow(R0,2) + 
                       2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                       pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2))) + 
                  pow(L2,2)*pow(Qr2(i,j,k),2)*
                   sqrt((1 + pow(R0,2))*
                     (1 + pow(R0,2) + 
                       2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                       pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2))))))/(16.*pow(rh,6))))/
      (L1*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*pow(L2,2)*pow(rh,4)*(1 - rgrid[i])*
        ((pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
             ((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt1(i,j,k),2))/pow(rh,4) - 
               pow(1 - pow(1 - rgrid[i],2),3)*
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
                   (2.*pow(rh,4)))*pow(Qt1(i,j,k),2) - 
               4*(1 - rgrid[i])*sqrt(1 + pow(R0,2))))/(4.*pow(rh,4)) \
+ (pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
             ((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt2(i,j,k),2))/pow(rh,4) - 
               pow(1 - pow(1 - rgrid[i],2),3)*
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
                   (2.*pow(rh,4)))*pow(Qt2(i,j,k),2) - 
               4*(1 - rgrid[i])*sqrt(1 + pow(R0,2))))/(4.*pow(rh,4)) \
- (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (pow(1 - pow(1 - rgrid[i],2),2)*
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
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
              (-Qrr(i,j,k) + Qtt(i,j,k))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             (pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                pow(-1 + rgrid[i],6)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),3)*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(-((1 - pow(1 - rgrid[i],2))*
                     (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                       L2*Qr2(i,j,k)*Qt2(i,j,k))) + Qtr(i,j,k),2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (16.*pow(rh,10)) - 
             (pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (pow(1 - pow(1 - rgrid[i],2),3)*
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
                   (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                   Qt1(i,j,k)*Qt2(i,j,k)*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
                  2*(1 - rgrid[i])*
                   ((1 - pow(1 - rgrid[i],2))*
                      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                      (-Qrr(i,j,k) + Qtt(i,j,k)) + 
                     2*R0*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                      (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (2.*pow(rh,4)) + 
             (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],4)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                (pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
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
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   pow(-((1 - pow(1 - rgrid[i],2))*
                        (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                        L2*Qr2(i,j,k)*Qt2(i,j,k))) + Qtr(i,j,k),2)*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) + 
                  8*(1 - rgrid[i])*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                   Qt1(i,j,k)*Qt2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) + 
                  4*pow(rh,2)*(1 - rgrid[i])*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   (2*L1*L2*R0*Qr1(i,j,k)*Qr2(i,j,k)*
                      sqrt(1 + pow(R0,2) + 
                        2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                        pow(1 - pow(1 - rgrid[i],2),2)*
                         pow(R(i,j,k),2)) + 
                     pow(L1,2)*pow(Qr1(i,j,k),2)*
                      sqrt((1 + pow(R0,2))*
                        (1 + pow(R0,2) + 
                         2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                         pow(1 - pow(1 - rgrid[i],2),2)*
                         pow(R(i,j,k),2))) + 
                     pow(L2,2)*pow(Qr2(i,j,k),2)*
                      sqrt((1 + pow(R0,2))*
                        (1 + pow(R0,2) + 
                         2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                         pow(1 - pow(1 - rgrid[i],2),2)*
                         pow(R(i,j,k),2))))))/(16.*pow(rh,8)))))/
      (L1*pow(1 + pow(R0,2),0.25)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     ((1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(Qt1(i,j,k),2)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        (L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  ar.diff(d001,i,j,k)*(-(((1 - pow(1 - rgrid[i],2))*
          pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          pow(Qt1(i,j,k),2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
        (L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  B.diff(d010,i,j,k)*((-16*a2.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*L2*sqrt(1 + pow(R0,2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (16*a1.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (16*a0.diff(d001,i,j,k)*pow(rh,4)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (ar(i,j,k)*(1 - rgrid[i])*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
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
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
              (2.*pow(rh,4)))*
           (pow(Qt1(i,j,k),2) - 
             pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              pow(Qt2(i,j,k),2)) + 
          4*(1 - rgrid[i])*(pow(L1 + 
                L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              ((pow(1 - pow(1 - rgrid[i],2),2)*
                   pow(-1 + rgrid[i],2)*
                   (pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                   pow(Qt2(i,j,k),2))/(4.*pow(rh,4)) - 
                sqrt(1 + pow(R0,2))) + 
             pow(L1,2)*(-(pow(1 - pow(1 - rgrid[i],2),2)*
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
                sqrt(1 + pow(R0,2)))))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,3)*pow(1 + pow(R0,2),0.25)*
        pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  a2.diff(d100,i,j,k)*(-(Qtr.diff(d010,i,j,k)*
         pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
         pow(-1 + rgrid[i],2)*
         (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
           pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
               pow(rgrid[i],3) + 
              4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                 4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
         (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
         (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (2.*L1*sqrt(1 + pow(R0,2))*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(Qt1(i,j,k),2)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (2.*sqrt(1 + pow(R0,2))*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (Qtr.diff(d001,i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (2.*L2*sqrt(1 + pow(R0,2))*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (2.*sqrt(1 + pow(R0,2))*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (Qt2.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*R0*Qr2(i,j,k) + L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
           R(i,j,k) + L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*sqrt(1 + 
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (2.*L1*sqrt(1 + pow(R0,2))*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (Qt1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*R0*Qr2(i,j,k) + L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
           R(i,j,k) + L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*sqrt(1 + 
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (2.*L2*sqrt(1 + pow(R0,2))*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  Qt1.diff(d001,i,j,k)*((-8*a0(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k))/
      (pow(1 + pow(R0,2),0.25)*
        (2 + 2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     a1.diff(d001,i,j,k)*((pow(1 - pow(1 - rgrid[i],2),5)*
           pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
               2)*pow(Qr1(i,j,k),2)*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             pow(L2,2)*pow(Qr2(i,j,k),2)*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              Qr1(i,j,k)*Qr2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (8*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           ((pow(L1,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr1(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr2(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
                (L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                   pow(-1 + rgrid[i],2)*
                   (pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   Qr1(i,j,k)*Qr2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                 (2.*pow(rh,2)))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L2,2)*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     a0.diff(d001,i,j,k)*((16*pow(rh,2)*
           pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
           (-(pow(L1,2)*pow(L2,2)*
                 pow(1 - pow(1 - rgrid[i],2),4)*
                 pow(-1 + rgrid[i],4)*
                 pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   pow(rh,2)*
                    (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                 pow(Qr2(i,j,k),2)*pow(Qt1(i,j,k),2)*
                 (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                 (1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
))/(16.*pow(rh,6)) - pow(L1 + 
                L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
                (L1*pow(1 - pow(1 - rgrid[i],2),3)*
                   pow(-1 + rgrid[i],4)*
                   pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                   Qr1(i,j,k)*Qt1(i,j,k)*
                   (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                      Qt1(i,j,k) - Qtr(i,j,k))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
                 (16.*pow(rh,6)))*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             (pow(L1,2)*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qr2(i,j,k)*(2*L2*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   Qr2(i,j,k) + 
                  ((1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                     (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                       pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                     Qt1(i,j,k)*
                     (2*L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                       Qt1(i,j,k) - Qtr(i,j,k))*
                     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                   (4.*pow(rh,4)))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2))))/
         (pow(L1,2)*pow(L2,2)*pow(1 + pow(R0,2),0.25)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (32*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],4)*
           (-(pow(L1,2)*pow(L2,2)*
                 pow(1 - pow(1 - rgrid[i],2),4)*
                 pow(-1 + rgrid[i],4)*
                 pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   pow(rh,2)*
                    (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                 (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                 pow(Qr2(i,j,k),2)*pow(Qt1(i,j,k),2)*
                 (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                 (1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
)/(16.*pow(rh,6)) - pow(L1 + 
                L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
                (L1*pow(1 - pow(1 - rgrid[i],2),3)*
                   pow(-1 + rgrid[i],4)*
                   pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                   Qr1(i,j,k)*Qt1(i,j,k)*
                   (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                      Qt1(i,j,k) - Qtr(i,j,k))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
                 (16.*pow(rh,6)))*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (pow(L2,2)*pow(rh,2)*
                   pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
                   pow(Qr2(i,j,k),2) + 
                  (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                   pow(Qt1(i,j,k),2)*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
                  (L2*(1 - pow(1 - rgrid[i],2))*
                     pow(-1 + rgrid[i],2)*
                     (pow(H,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                     Qr2(i,j,k)*Qt1(i,j,k)*
                     (2*L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                        Qt1(i,j,k) - Qtr(i,j,k))*
                     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                   (4.*pow(rh,2)))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,4))))/
         (pow(L1,2)*pow(L2,2)*pow(1 + pow(R0,2),0.25)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     (4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*
        (L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (-(L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                  Qr1(i,j,k)*
                  ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                     Qt2(i,j,k)*
                     (1 + pow(R0,2) + 
                       2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                       pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2)) - 
                    Qt1(i,j,k)*
                     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                     sqrt(1 + pow(R0,2) + 
                       2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                       pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2)))) + 
               L2*Qr2(i,j,k)*(Qt1(i,j,k)*
                   (1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) - 
                  (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
))))/(4.*pow(rh,4))))/
      (pow(1 + pow(R0,2),0.25)*
        (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (ar.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*R0*Qr2(i,j,k) + L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
           R(i,j,k) + L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*sqrt(1 + 
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (pow(L2,2)*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  Qt2.diff(d010,i,j,k)*((-8*L2*a0(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k))/
      (pow(1 + pow(R0,2),0.25)*
        (2*L1 + 2*L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (2*L2*pow(rh,4)*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        (((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           pow(rh,4) - (1 - pow(1 - rgrid[i],2))*
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
              (2.*pow(rh,4))))*Qr2(i,j,k)*
        (4*a0(i,j,k)*pow(1 - rgrid[i],2) + 
          (ar(i,j,k)*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             ((1 - pow(1 - rgrid[i],2))*
                (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
               Qtr(i,j,k)))/(4.*pow(rh,6))))/
      (L1*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (ar.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(L1,2)*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (2*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        (-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
             pow(rh,4)) + (1 - pow(1 - rgrid[i],2))*
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
        (4*L2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k) + 
          (ar(i,j,k)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             ((L2*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                  Qr2(i,j,k)*
                  (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                     Qt2(i,j,k) - Qtr(i,j,k)))/(4.*pow(rh,2)) - 
               Qt1(i,j,k)*(-(L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(-1 + rgrid[i],2)*
                      (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                      Qr1(i,j,k)*Qr2(i,j,k))/(4.*pow(rh,2)) + 
                  (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
               (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           pow(rh,4)))/
      (L1*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     a1.diff(d001,i,j,k)*(-((pow(1 - pow(1 - rgrid[i],2),5)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (pow(L1,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(Qr1(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
+ pow(L2,2)*pow(Qr2(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
+ 2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
                Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (L1*L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
) + (8*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           ((pow(L1,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr1(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr2(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
                (L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                   pow(-1 + rgrid[i],2)*
                   (pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   Qr1(i,j,k)*Qr2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                 (2.*pow(rh,2)))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     a2.diff(d010,i,j,k)*((pow(1 - pow(1 - rgrid[i],2),5)*
           pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
               2)*pow(Qr1(i,j,k),2)*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             pow(L2,2)*pow(Qr2(i,j,k),2)*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              Qr1(i,j,k)*Qr2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (sqrt(1 + pow(R0,2))*pow(rh,4)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (8*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           ((pow(L1,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr1(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr2(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
                (L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                   pow(-1 + rgrid[i],2)*
                   (pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   Qr1(i,j,k)*Qr2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                 (2.*pow(rh,2)))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L1,2)*sqrt(1 + pow(R0,2))*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     a0.diff(d010,i,j,k)*(-((L2*pow(1 - pow(1 - rgrid[i],2),5)*
             pow(1 - rgrid[i],4)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (pow(L1,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                (1 - pow(1 - rgrid[i],2))*pow(Qr1(i,j,k),2)*
                Qt2(i,j,k)*(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
+ L2*Qr2(i,j,k)*(L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                   Qt2(i,j,k) - Qtr(i,j,k))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
+ L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
                (2*L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                   Qt2(i,j,k) - Qtr(i,j,k))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (pow(1 + pow(R0,2),0.25)*pow(rh,4)*
             (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
) + (8*L2*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
           Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           ((pow(L1,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr1(i,j,k),2)*Qt2(i,j,k)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (L2*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
                (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                   Qt2(i,j,k) - Qtr(i,j,k))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k) + 
                (L1*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                   (pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   Qr1(i,j,k)*
                   (2*L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                      Qt2(i,j,k) - Qtr(i,j,k))*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                 (4.*pow(rh,2)))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L1,2)*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     a0.diff(d001,i,j,k)*((16*pow(rh,4)*
           pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
           ((pow(L1,2)*pow(L2,2)*
                pow(1 - pow(1 - rgrid[i],2),4)*
                pow(-1 + rgrid[i],4)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                pow(Qr2(i,j,k),2)*pow(Qt1(i,j,k),2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
)/(16.*pow(rh,6)) + pow(L1 + 
                L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              (-1 - (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
                (L1*pow(1 - pow(1 - rgrid[i],2),3)*
                   pow(-1 + rgrid[i],4)*
                   pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                   Qr1(i,j,k)*Qt1(i,j,k)*
                   (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                      Qt1(i,j,k) - Qtr(i,j,k))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
                 (16.*pow(rh,6)))*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (pow(L1,2)*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qr2(i,j,k)*(-2*L2*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   Qr2(i,j,k) + 
                  ((1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                     (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                       pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                     Qt1(i,j,k)*
                     (2*L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                       Qt1(i,j,k) - Qtr(i,j,k))*
                     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                   (4.*pow(rh,4)))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2))))/
         (pow(L1,3)*L2*pow(1 + pow(R0,2),0.25)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (32*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],4)*
           ((pow(L1,2)*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),4)*
                pow(-1 + rgrid[i],4)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr2(i,j,k),2)*pow(Qt1(i,j,k),2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (16.*pow(rh,6)) + 
             pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
              (-1 - (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
                (L1*pow(1 - pow(1 - rgrid[i],2),3)*
                   pow(-1 + rgrid[i],4)*
                   pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                   Qr1(i,j,k)*Qt1(i,j,k)*
                   (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                      Qt1(i,j,k) - Qtr(i,j,k))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
                 (16.*pow(rh,6)))*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (-(pow(L2,2)*pow(rh,2)*
                     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
                     pow(Qr2(i,j,k),2)) + 
                  (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                   pow(Qt1(i,j,k),2)*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
                  (L2*(1 - pow(1 - rgrid[i],2))*
                     pow(-1 + rgrid[i],2)*
                     (pow(H,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                     Qr2(i,j,k)*Qt1(i,j,k)*
                     (2*L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                        Qt1(i,j,k) - Qtr(i,j,k))*
                     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                   (4.*pow(rh,2)))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,4))))/
         (L2*pow(1 + pow(R0,2),0.25)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           pow(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     (4*L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*
        (L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
                ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                   (1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) - 
                  Qt1(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2))) - 
               L2*Qr2(i,j,k)*(Qt1(i,j,k)*
                   (1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) - 
                  (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
))))/(4.*pow(rh,4))))/
      (L1*pow(1 + pow(R0,2),0.25)*
        (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (ar.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*R0*Qr2(i,j,k) + L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
           R(i,j,k) + L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*sqrt(1 + 
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  a2.diff(d010,i,j,k)*((pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*sqrt(1 + pow(R0,2))*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(Qt1(i,j,k),2)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (sqrt(1 + pow(R0,2))*pow(rh,4)*
        (L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     Qt1.diff(d001,i,j,k)*(-((pow(1 - pow(1 - rgrid[i],2),5)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (pow(L1,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(Qr1(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
+ pow(L2,2)*pow(Qr2(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
+ 2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
                Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (L1*L2*sqrt(1 + pow(R0,2))*pow(rh,2)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
) + (8*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           ((pow(L1,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr1(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr2(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
                (L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                   pow(-1 + rgrid[i],2)*
                   (pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   Qr1(i,j,k)*Qr2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                 (2.*pow(rh,2)))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*sqrt(1 + pow(R0,2))*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     (8*Qrr.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*L2*sqrt(1 + pow(R0,2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (32*Q.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*L2*sqrt(1 + pow(R0,2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)) + 
     (8*Qtt.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*L2*sqrt(1 + pow(R0,2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
     (Qtr.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*R0*Qr2(i,j,k) + L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
           R(i,j,k) + L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*sqrt(1 + 
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (L1*L2*sqrt(1 + pow(R0,2))*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  H*((-8*Qrr.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (32*Q.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)) - 
     (8*Qtt.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
     (pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (Qt1.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (Qtr.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(Qt1(i,j,k),2)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     Qt2.diff(d010,i,j,k)*((pow(1 - pow(1 - rgrid[i],2),5)*
           pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
               2)*pow(Qr1(i,j,k),2)*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             pow(L2,2)*pow(Qr2(i,j,k),2)*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              Qr1(i,j,k)*Qr2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (8*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           ((pow(L1,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr1(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr2(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
                (L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                   pow(-1 + rgrid[i],2)*
                   (pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   Qr1(i,j,k)*Qr2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                 (2.*pow(rh,2)))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     Qt1.diff(d001,i,j,k)*(-((pow(1 - pow(1 - rgrid[i],2),5)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (pow(L1,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(Qr1(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
+ pow(L2,2)*pow(Qr2(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
+ 2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
                Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
             (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
) + (8*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           ((pow(L1,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr1(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr2(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
                (L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                   pow(-1 + rgrid[i],2)*
                   (pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   Qr1(i,j,k)*Qr2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                 (2.*pow(rh,2)))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     (8*Qrr.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (32*Q.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)) + 
     (8*Qtt.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
     (Qt2.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*R0*Qr2(i,j,k) + L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
           R(i,j,k) + L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*sqrt(1 + 
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (2.*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (Qtr.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*R0*Qr2(i,j,k) + L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
           R(i,j,k) + L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*sqrt(1 + 
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  Qt2.diff(d100,i,j,k)*((a1.diff(d100,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (4.*pow(rh,4)*pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (ar.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (2.*L1*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (a0.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],3)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          Qtr(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*L1*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (a0.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],3)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
             Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (a2.diff(d100,i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (4.*sqrt(1 + pow(R0,2))*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (ar.diff(d001,i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (2.*L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (a1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*R0*Qr2(i,j,k) + L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
           R(i,j,k) + L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*sqrt(1 + 
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (2.*L2*pow(rh,4)*pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) \
+ (a2.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*R0*Qr2(i,j,k) + L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
           R(i,j,k) + L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*sqrt(1 + 
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (2.*L1*sqrt(1 + pow(R0,2))*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0*Qt1(i,j,k) + (1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
           R(i,j,k) - (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qt2(i,j,k)*sqrt(1 + 
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (2.*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (a0.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0*Qt1(i,j,k) + (1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*R(i,j,k) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (4.*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  Qtr.diff(d001,i,j,k)*(-((a0.diff(d001,i,j,k)*
          pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
            (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
               Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (pow(L2,2)*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) - 
     (ar.diff(d001,i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (a1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*R0*Qr2(i,j,k) + L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
           R(i,j,k) + L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*sqrt(1 + 
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (pow(L2,2)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (-(R0*Qt1(i,j,k)) - (1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
           R(i,j,k) + (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qt2(i,j,k)*sqrt(1 + 
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  a0.diff(d100,i,j,k)*((-2*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],3)*Q(i,j,k)*Qt1(i,j,k))/
      (pow(1 + pow(R0,2),0.25)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (2*Q.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],3)*Qr1(i,j,k)*Qt1(i,j,k))/
      (pow(1 + pow(R0,2),0.25)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (2*Q.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],3)*Qr2(i,j,k)*Qt1(i,j,k))/
      (pow(1 + pow(R0,2),0.25)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],3)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (2*pow(L1,2)*L2*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],3)*
        Qt1(i,j,k)*(-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
          (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (Qt1(i,j,k)*(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
- (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (4.*pow(rh,4))))/
      (pow(1 + pow(R0,2),0.25)*
        (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     Qtr.diff(d010,i,j,k)*((-2*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],3))/
         (pow(1 + pow(R0,2),0.25)*
           (L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
        (2*L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*
           (-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) \
+ (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                (Qt1(i,j,k)*(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) - 
                  (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
)))/(4.*pow(rh,4))))/
         (L1*pow(1 + pow(R0,2),0.25)*
           (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     Qt2.diff(d010,i,j,k)*((8*L2*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],3)*Qr2(i,j,k))/
         (pow(1 + pow(R0,2),0.25)*
           (4*L1 + 4*L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
        (2*L2*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],3)*
           (L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                   Qr1(i,j,k)*
                   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                      Qt2(i,j,k)*
                      (1 + pow(R0,2) + 
                        2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                        pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2)) - 
                     Qt1(i,j,k)*
                      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                      sqrt(1 + pow(R0,2) + 
                        2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                        pow(1 - pow(1 - rgrid[i],2),2)*
                         pow(R(i,j,k),2))) - 
                  L2*Qr2(i,j,k)*
                   (Qt1(i,j,k)*
                      (1 + pow(R0,2) + 
                        2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                        pow(1 - pow(1 - rgrid[i],2),2)*
                         pow(R(i,j,k),2)) - 
                     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                      Qt2(i,j,k)*
                      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                      sqrt(1 + pow(R0,2) + 
                        2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                        pow(1 - pow(1 - rgrid[i],2),2)*
                         pow(R(i,j,k),2)))))/(4.*pow(rh,4))))/
         (L1*pow(1 + pow(R0,2),0.25)*
           (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     Qt1.diff(d001,i,j,k)*((8*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],3)*Qr2(i,j,k))/
         (pow(1 + pow(R0,2),0.25)*
           (4 + 4*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
        (2*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],3)*
           (L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                (-(L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                     Qr1(i,j,k)*
                     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                        Qt2(i,j,k)*
                        (1 + pow(R0,2) + 
                        2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                        pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2)) - 
                       Qt1(i,j,k)*
                        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                        sqrt(1 + pow(R0,2) + 
                         2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                         pow(1 - pow(1 - rgrid[i],2),2)*
                         pow(R(i,j,k),2)))) + 
                  L2*Qr2(i,j,k)*
                   (Qt1(i,j,k)*
                      (1 + pow(R0,2) + 
                        2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                        pow(1 - pow(1 - rgrid[i],2),2)*
                         pow(R(i,j,k),2)) - 
                     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                      Qt2(i,j,k)*
                      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                      sqrt(1 + pow(R0,2) + 
                        2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                        pow(1 - pow(1 - rgrid[i],2),2)*
                         pow(R(i,j,k),2)))))/(4.*pow(rh,4))))/
         (pow(1 + pow(R0,2),0.25)*
           (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) - 
     (Qtr.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],3)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (-(R0*Qt1(i,j,k)) - (1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
           R(i,j,k) + (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qt2(i,j,k)*sqrt(1 + 
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (2.*L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  Q.diff(d001,i,j,k)*((-4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*Qt1(i,j,k))/
      (pow(1 + pow(R0,2),0.25)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (32*a1.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)) + 
     a0.diff(d001,i,j,k)*((4*pow(1 - pow(1 - rgrid[i],2),5)*
           pow(1 - rgrid[i],4)*pow(Qr2(i,j,k),2)*Qt1(i,j,k))/
         (pow(1 + pow(R0,2),0.25)*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
        (8*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
           Qt1(i,j,k)*(pow(1 - pow(1 - rgrid[i],2),2)*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
              pow(Qr2(i,j,k),2) + 
             (4*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),
                   2)))/
              (pow(L2,2)*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))))/
         (pow(1 + pow(R0,2),0.25)*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)))) + 
  a0.diff(d010,i,j,k)*((4*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],4)*Q(i,j,k)*Qr1(i,j,k)*Qt1(i,j,k))/
      (pow(1 + pow(R0,2),0.25)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (16*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
        Qr1(i,j,k)*Qt1(i,j,k))/
      (pow(1 + pow(R0,2),0.25)*
        (4 + 4*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (16*R.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (8*Qrr.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
        (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (8*Qtt.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
     (Qtr.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],4)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          Qtr(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          Qtr(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(Qt1(i,j,k),2)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)*
           (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          Qtr(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        (L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     Qt1.diff(d001,i,j,k)*((pow(1 - pow(1 - rgrid[i],2),5)*
           pow(1 - rgrid[i],4)*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
               2)*(1 - pow(1 - rgrid[i],2))*pow(Qr1(i,j,k),2)*
              Qt2(i,j,k)*(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             L2*Qr2(i,j,k)*(L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                 Qt2(i,j,k) - Qtr(i,j,k))*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (2*L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                Qtr(i,j,k))*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
           (L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (8*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
           Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           ((pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
                 2)*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr1(i,j,k),2)*Qt2(i,j,k)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (L2*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
                (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                  Qtr(i,j,k))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2)) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k) + 
                (L1*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                   (pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                         (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                         4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                   (2*L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                      Qt2(i,j,k) - Qtr(i,j,k))*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                 (4.*pow(rh,2)))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L2*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) - 
     (16*Qt2.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     Q.diff(d001,i,j,k)*((4*pow(1 - pow(1 - rgrid[i],2),5)*
           pow(1 - rgrid[i],4)*Qr1(i,j,k)*Qr2(i,j,k)*Qt1(i,j,k))/
         (pow(1 + pow(R0,2),0.25)*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
        (8*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
           (L1*pow(1 - pow(1 - rgrid[i],2),2)*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
              Qr2(i,j,k)*Qt1(i,j,k) - 
             (4*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
                sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),
                   2)))/
              (L2*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))))/
         (L1*pow(1 + pow(R0,2),0.25)*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3))))
;

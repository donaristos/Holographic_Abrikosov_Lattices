
elem[149]+=
-4*Qtr.diff(d110,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]) + 
  8*Qtr.diff(d020,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],2)*Qr1(i,j,k) + 
  8*Qtr.diff(d011,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],2)*Qr2(i,j,k) + 
  Qt1.diff(d010,i,j,k)*(4*L1*(1 - pow(1 - rgrid[i],2))*
      (-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*Qr1(i,j,k) \
+ (8*L1*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*
        ((2*pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           pow(rh,4) + (1 - pow(1 - rgrid[i],2))*
           ((pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
              pow(rh,4) + 2*(1 - rgrid[i])*
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
                 (2.*pow(rh,4)))))*Qr1(i,j,k))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))) - 
     (8*L1*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k)*Qrr(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) \
- (2*L1*(-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
     Qt1(i,j,k))/(1 - pow(1 - rgrid[i],2)) - 
  (2*L1*pow(rh,4)*(-2*(1 - pow(1 - rgrid[i],2)) + 
       4*pow(1 - rgrid[i],2))*
     ((pow(-1 + rgrid[i],2)*
          (3*pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],2) + 
            3*pow(H,2)*pow(-2 + rgrid[i],2)*pow(rgrid[i],3) + 
            pow(rh,2)*(3*pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],2) + 
               3*pow(mu,2)*pow(-2 + rgrid[i],2)*pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (2 + 6*rgrid[i] - 12*pow(rgrid[i],2) + 
                  4*pow(rgrid[i],3)))))/(4.*pow(rh,4)) + 
       ((-1 + rgrid[i])*(pow(H,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
        (2.*pow(rh,4)))*Qt1(i,j,k))/
   ((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) - 
  (4*L1*pow(1 - rgrid[i],2)*Q(i,j,k)*Qt1(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
  (4*L1*pow(1 - rgrid[i],2)*Qrr(i,j,k)*Qt1(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
  Q.diff(d100,i,j,k)*(-((Qt1.diff(d100,i,j,k)*L1*
          pow(1 - pow(1 - rgrid[i],2),2))/
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (2*Qtr.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i]))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
     (2*L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qt1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  (pow(L1,2)*pow(1 - pow(1 - rgrid[i],2),2)*
     (-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
     pow(-1 + rgrid[i],4)*pow(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*Qr1(i,j,k)*
     pow(Qt1(i,j,k),2)*Qtr(i,j,k))/(8.*pow(rh,6)) + 
  (L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
     (-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
     pow(-1 + rgrid[i],4)*pow(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*Qr2(i,j,k)*
     Qt1(i,j,k)*Qt2(i,j,k)*Qtr(i,j,k))/(8.*pow(rh,6)) + 
  (L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/pow(rh,4) \
- (1 - pow(1 - rgrid[i],2))*
        ((pow(-1 + rgrid[i],2)*
             (3*pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],2) + 
               3*pow(H,2)*pow(-2 + rgrid[i],2)*pow(rgrid[i],3) + 
               pow(rh,2)*(3*pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],2) + 
                  3*pow(mu,2)*pow(-2 + rgrid[i],2)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (2 + 6*rgrid[i] - 12*pow(rgrid[i],2) + 
                     4*pow(rgrid[i],3)))))/(4.*pow(rh,4)) + 
          ((-1 + rgrid[i])*(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           (2.*pow(rh,4))))*Qr2(i,j,k)*Qt1(i,j,k)*Qt2(i,j,k)*Qtr(i,j,k))/
   (2.*pow(rh,2)) - (L1*L2*(-2*(1 - pow(1 - rgrid[i],2)) + 
       4*pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qr2(i,j,k)*
     (4*R0*(1 - rgrid[i]) + pow(1 - pow(1 - rgrid[i],2),3)*
        ((pow(-1 + rgrid[i],2)*
             (3*pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],2) + 
               3*pow(H,2)*pow(-2 + rgrid[i],2)*pow(rgrid[i],3) + 
               pow(rh,2)*(3*pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],2) + 
                  3*pow(mu,2)*pow(-2 + rgrid[i],2)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (2 + 6*rgrid[i] - 12*pow(rgrid[i],2) + 
                     4*pow(rgrid[i],3)))))/(4.*pow(rh,4)) + 
          ((-1 + rgrid[i])*(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           (2.*pow(rh,4)))*Qt1(i,j,k)*Qt2(i,j,k))*Qtr(i,j,k))/
   (8.*pow(rh,2)*(1 - rgrid[i])) + 
  (L1*(1 - pow(1 - rgrid[i],2))*
     (-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
     pow(-1 + rgrid[i],4)*pow(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*Qt1(i,j,k)*
     ((1 - pow(1 - rgrid[i],2))*
        (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)\
)*Qtr(i,j,k))/(8.*pow(rh,6)) + 
  Qr2.diff(d100,i,j,k)*(4*Qt2.diff(d010,i,j,k)*L2*
      pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i]) + 
     (L1*L2*pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
        pow(-1 + rgrid[i],4)*
        pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*Qt1(i,j,k)*
        Qt2(i,j,k)*Qtr(i,j,k))/(4.*pow(rh,6)) - 
     (L1*L2*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (4*R0*(1 - rgrid[i]) + 
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
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
              (2.*pow(rh,4)))*Qt1(i,j,k)*Qt2(i,j,k))*Qtr(i,j,k))/
      (4.*pow(rh,2))) + (4*L1*pow(rh,4)*Qt1(i,j,k)*
     (-2*(1 - rgrid[i])*((pow(-1 + rgrid[i],2)*
             (3*pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],2) + 
               3*pow(H,2)*pow(-2 + rgrid[i],2)*pow(rgrid[i],3) + 
               pow(rh,2)*(3*pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],2) + 
                  3*pow(mu,2)*pow(-2 + rgrid[i],2)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (2 + 6*rgrid[i] - 12*pow(rgrid[i],2) + 
                     4*pow(rgrid[i],3)))))/(4.*pow(rh,4)) + 
          ((-1 + rgrid[i])*(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           (2.*pow(rh,4)))*(1 + 
          (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             (2*L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) + 
               L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
               Qtr(i,j,k))*Qtr(i,j,k))/(16.*pow(rh,6))) + 
       (pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (-2 + ((1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
               pow(-1 + rgrid[i],4)*
               pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
               (2*L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                  Qt1(i,j,k) + 
                 L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                 Qtr(i,j,k))*Qtr(i,j,k))/(4.*pow(rh,6))))/
        (2.*pow(rh,4))))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
  (8*a0.diff(d010,i,j,k)*a2.diff(d100,i,j,k)*L2*pow(rh,2)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],3)*Qr2(i,j,k))/
   (pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
  (16*L1*pow(rh,2)*pow(a0(i,j,k),2)*pow(1 - pow(1 - rgrid[i],2),3)*
     pow(1 - rgrid[i],2)*Qt1(i,j,k))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
  (4*pow(a0.diff(d100,i,j,k),2)*L1*pow(rh,2)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*Qt1(i,j,k))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
  (4*Qtt.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
     pow(1 - rgrid[i],2)*Qr2(i,j,k)*Qt1(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
  (16*a0.diff(d001,i,j,k)*L1*pow(rh,2)*a0(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*Qr2(i,j,k)*
     Qt1(i,j,k))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
  (4*L1*pow(1 - rgrid[i],2)*Qt1(i,j,k)*Qtt(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
  ar.diff(d010,i,j,k)*((-16*pow(rh,2)*a0(i,j,k)*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
     (16*a0.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*Qr2(i,j,k)\
)/(pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))) + 
  Qtt.diff(d100,i,j,k)*((Qt1.diff(d100,i,j,k)*L1*
        pow(1 - pow(1 - rgrid[i],2),2))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
     (2*Qtr.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i]))/(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
     (2*L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qt1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
  (4*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(1 - rgrid[i],2)*
     Qt2(i,j,k)*R(i,j,k))/
   sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) - 
  (4*L1*B(i,j,k)*pow(1 - rgrid[i],2)*Qt1(i,j,k)*
     (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
  R.diff(d100,i,j,k)*((4*Qt2.diff(d010,i,j,k)*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*Qr1(i,j,k))/
      sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
     (4*Qt1.diff(d001,i,j,k)*pow(L1,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*Qr1(i,j,k))/
      (L2*sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
     (Qt2.diff(d100,i,j,k)*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
     (2*Qtr.diff(d001,i,j,k)*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i]))/
      (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qt2(i,j,k))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  R.diff(d001,i,j,k)*((-4*Qtr.diff(d001,i,j,k)*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k))/
      (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (4*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k)*Qt2(i,j,k))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
     (16*Qt1.diff(d001,i,j,k)*L1*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        ((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qr2(i,j,k))/(2.*pow(rh,2)) - 
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) - 
  (pow(L1,2)*(-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qr1(i,j,k)*
     Qtr(i,j,k)*(pow(1 - pow(1 - rgrid[i],2),3)*
        ((pow(-1 + rgrid[i],2)*
             (3*pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],2) + 
               3*pow(H,2)*pow(-2 + rgrid[i],2)*pow(rgrid[i],3) + 
               pow(rh,2)*(3*pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],2) + 
                  3*pow(mu,2)*pow(-2 + rgrid[i],2)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (2 + 6*rgrid[i] - 12*pow(rgrid[i],2) + 
                     4*pow(rgrid[i],3)))))/(4.*pow(rh,4)) + 
          ((-1 + rgrid[i])*(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           (2.*pow(rh,4)))*pow(Qt1(i,j,k),2) + 
       4*(1 - rgrid[i])*sqrt(1 + pow(R0,2))))/
   (8.*pow(rh,2)*(1 - rgrid[i])) + 
  Qr1.diff(d100,i,j,k)*(4*Qt1.diff(d010,i,j,k)*L1*
      pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i]) + 
     (pow(L1,2)*pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
        pow(-1 + rgrid[i],4)*
        pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
        pow(Qt1(i,j,k),2)*Qtr(i,j,k))/(4.*pow(rh,6)) - 
     (pow(L1,2)*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qtr(i,j,k)*
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
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
              (2.*pow(rh,4)))*pow(Qt1(i,j,k),2) + 
          4*(1 - rgrid[i])*sqrt(1 + pow(R0,2))))/(4.*pow(rh,2))) - 
  (L1*(-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qtr(i,j,k)*
     (pow(1 - pow(1 - rgrid[i],2),3)*
        ((pow(-1 + rgrid[i],2)*
             (3*pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],2) + 
               3*pow(H,2)*pow(-2 + rgrid[i],2)*pow(rgrid[i],3) + 
               pow(rh,2)*(3*pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],2) + 
                  3*pow(mu,2)*pow(-2 + rgrid[i],2)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (2 + 6*rgrid[i] - 12*pow(rgrid[i],2) + 
                     4*pow(rgrid[i],3)))))/(4.*pow(rh,4)) + 
          ((-1 + rgrid[i])*(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           (2.*pow(rh,4)))*Qt1(i,j,k)*
        (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
       pow(1 - pow(1 - rgrid[i],2),2)*
        ((pow(-1 + rgrid[i],2)*
             (3*pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],2) + 
               3*pow(H,2)*pow(-2 + rgrid[i],2)*pow(rgrid[i],3) + 
               pow(rh,2)*(3*pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],2) + 
                  3*pow(mu,2)*pow(-2 + rgrid[i],2)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (2 + 6*rgrid[i] - 12*pow(rgrid[i],2) + 
                     4*pow(rgrid[i],3)))))/(4.*pow(rh,4)) + 
          ((-1 + rgrid[i])*(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           (2.*pow(rh,4)))*Qt1(i,j,k)*Qtr(i,j,k) + 
       4*(1 - rgrid[i])*(L2*R0*Qr2(i,j,k) + 
          L1*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))))/
   (8.*pow(rh,2)*(1 - rgrid[i])) + 
  (8*L1*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qrr(i,j,k)*
     ((L1*pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],4)*
          pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
          (-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                 (pow(H,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
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
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
                (2.*pow(rh,4))))*Qr1(i,j,k)*pow(Qt1(i,j,k),2)*
          Qtr(i,j,k))/(16.*pow(rh,6)) + 
       Qt1(i,j,k)*(((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             pow(-1 + rgrid[i],6)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),3)*
             Qtr(i,j,k)*(-(L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                  Qt2(i,j,k)) + Qtr(i,j,k)))/(16.*pow(rh,10)) + 
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
           (1 + (pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],4)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                  Qtr(i,j,k))*Qtr(i,j,k))/(16.*pow(rh,6)))) + 
       ((1 - rgrid[i])*pow(-1 + rgrid[i],4)*
          pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*Qtr(i,j,k)*
          (L2*R0*Qr2(i,j,k) + L1*Qr1(i,j,k)*sqrt(1 + pow(R0,2))))/
        (4.*pow(rh,6))))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
  Qrr.diff(d100,i,j,k)*((Qt1.diff(d100,i,j,k)*L1*
        pow(1 - pow(1 - rgrid[i],2),2))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     (2*Qtr.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i]))/(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
     (4*Qt1.diff(d010,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
     (2*Qt1.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr2(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
     (2*Qt2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr2(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     (2*L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qt1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     (4*L1*pow(rh,4)*pow(1 - pow(1 - rgrid[i],2),2)*
        ((L1*pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             (-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                    (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
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
                     (pow(H,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
                   (2.*pow(rh,4))))*Qr1(i,j,k)*pow(Qt1(i,j,k),2)*
             Qtr(i,j,k))/(16.*pow(rh,6)) + 
          Qt1(i,j,k)*(((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                pow(-1 + rgrid[i],6)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),3)*
                Qtr(i,j,k)*(-(L2*(1 - pow(1 - rgrid[i],2))*
                     Qr2(i,j,k)*Qt2(i,j,k)) + Qtr(i,j,k)))/
              (16.*pow(rh,10)) + 
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
              (1 + (pow(1 - pow(1 - rgrid[i],2),2)*
                   pow(-1 + rgrid[i],4)*
                   pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                   (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                      Qt2(i,j,k) - Qtr(i,j,k))*Qtr(i,j,k))/
                 (16.*pow(rh,6)))) + 
          ((1 - rgrid[i])*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             Qtr(i,j,k)*(L2*R0*Qr2(i,j,k) + 
               L1*Qr1(i,j,k)*sqrt(1 + pow(R0,2))))/(4.*pow(rh,6))))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))) + 
  (32*pow(a0.diff(d010,i,j,k),2)*pow(rh,2)*
     pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],6)*Qr1(i,j,k)*
     Qtr(i,j,k)*sqrt(1 + pow(R0,2) + 
       2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  Qr2.diff(d001,i,j,k)*(-4*L1*(1 - pow(1 - rgrid[i],2))*
      pow(1 - rgrid[i],2)*Qt1(i,j,k)*
      pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2) + 
     (8*Qt1.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*
        (L2*Qr2(i,j,k)*pow(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2),1.5) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + pow(R0,3) + 
             (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
             pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3))))/
      (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*Qtr.diff(d001,i,j,k)*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/L2 + 
     4*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
      (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qt2(i,j,k)*
      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
      sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  Qtt.diff(d010,i,j,k)*((4*Qtr.diff(d001,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
     (8*L1*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k)*Qt1(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) \
- (4*L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k)*Qt2(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) \
+ (8*Qtr.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k)*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      ((1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  a1.diff(d100,i,j,k)*((8*L1*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
        a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i]))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
     (8*a0.diff(d001,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],3)*
        Qr2(i,j,k))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
     (16*a0.diff(d010,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],3)*
        Qr1(i,j,k)*sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  a0.diff(d010,i,j,k)*((-16*ar.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
        Qr2(i,j,k))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
     (32*a0.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],6)*
        Qr2(i,j,k)*Qtr(i,j,k))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
     (32*ar.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
        Qr1(i,j,k)*sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (16*pow(rh,2)*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],4)*((1 - pow(1 - rgrid[i],2))*
           (2*L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) + 
          Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  a0.diff(d100,i,j,k)*((-4*a1.diff(d100,i,j,k)*L1*
        pow(1 + pow(R0,2),0.25)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
     (8*ar.diff(d010,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],3))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
     (16*L1*pow(rh,2)*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],3)*Qt1(i,j,k))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
     (8*a0.diff(d001,i,j,k)*L1*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],5)*
        Qr2(i,j,k)*Qt1(i,j,k))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
     (8*a0.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],5)*
        ((1 - pow(1 - rgrid[i],2))*
           (2*L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) + 
          Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  Qrr.diff(d010,i,j,k)*((16*Qt2.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (L1*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (16*Qt1.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  Q.diff(d010,i,j,k)*((-8*Qtr.diff(d010,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
     (4*Qtr.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (8*L1*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k)*Qt1(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (4*L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k)*Qt2(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
     (16*Qt2.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (L1*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (16*Qt1.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  Qtr.diff(d010,i,j,k)*(-4*(1 - pow(1 - rgrid[i],2)) - 
     2*(-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2)) - 
     16*pow(1 - rgrid[i],2) - 
     (16*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        ((pow(-1 + rgrid[i],2)*
             (3*pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],2) + 
               3*pow(H,2)*pow(-2 + rgrid[i],2)*pow(rgrid[i],3) + 
               pow(rh,2)*(3*pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],2) + 
                  3*pow(mu,2)*pow(-2 + rgrid[i],2)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (2 + 6*rgrid[i] - 12*pow(rgrid[i],2) + 
                     4*pow(rgrid[i],3)))))/(4.*pow(rh,4)) + 
          ((-1 + rgrid[i])*(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           (2.*pow(rh,4))))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
     (4*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Q(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
     (4*Q.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (4*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qrr(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     (4*Qtt.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
     (4*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qtt(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
     (4*B(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) - 
     (4*B.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*
        (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
     (4*Qr1.diff(d001,i,j,k)*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/L2 + 
     (4*Qr2.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*pow(R0 + 
          (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  R.diff(d010,i,j,k)*((-8*Qtr.diff(d001,i,j,k)*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k))/
      (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*Qtr.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (4*L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k)*Qt1(i,j,k))/
      ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (8*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k)*Qt2(i,j,k))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
     (16*Qt2.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        ((-3*pow(L1,2)*pow(1 + 
               B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr1(i,j,k),2))/(4.*pow(rh,2)) + 
          (pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2))/(4.*pow(rh,2)) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (16*Qt1.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        ((3*pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
              2)*pow(1 - pow(1 - rgrid[i],2),2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr1(i,j,k),2))/(4.*pow(rh,2)) - 
          (pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2))/(4.*pow(rh,2)) + 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  Qr1.diff(d001,i,j,k)*((-4*Qtr.diff(d001,i,j,k)*pow(L1,2)*
        pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      pow(L2,2) + (4*pow(L1,2)*
        pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qt2(i,j,k)*
        (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/L2 - 
     (4*pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qt1(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/L2 + 
     (8*Qt1.diff(d001,i,j,k)*pow(L1,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      pow(L2,2)) + (Qtr.diff(d100,i,j,k)*L1*
     (1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
          ((1 - pow(1 - rgrid[i],2))*
             (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
            2*Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        pow(rh,4) - pow(1 - pow(1 - rgrid[i],2),2)*
        ((pow(-1 + rgrid[i],2)*
             (3*pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],2) + 
               3*pow(H,2)*pow(-2 + rgrid[i],2)*pow(rgrid[i],3) + 
               pow(rh,2)*(3*pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],2) + 
                  3*pow(mu,2)*pow(-2 + rgrid[i],2)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (2 + 6*rgrid[i] - 12*pow(rgrid[i],2) + 
                     4*pow(rgrid[i],3)))))/(4.*pow(rh,4)) + 
          ((-1 + rgrid[i])*(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           (2.*pow(rh,4)))*Qt1(i,j,k)*
        ((1 - pow(1 - rgrid[i],2))*
           (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
          2*Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
       4*(1 - rgrid[i])*(L2*R0*Qr2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*Qr1(i,j,k)*sqrt((1 + pow(R0,2))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))))/
   (4.*pow(rh,2)*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),
        2))) - (L1*(-2*(1 - pow(1 - rgrid[i],2)) + 
       4*pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qtr(i,j,k)*
     (pow(1 - pow(1 - rgrid[i],2),3)*
        ((pow(-1 + rgrid[i],2)*
             (3*pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],2) + 
               3*pow(H,2)*pow(-2 + rgrid[i],2)*pow(rgrid[i],3) + 
               pow(rh,2)*(3*pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],2) + 
                  3*pow(mu,2)*pow(-2 + rgrid[i],2)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (2 + 6*rgrid[i] - 12*pow(rgrid[i],2) + 
                     4*pow(rgrid[i],3)))))/(4.*pow(rh,4)) + 
          ((-1 + rgrid[i])*(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           (2.*pow(rh,4)))*Qt1(i,j,k)*
        (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
       (2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
          Qtr(i,j,k)*sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        pow(rh,4) - 2*pow(1 - pow(1 - rgrid[i],2),2)*Qt1(i,j,k)*
        (((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)))/
           (2.*pow(rh,4)) + 
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
              (2.*pow(rh,4)))*Qtr(i,j,k))*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
       4*(1 - rgrid[i])*(L2*R0*Qr2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*Qr1(i,j,k)*sqrt((1 + pow(R0,2))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))))/
   (8.*pow(rh,2)*(1 - rgrid[i])*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (2*L1*pow(rh,4)*((-2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
           ((pow(-1 + rgrid[i],2)*
                (6*pow(H,2)*pow(-2 + rgrid[i],3)*rgrid[i] + 
                  18*pow(H,2)*pow(-2 + rgrid[i],2)*
                   pow(rgrid[i],2) + 
                  6*pow(H,2)*(-2 + rgrid[i])*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (6*pow(mu,2)*pow(-2 + rgrid[i],3)*rgrid[i] + 
                     18*pow(mu,2)*pow(-2 + rgrid[i],2)*
                      pow(rgrid[i],2) + 
                     6*pow(mu,2)*(-2 + rgrid[i])*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (6 - 24*rgrid[i] + 12*pow(rgrid[i],2)))))/
              (4.*pow(rh,4)) + 
             ((-1 + rgrid[i])*
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
                        4*pow(rgrid[i],3)))))/pow(rh,4) + 
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   4*pow(rh,2)*
                    (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                      4*pow(rgrid[i],3) + pow(rgrid[i],4))))/
              (2.*pow(rh,4))) - 
          2*(2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
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
              (2.*pow(rh,4))))*Qt1(i,j,k)*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
       ((-6*(1 - pow(1 - rgrid[i],2)) - 4*pow(1 - rgrid[i],2))*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
          pow(-1 + rgrid[i],6)*
          pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))),3)*
          Qt1(i,j,k)*((1 - pow(1 - rgrid[i],2))*
             (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
            Qtr(i,j,k))*Qtr(i,j,k)*
          sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        (16.*pow(rh,10)) - ((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
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
             (2.*pow(rh,4)))*Qtr(i,j,k)*
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
                (2.*pow(rh,4)))*Qt1(i,j,k)*
             (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
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
                (2.*pow(rh,4)))*Qt1(i,j,k)*Qtr(i,j,k)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            4*(1 - rgrid[i])*(L2*R0*Qr2(i,j,k)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
+ L1*Qr1(i,j,k)*sqrt((1 + pow(R0,2))*
                  (1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
))))/pow(rh,2) + (pow(-1 + rgrid[i],4)*
          pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*Qtr(i,j,k)*
          (-(pow(1 - pow(1 - rgrid[i],2),4)*
               (2*(1 - rgrid[i])*
                  ((pow(-1 + rgrid[i],2)*
                       (6*pow(H,2)*pow(-2 + rgrid[i],3)*
                       rgrid[i] + 
                        18*pow(H,2)*pow(-2 + rgrid[i],2)*
                       pow(rgrid[i],2) + 
                        6*pow(H,2)*(-2 + rgrid[i])*
                       pow(rgrid[i],3) + 
                        pow(rh,2)*
                        (6*pow(mu,2)*pow(-2 + rgrid[i],3)*
                       rgrid[i] + 
                        18*pow(mu,2)*pow(-2 + rgrid[i],2)*
                       pow(rgrid[i],2) + 
                        6*pow(mu,2)*(-2 + rgrid[i])*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (6 - 24*rgrid[i] + 12*pow(rgrid[i],2)))))/
                     (4.*pow(rh,4)) + 
                    ((-1 + rgrid[i])*
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
                        4*pow(rgrid[i],3)))))/pow(rh,4) + 
                    (pow(H,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))/
                     (2.*pow(rh,4))) + 
                 8*((pow(-1 + rgrid[i],2)*
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
                     (2.*pow(rh,4))))*Qt1(i,j,k)*
               (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))*
               sqrt(1 + pow(R0,2) + 
                 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                 pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
            16*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
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
                (2.*pow(rh,4)))*Qt1(i,j,k)*Qtr(i,j,k)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            pow(1 - pow(1 - rgrid[i],2),3)*Qt1(i,j,k)*
             (16*L1*pow(1 - rgrid[i],2)*
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
                   (2.*pow(rh,4)))*Qr1(i,j,k)*Qt1(i,j,k) + 
               16*L2*pow(1 - rgrid[i],2)*
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
                   (2.*pow(rh,4)))*Qr2(i,j,k)*Qt2(i,j,k) + 
               (2*(1 - rgrid[i])*
                   ((pow(-1 + rgrid[i],2)*
                       (6*pow(H,2)*pow(-2 + rgrid[i],3)*
                       rgrid[i] + 
                       18*pow(H,2)*pow(-2 + rgrid[i],2)*
                       pow(rgrid[i],2) + 
                       6*pow(H,2)*(-2 + rgrid[i])*
                       pow(rgrid[i],3) + 
                       pow(rh,2)*
                       (6*pow(mu,2)*pow(-2 + rgrid[i],3)*
                       rgrid[i] + 
                       18*pow(mu,2)*pow(-2 + rgrid[i],2)*
                       pow(rgrid[i],2) + 
                       6*pow(mu,2)*(-2 + rgrid[i])*
                       pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                       (6 - 24*rgrid[i] + 12*pow(rgrid[i],2)))))/
                      (4.*pow(rh,4)) + 
                     ((-1 + rgrid[i])*
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
                       4*pow(rgrid[i],3)))))/pow(rh,4) + 
                     (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))/
                      (2.*pow(rh,4))) + 
                  8*((pow(-1 + rgrid[i],2)*
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
                      (2.*pow(rh,4))))*Qtr(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            24*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             (L2*R0*Qr2(i,j,k)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
+ L1*Qr1(i,j,k)*sqrt((1 + pow(R0,2))*
                  (1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
)) - 16*pow(1 - rgrid[i],3)*
             (L2*R0*Qr2(i,j,k)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
               L1*Qr1(i,j,k)*sqrt((1 + pow(R0,2))*
                  (1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))\
)))/(16.*pow(rh,6))))/
   ((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (16*Qrr.diff(d001,i,j,k)*Qt1.diff(d001,i,j,k)*L1*pow(rh,2)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (pow(L2,2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  4*L1*B(i,j,k)*pow(1 - rgrid[i],2)*Qt2(i,j,k)*
   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
   sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
  Q.diff(d001,i,j,k)*((4*L1*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*Qt1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
     (16*Qt1.diff(d001,i,j,k)*L1*pow(rh,2)*
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
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  Qtr.diff(d001,i,j,k)*((4*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*R(i,j,k))/
      (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/L2) + 
  Qr2.diff(d010,i,j,k)*(-4*L2*(1 - pow(1 - rgrid[i],2))*
      pow(1 - rgrid[i],2)*Qt2(i,j,k)*
      (-1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
     (4*Qtr.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        (-1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
     (8*Qt1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*
        (-(L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (-1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (8*Qt2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*
        (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (-1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*Qtr.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))) + 
     (4*L2*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qt1(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))) + 
  Qr1.diff(d010,i,j,k)*(4*L1*(1 - pow(1 - rgrid[i],2))*
      pow(1 - rgrid[i],2)*Qt1(i,j,k)*
      (2 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
     (8*Qt1.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*
        (-(L2*Qr2(i,j,k)*pow(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2),1.5)) \
- L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + pow(R0,3) + 
             (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
             pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3))))/
      (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (8*Qt2.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*
        (L2*Qr2(i,j,k)*pow(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2),1.5) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + pow(R0,3) + 
             (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
             pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3))))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) - 
     (4*Qtr.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        (2 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
     (4*Qtr.diff(d001,i,j,k)*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/L2 - 
     4*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
      (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qt2(i,j,k)*
      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  B.diff(d100,i,j,k)*(-((Qt1.diff(d100,i,j,k)*L1*
          pow(1 - pow(1 - rgrid[i],2),2)*
          (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))) + 
     (2*Qtr.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*(1 + pow(R0 + 
            (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) - 
     (2*L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qt1(i,j,k)*
        (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
     Qt2.diff(d100,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) - 
     (2*Qtr.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/L2 + 
     2*L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qt2(i,j,k)*
      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) - 
     (4*Qt2.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*(L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
     (4*Qt1.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*(L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))))) + 
  Qt1.diff(d001,i,j,k)*(2*L1*(1 - pow(1 - rgrid[i],2))*
      (-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*Qr2(i,j,k) \
+ (4*L1*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*
        ((2*pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           pow(rh,4) + (1 - pow(1 - rgrid[i],2))*
           ((pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
              pow(rh,4) + 2*(1 - rgrid[i])*
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
                 (2.*pow(rh,4)))))*Qr2(i,j,k))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))) - 
     (4*L1*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k)*Qrr(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) \
- (8*pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k)*R(i,j,k))/
      (L2*sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
     (8*L1*B(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(L1*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))) + 
  B.diff(d001,i,j,k)*((4*L1*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*Qt1(i,j,k)*
        (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
     (4*Qtr.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/L2 - 
     4*L1*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
      Qr2(i,j,k)*Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) - 
     (16*Qt1.diff(d001,i,j,k)*L1*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        ((L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
           (4.*pow(rh,2)) + 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
             (L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              (4.*pow(rh,2))) + 
          (pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (2.*pow(rh,2)))*sqrt(1 + 
          pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  Qt2.diff(d010,i,j,k)*(2*L2*(1 - pow(1 - rgrid[i],2))*
      (-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*Qr2(i,j,k) \
+ (4*L2*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*
        ((2*pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           pow(rh,4) + (1 - pow(1 - rgrid[i],2))*
           ((pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
              pow(rh,4) + 2*(1 - rgrid[i])*
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
                 (2.*pow(rh,4)))))*Qr2(i,j,k))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))) - 
     (4*L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k)*Qrr(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) \
+ (8*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k)*R(i,j,k))/
      sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
     (16*R.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (-(L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
              pow(-1 + rgrid[i],2)*
              (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   4*pow(rh,2)*
                    (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                      4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
              Qr2(i,j,k))/(2.*pow(rh,2)) + 
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))/
      (L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (8*Qr2.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*
        (-(L2*Qr2(i,j,k)*pow(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2),1.5)) \
- L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + pow(R0,3) + 
             (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
             pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3))))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) - 
     (8*Qr1.diff(d001,i,j,k)*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/L2 - 
     (16*Qrr.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (16*Q.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (8*B(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
     4*B.diff(d001,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
      pow(1 - rgrid[i],2)*sqrt(1 + 
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
      ((4*pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
         (pow(L2,2)*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
        (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*Qr1(i,j,k)*Qr2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/L2 + 
        (2*pow(1 - pow(1 - rgrid[i],2),2)*pow(Qr2(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))))) + 
  Qt1.diff(d100,i,j,k)*(-((L1*
          (-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2)))/
        (1 - rgrid[i])) - (2*L1*(1 - pow(1 - rgrid[i],2))*
        (1 - rgrid[i])*Q(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) \
+ (4*Q.diff(d010,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (2*Q.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr2(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (2*L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qrr(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     L1*(1 - pow(1 - rgrid[i],2))*
      (-2/(1 - rgrid[i]) - (4*pow(rh,4)*
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
              (2.*pow(rh,4))))/
         (pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
        ((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
           pow(-1 + rgrid[i],4)*
           pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
           ((1 - pow(1 - rgrid[i],2))*
              (2*L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
             Qtr(i,j,k))*Qtr(i,j,k))/(4.*pow(rh,6)) + 
        (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
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
              (2.*pow(rh,4)))*Qtr(i,j,k)*
           (-((1 - pow(1 - rgrid[i],2))*
                (2*L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))) \
+ Qtr(i,j,k)))/(4.*pow(rh,2))) - 
     (2*Qtt.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr2(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
     (2*L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtt(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
     (2*R.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr2(i,j,k))/
      ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (2*L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
     (2*B.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr2(i,j,k)*
        (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) - 
     (2*Qr1.diff(d001,i,j,k)*pow(L1,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/L2 - 
     (4*Qtt.diff(d010,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr1(i,j,k)*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      ((1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (2*Qr2.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
     (2*Qr1.diff(d010,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*(2 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
     (2*Qr2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
     B.diff(d010,i,j,k)*((2*L1*pow(1 - pow(1 - rgrid[i],2),3)*
           (1 - rgrid[i])*Qr1(i,j,k)*
           (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
        (2*pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2))) + 
  Qt2.diff(d100,i,j,k)*((2*Q.diff(d010,i,j,k)*L2*
        pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*Qr2(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (L1*L2*pow(1 - pow(1 - rgrid[i],2),3)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
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
              (2.*pow(rh,4))))*Qr2(i,j,k)*Qt1(i,j,k)*Qtr(i,j,k))/
      (4.*pow(rh,2)) - (2*Qtt.diff(d010,i,j,k)*L2*
        pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*Qr2(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
     (2*Qr1.diff(d001,i,j,k)*pow(L1,2)*
        pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/L2 + 
     (4*R.diff(d010,i,j,k)*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*Qr1(i,j,k))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
     (2*R.diff(d001,i,j,k)*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*Qr2(i,j,k))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) - 
     (2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*R(i,j,k))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
     2*Qr2.diff(d001,i,j,k)*L1*
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
      pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
      sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
     (2*Qr2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*(-1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
     2*L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) - 
     2*Qr1.diff(d010,i,j,k)*L1*
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
      pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) - 
     2*B.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
      (1 - rgrid[i])*Qr2(i,j,k)*
      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
     B.diff(d010,i,j,k)*(-2*L1*pow(1 - pow(1 - rgrid[i],2),3)*
         (1 - rgrid[i])*Qr1(i,j,k)*
         (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
         sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) - 
        (2*pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))))) + 
  B.diff(d010,i,j,k)*((4*L1*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qt1(i,j,k)*
        (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) - 
     4*L1*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
      Qr1(i,j,k)*Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
      sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
     (4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qt1(i,j,k)*(L2*Qr2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) - 
     (4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qt2(i,j,k)*(L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) - 
     (4*Qt1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*(3*pow(L1,2)*
           pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(Qr1(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(L2,2)*pow(Qr2(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          4*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L2*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) + 
     4*Qt2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),4)*
      pow(1 - rgrid[i],2)*sqrt(1 + 
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
      ((3*L1*pow(Qr1(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/L2 + 
        (L2*pow(Qr2(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
         (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) + 
        (4*Qr1(i,j,k)*Qr2(i,j,k)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))) + 
     Qtr.diff(d010,i,j,k)*((-4*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],2)*Qr1(i,j,k)*
           (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) - 
        (4*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2))) + 
     Qtr.diff(d001,i,j,k)*((4*L1*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],2)*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/L2 \
+ (4*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L2 + L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))))
;

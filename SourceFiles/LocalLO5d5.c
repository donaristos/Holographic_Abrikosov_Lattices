
elem[68]+=
-2/(1 - pow(1 - rgrid[i],2)) - 
  (32*R.diff(d011,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   (L1*L2*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
  (2*pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
     Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   (pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
  (8*pow(B(i,j,k),2)*pow(1 - rgrid[i],2)*
     pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
   pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) + 
  (8*pow(1 - rgrid[i],2)*pow(R(i,j,k),2)*
     pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
   pow(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2),2) - 
  (4*pow(1 - rgrid[i],2)*pow(R(i,j,k),2))/
   (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
  (8*pow(1 - rgrid[i],2)*R(i,j,k)*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   ((1 - pow(1 - rgrid[i],2))*
     (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
  (4*pow(B(i,j,k),2)*pow(1 - rgrid[i],2)*
     (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
   pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) - 
  4*pow(Qr1.diff(d010,i,j,k),2)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],2)*(1 + 3*pow(R0,2) + 
     6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
     3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
  4*pow(Qr2.diff(d001,i,j,k),2)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],2)*(1 + 3*pow(R0,2) + 
     6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
     3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
  (4*pow(Qr2.diff(d010,i,j,k),2)*pow(L2,2)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
     (1 + 3*pow(R0,2) + 6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
   (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) - 
  (4*pow(Qr1.diff(d001,i,j,k),2)*pow(L1,2)*
     pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
     (1 + 3*pow(R0,2) + 6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/pow(L2,2) \
- (2*Qr2.diff(d100,i,j,k)*Qrr.diff(d010,i,j,k)*L2*
     pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
  (2*Qr1.diff(d100,i,j,k)*Qrr.diff(d001,i,j,k)*L1*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   (L2*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (16*R.diff(d020,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (16*R.diff(d002,i,j,k)*pow(rh,2)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   (pow(L2,2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (4*Qtt.diff(d001,i,j,k)*L1*pow(rh,4)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
           (2.*pow(rh,4))))*Qr1(i,j,k)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k),2)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(Qt1(i,j,k),2)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
     pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
   (2.*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(1 - rgrid[i],2)*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(Qt2(i,j,k),2)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
     pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
   (2.*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  pow(R.diff(d100,i,j,k),2)*
   ((2*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      pow(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2),2) - 
     pow(1 - pow(1 - rgrid[i],2),2)/
      (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  pow(B.diff(d100,i,j,k),2)*
   ((-2*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) - 
     (pow(1 - pow(1 - rgrid[i],2),2)*
        (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) + 
  pow(B.diff(d010,i,j,k),2)*
   ((-8*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        pow(Qr1(i,j,k),2)*pow(R0 + 
          (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) - 
     (4*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        pow(Qr1(i,j,k),2)*(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) + 
  pow(B.diff(d001,i,j,k),2)*
   ((-8*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        pow(Qr2(i,j,k),2)*pow(R0 + 
          (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) - 
     (4*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        pow(Qr2(i,j,k),2)*(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) + 
  Qtt.diff(d010,i,j,k)*((16*Qtt.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k),2)) + 
     (4*L2*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  pow(Qtt.diff(d010,i,j,k),2)*
   ((4*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (1 + 3*pow(R0,2) + 6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  pow(Qtt.diff(d001,i,j,k),2)*
   ((4*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (1 + 3*pow(R0,2) + 6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  Qr2.diff(d001,i,j,k)*((-8*B(i,j,k)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*
        (1 + 3*pow(R0,2) + 6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) - 
     (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qr1(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
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
              (2.*pow(rh,4)))*Qt2(i,j,k)*
           (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) + 
          ((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt2(i,j,k)*Qtr(i,j,k))/pow(rh,4) - 
          pow(1 - pow(1 - rgrid[i],2),2)*Qt2(i,j,k)*
           (((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)))/
              pow(rh,4) + ((pow(-1 + rgrid[i],2)*
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
                 (2.*pow(rh,4)))*Qtr(i,j,k)) + 
          4*(1 - rgrid[i])*(L1*R0*Qr1(i,j,k) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2)))))/
      (2.*pow(rh,2)*sqrt(1 + 
          pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  Qrr.diff(d001,i,j,k)*((4*Qr1.diff(d001,i,j,k)*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L2*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*
        (-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
        Qr1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      ((L2 + L2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (((-4*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]) + 
               8*pow(1 - rgrid[i],3))*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           (2.*pow(rh,4)) - 4*(1 - pow(1 - rgrid[i],2))*
           pow(1 - rgrid[i],2)*
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
              (2.*pow(rh,4))) - 
          (pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*
             pow(-1 + rgrid[i],6)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),3)*
             pow(-((1 - pow(1 - rgrid[i],2))*
                  (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))\
) + Qtr(i,j,k),2))/(4.*pow(rh,10)) + 
          (pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
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
                       L2*Qr2(i,j,k)*Qt2(i,j,k))) + Qtr(i,j,k),2) + 
               4*(1 - rgrid[i])*
                (2*L1*L2*R0*Qr1(i,j,k)*Qr2(i,j,k) + 
                  pow(L1,2)*pow(Qr1(i,j,k),2)*sqrt(1 + pow(R0,2)) + 
                  pow(L2,2)*pow(Qr2(i,j,k),2)*sqrt(1 + pow(R0,2)))))/
           (4.*pow(rh,6))))/
      (L2*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  Qrr.diff(d010,i,j,k)*((16*Qrr.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
     (4*Qr2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (L2*(1 - pow(1 - rgrid[i],2))*
        (-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
        Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (4*Qr2.diff(d001,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (L2*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (((-4*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]) + 
               8*pow(1 - rgrid[i],3))*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           (2.*pow(rh,4)) - 4*(1 - pow(1 - rgrid[i],2))*
           pow(1 - rgrid[i],2)*
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
              (2.*pow(rh,4))) - 
          (pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*
             pow(-1 + rgrid[i],6)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),3)*
             pow(-((1 - pow(1 - rgrid[i],2))*
                  (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))\
) + Qtr(i,j,k),2))/(4.*pow(rh,10)) + 
          (pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
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
                       L2*Qr2(i,j,k)*Qt2(i,j,k))) + Qtr(i,j,k),2) + 
               4*(1 - rgrid[i])*
                (2*L1*L2*R0*Qr1(i,j,k)*Qr2(i,j,k) + 
                  pow(L1,2)*pow(Qr1(i,j,k),2)*sqrt(1 + pow(R0,2)) + 
                  pow(L2,2)*pow(Qr2(i,j,k),2)*sqrt(1 + pow(R0,2)))))/
           (4.*pow(rh,6))))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  (4*h2.diff(d001,i,j,k)*q*pow(rh,2)*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (8*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
       (4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
        (pow(1 + pow(R0,2),0.25)*
          sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
       4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (-2*Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
             (1 + 2*pow(R0,2) + 
               4*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               2*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
       (4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        pow(1 + pow(R0,2),0.25)))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  (4*h1.diff(d001,i,j,k)*q*pow(rh,2)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (-8*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
       (4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
        (pow(1 + pow(R0,2),0.25)*
          sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
       4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (2*Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
             (1 + 2*pow(R0,2) + 
               4*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               2*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
       (4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        pow(1 + pow(R0,2),0.25)))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
  (2*pow(a0(i,j,k),2)*pow(1 - pow(1 - rgrid[i],2),4)*
     pow(1 - rgrid[i],2)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
     (pow(Qt1(i,j,k),2)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
       pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        pow(Qt2(i,j,k),2)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
       2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt1(i,j,k)*Qt2(i,j,k)*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   (pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
  (2*pow(a0(i,j,k),2)*pow(1 - pow(1 - rgrid[i],2),4)*
     pow(1 - rgrid[i],2)*(-2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        Qt1(i,j,k)*Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) \
+ pow(Qt1(i,j,k),2)*sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
       pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        pow(Qt2(i,j,k),2)*sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   (pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  h2.diff(d010,i,j,k)*((32*h2.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (4*q*pow(rh,2)*h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        ((8*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
           pow(1 + pow(R0,2),0.25) - 
          (4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (-2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             (Qt1(i,j,k)*(1 + 2*pow(R0,2) + 
                  4*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  2*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
)/sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
          4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  h1.diff(d010,i,j,k)*((32*h1.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (4*q*pow(rh,2)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        ((-8*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
           pow(1 + pow(R0,2),0.25) + 
          (4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             (Qt1(i,j,k)*(1 + 2*pow(R0,2) + 
                  4*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  2*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
)/sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
          4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(a2.diff(d010,i,j,k),2)*
   ((-2*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(Qr1(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(L2,2)*pow(Qr2(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(L1,2)*sqrt(1 + pow(R0,2))*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (2*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        (2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(Qr1(i,j,k),2)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          pow(L2,2)*pow(Qr2(i,j,k),2)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(L1,2)*sqrt(1 + pow(R0,2))*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(Qt2.diff(d010,i,j,k),2)*
   ((pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(Qr1(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(L2,2)*pow(Qr2(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*pow(L1,2)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(Qr1(i,j,k),2)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          pow(L2,2)*pow(Qr2(i,j,k),2)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*pow(L1,2)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(Qt1.diff(d001,i,j,k),2)*
   ((pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(Qr1(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(L2,2)*pow(Qr2(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*pow(L2,2)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(Qr1(i,j,k),2)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          pow(L2,2)*pow(Qr2(i,j,k),2)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*pow(L2,2)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(a1.diff(d001,i,j,k),2)*
   ((-2*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*sqrt(1 + pow(R0,2))*
        (pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(Qr1(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(L2,2)*pow(Qr2(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(L2,2)*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (2*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        sqrt(1 + pow(R0,2))*(2*L1*L2*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(Qr1(i,j,k),2)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          pow(L2,2)*pow(Qr2(i,j,k),2)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(L2,2)*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(R.diff(d010,i,j,k),2)*
   (4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
      (-((pow(1 - pow(1 - rgrid[i],2),2)*pow(Qr1(i,j,k),2))/
           (1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (4*pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
         (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     (48*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)*
        ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
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
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        pow(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2),2.5)) \
- (16*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)*
        ((pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr1(i,j,k),2))/(4.*pow(rh,2)) + 
          2*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        pow(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2),2))) + 
  pow(R.diff(d001,i,j,k),2)*
   ((48*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)*
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
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
           (4.*pow(rh,2))))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        pow(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2),2.5)) \
- (16*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
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
           (4.*pow(rh,2))))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        pow(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2),1.5)) \
- (16*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)*
        ((pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2))/(4.*pow(rh,2)) + 
          2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        pow(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2),2))) + 
  ar.diff(d001,i,j,k)*((4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(rh,2)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(rh,2)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  Qt1.diff(d010,i,j,k)*((2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*(-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
              (2.*pow(rh,4))))*Qr2(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (-(L1*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                 (pow(H,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   pow(rh,2)*
                    (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                 (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                 ((1 - pow(1 - rgrid[i],2))*
                    (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                      L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)))/
              (4.*pow(rh,2)) + 
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) - 
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
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
              (2.*pow(rh,4))))*Qr2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (-(Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        (L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(a0.diff(d100,i,j,k),2)*
   (-(pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
         (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
         (pow(Qt1(i,j,k),2)*
            (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
           pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
            pow(Qt2(i,j,k),2)*
            (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
           2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt1(i,j,k)*
            Qt2(i,j,k)*sqrt(1 + pow(R0,2) + 
              2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
              pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
        (-2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt1(i,j,k)*
           Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(Qt1(i,j,k),2)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(Qt2(i,j,k),2)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) - 
  (2*Qt2.diff(d001,i,j,k)*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
     (-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
          pow(rh,4)) + (1 - pow(1 - rgrid[i],2))*
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
           (2.*pow(rh,4))))*Qr1(i,j,k)*
     ((L2*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
          ((1 - pow(1 - rgrid[i],2))*
             (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
            Qtr(i,j,k))*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
        (4.*pow(rh,2)) - (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (Qt1(i,j,k)*pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2) + 
          Qt1(i,j,k)*(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
   (L2*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  Qt1.diff(d001,i,j,k)*((pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(rh,4)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L2*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L2*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(rh,4)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
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
              (2.*pow(rh,4))))*Qr1(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (-(L1*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                 (pow(H,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   pow(rh,2)*
                    (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                 (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                 ((1 - pow(1 - rgrid[i],2))*
                    (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                      L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)))/
              (4.*pow(rh,2)) + 
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) - 
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      ((L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
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
              (2.*pow(rh,4))))*Qr1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     Qtr.diff(d001,i,j,k)*(-((pow(1 - pow(1 - rgrid[i],2),3)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
               L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (pow(L2,2)*pow(rh,4)*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) - 
        (pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L2,2)*pow(rh,4)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))) + 
  Qt2.diff(d010,i,j,k)*(-((pow(1 - pow(1 - rgrid[i],2),2)*
          pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
            L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (pow(rh,4)*(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     (pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(rh,4)*(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
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
              (2.*pow(rh,4))))*Qr2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
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
              (2.*pow(rh,4))))*Qr2(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        ((L2*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
             (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
               Qtr(i,j,k)))/(4.*pow(rh,2)) - 
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
                 (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                 Qr2(i,j,k))/(4.*pow(rh,2)) + 
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     Qtr.diff(d010,i,j,k)*(-((pow(1 - pow(1 - rgrid[i],2),3)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             (L2*Qr2(i,j,k)*(R0 + 
                  (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
               L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (pow(L1,2)*pow(rh,4)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
) - (pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L1,2)*pow(rh,4)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     Qtr.diff(d001,i,j,k)*((pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
        (pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     Qt1.diff(d001,i,j,k)*(-((pow(1 - pow(1 - rgrid[i],2),4)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             (pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(Qr1(i,j,k),2)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
               pow(L2,2)*pow(Qr2(i,j,k),2)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
               2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                Qr1(i,j,k)*Qr2(i,j,k)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (L1*L2*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) \
- (pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              pow(Qr1(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             pow(L2,2)*pow(Qr2(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))))) - 
  (pow(q,2)*pow(rh,2)*(pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (16*(1 + pow(R0,2))*pow(a1(i,j,k),2)*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
       16*pow(a2(i,j,k),2)*pow(1 + 
          B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
       32*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*a2(i,j,k)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
       16*pow(a0(i,j,k),2)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],4)*sqrt(1 + pow(R0,2))*
        (-2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt1(i,j,k)*
           Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(Qt1(i,j,k),2)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(Qt2(i,j,k),2)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
       8*a1(i,j,k)*(-4*a2(i,j,k)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2)) + 
          4*pow(1 + pow(R0,2),0.75)*a0(i,j,k)*
           (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))))/
   (2.*sqrt(1 + pow(R0,2))*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
  (pow(q,2)*pow(rh,2)*(pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
     ((16*pow(a2(i,j,k),2)*pow(1 + 
            B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
        sqrt(1 + pow(R0,2)) + 
       16*pow(a1(i,j,k),2)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0,2)) - 
       (32*a0(i,j,k)*a2(i,j,k)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
          ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
            Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        pow(1 + pow(R0,2),0.25) + 
       16*pow(a0(i,j,k),2)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],4)*
        (pow(Qt1(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(Qt2(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt1(i,j,k)*
           Qt2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
       8*a1(i,j,k)*(-4*a2(i,j,k)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
           (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (-(Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))))/
   (2.*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  a1.diff(d001,i,j,k)*(ar.diff(d001,i,j,k)*
      ((4*pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L2,2)*pow(rh,2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
        (4*pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L2,2)*pow(rh,2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) - 
     (4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        (L2*Qr2(i,j,k)*((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
      (L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
          L2*Qr2(i,j,k)*(-(Qt1(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
      (L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  (2*pow(rh,4)*R(i,j,k)*(-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
            (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
            (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
            (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
            (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                   pow(-1 + rgrid[i],2)*
                   (pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                   pow(Qt1(i,j,k),2))/pow(rh,4)) + 
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
                  (2.*pow(rh,4)))*pow(Qt1(i,j,k),2) + 
              4*(1 - rgrid[i])*sqrt(1 + pow(R0,2))))/pow(rh,4)) - 
       (pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*(1 - rgrid[i])*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                 pow(-1 + rgrid[i],2)*
                 (pow(H,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   pow(rh,2)*
                    (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                 pow(Qt2(i,j,k),2))/pow(rh,4)) + 
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
                (2.*pow(rh,4)))*pow(Qt2(i,j,k),2) + 
            4*(1 - rgrid[i])*sqrt(1 + pow(R0,2))))/pow(rh,4) + 
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        ((pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],6)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),3)*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(-((1 - pow(1 - rgrid[i],2))*
                  (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                    L2*Qr2(i,j,k)*Qt2(i,j,k))) + Qtr(i,j,k),2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
           (8.*pow(rh,10)*sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
          (2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
             (2 + (1 - pow(1 - rgrid[i],2))*(Qrr(i,j,k) + Qtt(i,j,k)))*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          (pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             ((1 - pow(1 - rgrid[i],2))*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                (2*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
                  2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
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
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
               4*pow(1 - rgrid[i],2)*
                ((1 - pow(1 - rgrid[i],2))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   (-Qrr(i,j,k) + Qtt(i,j,k)) + 
                  2*R0*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))))/
           (2.*pow(rh,4)*sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
          ((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
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
                     (pow(H,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
                   (2.*pow(rh,4)))*Qt1(i,j,k)*Qt2(i,j,k))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           pow(rh,4) - (pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             ((pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
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
                  (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
               (8*(1 - rgrid[i])*
                  (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                  Qt1(i,j,k)*Qt2(i,j,k)*
                  pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
               4*pow(rh,2)*(1 - rgrid[i])*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                ((2*L1*L2*R0*Qr1(i,j,k)*Qr2(i,j,k))/
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) + 
                  (pow(L1,2)*(1 + pow(R0,2))*pow(Qr1(i,j,k),2))/
                   sqrt((1 + pow(R0,2))*
                     (1 + pow(R0,2) + 
                       2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                       pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2))) + 
                  (pow(L2,2)*(1 + pow(R0,2))*pow(Qr2(i,j,k),2))/
                   sqrt((1 + pow(R0,2))*
                     (1 + pow(R0,2) + 
                       2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                       pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2)))) + 
               8*(1 - rgrid[i])*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
                Qt2(i,j,k)*sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (8.*pow(rh,8)))))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*(1 - pow(1 - rgrid[i],2))*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  R.diff(d010,i,j,k)*((16*Q.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (16*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k)*R(i,j,k)*pow(R0 + 
          (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      pow(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2),2) + 
     (8*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k)*R(i,j,k))/
      (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
     (8*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qr1(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
     R.diff(d001,i,j,k)*((32*pow(rh,2)*
           pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
         (L1*L2*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
        (64*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)*
           ((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                Qr2(i,j,k))/(4.*pow(rh,2)) - 
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))/
         (L1*L2*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           pow(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2),2)) + 
        (32*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (-(L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                 pow(-1 + rgrid[i],2)*
                 (pow(H,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   pow(rh,2)*
                    (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                 (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                 Qr2(i,j,k))/(4.*pow(rh,2)) + 
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))/
         (L1*L2*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))) + 
     ((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        (4*pow(L1,2)*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (4*R0*(1 - rgrid[i]) - 
             (pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qt1(i,j,k)*Qt2(i,j,k))/pow(rh,4) + 
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
                 (2.*pow(rh,4)))*Qt1(i,j,k)*Qt2(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          pow(L1,2)*L2*Qr2(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)*
           (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt1(i,j,k),2))/
                (pow(rh,4)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
))) + (pow(1 - pow(1 - rgrid[i],2),3)*
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
                   (2.*pow(rh,4)))*pow(Qt1(i,j,k),2))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (4*(1 + pow(R0,2))*(1 - rgrid[i]))/
              sqrt((1 + pow(R0,2))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))) \
- pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           ((8*L1*(1 - rgrid[i])*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (L2*pow(1 - pow(1 - rgrid[i],2),3)*
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
                   (2.*pow(rh,4)))*Qr2(i,j,k)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                pow(Qt2(i,j,k),2)*
                pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             4*L2*(1 - rgrid[i])*Qr2(i,j,k)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)*
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
                   pow(Qt2(i,j,k),2))/
                 (4.*pow(rh,4)*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2))) - 
                (1 + pow(R0,2))/
                 sqrt((1 + pow(R0,2))*
                   (1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)))) + 
             L2*pow(1 - pow(1 - rgrid[i],2),3)*
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
                 (2.*pow(rh,4)))*Qr2(i,j,k)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              pow(Qt2(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             4*L2*(1 - rgrid[i])*Qr2(i,j,k)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
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
                   pow(Qt2(i,j,k),2)*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)))/(4.*pow(rh,4)) - 
                sqrt((1 + pow(R0,2))*
                  (1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
))) - pow(L1,2)*L2*Qr2(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt1(i,j,k),2)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
))/pow(rh,4)) + pow(1 - pow(1 - rgrid[i],2),3)*
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
                 (2.*pow(rh,4)))*pow(Qt1(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             4*(1 - rgrid[i])*
              sqrt((1 + pow(R0,2))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))\
)/(pow(L1,3)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (2*pow(L1,2)*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (4*R0*(1 - rgrid[i]) - 
             (pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qt1(i,j,k)*Qt2(i,j,k))/pow(rh,4) + 
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
                 (2.*pow(rh,4)))*Qt1(i,j,k)*Qt2(i,j,k))*
           (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           (8*L1*(1 - rgrid[i])*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             L2*pow(1 - pow(1 - rgrid[i],2),3)*
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
                 (2.*pow(rh,4)))*Qr2(i,j,k)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              pow(Qt2(i,j,k),2)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             4*L2*(1 - rgrid[i])*Qr2(i,j,k)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
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
                   pow(Qt2(i,j,k),2)*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)))/(4.*pow(rh,4)) - 
                sqrt((1 + pow(R0,2))*
                  (1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))\
)) - pow(L1,2)*L2*Qr2(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                         4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt1(i,j,k),2)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
)/pow(rh,4)) + pow(1 - pow(1 - rgrid[i],2),3)*
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
                 (2.*pow(rh,4)))*pow(Qt1(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             4*(1 - rgrid[i])*sqrt((1 + pow(R0,2))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))))/
      (pow(L1,3)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        (L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        pow(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2),1.5))) \
+ R.diff(d001,i,j,k)*((-16*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*R(i,j,k)*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      pow(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2),2) + 
     (8*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k)*R(i,j,k))/
      (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
     (8*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qr2(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
     ((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        (4*pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (4*R0*(1 - rgrid[i]) - 
             (pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qt1(i,j,k)*Qt2(i,j,k))/pow(rh,4) + 
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
                 (2.*pow(rh,4)))*Qt1(i,j,k)*Qt2(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           Qr1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)*
           (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt2(i,j,k),2))/
                (pow(rh,4)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
))) + (pow(1 - pow(1 - rgrid[i],2),3)*
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
                   (2.*pow(rh,4)))*pow(Qt2(i,j,k),2))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (4*(1 + pow(R0,2))*(1 - rgrid[i]))/
              sqrt((1 + pow(R0,2))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))) \
- L1*((8*L2*(1 - rgrid[i])*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (L1*pow(1 - pow(1 - rgrid[i],2),3)*
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
                   (2.*pow(rh,4)))*Qr1(i,j,k)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                pow(Qt1(i,j,k),2)*
                pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             4*L1*(1 - rgrid[i])*Qr1(i,j,k)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)*
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
                   pow(Qt1(i,j,k),2))/
                 (4.*pow(rh,4)*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2))) - 
                (1 + pow(R0,2))/
                 sqrt((1 + pow(R0,2))*
                   (1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)))) + 
             L1*pow(1 - pow(1 - rgrid[i],2),3)*
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
                 (2.*pow(rh,4)))*Qr1(i,j,k)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              pow(Qt1(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             4*L1*(1 - rgrid[i])*Qr1(i,j,k)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
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
                   pow(Qt1(i,j,k),2)*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)))/(4.*pow(rh,4)) - 
                sqrt((1 + pow(R0,2))*
                  (1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
))) - pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*Qr1(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt2(i,j,k),2)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
))/pow(rh,4)) + pow(1 - pow(1 - rgrid[i],2),3)*
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
                 (2.*pow(rh,4)))*pow(Qt2(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             4*(1 - rgrid[i])*
              sqrt((1 + pow(R0,2))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))\
)/(L1*L2*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     ((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (2*pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qr1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (4*R0*(1 - rgrid[i]) - 
             (pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qt1(i,j,k)*Qt2(i,j,k))/pow(rh,4) + 
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
                 (2.*pow(rh,4)))*Qt1(i,j,k)*Qt2(i,j,k))*
           (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          L1*(8*L2*(1 - rgrid[i])*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             L1*pow(1 - pow(1 - rgrid[i],2),3)*
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
                 (2.*pow(rh,4)))*Qr1(i,j,k)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              pow(Qt1(i,j,k),2)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             4*L1*(1 - rgrid[i])*Qr1(i,j,k)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
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
                   pow(Qt1(i,j,k),2)*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)))/(4.*pow(rh,4)) - 
                sqrt((1 + pow(R0,2))*
                  (1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))\
)) - pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*Qr1(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                         4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt2(i,j,k),2)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
)/pow(rh,4)) + pow(1 - pow(1 - rgrid[i],2),3)*
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
                 (2.*pow(rh,4)))*pow(Qt2(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             4*(1 - rgrid[i])*sqrt((1 + pow(R0,2))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))))/
      (L2*(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        pow(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2),1.5))) \
- (2*pow(rh,4)*R(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
     (((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
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
            4*(1 - rgrid[i])*sqrt(1 + pow(R0,2))))/(2.*pow(rh,4)) + 
       (pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*(1 - rgrid[i])*
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
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
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
            4*(1 - rgrid[i])*sqrt(1 + pow(R0,2))))/(2.*pow(rh,4)) + 
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        ((pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],6)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
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
           (8.*pow(rh,10)) + 
          2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           (2 + (1 - pow(1 - rgrid[i],2))*(Qrr(i,j,k) + Qtt(i,j,k)))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          (pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             ((1 - pow(1 - rgrid[i],2))*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                (2*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
                  2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
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
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
               4*pow(1 - rgrid[i],2)*
                ((1 - pow(1 - rgrid[i],2))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   (-Qrr(i,j,k) + Qtt(i,j,k)) + 
                  2*R0*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (2.*pow(rh,4)) - (pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
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
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
               8*(1 - rgrid[i])*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
                Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
               4*pow(rh,2)*(1 - rgrid[i])*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                (2*L1*L2*R0*Qr1(i,j,k)*Qr2(i,j,k)*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
) + pow(L1,2)*pow(Qr1(i,j,k),2)*
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
                        pow(R(i,j,k),2))))))/(8.*pow(rh,8)))))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*(1 - pow(1 - rgrid[i],2))*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
     pow(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2),1.5)) + 
  (2*pow(rh,4)*(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
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
            4*(1 - rgrid[i])*sqrt(1 + pow(R0,2))))/(2.*pow(rh,4)) + 
       (pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*(1 - rgrid[i])*
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
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
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
            4*(1 - rgrid[i])*sqrt(1 + pow(R0,2))))/(2.*pow(rh,4)) + 
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        ((pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],6)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
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
           (8.*pow(rh,10)) + 
          2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           (2 + (1 - pow(1 - rgrid[i],2))*(Qrr(i,j,k) + Qtt(i,j,k)))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          (pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             ((1 - pow(1 - rgrid[i],2))*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                (2*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
                  2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
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
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
               4*pow(1 - rgrid[i],2)*
                ((1 - pow(1 - rgrid[i],2))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   (-Qrr(i,j,k) + Qtt(i,j,k)) + 
                  2*R0*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (2.*pow(rh,4)) - (pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
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
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
               8*(1 - rgrid[i])*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
                Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
               4*pow(rh,2)*(1 - rgrid[i])*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                (2*L1*L2*R0*Qr1(i,j,k)*Qr2(i,j,k)*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
) + pow(L1,2)*pow(Qr1(i,j,k),2)*
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
                        pow(R(i,j,k),2))))))/(8.*pow(rh,8)))))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  R.diff(d100,i,j,k)*((8*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        R(i,j,k)*pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      pow(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2),2) - 
     (4*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*R(i,j,k))/
      (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
     (4*(1 - rgrid[i])*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
     R.diff(d010,i,j,k)*((-8*pow(1 - pow(1 - rgrid[i],2),3)*
           (1 - rgrid[i])*Qr1(i,j,k)*
           pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
         pow(1 + pow(R0,2) + 
           2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
           pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2),2) + 
        (4*pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*Qr1(i,j,k))/
         (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
           pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
     R.diff(d001,i,j,k)*((-8*pow(1 - pow(1 - rgrid[i],2),3)*
           (1 - rgrid[i])*Qr2(i,j,k)*
           pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
         pow(1 + pow(R0,2) + 
           2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
           pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2),2) + 
        (4*pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*Qr2(i,j,k))/
         (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
           pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
     (pow(rh,4)*(-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
               (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
               (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
               (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                      pow(-1 + rgrid[i],2)*
                      (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                      pow(Qt1(i,j,k),2))/pow(rh,4)) + 
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
                     (2.*pow(rh,4)))*pow(Qt1(i,j,k),2) + 
                 4*(1 - rgrid[i])*sqrt(1 + pow(R0,2))))/pow(rh,4)) - 
          (pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
             (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                    pow(-1 + rgrid[i],2)*
                    (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      pow(rh,2)*
                       (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                    pow(Qt2(i,j,k),2))/pow(rh,4)) + 
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
                   (2.*pow(rh,4)))*pow(Qt2(i,j,k),2) + 
               4*(1 - rgrid[i])*sqrt(1 + pow(R0,2))))/pow(rh,4) + 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           ((pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
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
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              (8.*pow(rh,10)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) \
+ (2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
                (2 + (1 - pow(1 - rgrid[i],2))*
                   (Qrr(i,j,k) + Qtt(i,j,k)))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                ((1 - pow(1 - rgrid[i],2))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                   (2*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
                     2*pow(1 - pow(1 - rgrid[i],2),2)*
                      (1 - rgrid[i])*
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
                      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
                  4*pow(1 - rgrid[i],2)*
                   ((1 - pow(1 - rgrid[i],2))*
                      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                      (-Qrr(i,j,k) + Qtt(i,j,k)) + 
                     2*R0*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                      (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))))/
              (2.*pow(rh,4)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) \
+ ((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
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
                        (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
                      (2.*pow(rh,4)))*Qt1(i,j,k)*Qt2(i,j,k))*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              pow(rh,4) - (pow(1 - pow(1 - rgrid[i],2),2)*
                (1 - rgrid[i])*pow(-1 + rgrid[i],4)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                ((pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
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
                     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) + 
                  (8*(1 - rgrid[i])*
                     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                     Qt1(i,j,k)*Qt2(i,j,k)*
                     pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)\
)/sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) + 
                  4*pow(rh,2)*(1 - rgrid[i])*
                   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                   ((2*L1*L2*R0*Qr1(i,j,k)*Qr2(i,j,k))/
                      sqrt(1 + pow(R0,2) + 
                        2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                        pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2)) + 
                     (pow(L1,2)*(1 + pow(R0,2))*
                        pow(Qr1(i,j,k),2))/
                      sqrt((1 + pow(R0,2))*
                        (1 + pow(R0,2) + 
                        2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                        pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2))) + 
                     (pow(L2,2)*(1 + pow(R0,2))*
                        pow(Qr2(i,j,k),2))/
                      sqrt((1 + pow(R0,2))*
                        (1 + pow(R0,2) + 
                         2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                         pow(1 - pow(1 - rgrid[i],2),2)*
                         pow(R(i,j,k),2)))) + 
                  8*(1 - rgrid[i])*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                   Qt1(i,j,k)*Qt2(i,j,k)*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
)))/(8.*pow(rh,8)))))/
      ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (pow(rh,4)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
               4*(1 - rgrid[i])*sqrt(1 + pow(R0,2))))/(2.*pow(rh,4)) + 
          (pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
             (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
               4*(1 - rgrid[i])*sqrt(1 + pow(R0,2))))/(2.*pow(rh,4)) + 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           ((pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
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
              (8.*pow(rh,10)) + 
             2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
              (2 + (1 - pow(1 - rgrid[i],2))*(Qrr(i,j,k) + Qtt(i,j,k)))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                ((1 - pow(1 - rgrid[i],2))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                   (2*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
                     2*pow(1 - pow(1 - rgrid[i],2),2)*
                      (1 - rgrid[i])*
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
                      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
                  4*pow(1 - rgrid[i],2)*
                   ((1 - pow(1 - rgrid[i],2))*
                      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                      (-Qrr(i,j,k) + Qtt(i,j,k)) + 
                     2*R0*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                      (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                      (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (2.*pow(rh,4)) - 
             (pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                pow(-1 + rgrid[i],4)*
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
                     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
) + 8*(1 - rgrid[i])*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                   Qt1(i,j,k)*Qt2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
) + 4*pow(rh,2)*(1 - rgrid[i])*
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
                         pow(R(i,j,k),2))))))/(8.*pow(rh,8)))))/
      ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        pow(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2),1.5))) \
+ (pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(Qt1(i,j,k),2)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (2.*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(1 - rgrid[i],2)*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(Qt2(i,j,k),2)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (pow(rh,4)*(2 + 2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  pow(a2.diff(d100,i,j,k),2)*
   (-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
         pow(1 - pow(1 - rgrid[i],2),2)*
         pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (2.*sqrt(1 + pow(R0,2))*pow(rh,2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (2.*sqrt(1 + pow(R0,2))*pow(rh,2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(ar.diff(d010,i,j,k),2)*
   ((-2*pow(1 - rgrid[i],2)*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (pow(L1,2)*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (2*pow(1 - rgrid[i],2)*sqrt(1 + 
          pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,2)*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(ar.diff(d001,i,j,k),2)*
   ((-2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(1 - rgrid[i],2)*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (pow(L2,2)*pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(1 - rgrid[i],2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) \
+ pow(Qrr.diff(d010,i,j,k),2)*
   ((-4*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*sqrt(1 + 
          pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))) + 
  pow(Qrr.diff(d001,i,j,k),2)*
   ((-4*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))) + 
  pow(h1.diff(d010,i,j,k),2)*
   ((-8*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (8*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(h2.diff(d010,i,j,k),2)*
   ((-8*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (8*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(h1.diff(d001,i,j,k),2)*
   ((-8*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (8*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(h2.diff(d001,i,j,k),2)*
   ((-8*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (8*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(Qt1.diff(d100,i,j,k),2)*
   ((pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (8.*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (8.*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(Qtr.diff(d010,i,j,k),2)*
   ((pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (2.*pow(L1,2)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (2.*pow(L1,2)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(Qtr.diff(d001,i,j,k),2)*
   (((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (2.*pow(L2,2)*pow(rh,4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (2.*pow(L2,2)*pow(rh,4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(Qt2.diff(d100,i,j,k),2)*
   (((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (8.*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (4.*pow(rh,4)*(2 + 2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  Qtr.diff(d010,i,j,k)*((-2*Qtr.diff(d001,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (2*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(rh,4)*(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
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
              (2.*pow(rh,4))))*Qr2(i,j,k)*
        ((1 - pow(1 - rgrid[i],2))*
           (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
          Qtr(i,j,k))*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (2.*pow(rh,2)*(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     ((1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     Qt1.diff(d001,i,j,k)*((pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) - 
     ((1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  Qtr.diff(d001,i,j,k)*((2*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(rh,4)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
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
              (2.*pow(rh,4))))*Qr1(i,j,k)*
        ((1 - pow(1 - rgrid[i],2))*
           (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
          Qtr(i,j,k))*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (2.*L2*pow(rh,2)*sqrt(1 + 
          pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (pow(rh,4)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(rh,4)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  B.diff(d100,i,j,k)*((-8*B(i,j,k)*(1 - pow(1 - rgrid[i],2))*
        (1 - rgrid[i])*pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) - 
     (4*B(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) + 
     (4*Qr1.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*(1 + 3*pow(R0,2) + 
          6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) - 
     (4*Qr2.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*(1 + 3*pow(R0,2) + 
          6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
     (4*B.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr1(i,j,k)*
        (1 + 3*pow(R0,2) + 6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) + 
     (4*B.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr2(i,j,k)*
        (1 + 3*pow(R0,2) + 6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) - 
     (12*Qr1.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/L2 + 
     (12*Qr2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2))) + 
  Qr1.diff(d001,i,j,k)*((-2*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*(1 - rgrid[i])*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (4*(1 - rgrid[i]) - (L1*pow(1 - pow(1 - rgrid[i],2),3)*
             (1 - rgrid[i])*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             Qr1(i,j,k)*Qt1(i,j,k)*
             ((1 - pow(1 - rgrid[i],2))*
                (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
               Qtr(i,j,k)))/(4.*pow(rh,6)) + 
          (L1*pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qr1(i,j,k)*(pow(1 - pow(1 - rgrid[i],2),3)*
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
                (L2*R0*Qr2(i,j,k) + L1*Qr1(i,j,k)*sqrt(1 + pow(R0,2))))\
)/(4.*pow(rh,2))))/
      (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (24*L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/L2 - 
     (24*Qr2.diff(d001,i,j,k)*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/L2) + 
  B.diff(d010,i,j,k)*((-32*Q.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (8*B(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k)*(1 + 3*pow(R0,2) + 
          6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) - 
     (8*Qr1.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k)*
        (1 + 3*pow(R0,2) + 6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
     (8*Qr2.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k)*
        (1 + 3*pow(R0,2) + 6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
     (16*R.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     B.diff(d001,i,j,k)*((-8*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qr2(i,j,k)*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
         pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) - 
        (64*pow(rh,4)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
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
              (4.*pow(rh,2))))/
         (L1*L2*pow(rh + rh*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     (2*pow(L2,2)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
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
                sqrt(1 + pow(R0,2))))))/
      (pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),3)*
        (L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (24*Qr1.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/L2 - 
     (24*Qr2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2))) + 
  B.diff(d001,i,j,k)*((8*B(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*
        (1 + 3*pow(R0,2) + 6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) + 
     (8*Qr2.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*
        (1 + 3*pow(R0,2) + 6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) - 
     (16*R.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
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
                sqrt(1 + pow(R0,2))))))/
      (L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (24*Qr1.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/L2) + 
  pow(H,2)*((-8*pow(rh,2)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        ((pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
             pow(1 - pow(1 - rgrid[i],2),4)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             pow(Qr1(i,j,k),2)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
           (4.*pow(rh,4)) + 
          pow(L2,2)*((pow(1 - pow(1 - rgrid[i],2),4)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                pow(Qr2(i,j,k),2))/(4.*pow(rh,4)) + 
             pow(q,2)*(pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
              pow(1 - pow(1 - rgrid[i],2),2)*pow(ygrid[k],2)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          (L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),4)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qr1(i,j,k)*Qr2(i,j,k)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (2.*pow(rh,4))))/
      ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (8*pow(rh,2)*pow(1 - rgrid[i],2)*
        ((pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
             pow(1 - pow(1 - rgrid[i],2),4)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             pow(Qr1(i,j,k),2))/(4.*pow(rh,4)) + 
          pow(L2,2)*((pow(1 - pow(1 - rgrid[i],2),4)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                pow(Qr2(i,j,k),2))/(4.*pow(rh,4)) + 
             pow(q,2)*(pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
              pow(1 - pow(1 - rgrid[i],2),2)*pow(ygrid[k],2)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
          (L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),4)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qr1(i,j,k)*Qr2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
           (2.*pow(rh,4)*sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  pow(a0.diff(d001,i,j,k),2)*
   ((-32*pow(rh,4)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],6)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        ((pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),4)*
             pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             pow(Qr2(i,j,k),2)*pow(Qt1(i,j,k),2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
           (16.*pow(rh,6)) + 
          pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           (-1 - (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
             (pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],4)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                pow(-(L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                     Qt1(i,j,k)) + Qtr(i,j,k),2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
              (16.*pow(rh,6)))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          (L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),3)*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             Qr2(i,j,k)*Qt1(i,j,k)*
             (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
               Qtr(i,j,k))*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (8.*pow(rh,6))))/
      (pow(L2,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],4)*
        pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (32*pow(rh,4)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],6)*((pow(L2,2)*
             pow(1 - pow(1 - rgrid[i],2),4)*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             pow(Qr2(i,j,k),2)*pow(Qt1(i,j,k),2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
           (16.*pow(rh,6)) + 
          pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           (-1 - (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
             (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],4)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                pow(-(L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                     Qt1(i,j,k)) + Qtr(i,j,k),2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
              (16.*pow(rh,6))) + 
          (L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),3)*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             Qr2(i,j,k)*Qt1(i,j,k)*
             (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
               Qtr(i,j,k))*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
           (8.*pow(rh,6)*sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],4)*
        pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))) + 
  pow(a0.diff(d010,i,j,k),2)*
   ((-32*pow(rh,4)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],6)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        ((pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
             pow(1 - pow(1 - rgrid[i],2),4)*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             pow(Qr1(i,j,k),2)*pow(Qt2(i,j,k),2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
           (16.*pow(rh,6)) + 
          (-1 - (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
             (pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],4)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                pow(-(L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                     Qt2(i,j,k)) + Qtr(i,j,k),2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
              (16.*pow(rh,6)))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),3)*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             Qr1(i,j,k)*Qt2(i,j,k)*
             (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
               Qtr(i,j,k))*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (8.*pow(rh,6))))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],4)*
        pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (32*pow(rh,4)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],6)*(-1 - 
          (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
          (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
              2)*pow(1 - pow(1 - rgrid[i],2),4)*
             pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             pow(Qr1(i,j,k),2)*pow(Qt2(i,j,k),2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
           (16.*pow(rh,6)) + 
          (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             pow(-(L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                  Qt2(i,j,k)) + Qtr(i,j,k),2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
           (16.*pow(rh,6)) + 
          (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),3)*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             Qr1(i,j,k)*Qt2(i,j,k)*
             (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
               Qtr(i,j,k))*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
           (8.*pow(rh,6)*sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],4)*
        pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))) + 
  Q.diff(d001,i,j,k)*((-16*R.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     ((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (-2*pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (4*R0*(1 - rgrid[i]) - 
             (pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qt1(i,j,k)*Qt2(i,j,k))/pow(rh,4) + 
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
                 (2.*pow(rh,4)))*Qt1(i,j,k)*Qt2(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(L1,2)*(-((pow(1 - pow(1 - rgrid[i],2),2)*
                  (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt1(i,j,k),2)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*
                     pow(R(i,j,k),2)))/pow(rh,4)) + 
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
                 (2.*pow(rh,4)))*pow(Qt1(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             4*(1 - rgrid[i])*
              sqrt((1 + pow(R0,2))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))) \
+ pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt2(i,j,k),2)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
))/pow(rh,4)) + pow(1 - pow(1 - rgrid[i],2),3)*
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
                 (2.*pow(rh,4)))*pow(Qt2(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             4*(1 - rgrid[i])*
              sqrt((1 + pow(R0,2))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))\
)/(L1*L2*pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     ((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (-2*pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (4*R0*(1 - rgrid[i]) - 
             (pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qt1(i,j,k)*Qt2(i,j,k))/pow(rh,4) + 
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
                 (2.*pow(rh,4)))*Qt1(i,j,k)*Qt2(i,j,k)) + 
          pow(L1,2)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt1(i,j,k),2))/
                (pow(rh,4)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
))) + (pow(1 - pow(1 - rgrid[i],2),3)*
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
                   (2.*pow(rh,4)))*pow(Qt1(i,j,k),2))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (4*(1 + pow(R0,2))*(1 - rgrid[i]))/
              sqrt((1 + pow(R0,2))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))) \
+ pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt2(i,j,k),2))/
                (pow(rh,4)*sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
)) + (pow(1 - pow(1 - rgrid[i],2),3)*
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
                   (2.*pow(rh,4)))*pow(Qt2(i,j,k),2))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (4*(1 + pow(R0,2))*(1 - rgrid[i]))/
              sqrt((1 + pow(R0,2))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*L2*pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  Q.diff(d010,i,j,k)*((16*R.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (32*B.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (16*R.diff(d010,i,j,k)*pow(rh,2)*
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
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (-2*pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (4*R0*(1 - rgrid[i]) - 
             (pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qt1(i,j,k)*Qt2(i,j,k))/pow(rh,4) + 
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
                 (2.*pow(rh,4)))*Qt1(i,j,k)*Qt2(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(L1,2)*(-((pow(1 - pow(1 - rgrid[i],2),2)*
                  (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt1(i,j,k),2)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*
                     pow(R(i,j,k),2)))/pow(rh,4)) + 
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
                 (2.*pow(rh,4)))*pow(Qt1(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             4*(1 - rgrid[i])*
              sqrt((1 + pow(R0,2))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))) \
+ pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt2(i,j,k),2)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
))/pow(rh,4)) + pow(1 - pow(1 - rgrid[i],2),3)*
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
                 (2.*pow(rh,4)))*pow(Qt2(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             4*(1 - rgrid[i])*
              sqrt((1 + pow(R0,2))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))\
)/(pow(L1,3)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (-2*pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (4*R0*(1 - rgrid[i]) - 
             (pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                Qt1(i,j,k)*Qt2(i,j,k))/pow(rh,4) + 
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
                 (2.*pow(rh,4)))*Qt1(i,j,k)*Qt2(i,j,k)) + 
          pow(L1,2)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt1(i,j,k),2))/
                (pow(rh,4)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
))) + (pow(1 - pow(1 - rgrid[i],2),3)*
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
                   (2.*pow(rh,4)))*pow(Qt1(i,j,k),2))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (4*(1 + pow(R0,2))*(1 - rgrid[i]))/
              sqrt((1 + pow(R0,2))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))) \
+ pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (-((pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  pow(Qt2(i,j,k),2))/
                (pow(rh,4)*sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
)) + (pow(1 - pow(1 - rgrid[i],2),3)*
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
                   (2.*pow(rh,4)))*pow(Qt2(i,j,k),2))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (4*(1 + pow(R0,2))*(1 - rgrid[i]))/
              sqrt((1 + pow(R0,2))*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,3)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  pow(a1.diff(d100,i,j,k),2)*
   (-(pow(1 - pow(1 - rgrid[i],2),2)*
         pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)*
         sqrt(1 + pow(R0,2)))/
      (2.*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (pow(1 - pow(1 - rgrid[i],2),2)*sqrt(1 + pow(R0,2))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (2.*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  Qt1.diff(d100,i,j,k)*(-(Qt2.diff(d100,i,j,k)*
         pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
         (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
           pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
               pow(rgrid[i],3) + 
              4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                 4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
         (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
         (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (2.*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (Qtr.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L2*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     ((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     ((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (2.*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     Qt2.diff(d010,i,j,k)*((pow(1 - pow(1 - rgrid[i],2),3)*
           (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (2.*L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (2.*L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     Qt1.diff(d001,i,j,k)*(-(pow(1 - pow(1 - rgrid[i],2),3)*
            (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
            (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
            (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
            (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
              L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
               sqrt(1 + pow(R0,2) + 
                 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                 pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (2.*L2*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (2.*L2*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     ((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (2.*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     Qtr.diff(d010,i,j,k)*(-(pow(1 - pow(1 - rgrid[i],2),2)*
            (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
            (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
            pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
         (2.*L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (2.*L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))))) + 
  Qt2.diff(d100,i,j,k)*((Qtr.diff(d010,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     ((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
      (2.*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     Qt2.diff(d010,i,j,k)*(-(pow(1 - pow(1 - rgrid[i],2),3)*
            (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
            (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
            (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
              L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
               sqrt(1 + pow(R0,2) + 
                 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                 pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (2.*pow(rh,4)*(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
        (pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (2.*pow(rh,4)*(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     Qt1.diff(d001,i,j,k)*((pow(1 - pow(1 - rgrid[i],2),3)*
           (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (2.*pow(rh,4)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
        (pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (2.*pow(rh,4)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (2.*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     Qtr.diff(d001,i,j,k)*(-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
            pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
            pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
            (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
            pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
         (2.*pow(rh,4)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (2.*pow(rh,4)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))))) + 
  Qr1.diff(d010,i,j,k)*((8*B(i,j,k)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*
        (1 + 3*pow(R0,2) + 6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
     8*Qr2.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
      pow(1 - rgrid[i],2)*(1 + 3*pow(R0,2) + 
        6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
     (8*B.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*
        (1 + 3*pow(R0,2) + 6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
     (4*Qrr.diff(d001,i,j,k)*L1*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L2*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qr2(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
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
              (2.*pow(rh,4)))*Qt1(i,j,k)*
           (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) + 
          ((1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt1(i,j,k)*Qtr(i,j,k))/pow(rh,4) - 
          pow(1 - pow(1 - rgrid[i],2),2)*Qt1(i,j,k)*
           (((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)))/
              pow(rh,4) + ((pow(-1 + rgrid[i],2)*
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
                 (2.*pow(rh,4)))*Qtr(i,j,k)) + 
          4*(1 - rgrid[i])*(L2*R0*Qr2(i,j,k) + 
             L1*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))))/
      (2.*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (24*Qr2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))) + 
     Qr1.diff(d001,i,j,k)*((8*L1*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),3))/
         (L2*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) \
+ (16*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/L2)\
) + Qr2.diff(d010,i,j,k)*(8*Qr1.diff(d001,i,j,k)*
      pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
      (1 + 3*pow(R0,2) + 6*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
        3*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
     (2*L2*(1 - rgrid[i])*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (4*(1 - rgrid[i]) - (L2*pow(1 - pow(1 - rgrid[i],2),3)*
             (1 - rgrid[i])*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             Qr2(i,j,k)*Qt2(i,j,k)*
             ((1 - pow(1 - rgrid[i],2))*
                (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
               Qtr(i,j,k)))/(4.*pow(rh,6)) + 
          (L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
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
                   (2.*pow(rh,4)))*Qt2(i,j,k)*
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
                   (2.*pow(rh,4)))*Qt2(i,j,k)*Qtr(i,j,k) + 
               4*(1 - rgrid[i])*
                (L1*R0*Qr1(i,j,k) + L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2))))\
)/(4.*pow(rh,2))))/
      ((L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (24*L2*B(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) - 
     (24*B.diff(d001,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) + 
     Qr2.diff(d001,i,j,k)*((8*L2*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),3))/
         ((L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (16*L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2))))) + 
  ar.diff(d010,i,j,k)*((8*ar.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (-(Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     a1.diff(d001,i,j,k)*((-4*pow(1 + pow(R0,2),0.25)*
           pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (4*pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a0.diff(d001,i,j,k)*((4*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],4)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
                Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
           (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k) + 
             ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
                  Qtr(i,j,k))*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))\
)/sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L1*L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))))) + 
  a2.diff(d010,i,j,k)*(ar.diff(d010,i,j,k)*
      ((4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L1,2)*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L1,2)*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     ar.diff(d001,i,j,k)*((-4*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
        (4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     a1.diff(d001,i,j,k)*((4*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              pow(Qr1(i,j,k),2)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             pow(L2,2)*pow(Qr2(i,j,k),2)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              Qr1(i,j,k)*Qr2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (4*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
           (2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              Qr1(i,j,k)*Qr2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              pow(Qr1(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             pow(L2,2)*pow(Qr2(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) - 
     (4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*
        (-(L2*Qr2(i,j,k)*((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
               Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))) \
+ L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
      (L1*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
          L2*Qr2(i,j,k)*(-(Qt1(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
      (L1*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     a0.diff(d001,i,j,k)*((-4*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],4)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (pow(L2,2)*(1 - pow(1 - rgrid[i],2))*pow(Qr2(i,j,k),2)*
              Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              Qr1(i,j,k)*(L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                 Qt1(i,j,k) - Qtr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr2(i,j,k)*
              (2*L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
                Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (4*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
           (pow(L2,2)*(1 - pow(1 - rgrid[i],2))*pow(Qr2(i,j,k),2)*
              Qt1(i,j,k) + L1*
              pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              Qr1(i,j,k)*(L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                 Qt1(i,j,k) - Qtr(i,j,k)) + 
             (L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr2(i,j,k)*
                (2*L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                   Qt1(i,j,k) - Qtr(i,j,k))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L1*L2*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))))) + 
  a0.diff(d001,i,j,k)*((-4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],4)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
             Qtr(i,j,k))*(-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) \
+ Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
          L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
           (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
      (L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],4)*
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
             Qtr(i,j,k))*(-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                Qt2(i,j,k)) + 
             (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
          L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
           (Qt1(i,j,k) - ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     ar.diff(d001,i,j,k)*((-4*L1*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],4)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
                Qtr(i,j,k))*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) \
+ L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L2,2)*pow(rh,2)*
           (L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (4*L1*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
                Qtr(i,j,k)) + 
             (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (pow(L2,2)*pow(rh,2)*
           (L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a1.diff(d001,i,j,k)*((4*pow(1 + pow(R0,2),0.25)*
           pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (pow(L2,2)*(1 - pow(1 - rgrid[i],2))*pow(Qr2(i,j,k),2)*
              Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              Qr1(i,j,k)*(L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                 Qt1(i,j,k) - Qtr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr2(i,j,k)*
              (2*L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
                Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L2,2)*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (4*pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],4)*
           (pow(L2,2)*(1 - pow(1 - rgrid[i],2))*pow(Qr2(i,j,k),2)*
              Qt1(i,j,k) + L1*
              pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              Qr1(i,j,k)*(L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                 Qt1(i,j,k) - Qtr(i,j,k)) + 
             (L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr2(i,j,k)*
                (2*L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                   Qt1(i,j,k) - Qtr(i,j,k))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (pow(L2,2)*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))))) + 
  H*((-32*h2.diff(d001,i,j,k)*q*pow(rh,2)*h1(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*ygrid[k]*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (32*h1.diff(d001,i,j,k)*q*pow(rh,2)*h2(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*ygrid[k]*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     a1.diff(d100,i,j,k)*((-2*pow(1 + pow(R0,2),0.25)*
           pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (2*pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),3)*
           (1 - rgrid[i])*(L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              Qr1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     ar.diff(d010,i,j,k)*((4*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a2.diff(d100,i,j,k)*((2*pow(1 - pow(1 - rgrid[i],2),3)*
           (1 - rgrid[i])*(L2*Qr2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
        (2*pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     ar.diff(d001,i,j,k)*((-4*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(rh,2)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
        (4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(rh,2)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     a2.diff(d010,i,j,k)*((-4*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              pow(Qr1(i,j,k),2)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             pow(L2,2)*pow(Qr2(i,j,k),2)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              Qr1(i,j,k)*Qr2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (4*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
           (2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              Qr1(i,j,k)*Qr2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              pow(Qr1(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             pow(L2,2)*pow(Qr2(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a1.diff(d001,i,j,k)*((4*pow(1 + pow(R0,2),0.25)*
           pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              pow(Qr1(i,j,k),2)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             pow(L2,2)*pow(Qr2(i,j,k),2)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              Qr1(i,j,k)*Qr2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (4*pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],2)*
           (2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              Qr1(i,j,k)*Qr2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             pow(L1 + L1*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              pow(Qr1(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             pow(L2,2)*pow(Qr2(i,j,k),2)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     (4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*
        (L2*Qr2(i,j,k)*((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
      (pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
          L2*Qr2(i,j,k)*(-(Qt1(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
      (pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (4*L2*pow(q,2)*pow(rh,2)*
        (pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*ygrid[k]*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (4*a1(i,j,k)*(1 + 2*pow(R0,2) + 
             4*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             2*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
           sqrt(1 + pow(R0,2)) - 
          8*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
           (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (Qt1(i,j,k)*(1 + 2*pow(R0,2) + 
                4*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                2*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
      (pow(1 + pow(R0,2),0.25)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
     a0.diff(d100,i,j,k)*((-2*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],3)*
           (L2*Qr2(i,j,k)*((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                 Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) \
- Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
) - L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
                (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                 sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))\
))/(pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
        (2*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],3)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                 (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
                Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) \
+ L2*Qr2(i,j,k)*(-(Qt1(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
                (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                 sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))\
)/(pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     h2.diff(d010,i,j,k)*((16*L2*q*pow(rh,2)*h1(i,j,k)*
           pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           ygrid[k]*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
         (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (16*L2*q*pow(rh,2)*h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*ygrid[k]*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     h1.diff(d010,i,j,k)*((-16*L2*q*pow(rh,2)*h2(i,j,k)*
           pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           ygrid[k]*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
         (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (16*L2*q*pow(rh,2)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*ygrid[k]*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a0.diff(d001,i,j,k)*((-4*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],4)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (pow(L2,2)*(1 - pow(1 - rgrid[i],2))*pow(Qr2(i,j,k),2)*
              Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              Qr1(i,j,k)*(L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                 Qt1(i,j,k) - Qtr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr2(i,j,k)*
              (2*L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
                Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (4*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
           (pow(L2,2)*(1 - pow(1 - rgrid[i],2))*pow(Qr2(i,j,k),2)*
              Qt1(i,j,k) + L1*
              pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              Qr1(i,j,k)*(L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                 Qt1(i,j,k) - Qtr(i,j,k)) + 
             (L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr2(i,j,k)*
                (2*L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                   Qt1(i,j,k) - Qtr(i,j,k))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a0.diff(d010,i,j,k)*((4*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],4)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
               2)*(1 - pow(1 - rgrid[i],2))*pow(Qr1(i,j,k),2)*
              Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*(L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                 Qt2(i,j,k) - Qtr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (2*L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (4*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
           (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
               2)*(1 - pow(1 - rgrid[i],2))*pow(Qr1(i,j,k),2)*
              Qt2(i,j,k) + L2*Qr2(i,j,k)*
              (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                Qtr(i,j,k)) + 
             (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
                (2*L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                   Qt2(i,j,k) - Qtr(i,j,k))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))))) + 
  a0.diff(d100,i,j,k)*((2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],3)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (pow(Qt1(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(Qt2(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt1(i,j,k)*
           Qt2(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],3)*
        (-2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt1(i,j,k)*
           Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          pow(Qt1(i,j,k),2)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           pow(Qt2(i,j,k),2)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     a2.diff(d100,i,j,k)*((pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     ar.diff(d001,i,j,k)*((-2*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],3)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(rh,2)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (2*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*
           (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(rh,2)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a1.diff(d100,i,j,k)*(-((pow(1 + pow(R0,2),0.25)*
             pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
             ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
               Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
        (pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     ar.diff(d010,i,j,k)*((2*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],3)*
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
        (2*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (-(Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     a2.diff(d010,i,j,k)*((2*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],3)*
           (-(L2*Qr2(i,j,k)*((1 + 
                     B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
                  Qt1(i,j,k)*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)))) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
                (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                 sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))\
))/(L1*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
        (2*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],3)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                 (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
                Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) \
+ L2*Qr2(i,j,k)*(-(Qt1(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
                (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                 sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))\
)/(L1*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     a1.diff(d001,i,j,k)*((2*pow(1 + pow(R0,2),0.25)*
           pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],3)*
           (L2*Qr2(i,j,k)*((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                 Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) \
- Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
) - L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
                (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                 sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))\
))/(L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
        (2*pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],3)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                 (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
                Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) \
+ L2*Qr2(i,j,k)*(-(Qt1(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
                (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                 sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))\
)/(L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     a0.diff(d001,i,j,k)*((2*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],5)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
                Qtr(i,j,k))*(-((1 + 
                     B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
                Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
) + L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
              (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
                (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                 sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))\
))/(L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (2*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],5)*
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
                Qtr(i,j,k))*(-((1 + 
                     B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)) \
+ (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                 sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
) + L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
              (Qt1(i,j,k) - ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                   Qt2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                 sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))\
)*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a0.diff(d010,i,j,k)*((2*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],5)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*pow(Qt2(i,j,k),2)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             Qt1(i,j,k)*Qtr(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) + 
                Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)*
              (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
                (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                 sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))\
)/(L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (2*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],5)*
           (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*pow(Qt2(i,j,k),2) \
+ Qt1(i,j,k)*Qtr(i,j,k) - ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                Qt2(i,j,k)*(L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                   Qt1(i,j,k) + Qtr(i,j,k))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)*
              (Qt1(i,j,k) - ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                   Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))\
)/sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))))) + 
  a2.diff(d100,i,j,k)*((-4*ar.diff(d010,i,j,k)*
        (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
        (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     a2.diff(d010,i,j,k)*((2*pow(1 - pow(1 - rgrid[i],2),3)*
           (1 - rgrid[i])*(L2*Qr2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (sqrt(1 + pow(R0,2))*pow(rh,2)*
           (L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
        (2*pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (sqrt(1 + pow(R0,2))*pow(rh,2)*
           (L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     a1.diff(d001,i,j,k)*((-2*pow(1 - pow(1 - rgrid[i],2),3)*
           (1 - rgrid[i])*(L2*Qr2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(rh,2)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
        (2*pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(rh,2)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     ar.diff(d001,i,j,k)*((2*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
           pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
         (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a0.diff(d001,i,j,k)*((2*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],3)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
                Qtr(i,j,k))*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) \
+ L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (2*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],3)*
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
                Qtr(i,j,k)) + 
             (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a0.diff(d010,i,j,k)*((-2*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],3)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (2*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],3)*
           (L2*R0*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
             R0*Qtr(i,j,k) + L2*pow(1 - pow(1 - rgrid[i],2),2)*
              Qr2(i,j,k)*Qt2(i,j,k)*R(i,j,k) - 
             (1 - pow(1 - rgrid[i],2))*Qtr(i,j,k)*R(i,j,k) + 
             L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt2(i,j,k)*
              sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) \
+ L1*B(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*Qr1(i,j,k)*Qt2(i,j,k)*
              sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
)/(pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))))) + 
  a1.diff(d100,i,j,k)*((2*a2.diff(d100,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (4*ar.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
        (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L2*pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (2*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (2*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     a2.diff(d010,i,j,k)*((-2*pow(1 - pow(1 - rgrid[i],2),3)*
           (1 - rgrid[i])*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (2*pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a1.diff(d001,i,j,k)*((2*pow(1 - pow(1 - rgrid[i],2),3)*
           (1 - rgrid[i])*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2))*
           (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (2*pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*
           sqrt(1 + pow(R0,2))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     ar.diff(d010,i,j,k)*((2*pow(1 + pow(R0,2),0.25)*
           (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
           pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))/
         (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (2*pow(1 + pow(R0,2),0.25)*(1 - pow(1 - rgrid[i],2))*
           (1 - rgrid[i])*sqrt(1 + 
             pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a0.diff(d001,i,j,k)*((-2*pow(1 + pow(R0,2),0.25)*
           pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],3)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
                Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (2*pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],3)*
           (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt1(i,j,k) + 
             ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                   Qt1(i,j,k) - Qtr(i,j,k))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a0.diff(d010,i,j,k)*((2*pow(1 + pow(R0,2),0.25)*
           pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],3)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             Qtr(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (2*pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],3)*
           (L1*R0*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt2(i,j,k) + 
             L1*R0*B(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*Qr1(i,j,k)*
              Qt2(i,j,k) + L1*pow(1 - pow(1 - rgrid[i],2),2)*
              Qr1(i,j,k)*Qt2(i,j,k)*R(i,j,k) + 
             L1*B(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*Qr1(i,j,k)*
              Qt2(i,j,k)*R(i,j,k) + 
             L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)*
              sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) \
- Qtr(i,j,k)*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
         (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))))) + 
  a0.diff(d010,i,j,k)*((-4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],4)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*pow(Qt2(i,j,k),2)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          Qt1(i,j,k)*Qtr(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) + 
             Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)*
           (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
      (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
        (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
           (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*pow(Qt2(i,j,k),2) + 
          Qt1(i,j,k)*Qtr(i,j,k) - 
          ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
             (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) + 
               Qtr(i,j,k))*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)*
           (Qt1(i,j,k) - ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     a2.diff(d010,i,j,k)*((4*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],4)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
               2)*(1 - pow(1 - rgrid[i],2))*pow(Qr1(i,j,k),2)*
              Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*(L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                 Qt2(i,j,k) - Qtr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (2*L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L1,2)*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (4*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
           (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
               2)*(1 - pow(1 - rgrid[i],2))*pow(Qr1(i,j,k),2)*
              Qt2(i,j,k) + L2*Qr2(i,j,k)*
              (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                Qtr(i,j,k)) + 
             (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
                (2*L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                   Qt2(i,j,k) - Qtr(i,j,k))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (pow(L1,2)*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a1.diff(d001,i,j,k)*((-4*pow(1 + pow(R0,2),0.25)*
           pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
               2)*(1 - pow(1 - rgrid[i],2))*pow(Qr1(i,j,k),2)*
              Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L2*Qr2(i,j,k)*(L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                 Qt2(i,j,k) - Qtr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
              (2*L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (4*pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],4)*
           (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
               2)*(1 - pow(1 - rgrid[i],2))*pow(Qr1(i,j,k),2)*
              Qt2(i,j,k) + L2*Qr2(i,j,k)*
              (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                Qtr(i,j,k)) + 
             (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
                (2*L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                   Qt2(i,j,k) - Qtr(i,j,k))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L1*L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     a0.diff(d001,i,j,k)*((64*pow(rh,4)*
           pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],6)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           ((L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(1 - pow(1 - rgrid[i],2),3)*
                pow(-1 + rgrid[i],4)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                Qr1(i,j,k)*Qt2(i,j,k)*
                (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                   Qt1(i,j,k) - Qtr(i,j,k))*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              (16.*pow(rh,6)) + 
             (L2*pow(1 - pow(1 - rgrid[i],2),3)*
                pow(-1 + rgrid[i],4)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                Qr2(i,j,k)*Qt1(i,j,k)*
                (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                   Qt2(i,j,k) - Qtr(i,j,k))*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              (16.*pow(rh,6)) + 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (-1 - (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
                (pow(1 - pow(1 - rgrid[i],2),2)*
                   pow(-1 + rgrid[i],4)*
                   pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                   (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                      Qr1(i,j,k)*Qr2(i,j,k)*Qt1(i,j,k)*Qt2(i,j,k) - 
                     (1 - pow(1 - rgrid[i],2))*
                      (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                        L2*Qr2(i,j,k)*Qt2(i,j,k))*Qtr(i,j,k) + 
                     pow(Qtr(i,j,k),2))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
                 (16.*pow(rh,6)))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],4)*
           pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (64*pow(rh,4)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],6)*
           ((L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(1 - pow(1 - rgrid[i],2),3)*
                pow(-1 + rgrid[i],4)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                Qr1(i,j,k)*Qt2(i,j,k)*
                (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
                   Qt1(i,j,k) - Qtr(i,j,k))*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
              (16.*pow(rh,6)) + 
             (L2*pow(1 - pow(1 - rgrid[i],2),3)*
                pow(-1 + rgrid[i],4)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                Qr2(i,j,k)*Qt1(i,j,k)*
                (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
                   Qt2(i,j,k) - Qtr(i,j,k))*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
              (16.*pow(rh,6)) + 
             ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (-1 - (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) + 
                  (pow(1 - pow(1 - rgrid[i],2),2)*
                     pow(-1 + rgrid[i],4)*
                     pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                       pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                     (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                        Qr1(i,j,k)*Qr2(i,j,k)*Qt1(i,j,k)*Qt2(i,j,k) - 
                       (1 - pow(1 - rgrid[i],2))*
                        (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                       L2*Qr2(i,j,k)*Qt2(i,j,k))*Qtr(i,j,k) + 
                       pow(Qtr(i,j,k),2))*
                     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
                   (16.*pow(rh,6)))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],4)*
           pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))) + 
     ar.diff(d001,i,j,k)*((4*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],4)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (L2*pow(rh,2)*(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
        (4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
           (L2*R0*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
             R0*Qtr(i,j,k) + L2*pow(1 - pow(1 - rgrid[i],2),2)*
              Qr2(i,j,k)*Qt2(i,j,k)*R(i,j,k) - 
             (1 - pow(1 - rgrid[i],2))*Qtr(i,j,k)*R(i,j,k) + 
             L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt2(i,j,k)*
              sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) \
+ L1*B(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*Qr1(i,j,k)*Qt2(i,j,k)*
              sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
)/(L1*L2*pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     ar.diff(d010,i,j,k)*((-4*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],4)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             Qtr(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (pow(L1,2)*pow(rh,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
        (4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
           (L1*R0*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt2(i,j,k) + 
             L1*R0*B(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*Qr1(i,j,k)*
              Qt2(i,j,k) + L1*pow(1 - pow(1 - rgrid[i],2),2)*
              Qr1(i,j,k)*Qt2(i,j,k)*R(i,j,k) + 
             L1*B(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*Qr1(i,j,k)*
              Qt2(i,j,k)*R(i,j,k) + 
             L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)*
              sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) \
- Qtr(i,j,k)*sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
         (pow(L1,2)*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))
;

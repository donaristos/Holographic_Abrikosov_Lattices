
elem[103]+=
(Qtr.diff(d100,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
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
           (2.*pow(rh,4))))*Qr1(i,j,k)*Qr2(i,j,k))/(4.*pow(rh,2)) + 
  (16*a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   (pow(L2,2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
  (8*a0.diff(d001,i,j,k)*a1.diff(d100,i,j,k)*
     pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),2)*
     pow(1 - rgrid[i],3)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   (pow(L2,2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  (16*a0.diff(d001,i,j,k)*ar.diff(d010,i,j,k)*
     (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],4)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   (L1*pow(L2,2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  a2.diff(d010,i,j,k)*((-16*a0(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (32*a0.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],4)*Qr2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  (Qrr.diff(d100,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),3)*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
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
           (2.*pow(rh,4))))*Qr1(i,j,k)*Qr2(i,j,k)*
     ((1 - pow(1 - rgrid[i],2))*
        (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k))\
)/(4.*pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
  (L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
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
           (2.*pow(rh,4))))*Qr1(i,j,k)*Qr2(i,j,k)*Qrr(i,j,k)*
     ((1 - pow(1 - rgrid[i],2))*
        (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k))\
)/(2.*pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) - 
  (Qr1.diff(d100,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
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
           (2.*pow(rh,4))))*Qr2(i,j,k)*
     ((1 - pow(1 - rgrid[i],2))*
        (2*L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
       Qtr(i,j,k)))/(4.*pow(rh,2)) + 
  (L1*(1 - pow(1 - rgrid[i],2))*
     (-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
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
           (2.*pow(rh,4))))*Qr1(i,j,k)*Qr2(i,j,k)*
     ((1 - pow(1 - rgrid[i],2))*
        (2*L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
       Qtr(i,j,k)))/(8.*pow(rh,2)*(1 - rgrid[i])) + 
  (L1*(1 - pow(1 - rgrid[i],2))*
     (-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
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
           (2.*pow(rh,4))))*Qr1(i,j,k)*Qr2(i,j,k)*Qtr(i,j,k))/
   (8.*pow(rh,2)*(1 - rgrid[i])) + 
  (8*h2.diff(d100,i,j,k)*q*pow(rh,2)*a0(i,j,k)*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
  (8*h1.diff(d100,i,j,k)*q*pow(rh,2)*a0(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  (Qt2.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2))*
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
           (2.*pow(rh,4))))*Qr2(i,j,k)*
     ((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          Qr2(i,j,k))/(4.*pow(rh,2)) - 
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))/
   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
  (2*(1 - rgrid[i])*(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
           (2.*pow(rh,4))))*Qr2(i,j,k)*Qt2(i,j,k)*
     ((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          Qr2(i,j,k))/(4.*pow(rh,2)) - 
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))/
   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
  (2*Qtr.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
     ((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
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
   (L2*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  (2*Qr2.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
     (1 - rgrid[i])*(-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
           (2.*pow(rh,4))))*Qr1(i,j,k)*Qt2(i,j,k)*
     ((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
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
   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
  (2*Qrr.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
     (1 - rgrid[i])*(-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
     ((1 - pow(1 - rgrid[i],2))*
        (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)\
)*((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
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
   (L2*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
  (2*Qr1.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
           (2.*pow(rh,4))))*((1 - pow(1 - rgrid[i],2))*
        (2*L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
       Qtr(i,j,k))*((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
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
   (L2*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  H*((-16*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))\
)/(L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (8*a0.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],3)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))\
)/(L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (32*a0.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],4)*Qr2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
      (L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (32*L1*pow(q,2)*pow(rh,2)*a0(i,j,k)*
        (pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*ygrid[k]*
        Qr1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  Qt2.diff(d010,i,j,k)*((4*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
      (L1*L2*pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (8*pow(rh,4)*(1 - rgrid[i])*
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
              (2.*pow(rh,4))))*
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
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
        ((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
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
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) - 
  (32*pow(a0.diff(d001,i,j,k),2)*pow(1 - pow(1 - rgrid[i],2),4)*
     pow(1 - rgrid[i],6)*Qr2(i,j,k)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
   (pow(L2,2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
  (32*h2.diff(d010,i,j,k)*q*pow(rh,2)*a0(i,j,k)*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*Qr1(i,j,k)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (32*h1.diff(d010,i,j,k)*q*pow(rh,2)*a0(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*Qr1(i,j,k)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
  (16*h2.diff(d001,i,j,k)*q*pow(rh,2)*a0(i,j,k)*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
       L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   (pow(L2,2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  (16*h1.diff(d001,i,j,k)*q*pow(rh,2)*a0(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
       L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   (pow(L2,2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
  (Qt1.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2))*
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
           (2.*pow(rh,4))))*Qr2(i,j,k)*
     ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
       (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr1(i,j,k),2)*
          sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        (4.*pow(rh,2))))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
  (2*(1 - rgrid[i])*(-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
           (2.*pow(rh,4))))*Qr2(i,j,k)*Qt1(i,j,k)*
     ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
       (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr1(i,j,k),2)*
          sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        (4.*pow(rh,2))))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (8*Qt1.diff(d010,i,j,k)*pow(rh,2)*(1 - rgrid[i])*
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
           (2.*pow(rh,4))))*((L1*L2*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          Qr2(i,j,k))/(4.*pow(rh,2)) + 
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
     ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
       (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr1(i,j,k),2)*
          sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        (4.*pow(rh,2))))/
   (L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
  (2*Qr2.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
     (1 - rgrid[i])*(-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
           (2.*pow(rh,4))))*Qr1(i,j,k)*Qt2(i,j,k)*
     (-(pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr2(i,j,k),2)\
)/(4.*pow(rh,2)) + (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   (L2*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  (2*Qrr.diff(d001,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*
     (1 - rgrid[i])*(-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
     ((1 - pow(1 - rgrid[i],2))*
        (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)\
)*(-(pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr2(i,j,k),2)\
)/(4.*pow(rh,2)) + (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   (pow(L2,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) - 
  (R.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
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
           (2.*pow(rh,4))))*Qr2(i,j,k)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (8*R.diff(d010,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     (1 - rgrid[i])*(-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
           (2.*pow(rh,4))))*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     ((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          Qr2(i,j,k))/(4.*pow(rh,2)) + 
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
     (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   (L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
           (2.*pow(rh,4))))*Qr2(i,j,k)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*R(i,j,k)*
     (-(Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
           (2.*pow(rh,4))))*Q(i,j,k)*Qr2(i,j,k)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (Qt1(i,j,k)*(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
  (2*(1 - rgrid[i])*(-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
           (2.*pow(rh,4))))*Qr2(i,j,k)*
     ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
        (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (-(L1*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
              (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   4*pow(rh,2)*
                    (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                      4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
              ((1 - pow(1 - rgrid[i],2))*
                 (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) \
- Qtr(i,j,k)))/(4.*pow(rh,2)) + 
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (Qr2.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2))*
     ((((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
              (2.*pow(rh,4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
        (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
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
        (-(L1*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
              (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   4*pow(rh,2)*
                    (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                      4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
              ((1 - pow(1 - rgrid[i],2))*
                 (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                   2*L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)))/
           (4.*pow(rh,2)) + 
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  ((-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*Qr2(i,j,k)*
     ((((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
              (2.*pow(rh,4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
        (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
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
        (-(L1*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
              (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   4*pow(rh,2)*
                    (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                      4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
              ((1 - pow(1 - rgrid[i],2))*
                 (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                   2*L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)))/
           (4.*pow(rh,2)) + 
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   (2.*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*(1 - rgrid[i])*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  a0.diff(d001,i,j,k)*((32*a1.diff(d001,i,j,k)*
        pow(1 + pow(R0,2),0.25)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],4)*Qr2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (32*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],4)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  a0.diff(d100,i,j,k)*((8*a2.diff(d010,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
      (L1*L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (8*a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (8*a0.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],5)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (16*a0.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],5)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  a0.diff(d010,i,j,k)*((16*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (16*a0.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],6)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (-2*L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) + 
          Qtr(i,j,k))*sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
      (L1*pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  (pow(rh,2)*(-16*pow(L2,2)*pow(rh,2)*
        (-((pow(1 - rgrid[i],3)*pow(-1 + rgrid[i],4)*
               pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))),2))/
             pow(rh,8)) + ((1 - pow(1 - rgrid[i],2))*
             (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (-(pow(-1 + rgrid[i],2)*
                   (pow(H,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     pow(rh,2)*
                      (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
                (2.*pow(rh,4)) - 
               2*(1 - rgrid[i])*
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
                   (2.*pow(rh,4)))))/pow(rh,4) + 
          pow(1 - pow(1 - rgrid[i],2),2)*
           (((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                ((pow(-1 + rgrid[i],2)*
                     (6*pow(H,2)*pow(-2 + rgrid[i],3)*rgrid[i] + 
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
                   (2.*pow(rh,4))))/(2.*pow(rh,4)) + 
             (pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
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
                     (pow(H,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
                   (2.*pow(rh,4))))/pow(rh,4) + 
             2*(1 - rgrid[i])*
              pow((pow(-1 + rgrid[i],2)*
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
                 (2.*pow(rh,4)),2)))*Qr2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
        (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
       (32*L1*a0(i,j,k)*(pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
          pow(q + q*B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],5)*
          Qr1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          (-4*L2*a2(i,j,k) + 
            4*L2*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
             (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
             Qt2(i,j,k))*(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        pow(1 + pow(R0,2),0.25) + 
       L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        ((L1*L2*(1 - pow(1 - rgrid[i],2))*
             (6*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
             (1 - rgrid[i])*pow(-1 + rgrid[i],6)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),3)*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qr2(i,j,k)*((1 - pow(1 - rgrid[i],2))*
                (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
               Qtr(i,j,k)))/pow(rh,8) + 
          (L2*pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             Qr2(i,j,k)*(16*L1*pow(rh,2)*
                pow(1 - pow(1 - rgrid[i],2),2)*
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
                     (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                       pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
                   (2.*pow(rh,4)))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                ((1 - pow(1 - rgrid[i],2))*
                   (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                     L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)) - 
               L1*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),3)*
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
                      (2.*pow(rh,4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                ((1 - pow(1 - rgrid[i],2))*
                   (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                     L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)) - 
               8*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) \
- 16*pow(1 - rgrid[i],3)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))/
           pow(rh,6) - (4*L2*(1 - pow(1 - rgrid[i],2))*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qr2(i,j,k)*(4*L1*pow(rh,2)*
                pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                pow((pow(-1 + rgrid[i],2)*
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
                   (2.*pow(rh,4)),2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                ((1 - pow(1 - rgrid[i],2))*
                   (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                     L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)) + 
               8*pow(1 - rgrid[i],2)*
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
                Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) \
- (1 - pow(1 - rgrid[i],2))*
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
                  4*((pow(-1 + rgrid[i],2)*
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
                      (2.*pow(rh,4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))/pow(rh,2) \
+ (2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             (64*pow(q,2)*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
                ar(i,j,k)*(pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
                pow(1 - rgrid[i],4) + 
               16*L2*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
                (1 - pow(1 - rgrid[i],2))*
                pow((pow(-1 + rgrid[i],2)*
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
                   (2.*pow(rh,4)),2)*Qr2(i,j,k)*Qt2(i,j,k) + 
               64*pow(q,2)*pow(1 + pow(R0,2),0.25)*
                pow(a0(i,j,k),2)*
                (pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
                pow(1 - pow(1 - rgrid[i],2),3)*
                pow(1 - rgrid[i],6)*
                (2*L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) \
- 16*pow(q,2)*a0(i,j,k)*(pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
                pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
                (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
                  4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
                   pow(1 - rgrid[i],2)*Qtr(i,j,k) + 
                  8*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))))/
           pow(1 + pow(R0,2),0.25))*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   (8.*pow(L2,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (8*Qt2.diff(d001,i,j,k)*pow(rh,2)*(1 - rgrid[i])*
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
           (2.*pow(rh,4))))*((L1*pow(L2,3)*
          pow(1 - pow(1 - rgrid[i],2),4)*pow(-1 + rgrid[i],4)*
          pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*Qr1(i,j,k)*
          pow(Qr2(i,j,k),3))/(16.*pow(rh,4)) + 
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k),2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
       (L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          (L2*Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
            L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qr1(i,j,k)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (4.*pow(rh,2))))/
   (pow(L2,2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  (Q.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*Qr2(i,j,k)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     ((1 - pow(1 - rgrid[i],2))*
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
           (2.*pow(rh,4)))*(Qt1(i,j,k)*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
       ((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (-(Qt1(i,j,k)*(1 + pow(R0,2) + 
                 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                 pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
            (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        pow(rh,4)))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
  (8*Q.diff(d010,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     (1 - rgrid[i])*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     ((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          Qr2(i,j,k))/(4.*pow(rh,2)) + 
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
     ((1 - pow(1 - rgrid[i],2))*
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
           (2.*pow(rh,4)))*(Qt1(i,j,k)*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
       ((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (-(Qt1(i,j,k)*(1 + pow(R0,2) + 
                 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                 pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
            (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        pow(rh,4)))/
   (L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  Qt1.diff(d001,i,j,k)*((-4*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
      (pow(L2,2)*pow(rh,2)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (8*pow(rh,4)*(1 - rgrid[i])*
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
        (-(pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
              (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   4*pow(rh,2)*
                    (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                      4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
              pow(Qr1(i,j,k),2)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (4.*pow(rh,2)) + (pow(L2,2)*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (4.*pow(rh,2)) - (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
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
                 pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
                 pow(Qr1(i,j,k),2)*pow(Qr2(i,j,k),2))/
              (16.*pow(rh,4)) + 
             pow(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k),2)*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(L2,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  (8*R.diff(d001,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     (1 - rgrid[i])*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (-(pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
           pow(-1 + rgrid[i],4)*
           pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           pow(Qr2(i,j,k),2)*
           (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (4.*pow(rh,6)) - (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*
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
           (2.*pow(rh,4)))*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
          Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
       (pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (pow(L2,2)*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),3)*
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
             pow(Qr2(i,j,k),2)*
             (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
               (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) \
+ 4*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*(1 - rgrid[i])*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             (-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                  (1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
) + Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
        (4.*pow(rh,4))))/
   (pow(L2,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
  (8*Q.diff(d001,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     (1 - rgrid[i])*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
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
             (2.*pow(rh,4)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
          (-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
               (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
            Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))) - 
       (pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
          pow(-1 + rgrid[i],4)*
          pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr2(i,j,k),2)*
          (Qt1(i,j,k)*(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (4.*pow(rh,6)) + (pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (4*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*(1 - rgrid[i])*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
             (-((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                  (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
               Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) \
+ pow(L2,2)*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),3)*
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
             pow(Qr2(i,j,k),2)*
             (Qt1(i,j,k)*(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
               (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
        (4.*pow(rh,4))))/
   (pow(L2,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),3)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (B.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
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
           (2.*pow(rh,4))))*Qr2(i,j,k)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  (2*B(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
           (2.*pow(rh,4))))*Qr2(i,j,k)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
  (8*B.diff(d010,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     (1 - rgrid[i])*(-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
           (2.*pow(rh,4))))*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     Qt1(i,j,k)*((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
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
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (L1*L2*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  (8*B.diff(d001,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     (1 - rgrid[i])*(-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
           (2.*pow(rh,4))))*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     Qt1(i,j,k)*(-(pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
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
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (pow(L2,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  (2*Qtr.diff(d001,i,j,k)*L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
     ((pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr2(i,j,k),2))/
        (4.*pow(rh,2)) - (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
   (pow(L2,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  (2*Qr1.diff(d001,i,j,k)*L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
           (2.*pow(rh,4))))*((1 - pow(1 - rgrid[i],2))*
        (2*L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
       Qtr(i,j,k))*((pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr2(i,j,k),2))/
        (4.*pow(rh,2)) - (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
   (pow(L2,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))
;

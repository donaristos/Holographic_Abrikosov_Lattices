
elem[239]+=
(q*ar(i,j,k)*(2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2)))/
   pow(1 - pow(1 - rgrid[i],2),2) - 
  (2*ar.diff(d100,i,j,k)*q*(1 - rgrid[i]))/(1 - pow(1 - rgrid[i],2)) + 
  2*Qtr.diff(d100,i,j,k)*q*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
   pow(1 - rgrid[i],3) + (4*q*pow(rh,4)*ar(i,j,k)*
     ((2*pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/pow(rh,4) \
+ (1 - pow(1 - rgrid[i],2))*(-(pow(-1 + rgrid[i],2)*
              (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   4*pow(rh,2)*
                    (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                      4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           (2.*pow(rh,4)) - 2*(1 - rgrid[i])*
           ((pow(-1 + rgrid[i],2)*
                (3*pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],2) + 
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
              (2.*pow(rh,4))))))/
   (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
  2*a1.diff(d100,i,j,k)*L1*q*pow(1 + pow(R0,2),0.25)*
   (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) + 
  4*ar.diff(d010,i,j,k)*q*pow(1 - rgrid[i],2)*Qr1(i,j,k) - 
  2*Qt1.diff(d100,i,j,k)*L1*q*a0(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*Qr1(i,j,k) - 
  4*Qtr.diff(d010,i,j,k)*q*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],4)*Qr1(i,j,k) + 
  (4*L1*q*pow(1 + pow(R0,2),0.25)*pow(rh,4)*a1(i,j,k)*
     ((-2*pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/pow(rh,4) \
+ (1 - pow(1 - rgrid[i],2))*
        ((pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
           (2.*pow(rh,4)) + 
          2*(1 - rgrid[i])*((pow(-1 + rgrid[i],2)*
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
              (2.*pow(rh,4)))))*Qr1(i,j,k))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
  (2*L2*q*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k))/
   pow(1 + pow(R0,2),0.25) + 
  (2*a2.diff(d100,i,j,k)*L2*q*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
     Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
  4*ar.diff(d001,i,j,k)*q*pow(1 - rgrid[i],2)*Qr2(i,j,k) - 
  (8*L2*q*a2(i,j,k)*pow(1 - rgrid[i],2)*Qr2(i,j,k))/
   pow(1 + pow(R0,2),0.25) - 
  2*Qt2.diff(d100,i,j,k)*L2*q*a0(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*Qr2(i,j,k) - 
  4*Qtr.diff(d001,i,j,k)*q*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],4)*Qr2(i,j,k) + 
  (8*L2*q*pow(rh,4)*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
        (2.*pow(rh,4)))*Qr2(i,j,k))/
   (pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
  (768*pow(rh,4)*S*h1(i,j,k)*h2(i,j,k)*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) - 
  2*L1*q*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
   Qr1(i,j,k)*Qt1(i,j,k) + 4*L1*q*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
   pow(1 - rgrid[i],4)*Qr1(i,j,k)*Qt1(i,j,k) - 
  (8*L1*q*pow(rh,4)*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
     pow(1 - rgrid[i],3)*((pow(-1 + rgrid[i],2)*
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
        (2.*pow(rh,4)))*Qr1(i,j,k)*Qt1(i,j,k))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
  2*Qr1.diff(d100,i,j,k)*q*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
   (L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k) - 
     L1*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
      Qt1(i,j,k)) + q*(-2*(1 - pow(1 - rgrid[i],2)) + 
     4*pow(1 - rgrid[i],2))*Qr1(i,j,k)*
   (L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k) - 
     L1*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
      Qt1(i,j,k)) + 4*Qr1.diff(d001,i,j,k)*q*
   pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*Qr2(i,j,k)*
   (-(L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)) + 
     L1*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
      Qt1(i,j,k)) - 2*L2*q*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],2)*Qr2(i,j,k)*Qt2(i,j,k) + 
  4*L2*q*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],4)*
   Qr2(i,j,k)*Qt2(i,j,k) - (8*L2*q*pow(rh,4)*a0(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*
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
        (2.*pow(rh,4)))*Qr2(i,j,k)*Qt2(i,j,k))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
  2*Qr2.diff(d100,i,j,k)*q*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
   ((L2*a2(i,j,k))/pow(1 + pow(R0,2),0.25) - 
     L2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
      Qt2(i,j,k)) + q*(-2*(1 - pow(1 - rgrid[i],2)) + 
     4*pow(1 - rgrid[i],2))*Qr2(i,j,k)*
   ((L2*a2(i,j,k))/pow(1 + pow(R0,2),0.25) - 
     L2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
      Qt2(i,j,k)) + 4*Qr2.diff(d010,i,j,k)*q*
   pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*Qr1(i,j,k)*
   (-((L2*a2(i,j,k))/pow(1 + pow(R0,2),0.25)) + 
     L2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
      Qt2(i,j,k)) + (2*Q.diff(d010,i,j,k)*q*
     pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*Qr1(i,j,k)*
     ((2*ar(i,j,k)*(1 - rgrid[i]))/(1 - pow(1 - rgrid[i],2)) - 
       2*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) - 
       (2*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
          Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
       2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],3)*
        ((1 - pow(1 - rgrid[i],2))*
           (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
          Qtr(i,j,k))))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
  (2*Q.diff(d001,i,j,k)*q*pow(1 - pow(1 - rgrid[i],2),2)*
     (1 - rgrid[i])*Qr2(i,j,k)*
     ((2*ar(i,j,k)*(1 - rgrid[i]))/(1 - pow(1 - rgrid[i],2)) - 
       2*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) - 
       (2*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
          Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
       2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],3)*
        ((1 - pow(1 - rgrid[i],2))*
           (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
          Qtr(i,j,k))))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
  (Qrr.diff(d100,i,j,k)*q*(1 - pow(1 - rgrid[i],2))*
     ((2*ar(i,j,k)*(1 - rgrid[i]))/(1 - pow(1 - rgrid[i],2)) - 
       2*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) - 
       (2*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
          Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
       2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],3)*
        ((1 - pow(1 - rgrid[i],2))*
           (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
          Qtr(i,j,k))))/(2.*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
  2*Qr1.diff(d010,i,j,k)*q*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
   ((2*ar(i,j,k)*(1 - rgrid[i]))/(1 - pow(1 - rgrid[i],2)) - 
     4*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
      (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) - 
     (2*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
     2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],3)*
      ((1 - pow(1 - rgrid[i],2))*
         (2*L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
        Qtr(i,j,k))) + 2*Qr2.diff(d001,i,j,k)*q*
   (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
   ((2*ar(i,j,k)*(1 - rgrid[i]))/(1 - pow(1 - rgrid[i],2)) - 
     2*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
      (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) - 
     (4*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
     2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],3)*
      ((1 - pow(1 - rgrid[i],2))*
         (L1*Qr1(i,j,k)*Qt1(i,j,k) + 2*L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
        Qtr(i,j,k))) + 4*q*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
   pow(1 - rgrid[i],2)*((1 - pow(1 - rgrid[i],2))*
      (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)) + 
  2*q*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
   Qtr(i,j,k) + q*a0(i,j,k)*(-2*(1 - pow(1 - rgrid[i],2)) + 
     4*pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qtr(i,j,k) - 
  8*q*a0(i,j,k)*pow(1 - rgrid[i],4)*Qtr(i,j,k) + 
  (8*q*pow(rh,4)*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],3)*((pow(-1 + rgrid[i],2)*
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
        (2.*pow(rh,4)))*Qtr(i,j,k))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
  2*a0.diff(d100,i,j,k)*q*(1 - pow(1 - rgrid[i],2))*
   pow(1 - rgrid[i],3)*(-((1 - pow(1 - rgrid[i],2))*
        (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))) + Qtr(i,j,k)) \
+ 2*q*pow(1 - rgrid[i],2)*((-4*ar(i,j,k))/
      pow(1 - pow(1 - rgrid[i],2),2) + 
     4*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*Qr1(i,j,k) + 
     (4*L2*a2(i,j,k)*Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
     4*a0(i,j,k)*pow(1 - rgrid[i],2)*
      (-((1 - pow(1 - rgrid[i],2))*
           (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))) + 
        Qtr(i,j,k))) + (Q.diff(d100,i,j,k)*q*(1 - pow(1 - rgrid[i],2))*
     ((-2*ar(i,j,k)*(1 - rgrid[i]))/(1 - pow(1 - rgrid[i],2)) + 
       2*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) + 
       (2*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
          Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
       2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],3)*
        (-((1 - pow(1 - rgrid[i],2))*
             (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))) + 
          Qtr(i,j,k))))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
  (Qtt.diff(d100,i,j,k)*q*(1 - pow(1 - rgrid[i],2))*
     ((2*ar(i,j,k)*(1 - rgrid[i]))/(1 - pow(1 - rgrid[i],2)) - 
       2*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) - 
       (2*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
          Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
       2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],3)*
        ((1 - pow(1 - rgrid[i],2))*
           (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
          Qtr(i,j,k))))/(2.*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
  4*Qt1.diff(d001,i,j,k)*L1*q*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
   pow(1 - rgrid[i],4)*(pow(1 - pow(1 - rgrid[i],2),2)*Qr1(i,j,k)*
      Qr2(i,j,k) - (4*pow(rh,2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  4*Qt2.diff(d010,i,j,k)*L2*q*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
   pow(1 - rgrid[i],4)*(pow(1 - pow(1 - rgrid[i],2),2)*Qr1(i,j,k)*
      Qr2(i,j,k) - (4*pow(rh,2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  (4*a2.diff(d010,i,j,k)*L2*q*pow(1 - rgrid[i],2)*
     (-(pow(1 - pow(1 - rgrid[i],2),2)*Qr1(i,j,k)*Qr2(i,j,k)) + 
       (4*pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
        (L1*L2*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/
   pow(1 + pow(R0,2),0.25) + 
  4*a1.diff(d001,i,j,k)*L1*q*pow(1 + pow(R0,2),0.25)*
   pow(1 - rgrid[i],2)*(-(pow(1 - pow(1 - rgrid[i],2),2)*Qr1(i,j,k)*
        Qr2(i,j,k)) + (4*pow(rh,2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  (q*pow(1 - rgrid[i],2)*Qrr(i,j,k)*
     (4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*
        (-4*L2*a2(i,j,k)*Qr2(i,j,k) + 
          4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*pow(1 - rgrid[i],2)*
           (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) + 
             L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
             Qtr(i,j,k)) - 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))\
))/(2.*pow(1 + pow(R0,2),0.25)*(1 - pow(1 - rgrid[i],2))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
  (q*pow(1 - rgrid[i],2)*Q(i,j,k)*
     (-4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*
        (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
          4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*pow(1 - rgrid[i],2)*
           (-((1 - pow(1 - rgrid[i],2))*
                (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))) + 
             Qtr(i,j,k)) + 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))\
))/(pow(1 + pow(R0,2),0.25)*(1 - pow(1 - rgrid[i],2))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  (q*pow(1 - rgrid[i],2)*Qtt(i,j,k)*
     (-4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*
        (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
          4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*pow(1 - rgrid[i],2)*
           (-((1 - pow(1 - rgrid[i],2))*
                (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))) + 
             Qtr(i,j,k)) + 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))\
))/(2.*pow(1 + pow(R0,2),0.25)*(1 - pow(1 - rgrid[i],2))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
  (16*a1.diff(d010,i,j,k)*q*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
     pow(1 - rgrid[i],2)*((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
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
   (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (16*Qt1.diff(d010,i,j,k)*q*pow(rh,2)*a0(i,j,k)*
     (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],4)*
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
   (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
  (2*Qrr.diff(d001,i,j,k)*q*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*((L2*ar(i,j,k)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k))/
        pow(rh,2) - (pow(L2,2)*a2(i,j,k)*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr2(i,j,k),2))/
        (pow(1 + pow(R0,2),0.25)*pow(rh,2)) + 
       (L1*L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
          pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          Qr2(i,j,k)*Qt1(i,j,k))/pow(rh,2) + 
       4*R0*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) + 
       (pow(L2,2)*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
          pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr2(i,j,k),2)*
          Qt2(i,j,k))/pow(rh,2) - 
       (L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
          pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
          Qtr(i,j,k))/pow(rh,2) + 
       4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*R(i,j,k) - 
       4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
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
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
       (4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        pow(1 + pow(R0,2),0.25) - 
       4*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
  (2*Qtt.diff(d001,i,j,k)*q*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*((L2*ar(i,j,k)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k))/
        pow(rh,2) - (pow(L2,2)*a2(i,j,k)*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr2(i,j,k),2))/
        (pow(1 + pow(R0,2),0.25)*pow(rh,2)) + 
       (L1*L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
          pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          Qr2(i,j,k)*Qt1(i,j,k))/pow(rh,2) - 
       4*R0*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) + 
       (pow(L2,2)*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
          pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr2(i,j,k),2)*
          Qt2(i,j,k))/pow(rh,2) - 
       (L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
          pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
          Qtr(i,j,k))/pow(rh,2) - 
       4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*R(i,j,k) + 
       4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
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
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) - 
       (4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        pow(1 + pow(R0,2),0.25) + 
       4*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
  (16*a0.diff(d010,i,j,k)*q*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],4)*((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        Qt1(i,j,k)*(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
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
   (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (4*R.diff(d001,i,j,k)*q*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     ((-4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
        pow(1 + pow(R0,2),0.25) + 
       4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        sqrt(1 + pow(R0,2) + 
          2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
       4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (4*R.diff(d010,i,j,k)*q*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (-4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
       (4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        pow(1 + pow(R0,2),0.25) + 
       4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
   (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  4*a0.diff(d001,i,j,k)*q*(1 - pow(1 - rgrid[i],2))*
   pow(1 - rgrid[i],4)*((1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
      ((1 - pow(1 - rgrid[i],2))*
         (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)\
) - (4*pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  (16*B.diff(d010,i,j,k)*q*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k) - 
       L1*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        Qt1(i,j,k))*sqrt(1 + pow(R0 + 
         (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
  (16*B.diff(d001,i,j,k)*q*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     ((L2*a2(i,j,k))/pow(1 + pow(R0,2),0.25) - 
       L2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        Qt2(i,j,k))*sqrt(1 + pow(R0 + 
         (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (pow(L2,2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  (4*a2.diff(d001,i,j,k)*L2*q*pow(1 - rgrid[i],2)*
     (-(pow(1 - pow(1 - rgrid[i],2),2)*pow(Qr2(i,j,k),2)) - 
       (4*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (pow(L2,2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/
   pow(1 + pow(R0,2),0.25) + 
  4*Qt2.diff(d001,i,j,k)*L2*q*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
   pow(1 - rgrid[i],4)*(pow(1 - pow(1 - rgrid[i],2),2)*
      pow(Qr2(i,j,k),2) + (4*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  H*(-2*Qr1.diff(d100,i,j,k)*L1*L2*q*(1 - pow(1 - rgrid[i],2))*
      (1 - rgrid[i])*ygrid[k] - 
     L1*L2*q*(-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
      ygrid[k]*Qr1(i,j,k) - 8*L1*L2*q*pow(1 - rgrid[i],2)*ygrid[k]*
      Qr1(i,j,k) + 8*Qr1.diff(d010,i,j,k)*L1*L2*q*
      pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*ygrid[k]*
      Qr1(i,j,k) + 4*Qr2.diff(d001,i,j,k)*L1*L2*q*
      pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*ygrid[k]*
      Qr1(i,j,k) - (2*Q.diff(d100,i,j,k)*L1*L2*q*
        pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*ygrid[k]*
        Qr1(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
     (4*L1*L2*q*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        ygrid[k]*Q(i,j,k)*Qr1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (4*Q.diff(d010,i,j,k)*L1*L2*q*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*ygrid[k]*pow(Qr1(i,j,k),2))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     4*Qr1.diff(d001,i,j,k)*L1*L2*q*pow(1 - pow(1 - rgrid[i],2),2)*
      pow(1 - rgrid[i],2)*ygrid[k]*Qr2(i,j,k) + 
     (4*Q.diff(d001,i,j,k)*L1*L2*q*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*ygrid[k]*Qr1(i,j,k)*Qr2(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (Qrr.diff(d100,i,j,k)*L1*L2*q*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*ygrid[k]*Qr1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
     (16*R.diff(d001,i,j,k)*q*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*ygrid[k]*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (4*L1*L2*q*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        ygrid[k]*Qr1(i,j,k)*Qrr(i,j,k))/
      (2 + 2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
     (Qtt.diff(d100,i,j,k)*L1*L2*q*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*ygrid[k]*Qr1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
     (4*L1*L2*q*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        ygrid[k]*Qr1(i,j,k)*Qtt(i,j,k))/
      (2 + 2*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
     2*Qrr.diff(d001,i,j,k)*q*(1 - pow(1 - rgrid[i],2))*
      pow(1 - rgrid[i],2)*ygrid[k]*
      (-((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*Qr1(i,j,k)*Qr2(i,j,k))/
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) - 
        (4*pow(rh,2)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
         (pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     (8*Qtt.diff(d001,i,j,k)*q*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*ygrid[k]*
        ((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             pow(-1 + rgrid[i],2)*
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
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
     (16*R.diff(d010,i,j,k)*L2*q*pow(rh,2)*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*ygrid[k]*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     q*(L1*L2*Qr1(i,j,k)*(-2*(1 - pow(1 - rgrid[i],2))*ygrid[k] + 
           4*pow(1 - rgrid[i],2)*
            (2*ygrid[k] + pow(1 - pow(1 - rgrid[i],2),2)*Qr2(i,j,k))) \
- (8*pow(rh,2)*(1 - rgrid[i])*
           (2*R0*(1 - rgrid[i]) + 
             pow(1 - pow(1 - rgrid[i],2),2)*
              (L1*L2*pow(rh,2)*
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
                    (2.*pow(rh,4)))*ygrid[k]*Q(i,j,k)*Qr1(i,j,k) + 
                2*(1 - rgrid[i])*Qrr(i,j,k)*R(i,j,k)) + 
             (1 - pow(1 - rgrid[i],2))*
              (L1*L2*pow(rh,2)*
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
                    (2.*pow(rh,4)))*ygrid[k]*Qr1(i,j,k) + 
                2*(1 - rgrid[i])*(R0*Qrr(i,j,k) + R(i,j,k)))))/
         (pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     (8*Qrr.diff(d010,i,j,k)*L2*q*pow(rh,2)*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*ygrid[k]*
        ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
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
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (8*Qtt.diff(d010,i,j,k)*L2*q*pow(rh,2)*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*ygrid[k]*
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
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (16*B.diff(d010,i,j,k)*L2*q*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*ygrid[k]*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  (2*Qtt.diff(d010,i,j,k)*q*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*(-4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
       4*pow(R0,2)*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
       4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) + 
       4*pow(R0,2)*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) - 
       8*R0*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        (1 - pow(1 - rgrid[i],2))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*R(i,j,k) + 
       8*R0*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*R(i,j,k) - 
       4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*pow(R(i,j,k),2) + 
       4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
        pow(R(i,j,k),2) + (L1*ar(i,j,k)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(rh,2) - (pow(L1,2)*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr1(i,j,k),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(rh,2) - (L1*L2*a2(i,j,k)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          Qr2(i,j,k)*sqrt(1 + 
            pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (pow(1 + pow(R0,2),0.25)*pow(rh,2)) + 
       (4*R0*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(1 + pow(R0,2),0.25) + 
       (pow(L1,2)*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr1(i,j,k),2)*
          Qt1(i,j,k)*sqrt(1 + 
            pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(rh,2) + (L1*L2*a0(i,j,k)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          Qr2(i,j,k)*Qt2(i,j,k)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(rh,2) - 4*R0*a0(i,j,k)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) - 
       (L1*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          Qtr(i,j,k)*sqrt(1 + 
            pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(rh,2) + (4*a2(i,j,k)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 - pow(1 - rgrid[i],2))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*R(i,j,k)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(1 + pow(R0,2),0.25) - 
       4*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*R(i,j,k)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
   (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (2*Qrr.diff(d010,i,j,k)*q*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*(-4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
       4*pow(R0,2)*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
       4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) + 
       4*pow(R0,2)*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        Qt1(i,j,k) - 8*R0*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        (1 - pow(1 - rgrid[i],2))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*R(i,j,k) + 
       8*R0*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        Qt1(i,j,k)*R(i,j,k) - 4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*pow(R(i,j,k),2) + 
       4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        Qt1(i,j,k)*pow(R(i,j,k),2) - 
       (L1*ar(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(rh,2) + (pow(L1,2)*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr1(i,j,k),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(rh,2) + (L1*L2*a2(i,j,k)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          Qr2(i,j,k)*sqrt(1 + 
            pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (pow(1 + pow(R0,2),0.25)*pow(rh,2)) + 
       (4*R0*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(1 + pow(R0,2),0.25) - 
       (pow(L1,2)*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr1(i,j,k),2)*
          Qt1(i,j,k)*sqrt(1 + 
            pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(rh,2) - (L1*L2*a0(i,j,k)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          Qr2(i,j,k)*Qt2(i,j,k)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(rh,2) - 4*R0*a0(i,j,k)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
       (L1*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
          Qtr(i,j,k)*sqrt(1 + 
            pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(rh,2) + (4*a2(i,j,k)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 - pow(1 - rgrid[i],2))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*R(i,j,k)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(1 + pow(R0,2),0.25) - 
       4*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*R(i,j,k)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
   (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))
;

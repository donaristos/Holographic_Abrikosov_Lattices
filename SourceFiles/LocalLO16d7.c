
elem[246]+=
(-2*L2*q*a2(i,j,k)*h1(i,j,k)*(1 - pow(1 - rgrid[i],2)))/
   pow(1 + pow(R0,2),0.25) - 
  4*h2.diff(d101,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]) - 
  (2*a2.diff(d100,i,j,k)*L2*q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     (1 - rgrid[i]))/pow(1 + pow(R0,2),0.25) - 
  4*ar.diff(d001,i,j,k)*q*h1(i,j,k)*pow(1 - rgrid[i],2) + 
  (8*L2*q*a2(i,j,k)*h1(i,j,k)*pow(1 - rgrid[i],2))/
   pow(1 + pow(R0,2),0.25) + 
  2*Qt2.diff(d100,i,j,k)*L2*q*a0(i,j,k)*h1(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3) + 
  4*Qtr.diff(d001,i,j,k)*q*a0(i,j,k)*h1(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4) - 
  (8*L2*q*pow(rh,4)*a2(i,j,k)*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     (1 - rgrid[i])*((pow(-1 + rgrid[i],2)*
          (3*pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],2) + 
            3*pow(H,2)*pow(-2 + rgrid[i],2)*pow(rgrid[i],3) + 
            pow(rh,2)*(3*pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],2) + 
               3*pow(mu,2)*pow(-2 + rgrid[i],2)*pow(rgrid[i],3) + 
               4*pow(rh,2)*(2 + 6*rgrid[i] - 12*pow(rgrid[i],2) + 
                  4*pow(rgrid[i],3)))))/(4.*pow(rh,4)) + 
       ((-1 + rgrid[i])*(pow(H,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
        (2.*pow(rh,4))))/
   (pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
  8*h2.diff(d011,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],2)*Qr1(i,j,k) + 
  (4*a2.diff(d010,i,j,k)*L2*q*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*Qr1(i,j,k))/
   pow(1 + pow(R0,2),0.25) + 
  4*a1.diff(d001,i,j,k)*L1*q*pow(1 + pow(R0,2),0.25)*h1(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*Qr1(i,j,k) - 
  4*Qt1.diff(d001,i,j,k)*L1*q*a0(i,j,k)*h1(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*Qr1(i,j,k) - 
  4*Qt2.diff(d010,i,j,k)*L2*q*a0(i,j,k)*h1(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*Qr1(i,j,k) + 
  8*h2.diff(d002,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],2)*Qr2(i,j,k) + 
  (8*a2.diff(d001,i,j,k)*L2*q*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*Qr2(i,j,k))/
   pow(1 + pow(R0,2),0.25) - 
  (8*pow(L2,2)*pow(q,2)*pow(a2(i,j,k),2)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*Qr2(i,j,k))/
   sqrt(1 + pow(R0,2)) - 8*Qt2.diff(d001,i,j,k)*L2*q*a0(i,j,k)*
   h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
   Qr2(i,j,k) + (8*L1*L2*pow(q,2)*a0(i,j,k)*a2(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*Qr1(i,j,k)*
     Qt1(i,j,k))/pow(1 + pow(R0,2),0.25) + 
  4*Qr1.diff(d001,i,j,k)*q*h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],2)*(L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k) - 
     L1*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
      Qt1(i,j,k)) - 2*L2*q*a0(i,j,k)*h1(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*Qt2(i,j,k) + 
  2*a0.diff(d100,i,j,k)*L2*q*h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],3)*Qt2(i,j,k) - 
  4*L2*q*a0(i,j,k)*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*
   pow(1 - rgrid[i],4)*Qt2(i,j,k) + 
  (8*L2*q*pow(rh,4)*a0(i,j,k)*h1(i,j,k)*
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
        (2.*pow(rh,4)))*Qt2(i,j,k))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
  (16*pow(L2,2)*pow(q,2)*a0(i,j,k)*a2(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*Qr2(i,j,k)*
     Qt2(i,j,k))/pow(1 + pow(R0,2),0.25) - 
  8*L1*L2*pow(q,2)*pow(a0(i,j,k),2)*h2(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],6)*Qr1(i,j,k)*
   Qt1(i,j,k)*Qt2(i,j,k) - 8*pow(L2,2)*pow(q,2)*pow(a0(i,j,k),2)*
   h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],6)*
   Qr2(i,j,k)*pow(Qt2(i,j,k),2) + 
  h1.diff(d100,i,j,k)*L2*q*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
   ((-4*a2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
     4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
      Qt2(i,j,k)) + 2*L2*q*h1(i,j,k)*pow(1 - rgrid[i],2)*
   ((-4*a2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
     4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
      Qt2(i,j,k)) + (L2*q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*Q(i,j,k)*
     ((-4*a2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
       4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        Qt2(i,j,k)))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
  (L2*q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
     Qrr(i,j,k)*((-4*a2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
       4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        Qt2(i,j,k)))/(2.*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
  q*h1(i,j,k)*(-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
   (-((L2*a2(i,j,k))/pow(1 + pow(R0,2),0.25)) + 
     L2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
      Qt2(i,j,k)) - (2*L2*pow(q,2)*ar(i,j,k)*h2(i,j,k)*
     pow(1 - rgrid[i],2)*(-4*a2(i,j,k) + 
       4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*Qt2(i,j,k)))/pow(1 + pow(R0,2),0.25) - 
  (2*h1.diff(d010,i,j,k)*L2*q*pow(1 - pow(1 - rgrid[i],2),2)*
     pow(1 - rgrid[i],2)*Qr1(i,j,k)*
     (-4*a2(i,j,k) + 4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qt2(i,j,k)))/
   pow(1 + pow(R0,2),0.25) + 
  2*L1*L2*pow(q,2)*a1(i,j,k)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],2)*Qr1(i,j,k)*
   (-4*a2(i,j,k) + 4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
      (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qt2(i,j,k)) + 
  Qr2.diff(d001,i,j,k)*(8*h2.diff(d001,i,j,k)*
      pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2) + 
     2*L2*q*h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
      pow(1 - rgrid[i],2)*((4*a2(i,j,k))/pow(1 + pow(R0,2),0.25) - 
        4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
         Qt2(i,j,k))) + Q.diff(d100,i,j,k)*
   ((-2*h2.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i]))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (L2*q*h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        ((-4*a2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
          4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           Qt2(i,j,k)))/(2.*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  Qrr.diff(d100,i,j,k)*((h2.diff(d001,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i]))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
     (L2*q*h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        ((-4*a2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
          4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           Qt2(i,j,k)))/(4.*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))) + 
  Qrr.diff(d010,i,j,k)*((-2*h2.diff(d001,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     (L2*q*h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k)*
        (-4*a2(i,j,k) + 4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
           (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qt2(i,j,k)))/
      (2.*pow(1 + pow(R0,2),0.25)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))) + 
  Qr1.diff(d010,i,j,k)*(4*h2.diff(d001,i,j,k)*
      pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2) + 
     2*q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
      ((2*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]))/
         pow(1 + pow(R0,2),0.25) - 
        2*L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
         pow(1 - rgrid[i],3)*Qt2(i,j,k))) + 
  Q.diff(d010,i,j,k)*((4*h2.diff(d001,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (2*q*h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        Qr1(i,j,k)*((2*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*
             (1 - rgrid[i]))/pow(1 + pow(R0,2),0.25) - 
          2*L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],3)*Qt2(i,j,k)))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
  (8*L2*pow(q,2)*a0(i,j,k)*a2(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*Qtr(i,j,k))/
   pow(1 + pow(R0,2),0.25) + 
  8*L2*pow(q,2)*pow(a0(i,j,k),2)*h2(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],6)*Qt2(i,j,k)*
   Qtr(i,j,k) + 4*a0.diff(d001,i,j,k)*q*h1(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
   (-((1 - pow(1 - rgrid[i],2))*
        (L1*Qr1(i,j,k)*Qt1(i,j,k) + 2*L2*Qr2(i,j,k)*Qt2(i,j,k))) + 
     Qtr(i,j,k)) + (L2*q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*((-4*a2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
       4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        Qt2(i,j,k))*Qtt(i,j,k))/
   (2.*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
  h2.diff(d100,i,j,k)*((-2*Q.diff(d001,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i]))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (Qrr.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i]))/(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
     (Qtt.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i]))/(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
  h2.diff(d010,i,j,k)*(4*Qr1.diff(d001,i,j,k)*
      pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2) + 
     (4*Q.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
     (2*Qrr.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     (2*Qtt.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
  H*(-4*Qr1.diff(d001,i,j,k)*L1*L2*q*h1(i,j,k)*
      pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*ygrid[k] - 
     8*h1.diff(d001,i,j,k)*L1*L2*q*pow(1 - pow(1 - rgrid[i],2),2)*
      pow(1 - rgrid[i],2)*ygrid[k]*Qr1(i,j,k) - 
     (4*Q.diff(d001,i,j,k)*L1*L2*q*h1(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*ygrid[k]*
        Qr1(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (2*Qrr.diff(d001,i,j,k)*L1*L2*q*h1(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*ygrid[k]*
        Qr1(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
     (2*L1*L2*q*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k)*(2*pow(1 + pow(R0,2),0.25)*h1(i,j,k) + 
          L2*q*h2(i,j,k)*ygrid[k]*
           (-4*a2(i,j,k) + 4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
              (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qt2(i,j,k)\
)))/pow(1 + pow(R0,2),0.25) - 
     (2*Qtt.diff(d001,i,j,k)*L1*L2*q*h1(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*ygrid[k]*
        Qr1(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
  Qtt.diff(d100,i,j,k)*(-((h2.diff(d001,i,j,k)*
          pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i]))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
     (L2*q*h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        ((-4*a2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
          4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           Qt2(i,j,k)))/(4.*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))) + 
  Qtt.diff(d010,i,j,k)*((2*h2.diff(d001,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
     (L2*q*h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k)*
        (-4*a2(i,j,k) + 4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
           (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qt2(i,j,k)))/
      (2.*pow(1 + pow(R0,2),0.25)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))) + 
  h2.diff(d001,i,j,k)*(-4*pow(1 - rgrid[i],2) - 
     (8*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))) - 
     (4*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Q(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (4*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qrr(i,j,k))/
      (2 + 2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
     (4*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qtt(i,j,k))/
      (2 + 2*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
  (2*h1.diff(d001,i,j,k)*q*pow(1 - rgrid[i],2)*
     (4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*
        (-8*L2*a2(i,j,k)*Qr2(i,j,k) + 
          4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*pow(1 - rgrid[i],2)*
           (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) + 
             2*L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
             Qtr(i,j,k)) - 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))\
))/pow(1 + pow(R0,2),0.25) + 
  Qrr.diff(d001,i,j,k)*((-4*h2.diff(d001,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     (4*h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2))/
      (2 + 2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     (q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*
           (-8*L2*a2(i,j,k)*Qr2(i,j,k) + 
             4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
              pow(1 - rgrid[i],2)*
              (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) + 
                2*L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                Qtr(i,j,k)) - 4*L1*a1(i,j,k)*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2)))))/
      (2.*pow(1 + pow(R0,2),0.25)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))) + 
  Q.diff(d001,i,j,k)*((-4*h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (8*h2.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (2*q*h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
        Qr2(i,j,k)*((2*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*
             (1 - rgrid[i]))/pow(1 + pow(R0,2),0.25) - 
          2*L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],3)*Qt2(i,j,k)))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (-4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*
           (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
             4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
              pow(1 - rgrid[i],2)*
              (-((1 - pow(1 - rgrid[i],2))*
                   (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))) \
+ Qtr(i,j,k)) + 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))))/
      (pow(1 + pow(R0,2),0.25)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  Qtt.diff(d001,i,j,k)*((4*h2.diff(d001,i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
     (4*h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2))/
      (2 + 2*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
     (q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (-4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*
           (8*L2*a2(i,j,k)*Qr2(i,j,k) + 
             4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
              pow(1 - rgrid[i],2)*
              (-((1 - pow(1 - rgrid[i],2))*
                   (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                     2*L2*Qr2(i,j,k)*Qt2(i,j,k))) + Qtr(i,j,k)) + 
             4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))))/
      (2.*pow(1 + pow(R0,2),0.25)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))) - 
  (4*a0.diff(d010,i,j,k)*L2*q*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*Qr1(i,j,k)*
     Qt2(i,j,k)*sqrt(1 + pow(R0,2) + 
       2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
   sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))
;

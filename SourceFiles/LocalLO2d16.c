
elem[31]+=
4*h2.diff(d100,i,j,k)*(1 - rgrid[i])*
   (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
  (8*h2(i,j,k)*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   (1 - pow(1 - rgrid[i],2)) - 
  8*h2.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
   Qr1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
  8*h2.diff(d001,i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
   Qr2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
  H*(-4*h1.diff(d100,i,j,k)*L1*L2*q*pow(1 - pow(1 - rgrid[i],2),2)*
      (1 - rgrid[i])*ygrid[k]*Qr1(i,j,k)*
      (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     8*h1.diff(d010,i,j,k)*L1*L2*q*pow(1 - pow(1 - rgrid[i],2),3)*
      pow(1 - rgrid[i],2)*ygrid[k]*pow(Qr1(i,j,k),2)*
      (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     8*h1.diff(d001,i,j,k)*L1*L2*q*pow(1 - pow(1 - rgrid[i],2),3)*
      pow(1 - rgrid[i],2)*ygrid[k]*Qr1(i,j,k)*Qr2(i,j,k)*
      (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     8*L1*L2*pow(q,2)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
      (1 - rgrid[i])*ygrid[k]*Qr1(i,j,k)*
      (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
      ((2*ar(i,j,k)*(1 - rgrid[i]))/(1 - pow(1 - rgrid[i],2)) - 
        2*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
         (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) - 
        (2*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
           Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
        2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],3)*
         ((1 - pow(1 - rgrid[i],2))*
            (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
           Qtr(i,j,k)))) - 2*h1.diff(d100,i,j,k)*q*
   (1 - pow(1 - rgrid[i],2))*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
   ((2*ar(i,j,k)*(1 - rgrid[i]))/(1 - pow(1 - rgrid[i],2)) - 
     2*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
      (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) - 
     (2*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
     2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],3)*
      ((1 - pow(1 - rgrid[i],2))*
         (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k))\
) + 4*h1.diff(d010,i,j,k)*q*pow(1 - pow(1 - rgrid[i],2),2)*
   (1 - rgrid[i])*Qr1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
   ((2*ar(i,j,k)*(1 - rgrid[i]))/(1 - pow(1 - rgrid[i],2)) - 
     2*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
      (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) - 
     (2*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
     2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],3)*
      ((1 - pow(1 - rgrid[i],2))*
         (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k))\
) + 4*h1.diff(d001,i,j,k)*q*pow(1 - pow(1 - rgrid[i],2),2)*
   (1 - rgrid[i])*Qr2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
   ((2*ar(i,j,k)*(1 - rgrid[i]))/(1 - pow(1 - rgrid[i],2)) - 
     2*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
      (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) - 
     (2*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
     2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],3)*
      ((1 - pow(1 - rgrid[i],2))*
         (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k))\
) + (2*pow(rh,4)*h2(i,j,k)*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (-16*sqrt(1 + pow(R0,2)) - 
       192*S*(pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
        pow(1 - pow(1 - rgrid[i],2),2)*sqrt(1 + pow(R0,2)) + 
       (4*pow(q,2)*pow(ar(i,j,k),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          sqrt(1 + pow(R0,2)))/pow(rh,4) - 
       16*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)*sqrt(1 + pow(R0,2)) - 
       192*S*(pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
        pow(1 - pow(1 - rgrid[i],2),3)*Qrr(i,j,k)*
        sqrt(1 + pow(R0,2)) + 
       (4*pow(q,2)*pow(a0(i,j,k),2)*
          pow(1 - pow(1 - rgrid[i],2),6)*pow(1 - rgrid[i],4)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          pow(L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k),2)*
          sqrt(1 + pow(R0,2)))/pow(rh,4) - 
       (2*pow(q,2)*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
          pow(1 - pow(1 - rgrid[i],2),5)*pow(1 - rgrid[i],2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))*
          (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
            4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
             pow(1 - rgrid[i],2)*Qtr(i,j,k) + 
            4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2))))/pow(rh,4) \
- (2*pow(q,2)*pow(1 + pow(R0,2),0.25)*ar(i,j,k)*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
            4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
             pow(1 - rgrid[i],2)*
             (-((1 - pow(1 - rgrid[i],2))*
                  (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))\
) + Qtr(i,j,k)) + 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2))))/
        pow(rh,4) + (pow(q,2)*pow(1 - pow(1 - rgrid[i],2),4)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (16*pow(L1,2)*(1 + pow(R0,2))*pow(a1(i,j,k),2)*
             pow(Qr1(i,j,k),2) + 
            pow(4*L2*a2(i,j,k)*Qr2(i,j,k) + 
              4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
               pow(1 - rgrid[i],2)*Qtr(i,j,k),2) + 
            8*L1*a1(i,j,k)*Qr1(i,j,k)*
             (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
               4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
                pow(1 - rgrid[i],2)*Qtr(i,j,k))*sqrt(1 + pow(R0,2))))/
        (4.*pow(rh,4))))/
   (sqrt(1 + pow(R0,2))*(1 - pow(1 - rgrid[i],2))*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
  (8*pow(H,2)*pow(L1,2)*pow(L2,2)*pow(q,2)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
     pow(ygrid[k],2)*pow(Qr1(i,j,k),2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
   sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))
;

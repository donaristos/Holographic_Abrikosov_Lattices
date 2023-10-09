
elem[253]+=
(q*h1(i,j,k)*(-2*(1 - pow(1 - rgrid[i],2)) - 4*pow(1 - rgrid[i],2)))/
   pow(1 - pow(1 - rgrid[i],2),2) + 
  (4*h1.diff(d100,i,j,k)*q*(1 - rgrid[i]))/(1 - pow(1 - rgrid[i],2)) - 
  4*Qr1.diff(d010,i,j,k)*q*h1(i,j,k)*pow(1 - rgrid[i],2) - 
  4*Qr2.diff(d001,i,j,k)*q*h1(i,j,k)*pow(1 - rgrid[i],2) + 
  (8*q*h1(i,j,k)*pow(1 - rgrid[i],2))/pow(1 - pow(1 - rgrid[i],2),2) - 
  (8*pow(q,2)*ar(i,j,k)*h2(i,j,k)*pow(1 - rgrid[i],2))/
   pow(1 - pow(1 - rgrid[i],2),2) + 
  (2*Q.diff(d100,i,j,k)*q*h1(i,j,k)*(1 - rgrid[i]))/
   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
  (4*q*h1(i,j,k)*pow(1 - rgrid[i],2)*Q(i,j,k))/
   (1 - pow(1 - rgrid[i],2) + pow(1 - pow(1 - rgrid[i],2),2)*Q(i,j,k)) \
- 8*h1.diff(d010,i,j,k)*q*pow(1 - rgrid[i],2)*Qr1(i,j,k) - 
  8*H*L1*L2*pow(q,2)*h2(i,j,k)*pow(1 - rgrid[i],2)*ygrid[k]*Qr1(i,j,k) - 
  (4*Q.diff(d010,i,j,k)*q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*Qr1(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
  8*h1.diff(d001,i,j,k)*q*pow(1 - rgrid[i],2)*Qr2(i,j,k) - 
  (4*Q.diff(d001,i,j,k)*q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*Qr2(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
  (Qrr.diff(d100,i,j,k)*q*h1(i,j,k)*(1 - rgrid[i]))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
  (2*Qrr.diff(d010,i,j,k)*q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*Qr1(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
  (2*Qrr.diff(d001,i,j,k)*q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*Qr2(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
  (4*q*h1(i,j,k)*pow(1 - rgrid[i],2)*Qrr(i,j,k))/
   (2*(1 - pow(1 - rgrid[i],2)) + 
     2*pow(1 - pow(1 - rgrid[i],2),2)*Qrr(i,j,k)) + 
  (Qtt.diff(d100,i,j,k)*q*h1(i,j,k)*(1 - rgrid[i]))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
  (2*Qtt.diff(d010,i,j,k)*q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*Qr1(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
  (2*Qtt.diff(d001,i,j,k)*q*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*Qr2(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
  (4*q*h1(i,j,k)*pow(1 - rgrid[i],2)*Qtt(i,j,k))/
   (2*(1 - pow(1 - rgrid[i],2)) + 
     2*pow(1 - pow(1 - rgrid[i],2),2)*Qtt(i,j,k)) + 
  (q*((-16*h1(i,j,k)*pow(1 - rgrid[i],2))/
        pow(1 - pow(1 - rgrid[i],2),2) + 
       (8*pow(rh,4)*h1(i,j,k)*
          ((pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
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
                (2.*pow(rh,4)))))/
        ((1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))) - 
       16*q*a0(i,j,k)*h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],4)*(L1*Qr1(i,j,k)*Qt1(i,j,k) + 
          L2*Qr2(i,j,k)*Qt2(i,j,k)) + 
       (4*q*h2(i,j,k)*pow(1 - rgrid[i],2)*
          (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
            4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*pow(1 - rgrid[i],2)*
             Qtr(i,j,k) + 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2))))/
        pow(1 + pow(R0,2),0.25)))/2.
;

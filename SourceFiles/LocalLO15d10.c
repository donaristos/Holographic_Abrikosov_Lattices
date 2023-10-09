
elem[233]+=
-2*q*a0(i,j,k)*h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2) + 
  q*a0(i,j,k)*h2(i,j,k)*(-2*(1 - pow(1 - rgrid[i],2)) + 
     4*pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2) + 
  4*h2.diff(d100,i,j,k)*q*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
   pow(1 - rgrid[i],3) + 2*a0.diff(d100,i,j,k)*q*h2(i,j,k)*
   (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],3) + 
  8*pow(q,2)*a0(i,j,k)*ar(i,j,k)*h1(i,j,k)*pow(1 - rgrid[i],4) - 
  4*Qr1.diff(d010,i,j,k)*q*a0(i,j,k)*h2(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4) - 
  4*Qr2.diff(d001,i,j,k)*q*a0(i,j,k)*h2(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4) + 
  (8*q*pow(rh,4)*a0(i,j,k)*h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],3)*((pow(-1 + rgrid[i],2)*
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
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
  (2*Q.diff(d100,i,j,k)*q*a0(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3))/
   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
  (16*q*a0(i,j,k)*h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],4)*Q(i,j,k))/
   (4 + 4*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
  8*h2.diff(d010,i,j,k)*q*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],4)*Qr1(i,j,k) - 
  8*L1*pow(q,2)*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*a1(i,j,k)*h1(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*Qr1(i,j,k) + 
  8*H*L1*L2*pow(q,2)*a0(i,j,k)*h1(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*ygrid[k]*
   Qr1(i,j,k) - (4*Q.diff(d010,i,j,k)*q*a0(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*Qr1(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
  8*h2.diff(d001,i,j,k)*q*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],4)*Qr2(i,j,k) - 
  (8*L2*pow(q,2)*a0(i,j,k)*a2(i,j,k)*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*Qr2(i,j,k))/
   pow(1 + pow(R0,2),0.25) - 
  4*a0.diff(d001,i,j,k)*q*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
   pow(1 - rgrid[i],4)*Qr2(i,j,k) - 
  (4*Q.diff(d001,i,j,k)*q*a0(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*Qr2(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
  (Qrr.diff(d100,i,j,k)*q*a0(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
  (2*Qrr.diff(d010,i,j,k)*q*a0(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*Qr1(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
  (2*Qrr.diff(d001,i,j,k)*q*a0(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*Qr2(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
  (16*q*a0(i,j,k)*h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],4)*Qrr(i,j,k))/
   (8 + 8*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
  8*L1*pow(q,2)*pow(a0(i,j,k),2)*h1(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],6)*Qr1(i,j,k)*
   Qt1(i,j,k) + 8*L2*pow(q,2)*pow(a0(i,j,k),2)*h1(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],6)*Qr2(i,j,k)*
   Qt2(i,j,k) - 8*pow(q,2)*pow(a0(i,j,k),2)*h1(i,j,k)*
   pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],6)*Qtr(i,j,k) + 
  (Qtt.diff(d100,i,j,k)*q*a0(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
  (2*Qtt.diff(d010,i,j,k)*q*a0(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*Qr1(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
  (2*Qtt.diff(d001,i,j,k)*q*a0(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*Qr2(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
  (16*q*a0(i,j,k)*h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],4)*Qtt(i,j,k))/
   (8 + 8*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
  (4*a0.diff(d010,i,j,k)*q*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
     pow(1 - rgrid[i],4)*Qr1(i,j,k)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
   sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))
;

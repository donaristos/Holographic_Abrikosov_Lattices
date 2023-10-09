
if(rgrid.IsRightBoundary(i)){
elem[12]=
a2(i,j,k)-a2(i+1,j,k)
;
}
else if(rgrid.IsLeftBoundary(i)){
elem[12]=
a2.diff(d100,i-1,j,k)-a2.diff(d100,i,j,k)
;
}
else if(k==0){
elem[12]=
a2(i,j,0)-a2(i,j,ygrid.TotalLength()-1)
;
}
else if(k==ygrid.TotalLength()-1){
elem[12]=
a2.diff(d001,i,j,ygrid.TotalLength()-1)-a2.diff(d001,i,j,0)
;
}
else{
{elem[12]=
a2.diff(d200,i,j,k) + a2.diff(d100,i,j,k)*
   (1/(1 - rgrid[i]) + (4*pow(rh,4)*
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
                4*pow(rgrid[i],3) + pow(rgrid[i],4)))))) - 
  4*a2.diff(d110,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
   Qr1(i,j,k) + 2*a2.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*
   (1 - rgrid[i])*(-(1/(1 - rgrid[i])) - 
     (4*pow(rh,4)*((pow(-1 + rgrid[i],2)*
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
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))))*Qr1(i,j,k) - 
  4*a2.diff(d101,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
   Qr2(i,j,k) - 8*a2.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*Qr2(i,j,k) - 
  (16*h2.diff(d001,i,j,k)*q*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
     h1(i,j,k)*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
  (16*h1.diff(d001,i,j,k)*q*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
     h2(i,j,k)*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) - 
  (16*pow(q,2)*pow(rh,4)*a2(i,j,k)*
     (pow(h1(i,j,k),2)*pow(1 - pow(1 - rgrid[i],2),2) + 
       pow(h2(i,j,k),2)*pow(1 - pow(1 - rgrid[i],2),2))*
     pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
  (pow(1 + pow(R0,2),0.25)*
     (2*Qr1.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2))*
        (1 - rgrid[i]) + (-2*(1 - pow(1 - rgrid[i],2)) + 
          4*pow(1 - rgrid[i],2))*Qr1(i,j,k))*
     (-((a2.diff(d010,i,j,k)*L2)/pow(1 + pow(R0,2),0.25)) + 
       a1.diff(d001,i,j,k)*L1*pow(1 + pow(R0,2),0.25) - 
       a0.diff(d001,i,j,k)*L1*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*Qt1(i,j,k)))/L2 - 
  (2*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
     (-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
     (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
  (8*a0.diff(d010,i,j,k)*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
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
        (2.*pow(rh,4)))*Qr1(i,j,k)*Qt2(i,j,k))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))) - 
  a0.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
   (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
   (2*Qr2.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]) + 
     (-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*Qr2(i,j,k))*
   Qt2(i,j,k) + (pow(1 + pow(R0,2),0.25)*
     (Q.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)) + 
       2*(1 - rgrid[i])*Q(i,j,k))*
     ((L2*(-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qt2(i,j,k))/
        (2.*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
       (2*a0.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
          pow(1 - rgrid[i],3)*Qr1(i,j,k)*Qt2(i,j,k))/
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
       (2*a0.diff(d001,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
          pow(1 - rgrid[i],3)*Qr2(i,j,k)*Qt2(i,j,k))/
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))))/L2 + 
  (a0.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*pow(1 - rgrid[i],2)*
     (2*Qtr.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2))*
        (1 - rgrid[i]) + (-2*(1 - pow(1 - rgrid[i],2)) + 
          4*pow(1 - rgrid[i],2))*Qtr(i,j,k)))/L2 + 
  (a0.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*pow(1 - rgrid[i],2)*
     ((-4*pow(rh,4)*((pow(-1 + rgrid[i],2)*
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
          (4*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr1(i,j,k)*Qt1(i,j,k) + 
            2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr2(i,j,k)*Qt2(i,j,k) - 
            4*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)))/
        (pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
       ((-2*(1 - pow(1 - rgrid[i],2)) + 8*pow(1 - rgrid[i],2))*
          (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr1(i,j,k)*Qt1(i,j,k) + 
            2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr2(i,j,k)*Qt2(i,j,k) - 
            2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)))/
        (2.*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]))))/L2 + 
  (pow(1 + pow(R0,2),0.25)*
     (Qrr.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)) + 
       2*(1 - rgrid[i])*Qrr(i,j,k))*
     (-(a2.diff(d100,i,j,k)*L2)/
        (2.*pow(1 + pow(R0,2),0.25)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
       (ar.diff(d001,i,j,k)*(1 - rgrid[i]))/
        ((1 - pow(1 - rgrid[i],2))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
       (a2.diff(d010,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*Qr1(i,j,k))/
        (pow(1 + pow(R0,2),0.25)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) - 
       (a1.diff(d001,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
       (a0.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*
          (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr1(i,j,k)*Qt1(i,j,k) + 
            2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr2(i,j,k)*Qt2(i,j,k) - 
            2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)))/
        (2.*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))))/L2 + 
  (pow(1 + pow(R0,2),0.25)*
     (Qtt.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)) + 
       2*(1 - rgrid[i])*Qtt(i,j,k))*
     ((a2.diff(d100,i,j,k)*L2)/
        (2.*pow(1 + pow(R0,2),0.25)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
       (ar.diff(d001,i,j,k)*(1 - rgrid[i]))/
        ((1 - pow(1 - rgrid[i],2))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
       (a2.diff(d010,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*Qr1(i,j,k))/
        (pow(1 + pow(R0,2),0.25)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
       (a1.diff(d001,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
       (L2*(-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qt2(i,j,k))/
        (2.*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
       (2*a0.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
          pow(1 - rgrid[i],3)*Qr1(i,j,k)*Qt2(i,j,k))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
       (a0.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*
          (-2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr1(i,j,k)*Qt1(i,j,k) + 
            2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr2(i,j,k)*Qt2(i,j,k) + 
            2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)))/
        (2.*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))))/L2 + 
  a2.diff(d011,i,j,k)*(8*pow(1 - pow(1 - rgrid[i],2),2)*
      pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qr2(i,j,k) - 
     (32*pow(rh,2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  (pow(1 + pow(R0,2),0.25)*
     (R.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)) + 
       2*(1 - rgrid[i])*R(i,j,k))*
     (-((a1.diff(d100,i,j,k)*L2*pow(1 + pow(R0,2),0.25))/
          ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
            sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) \
+ (2*ar.diff(d010,i,j,k)*L2*(1 - rgrid[i]))/
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 - pow(1 - rgrid[i],2))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (2*a2.diff(d010,i,j,k)*pow(L2,2)*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*Qr2(i,j,k))/
        (L1*pow(1 + pow(R0,2),0.25)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (2*a1.diff(d001,i,j,k)*L2*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k))/
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (L2*(-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qt1(i,j,k))/
        (2.*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (2*a0.diff(d001,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
          pow(1 - rgrid[i],3)*Qr2(i,j,k)*Qt1(i,j,k))/
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a0.diff(d010,i,j,k)*L2*pow(1 - rgrid[i],2)*
          (2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr2(i,j,k)*Qt2(i,j,k) - 
            2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)))/
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))))/L2 \
+ (R.diff(d010,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*
     ((2*a1.diff(d100,i,j,k)*L2*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k))/
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (4*ar.diff(d010,i,j,k)*L2*pow(1 - rgrid[i],2)*Qr1(i,j,k))/
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (4*a2.diff(d010,i,j,k)*pow(L2,2)*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
          Qr1(i,j,k)*Qr2(i,j,k))/
        (L1*pow(1 + pow(R0,2),0.25)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (4*a1.diff(d001,i,j,k)*L2*pow(1 + pow(R0,2),0.25)*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
          Qr1(i,j,k)*Qr2(i,j,k))/
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (L2*(-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
          Qr1(i,j,k)*Qt1(i,j,k))/
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (4*a0.diff(d010,i,j,k)*pow(rh,2)*pow(1 - rgrid[i],2)*
          ((L1*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
               pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
               (-2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                  (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) + 
                 2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                  Qtr(i,j,k)))/(2.*pow(rh,2)) + 
            4*L2*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))/
        (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a0.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*
          ((-16*L2*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
               pow(1 - rgrid[i],2)*
               (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k))/
             (pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
            (4*L1*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),3)*
               pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qr2(i,j,k)*Qt1(i,j,k))/
             ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),
                  2)))))/(L1*L2)))/L2 + 
  (2*a1.diff(d001,i,j,k)*L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
     (1/(1 - rgrid[i]) - (4*(1 - rgrid[i]))/(1 - pow(1 - rgrid[i],2)) + 
       (4*pow(rh,4)*((pow(-1 + rgrid[i],2)*
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
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))))*Qr1(i,j,k)*
     sqrt(1 + pow(R0,2)))/L2 + 
  (Qtt.diff(d010,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*
     (-((a2.diff(d100,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
            (1 - rgrid[i])*Qr1(i,j,k))/
          (pow(1 + pow(R0,2),0.25)*
            (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))) + 
       (2*ar.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*Qr1(i,j,k))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
       (L2*(-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
          Qr1(i,j,k)*Qt2(i,j,k))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
       (2*a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
          (-4*pow(1 - rgrid[i],2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
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
             pow(rh,2)))/
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (2*a2.diff(d010,i,j,k)*L2*pow(rh,2)*
          (4*pow(1 - rgrid[i],2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
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
             pow(rh,2)))/
        (pow(L1,2)*pow(1 + pow(R0,2),0.25)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (4*a0.diff(d010,i,j,k)*L2*pow(rh,2)*
          (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qt2(i,j,k)*
          (4*pow(1 - rgrid[i],2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
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
             pow(rh,2)))/
        (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (2*a0.diff(d001,i,j,k)*pow(rh,2)*pow(1 - rgrid[i],2)*
          (4*L1*L2*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             ((L1*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
                  (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                     (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) - 
                    2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                     (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                    2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                     Qtr(i,j,k)))/(2.*pow(rh,2)) + 
               4*L2*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (pow(L1,2)*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))))/L2 \
+ (Qrr.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*
     (-((ar(i,j,k)*(-2*(1 - pow(1 - rgrid[i],2)) - 
               4*pow(1 - rgrid[i],2)) + 
            2*ar.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2))*
             (1 - rgrid[i]))/
          (pow(1 - pow(1 - rgrid[i],2),2)*
            (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))) + 
       (2*a1.diff(d100,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
       (4*ar.diff(d010,i,j,k)*pow(1 - rgrid[i],2)*Qr1(i,j,k))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
       (4*a1.diff(d010,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
          pow(Qr1(i,j,k),2))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
       (3*a2.diff(d100,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*Qr2(i,j,k))/
        (pow(1 + pow(R0,2),0.25)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
       (2*ar.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*Qr2(i,j,k))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
       (4*a2.diff(d001,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
          pow(1 - rgrid[i],2)*pow(Qr2(i,j,k),2))/
        (pow(1 + pow(R0,2),0.25)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
       (2*a0.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*
          pow(1 - rgrid[i],3)*Qr1(i,j,k)*
          (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr1(i,j,k)*Qt1(i,j,k) + 
            2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr2(i,j,k)*Qt2(i,j,k) - 
            2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
       ((-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - rgrid[i])*(-2*L1*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) - 
            2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr2(i,j,k)*Qt2(i,j,k) + 
            2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)))/
        (2.*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
       (a2.diff(d010,i,j,k)*L2*
          ((-6*pow(1 - pow(1 - rgrid[i],2),2)*
               pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qr2(i,j,k))/
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
            (8*pow(rh,2)*pow(1 - rgrid[i],2)*
               (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
             (L1*L2*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/
        pow(1 + pow(R0,2),0.25) + 
       (a1.diff(d001,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*
          ((-4*pow(1 - pow(1 - rgrid[i],2),2)*
               pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qr2(i,j,k))/
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
            (16*pow(rh,2)*pow(1 - rgrid[i],2)*
               (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
             (L1*L2*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/2. + 
       (pow(rh,4)*(8*pow(1 - rgrid[i],2)*
             (((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  (2*L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
                     (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                     Qr1(i,j,k) + 
                    (2*L2*a2(i,j,k)*(1 - pow(1 - rgrid[i],2))*
                       (1 - rgrid[i])*Qr2(i,j,k))/
                     pow(1 + pow(R0,2),0.25)))/pow(rh,4) + 
               a0(i,j,k)*pow(1 - rgrid[i],2)*
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
                      (2.*pow(rh,4))))*
                (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                   (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) + 
                  2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                   (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                  2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                   Qtr(i,j,k))) + 
            (2*ar(i,j,k)*(1 - rgrid[i])*
               (-((4*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]) + 
                       8*pow(1 - rgrid[i],3))*
                     pow(-1 + rgrid[i],2)*
                     (pow(H,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
                  (2.*pow(rh,4)) - 
                 4*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
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
                     (2.*pow(rh,4))) - 
                 ((1 - rgrid[i])*pow(-1 + rgrid[i],6)*
                    pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      pow(rh,2)*
                       (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))),3)*
                    pow(2*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                       (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) + 
                      2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                       (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                      2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                       Qtr(i,j,k),2))/(16.*pow(rh,10)) + 
                 (pow(-1 + rgrid[i],4)*
                    pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      pow(rh,2)*
                       (pow(mu,2)*pow(-2 + rgrid[i],3)*
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
                       pow(2*L1*
                       pow(1 - pow(1 - rgrid[i],2),2)*
                       (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) + 
                        2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                       (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                        2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                        Qtr(i,j,k),2) + 
                      4*(1 - rgrid[i])*
                       (8*L1*L2*R0*pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qr2(i,j,k) + 
                        4*pow(L1,2)*
                        pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(1 - rgrid[i],2)*pow(Qr1(i,j,k),2)*
                        sqrt(1 + pow(R0,2)) + 
                        4*pow(L2,2)*
                        pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(1 - rgrid[i],2)*pow(Qr2(i,j,k),2)*
                        sqrt(1 + pow(R0,2)))))/(16.*pow(rh,6))))/
             (1 - pow(1 - rgrid[i],2))))/
        (2.*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
       (a0.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*
          ((2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
               (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  Qr1(i,j,k)*Qt1(i,j,k) + 
                 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                  (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                 2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)\
))/(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
            (16*pow(rh,2)*pow(1 - rgrid[i],2)*
               (L1*L2*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
                  (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
                 L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                  (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
))/(L1*pow(L2,2)*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/2.))/L2 + 
  (Qtt.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*
     ((4*pow(rh,4)*ar(i,j,k)*(1 - rgrid[i])*
          (((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
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
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
                (2.*pow(rh,4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
        (pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k),2)) - 
       (a2.diff(d100,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*Qr2(i,j,k))/
        (pow(1 + pow(R0,2),0.25)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
       (2*ar.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*Qr2(i,j,k))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
       (L2*(-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
          Qr2(i,j,k)*Qt2(i,j,k))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
       (a2.diff(d010,i,j,k)*L2*
          (4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             Qr1(i,j,k)*Qr2(i,j,k) - 
            (16*pow(rh,2)*pow(1 - rgrid[i],2)*
               (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
               (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
             (L1*L2*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/
        (2.*pow(1 + pow(R0,2),0.25)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
       (a1.diff(d001,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*
          (-4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             Qr1(i,j,k)*Qr2(i,j,k) + 
            (16*pow(rh,2)*pow(1 - rgrid[i],2)*
               (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
               (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
             (L1*L2*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/
        (2.*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
       (a0.diff(d010,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          pow(1 - rgrid[i],2)*Qt2(i,j,k)*
          (-4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             Qr1(i,j,k)*Qr2(i,j,k) + 
            (16*pow(rh,2)*pow(1 - rgrid[i],2)*
               (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
               (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
             (L1*L2*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
       (a0.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*
          (2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (-2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                Qr1(i,j,k)*Qt1(i,j,k) + 
               2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                Qr2(i,j,k)*Qt2(i,j,k) + 
               2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)) + 
            (16*pow(rh,2)*pow(1 - rgrid[i],2)*
               (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
               (L1*L2*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
                  (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
                 L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                  (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
))/(L1*pow(L2,2)*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/
        (2.*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))))/L2 + 
  (R.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*
     ((2*a1.diff(d100,i,j,k)*L2*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k))/
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (4*ar.diff(d010,i,j,k)*L2*pow(1 - rgrid[i],2)*Qr2(i,j,k))/
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (L2*(-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
          Qr2(i,j,k)*Qt1(i,j,k))/
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (16*a1.diff(d010,i,j,k)*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
          pow(1 - rgrid[i],2)*
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
       (16*a2.diff(d001,i,j,k)*pow(rh,2)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - rgrid[i],2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
        (L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (4*a0.diff(d010,i,j,k)*pow(rh,2)*pow(1 - rgrid[i],2)*
          ((L1*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
               pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
               (-2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                  (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) + 
                 2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                  Qtr(i,j,k)))/(2.*pow(rh,2)) - 
            4*L1*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))/
        (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a2.diff(d010,i,j,k)*
          ((-16*pow(rh,2)*pow(1 - rgrid[i],2)*
               (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
             (pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
            (4*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
               pow(1 - rgrid[i],2)*pow(Qr2(i,j,k),2))/
             ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),
                  2)))))/(L1*pow(1 + pow(R0,2),0.25)) - 
       (a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
          ((16*pow(rh,2)*pow(1 - rgrid[i],2)*
               (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
             (pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
            (4*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
               pow(1 - rgrid[i],2)*pow(Qr2(i,j,k),2))/
             ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),
                  2)))))/L2 + (a0.diff(d001,i,j,k)*
          (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qt1(i,j,k)*
          ((16*pow(rh,2)*pow(1 - rgrid[i],2)*
               (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
             (pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
            (4*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
               pow(1 - rgrid[i],2)*pow(Qr2(i,j,k),2))/
             ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),
                  2)))))/L2 - (ar(i,j,k)*(1 - rgrid[i])*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
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
             (pow(L1,2)*pow(L2,2)*
                pow(1 - pow(1 - rgrid[i],2),2)*pow(Qt1(i,j,k),2)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
               pow(L1,2)*pow(L2,2)*
                pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(1 - pow(1 - rgrid[i],2),2)*pow(Qt2(i,j,k),2)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
               2*pow(L1,2)*pow(L2,2)*
                (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                pow(1 - pow(1 - rgrid[i],2),2)*Qt1(i,j,k)*Qt2(i,j,k)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) \
+ 4*(1 - rgrid[i])*(pow(L2,2)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                (-(pow(L1,2)*pow(1 - pow(1 - rgrid[i],2),2)*
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
                  pow(L1,2)*sqrt(1 + pow(R0,2))) + 
               pow(L1,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                (-(pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
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
                  pow(L2,2)*sqrt(1 + pow(R0,2))) + 
               2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (-(L1*L2*R0) + 
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
                     Qt1(i,j,k)*Qt2(i,j,k))/(4.*pow(rh,4)))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
        (pow(L1,2)*pow(L2,2)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))))/L2 \
+ (ar.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
     (-((pow(L2,2)*pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
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
                    (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      pow(rh,2)*
                       (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
                  (2.*pow(rh,4)))*pow(Qt1(i,j,k),2) + 
              4*(1 - rgrid[i])*
               (-(pow(L1,2)*pow(1 - pow(1 - rgrid[i],2),2)*
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
                 pow(L1,2)*sqrt(1 + pow(R0,2)))))/pow(rh,4)) - 
       (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
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
          (pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),3)*
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
            4*(1 - rgrid[i])*(-(pow(L2,2)*
                   pow(1 - pow(1 - rgrid[i],2),2)*
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
               pow(L2,2)*sqrt(1 + pow(R0,2)))))/pow(rh,4) - 
       L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (-(L1*L2*(1 - rgrid[i])*pow(-1 + rgrid[i],6)*
              pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   4*pow(rh,2)*
                    (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                      4*pow(rgrid[i],3) + pow(rgrid[i],4))),3)*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
              pow(2*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                 (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) + 
                2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                 (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k),
               2)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (16.*pow(rh,10)) + 
          4*L1*L2*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
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
           (-((1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
             (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          (2*pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (L1*L2*pow(1 - pow(1 - rgrid[i],2),3)*
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
               2*L1*L2*(1 - rgrid[i])*
                ((1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                   (-((1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
                     (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
                  2*R0*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                   (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           pow(rh,4) + (pow(-1 + rgrid[i],4)*
             pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (L1*L2*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
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
                pow(2*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                   (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) + 
                  2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                   (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                  2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                   Qtr(i,j,k),2)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
               32*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(1 - rgrid[i],3)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
                Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
               4*L1*L2*pow(rh,2)*(1 - rgrid[i])*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                (8*L1*L2*R0*pow(1 - pow(1 - rgrid[i],2),2)*
                   pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qr2(i,j,k)*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
) + 4*pow(L1,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                   pow(1 - rgrid[i],2)*pow(Qr1(i,j,k),2)*
                   sqrt((1 + pow(R0,2))*
                     (1 + pow(R0,2) + 
                       2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                       pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2))) + 
                  4*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                   pow(1 - rgrid[i],2)*pow(Qr2(i,j,k),2)*
                   sqrt((1 + pow(R0,2))*
                     (1 + pow(R0,2) + 
                       2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                       pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2))))))/(16.*pow(rh,8)))))/
   (pow(L1,2)*pow(L2,3)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
  (Q.diff(d010,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*
     (-((L2*(-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
            pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
            Qr1(i,j,k)*Qt2(i,j,k))/
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
       (4*a0.diff(d010,i,j,k)*L2*pow(rh,2)*
          (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qt2(i,j,k)*
          (4*pow(1 - rgrid[i],2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
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
             pow(rh,2)))/
        (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (4*a0.diff(d001,i,j,k)*pow(rh,2)*pow(1 - rgrid[i],2)*
          (-4*L1*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            (pow(L1,2)*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
               pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
               Qr2(i,j,k)*Qt2(i,j,k)*
               sqrt(1 + pow(R0,2) + 
                 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                 pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
             pow(rh,2)))/
        (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (16*a2.diff(d010,i,j,k)*L2*pow(rh,2)*pow(1 - rgrid[i],2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (pow(L1,2)*pow(1 + pow(R0,2),0.25)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
       (16*a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
          pow(1 - rgrid[i],2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))))/L2 + 
  a2.diff(d020,i,j,k)*(4*pow(1 - pow(1 - rgrid[i],2),2)*
      pow(1 - rgrid[i],2)*pow(Qr1(i,j,k),2) + 
     (16*pow(rh,2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  (a0.diff(d001,i,j,k)*Qt1.diff(d010,i,j,k)*L1*
     pow(1 + pow(R0,2),0.25)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*(4*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*pow(Qr1(i,j,k),2) + 
       (16*pow(rh,2)*pow(1 - rgrid[i],2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/L2 + 
  a2.diff(d002,i,j,k)*(4*pow(1 - pow(1 - rgrid[i],2),2)*
      pow(1 - rgrid[i],2)*pow(Qr2(i,j,k),2) + 
     (16*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  (L1*pow(1 + pow(R0,2),0.25)*
     (Qt1.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)) + 
       2*(1 - rgrid[i])*Qt1(i,j,k))*
     (-(a2.diff(d100,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
        (4.*L1*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
       (ar.diff(d001,i,j,k)*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
        (2.*L1*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
       (a2.diff(d010,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (4.*pow(L1,2)*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (4.*L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       ((-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (-(L1*L2*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
               (1 + pow(R0,2) + 
                 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                 pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) \
+ L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (8.*pow(L1,2)*pow(rh,4)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (a0.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*
          pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr1(i,j,k)*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (4.*pow(L1,2)*pow(rh,4)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a0.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*
          (-2*pow(L1,2)*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            (L2*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
               (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                  (1 - rgrid[i])*Qr2(i,j,k)*Qt1(i,j,k)*
                  (1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*
                     pow(R(i,j,k),2)) + 
                 L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                  (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                     (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) - 
                    2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                     Qtr(i,j,k))*
                  (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
)))/(4.*pow(rh,4))))/
        (pow(L1,2)*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a1.diff(d100,i,j,k)*L2*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (4.*L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
       (ar.diff(d010,i,j,k)*L2*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (2.*pow(L1,2)*pow(rh,4)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/L2 + 
  (2*Qtr.diff(d010,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
     ((a2.diff(d100,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
        (4.*L1*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
       (ar.diff(d001,i,j,k)*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
        (2.*L1*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
       (a2.diff(d010,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (4.*pow(L1,2)*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (4.*L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       ((-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (L1*L2*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (8.*pow(L1,2)*pow(rh,4)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a0.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*
          pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr1(i,j,k)*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (4.*pow(L1,2)*pow(rh,4)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a0.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*
          (-2*pow(L1,2)*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            (L2*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
               (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                  (1 - rgrid[i])*Qr2(i,j,k)*Qt1(i,j,k)*
                  (1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*
                     pow(R(i,j,k),2)) + 
                 L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                  (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                     (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) - 
                    2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                     Qtr(i,j,k))*
                  (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
)))/(4.*pow(rh,4))))/
        (pow(L1,2)*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (a1.diff(d100,i,j,k)*L2*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (4.*L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
       (ar.diff(d010,i,j,k)*L2*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (2.*pow(L1,2)*pow(rh,4)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/L2 + 
  (2*Qr2.diff(d010,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
     ((a1.diff(d100,i,j,k)*pow(L2,2)*pow(1 + pow(R0,2),0.25)*
          (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) - 
       (2*ar.diff(d010,i,j,k)*pow(L2,2)*(1 - rgrid[i])*
          (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          (1 - pow(1 - rgrid[i],2))) + 
       (a0.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*
          (2*L1*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*Qt1(i,j,k)*
             pow(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2),1.5) \
+ L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                Qr1(i,j,k)*Qt1(i,j,k) - 
               2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k))*
             (R0 + pow(R0,3) + 
               (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*
                R(i,j,k) + 3*R0*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(R(i,j,k),2) + 
               pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)) + 
            2*pow(L1,2)*L2*pow(1 + 
               B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
             pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr1(i,j,k)*Qt2(i,j,k)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a2.diff(d010,i,j,k)*pow(L2,2)*
          (2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (pow(L1,2)*pow(1 + pow(R0,2),0.25)*
          pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) - 
       (a1.diff(d001,i,j,k)*L2*pow(1 + pow(R0,2),0.25)*
          (2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) + 
       (L2*(-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - rgrid[i])*(-(L1*L2*(1 - pow(1 - rgrid[i],2))*
               Qt1(i,j,k)*(1 + pow(R0,2) + 
                 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                 pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) \
+ L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (2.*pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) \
- (a0.diff(d010,i,j,k)*L2*pow(1 - rgrid[i],2)*
          (2*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr1(i,j,k)*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) - 
       (a2.diff(d100,i,j,k)*pow(L2,2)*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L1*pow(1 + pow(R0,2),0.25)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))) + 
       (2*ar.diff(d001,i,j,k)*L2*(1 - rgrid[i])*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 - pow(1 - rgrid[i],2)))))/L2 + 
  (2*Qr1.diff(d010,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
     ((a2.diff(d100,i,j,k)*L2*
          (-1 - pow(R0,2) - 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) - 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        pow(1 + pow(R0,2),0.25) + 
       (2*ar.diff(d001,i,j,k)*(1 - rgrid[i])*
          (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        (1 - pow(1 - rgrid[i],2)) - 
       ((-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - rgrid[i])*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          (L1*L2*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (2.*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (a0.diff(d010,i,j,k)*pow(1 - rgrid[i],2)*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          (2*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr1(i,j,k)*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
          (-2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (R0 + pow(R0,3) + 
               (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*
                R(i,j,k) + 3*R0*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(R(i,j,k),2) + 
               pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)) - 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (2 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a2.diff(d010,i,j,k)*L2*
          (2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (R0 + pow(R0,3) + 
               (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*
                R(i,j,k) + 3*R0*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(R(i,j,k),2) + 
               pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (2 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (L1*pow(1 + pow(R0,2),0.25)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a0.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*
          (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (-2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                Qtr(i,j,k)*(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
) + 2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*Qr1(i,j,k)*
                Qt1(i,j,k)*(2 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))\
)*sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (L1*L2*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
                (R0 + pow(R0,3) + 
                  (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*
                   R(i,j,k) + 
                  3*R0*pow(1 - pow(1 - rgrid[i],2),2)*
                   pow(R(i,j,k),2) + 
                  pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)) \
+ L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))\
)/(L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a1.diff(d100,i,j,k)*L2*pow(1 + pow(R0,2),0.25)*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) - 
       (2*ar.diff(d010,i,j,k)*L2*(1 - rgrid[i])*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 - pow(1 - rgrid[i],2)))))/L2 + 
  (2*Qr2.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
     ((4*a2.diff(d001,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*Qr2(i,j,k))/pow(1 + pow(R0,2),0.25) + 
       (2*ar.diff(d001,i,j,k)*(1 - rgrid[i])*
          (-1 - pow(R0,2) - 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) - 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        (1 - pow(1 - rgrid[i],2)) + 
       (a2.diff(d100,i,j,k)*L2*
          (-1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        pow(1 + pow(R0,2),0.25) + 
       (4*pow(rh,4)*(-((L2*a2(i,j,k)*(1 - rgrid[i])*
                 pow(-1 + rgrid[i],2)*
                 (pow(H,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   pow(rh,2)*
                    (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
               (pow(1 + pow(R0,2),0.25)*pow(rh,4))) + 
            L2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
             pow(1 - rgrid[i],2)*
             (((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
                pow(rh,4) - 
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
                   (2.*pow(rh,4))))*Qt2(i,j,k) + 
            (ar(i,j,k)*pow(-1 + rgrid[i],4)*
               pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
               (-(L2*pow(1 - pow(1 - rgrid[i],2),2)*
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
                    (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                       (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) + 
                      2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                       (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                      2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                       Qtr(i,j,k))) - 
                 4*(1 - rgrid[i])*
                  (2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                     Qr1(i,j,k)*
                     (L1*L2*R0 - 
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
                        Qt1(i,j,k)*Qt2(i,j,k))/(4.*pow(rh,4))) + 
                    (L2*pow(1 - pow(1 - rgrid[i],2),2)*
                       (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                       (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                       Qt2(i,j,k)*Qtr(i,j,k))/(2.*pow(rh,4)) + 
                    2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                     Qr2(i,j,k)*
                     (-(pow(L2,2)*
                        pow(1 - pow(1 - rgrid[i],2),2)*
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
                       pow(L2,2)*sqrt(1 + pow(R0,2))))))/
             (32.*pow(rh,6)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]))\
))/((1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
       a1.diff(d001,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*
        (2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          (2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))))) + 
       a0.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*
        (-((2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                Qr1(i,j,k)*Qt1(i,j,k) - 
               2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
          (2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr2(i,j,k)*Qt1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))) + 
       ((-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - rgrid[i])*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          (L1*L2*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (2.*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a0.diff(d010,i,j,k)*pow(1 - rgrid[i],2)*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          (2*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr1(i,j,k)*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a2.diff(d010,i,j,k)*L2*
          (-2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (R0 + pow(R0,3) + 
               (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*
                R(i,j,k) + 3*R0*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(R(i,j,k),2) + 
               pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)) - 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (-1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (L1*pow(1 + pow(R0,2),0.25)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (a1.diff(d100,i,j,k)*L2*pow(1 + pow(R0,2),0.25)*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
       (2*ar.diff(d010,i,j,k)*L2*(1 - rgrid[i])*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 - pow(1 - rgrid[i],2)))))/L2 + 
  (Qrr.diff(d010,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*
     ((a2.diff(d100,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*Qr1(i,j,k))/
        (pow(1 + pow(R0,2),0.25)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) - 
       (2*ar.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*Qr1(i,j,k))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
       (2*a0.diff(d001,i,j,k)*pow(rh,2)*pow(1 - rgrid[i],2)*
          (4*L1*L2*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             ((L1*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
                  (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                     (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) + 
                    2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                     (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                    2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                     Qtr(i,j,k)))/(2.*pow(rh,2)) + 
               4*L2*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (pow(L1,2)*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a1.diff(d001,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*
          ((4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
               pow(Qr1(i,j,k),2))/
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
            (16*pow(rh,2)*pow(1 - rgrid[i],2)*
               sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),
                  2)))/
             (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/2. + 
       (a2.diff(d010,i,j,k)*L2*
          ((-4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
               pow(Qr1(i,j,k),2))/
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
            (16*pow(rh,2)*pow(1 - rgrid[i],2)*
               sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),
                  2)))/
             (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/
        (2.*pow(1 + pow(R0,2),0.25))))/L2 + 
  (2*Qtr.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
     ((pow(rh,4)*(-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
                (2.*pow(rh,4))))*
          (4*a0(i,j,k)*pow(1 - rgrid[i],4) + 
            (ar(i,j,k)*(1 - rgrid[i])*pow(-1 + rgrid[i],4)*
               pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
               (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  Qr1(i,j,k)*Qt1(i,j,k) + 
                 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                  (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                 2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)\
))/(8.*pow(rh,6)*(1 - pow(1 - rgrid[i],2)))))/
        ((1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
       (a1.diff(d100,i,j,k)*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
        (4.*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
       (ar.diff(d010,i,j,k)*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
        (2.*L1*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
       (a0.diff(d001,i,j,k)*(1 - pow(1 - rgrid[i],2))*
          pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr2(i,j,k)*Qt1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
            L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                Qr1(i,j,k)*Qt1(i,j,k) - 
               2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (4.*L1*L2*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) \
+ ((-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - rgrid[i])*(-1 + 
            ((1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
               (-(L1*L2*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
                    (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
                 L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                  (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
)))/(4.*L1*L2*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/2. \
+ a0.diff(d010,i,j,k)*pow(1 - rgrid[i],2)*
        (2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) - 
          ((1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (2*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
               2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                Qtr(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
               2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                Qr1(i,j,k)*Qt2(i,j,k)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (4.*L1*L2*pow(rh,4)*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) - 
       (a2.diff(d100,i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (4.*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
       (ar.diff(d001,i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (2.*L2*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
       (a2.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
))/(4.*L1*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
       (a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
)/(4.*L2*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/L2 + 
  pow(1 + pow(R0,2),0.25)*(Qt2.diff(d100,i,j,k)*
      (1 - pow(1 - rgrid[i],2)) + 2*(1 - rgrid[i])*Qt2(i,j,k))*
   (-(a1.diff(d100,i,j,k)*pow(1 + pow(R0,2),0.25)*
         (1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
         (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
           pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
               pow(rgrid[i],3) + 
              4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                 4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
         (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
         (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (4.*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (ar.diff(d010,i,j,k)*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (2.*L1*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (a0.diff(d001,i,j,k)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
           Qr2(i,j,k)*Qt1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
              Qr1(i,j,k)*Qt1(i,j,k) - 
             2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (4.*L1*L2*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     ((-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
        (1 - rgrid[i])*(-1 + ((1 - pow(1 - rgrid[i],2))*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (L1*L2*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
               L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           (4.*L1*L2*pow(rh,4)*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/2. + 
     a0.diff(d010,i,j,k)*pow(1 - rgrid[i],2)*
      (2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) + 
        ((1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (2*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
              (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
             2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
              Qr1(i,j,k)*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
         (4.*L1*L2*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) \
+ (a2.diff(d100,i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (4.*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (ar.diff(d001,i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (2.*L2*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (a2.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k) + 
          2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
           Qr2(i,j,k)*R(i,j,k) + 
          2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (4.*L1*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
        (1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k) + 
          2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
           Qr2(i,j,k)*R(i,j,k) + 
          2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (4.*L2*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  (Qt1.diff(d001,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*
     (-(a1.diff(d100,i,j,k)*pow(1 + pow(R0,2),0.25)*
           (1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           (2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (4.*L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (ar.diff(d010,i,j,k)*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (2.*pow(L1,2)*pow(rh,4)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (a2.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          ((pow(L1,2)*pow(1 + 
                 B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
               pow(1 - pow(1 - rgrid[i],2),2)*
               pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
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
             pow(rh,2) + (pow(L2,2)*
               pow(1 - pow(1 - rgrid[i],2),2)*
               pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
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
             pow(rh,2) + (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (4*pow(1 - rgrid[i],2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
               (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                  pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
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
                pow(rh,2))*sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (pow(L1,2)*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          ((pow(L1,2)*pow(1 + 
                 B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
               pow(1 - pow(1 - rgrid[i],2),2)*
               pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
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
             pow(rh,2) + (pow(L2,2)*
               pow(1 - pow(1 - rgrid[i],2),2)*
               pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
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
             pow(rh,2) + (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (4*pow(1 - rgrid[i],2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
               (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                  pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
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
                pow(rh,2))*sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (L1*L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (a0.diff(d001,i,j,k)*(1 - pow(1 - rgrid[i],2))*
          pow(1 - rgrid[i],2)*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          ((L1*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),3)*
               pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
               pow(Qr2(i,j,k),2)*Qt1(i,j,k)*
               (1 + pow(R0,2) + 
                 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                 pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
             pow(rh,2) + (pow(L1,2)*
               pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
               (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
               pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
               (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                  (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) - 
                 2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                  Qtr(i,j,k))*
               (1 + pow(R0,2) + 
                 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                 pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
             (2.*pow(rh,2)) + 
            (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (4*L1*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) \
+ (L1*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
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
                  Qr2(i,j,k)*
                  (4*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                     (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) - 
                    2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                     Qtr(i,j,k))*
                  (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                (2.*pow(rh,2)))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (pow(L1,2)*L2*pow(rh,2)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (4*a0.diff(d010,i,j,k)*pow(rh,2)*pow(1 - rgrid[i],2)*
          ((pow(L1,2)*pow(L2,2)*
               pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
               pow(1 - pow(1 - rgrid[i],2),4)*
               pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],4)*
               pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
               pow(Qr1(i,j,k),2)*pow(Qt2(i,j,k),2)*
               (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
               (1 + pow(R0,2) + 
                 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                 pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
             (4.*pow(rh,6)) + 
            pow(L2,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             (-4*pow(1 - rgrid[i],2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
               (L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  pow(-1 + rgrid[i],4)*
                  pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                  Qr2(i,j,k)*Qt2(i,j,k)*
                  (2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                     (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                    2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                     Qtr(i,j,k))*
                  (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
                (8.*pow(rh,6)))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (-4*pow(L1,2)*pow(L2,2)*pow(rh,2)*
                  pow(1 - pow(1 - rgrid[i],2),2)*
                  pow(1 - rgrid[i],2)*
                  pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
                  pow(Qr1(i,j,k),2) + 
                 4*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                  pow(1 - rgrid[i],2)*
                  (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                  pow(Qt2(i,j,k),2)*
                  (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
                 (L1*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                    (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                    (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      pow(rh,2)*
                       (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                    (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                    Qr1(i,j,k)*Qt2(i,j,k)*
                    (4*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                       (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                      2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                       Qtr(i,j,k))*
                    (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                    (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                  (2.*pow(rh,2)))*
               sqrt(1 + pow(R0,2) + 
                 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                 pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
             (4.*pow(rh,4))))/
        (pow(L1,2)*pow(L2,2)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (pow(rh,4)*(-(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
                (2.*pow(rh,4))))*
          (8*pow(L1,2)*L2*a0(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],5)*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            (ar(i,j,k)*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (4*L1*L2*(1 - pow(1 - rgrid[i],2))*
                  pow(1 - rgrid[i],2)*
                  (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                  Qt1(i,j,k)*
                  (1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*
                     pow(R(i,j,k),2)) - 
                 L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                  (-(L1*L2*(1 - pow(1 - rgrid[i],2))*
                        (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
                        (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                       (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) + 
                        2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                       (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                        2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                        Qtr(i,j,k)))/(2.*pow(rh,2)) + 
                    4*L2*(1 - pow(1 - rgrid[i],2))*
                     pow(1 - rgrid[i],2)*
                     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                     Qt2(i,j,k)*
                     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
)))/(2.*pow(rh,4)*(1 - pow(1 - rgrid[i],2)))))/
        (pow(L1,2)*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       ((-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - rgrid[i])*(2*pow(L1,2)*pow(L2,2)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            (L2*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
               (-2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                  (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                  Qr1(i,j,k)*
                  (L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                     (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
                     (1 + pow(R0,2) + 
                       2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                       pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2)) - 
                    L1*L2*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
                     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                     sqrt(1 + pow(R0,2) + 
                       2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                       pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2))) + 
                 2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                  Qr2(i,j,k)*
                  (L1*L2*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
                     (1 + pow(R0,2) + 
                       2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                       pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2)) - 
                    L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                     (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
                     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                     sqrt(1 + pow(R0,2) + 
                       2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                       pow(1 - pow(1 - rgrid[i],2),2)*
                        pow(R(i,j,k),2)))))/(4.*pow(rh,4))))/
        (2.*pow(L1,2)*pow(L2,2)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a2.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
))/(4.*L1*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
       (ar.diff(d001,i,j,k)*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
)/(2.*L1*L2*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/L2 + 
  Qt2.diff(d010,i,j,k)*pow(1 + pow(R0,2),0.25)*
   (1 - pow(1 - rgrid[i],2))*
   ((a1.diff(d100,i,j,k)*pow(1 + pow(R0,2),0.25)*
        (1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (4.*L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (ar.diff(d010,i,j,k)*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (2.*pow(L1,2)*pow(rh,4)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (a2.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        ((pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr1(i,j,k),2)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           pow(rh,2) + (pow(L2,2)*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           pow(rh,2) + (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (4*pow(1 - rgrid[i],2)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
             (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              pow(rh,2))*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(L1,2)*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
        (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        ((pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr1(i,j,k),2)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           pow(rh,2) + (pow(L2,2)*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           pow(rh,2) + (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (4*pow(1 - rgrid[i],2)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
             (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                Qr2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              pow(rh,2))*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*L2*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (a0.diff(d001,i,j,k)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        ((L1*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),3)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2)*Qt1(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           pow(rh,2) + (pow(L1,2)*
             pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                Qr1(i,j,k)*Qt1(i,j,k) - 
               2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (2.*pow(rh,2)) + 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (4*L1*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) + 
             (L1*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
                (4*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                   (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) - 
                  2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                   Qtr(i,j,k))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
              (2.*pow(rh,2)))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (pow(L1,2)*L2*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (4*a0.diff(d010,i,j,k)*pow(rh,2)*pow(1 - rgrid[i],2)*
        (-(pow(L1,2)*pow(L2,2)*
              pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
              pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
              pow(-1 + rgrid[i],4)*
              pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   4*pow(rh,2)*
                    (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                      4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
              (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
              pow(Qr1(i,j,k),2)*pow(Qt2(i,j,k),2)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (4.*pow(rh,6)) - 
          pow(L2,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           (4*pow(1 - rgrid[i],2)*
              (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
             (L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                pow(-1 + rgrid[i],4)*
                pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
                Qr2(i,j,k)*Qt2(i,j,k)*
                (2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                   (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                  2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                   Qtr(i,j,k))*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
              (8.*pow(rh,6)))*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (4*pow(L1,2)*pow(L2,2)*pow(rh,2)*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(1 - rgrid[i],2)*
                pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
                pow(Qr1(i,j,k),2) + 
               4*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(1 - rgrid[i],2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                pow(Qt2(i,j,k),2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
               (L1*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                  (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                        pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                  Qr1(i,j,k)*Qt2(i,j,k)*
                  (4*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                     (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                    2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                     Qtr(i,j,k))*
                  (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
                  (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
                (2.*pow(rh,2)))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (4.*pow(rh,4))))/
      (pow(L1,2)*pow(L2,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     ((-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
        (1 - rgrid[i])*(2*pow(L1,2)*pow(L2,2)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          (L2*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
             (2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
                (L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                   (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
                   (1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) - 
                  L1*L2*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2))) - 
               2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                Qr2(i,j,k)*(L1*L2*(1 - pow(1 - rgrid[i],2))*
                   Qt1(i,j,k)*
                   (1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*
                      pow(R(i,j,k),2)) - 
                  L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                   (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
))))/(4.*pow(rh,4))))/
      (2.*pow(L1,2)*pow(L2,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (a2.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k) + 
          2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
           Qr2(i,j,k)*R(i,j,k) + 
          2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (4.*L1*pow(1 + pow(R0,2),0.25)*pow(rh,4)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     (ar.diff(d001,i,j,k)*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k) + 
          2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
           Qr2(i,j,k)*R(i,j,k) + 
          2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (2.*L1*L2*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
  (2*Qr1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
     (4*a1.diff(d010,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*
        (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k) + 
       a1.diff(d100,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*
        (-1 - pow(R0,2) - 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) - 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
       (2*ar.diff(d010,i,j,k)*(1 - rgrid[i])*
          (-1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        (1 - pow(1 - rgrid[i],2)) + 
       (4*pow(rh,4)*(-((L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
                 (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                 (pow(H,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   pow(rh,2)*
                    (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
               pow(rh,4)) + 
            L1*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
             pow(1 - rgrid[i],2)*
             (((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
                pow(rh,4) - 
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
                   (2.*pow(rh,4))))*Qt1(i,j,k) + 
            (ar(i,j,k)*pow(-1 + rgrid[i],4)*
               pow(pow(H,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
               (-(L1*pow(1 - pow(1 - rgrid[i],2),2)*
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
                    (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                       (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) + 
                      2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                       (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                      2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                       Qtr(i,j,k))) - 
                 4*(1 - rgrid[i])*
                  (2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                     Qr2(i,j,k)*
                     (L1*L2*R0 - 
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
                        Qt1(i,j,k)*Qt2(i,j,k))/(4.*pow(rh,4))) + 
                    (L1*pow(1 - pow(1 - rgrid[i],2),2)*
                       (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
                       (pow(H,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        pow(rh,2)*
                        (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                        4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                       Qt1(i,j,k)*Qtr(i,j,k))/(2.*pow(rh,4)) + 
                    2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                     Qr1(i,j,k)*
                     (-(pow(L1,2)*
                        pow(1 - pow(1 - rgrid[i],2),2)*
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
                       pow(L1,2)*sqrt(1 + pow(R0,2))))))/
             (32.*pow(rh,6)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]))\
))/((1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))) - 
       (a0.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr2(i,j,k)*Qt1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
            L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                Qr1(i,j,k)*Qt1(i,j,k) - 
               2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/L2 \
+ (a2.diff(d010,i,j,k)*L2*(-2*(1 - pow(1 - rgrid[i],2))*
             (1 - rgrid[i])*Qr2(i,j,k)*
             (-2 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            (2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
               (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
               sqrt(1 + pow(R0,2) + 
                 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                 pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
             L2))/pow(1 + pow(R0,2),0.25) + 
       ((-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - rgrid[i])*(L1*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/2. \
+ a0.diff(d010,i,j,k)*pow(1 - rgrid[i],2)*
        ((2*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
              Qr2(i,j,k)*Qt2(i,j,k) - 
             2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k))*
           (-1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          (2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (-2*L1*L2*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k) + 
               L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
           L2) + (a2.diff(d100,i,j,k)*L1*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        pow(1 + pow(R0,2),0.25) - 
       (2*ar.diff(d001,i,j,k)*L1*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*(1 - rgrid[i])*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L2*(1 - pow(1 - rgrid[i],2))) + 
       (a1.diff(d001,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
)/L2))/L2 + (H*pow(1 + pow(R0,2),0.25)*
     (2*L1*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        (-(1/(1 - rgrid[i])) - 
          (4*pow(rh,4)*((pow(-1 + rgrid[i],2)*
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
                (2.*pow(rh,4))))/
           (pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))))*
        Qr1(i,j,k) - L1*L2*(2*Qr1.diff(d100,i,j,k)*
           (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]) + 
          (-2*(1 - pow(1 - rgrid[i],2)) + 4*pow(1 - rgrid[i],2))*
           Qr1(i,j,k)) + (L1*L2*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*Qr1(i,j,k)*
          (Qrr.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)) + 
            2*(1 - rgrid[i])*Qrr(i,j,k)))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
       (L1*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
          (Qtt.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)) + 
            2*(1 - rgrid[i])*Qtt(i,j,k)))/
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
       (16*Q.diff(d001,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
          pow(1 - rgrid[i],2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
        (pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
       Qrr.diff(d001,i,j,k)*(1 - pow(1 - rgrid[i],2))*
        ((-2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qr2(i,j,k))/
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
          (8*pow(rh,2)*pow(1 - rgrid[i],2)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
           (pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
       (Qtt.diff(d001,i,j,k)*(1 - pow(1 - rgrid[i],2))*
          (4*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qr2(i,j,k) - 
            (16*pow(rh,2)*pow(1 - rgrid[i],2)*
               (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
               (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
             (pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/
        (2.*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
       (4*R.diff(d010,i,j,k)*pow(L2,2)*
          pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
          Qr1(i,j,k)*Qr2(i,j,k))/
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (4*R.diff(d001,i,j,k)*pow(L2,2)*
          pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
          pow(Qr2(i,j,k),2))/
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (2*pow(L2,2)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
          Qr2(i,j,k)*(R.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)) + 
            2*(1 - rgrid[i])*R(i,j,k)))/
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (2*Qtt.diff(d010,i,j,k)*L2*pow(rh,2)*
          (1 - pow(1 - rgrid[i],2))*
          (4*pow(1 - rgrid[i],2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
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
             pow(rh,2)))/
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
       (2*Qr2.diff(d001,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          (2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
            2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) - 
       2*Qr1.diff(d001,i,j,k)*L1*(1 - pow(1 - rgrid[i],2))*
        (1 - rgrid[i])*(2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
           Qr2(i,j,k)*(-1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
       (2*Qr2.diff(d010,i,j,k)*pow(L2,2)*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*(2*L2*(1 - pow(1 - rgrid[i],2))*
             (1 - rgrid[i])*Qr2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) - 
       (Qtr.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
          (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (2.*L1*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (L2*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (Qt1.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)) + 
            2*(1 - rgrid[i])*Qt1(i,j,k))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (4.*pow(rh,4)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (2*Qr1.diff(d010,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*(2*L2*(1 - pow(1 - rgrid[i],2))*
             (1 - rgrid[i])*Qr2(i,j,k)*
             (R0 + pow(R0,3) + 
               (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*
                R(i,j,k) + 3*R0*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(R(i,j,k),2) + 
               pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             (2 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (Qt1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
          Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          ((pow(L1,2)*pow(1 + 
                 B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
               pow(1 - pow(1 - rgrid[i],2),2)*
               pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
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
             pow(rh,2) + (pow(L2,2)*
               pow(1 - pow(1 - rgrid[i],2),2)*
               pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
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
             pow(rh,2) + (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (4*pow(1 - rgrid[i],2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
               (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                  pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
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
                pow(rh,2))*sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (Qt2.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
          Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          ((pow(L1,2)*pow(1 + 
                 B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
               pow(1 - pow(1 - rgrid[i],2),2)*
               pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
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
             pow(rh,2) + (pow(L2,2)*
               pow(1 - pow(1 - rgrid[i],2),2)*
               pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
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
             pow(rh,2) + (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (4*pow(1 - rgrid[i],2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
               (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                  pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
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
                pow(rh,2))*sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (L1*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (16*Q.diff(d010,i,j,k)*L2*pow(rh,2)*
          (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
       (Qrr.diff(d010,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          ((-4*pow(L1,2)*pow(1 - pow(1 - rgrid[i],2),2)*
               pow(1 - rgrid[i],2)*pow(Qr1(i,j,k),2))/
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
            (16*pow(rh,2)*pow(1 - rgrid[i],2)*
               sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),
                  2)))/
             ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))))/(2.*L1) + 
       (Qtr.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
          (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
))/(2.*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
       (L2*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*
                (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*Qt2(i,j,k)*
          (Qt2.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)) + 
            2*(1 - rgrid[i])*Qt2(i,j,k))*
          (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
))/(4.*pow(rh,4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
       (L2*(B.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)) + 
            2*B(i,j,k)*(1 - rgrid[i]))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
))/pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) + 
       (2*B.diff(d010,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
          (1 - rgrid[i])*Qr1(i,j,k)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
))/pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) + 
       (2*B.diff(d001,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
          (1 - rgrid[i])*Qr2(i,j,k)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
)/pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)))/L2 + 
  (Q.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*
     (-((L2*(-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
            pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
            Qr2(i,j,k)*Qt2(i,j,k))/
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
       (32*a2.diff(d010,i,j,k)*pow(rh,2)*pow(1 - rgrid[i],2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
        (L1*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
       (4*a0.diff(d010,i,j,k)*pow(rh,2)*pow(1 - rgrid[i],2)*
          (4*L1*L2*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
            L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
             ((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                  pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
                  (pow(H,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    pow(rh,2)*
                     (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                       4*pow(rh,2)*
                        (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                  (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                  Qr1(i,j,k)*Qr2(i,j,k))/pow(rh,2) - 
               8*pow(1 - rgrid[i],2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
        (pow(L1,2)*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (ar(i,j,k)*(1 - rgrid[i])*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
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
             (pow(L1,2)*pow(L2,2)*
                pow(1 - pow(1 - rgrid[i],2),2)*pow(Qt1(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
+ pow(L1,2)*pow(L2,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(1 - pow(1 - rgrid[i],2),2)*pow(Qt2(i,j,k),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
- 2*pow(L1,2)*pow(L2,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                pow(1 - pow(1 - rgrid[i],2),2)*Qt1(i,j,k)*
                Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) \
+ 4*(1 - rgrid[i])*(pow(L2,2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
                (-(pow(L1,2)*pow(1 - pow(1 - rgrid[i],2),2)*
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
                  pow(L1,2)*sqrt(1 + pow(R0,2))) + 
               pow(L1,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                (1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
                (-(pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
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
                  pow(L2,2)*sqrt(1 + pow(R0,2))) - 
               2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                (L1*L2*R0 - (L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
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
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))\
)/(pow(L1,2)*pow(L2,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - pow(1 - rgrid[i],2),2)*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (16*a1.diff(d010,i,j,k)*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
          pow(1 - rgrid[i],2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
       (16*a2.diff(d001,i,j,k)*pow(rh,2)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          pow(1 - rgrid[i],2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
       (a0.diff(d001,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          pow(1 - rgrid[i],2)*Qt2(i,j,k)*
          (4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2) + 
            (16*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
               pow(1 - rgrid[i],2)*
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
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)))/L2 + 
  Qt2.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
   (1 - pow(1 - rgrid[i],2))*
   (2*(-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
      (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qr2(i,j,k) + 
     a0.diff(d010,i,j,k)*pow(1 - rgrid[i],2)*
      (-8*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
         Qr1(i,j,k)*Qr2(i,j,k) + 
        (32*pow(rh,2)*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
         (L1*L2*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))) + 
     (pow(rh,4)*(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
        (8*L1*pow(L2,2)*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
           pow(1 - rgrid[i],5)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k) + 
          (ar(i,j,k)*(1 - rgrid[i])*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             ((L1*pow(L2,2)*(1 - pow(1 - rgrid[i],2))*
                  (1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
                  (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*
                     (1 - rgrid[i])*Qr1(i,j,k)*Qt1(i,j,k) + 
                    2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                     (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k) - 
                    2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                     Qtr(i,j,k)))/(2.*pow(rh,2)) - 
               4*pow(1 - rgrid[i],2)*
                (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
                (L1*L2*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)*
                   (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
                  L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                   (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
                   sqrt(1 + pow(R0,2) + 
                     2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                     pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
))))/(2.*pow(rh,4)*(1 - pow(1 - rgrid[i],2)))))/
      (L1*pow(L2,2)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
     a0.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*
      (-4*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
         pow(Qr2(i,j,k),2) - 
        (16*pow(rh,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
         (pow(L2,2)*pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))))) + 
  (pow(1 + pow(R0,2),0.25)*
     (B.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)) + 
       2*B(i,j,k)*(1 - rgrid[i]))*
     ((a2.diff(d100,i,j,k)*L2*
          (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        (pow(1 + pow(R0,2),0.25)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))) - 
       (2*ar.diff(d001,i,j,k)*(1 - rgrid[i])*
          (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
          (1 - pow(1 - rgrid[i],2))) - 
       (a1.diff(d100,i,j,k)*L2*pow(1 + pow(R0,2),0.25)*
          (R0 + pow(R0,3) + 
            (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
            pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)))/
        (pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (2*ar.diff(d010,i,j,k)*L2*(1 - rgrid[i])*
          (R0 + pow(R0,3) + 
            (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
            pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)))/
        (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          (1 - pow(1 - rgrid[i],2))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (a0.diff(d010,i,j,k)*pow(1 - rgrid[i],2)*
          (2*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
            2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
            2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr1(i,j,k)*Qt2(i,j,k)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) - 
       (a0.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*
          (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr2(i,j,k)*Qt1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
            L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                Qr1(i,j,k)*Qt1(i,j,k) - 
               2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) - 
       (a2.diff(d010,i,j,k)*L2*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
))/(L1*pow(1 + pow(R0,2),0.25)*
          pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) + 
       (a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
))/pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) + 
       ((-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - rgrid[i])*sqrt(1 + 
            pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
          (L1*L2*R0*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k) + 
            L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*Qt1(i,j,k)*
             R(i,j,k) - L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
)/(2.*L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2))))/L2 + 
  (B.diff(d010,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*
     ((-2*a2.diff(d100,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*Qr1(i,j,k)*
          (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        (pow(1 + pow(R0,2),0.25)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))) + 
       (4*ar.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*Qr1(i,j,k)*
          (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
       (2*a1.diff(d100,i,j,k)*L2*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
          (R0 + pow(R0,3) + 
            (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
            pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)))/
        (pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (4*ar.diff(d010,i,j,k)*L2*pow(1 - rgrid[i],2)*Qr1(i,j,k)*
          (R0 + pow(R0,3) + 
            (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
            pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)))/
        (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (2*a0.diff(d001,i,j,k)*(1 - pow(1 - rgrid[i],2))*
          pow(1 - rgrid[i],3)*Qr1(i,j,k)*
          (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr2(i,j,k)*Qt1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
            L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                Qr1(i,j,k)*Qt1(i,j,k) - 
               2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) - 
       (4*a0.diff(d010,i,j,k)*pow(rh,2)*pow(1 - rgrid[i],2)*
          (4*L2*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k) + 
            (L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
               pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
               (2*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                  (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k)*
                  (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
                 2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                  Qtr(i,j,k)*
                  (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
                 2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                  pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  Qr1(i,j,k)*Qt2(i,j,k)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*
                     pow(R(i,j,k),2))))/(2.*pow(rh,2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
       (2*a2.diff(d010,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*Qr1(i,j,k)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
))/(L1*pow(1 + pow(R0,2),0.25)*
          pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) - 
       (2*a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
))/pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) + 
       ((-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qr1(i,j,k)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
          (-(L1*L2*R0*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)) - 
            L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*Qt1(i,j,k)*
             R(i,j,k) + L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
)/(L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2))))/L2 + 
  (B.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
     (1 - pow(1 - rgrid[i],2))*
     ((-2*a2.diff(d100,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*Qr2(i,j,k)*
          (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        (pow(1 + pow(R0,2),0.25)*
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))) + 
       (4*ar.diff(d001,i,j,k)*pow(1 - rgrid[i],2)*Qr2(i,j,k)*
          (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2))) + 
       (2*a1.diff(d100,i,j,k)*L2*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
          (R0 + pow(R0,3) + 
            (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
            pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)))/
        (pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
       (4*ar.diff(d010,i,j,k)*L2*pow(1 - rgrid[i],2)*Qr2(i,j,k)*
          (R0 + pow(R0,3) + 
            (1 + 3*pow(R0,2))*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            3*R0*pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2) + 
            pow(1 - pow(1 - rgrid[i],2),3)*pow(R(i,j,k),3)))/
        (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
       (16*a2.diff(d001,i,j,k)*pow(rh,2)*pow(1 - rgrid[i],2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
       (16*a1.diff(d010,i,j,k)*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
          pow(1 - rgrid[i],2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
       (ar(i,j,k)*(1 - rgrid[i])*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
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
             (pow(L1,2)*pow(L2,2)*
                pow(1 - pow(1 - rgrid[i],2),2)*pow(Qt1(i,j,k),2) \
- pow(L1,2)*pow(L2,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                pow(1 - pow(1 - rgrid[i],2),2)*pow(Qt2(i,j,k),2)) \
+ 4*(1 - rgrid[i])*(pow(L2,2)*
                (-(pow(L1,2)*pow(1 - pow(1 - rgrid[i],2),2)*
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
                  pow(L1,2)*sqrt(1 + pow(R0,2))) + 
               pow(L1,2)*pow(1 + 
                  B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
                ((pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
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
                  pow(L2,2)*sqrt(1 + pow(R0,2)))))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (pow(L1,2)*pow(L2,2)*
          pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          pow(1 - pow(1 - rgrid[i],2),2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
       (2*a0.diff(d001,i,j,k)*(1 - pow(1 - rgrid[i],2))*
          pow(1 - rgrid[i],3)*Qr2(i,j,k)*
          (2*L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
             Qr2(i,j,k)*Qt1(i,j,k)*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
            L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (2*L1*pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                Qr1(i,j,k)*Qt1(i,j,k) - 
               2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qtr(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) + 
       (4*a0.diff(d010,i,j,k)*pow(rh,2)*pow(1 - rgrid[i],2)*
          (4*L1*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) - 
            (L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
               pow(-1 + rgrid[i],2)*
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
               (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
               (2*pow(L2,2)*pow(1 - pow(1 - rgrid[i],2),2)*
                  (1 - rgrid[i])*Qr2(i,j,k)*Qt2(i,j,k)*
                  (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
                 2*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
                  Qtr(i,j,k)*
                  (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
                 2*L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                  pow(1 - pow(1 - rgrid[i],2),2)*(1 - rgrid[i])*
                  Qr1(i,j,k)*Qt2(i,j,k)*
                  sqrt(1 + pow(R0,2) + 
                    2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                    pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)\
)))/(2.*pow(rh,2)))*sqrt(1 + 
            pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
        (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
          (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
       (2*a2.diff(d010,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
          (1 - rgrid[i])*Qr2(i,j,k)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
)/(L1*pow(1 + pow(R0,2),0.25)*
          pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)) - 
       (2*a1.diff(d001,i,j,k)*pow(1 + pow(R0,2),0.25)*
          (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
          (2*L2*R0*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
             Qr2(i,j,k) + 2*L2*pow(1 - pow(1 - rgrid[i],2),2)*
             (1 - rgrid[i])*Qr2(i,j,k)*R(i,j,k) + 
            2*L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))\
)/pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2) + 
       ((-4*a0(i,j,k) + 2*a0.diff(d100,i,j,k)*(1 - rgrid[i]))*
          (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qr2(i,j,k)*
          sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
          (-(L1*L2*R0*(1 - pow(1 - rgrid[i],2))*Qt1(i,j,k)) - 
            L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*Qt1(i,j,k)*R(i,j,k) + 
            L1*L2*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*Qt2(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
        (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2))))/L2
;}
}


elem[58]+=
((2*a1.diff(d100,i,j,k)*pow(1 + pow(R0,2),0.25)*
       pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],3)*
       Qr2(i,j,k)*Qt1(i,j,k))/
     (pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
    (4*ar.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],4)*Qr2(i,j,k)*Qt1(i,j,k))/
     (pow(rh,2)*(L1 + L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (4*a1.diff(d001,i,j,k)*L1*pow(1 + pow(R0,2),0.25)*
       pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
       (-(pow(L2,2)*(1 - pow(1 - rgrid[i],2))*pow(Qr2(i,j,k),2)*
            Qt1(i,j,k)) + L1*pow(1 + 
            B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*Qr1(i,j,k)*
          (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
            Qtr(i,j,k))))/
     (pow(L2,2)*pow(rh,2)*(L1 + 
         L1*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
       (L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*pow(Qt1(i,j,k),2) + 
         pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*Qt2(i,j,k)*
          (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
            Qtr(i,j,k))))/
     (pow(rh,2)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
    (2*a0.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
       pow(1 - rgrid[i],5)*(L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*
          pow(Qt1(i,j,k),2) + 
         pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*Qt2(i,j,k)*
          (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - 
            Qtr(i,j,k))))/
     (pow(rh,2)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
    (4*ar.diff(d001,i,j,k)*pow(1 + 
         B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
       (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],4)*
       (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - Qtr(i,j,k))\
)/(pow(L2,2)*pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (2*a2.diff(d100,i,j,k)*pow(1 + 
         B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
       pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],3)*
       (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) - Qtr(i,j,k))\
)/(pow(1 + pow(R0,2),0.25)*pow(rh,2)*
       (L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (4*H*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
       (pow(L2,2)*(1 - pow(1 - rgrid[i],2))*pow(Qr2(i,j,k),2)*
          Qt1(i,j,k) + L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),
           2)*Qr1(i,j,k)*(-(L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
               Qt1(i,j,k)) + Qtr(i,j,k))))/
     (pow(rh,2)*(L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (4*a2.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
       pow(1 - rgrid[i],4)*(pow(L2,2)*(1 - pow(1 - rgrid[i],2))*
          pow(Qr2(i,j,k),2)*Qt1(i,j,k) + 
         L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*Qr1(i,j,k)*
          (-(L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k)) + 
            Qtr(i,j,k))))/
     (L1*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
       (L2 + L2*(1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (4*a0.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
       pow(1 - rgrid[i],6)*(L1*
          pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*Qr1(i,j,k)*
          Qt2(i,j,k)*(L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*
             Qt1(i,j,k) - Qtr(i,j,k)) + 
         L2*Qr2(i,j,k)*Qt1(i,j,k)*
          (-(L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k)) + 
            Qtr(i,j,k))))/
     (L1*L2*pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
    (64*a0.diff(d001,i,j,k)*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*
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
          (16.*pow(rh,6)) + pow(1 + 
            B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k) - 
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
             (16.*pow(rh,6)))))/
     (pow(L2,2)*pow(-1 + rgrid[i],4)*
       pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
       (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))*
  D[-1 + d001](i,j,k,i1,j1,k1)
;

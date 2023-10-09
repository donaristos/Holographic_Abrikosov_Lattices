
elem[16]+=
((-2*Qtt.diff(d100,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr1(i,j,k)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
     pow(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k),2) + 
    (4*Qtt.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
       pow(1 - rgrid[i],2)*pow(Qr1(i,j,k),2)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
     pow(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k),2) + 
    (4*Qtt.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
       pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qr2(i,j,k)*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
     pow(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k),2) + 
    (8*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
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
               (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                 pow(rh,2)*
                  (pow(mu,2)*pow(-2 + rgrid[i],3)*
                     pow(rgrid[i],3) + 
                    4*pow(rh,2)*
                     (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                       4*pow(rgrid[i],3) + pow(rgrid[i],4)))))/
             (2.*pow(rh,4)))*(Qrr(i,j,k) - Qtt(i,j,k)) + 
         ((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
            (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
              pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                  pow(rgrid[i],3) + 
                 4*pow(rh,2)*
                  (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                    4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
            (-2*Qrr(i,j,k) + Qtt(i,j,k)))/(2.*pow(rh,4))))/
     (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
          pow(rgrid[i],3) + pow(rh,2)*
          (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       pow(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k),2)))*
  D[-1 + d010](i,j,k,i1,j1,k1)
;

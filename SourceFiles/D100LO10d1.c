
elem[144]+=
((Qt1.diff(d100,i,j,k)*L1*pow(1 - pow(1 - rgrid[i],2),2)*Qr1(i,j,k))/
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
    (2*Qtr.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr1(i,j,k))/
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
    (Qt2.diff(d100,i,j,k)*L2*pow(1 - pow(1 - rgrid[i],2),2)*
       Qr2(i,j,k))/(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
    (2*Qtr.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
       (1 - rgrid[i])*Qr2(i,j,k))/
     (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
    ((4*pow(rh,4)*(((1 - rgrid[i])*pow(-1 + rgrid[i],2)*
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
          (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qtr(i,j,k))/
        (pow(-1 + rgrid[i],2)*
          (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                pow(rgrid[i],3) + 
               4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                  4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
       4*L1*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr1(i,j,k)*
        Qt1(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
       4*L2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*Qr2(i,j,k)*
        Qt2(i,j,k)*(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/
     (2.*pow(1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k),2)))*
  D[-1 + d100](i,j,k,i1,j1,k1)
;

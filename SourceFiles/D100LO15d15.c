
elem[238]+=
(1/(1 - rgrid[i]) - 2*Qr1.diff(d010,i,j,k)*(1 - pow(1 - rgrid[i],2))*
      (1 - rgrid[i]) - 2*Qr2.diff(d001,i,j,k)*
      (1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]) + 
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
      (pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))) + 
     (Q.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (2*(1 - rgrid[i])*Q(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
     (2*Q.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*Qr1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
     (2*Q.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*Qr2(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
     (Qrr.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)))/
      (2 + 2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     (2*Qrr.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*Qr1(i,j,k))/
      (2 + 2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     (2*Qrr.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*Qr2(i,j,k))/
      (2 + 2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
     (2*(1 - rgrid[i])*Qrr(i,j,k))/
      (2 + 2*(1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
     (Qtt.diff(d100,i,j,k)*(1 - pow(1 - rgrid[i],2)))/
      (2 + 2*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
     (2*Qtt.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*Qr1(i,j,k))/
      (2 + 2*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) - 
     (2*Qtt.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        (1 - rgrid[i])*Qr2(i,j,k))/
      (2 + 2*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)) + 
     (2*(1 - rgrid[i])*Qtt(i,j,k))/
      (2 + 2*(1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))*
   D[-1 + d100](i,j,k,i1,j1,k1) + D[-1 + d200](i,j,k,i1,j1,k1)
;

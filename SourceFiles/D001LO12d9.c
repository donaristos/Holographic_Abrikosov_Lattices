
elem[184]+=
(4*a0.diff(d010,i,j,k)*L2*(1 - pow(1 - rgrid[i],2))*
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
         (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))))*
    D[-1 + d001](i,j,k,i1,j1,k1))/(L1*pow(1 + pow(R0,2),0.25))
;

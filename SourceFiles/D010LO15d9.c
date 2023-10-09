
elem[232]+=
(16*q*pow(rh,2)*a0(i,j,k)*h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*
    pow(1 - rgrid[i],4)*((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
         pow(-1 + rgrid[i],2)*
         (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
           pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
               pow(rgrid[i],3) + 
              4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                 4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
         (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*Qr2(i,j,k)\
)/(4.*pow(rh,2)) - (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
    D[-1 + d010](i,j,k,i1,j1,k1))/
  (L1*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
       pow(rgrid[i],3) + pow(rh,2)*
       (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
            4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
    (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))
;

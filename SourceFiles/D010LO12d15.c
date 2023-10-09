
elem[190]+=
(16*q*pow(rh,4)*h2(i,j,k)*pow(1 - rgrid[i],2)*
    (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
    D[-1 + d010](i,j,k,i1,j1,k1))/
  (L1*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
    (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
      pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
            4*pow(rgrid[i],3) + pow(rgrid[i],4)))))
;

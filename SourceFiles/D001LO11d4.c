
elem[163]+=
(16*a0.diff(d001,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
    pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
    sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))*
    D[-1 + d001](i,j,k,i1,j1,k1))/
  (pow(L2,2)*pow(-1 + rgrid[i],2)*
    (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
      pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
            4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
    (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))
;

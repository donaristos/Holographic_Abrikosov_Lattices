
elem[14]+=
(-8*pow(rh,6)*h1(i,j,k)*pow(1 - rgrid[i],2)*
    (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
    (16*pow(q,2)*pow(a0(i,j,k),2)*pow(1 - pow(1 - rgrid[i],2),2)*
       pow(1 - rgrid[i],4) + 
      (4*(1 + 12*S*(pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
            pow(1 - pow(1 - rgrid[i],2),2))*pow(-1 + rgrid[i],2)*
         (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
           pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
               pow(rgrid[i],3) + 
              4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                 4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
         (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))/pow(rh,2))*
    sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
      pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
  ((1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],4)*
    pow(pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
      pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
            4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
    sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))
;

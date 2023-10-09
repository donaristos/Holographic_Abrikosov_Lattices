
elem[123]+=
(-64*pow(q,2)*pow(1 + pow(R0,2),0.25)*pow(rh,6)*a0(i,j,k)*
    (pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*(1 - pow(1 - rgrid[i],2))*
    pow(1 - rgrid[i],4)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
  (pow(-1 + rgrid[i],4)*pow(pow(H,2)*pow(-2 + rgrid[i],3)*
       pow(rgrid[i],3) + pow(rh,2)*
       (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
            4*pow(rgrid[i],3) + pow(rgrid[i],4))),2)*
    (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))
;

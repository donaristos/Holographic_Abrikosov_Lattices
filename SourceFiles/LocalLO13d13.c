
elem[204]+=
-8*Qr2.diff(d001,i,j,k)*pow(1 - rgrid[i],2) + 
  (8*Qrr.diff(d001,i,j,k)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*Qr2(i,j,k))/
   (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
  (16*pow(q,2)*pow(rh,4)*(pow(h1(i,j,k),2) + pow(h2(i,j,k),2))*
     pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4)))))
;

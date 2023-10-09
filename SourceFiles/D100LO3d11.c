
elem[42]+=
((-8*pow(rh,2)*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
       pow(1 - rgrid[i],3)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))/
     (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
          pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
    (4*a0.diff(d100,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
       pow(1 - rgrid[i],4)*(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)))/
     (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
          pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
    (8*a0.diff(d010,i,j,k)*pow(rh,2)*
       pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],5)*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k))/
     (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
          pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
    (8*a0.diff(d001,i,j,k)*pow(rh,2)*
       pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],5)*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k))/
     (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
          pow(rgrid[i],3) + pow(rh,2)*
          (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))))*
  D[-1 + d100](i,j,k,i1,j1,k1)
;

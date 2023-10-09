
elem[78]+=
((16*h1.diff(d001,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
       pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
         pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
     (L1*L2*pow(-1 + rgrid[i],2)*
       (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
    (4*q*pow(rh,2)*h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*
       pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       ((-4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
            (1 + pow(R0,2) + 
              2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
              pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
          pow(1 + pow(R0,2),0.25) + 
         4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
          (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
          sqrt(1 + pow(R0,2) + 
            2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
            pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
         4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
          ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
             (1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
            Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
     (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
       (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
    (16*h1.diff(d010,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
       pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
     (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
       pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
          pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
    (16*H*L2*q*pow(rh,2)*h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*
       pow(1 - rgrid[i],2)*ygrid[k]*
       (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
       (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
       sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
     (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
       (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
         pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
             pow(rgrid[i],3) + 
            4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
               4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
       (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))))*
  D[-1 + d010](i,j,k,i1,j1,k1)
;

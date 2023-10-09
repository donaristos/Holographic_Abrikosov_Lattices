
elem[226]+=
(-4*h1(i,j,k)*pow(1 - rgrid[i],2)*Q(i,j,k))/
   pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) + 
  (4*h1(i,j,k)*pow(1 - rgrid[i],2))/
   (1 - pow(1 - rgrid[i],2) + pow(1 - pow(1 - rgrid[i],2),2)*Q(i,j,k)) \
+ h1.diff(d100,i,j,k)*((-2*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        Q(i,j,k))/pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) + 
     (2*(1 - rgrid[i]))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (2*Q.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr1(i,j,k))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) + 
     (2*Q.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr2(i,j,k))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  h1.diff(d001,i,j,k)*((4*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*Q(i,j,k)*Qr2(i,j,k))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) - 
     (4*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qr2(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  (32*pow(q,2)*R0*pow(rh,4)*a0(i,j,k)*a2(i,j,k)*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k))/
   (pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
  (32*pow(q,2)*R0*pow(rh,4)*pow(a0(i,j,k),2)*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],6)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*Qt2(i,j,k))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  (32*pow(q,2)*pow(rh,4)*a0(i,j,k)*a2(i,j,k)*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*R(i,j,k))/
   (pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
  (32*pow(q,2)*pow(rh,4)*pow(a0(i,j,k),2)*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],6)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*Qt2(i,j,k)*
     R(i,j,k))/(pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  (32*h1.diff(d011,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   (L1*L2*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
  (16*a2.diff(d010,i,j,k)*q*pow(rh,2)*h2(i,j,k)*
     (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   (L1*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
  (16*a1.diff(d001,i,j,k)*q*pow(1 + pow(R0,2),0.25)*pow(rh,2)*
     h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  (16*Qt2.diff(d010,i,j,k)*q*pow(rh,2)*a0(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   (L1*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  (16*Qt1.diff(d001,i,j,k)*q*pow(rh,2)*a0(i,j,k)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  (16*pow(q,2)*pow(rh,4)*pow(a2(i,j,k),2)*h1(i,j,k)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   (sqrt(1 + pow(R0,2))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  (16*pow(q,2)*pow(R0,2)*pow(rh,4)*pow(a2(i,j,k),2)*h1(i,j,k)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
   (sqrt(1 + pow(R0,2))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  (16*pow(q,2)*pow(rh,4)*pow(a0(i,j,k),2)*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],6)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*pow(Qt1(i,j,k),2))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  (16*pow(q,2)*pow(R0,2)*pow(rh,4)*pow(a0(i,j,k),2)*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],6)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*pow(Qt1(i,j,k),2))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
  (32*pow(q,2)*pow(rh,4)*a0(i,j,k)*a2(i,j,k)*h1(i,j,k)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k))/
   (pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
  (32*pow(q,2)*pow(R0,2)*pow(rh,4)*a0(i,j,k)*a2(i,j,k)*h1(i,j,k)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k))/
   (pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  (16*pow(q,2)*pow(rh,4)*pow(a0(i,j,k),2)*h1(i,j,k)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],6)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*pow(Qt2(i,j,k),2))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  (16*pow(q,2)*pow(R0,2)*pow(rh,4)*pow(a0(i,j,k),2)*h1(i,j,k)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],6)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*pow(Qt2(i,j,k),2))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  (32*pow(q,2)*R0*pow(rh,4)*pow(a2(i,j,k),2)*h1(i,j,k)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*R(i,j,k))/
   (sqrt(1 + pow(R0,2))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  (32*pow(q,2)*R0*pow(rh,4)*pow(a0(i,j,k),2)*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],6)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*pow(Qt1(i,j,k),2)*
     R(i,j,k))/((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
  (64*pow(q,2)*R0*pow(rh,4)*a0(i,j,k)*a2(i,j,k)*h1(i,j,k)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*R(i,j,k))/
   (pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  (32*pow(q,2)*R0*pow(rh,4)*pow(a0(i,j,k),2)*h1(i,j,k)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],6)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*pow(Qt2(i,j,k),2)*
     R(i,j,k))/(pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  (16*pow(q,2)*pow(rh,4)*pow(a2(i,j,k),2)*h1(i,j,k)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*pow(R(i,j,k),2))/
   (sqrt(1 + pow(R0,2))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  (16*pow(q,2)*pow(rh,4)*pow(a0(i,j,k),2)*h1(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),5)*pow(1 - rgrid[i],6)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*pow(Qt1(i,j,k),2)*
     pow(R(i,j,k),2))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) - 
  (32*pow(q,2)*pow(rh,4)*a0(i,j,k)*a2(i,j,k)*h1(i,j,k)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
     pow(R(i,j,k),2))/
   (pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  (16*pow(q,2)*pow(rh,4)*pow(a0(i,j,k),2)*h1(i,j,k)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),5)*pow(1 - rgrid[i],6)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*pow(Qt2(i,j,k),2)*
     pow(R(i,j,k),2))/
   (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
  h1.diff(d010,i,j,k)*((4*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*Q(i,j,k)*Qr1(i,j,k))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) - 
     (4*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*Qr1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
     (4*Q.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qr2(i,j,k))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) + 
     (16*R.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (8*Qrr.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (8*Qtt.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))) + 
  (q*h2(i,j,k)*pow(1 - rgrid[i],2)*Q(i,j,k)*
     (4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*
        (-4*L2*a2(i,j,k)*Qr2(i,j,k) + 
          4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*pow(1 - rgrid[i],2)*
           (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) + 
             L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
             Qtr(i,j,k)) - 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))\
))/(pow(1 + pow(R0,2),0.25)*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  (q*h2(i,j,k)*pow(1 - rgrid[i],2)*
     (-4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*
        (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
          4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*pow(1 - rgrid[i],2)*
           (-((1 - pow(1 - rgrid[i],2))*
                (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))) + 
             Qtr(i,j,k)) + 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))\
))/(pow(1 + pow(R0,2),0.25)*(1 - pow(1 - rgrid[i],2))*
     (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) + 
  Q.diff(d100,i,j,k)*(-((h1.diff(d100,i,j,k)*
          pow(1 - pow(1 - rgrid[i],2),2))/
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (2*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i]))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) + 
     (2*h1.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr1(i,j,k))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) + 
     (2*h1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
        (1 - rgrid[i])*Qr2(i,j,k))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) + 
     (q*h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*(1 - rgrid[i])*
        (4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*
           (-4*L2*a2(i,j,k)*Qr2(i,j,k) + 
             4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
              pow(1 - rgrid[i],2)*
              (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) + 
                L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                Qtr(i,j,k)) - 4*L1*a1(i,j,k)*Qr1(i,j,k)*
              sqrt(1 + pow(R0,2)))))/
      (2.*pow(1 + pow(R0,2),0.25)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  Q.diff(d010,i,j,k)*((4*h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) - 
     (4*h1.diff(d010,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*pow(Qr1(i,j,k),2))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) - 
     (4*h1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k)*Qr2(i,j,k))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) + 
     (q*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k)*(-4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*
           (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
             4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
              pow(1 - rgrid[i],2)*
              (-((1 - pow(1 - rgrid[i],2))*
                   (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))) \
+ Qtr(i,j,k)) + 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))))/
      (pow(1 + pow(R0,2),0.25)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  Q.diff(d001,i,j,k)*((4*h1(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) - 
     (4*h1.diff(d001,i,j,k)*pow(1 - pow(1 - rgrid[i],2),4)*
        pow(1 - rgrid[i],2)*pow(Qr2(i,j,k),2))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) + 
     (q*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr2(i,j,k)*(-4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*
           (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
             4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
              pow(1 - rgrid[i],2)*
              (-((1 - pow(1 - rgrid[i],2))*
                   (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))) \
+ Qtr(i,j,k)) + 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))))/
      (pow(1 + pow(R0,2),0.25)*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  (16*pow(q,2)*pow(rh,4)*pow(a1(i,j,k),2)*h1(i,j,k)*
     (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*sqrt(1 + pow(R0,2))*
     sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
       pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  (16*a0.diff(d001,i,j,k)*q*pow(rh,2)*h2(i,j,k)*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
       (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  a1.diff(d010,i,j,k)*((-4*L1*q*pow(1 + pow(R0,2),0.25)*h2(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
        pow(Qr1(i,j,k),2))/(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) + 
     (16*q*pow(1 + pow(R0,2),0.25)*pow(rh,4)*h2(i,j,k)*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr1(i,j,k),2)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (4.*pow(rh,2))))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  Qt1.diff(d010,i,j,k)*((16*L1*q*a0(i,j,k)*h2(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],4)*
        pow(Qr1(i,j,k),2))/(4 + 4*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
     (16*q*pow(rh,4)*a0(i,j,k)*h2(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
        ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr1(i,j,k),2)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (4.*pow(rh,2))))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  (8*pow(q,2)*pow(rh,4)*a1(i,j,k)*h1(i,j,k)*
     (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     (-4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
       4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           Qt2(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
          Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
   ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  h2.diff(d001,i,j,k)*((-2*q*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*
        (-4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*
           (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
             4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
              pow(1 - rgrid[i],2)*
              (-((1 - pow(1 - rgrid[i],2))*
                   (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))\
) + Qtr(i,j,k)) + 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))))/
      (pow(1 + pow(R0,2),0.25)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (8*q*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        ((L2*ar(i,j,k)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k))/
           pow(rh,2) - (pow(L2,2)*a2(i,j,k)*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2))/
           (pow(1 + pow(R0,2),0.25)*pow(rh,2)) + 
          (L1*L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qr2(i,j,k)*Qt1(i,j,k))/pow(rh,2) - 
          4*R0*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) + 
          (pow(L2,2)*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2)*Qt2(i,j,k))/pow(rh,2) - 
          (L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
             Qtr(i,j,k))/pow(rh,2) - 
          4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
           R(i,j,k) + 4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           (-(L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                 pow(-1 + rgrid[i],2)*
                 (pow(H,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   pow(rh,2)*
                    (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                 (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                 Qr2(i,j,k))/(4.*pow(rh,2)) + 
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) - 
          (4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           pow(1 + pow(R0,2),0.25) + 
          4*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  a0.diff(d010,i,j,k)*((16*q*h2(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],4)*
        Qr1(i,j,k)*((1 - pow(1 - rgrid[i],2))*
           (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k)) - 
          Qtr(i,j,k)))/(4 + 4*(1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
     (16*q*pow(rh,4)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],4)*((1 + 
             (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
           (1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
          (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (-(L1*(1 - pow(1 - rgrid[i],2))*pow(-1 + rgrid[i],2)*
                 (pow(H,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   pow(rh,2)*
                    (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                 (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                 ((1 - pow(1 - rgrid[i],2))*
                    (L1*Qr1(i,j,k)*Qt1(i,j,k) + 
                      L2*Qr2(i,j,k)*Qt2(i,j,k)) - Qtr(i,j,k)))/
              (4.*pow(rh,2)) + 
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  R.diff(d001,i,j,k)*((-16*h1.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*q*pow(rh,2)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        ((-4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
           pow(1 + pow(R0,2),0.25) + 
          4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             Qt1(i,j,k)*sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
      (L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  R.diff(d010,i,j,k)*((16*h1.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (16*h1.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) - 
     (4*q*pow(rh,2)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (-4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) + 
          (4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           pow(1 + pow(R0,2),0.25) + 
          4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (Qt1(i,j,k)*(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)) - 
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) - 
  (16*h1.diff(d020,i,j,k)*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
     pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
  (16*h1.diff(d002,i,j,k)*pow(rh,2)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (pow(L2,2)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  (16*a2.diff(d001,i,j,k)*q*pow(rh,2)*h2(i,j,k)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (L2*pow(1 + pow(R0,2),0.25)*pow(-1 + rgrid[i],2)*
     (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
       pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
  (16*Qt2.diff(d001,i,j,k)*q*pow(rh,2)*a0(i,j,k)*h2(i,j,k)*
     (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
     pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],4)*
     (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
     sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
   (L2*pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
        pow(rgrid[i],3) + pow(rh,2)*
        (pow(mu,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
             4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
     pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  Qrr.diff(d001,i,j,k)*((q*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*
        (-4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*
           (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
             4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
              pow(1 - rgrid[i],2)*
              (-((1 - pow(1 - rgrid[i],2))*
                   (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))\
) + Qtr(i,j,k)) + 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))))/
      (2.*pow(1 + pow(R0,2),0.25)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) + 
     (2*q*pow(rh,4)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        ((L2*ar(i,j,k)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k))/
           pow(rh,2) - (pow(L2,2)*a2(i,j,k)*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2))/
           (pow(1 + pow(R0,2),0.25)*pow(rh,2)) + 
          (L1*L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qr2(i,j,k)*Qt1(i,j,k))/pow(rh,2) + 
          4*R0*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) + 
          (pow(L2,2)*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2)*Qt2(i,j,k))/pow(rh,2) - 
          (L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
             Qtr(i,j,k))/pow(rh,2) + 
          4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
           R(i,j,k) - 4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           ((L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                Qr2(i,j,k))/(4.*pow(rh,2)) + 
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) + 
          (4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           pow(1 + pow(R0,2),0.25) - 
          4*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) - 
     (8*h1.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  H*((2*Q.diff(d100,i,j,k)*L1*L2*q*h2(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),3)*(1 - rgrid[i])*ygrid[k]*
        Qr1(i,j,k))/pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) + 
     (4*L1*L2*q*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*ygrid[k]*Q(i,j,k)*Qr1(i,j,k))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) - 
     (4*L1*L2*q*h2(i,j,k)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*ygrid[k]*Qr1(i,j,k))/
      (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
     (4*Q.diff(d010,i,j,k)*L1*L2*q*h2(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*ygrid[k]*
        pow(Qr1(i,j,k),2))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) - 
     (4*Q.diff(d001,i,j,k)*L1*L2*q*h2(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*ygrid[k]*
        Qr1(i,j,k)*Qr2(i,j,k))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2) + 
     (16*R.diff(d001,i,j,k)*q*pow(rh,2)*h2(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*ygrid[k]*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (8*Qrr.diff(d001,i,j,k)*q*pow(rh,2)*h2(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*ygrid[k]*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (32*h2.diff(d001,i,j,k)*q*pow(rh,2)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*ygrid[k]*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (pow(-1 + rgrid[i],2)*(pow(H,2)*pow(-2 + rgrid[i],3)*
           pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (16*R.diff(d010,i,j,k)*L2*q*pow(rh,2)*h2(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*ygrid[k]*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     Qtt.diff(d001,i,j,k)*((2*L1*L2*q*h2(i,j,k)*
           pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
           ygrid[k]*Qr1(i,j,k)*Qr2(i,j,k))/
         ((1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
        (8*q*pow(rh,4)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*ygrid[k]*
           (-(L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                 pow(-1 + rgrid[i],2)*
                 (pow(H,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   pow(rh,2)*
                    (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                 (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                 Qr2(i,j,k))/(4.*pow(rh,2)) + 
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))))/
         (pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))) + 
     Qrr.diff(d010,i,j,k)*((-2*L1*L2*q*h2(i,j,k)*
           pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
           ygrid[k]*pow(Qr1(i,j,k),2))/
         ((1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) - 
        (8*L2*q*pow(rh,4)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*ygrid[k]*
           ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) - 
             (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr1(i,j,k),2)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2))))/
         (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     h2.diff(d010,i,j,k)*((8*L1*L2*q*
           pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           ygrid[k]*pow(Qr1(i,j,k),2))/
         (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k)) - 
        (32*L2*q*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*
           pow(1 - rgrid[i],2)*ygrid[k]*
           ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr1(i,j,k),2)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2))))/
         (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     Qtt.diff(d010,i,j,k)*((2*L1*L2*q*h2(i,j,k)*
           pow(1 - pow(1 - rgrid[i],2),4)*pow(1 - rgrid[i],2)*
           ygrid[k]*pow(Qr1(i,j,k),2))/
         ((1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
        (8*L2*q*pow(rh,4)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*ygrid[k]*
           ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr1(i,j,k),2)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2))))/
         (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
     (8*q*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (-4*L2*q*a1(i,j,k)*h1(i,j,k)*ygrid[k]*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))*
           sqrt(1 + pow(R0,2)) + 
          4*L2*q*a2(i,j,k)*h1(i,j,k)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*ygrid[k]*
           (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          pow(1 + pow(R0,2),0.25)*
           (2*h2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
              sqrt(1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             4*L2*q*a0(i,j,k)*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*
              pow(1 - rgrid[i],2)*ygrid[k]*
              (Qt1(i,j,k)*(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) \
- (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*Qt2(i,j,k)*
                 (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))*
                 sqrt(1 + pow(R0,2) + 
                   2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                   pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))\
))/(pow(1 + pow(R0,2),0.25)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0,2) + 2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))) + 
     (16*B.diff(d010,i,j,k)*L2*q*pow(rh,2)*h2(i,j,k)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*ygrid[k]*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (L1*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  B.diff(d010,i,j,k)*((16*h1.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (16*q*pow(rh,2)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (L1*pow(1 + pow(R0,2),0.25)*a1(i,j,k) - 
          L1*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           Qt1(i,j,k))*sqrt(1 + 
          pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,2)*pow(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)),2)*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  B.diff(d001,i,j,k)*((-16*h1.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (16*q*pow(rh,2)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        ((L2*a2(i,j,k))/pow(1 + pow(R0,2),0.25) - 
          L2*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           Qt2(i,j,k))*sqrt(1 + 
          pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2))) + 
  Qtt.diff(d001,i,j,k)*((q*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*Qr2(i,j,k)*
        (4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*
           (-4*L2*a2(i,j,k)*Qr2(i,j,k) + 
             4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
              pow(1 - rgrid[i],2)*
              (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) + 
                L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                Qtr(i,j,k)) - 
             4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))))/
      (2.*pow(1 + pow(R0,2),0.25)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
     (2*q*pow(rh,4)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*
        ((L2*ar(i,j,k)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k))/
           pow(rh,2) - (pow(L2,2)*a2(i,j,k)*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2))/
           (pow(1 + pow(R0,2),0.25)*pow(rh,2)) + 
          (L1*L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qr2(i,j,k)*Qt1(i,j,k))/pow(rh,2) - 
          4*R0*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) + 
          (pow(L2,2)*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr2(i,j,k),2)*Qt2(i,j,k))/pow(rh,2) - 
          (L2*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
             pow(1 - rgrid[i],2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr2(i,j,k)*
             Qtr(i,j,k))/pow(rh,2) - 
          4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
           R(i,j,k) + 4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           (-(L1*L2*pow(1 - pow(1 - rgrid[i],2),2)*
                 pow(-1 + rgrid[i],2)*
                 (pow(H,2)*pow(-2 + rgrid[i],3)*
                    pow(rgrid[i],3) + 
                   pow(rh,2)*
                    (pow(mu,2)*pow(-2 + rgrid[i],3)*
                       pow(rgrid[i],3) + 
                      4*pow(rh,2)*
                       (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                 (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
                 Qr2(i,j,k))/(4.*pow(rh,2)) + 
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k))) - 
          (4*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           pow(1 + pow(R0,2),0.25) + 
          4*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
           sqrt(1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2))))/
      (L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
     (8*h1.diff(d001,i,j,k)*pow(rh,2)*
        (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L2,2)*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k)))) + 
  pow(H,2)*((-16*pow(L2,2)*pow(q,2)*pow(rh,4)*h1(i,j,k)*
        (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        pow(ygrid[k],2)*((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
           (1 + pow(R0,2) + 
             2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
             pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
          (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr1(i,j,k),2)*
             sqrt(1 + pow(R0,2) + 
               2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
               pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
           (2.*pow(rh,2))))/
      ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))) + 
     (8*pow(L2,2)*pow(q,2)*h1(i,j,k)*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*pow(ygrid[k],2)*
        (pow(L1,2)*pow(1 - pow(1 - rgrid[i],2),2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*pow(Qr1(i,j,k),2) \
+ (4*pow(rh,2)*(1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
           ((1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4)))))))/
      pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
  h2.diff(d010,i,j,k)*((-2*q*(1 - pow(1 - rgrid[i],2))*
        pow(1 - rgrid[i],2)*Qr1(i,j,k)*
        (-4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*
           (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
             4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
              pow(1 - rgrid[i],2)*
              (-((1 - pow(1 - rgrid[i],2))*
                   (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))\
) + Qtr(i,j,k)) + 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))))/
      (pow(1 + pow(R0,2),0.25)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))) - 
     (8*q*pow(rh,4)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
        (-4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
          4*pow(R0,2)*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
          4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) + 
          4*pow(R0,2)*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) - 
          8*R0*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           (1 - pow(1 - rgrid[i],2))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*R(i,j,k) + 
          8*R0*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
           R(i,j,k) - 4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           pow(1 - pow(1 - rgrid[i],2),2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*pow(R(i,j,k),2) \
+ 4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
           pow(R(i,j,k),2) + 
          (L1*ar(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
)/pow(rh,2) - (pow(L1,2)*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr1(i,j,k),2)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
)/pow(rh,2) - (L1*L2*a2(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qr2(i,j,k)*sqrt(1 + 
               pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
           (pow(1 + pow(R0,2),0.25)*pow(rh,2)) + 
          (4*R0*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
)/pow(1 + pow(R0,2),0.25) + 
          (pow(L1,2)*a0(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr1(i,j,k),2)*Qt1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
)/pow(rh,2) + (L1*L2*a0(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qr2(i,j,k)*Qt2(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
)/pow(rh,2) - 4*R0*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) - 
          (L1*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qtr(i,j,k)*sqrt(1 + 
               pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
           pow(rh,2) + (4*a2(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*R(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
)/pow(1 + pow(R0,2),0.25) - 
          4*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
           R(i,j,k)*sqrt(1 + pow(R0 + 
               (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  Qtt.diff(d010,i,j,k)*((8*h1.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
     (q*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*Qr1(i,j,k)*
        (4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*
           (-4*L2*a2(i,j,k)*Qr2(i,j,k) + 
             4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
              pow(1 - rgrid[i],2)*
              (L1*(1 - pow(1 - rgrid[i],2))*Qr1(i,j,k)*Qt1(i,j,k) + 
                L2*(1 - pow(1 - rgrid[i],2))*Qr2(i,j,k)*Qt2(i,j,k) - 
                Qtr(i,j,k)) - 
             4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))))/
      (2.*pow(1 + pow(R0,2),0.25)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) + 
     h1.diff(d010,i,j,k)*((2*pow(1 - pow(1 - rgrid[i],2),4)*
           pow(1 - rgrid[i],2)*pow(Qr1(i,j,k),2))/
         ((1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))) - 
        (8*pow(rh,2)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           ((1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
              (1 + pow(R0,2) + 
                2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)) + 
             (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
                pow(1 - pow(1 - rgrid[i],2),2)*
                pow(-1 + rgrid[i],2)*
                (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
                  pow(rh,2)*
                   (pow(mu,2)*pow(-2 + rgrid[i],3)*
                      pow(rgrid[i],3) + 
                     4*pow(rh,2)*
                      (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                        4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
                (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
                pow(Qr1(i,j,k),2)*
                sqrt(1 + pow(R0,2) + 
                  2*R0*(1 - pow(1 - rgrid[i],2))*R(i,j,k) + 
                  pow(1 - pow(1 - rgrid[i],2),2)*pow(R(i,j,k),2)))/
              (4.*pow(rh,2))))/
         (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(-1 + rgrid[i],2)*
           (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
             pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                 pow(rgrid[i],3) + 
                4*pow(rh,2)*
                 (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                   4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
           pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) - 
     (2*q*pow(rh,4)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(-4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
          4*pow(R0,2)*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
          4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) + 
          4*pow(R0,2)*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) - 
          8*R0*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           (1 - pow(1 - rgrid[i],2))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*R(i,j,k) + 
          8*R0*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
           R(i,j,k) - 4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           pow(1 - pow(1 - rgrid[i],2),2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*pow(R(i,j,k),2) \
+ 4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
           pow(R(i,j,k),2) + 
          (L1*ar(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
)/pow(rh,2) - (pow(L1,2)*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr1(i,j,k),2)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
)/pow(rh,2) - (L1*L2*a2(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qr2(i,j,k)*sqrt(1 + 
               pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
           (pow(1 + pow(R0,2),0.25)*pow(rh,2)) + 
          (4*R0*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
)/pow(1 + pow(R0,2),0.25) + 
          (pow(L1,2)*a0(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr1(i,j,k),2)*Qt1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
)/pow(rh,2) + (L1*L2*a0(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qr2(i,j,k)*Qt2(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
)/pow(rh,2) - 4*R0*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) - 
          (L1*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qtr(i,j,k)*sqrt(1 + 
               pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
           pow(rh,2) + (4*a2(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*R(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))\
)/pow(1 + pow(R0,2),0.25) - 
          4*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
           R(i,j,k)*sqrt(1 + pow(R0 + 
               (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qtt(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))) + 
  Qrr.diff(d010,i,j,k)*((8*h1.diff(d001,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        (R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k)))/
      (L1*L2*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) + 
     (q*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        Qr1(i,j,k)*(-4*pow(1 + pow(R0,2),0.25)*ar(i,j,k) + 
          pow(1 - pow(1 - rgrid[i],2),2)*
           (4*L2*a2(i,j,k)*Qr2(i,j,k) + 
             4*pow(1 + pow(R0,2),0.25)*a0(i,j,k)*
              pow(1 - rgrid[i],2)*
              (-((1 - pow(1 - rgrid[i],2))*
                   (L1*Qr1(i,j,k)*Qt1(i,j,k) + L2*Qr2(i,j,k)*Qt2(i,j,k))) \
+ Qtr(i,j,k)) + 4*L1*a1(i,j,k)*Qr1(i,j,k)*sqrt(1 + pow(R0,2)))))/
      (2.*pow(1 + pow(R0,2),0.25)*
        (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))) - 
     (8*h1.diff(d010,i,j,k)*pow(rh,2)*
        pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
      (pow(L1,2)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
        pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)) - 
     (2*q*pow(rh,4)*h2(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
        pow(1 - rgrid[i],2)*(-4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) - 
          4*pow(R0,2)*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k)) + 
          4*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) + 
          4*pow(R0,2)*a0(i,j,k)*(1 - pow(1 - rgrid[i],2))*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k) - 
          8*R0*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           (1 - pow(1 - rgrid[i],2))*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*R(i,j,k) + 
          8*R0*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),2)*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*R(i,j,k) \
- 4*pow(1 + pow(R0,2),0.25)*a1(i,j,k)*
           pow(1 - pow(1 - rgrid[i],2),2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*pow(R(i,j,k),2) + 
          4*a0(i,j,k)*pow(1 - pow(1 - rgrid[i],2),3)*
           pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt1(i,j,k)*
           pow(R(i,j,k),2) - 
          (L1*ar(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
           pow(rh,2) + (pow(L1,2)*pow(1 + pow(R0,2),0.25)*
             a1(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr1(i,j,k),2)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
           pow(rh,2) + (L1*L2*a2(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qr2(i,j,k)*sqrt(1 + 
               pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
           (pow(1 + pow(R0,2),0.25)*pow(rh,2)) + 
          (4*R0*a2(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
           pow(1 + pow(R0,2),0.25) - 
          (pow(L1,2)*a0(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*
             pow(Qr1(i,j,k),2)*Qt1(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
           pow(rh,2) - (L1*L2*a0(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),3)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qr2(i,j,k)*Qt2(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
           pow(rh,2) - 4*R0*a0(i,j,k)*
           (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           (1 - pow(1 - rgrid[i],2))*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)) + 
          (L1*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
             pow(-1 + rgrid[i],2)*
             (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
               pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
                   pow(rgrid[i],3) + 
                  4*pow(rh,2)*
                   (1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                     4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
             (1 + (1 - pow(1 - rgrid[i],2))*Q(i,j,k))*Qr1(i,j,k)*
             Qtr(i,j,k)*sqrt(1 + 
               pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
           pow(rh,2) + (4*a2(i,j,k)*
             (1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
             (1 - pow(1 - rgrid[i],2))*
             (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*R(i,j,k)*
             sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2)))/
           pow(1 + pow(R0,2),0.25) - 
          4*a0(i,j,k)*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*
           pow(1 - pow(1 - rgrid[i],2),2)*pow(1 - rgrid[i],2)*
           (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*Qt2(i,j,k)*R(i,j,k)*
           sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))/
      (L1*(1 + B(i,j,k)*(1 - pow(1 - rgrid[i],2)))*pow(-1 + rgrid[i],2)*
        (pow(H,2)*pow(-2 + rgrid[i],3)*pow(rgrid[i],3) + 
          pow(rh,2)*(pow(mu,2)*pow(-2 + rgrid[i],3)*
              pow(rgrid[i],3) + 
             4*pow(rh,2)*(1 + 2*rgrid[i] + 3*pow(rgrid[i],2) - 
                4*pow(rgrid[i],3) + pow(rgrid[i],4))))*
        pow(rh + rh*(1 - pow(1 - rgrid[i],2))*Q(i,j,k),2)*
        (1 + (1 - pow(1 - rgrid[i],2))*Qrr(i,j,k))*
        sqrt(1 + pow(R0 + (1 - pow(1 - rgrid[i],2))*R(i,j,k),2))))
;

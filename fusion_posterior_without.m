T=0.5;N=1000;dt=T/N;D=5;L=500;
x0=0;y0=0;theta0=0;xend=1;yend=0;
x1=-1;y1=1;xend1=-0.35;yend1=1;
x2=-1;y2=-1;xend2=-0.35;yend2=-1;
x3=1;y3=0;xend3=1.65;yend3=0;
xi=0;yi=0;xendi=0.65;yendi=0;
v=1.3;l=0.200;r=0.033;

w1 = v/r;
w2 = v/r;


for i=1:L
   randn('state',i+1)
   dW1 = sqrt(dt) * randn(1,N);
   randn('state',i+10002)
   dW2 = sqrt(dt) * randn(1,N); %Wiener process
   xtemp=x1;
   ytemp=y1;
   thetatemp=theta0; %Initialization
   for j=1:N
      xtemp = xtemp+((r*cos(thetatemp)*(w1+w2)*dt)/2)+((sqrt(D)*r*cos(thetatemp)*(dW1(j)+dW2(j)))/2);
      ytemp = ytemp+((r*sin(thetatemp)*(w1+w2)*dt)/2)+((sqrt(D)*r*sin(thetatemp)*(dW1(j)+dW2(j)))/2);
      thetatemp = thetatemp+((r*(w1-w2)*dt)/l)+((sqrt(D)*r*(dW1(j)-dW2(j)))/l);    %kinematic equation with SDE
      x(i,j)=xtemp;
      y(i,j)=ytemp;
      t(i,j)=thetatemp; %store data 
  end
  xf1(i)=x(i,N);
  yf1(i)=y(i,N);
  tf1(i)=t(i,N);%assemble of final pose of each path
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% team member1
x1t=xf1(250);
y1t=yf1(250);
t1t=tf1(250);

for i=1:L
   randn('state',i+10001)
   dW1 = sqrt(dt) * randn(1,N);
   randn('state',i+20002)
   dW2 = sqrt(dt) * randn(1,N); %Wiener process
   xtemp=x2;
   ytemp=y2;
   thetatemp=theta0; %Initialization
   for j=1:N
      xtemp = xtemp+((r*cos(thetatemp)*(w1+w2)*dt)/2)+((sqrt(D)*r*cos(thetatemp)*(dW1(j)+dW2(j)))/2);
      ytemp = ytemp+((r*sin(thetatemp)*(w1+w2)*dt)/2)+((sqrt(D)*r*sin(thetatemp)*(dW1(j)+dW2(j)))/2);
      thetatemp = thetatemp+((r*(w1-w2)*dt)/l)+((sqrt(D)*r*(dW1(j)-dW2(j)))/l);    %kinematic equation with SDE
      x_2(i,j)=xtemp;
      y_2(i,j)=ytemp;
      t_2(i,j)=thetatemp; %store data 
  end
  xf2(i)=x_2(i,N);
  yf2(i)=y_2(i,N);
  tf2(i)=t_2(i,N);%assemble of final pose of each path
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%team member2
x2t=xf2(250);
y2t=yf2(250);
t2t=tf2(250);

for i=1:L
   randn('state',i+20001)
   dW1 = sqrt(dt) * randn(1,N);
   randn('state',i+30002)
   dW2 = sqrt(dt) * randn(1,N); %Wiener process
   xtemp=x3;
   ytemp=y3;
   thetatemp=theta0; %Initialization
   for j=1:N
      xtemp = xtemp+((r*cos(thetatemp)*(w1+w2)*dt)/2)+((sqrt(D)*r*cos(thetatemp)*(dW1(j)+dW2(j)))/2);
      ytemp = ytemp+((r*sin(thetatemp)*(w1+w2)*dt)/2)+((sqrt(D)*r*sin(thetatemp)*(dW1(j)+dW2(j)))/2);
      thetatemp = thetatemp+((r*(w1-w2)*dt)/l)+((sqrt(D)*r*(dW1(j)-dW2(j)))/l);    %kinematic equation with SDE
      x_3(i,j)=xtemp;
      y_3(i,j)=ytemp;
      t_3(i,j)=thetatemp; %store data 
  end
  xf3(i)=x_3(i,N);
  yf3(i)=y_3(i,N);
  tf3(i)=t_3(i,N);%assemble of final pose of each path
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%team member3
x3t=xf3(250);
y3t=yf3(250);
t3t=tf3(250);

for i=1:L
   randn('state',i+30001)
   dW1 = sqrt(dt) * randn(1,N);
   randn('state',i+40002)
   dW2 = sqrt(dt) * randn(1,N); %Wiener process
   xtemp=xi;
   ytemp=yi;
   thetatemp=theta0; %Initialization
   for j=1:N
      xtemp = xtemp+((r*cos(thetatemp)*(w1+w2)*dt)/2)+((sqrt(D)*r*cos(thetatemp)*(dW1(j)+dW2(j)))/2);
      ytemp = ytemp+((r*sin(thetatemp)*(w1+w2)*dt)/2)+((sqrt(D)*r*sin(thetatemp)*(dW1(j)+dW2(j)))/2);
      thetatemp = thetatemp+((r*(w1-w2)*dt)/l)+((sqrt(D)*r*(dW1(j)-dW2(j)))/l);    %kinematic equation with SDE
      x_i(i,j)=xtemp;
      y_i(i,j)=ytemp;
      t_i(i,j)=thetatemp; %store data 
  end
  xfi(i)=x_i(i,N);
  yfi(i)=y_i(i,N);
  tfi(i)=t_i(i,N);%assemble of final pose of each path
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%team core
xit=xfi(250);
yit=yfi(250);
tit=tfi(250);

me1=[1 0 v*T;0 1 0;0 0 1];
co1 = [D*r^2*T/2 0 0;0 2*D*w1^2*r^4*T^3/(3*l^2) D*w1*r^3*T^2/l^2;0 D*w1*r^3*T^2/l^2 2*D*r^2*T/l^2];

me2=[1 0 v*T;0 1 0;0 0 1];
co2=[D*r^2*T/2 0 0;0 2*D*w1^2*r^4*T^3/(3*l^2) D*w1*r^3*T^2/l^2;0 D*w1*r^3*T^2/l^2 2*D*r^2*T/l^2];

me3=[1 0 v*T;0 1 0;0 0 1];
co3=[D*r^2*T/2 0 0;0 2*D*w1^2*r^4*T^3/(3*l^2) D*w1*r^3*T^2/l^2;0 D*w1*r^3*T^2/l^2 2*D*r^2*T/l^2];

mei=[1 0 xi+v*T;0 1 yi;0 0 1];
coi=[D*r^2*T/2 0 0;0 2*D*w1^2*r^4*T^3/(3*l^2) D*w1*r^3*T^2/l^2;0 D*w1*r^3*T^2/l^2 2*D*r^2*T/l^2];

a1=[1 0 x1;0 1 y1;0 0 1];
a2=[1 0 x2;0 1 y2;0 0 1];
a3=[1 0 x3;0 1 y3;0 0 1];
ai=[1 0 xi;0 1 yi;0 0 1];

g1=[cos(t1t) -sin(t1t) x1t;sin(t1t) cos(t1t) y1t;0 0 1];
g2=[cos(t2t) -sin(t2t) x2t;sin(t2t) cos(t2t) y2t;0 0 1];
g3=[cos(t3t) -sin(t3t) x3t;sin(t3t) cos(t3t) y3t;0 0 1];
gi=[cos(tit) -sin(tit) xit;sin(tit) cos(tit) yit;0 0 1];

mi1=inv(gi)*g1;
mi2=inv(gi)*g2;
mi3=inv(gi)*g3;

E=[1 0 0;0 1 0;0 0 1];
Xi=logm(E);
adXi=[0 Xi(1,2) Xi(2,3);Xi(2,1) 0 -Xi(1,3);0 0 0];
Ti=E+adXi/2;
Si=Ti'*inv(coi)*Ti;%%%%%%%%%%%%%%%%%%%Si

q1=mi1*inv(me1)*inv(a1)*ai*mei;
X1=logm(q1);
xo1=[X1(1,3);X1(2,3);X1(2,1)];
adX1=[0 X1(1,2) X1(2,3);X1(2,1) 0 -X1(1,3);0 0 0];
Admi1=[mi1(1,1) mi1(1,2) mi1(2,3);mi1(2,1) mi1(2,2) -mi1(1,3);0 0 1];
Ad1=inv(Admi1);
T1=E+adX1/2;
S1=T1'*Ad1'*inv(co1)*Ad1*T1;%%%%%%%%%%%%%%%%%%%%%%%%S1

q2=mi2*inv(me2)*inv(a2)*ai*mei;
X2=logm(q2);
xo2=[X2(1,3);X2(2,3);X2(2,1)];
adX2=[0 X2(1,2) X2(2,3);X2(2,1) 0 -X2(1,3);0 0 0];
Admi2=[mi2(1,1) mi2(1,2) mi2(2,3);mi2(2,1) mi2(2,2) -mi2(1,3);0 0 1];
Ad2=inv(Admi2);
T2=E+adX2/2;
S2=T2'*Ad2'*inv(co2)*Ad2*T2;%%%%%%%%%%%%%%%%%%%%%%%%S2

q3=mi3*inv(me3)*inv(a3)*ai*mei;
X3=logm(q3);
xo3=[X3(1,3);X3(2,3);X3(2,1)];
adX3=[0 X3(1,2) X3(2,3);X3(2,1) 0 -X3(1,3);0 0 0];
Admi3=[mi3(1,1) mi3(1,2) mi3(2,3);mi3(2,1) mi3(2,2) -mi3(1,3);0 0 1];
Ad3=inv(Admi3);
T3=E+adX3/2;
S3=T3'*Ad3'*inv(co3)*Ad3*T3;%%%%%%%%%%%%%%%%%%%%%%%%S3

Sba=Si+S1+S2+S3;
xba=inv(Sba)*(S1*xo1+S2*xo2+S3*xo3);
%xba=[0.05;-0.02;0];
Xba=[0 -xba(3,1) xba(1,1);xba(3,1) 0 xba(2,1);0 0 0];

mui=mei*expm(-Xba);%fusion mean

Xdot=logm(inv(mui)*mei);
adXd=[0 Xdot(1,2) Xdot(2,3);Xdot(2,1) 0 -Xdot(1,3);0 0 0];
Tdot=E+adXd;
covf=Tdot*inv(Sba)*Tdot';%fusion covariance

plot([x1,xend1],[y1,yend1],'r--'),hold on 
grid on
scatter(x1,y1,'*','r'),hold on
axis([-1 2 -2 1.5]);
scatter(xf1,yf1,'.','g'),hold on
scatter(xend1,yend1,'*','r'),hold on
scatter(x1t,y1t,'x','k','LineWidth',1.5),hold on

plot([x2,xend2],[y2,yend2],'r--'),hold on 
scatter(x2,y2,'*','r'),hold on
scatter(xf2,yf2,'.'),hold on
scatter(xend2,yend2,'*','r'),hold on
scatter(x2t,y2t,'x','k','LineWidth',1.5),hold on

plot([x3,xend3],[y3,yend3],'r--'),hold on 
scatter(x3,y3,'*','r'),hold on
scatter(xf3,yf3,'.','c'),hold on
scatter(xend3,yend3,'*','r'),hold on
scatter(x3t,y3t,'x','k','LineWidth',1.5),hold on

plot([xi,xendi],[yi,yendi],'r--'),hold on 
scatter(xi,yi,'*','r'),hold on
scatter(xfi,yfi,'.','m'),hold on
scatter(xendi,yendi,'*','r'),hold on
scatter(xit,yit,'x','k','LineWidth',1.5),hold on
scatter(mui(1,3),mui(2,3),'+','g')
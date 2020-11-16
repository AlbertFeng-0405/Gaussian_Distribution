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

UE = [1 0 v*T;0 1 0;0 0 1];
cov2 = [D*r^2*T/2 0 0;0 2*D*w1^2*r^4*T^3/(3*l^2) D*w1*r^3*T^2/l^2;0 D*w1*r^3*T^2/l^2 2*D*r^2*T/l^2];

CsE = ((2*pi)^1.5)*sqrt(abs(det(cov2)));
RN=200;dT=0.01;LN=500;
 xe=0:dT:2;
 ye=-1:dT:1;
 [XE,YE]=meshgrid(xe,ye);
 fxy = zeros(RN+1,RN+1);
 stp = 2*pi/LN;
 oe = -pi:stp:pi;
 oe(find(oe==0))=[];
 for i=1:RN+1
     for j=1:RN+1
         s=0;
         fg = zeros(1,LN-1);
         fg1 = zeros(1,LN-1);
         for k=1:LN-1
         g = [cos(oe(k)) -sin(oe(k)) xe(j);sin(oe(k)) cos(oe(k)) ye(i);0 0 1];
         yT = logm(inv(UE)*g);
         EE = [yT(1,3) yT(2,3) yT(2,1)]*inv(cov2)*[yT(1,3);yT(2,3);yT(2,1)];
         fg(k) = exp((-1*EE)/2)/CsE;
         g1 = [cos(oe(k+1)) -sin(oe(k+1)) xe(j);sin(oe(k+1)) cos(oe(k+1)) ye(i);0 0 1];
         yT1 = logm(inv(UE)*g1);
         EE1 = [yT1(1,3) yT1(2,3) yT1(2,1)]*inv(cov2)*[yT1(1,3);yT1(2,3);yT1(2,1)];
         fg1(k) = exp((-1*EE1)/2)/CsE;
         s = s+(fg1(k)+fg(k))*stp/2;
         end
         fxy(i,j) = s;
     end
 end

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

zmax = fix(max(max(fxy)));
zmin = 0.05;
L = [zmin,(zmin+zmax)/3,zmax];
contour(XE,YE,fxy,L,'k','LineWidth',1.5);

figure
mesh(XE,YE,fxy),title('Epdf')


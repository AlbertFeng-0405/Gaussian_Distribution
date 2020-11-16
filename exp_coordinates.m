T=1;N=1000;dt=T/N;D=2;L=10000;
x0=0;y0=0;theta0=0;xend=1;yend=0;
v=1;l=0.200;r=0.033;

w1 = v/r;
w2 = v/r;


for i=1:L
   randn('state',i+1)
   dW1 = sqrt(dt) * randn(1,N);
   randn('state',i+10002)
   dW2 = sqrt(dt) * randn(1,N); %Wiener process
   xtemp=x0;
   ytemp=y0;
   thetatemp=theta0; %Initialization
   for j=1:N
      xtemp = xtemp+((r*cos(thetatemp)*(w1+w2)*dt)/2)+((sqrt(D)*r*cos(thetatemp)*(dW1(j)+dW2(j)))/2);
      ytemp = ytemp+((r*sin(thetatemp)*(w1+w2)*dt)/2)+((sqrt(D)*r*sin(thetatemp)*(dW1(j)+dW2(j)))/2);
      thetatemp = thetatemp+((r*(w1-w2)*dt)/l)+((sqrt(D)*r*(dW1(j)-dW2(j)))/l);    %kinematic equation with SDE
      x(i,j)=xtemp;
      y(i,j)=ytemp;
      t(i,j)=thetatemp; %store data 
  end
  xf(i)=x(i,N);
  yf(i)=y(i,N);
  tf(i)=t(i,N);%assemble of final pose of each path
end

for m = 1:L
    v1(m)=(tf(m)*sin(tf(m))*xf(m)+tf(m)*(1-cos(tf(m)))*yf(m))/(2-2*cos(tf(m)));
    v2(m)=(tf(m)*sin(tf(m))*yf(m)+tf(m)*(cos(tf(m))-1)*xf(m))/(2-2*cos(tf(m)));
end

xm=sum(v1)/L;ym=sum(v2)/L;

multi=[0 0;0 0];
for o=1:L
    multi=multi+([v1(o)-xm;v2(o)-ym]*[v1(o)-xm v2(o)-ym]);
end
cov1 = multi/L;

uc=[xm;ym];
COVC=[cov1(1,1) cov1(1,2);cov1(2,1) cov1(2,2)];
xc=0:0.01:2;
yc=-1:0.01:1;

[XC,YC]=meshgrid(xc,yc);
s2xc=COVC(1,1);
s2yc=COVC(2,2);
sxc=sqrt(s2xc);
syc=sqrt(s2yc);
cc=COVC(1,2);
rc=cc/(sxc*syc);
ac=1/(2*pi*sxc*syc*sqrt(1-rc^2));
b1c=-1/(2*(1-rc^2));
b2c=((XC-uc(1))./sxc).^2;
b3c=((YC-uc(2))./syc).^2;
b4c=2*rc.*(XC-uc(1)).*(YC-uc(2))./(sxc*syc)
ZC=ac*exp(b1c*(b2c+b3c-b4c));


plot([x0,xend],[y0,yend],'r--'),hold on 
grid on
scatter(x0,y0,'*','r'),hold on
axis([-0.5 1.5 -1 1]);
scatter(v1,v2,'.'),hold on
scatter(xend,yend,'*','r'),hold on
mmax = fix(max(max(ZC)));
mmin = 0.2;
M = [mmin,(mmin+mmax)/3,mmax];
contour(XC,YC,ZC,M,'m','LineWidth',1.5),hold on
xlabel('v1','FontSize',16);
ylabel('v2','FontSize',16,'Rotation',0,'HorizontalAlignment','right');
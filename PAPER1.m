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

xm=sum(xf)/L;ym=sum(yf)/L;tm=sum(tf)/L;
mean1 = [xm;ym;tm];

multi=[0 0 0;0 0 0;0 0 0];
for o=1:L
    multi=multi+([xf(o)-xm;yf(o)-ym;tf(o)-tm]*[xf(o)-xm yf(o)-ym tf(o)-tm]);
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
ZC=ac*exp(b1c*(b2c+b3c-b4c));%Gaussian Distribution in Cartesian coordinates

ue=cell(1,100);
g=cell(1,L);
er=cell(1,L);
g{1,1}=[0 0 0;0 0 0;0 0 0];
for r=1:100
    er{1,r}=[0 0 0;0 0 0;0 0 0];
end
for p=1:L
    g{1,p}=[cos(tf(p)) -sin(tf(p)) xf(p);sin(tf(p)) cos(tf(p)) yf(p);0 0 1];
end
ue{1,1}=g{1,1};
for q=1:99
    tempmat = [0 0 0;0 0 0;0 0 0];
    for s=1:L
        tempmat = tempmat + logm(g{1,s}*inv(ue{1,q}));
        er{1,q} = tempmat;
    end
    ue{1,q+1} = ue{1,q}*expm(er{1,q}/L);
end

UE=ue{1,100};%mean of exponential coordinates

multie = [0 0 0;0 0 0;0 0 0];
ydelta = cell(1,L);
yT = cell(1,L);
for w=1:L
    tempsp = [0 0 0;0 0 0;0 0 0];
    tempsp = logm(inv(UE)*g{1,w});
    ydelta{1,w} = [tempsp(1,3);tempsp(2,3);tempsp(2,1)];
    yT{1,w} = [tempsp(1,3) tempsp(2,3) tempsp(2,1)];
    multie = multie + ydelta{1,w}*yT{1,w};
end
cov2 = multie/L;%covariance of exponential coordinates
%
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

%
meanprop = [1 0 v*T;0 1 0;0 0 1];

covprop = [D*r^r*T/2 0 0;0 2*D*w1^2*r^4*T^3/3*l^2 D*w1*r^3*T^2/l^2;0 D*w1*r^3*T^2/l^2 2*D*r^2*T/l^2];


plot([x0,xend],[y0,yend],'r--'),hold on 
scatter(x0,y0,'*','r'),hold on
axis([-0.5 1.5 -1 1]);
scatter(xf,yf,'.'),hold on
scatter(xend,yend,'*','r'),hold on
mmax = fix(max(max(ZC)));
mmin = 0.05;
M = [mmin,(mmin+mmax)/3,mmax];
contour(XC,YC,ZC,M,'m','LineWidth',1.5),hold on
zmax = fix(max(max(fxy)));
zmin = 0.05;
L = [zmin,(zmin+zmax)/3,zmax];
contour(XE,YE,fxy,L,'k','LineWidth',1.5);
grid on
xlabel('x','FontSize',16);
ylabel('y','FontSize',16,'Rotation',0,'HorizontalAlignment','right');
figure
mesh(XC,YC,ZC),title('Cpdf')
figure
mesh(XE,YE,fxy),title('Epdf')


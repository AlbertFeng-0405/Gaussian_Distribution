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

xm=sum(v1)/L;ym=sum(v2)/L;tm=sum(tf)/L;

multie=[0 0 0;0 0 0;0 0 0];
for o=1:L
    multie=multie+([v1(o)-xm;v2(o)-ym;tf(o)-tm]*[v1(o)-xm v2(o)-ym tf(o)-tm]);
end
cov2 = multie/L;

d1=zeros(1,L);
for q=1:L
    d1(q)=[v1(q)-xm v2(q)-ym tf(q)-tm]*inv(cov2)*[v1(q)-xm;v2(q)-ym;tf(q)-tm];
end
d2=sort(d1);

x=0:30;
y=x;
pt=((1:L)-0.5)/L;  
x2=chi2inv(pt,3); 
scatter(x2,d2','*'),hold on
plot(x,y);
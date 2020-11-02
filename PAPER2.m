T=1;N=1000;dt=T/N;D=1.5;L=10000;
x0=0;y0=0;theta0=0;
dalpha=pi/2;l=0.200;r=0.033;a=8;
xend=8;yend=8;
b=0:-0.001:-pi/2;

w1 = (dalpha*(a+l/2))/r;
w2 = (dalpha*(a-l/2))/r;

xr = a*cos(b);
yr = a+a*sin(b);

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

scatter(x0,y0,'*','r'),hold on
axis([-2 13 -2 13]);
scatter(xf,yf,'.'),hold on
plot(xr,yr,'r--'),hold on
scatter(xend,yend,'*','r'),hold off
grid on
xlabel('x','FontSize',16);
ylabel('y','FontSize',16,'Rotation',0,'HorizontalAlignment','right');
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

ue=cell(1,100);
g=cell(1,L);
er=cell(1,L);
g{1,1}=[0 0 0;0 0 0;0 0 0];
for r=1:100
    er{1,r}=[0 0 0;0 0 0;0 0 0];
end
for p=1:L
    %g{1,p}=[cos(tf(p)) -sin(tf(p)) ((yf(p)*(-1+cos(tf(p)))+xf(p)*sin(tf(p)))/tf(p));sin(tf(p)) cos(tf(p)) ((xf(p)*(1-cos(tf(p)))+yf(p)*sin(tf(p)))/tf(p));0 0 1];
    g{1,p}=[cos(tf(p)) -sin(tf(p)) xf(p);sin(tf(p)) cos(tf(p)) yf(p);0 0 1];
end
%ue{1,1}=g{1,1};
ue{1,1}=[1 0 1;0 1 0;0 0 1];
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
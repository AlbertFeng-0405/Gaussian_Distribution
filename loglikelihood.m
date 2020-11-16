T=1;N=1000;dt=T/N;L=10000;
x0=0;y0=0;theta0=0;xend=1;yend=0;
v=1;l=0.200;r=0.033;

D=1:1:7;

w1 = v/r;
w2 = v/r;

for m=1:7
    for i=1:L
     randn('state',i+1)
     dW1 = sqrt(dt) * randn(1,N);
     randn('state',i+10002)
     dW2 = sqrt(dt) * randn(1,N); %Wiener process
     xtemp=x0;
     ytemp=y0;
     thetatemp=theta0; %Initialization
     for j=1:N
      xtemp = xtemp+((r*cos(thetatemp)*(w1+w2)*dt)/2)+((sqrt(D(m))*r*cos(thetatemp)*(dW1(j)+dW2(j)))/2);
      ytemp = ytemp+((r*sin(thetatemp)*(w1+w2)*dt)/2)+((sqrt(D(m))*r*sin(thetatemp)*(dW1(j)+dW2(j)))/2);
      thetatemp = thetatemp+((r*(w1-w2)*dt)/l)+((sqrt(D(m))*r*(dW1(j)-dW2(j)))/l);    %kinematic equation with SDE
     end
    xf(m,i)=xtemp;
    yf(m,i)=ytemp;
    tf(m,i)=thetatemp;%assemble of final pose of each path
    end
end
%%%%%%%%
xm=sum(xf,2)/L;ym=sum(yf,2)/L;
mean1=cell(1,7);
for q=1:7
   mean1{1,q} = [xm(q);ym(q)];
end
%%%%%%%%
cov1=cell(1,7);
multi=cell(1,7);
for n=1:7
    multi{1,n}=[0 0;0 0];
    for o=1:L
        multi{1,n}=multi{1,n}+([xf(n,o)-xm(n);yf(n,o)-ym(n)]*[xf(n,o)-xm(n) yf(n,o)-ym(n)]);
    end
    cov1{1,n} = multi{1,n}/L;
end
%%%%%
XT=zeros(1,7);
for l=1:7
    for k=1:L
       XT(l)=XT(l)-([xf(l,k)-xm(l) yf(l,k)-xm(l)]*inv(cov1{1,l})*[xf(l,k)-xm(l);yf(l,k)-xm(l)])/2;
    end
end
LL_CART=zeros(1,7);
for e=1:7
    LL_CART(e) = -L*log(2*pi)-L*log(det(cov1{1,e}))/2+XT(e);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ue=cell(7,100);
g=cell(7,L);
er=cell(7,L);
for o=1:7
  for p=1:L
    g{o,p}=[cos(tf(o,p)) -sin(tf(o,p)) xf(o,p);sin(tf(o,p)) cos(tf(o,p)) yf(o,p);0 0 1];
  end
end
for lm=1:7
    ue{lm,1}=g{lm,1};
end
tempmat=cell(1,7);
for e=1:7
    for q=1:99
      tempmat{1,e} = [0 0 0;0 0 0;0 0 0];
      for s=1:L
        tempmat{1,e} = tempmat{1,e} + logm(g{e,s}*inv(ue{e,q}));
        er{e,q} = tempmat{1,e};
      end
    ue{e,q+1} = ue{e,q}*expm(er{e,q}/L);
    end
end
UE=cell(1,7);
for ne=1:7
    UE{1,ne}=ue{ne,100};
end%mean of exponential coordinates

multie=cell(1,7);
ydelta = cell(7,L);
yT = cell(7,L);
cov2=cell(1,7);
for r=1:7
multie{1,r} = [0 0 0;0 0 0;0 0 0];

 for w=1:L
    tempsp = [0 0 0;0 0 0;0 0 0];
    tempsp = logm(inv(UE{1,r})*g{r,w});
    ydelta{r,w} = [tempsp(1,3);tempsp(2,3);tempsp(2,1)];
    yT{r,w} = [tempsp(1,3) tempsp(2,3) tempsp(2,1)];
    multie{1,r} = multie{1,r} + ydelta{r,w}*yT{r,w};
 end
cov2{1,r} = multie{1,r}/L;
end

YT=zeros(1,7);
ytri=cell(7,L);
for h=1:7
    for s=1:L
        ytri{h,s}=logm(inv(UE{1,h})*[cos(tf(h,s)) -sin(tf(h,s)) xf(h,s);sin(tf(h,s)) cos(tf(h,s)) yf(h,s);0 0 1]);
        YT(h)=YT(h)-([ytri{h,s}(1,3) ytri{h,s}(2,3) ytri{h,s}(2,1)]*inv(cov2{1,h})*[ytri{h,s}(1,3);ytri{h,s}(2,3);ytri{h,s}(2,1)])/2;
    end
end

LL_EXP=zeros(1,7);
for ow=1:7
    LL_EXP(ow) = -3*L*log(2*pi)/2-L*log(det(cov2{1,ow}))/2+YT(ow);
end

LLR=zeros(1,7);
for esc=1:7
    LLR(esc)=exp(LL_EXP(esc)-LL_CART(esc));
end

plot(D,LL_CART),hold on
plot(D,LL_EXP),hold on
plot(D,LLR),hold off

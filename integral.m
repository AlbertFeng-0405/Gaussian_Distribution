UE = [ 1.0000    0.0060    0.9908;-0.0060    1.0000   -0.0080;0         0    1.0000];
cov2 = [0.0007    0.0002    0.0002;0.0002    0.0574    0.0538;0.0002    0.0538    0.0547];

 CsE = 2*pi*sqrt(det(cov2));


% x=0.5:0.01:1.5;
% y=-0.5:0.01:0.5;

 RN=100;dT=1/RN;LN=10;
 xe=0.5:dT:1.5;
 ye=-0.5:dT:0.5;
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
         g = [cos(oe(k)) -sin(oe(k)) (ye(i)*(cos(oe(k))-1)+xe(j)*sin(oe(k)))/oe(k);sin(oe(k)) cos(oe(k)) (xe(j)*(1-cos(oe(k)))+ye(i)*sin(oe(k)))/oe(k);0 0 1];
         yT = logm(inv(UE)*g);
         EE = [yT(1,3) yT(2,3) yT(2,1)]*inv(cov2)*[yT(1,3);yT(2,3);yT(2,1)];
         fg(k) = exp((-1*EE)/2)/CsE;
         g1 = [cos(oe(k+1)) -sin(oe(k+1)) (ye(i)*(cos(oe(k+1))-1)+xe(j)*sin(oe(k+1)))/oe(k+1);sin(oe(k+1)) cos(oe(k+1)) (xe(j)*(1-cos(oe(k+1)))+ye(i)*sin(oe(k+1)))/oe(k+1);0 0 1];
         yT1 = logm(inv(UE)*g1);
         EE1 = [yT1(1,3) yT1(2,3) yT1(2,1)]*inv(cov2)*[yT1(1,3);yT1(2,3);yT1(2,1)];
         fg1(k) = exp((-1*EE1)/2)/CsE;
         s = s+(fg1(k)+fg(k))*stp/2;
         end
         fxy(i,j) = s;
     end
 end

 figure;
zmax = fix(max(max(fxy)));
%zmin = fix(min(min(fxy)));
zmin = 0.05;
L = [zmin,(zmin+zmax)/2,zmax];
contour(XE,YE,fxy,L,'k','LineWidth',1.5);
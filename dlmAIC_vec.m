function AIC = dlmAIC_vec(x_panel,y_panel,k,M,b)
% k is the reduced frequency
% M is the Mach number
% b is the characteristic half-chord (used for non-dimensionalization, if 
%   x_panel and y_panel are already non-dimensionlized, set b=1)

% AIC matrix coefficients
AIC = w_doublet(x_panel(:,3),y_panel(:,3),x_panel,y_panel,b,M,k,1);
AIC = AIC + w_doublet(x_panel(:,3),y_panel(:,3),x_panel(:,[2,1,3]),-y_panel(:,[2,1,3]),b,M,k,1);

end

function res = w_doublet(x,y,x_panel,y_panel,b,M,k,Cp)
% (x,y) are the coordinates of the  point where the velocity is evaluated
% (x_panel,y_panel) are the panel points where the source doublet is placed
% c is the characteristic chord (used for non-dimensionalization)
% M is the Mach number
% k is the reduced frequency
% Cp is the pressure difference coefficient
% The output velocity will be non-dimensional

% Compute the kernels
rx_0 = (x-(x_panel(:,1)+x_panel(:,2)).'/2)/b;
ry_0 = (y-(y_panel(:,1)+y_panel(:,2)).'/2)/b;
Kbar_0 = Kbar(rx_0,ry_0,M,k);
rx_1 = (x-x_panel(:,1).')/b;
ry_1 = (y-y_panel(:,1).')/b;
Kbar_1 = Kbar(rx_1,ry_1,M,k);
rx_2 = (x-x_panel(:,2).')/b;
ry_2 = (y-y_panel(:,2).')/b;
Kbar_2 = Kbar(rx_2,ry_2,M,k);

% Compute parameters
c_k = (2*x_panel(:,3)-x_panel(:,1)-x_panel(:,2)).'./b;
b_k = sqrt((x_panel(:,2)-x_panel(:,1)).^2 + (y_panel(:,2)-y_panel(:,1)).^2).'./b;
sinL_k = (y_panel(:,2)-y_panel(:,1)).'./b_k./b;

% Compute A coefficients
A_0 = Kbar_0;
A_1 = (Kbar_2-Kbar_1)./b_k;
A_2 = (2*Kbar_1-4*Kbar_0+2*Kbar_2)./b_k.^2;

% Compute B coefficients
fact0 = 4.*c_k.*b_k./(4.*ry_0.^2-b_k.^2.*sinL_k.^2);
fact1 = log((b_k.*sinL_k-2.*ry_0).^2./(b_k.*sinL_k+2.*ry_0).^2);
B_0 = fact0.*A_0;
B_1 = (c_k./(2.*sinL_k.^2)).*fact1.*A_1 ...
    + (ry_0./sinL_k).*fact0.*A_1;
B_2 = (c_k.*b_k./sinL_k.^2).*A_2 + ...
    + (ry_0.*c_k./sinL_k.^3).*fact1.*A_2 ...
    + (ry_0.^2./sinL_k.^2).*fact0.*A_2;

% Compute induced velocity
res = -Cp./(8*pi).*(B_0 + B_1 + B_2);

end

function res = Kbar(rx,ry,M,k)
% (rx, ry) are the components of the vector from a point in the doublet 
%   segment to a point (x,y) in the z=0 plane. The input is assumed
%   non-dimensionalized by the half-chord c/2.
% M is the Mach number
% k is the reduced frequency

% Tolerances
tol = 1e-6;
big = 1e20;

% Coefficients k1 and u1 -> Eqs. (280) and (281) in Ref. [1]
k1 = k*abs(ry);
u1 = zeros(size(ry));
ind1 = abs(ry)<tol;
ind2 = ind1 & rx>0;
ind3 = ind1 & ~ind2;
ind4 = ~ind1;
u1(ind2) = -big;
u1(ind3) = big;
u1(ind4) = (M.*sqrt(rx(ind4).^2+(1-M.^2)*ry(ind4).^2)-rx(ind4))./(abs(ry(ind4)).*(1-M.^2));

% Parameters for integral quadrature -> Page 89 in Ref. [1]
c_ = 0.372;
n_(1,1,:) = 1:11;
a_(1,1,:) = [
    0.24186198
   -2.7918027
    24.991079
   -111.59196
    271.43549
   -305.75288
   -41.183630
    545.98537
   -644.78155
    328.72755
   -64.279511
];

% Evaluation of J1 integral (only valid for u1>0) -> Eq. (274) in Ref. [1]
J1 = @(u1,k1) sum(a_.*exp(-n_.*c_.*u1)./(n_.^2.*c_.^2+k1.^2).*(n_.*c_-1i.*k1),3);

% Evaluation of I1 integral -> Eq. (271) in Ref. [1]
I1_fun = @(u1,k1) exp(-1i.*k1.*u1).*(1-u1./sqrt(1+u1.^2)+(-1i.*k1.*J1(u1,k1)));
ind1 = u1>=0;
ind2 = ind1 & u1>big;
ind3 = ind1 & ~ind2;
ind4 = ~ind1;
ind5 = ind4 & u1<-big;
ind6 = ind4 & ~ind5;
I1 = zeros(size(ry));
I1(ind2) = 0;
I1(ind3) = I1_fun(u1(ind3),k1(ind3));
I1(ind5) = 2*I1_fun(0,k1(ind5));
I1(ind6) = 2*real(I1_fun(0,k1(ind6))) - conj(I1_fun(-u1(ind6),k1(ind6)));

% Evaluation of the K1 parameter -> Eq. (278) in Ref. [1]
ind1 = u1>=big | u1<=-big;
ind2 = ~ind1;
K1 = zeros(size(ry));
K1(ind1) = -I1(ind1);
K1(ind2) = -I1(ind2)-M.*abs(ry(ind2))./sqrt(rx(ind2).^2+(1-M.^2).*ry(ind2).^2).*exp(-1i.*k1(ind2).*u1(ind2))./sqrt(1+u1(ind2).^2);

% Evaluation of Kbar -> Eq. (285) in Ref. [1] 
res = K1.*exp(-1i.*k.*rx);

en
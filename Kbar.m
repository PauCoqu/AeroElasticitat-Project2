function res = Kbar(rx,ry,M,k)

% Tolerances
tol = 1e-6;
big = 1e20;

% Coefficients k1 and u1 -> Eqs. (280) and (281) in Ref. [1]
k1 = k*abs(ry);
if abs(ry)<tol
    if rx>0
        u1 = -big;
    else
        u1 = big;
    end
else
    u1 = (M*sqrt(rx^2+(1-M^2)*ry^2)-rx)/(abs(ry)*(1-M^2));
end

% Parameters for integral quadrature -> Page 89 in Ref. [1]
c_ = 0.372;
n_ = (1:11)';
a_ = [
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
J1 = @(u1,k1) sum(a_.*exp(-n_*c_*u1)./(n_.^2*c_^2+k1^2).*(n_*c_-1i*k1));

% Evaluation of I1 integral -> Eq. (271) in Ref. [1]
I1_fun = @(u1,k1) exp(-1i*k1*u1)*(1-u1/sqrt(1+u1^2)+(-1i*k1*J1(u1,k1)));
if u1>=0
    if u1>big
        I1 = 0;
    else
        I1 = I1_fun(u1,k1);
    end
else % For u1<0 -> Eq. (275) in Ref. [1]
    if u1<-big
        I1 = 2*I1_fun(0,k1);
    else
        I1 = 2*real(I1_fun(0,k1)) - real(I1_fun(-u1,k1)) + 1i*imag(I1_fun(-u1,k1));
    end
end

% Evaluation of the K1 parameter -> Eq. (278) in Ref. [1]
if u1>=big || u1<=-big
    K1 = -I1;
else
    K1 = -I1-M*abs(ry)/sqrt(rx^2+(1-M^2)*ry^2)*exp(-1i*k1*u1)/sqrt(1+u1^2);
end

% Evaluation of Kbar -> Eq. (285) in Ref. [1] 
res = K1*exp(-1i*k*rx); 
end
function AIC = buildAIC(x_p, y_p, k, M_inf, b_semi)
% buildAIC
% Aerodynamic Influence Coefficient matrix using w_doublet
%
% INPUTS
%   x_p    : [N_panels x 3] x-coordinates (qc start, qc end, collocation)
%   y_p    : [N_panels x 3] y-coordinates (qc start, qc end, collocation)
%   k      : reduced frequency
%   M_inf  : Mach number
%   b_semi : semi-chord (c/2)
%
% OUTPUT
%   AIC    : [N_panels x N_panels]

N_panels = size(x_p,1);
AIC = zeros(N_panels, N_panels);

for i = 1:N_panels
    % Collocation point of panel i
    x_i = x_p(i,3);
    y_i = y_p(i,3);

    for j = 1:N_panels
        % Doublet contribution from panel j
        AIC(i,j) = AIC(i,j) + w_doublet(x_i,y_i,x_p(j,:),y_p(j,:),b_semi,M_inf,k,1);
        AIC(i,j) = AIC(i,j) + w_doublet(x_i,y_i,x_p(j,[2,1,3]),-y_p(j,[2,1,3]),b_semi,M_inf,k,1);
    end
end

end
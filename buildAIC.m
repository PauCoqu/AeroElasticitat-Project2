function AIC = buildAIC(x_p,y_p,c_root,lambda,y0,b,M_inf,k)
% Same AIC you already coded, but as a function so gust loop can call it with varying k.
N_panels = size(x_p,1);
AIC = zeros(N_panels,N_panels);
for i = 1:N_panels
    x_i = x_p(i,3);
    y_i = y_p(i,3);
    for j = 1:N_panels
        y_mid_j = y_p(j,3);
        c_local = c_root * (1 - (1 - lambda) * (y_mid_j - y0)/b);
        AIC(i,j) = AIC(i,j) + w_doublet(x_i,y_i,x_p(j,:),y_p(j,:),c_local/2,M_inf,k,1);
        AIC(i,j) = AIC(i,j) + w_doublet(x_i,y_i,x_p(j,[2,1,3]),-y_p(j,[2,1,3]),c_local/2,M_inf,k,1);
    end
end
end
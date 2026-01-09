function AIC = AIC_builder(k_, x_p, y_p, N_panels, c_root, taper, y0, b, M_inf)
% AIC_builder  Build DLM AIC matrix for a given reduced frequency k

AIC = complex(zeros(N_panels, N_panels));

for ii = 1:N_panels
    x_ci = x_p(ii,3);
    y_ci = y_p(ii,3);

    for jj = 1:N_panels
        y_mid_j = y_p(jj,3);
        c_loc   = c_root * (1 - (1 - taper) * (y_mid_j - y0)/b);

        AIC(ii,jj) = AIC(ii,jj) + w_doublet(x_ci, y_ci, x_p(jj,:),        y_p(jj,:),         c_loc/2, M_inf, k_, 1);
        AIC(ii,jj) = AIC(ii,jj) + w_doublet(x_ci, y_ci, x_p(jj,[2,1,3]), -y_p(jj,[2,1,3]),    c_loc/2, M_inf, k_, 1);
    end
end

end

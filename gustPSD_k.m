function Phi_k = gustPSD_k(k, sigma_g, k0)
    % Piecewise PSD given by the project statement in reduced frequency k
    % Output has units [(m/s)^2] per "unit k" (your subsequent dk/dw conversion handles omega-domain)
    Phi_k = zeros(size(k));
    mask_low = k < k0;
    mask_high = ~mask_low;
    
    Phi_k(mask_low)  = sigma_g^2 * 10^(8/3);
    Phi_k(mask_high) = sigma_g^2 * 10^(-8/3) * k(mask_high).^(-12/7);
end
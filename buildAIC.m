function AIC = buildAIC(x_panel, y_panel, c_root, M, k)
% buildAIC  Assemble DLM AIC matrix using Project 1 w_doublet/Kbar
% c_root is the root chord; we use c_root/2 as reference semi-chord (same as P1)

N_panels = size(x_panel,1);
AIC = zeros(N_panels,N_panels);

for i = 1:N_panels
    x_i = x_panel(i,3);
    y_i = y_panel(i,3);
    for j = 1:N_panels
        AIC(i,j) = AIC(i,j) + w_doublet(x_i,y_i,x_panel(j,:),y_panel(j,:),c_root/2,M,k,1);
        AIC(i,j) = AIC(i,j) + w_doublet(x_i,y_i,x_panel(j,[2,1,3]),-y_panel(j,[2,1,3]),c_root/2,M,k,1);
    end
end
end
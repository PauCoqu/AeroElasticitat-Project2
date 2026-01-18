function sigma_max = sigmaMax_fromResponse(q_full, nodes, span, material, x_sc_line)
    % Compute sigma_max over nodes using the project statement formula
    E = material.YoungModulus;
    nu = material.Poisson;
    G = E / (2 * (1 + nu));
    y = span.y(:);
    c_y = span.c_y(:);
    % Extract DOFs at nodes
    N_nodes = size(nodes,1);
    eta   = q_full(1:3:end);
    zeta  = q_full(2:3:end);
    theta = q_full(3:3:end);
    % For each span station (y(j)), average zeta and theta across chord nodes
    zeta_y  = zeros(length(y),1);
    theta_y = zeros(length(y),1);
    for j = 1:length(y)
        idx = find(abs(nodes(:,2) - y(j)) < 1e-12);
        zeta_y(j)  = mean(zeta(idx));
        theta_y(j) = mean(theta(idx));
    end
    % Spanwise derivatives (central differences or forward/backward at ends)
    dzeta_dy  = gradient(zeta_y, y);
    dtheta_dy = gradient(theta_y, y);
    % Evaluate sigma at each node
    sigma_node = zeros(N_nodes,1);
    for n = 1:N_nodes
        % station index in y-grid
        [~,j] = min(abs(y - nodes(n,2)));
        % thickness model (consistent with your plate element)
        h_local = 0.05 * c_y(j);  % or use your actual thickness distribution
        % distance to shear center in x
        r_sc = abs(nodes(n,1) - x_sc_line);
        term_bend = (E * h_local / 2) * dzeta_dy(j);
        term_tors = (sqrt(3) * G * r_sc) * dtheta_dy(j);
        sigma_node(n) = sqrt( term_bend^2 + term_tors^2 );
    end
    sigma_max = max(sigma_node);
end
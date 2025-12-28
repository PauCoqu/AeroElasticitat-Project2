%EXTRES

%PLOT DE NODES I ELEMENTS AMB ELS CORRESPONENTS NUMEROS (AL MODAL_ANALYSIS DESPRES DE NODES I ELEMENTS):

figure; hold on; axis equal; grid on

% Elements
for e = 1:size(elem,1)
    xy = nodes(elem(e,:),:);
    plot([xy(:,1); xy(1,1)], [xy(:,2); xy(1,2)], 'k-');

    % número d’element (al centre)
    xc = mean(xy(:,1));
    yc = mean(xy(:,2));
    text(xc, yc, num2str(e), 'Color','r','FontSize',8,...
        'HorizontalAlignment','center');
end

% Nodes
plot(nodes(:,1), nodes(:,2), 'b.', 'MarkerSize', 12)

for k = 1:size(nodes,1)
    text(nodes(k,1), nodes(k,2), [' ' num2str(k)], ...
        'Color','b','FontSize',8);
end

xlabel('x'); ylabel('y');
title('Wing structural mesh (nodes + elements)');

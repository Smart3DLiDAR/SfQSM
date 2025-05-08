function P = lineExtract(P)
SHOW_RESULTS = true;

options.collapse_order = 1;     %控制折叠的权重项，是否加上惩罚项,若点云点间距比较小，可以没有惩罚项，反之，加上惩罚项效果会好一些
[P.spls,P.corresp] = farthestSamplingBySphere(P.cpts, P.sample_radius);
P.spls_adj = connectByInheritNeighbor(P.cpts, P.spls, P.corresp, P.rings);    
[P.spls, P.spls_adj,P.corresp] = edgeCollapse(P.spls, P.corresp, P.spls_adj, options);

%%
if (SHOW_RESULTS)
    figure; movegui('northeast');set(gcf,'color','white');hold on;
    plot3( P.spls(:,1), P.spls(:,2), P.spls(:,3), '.r', 'markersize', 5);
    axis off; axis equal;set(gcf,'Renderer','OpenGL');
    plotConnectivity(P.spls, P.spls_adj);
end

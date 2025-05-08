function P = lineExtract(P)
SHOW_RESULTS = true;

options.collapse_order = 1;     %�����۵���Ȩ����Ƿ���ϳͷ���,�����Ƶ���Ƚ�С������û�гͷ����֮�����ϳͷ���Ч�����һЩ
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

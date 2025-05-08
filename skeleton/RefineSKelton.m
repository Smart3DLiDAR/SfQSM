function P = RefineSKelton(P,DebugShow)
%%
%SHOW_JOINTS = true;
%SHOW_ROOT_JOINT = true;
%SHOW_CYCLES = true;
%SHOW_IRRELEVANT_EXTRAMA = true;

%t1 = 0.1; % for inner branch nodes
%a1 = pi*5.0/7.0; % for inner branch nodes,  
%t2 = 0.2; % for irrelevant extrama;
%t3 = 20; % for small cycles;
%sprintf('thresholds for remove 1)inner nodes: %f angle, %f length \n 2)irrelevant extrama: %f\n3)small cycles: %d \n',t1, a1, t2, t3)

%load(sk_filename,'P');
%axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-90,0);
%% 0: 首先找到结合处
%[joints, roots,segments] = find_joints(P.pts, P.spls, P.corresp, P.spls_adj, SHOW_JOINTS);
%% 1: removing small cycles measured by topological length
%[P.spls, P.corresp, P.spls_adj, joints, segments] = remove_small_cycles(P.spls, P.corresp, P.spls_adj, joints, segments,t3, SHOW_CYCLES);
%% 2: 查找根节点，相对于root_id的全局距离，“骨架大小”
%[root_id, global_dist] = find_root_node(P.spls, P.spls_adj, joints, SHOW_ROOT_JOINT);
%% 3: 去除无关的分支点
%[P.spls, P.corresp, P.spls_adj, joints, segments] = remove_irrelevant_extrama(P.spls, P.corresp, P.spls_adj, joints, segments,global_dist, t2, SHOW_IRRELEVANT_EXTRAMA);
%% 4: 删除内部节点
%[P.spls, P.corresp, P.spls_adj, joints, segments] = remove_inner_nodes(P.spls, P.corresp, P.spls_adj, joints, segments,global_dist,t1,a1);
%[P.spls, P.spls_adj] = remove_upper_nodes(P.spls,P.spls_adj,true);% delate fake nodes;
%% 5: 更新连接图
%P.root_id = root_id;% 如果没有分支点，则第一个节点是根节点。
[~,P.root_id] = min(P.spls(:,3));
%[P, graph] = build_hierarchy(P,global_dist); % do not contain cycles
[P.spls, P.corresp, P.spls_adj, Graph] = build_graph(P.spls,P.corresp,P.spls_adj);% may contain cycles
%% 6: 将骨架节点移动到其对应局部点云领域的质心
%figure('Name','Build hierarchy','NumberTitle','off');
%set(gcf,'color','white');hold on; movegui('south');
%plot_skeleton(P.spls,P.spls_adj);
%scatter3(P.pts(:,1),P.pts(:,2),P.pts(:,3),60,'.','MarkerEdgeColor', GS.PC_COLOR);hold on;
for i = 1:size(P.spls,1)
    verts = P.pts(P.corresp==i,:);
    if size(verts,1) == 1
        P.spls(i,:) = verts;
    else
        P.spls(i,:) = mean(verts);
    end
end
for i = 1:size(P.spls_adj,1)
    links = find( P.spls_adj(i,:)==1 );
    if length(links) == 2
        verts = P.spls(links,:);
        P.spls(i,:) = mean(verts);
    end
end
%figure;set(gcf,'color','white');
%plot_embedded_graph(P.spls, Graph, 'b', 'linewidth', 2);
%axis off; axis equal; camorbit(0,0,'camera');  view(-90,0);
%% 最终骨架
%% finally
if(DebugShow)
figure;set(gcf,'color','white');movegui('northwest');set(gcf,'Renderer','OpenGL');
%scatter3(P.pts(:,1),P.pts(:,2),P.pts(:,3),20,'.');hold on;
scatter3(P.spls(:,1),P.spls(:,2),P.spls(:,3),300,'.');hold on;
plot_skeleton(P.spls, P.spls_adj);hold on;
axis off; axis equal; camorbit(0,0,'camera');view(0,90);
end

end
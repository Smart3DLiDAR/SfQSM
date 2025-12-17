function [ P_S ] = skeleton_laplacian( P,Parameters )
%SKELETON_LAPLACIAN 此处显示有关此函数的摘要
%   此处显示详细说明
options.USING_POINT_RING = GS.USING_POINT_RING;
P.npts = size(P.pts,1);
[P.bbox, P.diameter] = GS.compute_bbox(P.pts);
P.k_knn =Parameters.KnnNum;
if options.USING_POINT_RING
    P.rings = compute_point_point_ring(P.pts, P.k_knn, []);
else    
    P.frings = compute_vertex_face_ring(P.faces);
    P.rings = compute_vertex_ring(P.faces, P.frings);
end
[P.cpts, t, initWL, WC, sl] = contraction_by_mesh_laplacian(P, options);
P.sample_radius = P.diameter*Parameters.sampleScale;
P = rosa_lineextract(P,P.sample_radius, 1);


SHOW_JOINTS = false;
SHOW_ROOT_JOINT = false;
SHOW_CYCLES = false;
SHOW_IRRELEVANT_EXTRAMA = false;

t1 = Parameters.t1; % for inner branch nodes
a1 = Parameters.a1; % for inner branch nodes,  
t2 = Parameters.t2; % for irrelevant extrama;
t3 = Parameters.t3; % for small cycles;

[joints, roots,segments] = find_joints(P.pts, P.spls, P.corresp, P.spls_adj, SHOW_JOINTS);
%% 1: removing small cycles measured by topological length
[P.spls, P.corresp, P.spls_adj, joints, segments] = remove_small_cycles_mt(P.spls, P.corresp, P.spls_adj, joints,roots, segments,t3, SHOW_CYCLES);
%% 2: find root joint, global distance relative to root_id, "size of skeleton"
[root_id, global_dist] = find_root_node(P.spls, P.spls_adj, joints, SHOW_ROOT_JOINT);
%% 3: unify joints if there are no branch node inbetween of them. The effect is bad for some cases.
if ~isempty(joints)
    [P.spls, P.corresp, P.spls_adj, joints, root_id] = merge_nearby_joints(P.spls, P.corresp, P.spls_adj, joints, root_id);
end
%linshi
%
%% 4: removing irrelevant extrama
[P.spls, P.corresp, P.spls_adj, joints, segments] = remove_irrelevant_extrama(P.spls, P.corresp, P.spls_adj, joints, segments,global_dist, t2, SHOW_IRRELEVANT_EXTRAMA);
% P.root_id = root_id;% The first node is root node if there is no joints.
[P.spls,P.corresp,P.spls_adj, graph] = build_graph(P.spls,P.corresp,P.spls_adj);% may contain cycles
% 
% 
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
[joints ,roots, branches]=find_Joints_mt(P.spls, P.spls_adj,P.pts, true); 
P.joints=joints;
P.roots=roots;
P.branches=branches;
P_S=P;

end


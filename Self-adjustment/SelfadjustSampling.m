function P = SelfadjustSampling(P, R)
%% Skeleton points obtained with adaptive farthest distance sphere sampling
% The sampling function based on the linear characteristics of the contracted point cloud 
% is used to adaptively adjust the radius of the sampling sphere
% --input------------------------------------------------------------
% P: The information of individual tree point clouds 
%   P.cps : The contraction point clouds
%   P.rings : The one-ring of each individual tree point cloud
% R: Initial radius
% --output------------------------------------------------------------
% P.spls : Skeleton points
% P.corresp : The correspondence between point clouds and skeleton points
% P.seg : The set of contraction points corresponding to every skeleton points
%%
[~,Linearity] = define_contractionnodeshape(P.cpts,P.rings);
[P.spls, P.corresp] = FarthestSamplingbysphere(P,Linearity,R);
P.spls_adj = connect_by_inherit_neigh(P.cpts, P.spls, P.corresp, P.rings);

% Each raw point cloud corresponds to a skeleton point
% According to the row indexes with the same corresponding relationship, 
% the contracted point sets corresponding to each skeleton point can be found 
t = accumarray(P.corresp,[1:length(P.corresp)]',[],@(x) {x});
Seg = cellfun(@(x) P.pts(x,:),t,'UniformOutput',0);
P.seg = Seg;

end

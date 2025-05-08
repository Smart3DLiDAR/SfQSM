function [nodelabel,Linearity] = define_contractionnodeshape(cpts,rings)
%% Calculate the geometric features of each point and its Delaunay neighbors 
% to determine whether to distinguish the points at bifurcation or branch parts
% --input------------------------------------------------------------
% cps : The contraction point clouds
% rings : The one-ring of each individual tree point cloud
% --output------------------------------------------------------------
% nodelabel : the points at bifurcation or branch parts
%      nodelabel==1 : the points at branch parts
%      nodelabel==0 : the points at bifurcation parts
% Linearity : the linearity of each contraction point cloud
%%
nodelabel = zeros(size(cpts,1),1);
Linearity = zeros(size(cpts,1),1);
for i= 1:size(cpts,1)
    curring = rings{i,1};
    if length(curring)>=3
        neibor =cpts(curring,:);
        [~,~,latent] = pca(neibor, 'Economy', false); 
        Linearity = latent(1)/sum(latent);
        if Linearity >= 0.6   
            nodelabel(i) = 1;
        end
    else
        kdtree = KDTreeSearcher(cpts);
        id = knnsearch(kdtree,cpts(i,:),'k',5);
        neiber = cpts(id,:);
        [~,~,latent] = pca(neiber, 'Economy', false); 
        Linearity = latent(1)/sum(latent);
        if Linearity >= 0.6 
            nodelabel(i) = 1;
        end
    end
    Linearity(i) = Linearity;
end

end
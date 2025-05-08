function [Gpdis,adjdis] = BuildGraph(spls,Seg,R)
% --input------------------------------------------------------------------
% spls : skeleton points
% Seg : the contracted point sets of each skeleton point
% R : radius
% --output-----------------------------------------------------------------
% Gpdis : graph
% adjdis : adjacency matrix

% The code is modified from Wangdi

%% Obtain adjacency information
npts = size(spls,1);
kdtree = KDTreeSearcher(spls);
[neighbors1,~] = rangesearch(kdtree,spls,R);

% Calculate the shortest distance of the contraction point clusters 
% corresponding to adjacent skeleton points
dis = cell(npts,1);
for i = 1:length(neighbors1)
    curring = neighbors1{i,1};
    tempdis = zeros(1,length(curring));
    for j = 1:length(curring)
        temp = pdist2(Seg{i,1},Seg{curring(j),1});
        tempdis(j) = min(temp(:));
    end
    dis{i,1} = tempdis;
end

%% Obtain the source and target of each edge
source = [];
target = [];
weight = [];
for i = 1:size(neighbors1,1)%two
    curring = neighbors1{i,1}(2:end);
    curdistance1 = dis{i,1}(2:end);
    
    numring = length(curring);
    tempsourcering = repmat(i, [1 numring])';
    temptargetring = curring';
    

    tempweightring = zeros(numring,1);
    for j = 1:numring
        tempweightring(j,1) = curdistance1(j);
    end
    
    source = [source;tempsourcering];
    target = [target;temptargetring]; 
    weight = [weight;tempweightring];
end

%% Pruning Edge
tempdis = [];
pruned = [];
mediandis = [];
maxdis = [];
for i = 1:size(neighbors1,1)
    curring = neighbors1{i,1}(2:end);
    curdistance1 = dis{i,1}(2:end);
    tempdis = [tempdis,curdistance1];
    
    % Condition1£ºRemove the edges whose distance is longer than (average value + std)
    dt = mean(curdistance1) + std(curdistance1);
    temppruned = curdistance1' > dt;
    pruned = logical([pruned;temppruned]);
    
    tempmax = max(curdistance1);
    maxdis = [maxdis,tempmax];
end
% Condition2£ºDefine the maximum distance£¬
% and remove the edges whose distance is longer than (maxd+std)
maxd = mean(maxdis) + std(maxdis);
pruned2 = tempdis' > maxd;

% Pruning Edge
selfedge = source==target;% self edge
to_remove = selfedge + pruned + pruned2;% All edges that need to be removed
source = source(~to_remove);
target = target(~to_remove);
weight = weight(~to_remove);

%% Build Graph
Edge = [source,target];

A = zeros(size(spls,1),size(spls,1));
for i = 1:size(Edge,1)
    A(Edge(i,1),Edge(i,2)) = 1;
    A(Edge(i,2),Edge(i,1)) = 1;
end

% Judge the graph connectivity
[bin,binsize] = conncomp(digraph(A),'Type','weak','OutputForm','cell');
[~,rootid] = min(spls(:,3));
for i = 1:length(binsize)
    curbin = bin{1,i};
    % Find the subgraph where the root node is located
    [~,col] = find(curbin==rootid);
    if ~isempty(col)
        level = i;
        MainSkelton = bin{1,i};
    end
end

% Define the importance of subgraphs based on the number of nodes in the subgraph 
% and sort them according to their importance
[binsize,ia]=sort(binsize,'descend');
newbin = {};
newbin{1,1} = bin{1,level};
for i= 1:length(ia)
    if ia(i)~= level
        newbin{end+1,1} = bin{1,ia(i)};
    end
end

% Connect the subgraphs to the main skeleton graph in sequence
for i = 2:length(binsize)
    cursubgraph = newbin{i,1};
    [~,minid] = min(spls(cursubgraph,3));

    if length(MainSkelton)>1
        [neik,~] = knnsearch(spls(MainSkelton,1:3),spls(cursubgraph(minid),1:3),'k',3);
        neik = MainSkelton(neik);
    else
        neik = MainSkelton;
    end

    for j = 1:length(neik)
        Edge(end+1,1:2) = [cursubgraph(minid) neik(j)];
        tempdis = pdist2(Seg{cursubgraph(minid),1},Seg{neik(j),1});
        weight(end+1) = min(tempdis(:));
    end
    MainSkelton = unique([MainSkelton,cursubgraph],'stable');
end

% Remove repetition
[Edge,ia,~] = unique(Edge,'rows');
Edge2 = [Edge(:,2),Edge(:,1)];
finalEdge = [Edge;Edge2];
weight = weight(ia,:);
weight = [weight;weight];

% build graph
adjdis = sparse(finalEdge(:,1),finalEdge(:,2),weight,size(spls,1),size(spls,1));
a = table((1:size(spls,1))');
Gpdis = graph(adjdis,a);

end
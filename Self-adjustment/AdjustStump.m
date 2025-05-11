function P = AdjustStump(P)
%% ---1.1 部分树种基部存在很多分枝,进行树桩调整----------------------------------------
joints = zeros(0,1);
roots = zeros(0,1);
for i=1:size(P.spls,1)
    P.spls_adj(i,i) = 0; 
    links = find( P.spls_adj(i,:)==1 );
    if length(links) >2 
        joints(end+1,1)=i;
    elseif length(links) ==1 
        roots(end+1,1)=i;
    end
end

[jzmin,~] = min(P.spls(joints(:),3));
[~,rootmin] = min(P.spls(roots(:),3));

%%--Calculate the tree height to define the reference plane----------------------------------------------------------------------------
[zmax,~] = max(P.pts(:,3));
[zmin,zminid] = min(P.spls(:,3));
treeH = zmax - zmin;
if treeH >= 5 
    TH = zmin + 1.3;
else
    TH = 1.3 * (treeH/5);
end

% Find the position of the root node
% Judge whether the root needs to be adjusted
if jzmin <= TH 
    [lowid,~] = find(P.spls(joints(:),3) < TH );
    [~,maxjid] = max(P.spls(joints(lowid),3));
    curjoint = P.spls(joints(lowid(maxjid)),:);
    curjointid = joints(lowid(maxjid));
    [lowjid,~] = find(P.spls(:,3)< curjoint(:,3)+1);
    [~,sortid] = find(P.spls_adj(curjointid,:) == 1);
    id1 = unique([lowjid;sortid']);
    
    for i= 1:size(id1,1)
        [ROW,~] = find(P.corresp(:)==id1(i));
        points = P.pts(ROW,:);
        figure(2);set(gcf,'color','white');
        set(gca,'position',[0.1,0.1,0.8,0.8]);
        scatter3(points(:,1),points(:,2),points(:,3),3,'.b'); hold on;
        axis off;axis equal;
    end
    
    %%----------------------------------------------------------------------------------------
    % According to the physiological characteristics of plants, 
    % the main trunk is responsible for the upward transport of water and nutrients from the roots.
    % The roots actively absorb water and nutrients in a downward-oriented direction
    % Define the direction vectors V(z) that denote the direction 
    % from the bifurcation skeleton points at root towards other skeleton points(pi)
    
    % Condition 1 : V(z)>0 means this skeleton points(pi) is trunk skeleton poins
    % Condition 2 : V(z)<0 means this skeleton points(pi) is root skeleton poins
   
    v = zeros(size(id1,1),3);
    for i = 1:length(id1)
        v(i,1:3) = P.spls(id1(i),:) - P.spls(curjointid,:);
    end
    
    if ~isempty(id1)
        rootsection = [];
        for j = 1:length(id1)
            [temprow,~] = find(P.corresp(:) == id1(j));
            tempsection = P.pts(temprow,:); 
            rootsection = [rootsection;tempsection];
        end
        XY = rootsection(:,1:2);
        [THETA,RHO] = pertransform(XY);
        [W, ~, ptsC] = dbscan2105([THETA,RHO]', 0.25, 3);
        gapx = zeros(length(W),1);
        for j = 1:max(ptsC)
            curdbs = W{1,j};
            gapx(j) = max(curdbs(:,1)) - min(curdbs(:,1));
        end
        ratio = (sum(gapx))/(2*pi);
    
        if ratio >= 2/3
            [tempxc,tempyc,~] = circlefit(rootsection);
            if (tempxc > min(XY(:,1)) && tempxc < max(XY(:,1))) && (tempyc > min(XY(:,2)) && tempyc < max(XY(:,2)))
                P.spls(id1(1),1:3) = [tempxc,tempyc,min(P.spls(id1,3))];
            else 
                P.spls(id1(1),1:3) = [mean(P.spls(id1,1:2)),min(P.spls(id1,3))];
            end
        else
            Par = CircleFitByTaubin(XY);
            if (Par(:,1) > min(XY(:,1)) && Par(:,1) < max(XY(:,1))) && (Par(:,2) > min(XY(:,2)) && Par(:,2) < max(XY(:,2)))
                P.spls(id1(1),1:3) = [Par(:,1),Par(:,2),min(P.spls(id1,3))];
            else 
                P.spls(id1(1),1:3) = [mean(P.spls(id1,1:2)),min(P.spls(id1,3))];
            end
        end
        
        for j = 2:length(id1)
            P.spls( id1(j),: ) = NaN;
            P.corresp( P.corresp==id1(j) ) = id1(1);
            for k=1:size(P.spls_adj,1)
                if P.spls_adj(id1(j),k) == 1
                    P.spls_adj(id1(1),k)=1;
                    P.spls_adj(k,id1(1))=1;
                end
            end
            P.spls_adj(id1(j),:) = 0;
            P.spls_adj(:,id1(j)) = 0;
        end
        Rootpoint = P.spls(id1(1),1:3);
    else 
        Rootpoint = curjoint;
    end
else 
    roota = P.spls(roots(rootmin),:);
    rootb = P.spls(zminid,:);
    if roota == rootb
        Rootpoint = rootb;
    else
        Rootpoint = roota;
    end
end
%% 
[~,rootid] = min(P.spls(:,3));
[nei,~] = knnsearch(P.spls,P.spls(rootid,1:3),'k',3);
spls = P.spls(nei,1:3);
[x1, y1] = bezir_n(spls(:,1:2)', 7);
[~, z2] = bezir_n(spls(:,[1 3])', 7);
newP1 = [x1' y1' z2'];
newP = newP1(2:end-1,1:3);
for j = 1:size(newP,1)
    P.spls(end+1,1:3) = newP(j,:);
    idex = size(P.spls,1);
    if j == 1
        P.spls_adj(idex,nei(1)) = 1;
        P.spls_adj(nei(1),idex) = 1;
    elseif j ~= 1 && j ~= size(newP,1)
        [~,ia,~] = intersect(P.spls,newP(j-1,:),'rows');
        P.spls_adj(idex,ia) = 1;
        P.spls_adj(ia,idex) = 1;
    else
        [~,ia,~] = intersect(P.spls,newP(j-1,:),'rows');
        P.spls_adj(idex,[ia;nei(3)]) = 1;
        P.spls_adj([ia;nei(3)],idex) = 1;
    end
    [row,~] = find(P.pts(:,3)<newP(j,3)&P.pts(:,3)>=newP1(j,3));
    P.corresp(row)= idex;
end

Graph = zeros(0,2);
for i=1:size(P.spls_adj,1)
    for j=i+1:size(P.spls_adj,2)
        if( P.spls_adj(i,j)==1 )
            Graph(end+1,:) = [i, j];
        end
    end
end
tmp = find(isnan(P.spls(:,1)))';
Tmp = tmp(length(tmp):-1:1); 
for i=Tmp
    Graph(Graph>i) = Graph(Graph>i) - 1;
    P.corresp(P.corresp>i) = P.corresp(P.corresp>i) -1;
end
P.spls(tmp,:) = [];
spls_adj = zeros(size(P.spls,1),size(P.spls,1));
for a=1:size(Graph,1)
    i = Graph(a,1); j = Graph(a,2);
    spls_adj(i,j) = 1;
    spls_adj(j,i) = 1;
end
P.spls_adj = spls_adj;

%%
t = accumarray(P.corresp,[1:length(P.corresp)]',[],@(x) {x});
Seg = cellfun(@(x) P.pts(x,:),t,'UniformOutput',0);
P.seg = Seg;

end
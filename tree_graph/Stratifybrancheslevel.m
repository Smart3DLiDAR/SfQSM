function [Branch,branch] = Stratifybrancheslevel(spls,adj,Cid,pli,rootid)
%%----Divide the hierarchy of branches ----------------------------------------------------------
%---input -----------------------------------------------------------------
%  spls: Skeleton point
%  adj : Skeleton point connection matrix
%  Cid : The child nodes of each skeleton point
%  pli : The occurrence times of each skeleton point in the shortest path
%---output -----------------------------------------------------------------
%  Branch : Branches at all levels
%  branch :
%% --------------------------------------------------------------------------
seedid = rootid;
Candidate = zeros(size(spls,1),1); 
mainsteam = zeros(0,1); 
for i = 1:inf 
    id = seedid;
    col = Cid(id,:);
    col(:,all(col==0,1))=[];
    newroi = intersect(mainsteam,id);
    if ~isempty(newroi)
        break
    end
    mainsteam(end+1,1) = id;

    if length(col)==1
        seedid = col;
    elseif length(col)>1
        add = [];
        for j = 1:length(col)
            [row,~] = find(pli(:,1) == col(j));
            add = [add,row];
        end
        % Compare the frequency of occurrence of each child node in the path
        % The one with the highest occurrence frequency maybe is the main branch
        [~,maxid] = max(pli(add,2)); 
        if length(maxid) == 1
            seedid = col(maxid);
        else
            seedid = col(maxid(1));
        end
        otherid = setdiff(col,seedid);
        for k = 1:length(otherid)
            Candidate(otherid(k)) = 1;
        end
    else
        break
    end
end

%%  
branch{1} = mainsteam;
Branch{1} = branch;
lastlevel = mainsteam;
I = 2;
for i = 1:inf
    cand=find(Candidate==1);
    [~,id,~] = intersect(cand,lastlevel);
    cand(id) = [];
    branch = [];
    if ~isempty(cand)
        tempbranch =[];
        for j = 1:length(cand)
            Seedid = cand(j);
            [curbranch,Candidate] = findbranch(Seedid,Cid,pli,Candidate);
            [~,roiid,~] = intersect(curbranch,lastlevel);
            curbranch(roiid) = [];
            if ~isempty(curbranch)
                branch{end+1,1}= curbranch;
                tempbranch = [tempbranch;curbranch];
            else
                branch{end+1,1}= cand(j);
                tempbranch = [tempbranch;cand(j)];
            end
        end
        Branch{I,1} = branch;
        lastlevel = unique([lastlevel;tempbranch]);
    else
        break
    end
   I = I+1;
end


joints = zeros(0,1); % Storage bifurcate node
roots = zeros(0,1); % Storage end node
for i=1:size(spls,1)
    links = find( adj(i,:)~=0 );
    if length(links) >2 
        joints(end+1,1)=i;
    elseif length(links)==1
        roots(end+1,1) = i;
    end
end
Gp1 = graph(adj);
[tR,~] = shortestpathtree(Gp1,rootid,[1:size(spls,1)],'OutputForm','cell');

% branchorder is hierarchically divided based on the number of key nodes(bifurcate node) 
% that appear in the shortest path
label = zeros(size(spls,1),2); 
branch = {};
for i= 1:length(joints)
    curtR = tR{joints(i),1};
    [roi,ia,~] = intersect(curtR,joints);
    ia = sort(ia,'ascend');
    num = length(roi);
    for j = 1:num
        if j ==1
            id1 = 1;
            tmp = logical(label(curtR(id1:ia(j))));
            if sum(tmp)==0
                curbranch = curtR(id1:ia(j));
                label(curbranch,1) = j;
                branch{end+1,1} = curbranch;
            end
        else
            id1 = ia(j-1)+1;
            tmp = logical(label(curtR(id1:ia(j))));
            if sum(tmp)==0
                curbranch = curtR(id1:ia(j));
                label(curbranch,1) = j;
                branch{end+1,1} = curbranch;
            end
        end
    end
end

roots = setdiff(roots,rootid);
for i = 1:length(roots)
    curtR = tR{roots(i),1};
    [roi,ia,~] = intersect(curtR,joints);
    if ~isempty(roi)
        ia = sort(ia,'ascend');
        num = length(roi);
        curbranch = curtR(ia(num)+1:end);
        tmp = logical(label(curbranch));
        if sum(tmp)==0
            label(curbranch,1) = j+1;
            branch{end+1,1} = curbranch;
        end
    end
end


end
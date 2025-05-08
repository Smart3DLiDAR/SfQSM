function [nodelabel,spls_adj,pli,Pid,Cid,Branch,BranchL] = FindneedAdjustNodes(spls,Seg,rootid)
%%--------------------------------------------------------------------------------------------------
[~,dis] = knnsearch(spls,spls,'k',6);
OptimalR = max(dis(:,end));

[Gp,~] = BuildGraph(spls,Seg,OptimalR);
[tR,D] = shortestpathtree(Gp,rootid,[1:size(spls,1)],'OutputForm','cell');
adj = zeros(size(spls,1),size(spls,1));
for i = 1:size(spls,1)
    curpath = tR{i,:};
    if length(curpath)>1
        front = curpath(end-1);
        adj(i,front) = 1;
        adj(front,i) =1;
    end
end
C = cell2mat(tR');
pli = tabulate(C);
[Cid,Pid] = define_topologicalrelations(spls,adj,rootid);
%%
[Branch,~] = Stratifybrancheslevel(spls,adj,Cid,pli,rootid);
BranchL = CalculatBranchLength(spls,Pid,adj,rootid,Branch);

delateb = [];
for i = 2:length(Branch)
    curBL = BranchL{i,1};
    curBranchid = Branch{i,1};
    lastBL = BranchL{i-1,1};
    lastBranchid = Branch{i-1,1};

    curlratio = zeros(1,length(curBL));
    if length(lastBL) >1
        for j = 1:length(curBL)
            curlength = curBL(j);
            curbenginid = curBranchid{j,1}(1);
            p_curbenginid = Pid(curbenginid);
            for n = 1:length(lastBranchid)
                roi = intersect(lastBranchid{n,1},p_curbenginid);
                if ~isempty(roi)
                    lastBranchL = lastBL(n);
                end
            end          
            curlratio(j) = curlength/lastBranchL;
        end
    else
        for j = 1:length(curBL)
            curlength = curBL(j);
            curlratio(j) = curlength/lastBL;
        end       
    end
    
    th = mean(curlratio)+std(curlratio);
    lastorglabel = zeros(1,length(curlratio));
    tempdelete = [];
    for m = 1:length(curlratio)
        curbranch = Branch{i,1}{m,1};
        if curlratio(m) < th 
            tempdelete = curbranch;
            lastorglabel(m) = 1;
        end
        delateb = unique([delateb;tempdelete]);
    end 
end

splsssss = spls;
splsssss(delateb,:) = NaN;
spls_adj = adj;
spls_adj(delateb,:) = 0;
spls_adj(:,delateb) = 0;

joints = zeros(0,1);
for i=1:size(splsssss,1)
    spls_adj(i,i) = 0; 
    links = find( spls_adj(i,:)==1 );
    if length(links) >2 % bifurcation points
        joints(end+1,1)=i;
    end
end

% Hierarchical division based on 
% the number of bifurcation points that appear in the shortest path
nodelabel = zeros(size(spls,1),1); 
for i= 1:length(tR)
    curtR = tR{i,1};
    roi = intersect(curtR,joints);
    num = length(roi);
    nodelabel(i) = num +1;
end


end

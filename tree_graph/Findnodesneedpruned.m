function [spls4,adj4,newPid,newCid,link,Radius2] = Findnodesneedpruned(spls,adj,Seg,Radius,Branch,BranchL,Pid,rootid)
%% ------------------------------------------------------------------------
delateb = [];
for i = 2:3
    % There are many fine branches in the canopy area, 
    % The main task is to prune the artificial fine branches of the first few important levels
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
    
    th = mean(curlratio) + std(curlratio);
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

%% Delete the skeleton points and their connection conditions
spls2 = spls;
spls2(delateb,:) = NaN;
adj2 = adj;
adj2(delateb,:) = 0; 
adj2(:,delateb) = 0;
% Remove the redundant broken branches through the connectivity of the graph
[newbin,~] = conncomp(digraph(adj2),'Type','weak','OutputForm','cell');
for i = 1:length(newbin)
    roi = intersect(newbin{1,i},rootid);
    if ~isempty(roi)
        newMainSkelton = newbin{1,i};
    end
end
all = [1:size(spls2,1)];
gap = setdiff(all,newMainSkelton);
newdelateb = setdiff(gap,delateb);
spls2(newdelateb,:) = NaN;
adj2(newdelateb,:) = 0; 
adj2(:,newdelateb) = 0;

newSeg = Simplifyskeleton(spls2,spls,newdelateb,Seg);

% Reorder the skeleton points 
% and remove the eliminated skeleton points from the sequence
[spls3,adj3,newSeg,Radius1,~] = Reconstructedtopology(spls2,adj2,newSeg,Radius);
[spls4,Radius2] = smoothspls2(spls3,adj3,Radius1,newSeg);

[~,rootid]=min(spls4(:,3)); 
[newCid,newPid ] = define_topologicalrelations(spls4,adj3,rootid);
link = zeros(0,2);
for i = 1:length(newPid)
    if i ~= rootid
        link(end+1,:) = [newPid(i) i];
    end
end

adj4 = zeros(size(adj3,1),size(adj3,1));
for i = 1:size(adj3,1)
    for j = 1:size(adj3,2)
        if adj3(i,j) ~= 0
            % Replace the weights with Euclidean distances to calculate the branch lengths
            adj4(i,j) = pdist2(spls4(i,:),spls4(j,:),'euclidean');
        end
    end
end

end
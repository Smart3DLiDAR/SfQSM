function [curbranch,Candidate] = findbranch(seedid,Cid,pli,Candidate)
curbranch = zeros(0,1); %储存当前访问的点
for i = 1:inf
    %id = knnsearch(spls,seedpoint);% 找到当前点在当前枝干骨架点中的位置
    id = seedid;
    %[~,col] = find(adj(id,:)==1);% 找到id在原骨架中的连接情况搜索下一个需要调整的点
    col = Cid(id,:);
    col(:,all(col==0,1))=[];%找到id在原骨架中下一个可能需要调整的点
    newroi = intersect(curbranch,id);
    if ~isempty(newroi)
        break
    end
    curbranch(end+1,1) = id;%储存当前处理的点索引         
    
    if length(col)==1
        seedid = col;
    elseif length(col)>1
        add = [];
        for j = 1:length(col)
            [row,~] = find(pli(:,1) == col(j));
            add = [add,row];
        end
        [~,maxid] = max(pli(add,2)); %比较连接点各自在路径中出现的频数 出现次数多的为主要的分支
        seedid = col(maxid(1));
        otherid = setdiff([1:length(col)],col(maxid(1)));
        for k = 1:length(otherid)
            Candidate(col(otherid(k))) = 1;%遇到分叉处的其他分支点储存作为新的起始点的候选点
        end
    else
        break
    end
end


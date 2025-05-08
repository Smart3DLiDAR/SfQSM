function [Cid,Pid ] = DefineTopologicalrelations(spls,adj,rootid)

joints = zeros(0,1);%存储分叉点
for i=1:size(spls,1)
    adj(i,i) = 0; %点和点自身的连接设置为0
    links = find(adj(i,:)==1 );% 找到和i连接的点
    if length(links) >2 %分叉点
        joints(end+1,1)=i;
    end
end

%% --- 根据具有连接关系的获取每个点的父节点和子节点------------------------------------
Num = zeros(length(joints),1);
% 计算每个分叉点的子节点个数 判断出最大子节点数
for i = 1:size(joints,1)
    [~,col] = find(adj(joints(i),:)==1);
    Num(i,1) = length(col) -1;
end
N = max(Num);%最大子节点数

Cid = zeros(size(spls,1),N);%储存子节点
beginid = rootid;
curvisit = zeros(0,1);
list = zeros(size(spls,1),1);
otherb = zeros(0,1);
for i = 1:inf
    [row,~] = find(list(:,1)==0);
    if ~isempty(row)
        id = beginid;
        [~,col] = find(adj(id,:)==1);
        curvisit(end+1,1) = id;
        if length(curvisit) ==1
            roiid = [];
        else
            [alreadyvisit,~] = find(list(:,1)== -1);
            [~,roiid,~] = intersect(col,alreadyvisit);
        end
        col(roiid) = [];% 将上一个访问过的连接点删掉
        num = length(col);
        if num == 1
            Cid(id,1:num) = col;
            beginid = col;
        elseif num>1
            Cid(id,1:num) = col;
            beginid = col(1);
            otherb(end+1:end+num-1,1) = col(2:end);%储存其他的子节点
        else%若查询到末端枝干点后无子节点 从之前的未查询的枝干子节点重新开始
            [row,~] = find(list(otherb,1)==0);
            if ~isempty(row)
                beginid = otherb(row(1));
            else
                break
            end
        end
        list(id,1) = -1; %当前访问的骨架点标签改为-1
    else
        break
    end
end

%获取每个点的父节点（根节点无父节点）
Pid = zeros(size(spls,1),1);%储存父节点
for i = 1:size(spls,1)
    curid = i;
    for j  = 1 :size(Cid,1)
        [~,col] = find(Cid(j,:)== curid);
        if ~isempty(col)
            Pid(i,1) = j; 
        end
    end
end

end
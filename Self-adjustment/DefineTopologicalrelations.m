function [Cid,Pid ] = DefineTopologicalrelations(spls,adj,rootid)

joints = zeros(0,1);
for i=1:size(spls,1)
    adj(i,i) = 0;
    links = find(adj(i,:)==1 );
    if length(links) >2 
        joints(end+1,1)=i;
    end
end


Num = zeros(length(joints),1);
for i = 1:size(joints,1)
    [~,col] = find(adj(joints(i),:)==1);
    Num(i,1) = length(col) -1;
end
N = max(Num);

Cid = zeros(size(spls,1),N);
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
        col(roiid) = [];
        num = length(col);
        if num == 1
            Cid(id,1:num) = col;
            beginid = col;
        elseif num>1
            Cid(id,1:num) = col;
            beginid = col(1);
            otherb(end+1:end+num-1,1) = col(2:end);
        else
            [row,~] = find(list(otherb,1)==0);
            if ~isempty(row)
                beginid = otherb(row(1));
            else
                break
            end
        end
        list(id,1) = -1; 
    else
        break
    end
end


Pid = zeros(size(spls,1),1);
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

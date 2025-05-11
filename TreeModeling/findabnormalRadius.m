function [normallabel,link,Radius,TH] = findabnormalRadius(Radius,finnalBranch,rootid,newCid,newPid,subbranchlength,link)
%%
normallabel = zeros(length(Radius),1);

for i = 1:length(finnalBranch)
    for j = 1:length(finnalBranch{i,1})
        curbranch = finnalBranch{i,1}{j,1};
        for k = 1:length(curbranch)
            id = curbranch(k);
            cursubbranchL = subbranchlength(id);
            
            if id ~= rootid
                pid = newPid(id,1);%
                cursubbranchLP = subbranchlength(pid);
                % compared the sub-branch length
                minR = Radius(pid)*((cursubbranchL/cursubbranchLP)^2);
                if Radius(id)<= Radius(pid) && Radius(id) >= minR
                    normallabel(id,:)=1;% 1：normal；0：abnormal
                end
            else
                normallabel(id,1) = 1;
            end
        end
    end
end

for i = 1:size(link,1)
    id1 = link(i,2);
    cursubbranchL = subbranchlength(id1);
    
    id2 = link(i,1);
    cursubbranchLP = subbranchlength(id2);
    
    son = newCid(id2,:);
    son(son==0) = [];

    th1 = log(Radius(id1)/Radius(id2));
    th2 = log(cursubbranchL/cursubbranchLP);

    % If the current point is the only child of the parent node
    % compared the sub-branch length
    if th1~=0 && th2 ~=0
        th = th1/th2;
    else
        th = 1.5;
    end
    
    if normallabel(id1)==1 && Radius(id2) >= Radius(id1)
        link(i,3:6) = [Radius(id2) Radius(id1) 1 th];
    else
        link(i,3:6) = [Radius(id2) Radius(id1) 0 th];
    end
end

TH = cell(length(finnalBranch),1);
for i = 1:length(finnalBranch)
    tempth = {};
    for j = 1:length(finnalBranch{i,1})
        curbranch = finnalBranch{i,1}{j,1};
        temp = [];
        for k = 1:length(curbranch)
            id = curbranch(k);%当前子节点
            if id ~= rootid
                [row,~] = find(link(:,2)==id);
                temp(end+1,:) = link(row,5:6);
            end
        end
        tempth{j,1} = temp;
    end
    TH{i,1} = tempth;
end

end

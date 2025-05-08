function cyl = branchmodeling(spls,finnalBranch,newRadius1,fixedR,branch2,newPid)
%% Build individual tree model
link2 = [];
for i = 1:length(finnalBranch)
    for j = 1:length(finnalBranch{i,1})
        curbranch = finnalBranch{i,1}{j,1};
        curRadius = fixedR{i,1}{j,1};
        if i ==1
            cur1 = curbranch(1:end-1);
            cur1Radius = curRadius(1:end-1);
            cur2 = curbranch(2:end);
            cur2Radius = curRadius(2:end);
            temp = [cur1,cur2,cur1Radius,cur2Radius];
            link2 = [link2;temp];
        else
            ppid = newPid(curbranch(1));
            if length(curbranch)>1
                cur1 = [ppid;curbranch(1:end-1)];
                cur1Radius = [newRadius1(ppid);curRadius(1:end-1)];
                cur2 = curbranch;
                cur2Radius = curRadius;
                temp = [cur1,cur2,cur1Radius,cur2Radius];
                link2 = [link2;temp];
            else
                cur1 = ppid;
                cur2 = curbranch;
                temp = [cur1,cur2,newRadius1(cur1),curRadius];
                link2 = [link2;temp];
            end
        end
    end
end

for i= 1:size(link2,1)
    idex1 = link2(i,1);
    for j = 1:length(branch2)
        tmp = ismember(branch2{j,1},idex1);
        if sum(tmp) ~= 0
            label1 = j;
        end
    end
    clear j 
    
    idex2 = link2(i,1);
    for j = 1:length(branch2)
        tmp = ismember(branch2{j,1},idex2);
        if sum(tmp) ~= 0
            label2 = j;
        end
    end
    
    branchid = max([label1 label2]);
    link2(i,5) = branchid;
end
%%
radius = zeros(size(link2,1),1);
long = zeros(size(link2,1),1);
start = zeros(size(link2,1),3);
Axe = zeros(size(link2,1),3);
bran = zeros(size(link2,1),1);
BranchOrder = zeros(size(link2,1),1);
for i = 1:length(link2)
    radius(i) = mean([link2(i,3) link2(i,4)]);%*1.1
    long(i) = pdist2(spls(link2(i,1),1:3),spls(link2(i,2),1:3));
    start(i,:) = spls(link2(i,1),1:3);
    Axe(i,:) = (spls(link2(i,2),1:3)-spls(link2(i,1),1:3))./norm(spls(link2(i,2),1:3)-spls(link2(i,1),1:3));
    bran(i) = link2(i,5);
    [~,~,~,id1] = searchlocationinbranch(link2(i,1),finnalBranch);
    [~,~,~,id2] = searchlocationinbranch(link2(i,2),finnalBranch);
    BranchOrder(i) = max([id1(1) id2(1)])-1;
end

cyl.radius = radius;
cyl.length = long;
cyl.start = start;
cyl.axis = Axe;
cyl.branch = bran;
cyl.BranchOrder = BranchOrder;
end
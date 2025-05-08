function Branch = ReOdertoSearchBranch(P,joints,MainSkeleton)

for i=1:size(P.spls,1)
   newskl(i,:)=[P.spls(i,1:3),i];
end
leavesklp = newskl;



%% Setp1:search the first branch
% the first seach to get the second point
for i = 1:size(MainSkeleton)
    [sortIndex,~] = find(leavesklp(:,4)==MainSkeleton(i,1));
    curPoint = 
end

[sortIndex,~] = find(newskl(:,4)==zIndex);
curPoint = newskl(sortIndex,:);
tree = KDTreeSearcher(newskl(:,1:3));
[kIndex,~] = knnsearch(tree,curPoint(:,1:3),'k',2);
%kIndex = cell2mat(kIndex);
kIndex = kIndex(2:end);
kPoint = newskl(kIndex,:);%No.2
first(2,:) = [kPoint(:,1:4),1,1,2];%save No.2
% oriIndex = zIndex; % get the original Index of the first point
leavesklp = newskl;
leavesklp(sortIndex,:) = [];%delate No.1
clear sortIndex

% Loop search
for i=2:size(leavesklp,1)
   [sortIndex,~] = find(leavesklp(:,4)==first(i,4));
   curPoint = leavesklp(sortIndex,:);
   tree = KDTreeSearcher(leavesklp(:,1:3));
   [kIndex,~] = knnsearch(tree,curPoint(:,1:3),'k',2);
%  kIndex = cell2mat(kIndex);
   kIndex = kIndex(2:end);
   kPoint = leavesklp(kIndex,:);
   first(i+1,:) = [kPoint(:,1:4),1,1,i+1];% save
   leavesklp(sortIndex,:) = [];%delate 
   for j = 1:size(joints,1)
       if (first(end,4) == joints(j,4)) > 0
           break
       end
   end
   
   % If there is no search point within the search radius or the root node
   % reached means the branch ended
   tree = KDTreeSearcher(leavesklp(:,1:3));
   [kIndex,~] = rangesearch(tree,kPoint(:,1:3),0.045);
   kIndex = cell2mat(kIndex);
   kIndex = kIndex(2:end);
        
   m = joints(j,4);
   if (first(end,4) == m)>0 || isempty (kIndex)
      break
   end
end

%% Setp2:search the second branch
% Distinguish the different branches from the first joint
[sortIndex,~] = find(leavesklp(:,4)==first(end,4));
curPoint = leavesklp(sortIndex,:);
tree = KDTreeSearcher(leavesklp(:,1:3));
[kIndex,~] = rangesearch(tree,curPoint(:,1:3),0.045);
kIndex = cell2mat(kIndex);
kIndex = kIndex(2:end);
Num1 = length(kIndex); %the son number of the first joint(main branch)
second = cell(Num1,1);
for i = 1:Num1
    kPoint(i,:) = leavesklp(kIndex(1,i),:);
    second{i,1}(1,:) = [curPoint(:,1:4),2,i,1];
    second{i,1}(2,:) = [kPoint(i,:),2,i,2];
end
leavesklp(sortIndex,:) = [];
% Loop search the first son branch(No.1 joint)
k = 2;
for i = 1:Num1
    [sortIndex,~] = find(leavesklp(:,4) == second{i,1}(end,4));
    curPoint = leavesklp(sortIndex,1:3);
    curvector = curPoint - second{i,1}(end-1,1:3);
    [kIndex,~] = knnsearch(leavesklp(:,1:3),curPoint(:,1:3),'k',k+1);
    kIndex = kIndex(2:end);
    Num2 = length(kIndex);
    next = [];
    vector = [];
    theta = [];
    for ii = 1:Num2
        next(ii,:) = leavesklp(kIndex(ii),1:3);
        vector(ii,:) = (next(ii,:) - curPoint);
        theta(ii,:) = acos((vector(ii,:)*curvector')/(norm(vector(ii,:))*norm(curvector)));
    end
    [~,Ind1] = min(theta);
    second{i,1}(end+1,:) = [leavesklp(kIndex(Ind1),1:4),2,i,3];
end
leavesklp = delate(leavesklp,second);
clear i ii Ind1 next vector theta Num2
% Loop search
[leavesklp,second] = Loopsearch_branch(leavesklp,second,joints,Num1,2);

%% Setp3:search the third branch
m = 1;
for i = 1:Num1
    cutend = second{i,1}(end,4);
    [sortIndex,~] = find(leavesklp(:,4) == cutend );
    curPoint = leavesklp(sortIndex,:);
    for j = 1:size(joints,1)
        if cutend == joints(j,4)
            key(i,:) = curPoint;
            [kIndex,~] = rangesearch(leavesklp(:,1:3),curPoint(:,1:3),0.045);
            kIndex = cell2mat(kIndex);
            kIndex = kIndex(2:end);
            Num2(1,m) = length(kIndex);
            
            TMP{m} = kIndex;
        end
    end
    kIndex = cat(2,TMP{:});
    m = m+1;
end
clear i j m 
Num2 = cumsum(Num2);
Num3 = Num2(end);
%third = cell(Num3,1);
third{1,1}(1,:) = [key(1,1:4),3,1,1];
for i =2:Num3
    if i <= Num2(1)
        third{end+1,1}(1,:) = [key(1,1:4),3,i,1];
    else
        third{end+1,1}(1,:) = [key(2,1:4),3,i,1];
    end 
end

for i =1:Num3
    third{i,1}(2,:) = [leavesklp(kIndex(i),1:4),3,i,2];
end

% Loop search the first son branch(No.1 joint)
k = 2;
for i = 1:Num3
    [sortIndex,~] = find(leavesklp(:,4) == third{i,1}(end,4));
    curPoint = leavesklp(sortIndex,1:3);
    curvector = curPoint - third{i,1}(end-1,1:3);
    [kIndex,~] = knnsearch(leavesklp(:,1:3),curPoint(:,1:3),'k',k+1);
    kIndex = kIndex(2:end);
    Num4 = length(kIndex);
    next = [];
    vector = [];
    theta = [];
    for ii = 1:Num4
        next(ii,:) = leavesklp(kIndex(ii),1:3);
        vector(ii,:) = next(ii,:) - curPoint;
        theta(ii,:) = acos((vector(ii,:)*curvector')/(norm(vector(ii,:))*norm(curvector)));
    end
    [~,Ind1] = min(theta);
    third{i,1}(end+1,:) = [leavesklp(kIndex(Ind1),1:4),3,i,3];
end
% Loop search
[leavesklp,third] = Loopsearch_branch(leavesklp,third,joints,Num3,3);

%% Setp4:search the forth branch
ii = 1;
for i = 1:Num2
    for j = 1:size(joints,1)
        if third{i,1}(end,4) == joints(j,4)
            keyp(ii) = third{i,1}(end,4);
            ii = ii+1;
        end
    end
end
clear ii i j

[sortIndex,~] = find(leavesklp(:,4) == keyp(1));
curPoint = leavesklp(sortIndex,:);
tree = KDTreeSearcher(leavesklp(:,1:3));
[kIndex,~] = rangesearch(tree,curPoint(:,1:3),0.035);
kIndex = cell2mat(kIndex);
kIndex = kIndex(2:end);
Num3 = length(kIndex); %the son number of the second joints
forth = cell(Num3,1);
for m = 1:Num3
    kPoint(m,:) = leavesklp(kIndex(1,m),:);
    forth{m,1}(1,:) = [curPoint(:,1:4),4,m,1];
    forth{m,1}(2,:) = [kPoint(m,:),4,m,2];
end
leavesklp(sortIndex ,:) = [];%delate
clear j m sortIndex kIndex kPoint 

[sortIndex,~] = find(leavesklp(:,4) == keyp(2));
curPoint = leavesklp(sortIndex,:);
tree = KDTreeSearcher(leavesklp(:,1:3));
[kIndex,~] = rangesearch(tree,curPoint(:,1:3),0.047);
kIndex = cell2mat(kIndex);
kIndex = kIndex(2:end);
Num4 = length(kIndex); %the son number of the second joints
%forth = cell(Num4,1);
for m = 1:Num4
    kPoint(m,:) = leavesklp(kIndex(1,m),:);
    forth{Num3+m,1}(1,:) = [curPoint(:,1:4),4,Num3+m,1];
    forth{Num3+m,1}(2,:) = [kPoint(m,:),4,Num3+m,2];
end
leavesklp(sortIndex ,:) = [];%delate
clear j m sortIndex kIndex kPoint 

[sortIndex,~] = find(leavesklp(:,4) == keyp(3));
curPoint = leavesklp(sortIndex,:);
tree = KDTreeSearcher(leavesklp(:,1:3));
[kIndex,~] = rangesearch(tree,curPoint(:,1:3),0.045);
kIndex = cell2mat(kIndex);
kIndex = kIndex(2:end);
Num5 = length(kIndex); %the son number of the second joints
%forth = cell(Num5,1);
for m = 1:Num5
    kPoint(m,:) = leavesklp(kIndex(1,m),:);
    forth{Num3+Num4+m,1}(1,:) = [curPoint(:,1:4),4,Num3+Num4+m,1];
    forth{Num3+Num4+m,1}(2,:) = [kPoint(m,:),4,Num3+Num4+m,2];
end
leavesklp(sortIndex ,:) = [];%delate
clear j m sortIndex kIndex kPoint Num3 Num4 Num5
Num3 = size(forth,1);
% Loop search 
for i = 2:Num3
    [leavesklp,forth{i,1}] = search(leavesklp,forth{i,1},joints,4,i);
end
[leavesklp,forth{1,1}] = search(leavesklp,forth{1,1},joints,4,1);

%% Setp5:search the fifth branch
ii = 1;
for i = 1:Num3
    for j = 1:size(joints,1)
        if forth{i,1}(end,4) == joints(j,4)
            keyp(ii) = forth{i,1}(end,4);
            ii = ii+1;
        end
    end
end
clear ii i j

[sortIndex,~] = find(leavesklp(:,4) == keyp(1));
curPoint = leavesklp(sortIndex,:);
tree = KDTreeSearcher(leavesklp(:,1:3));
[kIndex,~] = rangesearch(tree,curPoint(:,1:3),0.045);
kIndex = cell2mat(kIndex);
kIndex = kIndex(2:end);
Num4 = length(kIndex); %the son number of the second joints
fifth = cell(Num4,1);
for m = 1:Num4
    kPoint(m,:) = leavesklp(kIndex(1,m),:);
    fifth{m,1}(1,:) = [curPoint(:,1:4),5,m,1];
    fifth{m,1}(2,:) = [kPoint(m,:),5,m,2];
end
leavesklp(sortIndex ,:) = [];%delate
clear j m sortIndex kIndex kPoint curPoint

[sortIndex,~] = find(leavesklp(:,4) == keyp(2));
curPoint = leavesklp(sortIndex,:);
tree = KDTreeSearcher(leavesklp(:,1:3));
[kIndex,~] = rangesearch(tree,curPoint(:,1:3),0.045);
kIndex = cell2mat(kIndex);
kIndex = kIndex(2:end);
Num5 = length(kIndex); %the son number of the second joints
%forth = cell(Num4,1);
for m = 1:Num5
    kPoint(m,:) = leavesklp(kIndex(1,m),:);
    fifth{Num4+m,1}(1,:) = [curPoint(:,1:4),5,Num4+m,1];
    fifth{Num4+m,1}(2,:) = [kPoint(m,:),5,Num4+m,2];
end
leavesklp(sortIndex ,:) = [];%delate
clear j m sortIndex kIndex kPoint curPoint

[sortIndex,~] = find(leavesklp(:,4) == keyp(3));
curPoint = leavesklp(sortIndex,:);
tree = KDTreeSearcher(leavesklp(:,1:3));
[kIndex,~] = rangesearch(tree,curPoint(:,1:3),0.045);
kIndex = cell2mat(kIndex);
kIndex = kIndex(2:end);
Num6 = length(kIndex); %the son number of the second joints
%forth = cell(Num5,1);
for m = 1:Num6
    kPoint(m,:) = leavesklp(kIndex(1,m),:);
    fifth{Num4+Num5+m,1}(1,:) = [curPoint(:,1:4),5,Num4+Num5+m,1];
    fifth{Num4+Num5+m,1}(2,:) = [kPoint(m,:),5,Num4+Num5+m,2];
end
leavesklp(sortIndex ,:) = [];%delate
clear  m sortIndex kIndex kPoint curPoint

[sortIndex,~] = find(leavesklp(:,4) == keyp(4));
curPoint = leavesklp(sortIndex,:);
tree = KDTreeSearcher(leavesklp(:,1:3));
[kIndex,~] = rangesearch(tree,curPoint(:,1:3),0.045);
kIndex = cell2mat(kIndex);
kIndex = kIndex(2:end);
Num7 = length(kIndex); %the son number of the second joints
%forth = cell(Num5,1);
for m = 1:Num7
    kPoint(m,:) = leavesklp(kIndex(1,m),:);
    fifth{Num4+Num5+Num6+m,1}(1,:) = [curPoint(:,1:4),5,Num4+Num5+Num6+m,1];
    fifth{Num4+Num5+Num6+m,1}(2,:) = [kPoint(m,:),5,Num4+Num5+Num6+m,2];
end
leavesklp(sortIndex ,:) = [];%delate
clear j m sortIndex kIndex kPoint Num4 Num5 Num6 Num7
Num4 = size(fifth,1);
%[leavesklp,fifth] = Loopsearch_branch(leavesklp,fifth,joints,Num4,5);
for n = 2:Num4
    for j = 1:size(joints,1)
        [kIndex,~] = rangesearch(leavesklp(:,1:3),kPoint(:,1:3),0.05);
        kIndex = cell2mat(kIndex);
        kIndex = kIndex(2:end);
        if fifth{n,1}(end,4) ~= joints(:,4) || isempty (kIndex)
            [leavesklp,fifth{n,1}] = search(leavesklp,fifth{n,1},joints,5,n);
        end
    end
end
clear n j kIndex
[leavesklp,fifth{1,1}] = search(leavesklp,fifth{1,1},joints,5,1);


%define output
Branch.first = first;
Branch.second = second;
Branch.third = third;
Branch.forth = forth;
Branch.fifth = fifth;

end
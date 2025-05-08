function [spls1,Radius] = smoothspls2(spls,adj,Radius,Seg)
%%
joints = zeros(0,1);
roots = zeros(0,1);
for i=1:size(spls,1)
    links = find( adj(i,:)==1 );
    if length(links) >2 
        joints(end+1,1)=i;
    elseif length(links)==1
        roots(end+1,1) = i;
    end
end
[~,rootid] = min(spls(:,3));
Gp1 = graph(adj);
[tR,~] = shortestpathtree(Gp1,rootid,[1:size(spls,1)],'OutputForm','cell');

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
    %end
end

spls1 = spls;
for i = 1:length(branch)
    curbranch = branch{i,1};
    if length(curbranch)>=3
        points = spls(curbranch,:);
        [x1, y1] = bezir_n(points(:,1:2)', size(points,1));%8
        for j = 1:length(curbranch)
            spls1(curbranch(j),1:2) = [x1(j) y1(j)];
            temp = pdist2(Seg{curbranch(j),1},spls1(curbranch(j),:));
            Radius(curbranch(j)) = mean(temp(:));
        end
    end
end
 
branch1 = {};%´¢´æbranch
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
            branch1{end+1,1} = curbranch;
        end
        spls1(roots(i),3) = max(Seg{roots(i),1}(:,3));
    end
end

for i = 1:length(branch1)
    curbranch = branch1{i,1};
    if length(curbranch)>=3
        x = spls1(curbranch,1);
        y = spls1(curbranch,2);
        z = spls1(curbranch,3);
            
        interx = smooth(z,x,0.6,'rloess');
        intery = smooth(z,y,0.6,'rloess');
        for k = 1:length(curbranch)
           spls1(curbranch(k),1:2) = [interx(k) intery(k)];
           tempdis = pdist2(Seg{curbranch(k),1},spls1(curbranch(k),:));
           Radius(curbranch(k)) = mean(tempdis(:));
        end
    end
end




end

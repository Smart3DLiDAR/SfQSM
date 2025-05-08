function Cylinders = fitCylinder(skl,pts,spls,corresp,graph,Pid,Cid,Path1,joints)
%% wood2,leaf2,data
z = spls(:,3);
zmin = min(z);
zmax = max(z);
treeH = zmax - zmin;

points =  spls(:,1:3);
m = 1;
seg = [];
radius =[];
for i = 1:size(spls,1) 
    seg{i,1} = pts(corresp==i,:);
    if size(seg{i,1},1) == 1
        radius(m,:) = pdist2(spls(i,1:3), seg{i,1});
        points(i,4) = radius(m,:);
    else
        [~,~,radius(m,:)] = circlefit(seg{i,1});
        points(i,4) = radius(m,:);
        m=m+1;
    end
end
clear i m

%Radius anomaly detection
[radius,Pr] = RadiusAnomalydetection(skl,spls,Pid,Cid,radius,joints,Path1);
points(:,4) = radius;

m = 1;
Cylinders = {};
cly=[];
for i = 1:size(graph,1) 
    index1 = graph(i,1);
    cly.Parameters(m,1:3) =  points(index1,1:3);%终点（顶）
    curseg1 = pts(corresp==index1,:);
    [sort1,~] = find(curseg1(:,3) <= cly.Parameters(m,3));
    seg1 = curseg1(sort1,:);
    
    index2 = graph(i,2);
    cly.Parameters(m,4:6) =  points(index2,1:3);%起点（底）
    curseg2 = pts(corresp==index2,:);
    [sort2,~] = find(curseg2(:,3) >= cly.Parameters(m,6));
    seg2 = curseg2(sort2,:);
    cly.fitpoints = cat(1,seg1,seg2);
    
    cly.Parameters(m,7) =  (points(index1,4) + points(index2,4))/2;
    cly.Radius(m,:) = cly.Parameters(m,7);
    cly.Center(m,1:3) = (cly.Parameters(m,1:3) + cly.Parameters(m,4:6))/2;
    cly.Orientation(m,1:3) = cly.Parameters(m,1:3) - cly.Parameters(m,4:6);
    cly.Orientation(m,1:3) = cly.Orientation(m,1:3)./ norm(cly.Orientation(m,1:3));
    cly.Height(m,:) = norm(cly.Parameters(m,1:3) - cly.Parameters(m,4:6));
    Cylinders{i,1} = cly; 
    clear cly;
end
clear i

%figure(14);
%pcshow(wood2);

%thresholdR = min(radius);
for i=1:size(Cylinders,1)
    model = Cylinders{i,1};
%   pcwrite(result{i},['cycliner',num2str(i),'.pcd'])        %保存每个圆柱
    %hold on;
    figure(14);
    %subplot(1,3,1)
    %scatter3(branch.pts(:,1),branch.pts(:,2),branch.pts(:,3),30,'.');
    %scatter3(data(:,1),data(:,2),data(:,3),30,'.');
    %pcshow(pc);
    %axis off; axis equal;
    %title('Primitive tree')
    %subplot(1,3,2)
    %scatter3(wood2(:,1),wood2(:,2),wood2(:,3),30,[0.4471, 0.3216, 0.1647],'.');hold on;
    %scatter3(leaf2(:,1),leaf2(:,2),leaf2(:,3),30,[0.2667, 0.5686, 0.1961],'.');
    %axis off; axis equal;
    %title('After separation of branches')
    % if model.Center(:,3) < 0.75*treeH
    %subplot(1,3,3)
    color = rand(1,3);
    drawCylinder2(model.Parameters,'FaceColor',color);
    axis off; axis equal; camorbit(0,0,'camera');hold on
    title('Cylinderfit tree')
%    else
%        DrawCylinderV1(model.Center,model.Orientation,...
%        model.Height,thresholdR,[1,0,0]);
%        axis off; axis equal; camorbit(0,0,'camera');view(0,90);hold on
%    end
%   fprintf("第%d个圆柱，点数%d个\n",i,size(location,1));
end

end
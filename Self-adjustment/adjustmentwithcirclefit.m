function sklpoints = adjustmentwithcirclefit(spls,Seg,rootid,th1)
%% ----adjust Skeleton points with circlefit-----------------------------------
%----input --------------------------------------------------------------------
% spls : Skeleton points
% seg : The set of contraction points corresponding to every skeleton points
% rootid : The index of root node
% th1 : The Threshold of the effective domain covered the angular domain
%----output--------------------------------------------------------------------
% sklpoints : The information of tree after adjustment based on circlefit
%     sklpoints.spls = Skeleton points after adjusting
%     sklpoints.Radius : The radius of Skeleton points
%     sklpoints.label : The Reliability of skeleton point
%                       (1:The reliable skeleton point;0:The indeterminate skeleton point)
%     sklpoints.Seg : The set of contraction points corresponding to every skeleton points
%------------------------------------------------------------------------------
newspls = zeros(size(spls,1),3);
Radius = zeros(size(spls,1),1);
nodeslabel = zeros(size(newspls,1),1);
for i = 1:size(spls,1)
    curid = i;
    curspls = spls(i,:);
    cur = Seg{i,1};
    if isempty(cur)
        nodeslabel(i)=2;
        continue
    end
    centroid=mean(cur,1);
    % Local coordinate system transformation
    [Axis,newP,tempseg] = transferofaxes(curid,centroid,spls,Seg);
    if size(tempseg,1) >=20
        % the transformed polar coordinates are clustered using a density-based spatial clustering method
        % and calculated the effective domain by summing together the connected parts
        % Convert to polar coordinates (from Wangdi)
        [THETA,RHO] = pertransform(newP(:,1:2));
        [W, ~, ptsC] = dbscan2105([THETA,RHO]', 0.25, 3);%DBSCAN

        gapx = zeros(length(W),1);
        for j = 1:max(ptsC)
            curdbs = W{1,j};
            gapx(j) = max(curdbs(:,1)) - min(curdbs(:,1));
        end
        ratio = (sum(gapx))/(2*pi);

        Par = CircleFitByTaubin(newP(:,1:2));
        xc = Par(:,1);
        yc = Par(:,2);
        zc = mean(newP(:,3));
        R = Par(:,3);

        if ratio >= th1%2/3
            % If the effective domain covers more than two thirds of the angular domain
            % then the corresponding skeleton point is considered a reliable skeleton point
            nodeslabel(curid,1) = 1;
        end

        % Convert the fitted center coordinates to the original coordinate system
        Coords = [xc yc zc];
        rmatrix=[Axis(2,:);Axis(3,:);Axis(4,:)];
        xyz=rmatrix'*Coords';
        xyz=xyz'+Axis(1,:);

        % The distance between the adjusted skeleton point and the center point
        gap1 = sqrt(sum((centroid - xyz).^2));
        % The distance between the point before adjustment and the center point
        gap2 = sqrt(sum((centroid - curspls).^2));
        if gap1 <= gap2
            newspls(i,1:3)= xyz;%储存调整后的坐标
        else
            % Cancel the label of the reliable point and change it to 0
            newspls(i,1:3)= curspls;
            nodeslabel(curid,1) = 0;
        end
        Radius(i) = R;
    else
        % It is not convenient to fit the circle
        newspls(i,1:3)= centroid;
        dis = zeros(size(cur,1),1);
        for m = 1:size(cur,1)
            dis(m) = pdist2(cur(m,:),centroid);
        end
        Radius(i) = mean(dis);
    end

    if i == rootid
        newspls(i,3)= min(cur(:,3));
    end
end

sklpoints.spls = newspls;
sklpoints.Radius = Radius;
sklpoints.label = nodeslabel;
sklpoints.Seg = Seg;

end
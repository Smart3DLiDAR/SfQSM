function [averageDistance] = PointsAverageSpacing(data)
kdtree = KDTreeSearcher(data);
[~,distances]=knnsearch(kdtree,data,'K',2);
totalDistance = sum(sum(distances,'omitnan'));
averageDistance = totalDistance/size(data,1);
fprintf("ƽ������?%f\n",averageDistance)
end

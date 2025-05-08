function [averageDistance] = PointsAverageSpacing(data)
kdtree = KDTreeSearcher(data);
[~,distances]=knnsearch(kdtree,data,'K',2);
totalDistance = sum(sum(distances,'omitnan'));
averageDistance = totalDistance/size(data,1);
fprintf("平均点间距?%f\n",averageDistance)
end

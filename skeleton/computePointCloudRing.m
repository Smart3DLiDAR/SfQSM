function ring = computePointCloudRing(pts, k, index)
%计算点云中的一环邻域

%pts:n*3矩阵，用于计算1-环领域的坐标。
%k:kNN的k
%索引：kNN的索引

npts = size(pts,1);     %获取点数
ring = cell(npts,1);    %存储环

if nargin < 3 || isempty(index) %函数参数小于3且index为空
    kdtree = KDTreeSearcher(pts);
    index = knnsearch(kdtree,pts,'K',k);
end

for i = 1:npts  
    neighbor = pts(index(i,:),:);       %邻近点
    coefs = pca(neighbor, 'Economy', false);    %使用PCA算法计算局部的特征向量
    x = [neighbor * coefs(:, 1), neighbor * coefs(:, 2)]; %使用两个主方向构建平面，并将点数据投影到该平面
    tempx = unique(x,'rows');
    if size(tempx,1)<=3
        ring{i,:} = index(i,2:7);
        continue
    end

    TRI = delaunayn(x);         %构建局部三角网
    
    [row,~] = find(TRI == 1);
    temp = TRI(row,:);      %找寻与该点相关的点
    
    temp = sort(temp,2);    %每行按照升序进行排序
    temp = temp(:,2:end);
    if isempty(temp)
        ring{i,:} = index(i,2:7);
        continue
    end
%     
    x=temp(:);
    x=sort(x);      %对边界点进行排序
    d=diff([x;max(x)+1]);   %计算相邻元素之差
    count = diff(find([1;d]));
    %d2 = diff([min(x)-1;x]);
    %count2 = diff(find([d2;1]));
    
    y =[x(find(d)) count];
    n_sorted_index = size(y,1);
    start = find(count==1);     %找到连续的索引
    if ~isempty(start)
        want_to_find = y(start(1),1);
    else
        want_to_find = temp(1,1);
        n_sorted_index = n_sorted_index+1;
    end
    
    %对边界索引进行排序
    j = 0;    
    sorted_index = zeros(1,n_sorted_index);
    while j < n_sorted_index
        j = j+1;
        sorted_index(j) = want_to_find;
        [row,col] = find(temp == want_to_find);     %找到都有那个边与该点相连接
        if ~isempty(col)
            if col(1) == 1      %判断是否为第一列
                want_to_find = temp(row(1),2);
                temp(row(1),2) = -1;
            else
                want_to_find = temp(row(1),1);
                temp(row(1),1) = -1;
            end    
        end
    end
    
    neighbor_index = index(i,sorted_index);
    
    ring{i} = neighbor_index;
end

end


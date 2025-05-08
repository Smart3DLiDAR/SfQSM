function [samplePoint,corresp] = farthestSamplingBySphere(pts,radius)

%使用大小为RADIUS的球以最远的采样方式对点云进行降采样
%加速：我们只更新距离超过%90！
%pts：要采样的点
%嘶声：嘶声
%RADIUS：采样半径
%spls：采样点
%corresp：pts和spls之间的对应关系，|pts|*1的数组

%% **********参数设置************
SHOW_SAMPLING_PROGRESS = true;
SHOW_RESULTS = true;

%% **********采样过程************
if SHOW_SAMPLING_PROGRESS || SHOW_RESULTS
    close all;
    figure; movegui('northwest');set(gcf,'color','white');hold on;
    plot3( pts(:,1), pts(:,2), pts(:,3), '.r', 'markersize', 1);
    axis off; axis equal;set(gcf,'Renderer','OpenGL');
end

tic
kdtree = KDTreeSearcher(pts);
samplePoint = [];
corresp = zeros(length(pts), 1);
mindst = nan(length(pts), 1);

for k=1:length(pts)
    if corresp(k)~=0, continue, end
        
    %根据距离查询所有点
    mindst(k) = inf;
    
    % 初始化优先队列
    while ~all(corresp~=0)
        [~, maxIdx] = max( mindst );
        if mindst(maxIdx) == 0
            break
        end

        [nIdxs, nDsts] = rangesearch( kdtree,pts(maxIdx,:), radius );
        nIdxs = nIdxs{1};
        nDsts = nDsts{1};
        % 如果邻近点都被标记，则跳过
        if all( corresp(nIdxs) ~= 0 )
            mindst(maxIdx) = 0; 
            continue;
        end

        % 新建节点，更新距离
        samplePoint(end+1,:) = pts(maxIdx,:);
        for i=1:length(nIdxs)
            if mindst(nIdxs(i))>nDsts(i) || isnan(mindst(nIdxs(i)))
               mindst(nIdxs(i)) = nDsts(i);
               corresp(nIdxs(i)) = size(samplePoint,1);
            end
        end

        if SHOW_SAMPLING_PROGRESS == true
            figure(1); plot3( pts(maxIdx,1), pts(maxIdx,2), pts(maxIdx,3), '*g');
        end
    end
end
toc

end


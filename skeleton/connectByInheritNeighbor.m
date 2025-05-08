function A = connectByInheritNeighbor(pts, samplePoint, corresp, neighbor)
%连接连续的邻近点:通过继承它们对应的样本（pts）的邻居来构建下样本（spls）的连接矩阵。
%
%pts：采样点
%samplePoint：下采样点
%corresp：pts和spls之间的对应关系，|pts|*1的数组
%neighbor：pts，a |pts|*1cell，邻居按距离升序排序
%A：下采样的连接矩阵（spls）
%

%% 是否显示过程
SHOW_CONNECT_PROGRESS = false;
SHOW_RESULTS = false;

if SHOW_CONNECT_PROGRESS || SHOW_RESULTS
    close all;
    figure; movegui('northwest');set(gcf,'color','white');hold on;
    plot3( pts(:,1), pts(:,2), pts(:,3), '.r', 'markersize', 1);
    plot3( samplePoint(:,1), samplePoint(:,2), samplePoint(:,3), '.g', 'markersize', 20);
    axis off; axis equal;set(gcf,'Renderer','OpenGL');
end
if ~iscell(neighbor)
    tmp = cell(size(neighbor,1),1);
    for i = 1:size(neighbor,1)
        tmp{i}=neighbor(i,:);
    end
    neighbor = tmp;  clear tmp;
end

A = zeros( length(samplePoint), length(samplePoint) );
%遍历每一个点
for pIdx=1:length(pts)    
    ns =neighbor{pIdx};     %局部多边形
    pc = corresp(pIdx);     %匹配点
    if pc == 0
        warning('该点没有匹配点');continue;
    end
    
    %遍历局部点
    for nIdx=1:length(ns)
        nc = corresp(ns(nIdx));     %邻近点的匹配点
        if nc == 0      
            warning('该点没有匹配点');continue;
        end
        if nc~=pc       %判断该点的匹配点与邻近点的匹配点是否相同
            A(pc,nc) = A(pc,nc) + 1;        %累加
            A(nc,pc) = A(nc,pc) + 1;
            if SHOW_CONNECT_PROGRESS
                figure(1); idx = [pc, nc];
                line( samplePoint(idx,1),samplePoint(idx,2),samplePoint(idx,3), 'LineWidth', 2,'Color','b');
            end
            break;
        end
    end
end

%% 如果存在孤立点，连接它的最近点
isopts = zeros(1,0);
for i=1:size(A,1)
    A(i,i) = 1;     %将对角元素置为1
    if length( find( A(i,:)>0) ) == 1       %判断是否为孤立点
        isopts(1, end+1) = i;
        warning('there are isolate points: %d', i);
    end
end
if isempty(isopts)
    spls_kdtree = KDTreeSearcher(samplePoint);
    for i=isopts
        [neighs,distance] = knnsearch(spls_kdtree, samplePoint(i,:), 2);
        for j = neighs
            A(i,j) = 1;
            A(j,i) = 1;
        end
    end
end

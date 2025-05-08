function [samplePoint, A, corresp] = edgeCollapse(samplePoint, corresp, spls_adj, options)
%对采样点进行简化：折叠邻接矩阵spls_adj，直到所描述的图是一维的。
%
%pts：采样点
%samplePoint：下采样
%corresp：pts和spls之间的对应关系，|pts|*1的数组
%spls_adj：下采样的连接矩阵（spls）
%options.collapse_order：控制折叠的权重项,0仅适用于欧式距离；1表示度，然后表示距离。
%A：边缘折叠后下采样（spls）的新连接矩阵

options.null = 0;
collapse_order = options.collapse_order;
SHOW_COLLAPSE_PROGRESS = true;
SHOW_RESULTS = true;

if SHOW_COLLAPSE_PROGRESS || SHOW_RESULTS
    close all;
    figure; movegui('northwest');set(gcf,'color','white');hold on;
    plot3( samplePoint(:,1), samplePoint(:,2), samplePoint(:,3), '.r', 'markersize', 5);
    axis off; axis equal;set(gcf,'Renderer','OpenGL');
end

A = spls_adj;       %邻接矩阵
A(A>0) = 1;

degrees = ones(size(samplePoint,1),1);     
for i=1:length(samplePoint)
    ns = find(A(i,:)==1);
    degrees(i) = length(ns)-1;      %统计该点与几个点相连
end

tricount = 0;
skeds = [];     %记录端点的平均度以及距离
for i=1:length(samplePoint)
    ns = find(A(i,:)==1);   %查找连接点
    ns = ns( ns>i );        %每个三角形只一次(如果一条边属于两个三角形，它会出现两次!)
    lns = length(ns);       %相连点数
    for j=1:lns
        for k=j+1:lns
            if A(ns(j),ns(k)) == 1      %判断这两个点是否连接，如果连接，则组成三角形
                tricount = tricount+1;
                skeds(end+1,1:3) = [i,ns(j), 0.5*(lns+degrees(ns(j))) ];                                                        %#ok<AGROW>
                skeds(end,4) = euclidean_distance(samplePoint(i,:), samplePoint(ns(j),:) ); 
                skeds(end+1,1:3) = [ns(j),ns(k), 0.5*(degrees(ns(j)) +degrees(ns(k)) )];                                                %#ok<AGROW>
                skeds(end,4) = euclidean_distance( samplePoint(ns(j),:), samplePoint(ns(k),:) );
                skeds(end+1,1:3) = [i,ns(k), 0.5*(lns+degrees(ns(k))) ];                                                     %#ok<AGROW>
                skeds(end,4) = euclidean_distance( samplePoint(ns(k),:), samplePoint(i,:) );
            end
        end
    end
end
    
%% 移除边缘
while true
    if SHOW_COLLAPSE_PROGRESS
        heds = [];
        for i=1:size(A,1)
            for j=1:size(A,2)
                if( A(i,j)==1 )
                    idx = [i;j];
                    heds(end+1) = line( samplePoint(idx,1),samplePoint(idx,2),samplePoint(idx,3), 'LineWidth', 2, 'Color', 'b');
                end
            end
        end
        drawnow update
        pause(0.2);     %暂停0.2s
        delete(heds);
    end

    fprintf('抽取骨架图，剩余%d连接边\n', size(skeds,1));
    
    % 不存在三角形，停止计算
    if size(skeds,1) == 0, break, end
    
    % 以最小代价移除边缘，删除第二个顶点
    if collapse_order == 1  % 代价是degree + distance
        mind = min( skeds(:,3) );       %最小的平均度
        tmpIdx = find(skeds(:,3)==mind);
        tmpSkeds = skeds(tmpIdx,4);

        [~, idx] = min( tmpSkeds );     %找到最小的距离边
        edge = skeds(tmpIdx(idx),1:2);
        skeds(tmpIdx(idx),:)=[];        %将该边进行删除
    else % 代价是距离
        [~, idx] = min( skeds(:,4) );
        edge = skeds(idx,1:2);
        skeds(idx,:)=[];
    end
    fprintf( '删除边: %d, %d\n', edge(1),edge(2));

    % 更新位置
    samplePoint( edge(2),: ) = mean( samplePoint( edge,: ) );
    samplePoint( edge(1),: ) = NaN;
    % 更新邻接矩阵
    for k=1:size(A,1)
        if A(edge(1),k) == 1
            A(edge(2),k)=1; 
            A(k,edge(2))=1; 
        end
    end
    % 移除该行
    A(edge(1),:) = 0;
    A(:,edge(1)) = 0;
    % 更新匹配点
    corresp(corresp==edge(1) ) = edge(2);
    
    %%
    tmpIdx = skeds( skeds(:,1)==edge(2), 2);
    tmpIdx = [tmpIdx; skeds( skeds(:,2)==edge(2), 1)];
    
    [rows,cols] = find(skeds(:,1:2)==edge(1));
    toBeRemoved =  zeros(0,1);   
    for i = 1:length(rows)
        col = 1 + mod(cols(i),2);
        if ismember( skeds(rows(i), col), tmpIdx )%remove
            toBeRemoved(end+1) = rows(i);
        else
            skeds(rows(i), cols(i)) = edge(2);
        end
    end
    if ~isempty(toBeRemoved)
        skeds(toBeRemoved,: ) = [];
    end
    
    %% 2) 移除包含edge(2)的骨架边，使其不再是一个三角形
    % 之后添加一个新的三角形包含edge(2).
    % 找到所有包含edge(2)的三角形
    ns = find( A(edge(2),:)==1 );
    ns = ns( ns~=edge(2) );
    lns = length(ns);

    tmpEdges = zeros(0,2); % 包含edge(2)
    tmpEdges1 = []; % 不包含edge(2)
    for j=1:lns
        for k=j+1:lns
            if A(ns(j),ns(k)) == 1
                tmpEdges(end+1,:) = [edge(2),ns(j)];  
                tmpEdges(end+1,:) = [edge(2),ns(k)];
                tmpEdges1(end+1,:) = [ns(j),ns(k)];                
            end
        end
    end
    
    % 移除所有包含edge(2)的边
    [rows,cols] = find(skeds(:,1:2)==edge(2));
    toBeRemoved =  zeros(0,1);
    tobedel = []; 
    for j = 1:length(rows)
        col = 1 + mod(cols(j),2);
        tmp = find( tmpEdges(:,2) == skeds(rows(j),col) );
        if tmp % 移除边    
            for k = 1:length(tmp)
                if ~ismember (tmp(k), tobedel)
                    tobedel = [tobedel; tmp(k)];
                end
            end
        else
            toBeRemoved(end+1) = rows(j);
        end
    end
    if ~isempty(toBeRemoved)
        skeds( toBeRemoved,: ) = [];
    end
    if ~isempty(tobedel)
        tmpEdges(tobedel,:) = [];
    end
    
    % 添加新的三角形
    tmpEdges = [tmpEdges; tmpEdges1];
    for j = 1:size(tmpEdges, 1)
        tedge = tmpEdges(j,:);
        [rows,cols] = find(skeds(:,1:2)==tedge(1));
        bin = false;
        for k = 1:length(rows)
            col = 1 + mod(cols(k),2);
            if skeds(rows(k),col)==tedge(2) %这个边已经在skeds中存在
                bin = true;
                break;
            end
        end
        if ~bin % 在骨架中添加边
            ns = find(A(tedge(1),:)==1);
            degrees(tedge(1)) = (length(ns)-1)*0.5;
            ns = find(A(tedge(2),:)==1);
            degrees(tedge(2)) = (length(ns)-1)*0.5;
            skeds(end+1,1:2) = tedge;
            skeds(end,3) = 0.5*(degrees(tedge(1))+degrees(tedge(2)));
            skeds(end,4) = euclidean_distance(samplePoint(tedge(1),:), samplePoint(tedge(2),:) );             
        end
    end
    
    %% 3) 更新边的距离与自由度
    [rows,cols] = find(skeds(:,1:2)==edge(2)); 
    ns = find(A(edge(2),:)==1);
    degrees(edge(2)) = (length(ns)-1)*0.5;
    
    for j = 1:length(rows)
        col = 1 + mod(cols(j),2);   
        k = skeds(rows(j),col);
        ns = find(A(k,:)==1);
        degrees(k) = (length(ns)-1)*0.5;
        
        skeds(rows(j),3) = 0.5*(degrees(edge(2))+degrees(k));
        skeds(rows(j),4) = euclidean_distance(samplePoint(edge(2),:), samplePoint(k,:) ); 
    end
end
function dist = euclidean_distance(p1, p2)
v=p1-p2;
dist = sqrt(dot(v,v));





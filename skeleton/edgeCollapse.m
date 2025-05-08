function [samplePoint, A, corresp] = edgeCollapse(samplePoint, corresp, spls_adj, options)
%�Բ�������м򻯣��۵��ڽӾ���spls_adj��ֱ����������ͼ��һά�ġ�
%
%pts��������
%samplePoint���²���
%corresp��pts��spls֮��Ķ�Ӧ��ϵ��|pts|*1������
%spls_adj���²��������Ӿ���spls��
%options.collapse_order�������۵���Ȩ����,0��������ŷʽ���룻1��ʾ�ȣ�Ȼ���ʾ���롣
%A����Ե�۵����²�����spls���������Ӿ���

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

A = spls_adj;       %�ڽӾ���
A(A>0) = 1;

degrees = ones(size(samplePoint,1),1);     
for i=1:length(samplePoint)
    ns = find(A(i,:)==1);
    degrees(i) = length(ns)-1;      %ͳ�Ƹõ��뼸��������
end

tricount = 0;
skeds = [];     %��¼�˵��ƽ�����Լ�����
for i=1:length(samplePoint)
    ns = find(A(i,:)==1);   %�������ӵ�
    ns = ns( ns>i );        %ÿ��������ֻһ��(���һ�����������������Σ������������!)
    lns = length(ns);       %��������
    for j=1:lns
        for k=j+1:lns
            if A(ns(j),ns(k)) == 1      %�ж����������Ƿ����ӣ�������ӣ������������
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
    
%% �Ƴ���Ե
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
        pause(0.2);     %��ͣ0.2s
        delete(heds);
    end

    fprintf('��ȡ�Ǽ�ͼ��ʣ��%d���ӱ�\n', size(skeds,1));
    
    % �����������Σ�ֹͣ����
    if size(skeds,1) == 0, break, end
    
    % ����С�����Ƴ���Ե��ɾ���ڶ�������
    if collapse_order == 1  % ������degree + distance
        mind = min( skeds(:,3) );       %��С��ƽ����
        tmpIdx = find(skeds(:,3)==mind);
        tmpSkeds = skeds(tmpIdx,4);

        [~, idx] = min( tmpSkeds );     %�ҵ���С�ľ����
        edge = skeds(tmpIdx(idx),1:2);
        skeds(tmpIdx(idx),:)=[];        %���ñ߽���ɾ��
    else % �����Ǿ���
        [~, idx] = min( skeds(:,4) );
        edge = skeds(idx,1:2);
        skeds(idx,:)=[];
    end
    fprintf( 'ɾ����: %d, %d\n', edge(1),edge(2));

    % ����λ��
    samplePoint( edge(2),: ) = mean( samplePoint( edge,: ) );
    samplePoint( edge(1),: ) = NaN;
    % �����ڽӾ���
    for k=1:size(A,1)
        if A(edge(1),k) == 1
            A(edge(2),k)=1; 
            A(k,edge(2))=1; 
        end
    end
    % �Ƴ�����
    A(edge(1),:) = 0;
    A(:,edge(1)) = 0;
    % ����ƥ���
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
    
    %% 2) �Ƴ�����edge(2)�ĹǼܱߣ�ʹ�䲻����һ��������
    % ֮�����һ���µ������ΰ���edge(2).
    % �ҵ����а���edge(2)��������
    ns = find( A(edge(2),:)==1 );
    ns = ns( ns~=edge(2) );
    lns = length(ns);

    tmpEdges = zeros(0,2); % ����edge(2)
    tmpEdges1 = []; % ������edge(2)
    for j=1:lns
        for k=j+1:lns
            if A(ns(j),ns(k)) == 1
                tmpEdges(end+1,:) = [edge(2),ns(j)];  
                tmpEdges(end+1,:) = [edge(2),ns(k)];
                tmpEdges1(end+1,:) = [ns(j),ns(k)];                
            end
        end
    end
    
    % �Ƴ����а���edge(2)�ı�
    [rows,cols] = find(skeds(:,1:2)==edge(2));
    toBeRemoved =  zeros(0,1);
    tobedel = []; 
    for j = 1:length(rows)
        col = 1 + mod(cols(j),2);
        tmp = find( tmpEdges(:,2) == skeds(rows(j),col) );
        if tmp % �Ƴ���    
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
    
    % ����µ�������
    tmpEdges = [tmpEdges; tmpEdges1];
    for j = 1:size(tmpEdges, 1)
        tedge = tmpEdges(j,:);
        [rows,cols] = find(skeds(:,1:2)==tedge(1));
        bin = false;
        for k = 1:length(rows)
            col = 1 + mod(cols(k),2);
            if skeds(rows(k),col)==tedge(2) %������Ѿ���skeds�д���
                bin = true;
                break;
            end
        end
        if ~bin % �ڹǼ�����ӱ�
            ns = find(A(tedge(1),:)==1);
            degrees(tedge(1)) = (length(ns)-1)*0.5;
            ns = find(A(tedge(2),:)==1);
            degrees(tedge(2)) = (length(ns)-1)*0.5;
            skeds(end+1,1:2) = tedge;
            skeds(end,3) = 0.5*(degrees(tedge(1))+degrees(tedge(2)));
            skeds(end,4) = euclidean_distance(samplePoint(tedge(1),:), samplePoint(tedge(2),:) );             
        end
    end
    
    %% 3) ���±ߵľ��������ɶ�
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





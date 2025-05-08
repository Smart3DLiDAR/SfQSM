function [samplePoint,corresp] = farthestSamplingBySphere(pts,radius)

%ʹ�ô�СΪRADIUS��������Զ�Ĳ�����ʽ�Ե��ƽ��н�����
%���٣�����ֻ���¾��볬��%90��
%pts��Ҫ�����ĵ�
%˻����˻��
%RADIUS�������뾶
%spls��������
%corresp��pts��spls֮��Ķ�Ӧ��ϵ��|pts|*1������

%% **********��������************
SHOW_SAMPLING_PROGRESS = true;
SHOW_RESULTS = true;

%% **********��������************
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
        
    %���ݾ����ѯ���е�
    mindst(k) = inf;
    
    % ��ʼ�����ȶ���
    while ~all(corresp~=0)
        [~, maxIdx] = max( mindst );
        if mindst(maxIdx) == 0
            break
        end

        [nIdxs, nDsts] = rangesearch( kdtree,pts(maxIdx,:), radius );
        nIdxs = nIdxs{1};
        nDsts = nDsts{1};
        % ����ڽ��㶼����ǣ�������
        if all( corresp(nIdxs) ~= 0 )
            mindst(maxIdx) = 0; 
            continue;
        end

        % �½��ڵ㣬���¾���
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


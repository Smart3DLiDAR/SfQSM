function A = connectByInheritNeighbor(pts, samplePoint, corresp, neighbor)
%�����������ڽ���:ͨ���̳����Ƕ�Ӧ��������pts�����ھ���������������spls�������Ӿ���
%
%pts��������
%samplePoint���²�����
%corresp��pts��spls֮��Ķ�Ӧ��ϵ��|pts|*1������
%neighbor��pts��a |pts|*1cell���ھӰ�������������
%A���²��������Ӿ���spls��
%

%% �Ƿ���ʾ����
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
%����ÿһ����
for pIdx=1:length(pts)    
    ns =neighbor{pIdx};     %�ֲ������
    pc = corresp(pIdx);     %ƥ���
    if pc == 0
        warning('�õ�û��ƥ���');continue;
    end
    
    %�����ֲ���
    for nIdx=1:length(ns)
        nc = corresp(ns(nIdx));     %�ڽ����ƥ���
        if nc == 0      
            warning('�õ�û��ƥ���');continue;
        end
        if nc~=pc       %�жϸõ��ƥ������ڽ����ƥ����Ƿ���ͬ
            A(pc,nc) = A(pc,nc) + 1;        %�ۼ�
            A(nc,pc) = A(nc,pc) + 1;
            if SHOW_CONNECT_PROGRESS
                figure(1); idx = [pc, nc];
                line( samplePoint(idx,1),samplePoint(idx,2),samplePoint(idx,3), 'LineWidth', 2,'Color','b');
            end
            break;
        end
    end
end

%% ������ڹ����㣬�������������
isopts = zeros(1,0);
for i=1:size(A,1)
    A(i,i) = 1;     %���Խ�Ԫ����Ϊ1
    if length( find( A(i,:)>0) ) == 1       %�ж��Ƿ�Ϊ������
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

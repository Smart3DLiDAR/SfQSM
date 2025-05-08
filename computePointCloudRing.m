function ring = computePointCloudRing(pts, k, index)
%��������е�һ������

%pts:n*3�������ڼ���1-����������ꡣ
%k:kNN��k
%������kNN������

npts = size(pts,1);     %��ȡ����
ring = cell(npts,1);    %�洢��

if nargin < 3 || isempty(index) %��������С��3��indexΪ��
    kdtree = KDTreeSearcher(pts);
    index = knnsearch(kdtree,pts,'K',k);
end

for i = 1:npts  
    neighbor = pts(index(i,:),:);       %�ڽ���
    coefs = pca(neighbor, 'Economy', false);    %ʹ��PCA�㷨����ֲ�����������
    x = [neighbor * coefs(:, 1), neighbor * coefs(:, 2)]; %ʹ�����������򹹽�ƽ�棬����������ͶӰ����ƽ��
    tempx = unique(x,'rows');
    if size(tempx,1)<=3
        ring{i,:} = index(i,2:7);
        continue
    end

    TRI = delaunayn(x);         %�����ֲ�������
    
    [row,~] = find(TRI == 1);
    temp = TRI(row,:);      %��Ѱ��õ���صĵ�
    
    temp = sort(temp,2);    %ÿ�а��������������
    temp = temp(:,2:end);
    if isempty(temp)
        ring{i,:} = index(i,2:7);
        continue
    end
%     
    x=temp(:);
    x=sort(x);      %�Ա߽���������
    d=diff([x;max(x)+1]);   %��������Ԫ��֮��
    count = diff(find([1;d]));
    %d2 = diff([min(x)-1;x]);
    %count2 = diff(find([d2;1]));
    
    y =[x(find(d)) count];
    n_sorted_index = size(y,1);
    start = find(count==1);     %�ҵ�����������
    if ~isempty(start)
        want_to_find = y(start(1),1);
    else
        want_to_find = temp(1,1);
        n_sorted_index = n_sorted_index+1;
    end
    
    %�Ա߽�������������
    j = 0;    
    sorted_index = zeros(1,n_sorted_index);
    while j < n_sorted_index
        j = j+1;
        sorted_index(j) = want_to_find;
        [row,col] = find(temp == want_to_find);     %�ҵ������Ǹ�����õ�������
        if ~isempty(col)
            if col(1) == 1      %�ж��Ƿ�Ϊ��һ��
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


function [curbranch,Candidate] = findbranch(seedid,Cid,pli,Candidate)
curbranch = zeros(0,1); %���浱ǰ���ʵĵ�
for i = 1:inf
    %id = knnsearch(spls,seedpoint);% �ҵ���ǰ���ڵ�ǰ֦�ɹǼܵ��е�λ��
    id = seedid;
    %[~,col] = find(adj(id,:)==1);% �ҵ�id��ԭ�Ǽ��е��������������һ����Ҫ�����ĵ�
    col = Cid(id,:);
    col(:,all(col==0,1))=[];%�ҵ�id��ԭ�Ǽ�����һ��������Ҫ�����ĵ�
    newroi = intersect(curbranch,id);
    if ~isempty(newroi)
        break
    end
    curbranch(end+1,1) = id;%���浱ǰ����ĵ�����         
    
    if length(col)==1
        seedid = col;
    elseif length(col)>1
        add = [];
        for j = 1:length(col)
            [row,~] = find(pli(:,1) == col(j));
            add = [add,row];
        end
        [~,maxid] = max(pli(add,2)); %�Ƚ����ӵ������·���г��ֵ�Ƶ�� ���ִ������Ϊ��Ҫ�ķ�֧
        seedid = col(maxid(1));
        otherid = setdiff([1:length(col)],col(maxid(1)));
        for k = 1:length(otherid)
            Candidate(col(otherid(k))) = 1;%�����ֲ洦��������֧�㴢����Ϊ�µ���ʼ��ĺ�ѡ��
        end
    else
        break
    end
end


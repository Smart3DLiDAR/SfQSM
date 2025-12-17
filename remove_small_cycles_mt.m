function [spls, corresp, spls_adj, joints, segments] = remove_small_cycles_mt(pts,spls, corresp, spls_adj, joints,roots,segments, threshold, show_cycles)
% removing small cycles measured by topological length <= threshold
G = spls_adj;
cycles = findcycles(sparse(G));% can be speed up by just using joints, not every nodes.
kept_joints = zeros(0,1);
for i = 1:length(cycles)
    cycle = cycles{i};
    if length(cycle)<3
        continue
    end
    if show_cycles
        tmp = [cycle, cycle(1)];
        for j = 1:(length(tmp)-1)
            myedge3(spls(tmp(j),:), spls(tmp(j+1),:), 'Color',[1 0 0]);
        end
    end
 %%%%%% mt delete big cycles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if length(cycle)>threshold %删除长的回路
        ja=[];
        for j = 1:length(cycle) 
            if ismember(cycle(j), joints(:,1))% a joint
                id = cycle(j);
                ja=[id;ja]; %找到关键点索引
            end
        end
        dis=[];
        Va=spls(ja,:);%找到关键点xyz
        Dis= pdist(Va);%计算关键点距离
        Dis_square=squareform(Dis); 
        SumDis=sum(Dis_square);
        [~,idx]=max(SumDis);
        del_joint_id=ja(idx); 
        del_joint=spls(del_joint_id,:);%找到需要删除的关键点
        nei=G(del_joint_id,:);
        nei(del_joint_id)=0;
        indices=find(nei==1);%找到该关键点的连接的点索引
        va=spls(indices,:);
        num=size(va,1);%该关键点的连接的点数量
        if(num>3)
            continue;
        end
        %计算该关键点与连接点的方向向量
        v1=va(1,:)-del_joint;
        v2=va(2,:)-del_joint;
        v3=va(3,:)-del_joint;
        v1=v1./norm(v1);
        v2=v2./norm(v2);
        v3=v3./norm(v3);
        theta12=abs(dot(v1,v2));
        theta13=abs(dot(v1,v3));
        theta23=abs(dot(v2,v3));
        
        if(theta12>=theta13&&theta12>=theta23)
          %% delete edge indices(3)-del_joint_id
          edgeV=indices(3);
          spls_adj(edgeV,del_joint_id) = 0;
          spls_adj(del_joint_id,edgeV) = 0;  
       
       elseif(theta13>=theta12&&theta13>=theta23)
           %% delete edge indices(2)-del_joint_id
          edgeV=indices(2);
          spls_adj(edgeV,del_joint_id) = 0;
          spls_adj(del_joint_id,edgeV) = 0;  
       elseif(theta23>=theta13&&theta23>=theta12)
            %% delete edge indices(1)-del_joint_id
          edgeV=indices(1);
          spls_adj(edgeV,del_joint_id) = 0;
          spls_adj(del_joint_id,edgeV) = 0;  
       end
      %idx= find(joints(:,1)==del_joint_id);
      %joints(idx,:)=[];
      continue;
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % delete the small cycle.  
    for j = 1:length(cycle)  
        if ismember(cycle(j), joints(:,1))% a joint
            id = cycle(j);%找到关键点索引
            cycle(j)=[];
            cycle = [id, cycle];
            break;
        end
    end
    %将小的回路合并成一个点
    spls( cycle(1),: ) = mean( spls( cycle,: ) ); 
    kept_joints(end+1,:) = cycle(1);
    
    for j = 2:length(cycle)      
        spls( cycle(j),: ) = NaN; 
        segments( cycle(j) ) = 0; 
        % update the correspondents
        corresp( corresp==cycle(j) ) = cycle(1);
        
        tmp = find( joints(:,1)==cycle(j) );
        if ~isempty(tmp) % a joint  
            % update the A matrix
            for k=1:size(spls_adj,1)
                if spls_adj(cycle(j),k) == 1
                   spls_adj(cycle(1),k)=1; 
                   spls_adj(k,cycle(1))=1; 
                end
            end
            joints(tmp,:) = [];
        end    
        % remove the row
        spls_adj(cycle(j),:) = 0;
        spls_adj(:,cycle(j)) = 0;  
    end    
end

if show_cycles && ~isempty(kept_joints)
    figure('Name','Remove small cycles','NumberTitle','off');set(gcf,'color','white');movegui('west');
    subplot(1,2,1)
    scatter3(pts(:,1),pts(:,2),pts(:,3),30,'.');
    axis off; axis equal;
    title('Primitive tree')
    subplot(1,2,2)
    plot_skeleton(spls, spls_adj);
    scatter3(spls(kept_joints,1),spls(kept_joints,2),spls(kept_joints,3),200,'.b');hold on;    
    axis off; axis equal; 
    title('skeleton without cycles')
end


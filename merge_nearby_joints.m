function [spls, spls_adj, joints,roots] = merge_nearby_joints(spls,spls_adj, joints)
while true
    for i=joints'    
        edges = find( spls_adj(i,:)==1 );
        edges = edges(edges~=i);
        tmp = ismember(edges, joints');
        if sum(tmp)~=1, continue, end;

        j = edges(tmp);
        edges = find( spls_adj(j,:)==1 );
        edges = edges(edges~=j);
        tmp = ismember(edges, joints');
        if sum(tmp)~=1, continue, end;

      % update the location
        spls( i,: ) = mean( spls( [i,j],: ) );
        spls( j,: ) = NaN;
        % update the correspondents
       %corresp(corresp==j ) = i;

        % update the A matrix
        for k=1:size(spls_adj,1)
            if spls_adj( j,k ) == 1, 
                spls_adj( i,k )=1; 
                spls_adj( k,i)=1; 
            end
        end
        % remove the row
        spls_adj( j,: ) = 0;
        spls_adj( :,j ) = 0;

        %segments(j) = 0;
        j = find( joints(:,1)==j );
        %if root_id == j
        %    root_id = find( joints(:,1)==i );            
        %end
        joints(j,:) = [];
%         break;
    end
    if i == joints(end,1), break, end;
end

joints = zeros(0,1);%存储分叉点
roots = zeros(0,1);%存储分叉点
for i=1:size(spls,1)
    % 边折叠后的连接矩阵存在点和点本身连接
    %adj(i,i) = 0; %点和点自身的连接设置为0
    links = find( spls_adj(i,:)==1 );% 找到和i连接的点
    if length(links) >2 % 分叉点
        joints(end+1,1)=i;
    elseif length(links)==1
        roots(end+1,1) = i;
    end
end

%%
% figure('Name','Merge nearby joints','NumberTitle','off');set(gcf,'color','white');view3d rot;
% movegui('north');
% plot_skeleton(spls, spls_adj);
% axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-90,0);view3d rot;
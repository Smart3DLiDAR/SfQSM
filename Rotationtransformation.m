function [final_slices] = Rotationtransformation(slices)
%对点云数据进行旋转平移处理
for ii = 1:inf
    TMP = cell(length(slices),1);
    indx = false(length(slices),1);
    for i = 1:length(slices)
        P = slices{i};
        [m,~] = size(P);
        if m>=5
            P = P-ones(m,1)*(sum(P,1)/m);
            C = P.'*P./(m-1);
            [V,D] = eig(C);
            [~,ind] = sort(diag(D),'descend');
            Vs = V(:,ind);
            % 计算分层点云主方向
            MainD = Vs(:,1);
            OY = [0,1,0];
            theta1 = acos(dot(MainD,OY)/(norm(MainD)*norm(OY)));
            %旋转变换
            rot1 = [cos(theta1) sin(theta1) 0;...
                   -sin(theta1) cos(theta1) 0;...
                              0           0 0];
            proj_pnts_1=rot1 .* P(i,:); 
            TMP{i} = proj_pnts_1;
        end
            
        tslices = TMP{:};
        
        for j= 1:length(tslices)
            M = slices{i};
            dx=M(:,1);
            dy=M(:,2);
            dz=M(:,3);
            dtheta = rad2deg(dz./sqrt(dx.^2+dy.^2+dz.^2));
            if dtheta < 30
                indx(i) = true;
            else
                t = accumarray(dtheta,[1:length(dtheta)]',[],@(x) {x});
                slice = cellfun(@(x) M(x,:),t,'UniformOutput',0);
                
                rot2 = [cos(dtheta) 0 sin(dtheta);...
                                  0 0           0;...
                       -sin(dtheta) 0 cos(dtheta)];
                trans_slices =rot2 .* slice(j,:);
            end
            final_slices = [tslices;trans_slices];
    end
end
end
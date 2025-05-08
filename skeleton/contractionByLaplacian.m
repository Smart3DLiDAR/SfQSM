function [cpts, t, initWL, WC, sl] = contractionByLaplacian(P, options)
% 通过Laplacian算子来收缩点云

%% 参数设置
RING_SIZE_TYPE = 1;     %环多边形的类型
SHOW_CONTRACTION_PROGRESS = true;       %是否显示收缩过程
Laplace_type = 'mvc';     %位移类型（等角投影） 'curvaturenormal''conformal'
tc = getoptions(options, 'tc', Method.CONTRACT_TERMINATION_CONDITION);      %收敛阈值
iterate_time = getoptions(options, '迭代次数', Method.MAX_CONTRACT_NUM);     %最大迭代次数
initWL = getoptions(options, 'WL', Method.compute_init_laplacian_constraint_weight(P,Laplace_type)); %初始收缩权重
% 根据类型类初始化WC数组
if strcmp(Laplace_type,'mvc')
    WC = getoptions(options, 'WC', 1);%*10;
elseif strcmp(Laplace_type,'conformal')
    WC = getoptions(options, 'WC', 1);
else
    WC = getoptions(options, 'WC', 1);
end
WH = ones(P.npts, 1)*WC; % 初始约束权
sl = getoptions(options, 'sl', Method.LAPLACIAN_CONSTRAINT_SCALE); % 初始化尺度因子，论文中是2
WL = 3;%;initWL

sprintf(['1) k of knn: %d\n 2) termination condition: %f \n 3)' ...
    'Init Contract weight: %f\n 4) Init handle weight: %f\n 5) Contract scalar: %f\n' ...
    '6) Max iter steps: %d'], P.k_knn, tc, initWL, WC, sl,iterate_time)

%% 开始迭代
t = 1; 
tic
%左边的等式
L = -computePointLaplacian(P.pts,Laplace_type,P.rings, options);
%th = computeringfeature(P.pts,P.rings);
%A = [th.*L.*WL;sparse(1:P.npts,1:P.npts, WH)];
A = [L.*WL;sparse(1:P.npts,1:P.npts, WH)];
%右边的等式
b = [zeros(P.npts,3);sparse(1:P.npts,1:P.npts, WH)*P.pts];
cpts = (A'*A)\(A'*b);       %最小二乘法
disp('最小二乘法求解完成！');
toc

if SHOW_CONTRACTION_PROGRESS
    tic
    figure;
    axis off;axis equal;
    set(gcf,'Renderer','OpenGL');
    hold on;set(gcf,'color','white');
    camorbit(0,0,'camera'); axis vis3d; view(-90,0);    
    h1 = scatter3(P.pts(:,1),P.pts(:,2), P.pts(:,3),10,'b','filled');
    h2 = scatter3(cpts(:,1),cpts(:,2), cpts(:,3),10,'r','filled');
    title(['iterate ',num2str(t),' time(s)'])    
    disp('绘制点云');
    toc
end
%%
tic
sizes = Method.one_ring_size(P.pts, P.rings, RING_SIZE_TYPE);   % 使用最小半径
size_new = Method.one_ring_size(cpts, P.rings, RING_SIZE_TYPE);
a(t) = sum(size_new)/sum(sizes);

disp('计算局部范围完成！');
toc

%开始迭代
while t<iterate_time
    L = -computePointLaplacian(cpts,Laplace_type, P.rings, options);
    %th = computeringfeature(cpts,P.rings);
    WL = sl*WL;   
    %WL = WL; %收缩权WL sl*WL   
    if WL>Method.MAX_LAPLACIAN_CONSTRAINT_WEIGHT
        WL=Method.MAX_LAPLACIAN_CONSTRAINT_WEIGHT;
    end

    if strcmp(Laplace_type,'mvc')
        %WH = WC.*(sizes./size_new)*10;      % 更新约束权
        WH = WH.*(sizes./size_new);%
    else
        %WH = WC.*(sizes./size_new);
        WH = WC.*WH.*(sizes./size_new);
    end
    
    %WH(WH>Method.MAX_POSITION_CONSTRAINT_WEIGHT) = Method.MAX_POSITION_CONSTRAINT_WEIGHT;
    %A = real([th.*L.*WL;sparse(1:P.npts,1:P.npts, WH)]);
    A = real([L.*WL;sparse(1:P.npts,1:P.npts, WH)]);

    % 更新右边的等式
    b(P.npts+1:end, :) = sparse(1:P.npts,1:P.npts, WH)*cpts;
    tmp = (A'*A)\(A'*b);
    size_new = Method.one_ring_size(tmp, P.rings, RING_SIZE_TYPE);  
    a(end+1) = sum(size_new)/sum(sizes);

    tmpbox = Method.compute_bbox(tmp);
    if sum( (tmpbox(4:6)-tmpbox(1:3))> ((P.bbox(4:6)-P.bbox(1:3))*1.2) ) > 0
        break;
    end
    if a(t)-a(end)<tc || isnan(a(end))
        break;
    else 
        cpts = tmp;
    end
    
    t = t+1;    %更新迭代次数
    if SHOW_CONTRACTION_PROGRESS
        % 显示前后点云     
        delete(h1);delete(h2);      %释放资源
        figure;
        h1 = scatter3(P.pts(:,1),P.pts(:,2), P.pts(:,3),10,WH,'b','filled');
        h2 = scatter3(cpts(:,1),cpts(:,2), cpts(:,3),10,ones(P.npts,1)*WL,'r','filled');
        title(['迭代次数：',num2str(t)]); 
        drawnow; 
    end
end
clear tmp;
%sprintf('%t')

if SHOW_CONTRACTION_PROGRESS
    figure;
    plot(a);
    xlabel('迭代次数');
    ylabel('原始体积与当前体积之比');
    if t<iterate_time
        if a(end) < tc
            sprintf('迭代满足收敛条件!')
        else
            warning('存在迭代过程不正常!');
        end
    end
end


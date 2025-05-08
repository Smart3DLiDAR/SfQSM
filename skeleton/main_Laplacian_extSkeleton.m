%% ************清理环境**************
clear
close all;
clc

%% **********获取点云数据************
%获取点云数据
[fileName,pathName]=uigetfile('*.txt','Input Data-File');   %选择要进行计算的三维点云数据文件路径

if isempty(fileName) || length(fileName) == 1
    fprintf("未选择点云文件！\n");
    return;
end

data=importdata([pathName,fileName]);   %加载点云数据
pc=pointCloud(data(:,1:3));  %加载点云数据
figure;set(gcf,'color','white');hold on
pcshow(pc)
%pcshow(P.cpts)
%colorbar('west')
axis off;axis equal;set(gcf,'Renderer','OpenGL');
%% *************下采样*************
%该步骤主要是为了节省计算时间，如果要求精度，options.USING_UNDOWNSAMPLE需设置为false
options.USING_UNDOWNSAMPLE =1;  %是否进行下采样

if(options.USING_UNDOWNSAMPLE)
    gridStep = 5*PointsAverageSpacing(pc.Location);     %基于点密度进行下采样
    pc = pcdownsample(pc,'gridAverage',gridStep);
    P.pts = double(pc.Location);
    fprintf("下采样后的点数：%d\n",size(P.pts,1));
else
    P.pts = double(pc.Location);
end

%% *************计算过程*************
tic
P.npts = size(P.pts,1);
P.radis = ones(P.npts,1);
%P.pts = Method.normalize(P.pts); %缩放至单位包围盒并移动至原点
[P.bbox, P.diameter] = Method.compute_bbox(P.pts); %计算包围盒
disp('点云读取成功!');
toc

P.k_knn = Method.compute_k_knn(P.npts);
P.rings = computePointCloudRing(P.pts, P.k_knn, []);%计算点云中的一环邻域
disp('计算网格成功');
toc

%% **********通过Laplacian算子来收缩点云**********
%Features = ComputeFeature_mt(P.pts,32,false);
SHOW_JOINTS = true;
SHOW_ROOT_JOINT = true;
SHOW_CYCLES = true;
SHOW_IRRELEVANT_EXTRAMA = true;

Parameters.KnnNum=30; %%%% k parameter for skeleton 
Parameters.t2=0.2;
Parameters.alpha=1;%%% alpha parameter for stem 
Parameters.K1=5;%%%% K1 for fine segmentation 
Parameters.K2=Parameters.K1;
Parameters.beta=0.0;%%%%% beta for fine segmentation
Parameters.sigma=0.2; %%%%% sigma for fine segmentation
Parameters.idstop=0;
%%%% default parameter from paper(Cao et al 2010) 
sampleScale=0.015;%0.015
t1 = 0.1; % for inner branch nodes
a1 = pi*5.0/7.0; % for inner branch nodes, 
t2 = Parameters.t2;    
t3 = 8; % for small cycles;
DebugShow = 1;

[P.cpts, t, initWL, WC, sl] = contractionByLaplacian(P, options);

P.sample_radius = P.diameter*sampleScale;
P = rosa_lineextract(P,P.sample_radius, 1);
%P = RefineSKelton(P,DebugShow);
%% 计算点簇几何特征分离树叶和枝干
t = accumarray(P.corresp,[1:length(P.corresp)]',[],@(x) {x});
Seg = cellfun(@(x) P.pts(x,:),t,'UniformOutput',0);

%%  
Anisotropy = Cal_Anisotropy(Seg);
Num = cellfun(@(x) size(x,1),Seg);

AL = cell(length(Seg),1);%各向异性
NL = cell(length(Seg),1);%数量
for i = 1:length(Seg)
    PP = Seg{i};
    AL{i} = repmat(Anisotropy(i,:),size(PP,1),1);
    NL{i} = repmat(Num(i,:),size(PP,1),1);
end
AePts = cell2mat(AL);
NuPts = cell2mat(NL);
Pts = cell2mat(Seg);

% test a range of thresholds
Athres_list = 0.50:0.02:0.8;
Nthres_list = 350:50:800;%250-500
% all combinations of two thresholds
allc = combvec(Athres_list,Nthres_list)';

Freq = zeros(size(Pts,1),1);
for i = 1:size(allc,1)
    % find wood points based on two thresholds
    ia = AePts >= allc(i,1) & NuPts >= allc(i,2);
    % count the frequency of being identified as wood
    Freq = Freq + ia;
end
Pli = Freq/size(allc,1);

% build adjacent graph (similar to the one in "GraphRG")
Graph = build_graph_structure(Pts,30,0);
% initial class probability from our "Pli"
initial_classif = single([Pli,1-Pli]);
% alpha expansion method, output is the regularized label per point
[l_lin_potts, ~, ~, ~] = alpha_expansion(initial_classif, Graph, 0, 1, 5);
% we also record the label without regularization (directly from "Pli")
[~, l_baseline] = max(initial_classif,[],2);

%% restore original order (the original order was lost during wrapping labels to segments)
idx = knnsearch(Pts,P.cpts);

% wood label ->1 , leaf label -> 0
% BiLabel refers to point label without regularization !!
BiLabel = l_baseline(idx);
BiLabel(BiLabel~=1) = 0;
% BiLabel_Regu refers to point label with regularization !!
BiLabel_Regu = l_lin_potts(idx);
BiLabel_Regu(BiLabel_Regu~=1) = 0; 

%% if visualize
if plot == 1
    leaf = Pts(l_baseline==1,:);
    wood = Pts(l_baseline~=1,:);
    %wood = Pts(l_lin_potts==1,:);
    %leaf = Pts(l_lin_potts~=1,:);
    
    leaf2 = Pts(l_lin_potts==1,:);
    wood2 = Pts(l_lin_potts~=1,:);
    %wood2 = P.pts(BiLabel_Regu==1,:);
    %leaf2 = P.pts(BiLabel_Regu~=1,:);
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,3,1)
    pcshow(pointCloud(P.pts));view(110,20);
    grid off
    xlabel('X(m)');ylabel('Y(m)');zlabel('Z(m)');
    title('原始点云')
    subplot(1,3,2)
    pcshow(wood, repmat([0.4471, 0.3216, 0.1647],size(wood,1),1));
    hold on
    pcshow(leaf, repmat([0.2667, 0.5686, 0.1961],size(leaf,1),1));
    hold off
    grid off
    xlabel('X(m)');ylabel('Y(m)');zlabel('Z(m)');view(110,20);
    title('对点云枝叶分离')
    %figure;set(gcf,'color','white');
    subplot(1,3,3)
    pcshow(wood2, repmat([0.4471, 0.3216, 0.1647],size(wood2,1),1));
    hold on
    if ~isempty(leaf2)
        pcshow(leaf2, repmat([0.2667, 0.5686, 0.1961],size(leaf2,1),1));
    end
    hold off
    grid off
    xlabel('X(m)');ylabel('Y(m)');zlabel('Z(m)');view(110,20);
    %axis off;axis equal;
    title('检查标签')
end

%% 对于几何特征为树叶的骨架点进行删除
zmax = max(P.pts(:,3));
zmin = min(P.pts(:,3));
treeH = zmax - zmin;
if treeH > 5
    [sortid,~] = find(leaf2(:,3)<=1.3);
    leaf2(sortid,:) = [];
else
    [sortid,~] = find(leaf2(:,3)<=(1.3*treeH)/5);
    leaf2(sortid,:) = [];
end
idx1 = knnsearch(P.pts,leaf2);
fakenode = unique(P.corresp(idx1,:));

%% 根据最短路径算法连接断枝
branch = P;
for i=1:length(Seg)
    cur = Seg{i};
    m = size(cur,1);
    if m <10
        branch.spls(i,:)=NaN;
        branch.spls_adj(i,:)=0;
        branch.spls_adj(:,i)=0;
        [row,~] = find(branch.corresp(:)==i);
        branch.corresp(row,:) = NaN;
        branch.cpts(row,:) = NaN;
        branch.pts(row,:) = NaN;
    end
end

for i = 1:length(fakenode)
    ind = fakenode(i);
    branch.spls(ind,:)=NaN;
    branch.spls_adj(ind,:)=0;
    branch.spls_adj(:,ind)=0;
    [row,~] = find(branch.corresp(:)==ind);
    branch.cpts(row,:) = NaN;
    branch.pts(row,:) = NaN;
    branch.corresp(row,:) = NaN;
end
[m,~] = find(isnan(branch.cpts)==1);
branch.cpts(m,:) = [];
[n,~]=find(isnan(branch.pts)==1);
branch.pts(n,:) = [];
[n,~]=find(isnan(branch.corresp)==1);
branch.corresp(n,:) = [];
%[joints, roots,segments] = find_joints(branch.pts, branch.spls, branch.corresp, branch.spls_adj, SHOW_JOINTS);
[branch.spls, branch.corresp, branch.spls_adj, GP] = build_graph(branch.spls,branch.corresp,branch.spls_adj);
%[fakebranch,branch.spls,B] = remove_upper_nodes(branch.pts,branch.spls,B,true);% delate fake nodes;
[joints ,roots, branches]=find_Joints_mt(branch.spls, branch.spls_adj,branch.pts, false);
[branch.spls,branch.pts,branch.cpts,branch.corresp,branch.spls_adj] = adjust_stump(branch.pts,branch.cpts, branch.spls,branch.spls_adj , branch.corresp,joints);
[branch.spls, branch.corresp, branch.spls_adj, GP] = build_graph(branch.spls,branch.corresp,branch.spls_adj);
[joints ,roots, branches]=find_Joints_mt(branch.spls, branch.spls_adj,branch.pts, 1);

% 根据图的连通性确定断枝
[bin,binsize] = conncomp(digraph(branch.spls_adj),'Type','weak','OutputForm','cell');
if length(bin)~=1
    [newB,Gp,Link,branch.spls,branch.corresp,joints,roots,rootid] = findshorpath(branch.pts,branch.spls, ...
        branch.spls_adj,branch.corresp,bin,joints,roots,branch.diameter);
else
    [newB,Gp,Link,branch.spls, branch.corresp,joints,roots,rootid] = findshorpath2(branch.pts,branch.spls, ...
        branch.spls_adj,branch.corresp,roots,branch.diameter);
end

%% 重新寻找分支点 定义分支顺序
[skl,Path] = defineBranchOrder(Gp,branch.spls,joints,roots,rootid);
[branch.spls, branch.corresp, newB, newGraph] = build_graph(branch.spls,branch.corresp,newB);
[joints ,roots, branches]=find_Joints_mt(branch.spls, newB,branch.pts, true);
[Pid,Cid] = search_id2(rootid,branch.spls,joints,Link);
%[skl,Path] = defineBranchOrder(Gp,branch.spls,joints,roots);
[joints ,roots, branches]=find_Joints_mt(branch.spls, newB,branch.pts, true);
%% 分叉优化
spls = branch.spls;
[Similarity,FeasibleA,movej] = cosineSimilarity(branch.pts,spls,joints,Pid,Cid,branch.diameter);
%D = zeros(size(FeasibleA,1),3);%可能区域
for i = 1:size(Similarity,1)
    parentid = Pid(joints(i));
    parent = spls(parentid,:);%其父节点
    dis1 = pdist2(spls(joints(i),:),parent);
    
    CurSimilarity = Similarity{i,1};
    num = size(CurSimilarity,2);
    if size(CurSimilarity,1) > 1
        %CSimilarity = zeros(size(Similarity{i,1},1),2);
        for ii = 1:num
            [~,maxS1] = max(CurSimilarity(:,ii));
            FeasibleB(ii,:) = FeasibleA{i,1}(maxS1,:);
        end
        d = mean(FeasibleB);
        dis2 = pdist2(d,parent);
        
        if dis2 < dis1
            spls(joints(i),:) = d;
        else
            spls(joints(i),:) = branch.spls(joints(i),:);
        end
    else
        spls(joints(i),:) = branch.spls(joints(i),:);
    end
end
%spls(joints,:) = D;
figure('Name','初始骨架点');set(gcf,'color','white');
set(gca,'position',[0.1,0.1,0.8,0.8]);
%scatter3(branch.pts(:,1),branch.pts(:,2),branch.pts(:,3),30,'.');  hold on;
subplot(1,2,1)
showoptions.sizep=100;showoptions.sizee=2;
plotSkeleton(branch.spls, newB, showoptions);
axis off;axis equal;set(gcf,'Renderer','OpenGL');
title("原骨架")
subplot(1,2,2)
showoption.sizep=100;showoption.sizee=2;
showoption.colorp=[0,0,1];showoption.colore=[1, .0, .0];
plotSkeleton(spls, newB, showoption);
axis off;axis equal;set(gcf,'Renderer','OpenGL');
title("修改后骨架")

%% *************曲线拟合*************
% 分段平滑处理
for i = 1:max(skl(:,4))
    id1 = i;
    [sort1,~] = find(skl(:,4)==id1);
    cur = skl(sort1,:);
    for ii = 1:max(skl(sort1,5))
        id2 = ii;
        [sort2,~] = find(cur(:,5)==id2);
        if length(sort2)>2
            points = cur(sort2,1:3);
            idx = knnsearch(spls,points);
            centerX = spls(idx,1);
            centerY = spls(idx,2);
            centerZ = spls(idx,3);
            interX = smooth(centerZ,centerX,0.9,'rloess');%
            interY = smooth(centerZ,centerY,0.9,'rloess');%
            spls(idx,1) = interX;
            spls(idx,2) = interY;
        end
    end
end
figure('Name','初始骨架点');set(gcf,'color','white');
set(gca,'position',[0.1,0.1,0.8,0.8]);
%scatter3(branch.pts(:,1),branch.pts(:,2),branch.pts(:,3),20,'.k'); hold on;
showoptions.sizep=100;showoptions.sizee=2;
plotSkeleton(spls, newB, showoptions);hold on;
%showoption.sizep=100;showoption.sizee=2;
%showoption.colorp=[.0 .0 1];showoption.colore=[.0 .0 1];
%plotSkeleton(spls, newB, showoption);
axis off;axis equal;set(gcf,'Renderer','OpenGL');
title("分片平滑拟合")

%% *************圆柱拟合*************
Cylinders = fitCylinder(skl,branch.pts,spls,branch.corresp,Link,Pid,Cid,Path,joints);
%% 整体精度评定
[acc,sd,PP] =  Quantitativeevaluation(Cylinders,branch.pts);
[m,~] = find(isnan(acc)==1);
[n,~]=find(isnan(sd)==1);


%acc = QSM(1).pmdistance.CylDist;
%sd = QSM(1).pmdistance.sd;
[m,~] = find(isnan(acc)==1);
[n,~]=find(isnan(sd)==1);

clylabel = zeros(0,1);
label = QSM(1).cylinder.BranchOrder;
for i = 1:size(label,1)
    cur = label(i,1);
    clylabel(end+1,1) = cur+1;
end

clylabel = zeros(0,1);
for i = 1:size(Cylinders,1)
    curcyl = Cylinders{i,1};
    point1 = curcyl.Parameters(:,1:3);
    idx1 = knnsearch(skl(:,1:3),point1);
    label1 = skl(idx1,4);
    point2 = curcyl.Parameters(:,4:6);
    idx2 = knnsearch(skl(:,1:3),point2);
    label2 = skl(idx1,4);
    clylabel(end+1,1) = max([label1;label2]);
end

if isempty(m) && isempty(n)
    ACC1 = (mean(acc))*100;
    ACC2 = (mean(acc)+median(acc))/2*100;
    SD1 = (mean(sd))*100; 
    SD2 = (mean(sd)+median(sd))/2*100;
    
    accw = zeros(0,1);
    sdw = zeros(0,1);
    for i = 1:size(clylabel)
        curacc = acc(i,1);
        cursd = sd(i,1);
        curclylabel = clylabel(i,1);
        weight =1-(curclylabel/max(clylabel));
        accw(end+1,1) = curacc*weight;
        sdw(end+1,1) = cursd*weight;
    end
    if max(clylabel)==1
        ACC3 = ACC1;
        SD3 = SD1;
    else
        ACC3 = (sum(accw)/size(accw,1))*100;
        SD3 = (sum(sdw)/size(sdw,1))*100;
    end
else
    acc(m,:) = [];
    sd(n,:) = [];
    ACC1 = (mean(acc))*100;
    ACC2 = (mean(acc)+median(acc))/2*100;
    SD1 = (mean(sd))*100; 
    SD2 = (mean(sd)+median(sd))/2*100;
    
    clylabel1 = clylabel;
    clylabel1(m,:) = [];
    clylabel2 = clylabel;
    clylabel2(n,:) = [];
    accw = zeros(0,1);
    for i = 1:size(clylabel1)
        curacc = acc(i,1);
        curclylabel = clylabel1(i,1);
        weight =1-(curclylabel/max(clylabel1));
        accw(end+1,1) = curacc*weight;
    end
    
    sdw = zeros(0,1);
    for j = 1:size(clylabel2)
        cursd = sd(j,1);
        curclylabel = clylabel2(j,1);
        weight =1-(curclylabel/max(clylabel2));
        sdw(end+1,1) = cursd*weight;
    end
    if max(clylabel)==1
        ACC3 = ACC1;
        SD3 = SD1;
    else
        ACC3 = (sum(accw)/size(accw,1))*100;
        SD3 = (sum(sdw)/size(sdw,1))*100;
    end
end

%% 主干精度评定
[acc,sd] =  Quantitativeevaluation(Cylinders);
[m,~] = find(isnan(acc)==1);
[n,~]=find(isnan(sd)==1);

clylabel = zeros(0,1);
for i = 1:size(Cylinders,1)
    curcyl = Cylinders{i,1};
    point1 = curcyl.Parameters(:,1:3);
    idx1 = knnsearch(skl(:,1:3),point1);
    label1 = skl(idx1,4);
    point2 = curcyl.Parameters(:,4:6);
    idx2 = knnsearch(skl(:,1:3),point2);
    label2 = skl(idx1,4);
    clylabel(end+1,1) = max([label1;label2]);
end

if max(clylabel)<=3
    if isempty(m) && isempty(n)
        ACC1 = (mean(acc))*100;
        ACC2 = (mean(acc)+median(acc))/2*100;
        SD1 = (mean(sd))*100;
        SD2 = (mean(sd)+median(sd))/2*100;
    
        accw = zeros(0,1);
        sdw = zeros(0,1);
        for i = 1:size(clylabel)
            curacc = acc(i,1);
            cursd = sd(i,1);
            curclylabel = clylabel(i,1);
            weight =abs(1-(curclylabel/max(clylabel)));
            accw(end+1,1) = curacc*weight;
            sdw(end+1,1) = cursd*weight;
        end
        if max(clylabel)==1
            ACC3 = ACC1;
            SD3 = SD1;
        else
            ACC3 = (sum(accw)/size(accw,1))*100;
            SD3 = (sum(sdw)/size(sdw,1))*100;
        end
    else
        acc(m,:) = [];
        sd(n,:) = [];
        ACC1 = (mean(acc))*100;
        ACC2 = (mean(acc)+median(acc))/2*100;
        SD1 = (mean(sd))*100;
        SD2 = (mean(sd)+median(sd))/2*100;
    
        clylabel1 = clylabel;
        clylabel1(m,:) = [];
        clylabel2 = clylabel;
        clylabel2(n,:) = [];
        accw = zeros(0,1);
        for i = 1:size(clylabel1)
            curacc = acc(i,1);
            curclylabel = clylabel1(i,1);
            weight =abs(1-(curclylabel/max(clylabel1)));
            accw(end+1,1) = curacc*weight;
        end
    
        sdw = zeros(0,1);
        for j = 1:size(clylabel2)
            cursd = sd(j,1);
            curclylabel = clylabel2(j,1);
            weight =abs(1-(curclylabel/max(clylabel2)));
            sdw(end+1,1) = cursd*weight;
        end
        if max(clylabel)==1
            ACC3 = ACC1;
            SD3 = SD1;
        else
            ACC3 = (sum(accw)/size(accw,1))*100;
            SD3 = (sum(sdw)/size(sdw,1))*100;
        end
    end
else%计算主干部分
    [sortid,~] = find(clylabel(:,1)<=3);
    if isempty(m) && isempty(n)
        ACC1 = (mean(acc(sortid,1)))*100;
        ACC2 = (mean(acc(sortid,1))+median(acc(sortid,1)))/2*100;
        SD1 = (mean(sd(sortid,1)))*100;
        SD2 = (mean(sd(sortid,1))+median(sd(sortid,1)))/2*100;
    
        accw = zeros(0,1);
        sdw = zeros(0,1);
        for i = 1:size(sortid,1)
            curacc = acc(sortid(i,1),1);
            cursd = sd(sortid(i,1),1);
            curclylabel = clylabel(sortid(i,1),1);
            weight =abs(1-(curclylabel/3));
            accw(end+1,1) = curacc*weight;
            sdw(end+1,1) = cursd*weight;
        end
        if max(clylabel)==1
            ACC3 = ACC1;
            SD3 = SD1;
        else
            ACC3 = (sum(accw)/size(accw,1))*100;
            SD3 = (sum(sdw)/size(sdw,1))*100;
        end
    else
        acc(m,:) = [];
        sd(n,:) = [];
        ACC1 = (mean(acc(sortid,1)))*100;
        ACC2 = (mean(acc(sortid,1))+median(acc(sortid,1)))/2*100;
        SD1 = (mean(sd(sortid,1)))*100;
        SD2 = (mean(sd(sortid,1))+median(sd(sortid,1)))/2*100;
    
        clylabel1 = clylabel(sortid,1);
        clylabel1(m,:) = [];
        clylabel2 = clylabel(sortid,1);
        clylabel2(n,:) = [];
        accw = zeros(0,1);
        for i = 1:size(clylabel1)
            curacc = acc(i,1);
            curclylabel = clylabel1(i,1);
            weight =abs(1-(curclylabel/3));
            accw(end+1,1) = curacc*weight;
        end
    
        sdw = zeros(0,1);
        for j = 1:size(clylabel2)
            cursd = sd(j,1);
            curclylabel = clylabel2(j,1);
            weight =abs(1-(curclylabel/3));
            sdw(end+1,1) = cursd*weight;
        end
        if max(clylabel)==1
            ACC3 = ACC1;
            SD3 = SD1;
        else
            ACC3 = (sum(accw)/size(accw,1))*100;
            SD3 = (sum(sdw)/size(sdw,1))*100;
        end
    end
end

%%

for i = 1:length(Cylinders)
    seg = Cylinders{i,1};
    points = seg.fitpoints;
    figure(10);set(gcf,'color','white');
    pcshow(points);hold on;
    %scatter3(P.spls(curid,1),P.spls(curid,2),P.spls(curid,3),100,'.b');hold on;
    axis off;axis equal;set(gcf,'Renderer','OpenGL');
end

%% 可视化
%原始点云与初始的骨架点
figure('Name','初始骨架点');set(gcf,'color','white');
%set(gca,'position',[0.1,0.1,0.8,0.8]);
%subplot(1,2,1)
hold on;
pcshow(branch.pts);
showoptions.sizep=150;showoptions.sizee=2;
showoption.colorp=[0,0,1];showoption.colore=[1, .0, .0];
plotSkeleton(branch.spls, branch.spls_adj, showoptions);
axis off;axis equal;
figure('Name','初始骨架点');set(gcf,'color','white');
%set(gca,'position',[0.1,0.1,0.8,0.8]);
%showoptions.sizep=100;showoptions.sizee=1;
%pcshow(branch.pts);hold on
showoptions.sizep=150;showoptions.sizee=2;
showoption.colorp=[0,0,1];showoption.colore=[1, .0, .0];
plotSkeleton(branch.spls, newB, showoptions);
axis off;axis equal;

%%
clylabel = zeros(0,1);
for i = 1:size(Cylinders,1)
    curcyl = Cylinders{i,1};
    point1 = curcyl.Parameters(:,1:3);
    idx1 = knnsearch(skl(:,1:3),point1);
    label1 = skl(idx1,4);
    point2 = curcyl.Parameters(:,4:6);
    idx2 = knnsearch(skl(:,1:3),point2);
    label2 = skl(idx1,4);
    clylabel(end+1,1) = max([label1;label2]);
end

for i=1:max(clylabel)
    curlayer = i;
    [sort1,~] = find(clylabel(:,1)== curlayer);
    color = rand(1,3);
    for ii = 1:size(sort1)
        model = Cylinders{sort1(ii),1};
        figure(13);
        drawCylinder2(model.Parameters,'FaceColor',color);
        axis off; axis equal; camorbit(0,0,'camera');hold on
        title('Cylinderfit tree')
    end
end

%% *************可视化*************
for i= 1:size(branch.spls,1)
    curid = i;
    [id,~] = find(branch.corresp ==curid);
    seg = branch.pts(id,:);
    figure(10);set(gcf,'color','white');
    scatter3(seg(:,1),seg(:,2),seg(:,3),30,'.');hold on;
    scatter3(P.spls(curid,1),P.spls(curid,2),P.spls(curid,3),100,'.b');hold on;
    axis off;axis equal;set(gcf,'Renderer','OpenGL');
end

%分支可视化
sklmax = max(skl(:,4));
seg = zeros(0,3);
for i= 1:sklmax
    curid = i;
    [id,~] = find(skl(:,4) ==curid);
    curseg = skl(id,1:3);
    idx = knnsearch(branch.spls,curseg);
    for ii = 1:size(idx,1)
        [idc,~] = find(branch.corresp(:,1) == idx(ii));
        addseg = branch.pts(idc,:);
        seg = cat(1,seg,addseg);
    end
    figure(10);set(gcf,'color','white');
    scatter3(seg(:,1),seg(:,2),seg(:,3),30,'.');hold on;
    %scatter3(P.spls(curid,1),P.spls(curid,2),P.spls(curid,3),100,'.b');hold on;
    axis off;axis equal;set(gcf,'Renderer','OpenGL');
end


%原始点云与初始的骨架点
figure;set(gcf,'color','white');
%pcshow(pointCloud(P.pts))
pcshow(P.cpts)
colorbar('west')
axis off;axis equal;set(gcf,'Renderer','OpenGL');

minz = min(branch.spls(:,3));
[stumpid,~] = find(branch.pts(:,3)<= (minz+1.3));
%[~,maxid] = max(branch.spls(stumpid(:),3));
figure('Name','初始骨架点');set(gcf,'color','white');
%set(gca,'position',[0.1,0.1,0.8,0.8]);
%pcshow(pc);  hold on;
%scatter3(branch.spls(57,1),branch.spls(57,2),branch.spls(57,3),100,'.y'); 
pcshow(branch.pts); %hold on;
%showoptions.sizep=150;showoptions.sizee=2;
%showoption.colorp=[0,0,1];showoption.colore=[1, .0, .0];
%plotSkeleton(branch.spls, branch.spls_adj, showoptions);hold on;
%scatter3(branch.spls(stumpid(maxid),1),branch.spls(stumpid(maxid),2),branch.spls(stumpid(maxid),3),150,'.y'); 
axis off;axis equal;
%title('初始骨架')


figure;set(gcf,'color','white');
%pcshow(wood2, repmat([0.4471, 0.3216, 0.1647],size(wood2,1),1));hold on
%pcshow(branch.pts);
pcshow(PP(:,1:3),PP(:,4));
colormap('jet');
axis off;axis equal;

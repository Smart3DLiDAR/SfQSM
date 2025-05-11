% ************clean**************
clear
close all;
clc

%% **********Obtain point cloud data************
[fileName,pathName]=uigetfile('*.txt','Input Data-File');

if isempty(fileName) || length(fileName) == 1
    fprintf("The point cloud file was not selected！\n");
    return;
end
data=importdata([pathName,fileName]);

pc=pointCloud(data(:,1:3));    
pc = pcdenoise(pc,'NumNeighbors',10);

%% *************downsampling*************
%This step is mainly to save calculation time if accuracy is required，
% options.USING_UNDOWNSAMPLE needs to be set to false
options.USING_UNDOWNSAMPLE = 1;

if(options.USING_UNDOWNSAMPLE)
    gridStep = 5*PointsAverageSpacing(pc.Location);
    pc = pcdownsample(pc,'gridAverage',gridStep);
    P.pts = double(pc.Location);
    fprintf("The number of points after downsampling：%d\n",size(P.pts,1));
else
    P.pts = double(pc.Location);
end

%% 1. Skeleton points self-adjusting based on geometric features
% -----------1.1 Shrinking by Laplacian-based contraction----------------------------
P.npts = size(P.pts,1);
P.radis = ones(P.npts,1);
%P.pts = Method.normalize(P.pts); 
[P.bbox, P.diameter] = Method.compute_bbox(P.pts); 
P.k_knn = Method.compute_k_knn(P.npts);
P.rings = computePointCloudRing(P.pts, P.k_knn, []);
[P.cpts, ~, ~, ~, ~] = contractionByLaplacian(P, options);

% -----------1.2 Fasrthest distance spherical sampling -------------------------------
% To obtain more details, the initial sampling scale(R1) in the experiment 
% is not set large and can be changed as needed
R1 = 0.2;
SP = SelfadjustSampling(P,R1);
SP = AdjustStump(SP);
% -----------1.3 Skeleton points self-adjusting --------------------------------------
[~,rootid] = min(SP.spls(:,3));
[nodelabel,~,pli,Pid,Cid,Branch,~] = FindneedAdjustNodes(SP.spls,SP.seg,rootid);
newspls = adjustmentwithcirclefit(SP.spls,SP.seg,rootid,2/3);
[newspls1,Radius] = adjustmentwithangle(newspls,SP.pts,rootid,nodelabel,Branch,Pid,Cid);

%% 2. Edge weight definition for the tree graph
[~,rootid] = min(newspls1(:,3));
[~,dis] = knnsearch(newspls1,newspls1,'k',6);
OptimalR = mean(dis(:,end))+std(dis(:,end));
[Gpdis,adjdis] = BuildGraph(newspls1,newspls.Seg,OptimalR);
%%----Obtain the skeleton lines according to the shortest path algorithm----------
[newtR,newD] = shortestpathtree(Gpdis,rootid,[1:size(newspls1,1)],'OutputForm','cell');
newadj = zeros(size(newspls1,1),size(newspls1,1));
for i = 1:size(newspls1,1)
    curpath = newtR{i,:};
    if length(curpath)>1
        front = curpath(end-1);
        newadj(i,front) = 1;
        newadj(front,i) =1;
    end
end

newC = cell2mat(newtR');
newpli = tabulate(newC);
[newCid,newPid ] = define_topologicalrelations(newspls1,newadj,rootid);
%%----Define the order of all branches and trunks---------------------------------
[newBranch,~] = Stratifybrancheslevel(newspls1,newadj,newCid,newpli,rootid);
%%----Calculate the branch length at each level-----------------------------------
newBranchL = CalculatBranchLength(newspls1,newPid,newadj,rootid,newBranch);
%%----Pruned Branch---------------------------------------------------------------
[finalspls,finaladj,newPid,newCid,link,Radius2] = Findnodesneedpruned(newspls1, ...
    newadj,newspls.Seg,newspls.Radius,newBranch,newBranchL,newPid,rootid);
finnalGp = graph(finaladj);
[~,rootid] = min(finalspls(:,3));
[finnaltR,finnalD] = shortestpathtree(finnalGp,rootid,[1:size(finalspls,1)],'OutputForm','cell');
finnalC = cell2mat(finnaltR');
finnalpli = tabulate(finnalC);
%%----Update the branch order and length after pruning ---------------------------
[finnalBranch,branch2] = Stratifybrancheslevel(finalspls,finaladj,newCid,finnalpli,rootid);
finnalBranchL = CalculatBranchLength(finalspls,newPid,finaladj,rootid,finnalBranch);
%% 3. Fractal self-similarity optimization for individual tree modeling
subbranchlength = findsubbranchlength(finalspls,finnalBranch,finnalGp,newPid);
[normallabel,newlink,Radius21,TH] = findabnormalRadius(Radius2,...
    finnalBranch,rootid,newCid,newPid,subbranchlength,link);
[newRadius1,newlink1,fixedR] = fixabnormalRadius(Radius21,finnalBranch,newlink,...
    rootid,newCid,newPid,subbranchlength,normallabel,TH);
%----build individual tree model---------------------------------------------------
cyl = branchmodeling(finalspls,finnalBranch,newRadius1,fixedR,branch2,newPid);
figure(3);set(gcf,'color','white');
set(gcf,'unit','centimeters','position',[6 6 8 10]);
plot_cylinder_model(cyl,'order',3,10,0.8);
camlight('headlight');material('dull');
axis off;axis equal;




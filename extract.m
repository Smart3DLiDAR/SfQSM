%��ȡ����
clear
close all
clc

%��ȡ��������
[fileName,pathName]=uigetfile('*.txt','Input Data-File');   %ѡ��Ҫ���м������ά���������ļ�·��

if isempty(fileName) || length(fileName) == 1
    fprintf("δѡ������ļ���\n");
    return;
end

Data=importdata([pathName,fileName]);   %���ص�������
pc=pointCloud(Data(:,1:3));
figure;
pcshow(pc);
location = pc.Location;
n = size(location,1);
e1 = zeros(n,1);
e2 = zeros(n,2);
e3 = zeros(n,3);
kdtree = KDTreeSearcher(location);
nearPoints = knnsearch(kdtree,location,"K",100);
flags = zeros(n,1);

%���ò���
lambda3_upper = 0.1;
lambda2_lower = 0.35;

%����ÿ�����Э����
for i=1:n
    nearIndex = nearPoints(i,:);
    nPoints = location(nearIndex,:);
    nPoints = nPoints - mean(nPoints);     %ȥ���Ļ�
    sysMatrix = nPoints'*nPoints;       %����Э�������
    s = svd(sysMatrix);
    tol = sum(s);
    lambda1 = s(1)/tol;
    lambda2 = s(2)/tol;
    lambda3 = s(3)/tol;
    if lambda3 <= lambda3_upper && lambda2 >= lambda2_lower
        flags(i) = 1;
    end
end
stem1 = location(flags==1,:);

%��Ҫ�������÷�����(ˮƽ��)
pcStem1 = pointCloud(stem1);
normals = pcnormals(pcStem1,100);       %��ȡÿһ����ķ�����
nx = normals(:,1);
ny = normals(:,2);
nz = normals(:,3);
angle = abs(asin(nz./sqrt(nx.^2+ny.^2+nz.^2)));
angleToRad = pi/180;
aThresh = 10;       %�Ƕȵ�����

stem2 = stem1(angle< aThresh * angleToRad,:);
% offsetm2 = location(angle > aThresh * angleToRad,:);
pcStem2 = pointCloud(stem2);
figure
pcshow(pcStem2)
title('������ȡ������')

%ͨ��ŷʽ�����һ������ȡ���ɲ���
minDistance = 0.2;
minPoints = 1000;
[labels,numClusters] = pcsegdist(pcStem2,minDistance);

idxValidPoints = find(labels);      %������Ϊ����ʾ���
labelColorIndex = labels(idxValidPoints);
segmentedPtCloud = select(pcStem2,idxValidPoints);
figure
colormap(hsv(numClusters))
pcshow(segmentedPtCloud.Location,labelColorIndex)
title('���ƾ���')

%ʹ��RANSAC�㷨�������Բ��
labelIndexs = 1:length(labels);
classNum = length(unique(labels));
referenceVector = [0,0,1];
maxDist = 0.01;

figure
hold on
grid on
rotate3d on
for i=1:classNum-1
    classIndex = labelIndexs(labels==i); 
    classPC = select(pcStem2,classIndex);
    [model,inlierIndices,outlierIndices] = pcfitcylinder(classPC,maxDist...
        ,referenceVector,25,"MaxNumTrials",10000,"Confidence",99);
    normalOri = model.Orientation./norm(model.Orientation);
    plot3(classPC.Location(:,1),...
        classPC.Location(:,2),...
        classPC.Location(:,3),'b.')
    DrawCylinderV1(model.Center,normalOri,...
        model.Height,model.Radius,[1,0,0])

    if(model.Radius > 0.5 || model.Height < 0.5)
       labels(classIndex) = 0;
    end
end
fprintf("RANSAC���Բ��\n");

result = labelIndexs(labels>0);
pcStem3 = select(pcStem2,result);
figure
pcshow(pcStem3)
title("���ɵ���")

% plot3(stem2(:,1),stem2(:,2),stem2(:,3),'b.')
% plot3(offsetm1(:,1),offsetm1(:,2),offsetm1(:,3),'r.')

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
Data=pc.Location(:,1:3);     %ȡ���ݵ�һ������
figure;
pcshow(pc);
n = size(Data,1);

x = Data(:,1);
y = Data(:,2);
z = Data(:,3);

xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

% �����ɽ�����Ƭ
zMin = min(z);
zMax = max(z);
h=zMax-zMin;    %�߳����ֵ����¸�֮��
d=0.1;       %��Ƭ�ĺ��Ϊ2d
cycleNum=ceil(h/d); %ѭ������
S=zeros(1,cycleNum);    %�����洢��Ƭ���������
V=zeros(1,cycleNum);    %�����洢����̨�����������

%�Ե������ݽ�����Ƭ����
slices=cell(1,cycleNum);
for i=1:size(Data, 1)
    tx=Data(i,1);
    ty=Data(i,2);
    tz=Data(i,3);
    index=floor((tz-zMin)/d) + 1;%��Ƭ���
    zUnder=zMin+(index -1) * d;%��һ��Ƭ��Сֵ
    zUp=zUnder + 2*d;%��Ƭ���ֵ
    if tz > zUnder && tz < zUp 
       slices{1,index}=[slices{1,index};tx ty tz];  %�洢ÿһ����Ƭ�еĵ�
    end
end

%% ���ӻ�����Ƭ����
cmap = hsv(length(slices));
pz = randperm(size(cmap,1),size(cmap,1));
cmap = cmap(pz,:);
col = cell(length(slices),1);
 for i =1:length(slices)
     ww = slices{i};
     col(i) = {repmat(cmap(i,:),size(ww,1),1)};
 end
 figure;pcshow(cell2mat(slices'),cell2mat(col));grid off;
 
centerPoint = [];
%ʹ��RANSAC�㷨�����Բ
for i=1:size(slices,2)
    points = slices{i};
    if size(points,1) < 3
        continue;
    end
    zVal = zeros(size(points,1),1);
    tmp = zMin + i*d - d/2;
    zVal(:) = tmp;
    [model,inlierIndices] = RANSACCircle(points,0.02,1000);
    plot3(points(:,1),points(:,2),zVal,".");
    DrawCircle(model.circleCenter,model.r,tmp)
    centerPoint = [centerPoint;model.circleCenter tmp model.r];
end
figure;
hold on
grid on
rotate3d on
pcshow(pc)
plot3(centerPoint(:,1),centerPoint(:,2),centerPoint(:,3), 'r-')

%ʹ�����������������߽���ƽ������
centerX = centerPoint(:,1);
centerY = centerPoint(:,2);
centerZ = centerPoint(:,3);

interX = smooth(centerZ,centerX,0.8,'rloess');
interY = smooth(centerZ,centerY,0.8,'rloess');

figure
hold on
grid on
rotate3d on
pcshow(pc)
plot3(interX,interY,centerZ,'r.')

SkeletonPoint = [-23.22 25.54 8.03];
R_search = 0.12;%���������뾶
idx = rangesearch(SkeletonPoint,Data,1.2);
ImpactPoint=[];
for i = 1:length(idx)
    if idx{i}==1
       ImpactPoint = Data(i,:);%��ȡ��������
    end
    ImpactPoint=ImpactPoint;
end
[newSkeletonPoint] = SpaceColonization(SkeletonPoint,ImpactPoint,R_search,d);


clear
close all
clc

%获取点云数据
[fileName,pathName]=uigetfile('*.txt','Input Data-File');   %选择要进行计算的三维点云数据文件路径

if isempty(fileName) || length(fileName) == 1
    fprintf("未选择点云文件！\n");
    return;
end

Data=importdata([pathName,fileName]);   %加载点云数据
pc=pointCloud(Data(:,1:3));
Data=pc.Location(:,1:3);     %取数据的一到三列
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

% 对树干进行切片
zMin = min(z);
zMax = max(z);
h=zMax-zMin;    %高程最大值与冠下高之差
d=0.1;       %切片的厚度为2d
cycleNum=ceil(h/d); %循环次数
S=zeros(1,cycleNum);    %创建存储切片面积的数组
V=zeros(1,cycleNum);    %创建存储各个台体体积的数组

%对点云数据进行切片处理
slices=cell(1,cycleNum);
for i=1:size(Data, 1)
    tx=Data(i,1);
    ty=Data(i,2);
    tz=Data(i,3);
    index=floor((tz-zMin)/d) + 1;%切片间隔
    zUnder=zMin+(index -1) * d;%下一切片最小值
    zUp=zUnder + 2*d;%切片最大值
    if tz > zUnder && tz < zUp 
       slices{1,index}=[slices{1,index};tx ty tz];  %存储每一个切片中的点
    end
end

%% 可视化各切片点云
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
%使用RANSAC算法来拟合圆
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

%使用样条函数对中心线进行平滑处理
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
R_search = 0.12;%设置搜索半径
idx = rangesearch(SkeletonPoint,Data,1.2);
ImpactPoint=[];
for i = 1:length(idx)
    if idx{i}==1
       ImpactPoint = Data(i,:);%获取待搜索点
    end
    ImpactPoint=ImpactPoint;
end
[newSkeletonPoint] = SpaceColonization(SkeletonPoint,ImpactPoint,R_search,d);


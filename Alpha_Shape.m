clear
close all;
clc


%获取点云数据


Data=importdata('3612_leaf.txt');   %加载点云数据
Data=Data(:,1:2);

if size(Data,1)>10000  %抽稀
    rnum = 10000;
else
    rnum = size(Data,1);
end

rng('default')  % For reproducibility
randindx = randperm(size(Data,1),rnum);
Data = Data(randindx,:);


Lr=zeros(1, 4);
lm=1;

tic
X=Data(:,1);        
Y=Data(:,2);
% S=[X Y];        
plot(X,Y,'.')
hold on;

alpha=1;
% adjmatrix=zeros(n);
n=size(X,1);
kdtree=KDTreeSearcher(Data);
%p0为中心2r为半径画圆，遍历查询小于2倍alpha值的点，并将符合条件的点与当前点关联起来   
[points,distance]=rangesearch(kdtree,Data,2*alpha);

m=size(points,1);
for i=1:m
    adjPoints=points{i};
    adjDistance = distance{i};
    cellLength=size(adjPoints,2);
   
    for j=2:cellLength
        index=adjPoints(j);
        sqDist=adjDistance(j)^2;
%       sqDist=sum((S(i,:)-S(index,:)).^2);
        H=sqrt(alpha*alpha/sqDist-0.25);        %圆心到弦之间的距离除以弦的距离
        cirCenter_X1=X(i)+0.5*(X(index)-X(i))-H*(Y(index)-Y(i));
        cirCenter_Y1=Y(i)+0.5*(Y(index)-Y(i))-H*(X(i)-X(index));
        cirCenter_X2=X(i)+0.5*(X(index)-X(i))+H*(Y(index)-Y(i));
        cirCenter_Y2=Y(i)+0.5*(Y(index)-Y(i))+H*(X(i)-X(index));
        adjPointsXY = Data(adjPoints(2:end),:);
        D1=sqrt(sum((adjPointsXY - [cirCenter_X1,cirCenter_Y1]).^2,2));
        D2=sqrt(sum((adjPointsXY - [cirCenter_X2,cirCenter_Y2]).^2,2));
        D1minVal=min(D1);
        D2minVal=min(D2);
        if D1minVal>=alpha || abs(D1minVal-alpha)<0.000001 || D2minVal>=alpha || abs(D2minVal-alpha)<0.000001
            line([X(i),X(index)],[Y(i),Y(index)]);
            Lr(lm, 1:4)=[X(i),Y(i), X(index),Y(index)];  %%一行四列，相邻点坐标位于一行
            lm=lm+1;
        end
    end
end
clear i j
m=1;
axis equal;
toc

  %%  display画树冠点；边界点连线 ；Lr顺序排列可以看出有一半是重复行
% scatter(Lr(1, 3), Lr(1, 4), '.', 'r')

% for i=1:size(Lr, 1)
%     line([Lr(i, 1),Lr(i, 3)],[Lr(i, 2),Lr(i, 4)]);
% end
% clear i 

% [~,bb] = sort(Lr(:, 1));
% Lr=Lr(bb, :); 

  %%  边界点连续排列
for i=1:size(Lr, 1)-1    
    for j=i+1:size(Lr, 1)        
        if double(isequal(Lr(i, 3:4), Lr(j, 1:2)))==1 && double(isequal(Lr(i, 1:2), Lr(j, 3:4)))==1
            Lr(j, :)=0;  %重复行置0
        elseif double(isequal(Lr(i, 3:4), Lr(j, 1:2)))==1 && double(isequal(Lr(i, 1:2), Lr(j, 3:4)))==0
            cc=Lr(j, :);
            Lr(j, 1:4)=Lr(i+1, :);
            Lr(i+1, 1:4)=cc;  %一行四列，如1122,2233将2233的行放到1122行下形成连续
        end
    end
end
clear i j

  %%  去零行
idx2=find(Lr(:, 1)~=0);
Lr=Lr(idx2, :);


  %%  display连续点挨个标记（检查代码成果）
% for i=1:size(Lr, 1)
%     scatter(Lr(i, 1), Lr(i, 2), '.', 'r')
%     pause
% end

  %%  画圆
for k=1:size(Lr, 1)
    S=(Lr(k, 1)-Lr(k, 3))^2+(Lr(k, 2)-Lr(k, 4))^2;
    H=sqrt(alpha^2 / S -1/4);

    cx=Lr(k, 1)+ (Lr(k, 3)-Lr(k, 1))/2 -H*(Lr(k, 4)-Lr(k, 2));
    cy=Lr(k, 2)+ (Lr(k, 4)-Lr(k, 2))/2 -H*(Lr(k, 1)-Lr(k, 3));

    th=0:pi/50:2*pi;
    Cx=alpha*cos(th)+cx;
    Cy=alpha*sin(th)+cy;

    plot(Cx, Cy, 'Color', 'g');
end
clear k

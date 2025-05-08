function [W,S] = computePointWeight(pts, type, rings, options)

% 计算权重矩阵
options.null = 0;

if nargin<2
    type = 'conformal';
end
        
switch lower(type)
    case 'combinatorial'
        W = compute_point_weight_combinatorial(pts, rings);
    case 'distance'
       warning('not implemented!'); 
    case 'spring'
        W = computePointPeightSpring(pts, rings);        
    case {'conformal','dcp'} % conformal laplacian  
        W = computePointWeightDcp(pts, rings);
    case 'laplace-beltrami' %
        W = computePointWeightDcp(pts, rings)*0.5;
    case 'mvc'% mvc laplacian
        W = computePointWeightMvc(pts, rings);  
    case 'curvaturenormal'
        W = computePointcurvaturenormaloperator1(pts, rings);
    otherwise
        error('Unknown type!!')
end

%#########################################################################
function W = compute_point_weight_combinatorial(points, rings)
n = length(points);
% W = size(n,n);  %这个更�?�用与小数据
W = sparse(n,n);	%创建�?个稀疏矩阵，对于大数据�?�言可以节省内存空间	
for i = 1:n
    ring = rings{i};
    if ring(1) == ring(end)
        ring = ring(1,1:(end-1));
    end
    for j = ring
        W(i,j) = 1.0;
    end
end
function W = computePointPeightSpring(points, rings)
n = length(points);
% W = size(n,n);
W = sparse(n,n);
for i = 1:n
    vi = points(i,:);
    ring = rings{i};
    if ring(1) == ring(end)
        ring = ring(1,1:(end-1));
    end
    for j = ring                
        vj = points(j,:);        
        W(i,j) = 1./sqrt(sum((vi-vj).^2));
    end
    tmp = sum(W(i,:));
    if tmp>10000
        W(i,:) = W(i,:)*10000/tmp;
    end
end
%#########################################################################
function W = computePointWeightDcp(points, rings)
%通过�?部网格多边形计算权重
n = length(points);
% W = size(n,n);
W = sparse(n,n);
for i = 1:n
    ring = rings{i};        %获取�?部多边形

    tmp = size(ring,2)-1;
    for ii = 1: tmp
        j = ring(ii);   %获取端点
        k = ring(ii+1);
        vi = points(i,:);   %���������
        vj = points(j,:);
        vk = points(k,:);
        
        % 
        u = vk-vi; v = vk-vj;
        cot1 = dot(u,v)/norm(cross(u,v));
        W(i,j) = W(i,j) + cot1;
        u = vj-vi; v = vj-vk;
        cot2 = dot(u,v)/norm(cross(u,v));
        W(i,k) = W(i,k) + cot2;        
    end
    
    tmp = abs(sum(W(i,:)));
    if tmp>10000
        W(i,:) = W(i,:)*10000/tmp;
    end
end
function W = computePointWeightMvc(points, rings)
n = length(points);
% W = size(n,n);
W = sparse(n,n);
for i = 1:n
    ring = rings{i};
    curpoints = points(ring,1:3);
    curpoints = unique(curpoints,'rows');
    [~,~,latent] = pca(curpoints(:,1:3), 'Economy', false); %ʹ��PCA�㷨����ֲ�����������
    L = (latent(1)-latent(2))/latent(1);

    tmp = size(ring,2)-1;
    for ii = 1: tmp
        j = ring(ii); k = ring(ii+1);
        vi = points(i,:);
        vj = points(j,:);
        vk = points(k,:);
        
        % angles
        alpha = myangle(vi-vk,vi-vj);
        % add weight
        W(i,j) = W(i,j) + tan( 0.5*alpha )/sqrt(sum((vi-vj).^2));
        W(i,k) = W(i,k) + tan( 0.5*alpha )/sqrt(sum((vi-vk).^2));
        %W(i,j) = W(i,j) + L*(tan( 0.5*alpha )/sqrt(sum((vi-vj).^2)));
        %W(i,k) = W(i,k) + L*(tan( 0.5*alpha )/sqrt(sum((vi-vk).^2)));
    end
    
    tmp = sum(W(i,:));
    if tmp>10000
        W(i,:) = W(i,:)*10000/tmp;
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = computePointcurvaturenormaloperator1(points, rings)
n = length(points);
% W = size(n,n);
W = sparse(n,n);
%L = zeros(n,1);
for i = 1:n
    ring = rings{i};
    curpoints = points(ring,1:3);
    curpoints = unique(curpoints,'rows');
    [~,~,latent] = pca(curpoints(:,1:3), 'Economy', false); %ʹ��PCA�㷨����ֲ�����������
    L = (latent(1)-latent(2))/latent(1);
    
    tmp = size(ring,2)-1;
    area = zeros(tmp,1);
    %L = zeros(1,n);
    for ii = 1: tmp
        j = ring(ii);   %获取端点
        k = ring(ii+1);
        vi = points(i,:);   %���������
        vj = points(j,:);
        vk = points(k,:);
        %L(j)= norm(vi-vj);
        
        % 
        alpha = myangle(vi-vk,vi-vj);
        
        Lij = norm(vi-vj);
        W(i,j) = W(i,j) + (tan( 0.5*alpha )/sqrt(sum((vi-vj).^2)))*Lij;%*()
        Lik = norm(vi-vk);
        W(i,k) = W(i,k) + (tan( 0.5*alpha )/sqrt(sum((vi-vk).^2)))*Lik;%*;(norm(vi-vk))
        area(ii) = mianji(vi,vj,vk);
    end
    S = sum(area); 
    
    all = (tmp/(4*S));
    for m = 1:length(ring)
        W(i,ring(m)) = all;
    end
    tmp = sum(W(i,:));
    
    if tmp>100%10000L*
        W(i,:) = W(i,:)*100/tmp;
    end
    %A = L*((1/(4*S))*tmp);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = computePointcurvaturenormaloperator2(points, rings)
n = length(points);
% W = size(n,n);
W = sparse(n,n);
%L = zeros(n,1);
for i = 1:n
    ring = rings{i};
    curpoints = points(ring,1:3);
    curpoints = unique(curpoints,'rows');
    [~,~,latent] = pca(curpoints(:,1:3), 'Economy', false); %ʹ��PCA�㷨����ֲ�����������
    L = (latent(1)-latent(2))/latent(1);
    
    tmp = size(ring,2)-1;
    area = zeros(tmp,1);
    %L = zeros(1,n);
    for ii = 1: tmp
        j = ring(ii);   %获取端点
        k = ring(ii+1);
        vi = points(i,:);   %���������
        vj = points(j,:);
        vk = points(k,:);
        %L(j)= norm(vi-vj);
        
        % 
        u = vk-vi; v = vk-vj;
        cot1 = dot(u,v)/norm(cross(u,v));
        Lij = norm(vi-vj);
        W(i,j) = W(i,j) + cot1*Lij;%*()
        u = vj-vi; v = vj-vk;
        cot2 = dot(u,v)/norm(cross(u,v));
        Lik = norm(vi-vk);
        W(i,k) = W(i,k) + cot2*Lik;%*;(norm(vi-vk))
        area(ii) = mianji(vi,vj,vk);
    end
    S = sum(area); 
    
    tmp = sum(W(i,:));
    all = L*((1/(4*S))*tmp);
    for m = 1:length(ring)
        W(i,ring(m)) = all;
    end
    if tmp>10000
        W(i,:) = W(i,:)*10000/tmp;
    end
    %A = ((1/(4*S))*tmp);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta = myangle(u,v)

du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
du = max(du,eps); dv = max(dv,eps);
beta = acos( sum(u.*v) / (du*dv) ); 



function s=mianji(A,B,C)
if length(A)==2  %���������Ƕ�άƽ�����꣬�����ά
    AB=[B-A 0];
    BC=[C-B 0];
elseif length(A)==3  %������������ά�ռ�����
    AB=B-A;
    BC=C-B;
end
Z=cross(AB,BC);  %���
s=1/2*norm(Z);  %ȡģ


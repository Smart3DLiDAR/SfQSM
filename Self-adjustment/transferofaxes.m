function [Axis,newP,tempseg] = transferofaxes(curid,centroid,spls,Seg)
%----input-----------------------------------------------------------------------
% spls : Skeleton points in the original coordinate system
% Seg : The set of contraction points corresponding to the skeleton points 
%       in the original coordinate system
% Branch :  All levels Branches  
%----output-----------------------------------------------------------------------
% Axis : Rotation matrix
% newP : The tree points in the new coordinate system
% tempseg: The set of contraction points corresponding to the skeleton points 
%          in the new coordinate system
%% -------------------------------------------------------------------------------
neispls = knnsearch(spls,spls(curid,:),'k',6);
tempseg = [];
for l = 1:length(neispls)
    cursplsseg = Seg{neispls(l),1};
    tempseg = [tempseg;cursplsseg];
end
[~,maxid] = max(spls(neispls,3));
[~,minid] = min(spls(neispls,3));
p1 = spls(neispls(minid),:);
p2 = spls(neispls(maxid),:);
    
[V,~,S] = pca(tempseg);
feature = S(2)/S(3);
if feature>2
    [~,maxid] = max(V(3,:));
    dir1 = V(:,maxid)/norm(V(:,maxid));
    dir1 = dir1';
else
    dir1=(p2-p1)./norm(p2-p1);
end
    
pCoords=find_ProjCoord_mt(dir1,centroid,tempseg);
coeff = pca(pCoords);
if(size(coeff,2)==3)
    dir2=coeff(:,1)';
    dir3=coeff(:,2)';
elseif(size(coeff,2)==2)
    dir2=coeff(:,1)';
    dir3=[0 0 0];
else
    dir2=[0 0 0];
    dir3=[0 0 0];
end
Axis=[centroid;dir2;dir3;dir1];
    
PNum=size(tempseg,1);
Dir1=repmat(Axis(2,:),[PNum 1]);
Dir2=repmat(Axis(3,:),[PNum 1]);
Dir3=repmat(Axis(4,:),[PNum 1]);
Displace=tempseg-Axis(1,:);
V1=dot(Displace',Dir1');
V2=dot(Displace',Dir2');
V3=dot(Displace',Dir3');
newP=[ V1' V2' V3'];

end

function [THETA,RHO]=pertransform(XY)

%% downsample to ~2k points to increase speed
if size(XY,1)>2000
    rnum = 2000;
else
    rnum = size(XY,1);
end

rng('default')  % For reproducibility
randindx = randperm(size(XY,1),rnum);
XY = XY(randindx,:);
    
xx = min(XY(:,1)):0.02:max(XY(:,1));
yy = min(XY(:,2)):0.02:max(XY(:,2));
    
[xx,yy] = meshgrid(xx,yy);%Generate a 2D grid
pot_loc = [xx(:),yy(:)];
    
N = nan(size(pot_loc,1),1);
for ii = 1:size(pot_loc,1)
    dist = sqrt((XY(:,1) - pot_loc(ii,1)).^2 + (XY(:,2) - pot_loc(ii,2)).^2);
    N(ii) = sum(abs(dist - prctile(dist,5))<0.02);
end
ia = find(N == max(N),1,'first');
L1M = pot_loc(ia,:);
    
r_init = median(sqrt(sum((XY - repmat(L1M,size(XY,1),1)).^2,2)));%radius
Tp = [L1M,r_init];%Obtain the center of the Probably circle¡¾x,y,r¡¿
    
%% centered by rough center 
x = XY(:,1) - Tp(1);
y = XY(:,2) - Tp(2);

[THETA,RHO] = cart2pol(x,y);
end
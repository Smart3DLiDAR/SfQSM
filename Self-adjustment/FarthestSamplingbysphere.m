function [spls, corresp] = FarthestSamplingbysphere(P,Linearity,R)
%% Skeleton points obtained with adaptive farthest distance sphere sampling
% The sampling function based on the linear characteristics of the contracted point cloud 
% is used to adaptively adjust the radius of the sampling sphere
% --input------------------------------------------------------------
% P: The contraction point clouds
% Linearity: The linear characteristics of the contraction point clouds
% R: Initial radius
% --output------------------------------------------------------------
% spls: Skeleton points
% corresp: The correspondence between point clouds and skeleton points
%
% @author: JJCAO
% Changed from Andrea Tagliasacchi's code

%% visual debug conditions
SHOW_SAMPLING_PROGRESS = false;
SHOW_RESULTS =false;

%%
if SHOW_SAMPLING_PROGRESS || SHOW_RESULTS
    close all;
    figure(1); movegui('northwest');set(gcf,'color','white');hold on;
    plot3( cpts(:,1), cpts(:,2), cpts(:,3), '.r', 'markersize', 1);
    axis off; axis equal;set(gcf,'Renderer','OpenGL');
end

cpts = P.cpts;
pts = P.pts;
% pause(5);
%%--- FURTHEST POINT DOWNSAMPLE THE CLOUD
%tic
kdtree = kdtree_build( cpts );
spls = zeros( 0, 3 );
corresp = zeros( length(cpts), 1 );
% mindst(i) is the min distance of pts(i) to the sample piont corresp(i) 
mindst = nan( length(cpts), 1 ); 


for k=1:length(cpts)
    if corresp(k)~=0 
        continue 
    end
    
    %--- query all the points for distances
    mindst(k) = inf; % make sure picked first
    
    %--- initialize the priority queue
    while ~all(corresp~=0) %~isempty( find(corresp==0, 1) )
        [maxValue, maxIdx] = max( mindst );
        if mindst(maxIdx) == 0
            break
        end

        % The radius of the farthest sampling sphere is determined based on linearity
        Rth = 1./(1+ exp(-(Linearity(maxIdx))));
        [nIdxs, nDsts] = kdtree_ball_query( kdtree,cpts(maxIdx,:), Rth*R);

        % if maxIdx and all its neighborhood has been marked, skip ahead
        if all( corresp(nIdxs) ~= 0 )
            mindst(maxIdx) = 0; 
            continue;
        end

        % create new node and update (closest) distances
        spls(end+1,:) = cpts(maxIdx,:); 
        for i=1:length(nIdxs)
            if mindst(nIdxs(i))>nDsts(i) || isnan(mindst(nIdxs(i)))
               mindst(nIdxs(i)) = nDsts(i);
               corresp(nIdxs(i)) = size(spls,1);
            end
        end

        if SHOW_SAMPLING_PROGRESS == true
            figure(1); plot3( cpts(maxIdx,1), cpts(maxIdx,2), cpts(maxIdx,3), '*g');
        end
    end
end
%toc
kdtree_delete( kdtree );

% if SHOW_RESULTS && ~SHOW_SAMPLING_PROGRESS
%     plot3( spls(:,1), spls(:,2), spls(:,3), '*g');
%     hold on;axis off; axis equal;set(gcf,'Renderer','OpenGL');
% end
end

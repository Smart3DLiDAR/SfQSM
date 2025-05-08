function [W, C, ptsC] = dbscan2105(P, E, minPts)

[~, Npts] = size(P);

%neighbours=dis(P',E);
MdlKDT = KDTreeSearcher(P');
neighbours= rangesearch(MdlKDT,P',E);

W={};
ptsC  = zeros(Npts,1);
C     = {};
Nc    = 0;               % Cluster counter.
Pvisit = zeros(Npts,1);  % Array to keep track of points that have been visited.

for n = 1:Npts
    if ~Pvisit(n)                            % If this point not visited yet
        Pvisit(n) = 1;                       % mark as visited
        
        neighbourPts=cell2mat(neighbours(n));
        
        if length(neighbourPts) <= minPts-1  % Not enough points to form a cluster
            ptsC(n) = 0;                    % Mark point n as noise.
            
        else                % Form a cluster...
            Nc = Nc + 1;    % Increment number of clusters and process
            % neighbourhood.
            
            C{Nc} = [n];    % Initialise cluster Nc with point n
            ptsC(n) = Nc;   % and mark point n as being a member of cluster Nc.
            
            ind = 1;        % Initialise index into neighbourPts array.
            
            % For each point P' in neighbourPts ...
            while ind <= length(neighbourPts)
                
                nb = neighbourPts(ind);
                
                if ~Pvisit(nb)        % If this neighbour has not been visited
                    Pvisit(nb) = 1;   % mark it as visited.
                    
                    % Find the neighbours of this neighbour and if it has
                    % enough neighbours add them to the neighbourPts list
                    %neighbourPtsP = regionQuery(P, nb, E);
                    %nei2= rangesearch(KDtree,P(:,nb)',E);
                    neighbourPtsP = cell2mat(neighbours(nb));
                    
                    if length(neighbourPtsP) >= minPts
                        neighbourPts = [neighbourPts  neighbourPtsP];
                        %neighbourPts = unique(neighbourPts);
                    end
                end
                
                % If this neighbour nb not yet a member of any cluster add it
                % to this cluster.
                if ~ptsC(nb)
                    C{Nc} = [C{Nc} nb];
                    ptsC(nb) = Nc;
                end
                
                ind = ind + 1;  % Increment neighbour point index and process
                % next neighbour
            end
        end
    end
end

for i=1:length(C)
    
    W(i) ={P(:,cell2mat(C(i)))'};
end
end % of dbscan
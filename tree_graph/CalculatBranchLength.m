function BranchL = CalculatBranchLength(spls,Pid,spls_adj,rootid,Branch)
% --input------------------------------------------------------------------ 
% spls : skeleton points
% spls_adj : Connection matrix of skeleton points
% Pid : the parent node of all skeleton points
% rootid : the root node
% Branch : branch
% --output-----------------------------------------------------------------
%  BranchL : the length of each branch at every level 

%% --------------------------------------------------------------------------------------------------
Adj = zeros(size(spls_adj,1),size(spls_adj,1));
for i = 1:size(spls_adj,1)
    for j = 1:size(spls_adj,2)       
        if spls_adj(i,j)==1
            dis = pdist2(spls(i,:),spls(j,:),'euclidean');
            Adj(i,j) = dis;
            Adj(j,i) = dis;
        end
    end
end
Gp = graph(Adj);
% [~,D] = shortestpathtree(Gp,rootid,[1:size(spls,1)],'OutputForm','cell');

BranchL = cell(length(Branch),1);
for i = 1:length(Branch)
    curlevel = Branch{i,1};
    if i ~=1
        tempbranchl = zeros(length(curlevel),1);
        for j = 1:length(curlevel)
            curbranch = curlevel{j};
            % Find the parent node of the current branch starting point in the upper level
            roi = Pid(curbranch(1));
            [~,dis] = shortestpath(Gp,curbranch(end),roi);
            tempbranchl(j,1) = dis;
        end
        BranchL{i} = tempbranchl;
    else 
        % Branch length of the first level: 
        % The length from the end point of the first level to the root node
        curlevel = cell2mat(curlevel);
        [~,dis] = shortestpath(Gp,rootid,curlevel(end));
        BranchL{i} = dis;
    end
end

end
function cursontreelength = findsubbranchlength(spls,finnalBranch,finnalGp,newPid)
%%----Calculate the length of the sub-branches-----------------------------
%---input -----------------------------------------------------------------
%  spls : Skeleton points
%  finnalBranch : Branch
%  finnalGp : Skeleton graph
%  newPid : The parent node of Skeleton points
%---output ----------------------------------------------------------------
%  cursontreelength : the length of the sub-branches
%% 
cursontreelength = zeros(1,size(spls,1));
for j = 1:size(spls,1)
    curid = j;
    for i = 1:length(finnalBranch)
        for k = 1:length(finnalBranch{i,1})
            curbranch = finnalBranch{i,1}{k,1};
            roi = intersect(curbranch,curid);
            if ~isempty(roi)
                [~,dis] = shortestpath(finnalGp,curid,curbranch(end));
            end
        end
    end
    if ~isempty(dis)
        cursontreelength(j) = dis;
    else
        cursontreelength(j) = 0;
    end
end

[~,col] = find(cursontreelength==0);
for i = 1:length(col)
    curnormal = col(i);
    pid = newPid(curnormal);
    cursontreelength(curnormal) = cursontreelength(pid);
end

end

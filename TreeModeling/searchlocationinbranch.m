function [curbranch,curlevel,lastevel,curlevelid] = searchlocationinbranch(curid,finnalBranch)
%%-------------------------------------------------------------------------------------------------
%%input 
%  curid           
%  finnalBranch    
%%output 
%  curbranch       
%  curlevel        
%  lastevel        
%  curlevelid      
%--------------------------------------------------------------------------------------------------
for i = 1:length(finnalBranch)
    for j = 1:length(finnalBranch{i,1})
        [row,~] = find(finnalBranch{i,1}{j,1} == curid);
        if ~isempty(row)
            curbranch = finnalBranch{i,1}{j,1};
            curlevelid = [i j];
            level = finnalBranch{i,1};
            curlevel = [];
            for k = 1:length(level)
                Temp = level{k,1};
                curlevel = [curlevel;Temp];
            end
        end
    end
end
clear i j 

if curlevelid(1) ~= 1
    last = finnalBranch{curlevelid(1)-1,1};
    lastevel = [];
    for j = 1:length(last)
        Temp = last{j,1};
        lastevel = [lastevel;Temp];
    end
else
    lastevel = [];
end
%curblength = finnalBranchL{curbranchid(1),1}(curbranchid(2),1);
%curblength = finnalD(curid);
end

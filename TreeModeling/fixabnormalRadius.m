function [newRadius,link,BranchRadius] = fixabnormalRadius(Radius,finnalBranch,link,rootid,newCid,newPid,subbranchlength,normallabel,TH)
%%
newRadius = zeros(length(Radius),3);
BranchRadius = {};
for i = 1:length(finnalBranch)
    curTH = cell2mat(TH{i,:});
    [row,~] = find(curTH(:,1)==1 & curTH(:,2)>=0);
    if ~isempty(row)
        curth(:,1) = mean(curTH(row,2),"omitnan");% Mean value
        curth(:,2) = min(curTH(row,2),[],'all',"omitnan");% Minimum value
    else
        curth = [1.5  1/2.49];
    end

    for j = 1:length(finnalBranch{i,1})
        curbranch = finnalBranch{i,1}{j,1};
        for k = 1:length(curbranch)
            id = curbranch(k);
            curson = newCid(id,:);
            curson(curson==0) = [];
            cursubbranchL = subbranchlength(id);
            
            if id ~= curbranch(1)
                pid = newPid(id,1);
                cursubbranchLP = subbranchlength(pid);
                
                son = newCid(pid,:);
                son(son==0) = [];
                if normallabel(id) == 0
                    % If the radius of this node is abnormal 
                    % it needs to be corrected
                    if length(son)==1
                        th = curth(1);% Mean value
                    else
                        th = curth(2);% Minimum value
                    end
                    tempth = cursubbranchL/cursubbranchLP;

                    tempR = newRadius(pid)*(tempth^th);
                    newRadius(id,:) = [tempR i th];
                    
                    [curfix1,~] = find(link(:,2)==id);
                    [curfix2,~] = find(link(:,1)==id);
                    link(curfix1,4:6) = [tempR 1 th];
                    link(curfix2,3) = tempR;

                    TH{i,1}{j,1}(k-1,:) = [1 th];
                else
                    [ind,~] = find(link(:,2)==id);
                    tempR = Radius(id);
                    newRadius(id,:) = [tempR i link(ind,6)];
                end
                
                if id ~= curbranch(end)
                    cid = curbranch(k+1);
                    [sonrow,~] = find(link(:,2) == cid);
                    Lc = subbranchlength(cid);
                    minR = tempR*((Lc/cursubbranchL)^2);
                    newth = (log(Radius(cid)/tempR))/(log(Lc/cursubbranchL));
                    if Radius(cid)<= tempR && Radius(cid) >= minR
                        link(sonrow,5:6) = [1 newth];
                        normallabel(cid) = 1;
                    else
                        link(sonrow,5) = 0;
                        normallabel(cid) = 0;
                    end
                end
            elseif id == curbranch(1) && id == rootid
                tempR = Radius(rootid);
                newRadius(id,:) = [tempR i 0];
            else
                pid = newPid(id,1);
                cursubbranchLP = subbranchlength(pid);
                tempth = cursubbranchL/cursubbranchLP;
                if Radius(id)<= newRadius(pid) 
                    [ind,~] = find(link(:,2)==id);
                    tempR = Radius(id);
                    newRadius(id,:) = [tempR i link(ind,6)];
                else    
                    tempR = newRadius(pid)*(tempth^curth(2));
                    newRadius(id,:) = [tempR i curth(2)];
                    [curfix1,~] = find(link(:,2)==id);
                    [curfix2,~] = find(link(:,1)==id);
                    link(curfix1,4:6) = [tempR 1 curth(2)];
                    link(curfix2,3) = tempR;

                    TH{i,1}{j,1}(1,:) = [1 curth(2)];
                end

                if id ~= curbranch(end)
                    cid = curbranch(k+1);
                    [sonrow,~] = find(link(:,2) == cid);
                    Lc = subbranchlength(cid);
                    minR = tempR*((Lc/cursubbranchL)^2);
                    newth = (log(Radius(cid)/tempR))/(log(Lc/cursubbranchL));
                    if Radius(cid)<= tempR && Radius(cid) >= minR
                        link(sonrow,5:6) = [1 newth];
                        normallabel(cid) = 1;
                    else
                        link(sonrow,5) = 0;
                        normallabel(cid) = 0;
                    end
                end
            end
            normallabel(id) = 1;
            if tempR<0.001
                tempR = 0.001;
            end
            BranchRadius{i,1}{j,1}(k,1) = tempR;
        end
    end
end

newRadius(newRadius(:,1)<0.001,:) = 0.001;
end

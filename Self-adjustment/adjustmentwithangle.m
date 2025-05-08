function [newspls2,Radius] = adjustmentwithangle(newspls,pts,rootid,nodelabel,Branch,Pid,Cid)
%%
Radius =  newspls.Radius;
nodeslabels1 = newspls.label;
newspls2 = newspls.spls;
Seg = newspls.Seg;
for i = 1:length(Branch)
    for j = 1:length(Branch{i,1})
        curbranch = Branch{i,1}{j,1};
        for k = 1:length(curbranch)
            id = curbranch(k);
            selfson = Cid(id,:);
            selfson(selfson==0) = [];
            if id ~= curbranch(1)
                pid = Pid(id,1);
                ppid = Pid(pid,1);
                son = Cid(pid,:);
                son(son==0) = [];
                
                if length(son)==1
                    % Only one child node
                    theta = 1/36;
                else
                    theta = 1/6;
                end
                
                curlabel = nodeslabels1(id);
                if curlabel ~=1 && nodelabel(id)<=5
                % Angle adjustment is only required if the point is not a reliable one 
                % and belongs to the first few important levels
                    adjustp = newspls2(id,1:3);
                    BeginP = newspls2(pid,1:3);
                    vert2 = (adjustp - BeginP)/norm(adjustp - BeginP);
                    if ppid ~=0
                        ppidP = newspls2(ppid,1:3);
                        vert1 = (BeginP - ppidP)/norm(BeginP - ppidP);
                        stopP = ppidP;
                    else
                        [reliabson,~]=find(nodeslabels1(selfson)==1);
                        if isempty(reliabson) && ppid==0
                            [row,~] = find(pts(:,3)<=(min(pts(:,3))+1.3));
                            neibor = pts(row,1:3);
                            [V,~,~] = pca(neibor, 'Economy', false); 
                            [~,maxid] = max(V(3,:));
                            vert1 = V(:,maxid)/norm(V(:,maxid));
                            vert1 = vert1';
                            stopP = BeginP + 1*vert1;
                        elseif isempty(reliabson) && ppid~=0
                            neibor = [Seg{id,1};Seg{pid,1}];
                            [V,~,~] = pca(neibor, 'Economy', false); 
                            [~,maxid] = max(V(3,:));
                            vert1 = V(:,maxid)/norm(V(:,maxid));
                            vert1 = vert1';
                            stopP = BeginP + 1*vert1;
                        elseif ~isempty(reliabson) && length(reliabson)==1
                            stopP = newspls2(selfson(reliabson),1:3);
                            vert1 = (stopP-BeginP)/norm(stopP-BeginP);
                        else
                            idex = selfson(nodeslabels1(selfson)==1);
                            v = zeros(length(idex),3);
                            for m = 1:length(idex)
                                tempstopP = newspls2(selfson(m),1:3);
                                v(m,:) = (tempstopP-BeginP)/norm(tempstopP-BeginP);
                            end
                            V = sum(v);
                            vert1 = V/norm(V);
                            stopP = BeginP + 1*vert1;
                        end
                    end
                    angle = acos(dot(vert1,vert2)/(norm(vert1)*norm(vert2)));
                    Th = angle/pi;
                    
                    [~,pedal] = footseeking(BeginP,stopP,adjustp);
                    if Th > theta && Th ~= 1
                        for l = 1:10
                            tempth = Th;
                            if tempth > theta  && tempth ~= 1
                                %If the threshold is not met, this point requires a secondary adjustment
                                adjustp = mean([pedal;adjustp]);
                                vert2 = (adjustp - BeginP)/norm(adjustp - BeginP);
                                if ppid ~=0
                                    ppidP = newspls2(ppid,1:3);
                                    vert1 = (BeginP - ppidP)/norm(BeginP - ppidP);
                                else
                                    [reliabson,~]=find(nodeslabels1(selfson)==1);
                                    if isempty(reliabson) && ppid==0
                                        [row,~] = find(pts(:,3)<=(min(pts(:,3))+1.3));
                                        neibor = pts(row,1:3);
                                        [V,~,S] = pca(neibor, 'Economy', false); 
                                        [~,maxid] = max(V(3,:));
                                        vert1 = V(:,maxid)/norm(V(:,maxid));
                                        vert1 = vert1';
                                    elseif isempty(reliabson) && ppid~=0
                                        neibor = [Seg{id,1};Seg{pid,1}];
                                        [V,~,S] = pca(neibor, 'Economy', false); 
                                        [~,maxid] = max(V(3,:));
                                        vert1 = V(:,maxid)/norm(V(:,maxid));
                                        vert1 = vert1';
                                    elseif ~isempty(reliabson) && length(reliabson)==1
                                        tempstopP = newspls2(selfson(reliabson),1:3);
                                        vert1 = (tempstopP-BeginP)/norm(tempstopP-BeginP);
                                    else
                                        idex = selfson(nodeslabels1(selfson)==1);
                                        v = zeros(length(idex),3);
                                        for m = 1:length(idex)
                                            tempstopP = newspls2(idex(m),1:3);
                                            v(m,:) = (tempstopP-BeginP)/norm(tempstopP-BeginP);
                                        end
                                        V = sum(v);
                                        vert1 = V/norm(V);
                                    end
                                end
                                angle = acos(dot(vert1,vert2)/(norm(vert1)*norm(vert2)));
                                Th = angle/pi;
                                if Th == 1
                                    break
                                end
                            else
                                break
                            end
                        end
                    end
                    curseg = Seg{id,1};
                    dis = zeros(size(curseg,1),1);
                    for m = 1:size(curseg,1)
                        dis(m) = pdist2(curseg(m,:),newspls2(id,:),'euclidean');
                    end
                    tempRadius = mean(dis);
                    
                    % Adjust the displacement before and after
                    gap = sqrt(sum((adjustp - newspls.spls(id,1:3)).^2));
                    if gap >= 0.5*tempRadius
                        adjustp = newspls.spls(id,1:3);
                    end
                    
                    newspls2(id,:) = adjustp;
                    
                    if size(unique([adjustp;newspls.spls(id,1:3)],'row'),1) ~=1
                        curseg = Seg{id,1};
                        dis = zeros(size(curseg,1),1);
                        for m = 1:size(curseg,1)
                            dis(m) = sqrt(sum((curseg(m,:) - adjustp).^2));
                        end
                        Radius(id) = mean(dis);
                    end
                    
                    nodeslabels1(id) = 1;
                end
            else%The current point is the first point on the branch
                if id ~= rootid && nodeslabels1(id)~=1 && nodelabel(id)<=5
                    pid = Pid(id,1);
                    ppid = Pid(pid,1);
                    theta = 1/6;
                    
                    adjustp = newspls2(id,1:3);
                    BeginP = newspls2(pid,1:3);
                    vert2 = (adjustp - BeginP)/norm(adjustp - BeginP);
                    if ppid ~=0
                        ppidP = newspls2(ppid,1:3);
                        vert1 = (BeginP - ppidP)/norm(BeginP - ppidP);
                        stopP = ppidP;
                    else
                        [reliabson,~]=find(nodeslabels1(selfson)==1);
                        if isempty(reliabson) && ppid==0
                            [row,~] = find(pts(:,3)<=(min(pts(:,3))+1.3));
                            neibor = pts(row,1:3);
                            [V,~,~] = pca(neibor, 'Economy', false); 
                            [~,maxid] = max(V(3,:));
                            vert1 = V(:,maxid)/norm(V(:,maxid));
                            vert1 = vert1';
                            stopP = BeginP + 1*vert1;
                        elseif isempty(reliabson) && ppid~=0
                            neibor = [Seg{id,1};Seg{pid,1}];
                            [V,~,~] = pca(neibor, 'Economy', false); 
                            [~,maxid] = max(V(3,:));
                            vert1 = V(:,maxid)/norm(V(:,maxid));
                            vert1 = vert1';
                            stopP = BeginP + 1*vert1;
                        elseif ~isempty(reliabson) && length(reliabson)==1
                            stopP = newspls2(selfson(reliabson),1:3);
                            vert1 = (stopP-BeginP)/norm(stopP-BeginP);
                        else
                            idex = selfson(nodeslabels1(selfson)==1);
                            v = zeros(length(idex),3);
                            for m = 1:length(idex)
                                tempstopP = newspls2(idex(m),1:3);
                                v(m,:) = (tempstopP-BeginP)/norm(tempstopP-BeginP);
                            end
                            V = sum(v);
                            vert1 = V/norm(V);
                            stopP = BeginP + 1*vert1;
                        end
                    end
                    angle = abs(acos(dot(vert1,vert2)/(norm(vert1)*norm(vert2))));
                    Th = angle/pi;
                    
                    [~,pedal] = footseeking(BeginP,stopP,adjustp);
                    if Th > theta && Th ~= 1
                        for l = 1:10
                            tempth = Th;
                            if tempth > theta && tempth ~= 1
                                adjustp = mean([pedal;adjustp]);
                                vert2 = (adjustp - BeginP)/norm(adjustp - BeginP);
                                if ppid ~=0
                                    ppidP = newspls2(ppid,1:3);
                                    vert1 = (BeginP - ppidP)/norm(BeginP - ppidP);
                                else
                                    [reliabson,~]=find(nodeslabels1(selfson)==1);
                                    if isempty(reliabson) && ppid==0
                                        [row,~] = find(pts(:,3)<=(min(pts(:,3))+1.3));
                                        neibor = pts(row,1:3);
                                        [V,~,~] = pca(neibor, 'Economy', false); 
                                        [~,maxid] = max(V(3,:));
                                        vert1 = V(:,maxid)/norm(V(:,maxid));
                                        vert1 = vert1';
                                    elseif isempty(reliabson) && ppid~=0
                                        neibor = [Seg{id,1};Seg{pid,1}];
                                        [V,~,~] = pca(neibor, 'Economy', false);
                                        [~,maxid] = max(V(3,:));
                                        vert1 = V(:,maxid)/norm(V(:,maxid));
                                        vert1 = vert1';
                                    elseif ~isempty(reliabson) && length(reliabson)==1
                                        tempstopP = newspls2(selfson(reliabson),1:3);
                                        vert1 = (tempstopP-BeginP)/norm(tempstopP-BeginP);
                                    else
                                        idex = selfson(nodeslabels1(selfson)==1);
                                        v = zeros(length(idex),3);
                                        for m = 1:length(idex)
                                            tempstopP = newspls2(selfson(idex(m)),1:3);
                                            v(m,:) = (tempstopP-BeginP)/norm(tempstopP-BeginP);
                                        end
                                        V = sum(v);
                                        vert1 = V/norm(V);
                                    end
                                end
                                angle = abs(acos(dot(vert1,vert2)/(norm(vert1)*norm(vert2))));
                                Th = angle/pi;
                                if Th ==1%If the lines are collinear, stop adjusting
                                    break
                                end
                            else
                                break
                            end
                        end
                    end
                    curseg = Seg{id,1};
                    dis = zeros(size(curseg,1),1);
                    for m = 1:size(curseg,1)
                        dis(m) = pdist2(curseg(m,:),newspls2(id,:),'euclidean');
                    end
                    tempRadius = mean(dis);
                    
                    gap = sqrt(sum((adjustp - newspls.spls(id,1:3)).^2));
                    if gap >= 0.5*tempRadius
                        adjustp = newspls.spls(id,1:3);
                    end
                    
                    newspls2(id,:) = adjustp;
                    
                    if size(unique([adjustp;newspls.spls(id,1:3)],'row'),1) ~=1
                        curseg = Seg{id,1};
                        dis = zeros(size(curseg,1),1);
                        for m = 1:size(curseg,1)
                            dis(m) = sqrt(sum((curseg(m,:) - adjustp).^2));
                        end
                        Radius(id) = mean(dis);
                    end
                    nodeslabels1(id) = 1;
                end
            end
        end
    end
end

end

function [fakebranch,spls,spls_adj] = remove_upper_nodes(pts,spls,spls_adj,SHOW_INNEREDGECOLLAPSE_PROGRESS)
[num, ~]=size(spls);
joints=zeros(0,1);
roots=zeros(0,1);
branches=zeros(0,1);
fakebranch=zeros(0,1);
for i=1:num   
    links = find( spls_adj(i,:)==1 );
    if length(links) == 1 % root
       roots(end+1)=i;
    end
    if length(links) == 2 % branch
        branches(end+1)=i;
    end
     if length(links) >2 % joints
        joints(end+1)=i;
    end
end

joints=joints';
roots=roots';
branches=branches';

for i=1:num   
    links = find( spls_adj(i,:)==1 );%row
    if length(links) == 0 % Ðé¼ÙÖ¦¸É/Ê÷Ò¶
       fakebranch(end+1)=i;
    elseif length(links) == 1
        [ind1,~] = find( roots(:,1)== i );
        [ind2,~] = find( roots(:,1)== links );
        if ~isempty(ind1) && ~isempty(ind2)
            fakebranch(end+1:end+2)=[i links];
        end
        clear ind1 ind2
    elseif length(links) == 2
        [ind1,~] = find( roots(:,1)== links(1));
        [ind2,~] = find( roots(:,1)== links(2));
        if ~isempty(ind1) && ~isempty(ind2)
            fakebranch(end+1:end+3)=[links(1) i links(2)];
        end
    end
end
clear ind

fakebranch = unique(fakebranch);
%for i = 1:length(fakebranch)
%    ind = fakebranch(i);
%    spls(ind,:)=NaN;
%    spls_adj(ind,:)=0;
%    spls_adj(:,ind)=0;
%    [row,~] = find(corresp(:)==ind);
%    cpts(row,:) = NaN;
%    pts(row,:) = NaN;
%    corresp(row,:) = NaN;
%end
%[m,~] = find(isnan(cpts)==1);
%cpts(m,:) = [];
%[n,~]=find(isnan(pts)==1);
%pts(n,:) = [];
%[n,~]=find(isnan(corresp)==1);
%corresp(n,:) = [];

spls(fakebranch(:),:) = NaN;
spls_adj(fakebranch(:),:) = 0;
spls_adj(:,fakebranch(:)) = 0;


if SHOW_INNEREDGECOLLAPSE_PROGRESS
    figure('Name','Remove upper nodes','NumberTitle','off');set(gcf,'color','white');
    subplot(1,2,1)
    scatter3(pts(:,1),pts(:,2),pts(:,3),30,'.');
    axis off; axis equal;
    title('Primitive tree')
    subplot(1,2,2)
    plot_skeleton(spls, spls_adj);
    axis off; axis equal;
    title('skeleton without fake skleton points')
end

end
function [newsplsssss,newspls_adj,newSeg,Radius,Graph] = Reconstructedtopology(newsplsssss,newspls_adj,newSeg,Radius)
%%
Graph = zeros(0,2);
for i=1:size(newspls_adj,1)
    for j=i+1:size(newspls_adj,2)
        if( newspls_adj(i,j)==1 )
            Graph(end+1,:) = [i, j];
        end
    end
end
% Find the positions of the deleted skeleton points 
tmp = find(isnan(newsplsssss(:,1)))';
Tmp = tmp(length(tmp):-1:1); 
for i=Tmp
    Graph(Graph>i) = Graph(Graph>i) - 1;
end
newsplsssss(tmp,:) = [];
newSeg(tmp,:) = [];
Radius(tmp,:) = [];
newspls_adj = zeros(size(newsplsssss,1),size(newsplsssss,1));
for a=1:size(Graph,1)
    i = Graph(a,1); j = Graph(a,2);
    newspls_adj(i,j) = 1;
    newspls_adj(j,i) = 1;
end

end

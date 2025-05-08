function Seg = Simplifyskeleton(spls1,spls,newdelateb,Seg)
%%------------------------------------------------------------------------------------------------
%%input 
%  spls1 : The skeleton points after pruning
%  spls  : The skeleton points before pruning
%  newdelateb  :    The skeleton points that need pruning
%  Seg  :  Contraction point set           
%%output 
%  Seg  : Contraction point set  after pruning   
%--------------------------------------------------------------------------------------------------

for i = 1:length(newdelateb)
    curd = newdelateb(i);
    nei = knnsearch(spls1,spls(curd,:),'k',2);
    [~,pedal] = footseeking(spls1(nei(1),:),spls1(nei(2),:),spls(curd,:));
    dis = zeros(length(nei),1);
    for j = 1:length(nei)
        dis(j) = pdist2(pedal,spls1(nei(j),:));
    end
    [~,minid] =min(dis);
    Seg{nei(minid),:} = union(Seg{nei(minid),:},Seg{curd,:},'row');
end

end

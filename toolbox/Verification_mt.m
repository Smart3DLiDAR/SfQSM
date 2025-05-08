function [Precison,Recall,Micro_F1,Macro_F1,OA,Precisons,Recalls,F1s,RegionsCorresp] = Verification_mt(pts, trueRegions, autoRegions,IdName,isWrite)
%VERIFICATION_MT 此处显示有关此函数的摘要
%   此处显示详细说明
 t_size=length(trueRegions)-1;
 a_size=length(autoRegions)-1;
 CorrespIndex=zeros(t_size,1);
 CorrespMax=zeros(t_size,1);
 for i=1: t_size
   t_indices=trueRegions{i+1};
   interNum=zeros(a_size,1);
   for j=1:a_size
     a_indices= autoRegions{j+1};  
     inter=intersect(t_indices,a_indices);
     interNum(j)=length(inter);
   end
   [maxNum,idx]=sort(interNum);
   CorrespIndex(i)=idx(end);
   CorrespMax(i)=maxNum(end);
 end
 %%%%%%%%%如果出现多个trueRegions对应一个autoRegions的情况
 %%%%%%%%%选取交集最大的作为对应，剩下的认为是FN
 for i=1:t_size
      if(CorrespIndex(i)==-1)
          continue;
      end
   for j=i+1:t_size 
      if(CorrespIndex(j)==-1)
          continue;
      end
      if(CorrespIndex(i)==CorrespIndex(j))
         if(CorrespMax(i)>CorrespMax(j))
           CorrespIndex(j)=-1;
           CorrespMax(j)=-1;
         else
           CorrespIndex(i)=-1;  
          CorrespMax(i)=-1;
         end
      end
   end
 end
 RegionsCorresp=CorrespIndex;
%  autoLabels=zeros(ptnum,1);
%  trueLabels=zeros(ptnum,1);
% for i=1:length(autoRegions)
%    label=find(CorrespIndex==i);
%    if(isempty(label))
%        label=-1;
%    end
%    ids=autoRegions{i};
%    autoLabels(ids)=label;
% end
% for i=1:length(trueRegions)
%    ids=trueRegions{i};
%    trueLabels(ids)=i;  
% end
  ptnum=0;
precison=zeros(t_size,1);recall=zeros(t_size,1);f1_score=zeros(t_size,1);
for i=2:length(trueRegions)
  tr=trueRegions{i};
  ptnum=ptnum+length(tr);
  j=CorrespIndex(i-1,1);
  if(j==-1)
     precison(i)=0; recall(i)=0;f1_score(i)=0;
     continue;
  end
  ar=autoRegions{j+1};
  precison(i)=CorrespMax(i-1)/length(ar);
  recall(i)=CorrespMax(i-1)/length(tr);
  if(recall(i)==0||precison(i)==0)
      f1_score(i)=0;
  else
  f1_score(i)=2*precison(i)*recall(i)/(recall(i)+precison(i));
  end
end
Precison=sum(precison)/(length(precison)-1);
Recall=sum(recall)/(length(recall)-1);
Micro_F1=2*(Precison*Recall)/(Precison+Recall);
Macro_F1=sum(f1_score)/(length(f1_score)-1);

Precisons=precison;
Recalls=recall;
F1s=f1_score;

  OA=sum(CorrespMax)/ptnum;
  txtName=[IdName '_p.txt']; 
  fid=fopen(txtName,'w');
  fprintf(fid,'%d Precison %f\r\n',0,Precison);
  fprintf(fid,'%d Recall %f\r\n',0,Recall);
  fprintf(fid,'%d Micro_F1 %f\r\n',0,Micro_F1);
  fprintf(fid,'%d Macro_F1 %f\r\n',0,Macro_F1);
  fprintf(fid,'%d OA %f\r\n',0,OA);
  for i=2:length(precison)
  fprintf(fid,'%d Precison %f\r\n',i-1,precison(i));
  fprintf(fid,'%d Recall %f\r\n',i-1,recall(i));
  fprintf(fid,'%d f1score %f\r\n',i-1,f1_score(i));    
  end
  fclose(fid);


 if(true)
    figure('Name','True segment','NumberTitle','off');set(gcf,'color','white');movegui('southwest');
    hold on;
    color=MyGS.MYCOLOR;
    scatter3(pts(trueRegions{1},1),pts(trueRegions{1},2),pts(trueRegions{1},3),5,[0 0 0], 'filled');
    for i=2:length(trueRegions)
    ids=trueRegions{i};
    scatter3(pts(ids,1),pts(ids,2),pts(ids,3),5,color(i-1,:), 'filled');
    hold on;
    end
    axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(180,-65);view3d ZOOM;
%     if(isWrite)
%     saveas(gcf,['.\result\',IdName,'_T'],'png');
%     end
 end
  
 if(true)
    figure('Name','True ear','NumberTitle','off');set(gcf,'color','white');movegui('southwest');
    hold on;
    color=MyGS.MYCOLOR;
 %   scatter3(pts(trueRegions{1},1),pts(trueRegions{1},2),pts(trueRegions{1},3),5,[0 0 0], 'filled');
    for i=2:length(trueRegions)
    ids=trueRegions{i};
    scatter3(pts(ids,1),pts(ids,2),pts(ids,3),5,color(i-1,:), 'filled');
    hold on;
    end
    axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(180,-65);view3d ZOOM;
%     if(isWrite)
%     saveas(gcf,['.\result\',IdName,'_T'],'png');
%     end
 end
 
 if(true)
    figure('Name','auto segment','NumberTitle','off');set(gcf,'color','white');movegui('southwest');
    hold on;
    color=MyGS.MYCOLOR; 
    scatter3(pts(autoRegions{1},1),pts(autoRegions{1},2),pts(autoRegions{1},3),5,[0 0 0], 'filled');
    for i=2:length(autoRegions)
    region=autoRegions{i};
    trId= find(CorrespIndex(:,1)+1==i);
    if(isempty(trId))
        trId=length(color)-1;
    end
    trId=trId(1);
    scatter3(pts(region,1),pts(region,2),pts(region,3),5,color(trId,:), 'filled');
    hold on;
    end
     
    
    axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(180,-65);view3d ZOOM;
%     if(isWrite)
%     saveas(gcf,['.\result\',IdName,'_A'],'png');
%     end
end
if(true)
    figure('Name','auto ear','NumberTitle','off');set(gcf,'color','white');movegui('southwest');
    hold on;
    color=MyGS.MYCOLOR; 
   % scatter3(pts(autoRegions{1},1),pts(autoRegions{1},2),pts(autoRegions{1},3),5,[0 0 0], 'filled');
    for i=2:length(autoRegions)
    region=autoRegions{i};
    trId= find(CorrespIndex(:,1)+1==i);
    if(isempty(trId))
        trId=length(color)-1;
    end
    trId=trId(1);
    scatter3(pts(region,1),pts(region,2),pts(region,3),5,color(trId,:), 'filled');
    hold on;
    end
     
    
    axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(180,-65);view3d ZOOM;
%     if(isWrite)
%     saveas(gcf,['.\result\',IdName,'_A'],'png');
%     end
end
 
 
 
end


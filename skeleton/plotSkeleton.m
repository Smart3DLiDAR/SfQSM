function plotSkeleton(pts, A, options)
%根据骨架点以及其邻接矩阵对其进行连接

if nargin<3
    options.sizep=400;
end
sizep = getoptions(options, 'sizep', 0.5);% 200
sizee = getoptions(options, 'sizee', 1); %2.5
%colorp = getoptions(options, 'colorp', [1, 0, 0]);%[1, .75, .79 ]
colore = getoptions(options, 'colore', [1, 0, 0]);%
colorp = getoptions(options, 'colorp', [1, 0, 0]);%[1, .75, .79 ]
%Color = getoptions(options, 'color'); 
%colorp = Color(50,:);
%colore = Color(200,:);%Color
scatter3(pts(:,1),pts(:,2),pts(:,3),sizep,'.','MarkerEdgeColor', colorp);  hold on;
plotConnectivity(pts, A, sizee, colore);


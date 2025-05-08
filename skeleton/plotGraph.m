function plotGraph(pts, GP, options)
%���ݹǼܵ��Լ����ڽӾ�������������

if nargin<3
    options.sizep=400;
end
sizep = getoptions(options, 'sizep', 200); 
sizee = getoptions(options, 'sizee', 2.5); 
colorp = getoptions(options, 'colorp', [1, .75, .79]); 
colore = getoptions(options, 'colore', [1, .0, .0]); 

scatter3(pts(:,1),pts(:,2),pts(:,3),sizep,'.','MarkerEdgeColor', colorp);  hold on;
plotedge(pts, GP, sizee, colore);
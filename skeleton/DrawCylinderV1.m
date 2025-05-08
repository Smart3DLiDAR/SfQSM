function DrawCylinderV1(point,a,h,r,color)
xlabel("X��")
ylabel("Y��")
zlabel("Z��")
%quiver3(point(1),point(2),point(3),a(1),a(2),a(3));
%����Բ����
[x,y,z]=cylinder(r,100);
z=z*h;      %�߶�
z=z-h/2;    %���������������ƶ�
x = x+point(1);     %����ƽ��,���ĵ�
y = y+point(2);
z = z+point(3);
Cylinder=mesh(z);
%����Բ������ת�ĽǶ�
V=[0 0 1];
angle=acos(dot(V,a)/(norm(V)*norm(a)))*180/pi;
axis_rot=cross(V,a);
%��Բ������ת����������
if angle~=0
    rotate(Cylinder,axis_rot,angle,point)
end
% ����Բ�������ɫ
set(Cylinder,'FaceColor',color,...
    'EdgeColor','none',...        %���ñߵ���ɫ
    'FaceAlpha',1)   
camlight
lighting gouraud
end

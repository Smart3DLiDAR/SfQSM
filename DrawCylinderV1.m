function DrawCylinderV1(point,a,h,r,color)
xlabel("X轴")
ylabel("Y轴")
zlabel("Z轴")
%quiver3(point(1),point(2),point(3),a(1),a(2),a(3));
%绘制圆柱体
[x,y,z]=cylinder(r,100);
z=z*h;      %尺度
z=z-h/2;    %将点云整体往下移动
x = x+point(1);     %进行平移,中心点
y = y+point(2);
z = z+point(3);
Cylinder=mesh(z);
%计算圆柱体旋转的角度
V=[0 0 1];
angle=acos(dot(V,a)/(norm(V)*norm(a)))*180/pi;
axis_rot=cross(V,a);
%将圆柱体旋转到期望方向
if angle~=0
    rotate(Cylinder,axis_rot,angle,point)
end
% 设置圆柱体的颜色
set(Cylinder,'FaceColor',color,...
    'EdgeColor','none',...        %禁用边的颜色
    'FaceAlpha',1)   
camlight
lighting gouraud
end

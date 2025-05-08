function [xc,yc,R]=circlefit(points)
% CIRCLEFIT fits a circle in x,y plane
% (x-a)^2 + (y-b)^2 = R^2
% x^2+y^2+a(1)*x+a(2)*y+a(3)=0

x=points(:,1);
y=points(:,2);

n=length(x);
xx=x.*x;
yy=y.*y;
xy=x.*y;

A=[sum(x) sum(y) n;sum(xy) sum(yy) sum(y);sum(xx) sum(xy) sum(x)];
B=[-sum(xx+yy);-sum(xx.*y+yy.*y);-sum(xx.*x+xy.*y)];
a=A\B;
xc = -0.5*a(1);
yc = -0.5*a(2);
R = sqrt(-(a(3)-xc^2-yc^2));

th = 0:pi/50:2*pi;
cx=R*cos(th)+xc;
cy=R*sin(th)+yc;

%figure
%scatter(x, y, '.b');hold on
%scatter(xc,yc,50, '.r');hold on
%scatter(curspls(:,1),curspls(:,2),50, '.k');hold on
%plot(cx, cy, 'Color', 'g');
%axis equal

end


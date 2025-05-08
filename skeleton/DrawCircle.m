function DrawCircle(center,r,z)
%ªÊ÷∆‘≤
theta = 0:pi/100:2*pi;
x=r*cos(theta)+center(1); 
y=r*sin(theta)+center(2);
tz = zeros(length(x),1);
tz(:) = z;
plot3(x,y,tz,'g-','LineWidth',2)
end

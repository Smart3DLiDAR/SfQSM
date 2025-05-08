function plotcircle(x,y,r)
    theta=0:0.01:2*pi+0.01;  
    Circle1=x+r*cos(theta);  
    Circle2=y+r*sin(theta);   
    plot(Circle1,Circle2,'g','linewidth',1);  
    axis equal
end


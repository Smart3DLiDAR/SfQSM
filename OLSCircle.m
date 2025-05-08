function r = OLSCircle(data,center)
%OLSCIRCLE 最小二乘法

df = data - center;
b = sqrt(sum(df.^2,2));     %这是所有点到圆心的距离
r = sum(b) / length(b);     %按照最小二乘法理论可以推出r等于其均值

end

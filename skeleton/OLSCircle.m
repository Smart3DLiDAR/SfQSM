function r = OLSCircle(data,center)
%OLSCIRCLE ��С���˷�

df = data - center;
b = sqrt(sum(df.^2,2));     %�������е㵽Բ�ĵľ���
r = sum(b) / length(b);     %������С���˷����ۿ����Ƴ�r�������ֵ

end

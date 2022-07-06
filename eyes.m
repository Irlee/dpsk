function eyes(y,nd,num)
%eyes 绘制眼图
%   y:低通滤波后的序列
%   nd:一个码元的采样点数
%   num:眼睛数量
L = length(y);
yls = y(2*nd+1:L);
av = (max(yls)+min(yls))/2;
av = abs(av);
if(av<10e-2)
    av=10e-4;
end
delta = (max(y)-min(y))/nd;
k = find(yls>av & yls<av+delta/5, 1);
if(isempty(k))
    disp('Warning:The input vector length is too short\r\n')
    k=0;
end
start = 2*nd+k;
gap = nd*num;
ii = start:gap:L-gap;
if(isempty(ii))
    disp('Warning:The input vector length is too short\r\n')
    start = 1;
end   
hold on
for ii = start:gap:L-gap
    plot(y(ii+1:ii+gap))
end
hold off

end

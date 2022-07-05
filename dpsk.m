clc
clear
close all
%% 采样率设置
Fs = 10e6;          %10MHz采样率
T = 1/Fs;           %采样周期
L = 0.1e6;          %采样点数
t = (1:L)*T;        %采样的时间窗：10ms
SNR = 10;           %分贝
%% 载波生成
f1 = 100e3;
f3 = 200e3;

fc1 = sin(2*pi*f1*t);
fc2 = sin(2*pi*f1*t + pi);
fc3 = sin(2*pi*f3*t);
%% 基带数据数据产生
f2 = 10e3;   T2 = 1/f2;     %10kbps信息速率
N = round(max(t)*f2);       %原始基带数据的长度
nd = length(t)/N;           %1bit对应的采样点数
P_data=randi([0 1],1,N);
%% 差分编码
diff_code = false(1,N+1);
for ii = 1:N
    if(P_data(ii)==1)
        diff_code(ii+1) = ~diff_code(ii);
    else
        diff_code(ii+1) = diff_code(ii);
    end
end
diff_code = diff_code(2:N+1);
%% manchester编码
manchester_code = zeros(1,2*N);
signal_Pdata = zeros(1,L);          %扩充原始信号
signal_diff = zeros(1,L);           %扩充差分信号
for ii = 1:N
    jj = 2*ii-1;
    zz = (ii-1)*nd + 1;
    if(diff_code(ii)==1)
        manchester_code(jj:jj+1) = [-1, 1];
    else
        manchester_code(jj:jj+1) = [1, -1];
    end
    signal_Pdata(zz:zz+nd-1) = P_data(ii);
    signal_diff(zz:zz+nd-1) = diff_code(ii);
end
signal_2polar = zeros(1,L);        %扩充基带manchester信号
ns = nd/2;
for ii = 1:2*N
    jj=(ii-1)*ns+1;
    signal_2polar(jj:jj+ns-1) = manchester_code(ii);
end
signal_spolar = (signal_2polar+1)/2;    %单极性信号
%% 调制
ask_signal = signal_spolar.*fc1;
fsk_signal = zeros(1,L);
for ii = 1:L
    if(signal_2polar(ii)==1)
        fsk_signal(ii) = fc1(ii);
    else
        fsk_signal(ii) = fc3(ii);
    end
end
% for ii = 1:L
%     if(signal_2polar(ii)==1)
%         dpsk_signal(ii) = fc1(ii);
%     else
%         dpsk_signal(ii) = fc2(ii);
%     end
% end
dpsk_signal = signal_2polar.*fc1;
%% 绘图
%10个周期的载波
figure(1)
subplot(3,1,1)
plot(t(1:1000), fc1(1:1000))
title('100kHz同相载波')

subplot(3,1,2)
plot(t(1:1000), fc2(1:1000))
title('反相载波')

subplot(3,1,3)
plot(t(1:1000), fc3(1:1000))
title('200kHz载波')

%%基带信号和编码信号
figure(2)
subplot(3,1,1)
x = 1:N;
stairs(x,P_data(1:N))
title('原始二进制信号')
ylim([-0.2 1.2])
subplot(3,1,2)
x = 1:N;
stairs(x,diff_code(1:N))
title('差分编码信号')
ylim([-0.2 1.2])
subplot(3,1,3)
x = 1:2*N;
stairs(x,manchester_code(1:2*N))
title('Manchester信号')
ylim([-1.2 1.2])
%%已调信号
figure(3)
x = t(1:5000);
Nx = length(x);
plot(x, ask_signal(1:Nx))
hold on
plot(x, fsk_signal(1:Nx)-3)
plot(x, dpsk_signal(1:Nx)-6)
plot(x, signal_spolar(1:Nx)-9)
ylim([-10,3])
legend({'2ask','2fsk','2dpsk','调制信号'},...
    'Location','best','NumColumns',2)
hold off

%% 加噪
% r1 = ask_signal + 0.1*randn(size(t));
% r2 = fsk_signal + 0.1*randn(size(t));
% r3 = dpsk_signal + 0.1*randn(size(t));
r1 = awgn(ask_signal,SNR);
r2 = awgn(fsk_signal,SNR);
r3 = awgn(dpsk_signal,SNR);
%% 2ask解调
hb = bandpass100k;
hl = lowpass20k;
hl10k = lowpass10k;
y1 = filter(hb,r1);

fc4 = sin(2*pi*f1*t + 3*pi/2);         %100kHz相干载波
y2 = filter(hl, y1.*fc4);
figure(4)
plot(t,y1)
hold on
plot(t,2*y2-3)
plot(t,signal_spolar-5)
ylim([-5.5 2.5])
title('2ask解调')
legend({'r1经带通滤波', 'r1解调后经低通滤波',...
    '差分Manchester信号'}, 'Location','northwest',...
    'NumColumns',2)
hold off
%% 2fsk
hb2 = bandpass200k;
y1 = filter(hb,r2);
fc5 = sin(2*pi*f3*t + 3*pi/4);  %200kHz相干载波
y2 = filter(hl, y1.*fc4);
y3 = filter(hb2, r2);
y4 = filter(hl, y3.*fc5);
figure(5)
plot(t,y1)
hold on
plot(t,2*y2-3)
plot(t,y3-6)
plot(t,2*y4-9)
plot(t,signal_spolar-12)
ylim([-12.5 3.5])
title('2fsk解调')
legend({'r2经带通1','低通1',...
    'r2经带通2','低通2','差分Manchester信号'},...
    'Location','northwest','NumColumns',3)
hold off

%% 2dpsk非相干解调
y1 = filter(hb,r3);
% y1 = dpsk_signal;
y1_delay = [zeros(1,nd), y1(1:end-nd)];
y3 = y1.*y1_delay;
y2 = filter(hl10k, y3);

figure(6)
plot(t,y1)
hold on
plot(t,y1_delay-3);
plot(t,y3-6)
plot(t,y2-9)
plot(t,signal_Pdata-12)
ylim([-12.5 3.5])
title('2dpsk非相干解调')
legend({'r3经带通滤波','r3延时一个码元',...
    '非相干解调','10k低通滤波','原始基带信号'},...
    'Location','best','NumColumns',3)
hold off
%注：2dpsk信号差分相干解调后直接得到原始基带信号
%   对差分Manchester信号解调时也可以直接得到原始信号
%   不需要码反变换，低通滤波得到的信号相位与原始信号反相
%   这里的差分编码规则是：1-相位跳变， 0-保持
%% 眼图
figure(7)
eyes(y2,nd,2)
title('2dpsk解调的眼图')

%% 频谱
figure(8)
subplot(3,1,1)
f = Fs*(0:(L/2))/L;
Y = fft(y1);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
plot(f,P1) 
title('r3经过带通滤波后的y1频谱图（幅度）')
xlabel('f (Hz)')
ylabel('y1')
xlim([0 350e3])

subplot(3,1,2)
Y = fft(y3);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
plot(f,P1) 
title('y1相乘后的频谱图（幅度）')
xlabel('f (Hz)')
ylabel('y3=y1*y1_delay')
xlim([0 350e3])

subplot(3,1,3)
Y = fft(y2);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
plot(f,P1) 
title('y3低通滤波后的频谱图（幅度）')
xlabel('f (Hz)')
ylabel('y2')
xlim([0 350e3])

clc
clear
close all

%% 载波生成
f1 = 100e3; T1 = 1/f1;
f3 = 300e3;
nc = 50;
TN = 10000;
dt = T1/nc;
t = dt:dt:TN*T1;
fc1 = sin(2*pi*f1*t);
fc2 = sin(2*pi*f1*t + pi);
fc3 = sin(2*pi*f3*t);
%% 基带数据数据产生
f2 = 10e3;   T2 = 1/f2;
N = max(t)/T2;
nd = length(t)/N;
P_data=randi([0 1],1,N);
%% 差分编码
diff_code = zeros(1,N);
for ii = 2:N
    if(P_data(ii)==P_data(ii-1))
        diff_code(ii-1)=1;
    else
        diff_code(ii-1)=0;
    end
end
%% manchester编码
manchester_code = zeros(1,2*N);
for ii = 1:N
    jj = 2*ii-1;
    if(diff_code(ii)==1)
        manchester_code(jj:jj+1) = [-1, 1];
    else
        manchester_code(jj:jj+1) = [1, -1];
    end
end
signal_2polar = zeros(1,nc*TN);        %扩充信号
ns = nd/2;
for ii = 1:2*N
    jj=(ii-1)*ns+1;
    signal_2polar(jj:jj+ns-1) = manchester_code(ii);
end
signal_spolar = (signal_2polar+1)/2;
%% 调制
ask_signal = signal_spolar.*fc1;
fsk_signal = zeros(1,nc*TN);
for ii = 1:nc*TN
    if(signal_2polar(ii)==1)
        fsk_signal(ii) = fc1(ii);
    else
        fsk_signal(ii) = fc3(ii);
    end
end
dpsk_signal = zeros(1,nc*TN);
for ii = 1:nc*TN
    if(signal_2polar(ii)==1)
        dpsk_signal(ii) = fc1(ii);
    else
        dpsk_signal(ii) = fc2(ii);
    end
end
%% 绘图
%%10个周期的载波
% figure(1)
% subplot(2,1,1)
% plot(t(1:10*nc), fc1(1:10*nc));
% title('同相载波');
% xlim([0 T1*10]);
% subplot(2,1,2)
% plot(t(1:10*nc), fc2(1:10*nc));
% title('反相载波');
% xlim([0 T1*10]);
% %%基带信号和编码信号
% figure(2)
% subplot(3,1,1)
% x = 1:100;
% stairs(x,P_data(1:100));
% title('原始二进制信号');
% ylim([-0.2 1.2]);
% subplot(3,1,2)
% x = 1:100;
% stairs(x,diff_code(1:100));
% title('差分编码信号');
% ylim([-0.2 1.2]);
% subplot(3,1,3)
% x = 1:200;
% stairs(x,manchester_code(1:200));
% title('Manchester信号');
% ylim([-1.2 1.2]);
% %%2ask已调信号
% figure(3)
% subplot(2,1,1);
% plot(t, ask_signal);
% hold on
% plot(t, signal_spolar-2.5);
% hold off
% title('2ask已调信号');
% xlim([0 T1*30]);
% ylim([-3 1.2]);
% subplot(2,1,2);
% plot(t,fc1);
% title('载波信号');
% xlim([0 T1*30])
% %%2fsk已调信号
% figure(4)
% subplot(3,1,1);
% plot(t, fsk_signal);
% xlim([0 T1*30])
% title('2fsk已调信号');
% subplot(3,1,2);
% plot(t,fc1);
% xlim([0 T1*30])
% title('载波信号1');
% subplot(3,1,3);
% plot(t,fc3);
% title('载波信号3');
% xlim([0 T1*30])
% %%2dpsk已调信号
% figure(5)
% subplot(3,1,1);
% plot(t, dpsk_signal);
% xlim([0 T1*30])
% title('2dpsk已调信号');
% subplot(3,1,2);
% plot(t,fc1);
% xlim([0 T1*30])
% title('载波信号1');
% subplot(3,1,3);
% plot(t,fc2);
% title('载波信号2');
% xlim([0 T1*30])




%% 解调
h = bandpass100k;
y1 = filter(h,ask_signal);
figure(6);
plot(t,y1*10);
hold on
plot(t,ask_signal-2);
xlim([0 T1*50]);
hold off

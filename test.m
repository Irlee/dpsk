clc
clear
close all
%% Noisy Signal
% Use Fourier transforms to find the frequency components of a signal buried in noise.
% Specify the parameters of a signal with a sampling frequency of 1 kHz and a signal duration of 1.5 seconds.
Fs = 10e6;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1e5;             % Length of signal
t = (0:L-1)*T;        % Time vector

%% Form a signal containing a 50 Hz sinusoid of amplitude 0.7 and a 120 Hz sinusoid of amplitude 1.
S = sin(2*pi*20e3*t) + sin(2*pi*100e3*t) + sin(2*pi*110e3*t) + sin(2*pi*120e3*t) + sin(2*pi*150e3*t);
S1 = square(2*pi*100e3*t);

%% Corrupt the signal with zero-mean white noise with a variance of 4.
X = S + 10*randn(size(t));

%% Plot the noisy signal in the time domain. It is difficult to identify the frequency components by looking at the signal X(t). 
figure(1)
subplot(3,1,1);
plot(1000*t(1:1000),X(1:1000))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')
subplot(3,1,2)
plot(1000*t(1:1000),S(1:1000))
title('No noise signal')
xlabel('t (milliseconds)')
ylabel('X(t)')
subplot(3,1,3)
plot(1000*t(1:1000),S1(1:1000))
title('square')
xlabel('t (milliseconds)')
ylabel('S1(t)')

%% Compute the Fourier transform of the signal. 
Y = fft(X);

%% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

%% Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.
f = Fs*(0:(L/2))/L;
figure(2)
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Now, take the Fourier transform of the original, uncorrupted signal and retrieve the exact amplitudes, 0.7 and 1.0.
Y = fft(S);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure(3)
subplot(2,1,1)
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of S(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

subplot(2,1,2)
Y = fft(S1);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(f,P1) 
title('Single-Sided Amplitude Spectrum of S1(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% 白噪声滤波
h = bandpass100k;
y = filter(h,(X-S)*100);
Fy = fft(y);
P2 = abs(Fy/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure(4)
subplot(2,1,1)
plot(f,P1) 
title('白噪声滤波后输出频谱')
xlabel('f (Hz)')
ylabel('|P1(f)|')

subplot(2,1,2)
Fy = fft((X-S)*100);
P2 = abs(Fy/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
plot(f,P1) 
title('白噪声频谱')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% 方波通过滤波器的波形
figure(5)
subplot(3,1,1)
plot(1000*t(1:10000),S1(1:10000))
title('square')
xlabel('t (milliseconds)')
ylabel('S1(t)')

y = filter(h,S1);
Fy = fft(y);
P2 = abs(Fy/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

subplot(3,1,2)
plot(1000*t(1:10000),y(1:10000))
title('方波滤波后')
xlabel('t (milliseconds)')
ylabel('S1(t)')

subplot(3,1,3)
plot(f,P1) 
title('方波通过滤波器输出频谱')
xlabel('f (Hz)')
ylabel('|P1(f)|')

figure(6)
[p f]= pspectrum(S);
plot(f,p)

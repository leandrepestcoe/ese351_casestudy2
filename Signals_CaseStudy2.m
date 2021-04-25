% Case Study 2
% Leandre Pestcoe and Julianne Wegmann
% ESE 351: Signals and Systems
% Created on: 4/20/21, Last Edited on: 4/23/21

%% Define Pulse Shape p(t)
Ts = 0.1; %symbol period (rate 1/Ts)
dt = 0.01; %sample period
w = 5*Ts; %width
t = -w:dt:w; %time vector
fs = 1/dt; %sample frequency

%sinc
sinc_p_t = sinc(t/Ts); %define sinc

%triangular pulse
triang_p_t = tripuls(t,w*2);

%decide which signal to use for simulation
p_t = sinc_p_t;

figure, plot(t,p_t), grid on;
xlabel('time (s)'), ylabel('p(t)'), title('Truncated Signal p(t)')

%% Nyquist Criteria

Nfft = 1024; %length of fft
f1 = (-fs/2:fs/Nfft:(fs/2));
f2 = (fs-2*pi/Ts:fs/Nfft:fs+(2*pi/Ts));
fft_signal = abs(fft(p_t,Nfft));
fft_signal1 = [fft_signal((length(fft_signal)/2)+1:end), fft_signal(1:(length(fft_signal)/2))];

figure
plot(f1(1:length(fft_signal1)),fft_signal1), grid on;
hold on
plot(f2(1:length(fft_signal1)),fft_signal1), grid on;
xlabel('frequency (Hz)'), ylabel('|Y(j\omega)|');

%% Random Bit Message/Noise-free PAM Signal y(t)

%define bit messages
N = 10; %number of bits
bits1 = 2*((rand(1,N)>0.5)-0.5);
bits2 = 2*((rand(1,N)>0.5)-0.5);
bits3 = 2*((rand(1,N)>0.5)-0.5);

%spacing out bits
xn1 = zeros(1,N*(Ts/dt));
xn2 = zeros(1,N*(Ts/dt));
xn3 = zeros(1,N*(Ts/dt));

for i=1:length(bits1)
    xn1((i-1)*(Ts/dt)+1)=bits1(i);
end
for i=1:length(bits2)
    xn2((i-1)*(Ts/dt)+1)=bits2(i);
end
for i=1:length(bits3)
    xn3((i-1)*(Ts/dt)+1)=bits3(i);
end

%compute individual CT transmitted signals
y_t1 = conv(xn1,p_t);
y_t2 = conv(xn2,p_t);
y_t3 = conv(xn3,p_t);

t = (0:length(y_t1)-1)*dt; %define new time vector

%shift xn so lines up with y(t) on graph
x_new1 = zeros(1,length(y_t1));
x_new1((length(p_t)+1)/2:(length(xn1)+(length(p_t)+1)/2)-1) = xn1;
x_new2 = zeros(1,length(y_t2));
x_new2((length(p_t)+1)/2:(length(xn2)+(length(p_t)+1)/2)-1) = xn2;
x_new3 = zeros(1,length(y_t3));
x_new3((length(p_t)+1)/2:(length(xn3)+(length(p_t)+1)/2)-1) = xn3;

figure();
subplot(3,1,1), plot(t,y_t1), grid on;
hold on
stem(t,x_new1);
title('Transmitted Signal 1'),xlabel('Time[s]'),ylabel('y_1(t)');
subplot(3,1,2), plot(t,y_t2), grid on;
hold on
stem(t,x_new2);
title('Transmitted Signal 2'),xlabel('Time[s]'),ylabel('y_2(t)');
subplot(3,1,3), plot(t,y_t3), grid on;
hold on
stem(t,x_new3);
title('Transmitted Signal 3'),xlabel('Time[s]'),ylabel('y_3(t)');


%y_t = [y_t1,y_t2,y_t3];

% figure(),plot(y_t),title('Transmitted Signal y(t)'),xlabel('Time[s]'),ylabel('y(t)');
% hold on;
% stem(x_new);


%% Text Message and Noise-free PAM Signal y(t)

message1 = 'Case study 2 rocks!';
bits1 = str2num(reshape(dec2bin(message1)',1,[])');
for i = 1:length(bits1)
    if bits1(i)==0
        bits1(i)=-1;
    end
end

message2 = 'Signals and systems';
bits2 = str2num(reshape(dec2bin(message2)',1,[])');
for i = 1:length(bits2)
    if bits2(i)==0
        bits2(i)=-1;
    end
end

message3 = 'Systems engineering';
bits3 = str2num(reshape(dec2bin(message3)',1,[])');
for i = 1:length(bits3)
    if bits3(i)==0
        bits3(i)=-1;
    end
end

xn1 = zeros(1,length(bits1)*(Ts/dt));
xn2 = zeros(1,length(bits2)*(Ts/dt));
xn3 = zeros(1,length(bits3)*(Ts/dt));

for i=1:length(bits1)
    xn1((i-1)*(Ts/dt)+1)=bits1(i);
end
for i=1:length(bits2)
    xn2((i-1)*(Ts/dt)+1)=bits2(i);
end
for i=1:length(bits3)
    xn3((i-1)*(Ts/dt)+1)=bits3(i);
end

y_t1 = conv(xn1, p_t);
y_t2 = conv(xn2, p_t);
y_t3 = conv(xn3, p_t);

%make sure all of the signals are same length
if length(y_t1)>length(y_t2)
    y_t2 = [y_t2,zeros(1,length(y_t1)-length(y_t2))];
    if length(y_t2)>length(y_t3)
        y_t3 = [y_t3,zeros(1,length(y_t2)-length(y_t3))];
    else
        y_t2 = [y_t2,zeros(1,length(y_t3)-length(y_t2))];
    end
else
    y_t1 = [y_t1,zeros(1,length(y_t2)-length(y_t1))];
    if length(y_t1)>length(y_t3)
        y_t3 = [y_t3,zeros(1,length(y_t1)-length(y_t3))];
    else
        y_t1 = [y_t1,zeros(1,length(y_t3)-length(y_t1))];
    end
end

N = length(y_t1)/(Ts/dt)-(Ts/dt);
t = (0:length(y_t1)-1)*dt; %define new time vector

%shift x so lines up with y(t) on graph
x_new1 = zeros(1,length(y_t1));
x_new1((length(p_t)+1)/2:(length(xn1)+(length(p_t)+1)/2)-1) = xn1;
x_new2 = zeros(1,length(y_t2));
x_new2((length(p_t)+1)/2:(length(xn2)+(length(p_t)+1)/2)-1) = xn2;
x_new3 = zeros(1,length(y_t3));
x_new3((length(p_t)+1)/2:(length(xn3)+(length(p_t)+1)/2)-1) = xn3;

figure();
subplot(3,1,1), plot(t,y_t1), grid on;
hold on
stem(t,x_message_new1);
title('Transmitted Signal 1'),xlabel('Time[s]'),ylabel('y_1(t)');
subplot(3,1,2), plot(t,y_t2), grid on;
hold on
stem(t,x_message_new2);
title('Transmitted Signal 2'),xlabel('Time[s]'),ylabel('y_2(t)');
subplot(3,1,3), plot(t,y_t3), grid on;
hold on
stem(t,x_message_new3);
title('Transmitted Signal 3'),xlabel('Time[s]'),ylabel('y_3(t)');

%% Signal Modulation (Up-Conversion)
wc20 = 2*pi*20; %20 Hz modulation
mod_signal20 = y_t1.*cos(wc20*t);

wc30 = 2*pi*30; %30 Hz modulation
mod_signal30 = y_t2.*cos(wc30*t);

wc40 = 2*pi*40; %40 Hz modulation
mod_signal40 = y_t3.*cos(wc40*t);

mod_signal = mod_signal20+mod_signal30+mod_signal40;

%plot the combined modulated signal
figure
subplot(2,1,1), plot(mod_signal), grid on;
xlabel('time (s)'), ylabel('y(t)cos(wc*t)'), title('Modulated Signal')

f = (0:length(mod_signal)/2)*fs/length(mod_signal);
XV = fft(mod_signal);
P2 = abs(XV/length(mod_signal));
P1 = P2(1:length(mod_signal)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,1,2),plot(f,P1);
xlabel('frequency[Hz]')
ylabel('|Y(j\omega)|')

%plot three signals separately
figure
subplot(2,3,1),plot(mod_signal20),grid on;
f = (0:length(mod_signal20)/2)*fs/length(mod_signal20);
XV = fft(mod_signal20);
P2 = abs(XV/length(mod_signal20));
P1 = P2(1:length(mod_signal20)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,3,4),plot(f,P1),xlabel('frequency[Hz]'),ylabel('|Y(j\omega)|'),grid on;

subplot(2,3,2),plot(mod_signal30),grid on;
f = (0:length(mod_signal30)/2)*fs/length(mod_signal30);
XV = fft(mod_signal30);
P2 = abs(XV/length(mod_signal30));
P1 = P2(1:length(mod_signal30)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,3,5),plot(f,P1),xlabel('frequency[Hz]'),ylabel('|Y(j\omega)|'),grid on;

subplot(2,3,3),plot(mod_signal40),grid on;
f = (0:length(mod_signal40)/2)*fs/length(mod_signal40);
XV = fft(mod_signal40);
P2 = abs(XV/length(mod_signal40));
P1 = P2(1:length(mod_signal40)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,3,6),plot(f,P1),xlabel('frequency[Hz]'),ylabel('|Y(j\omega)|'),grid on;


%% Noisy Recieved Signal r(t)
sigma = 0.5;
n_t = sigma*randn(1,length(mod_signal));
r_t = mod_signal + n_t;

figure();
plot(r_t);
title('Received Signal r(t)');
xlabel('Time[s]');
ylabel('Amplitude');

%% Signal Demodulation (Down-Conversion)

demod_signal20 = r_t.*cos(wc20*t);
demod_signal30 = r_t.*cos(wc30*t);
demod_signal40 = r_t.*cos(wc40*t);

%plot demod_signal20
figure
subplot(2,3,1), plot(demod_signal20), grid on;
xlabel('time (s)'), ylabel('x_r(t)'), title('Demod (without LPF)')

f = (0:length(demod_signal20)/2)*fs/length(demod_signal20);
XV = fft(demod_signal20);
P2 = abs(XV/length(demod_signal20));
P1 = P2(1:length(demod_signal20)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,3,4),plot(f,P1),xlabel('f, Hz'),ylabel('|X_r(j\omega)|'),grid on;

%plot demod_signal30
subplot(2,3,2), plot(demod_signal30), grid on;
xlabel('time (s)'), ylabel('x_r(t)'), title('Demod (without LPF)')

f = (0:length(demod_signal30)/2)*fs/length(demod_signal30);
XV = fft(demod_signal30);
P2 = abs(XV/length(demod_signal30));
P1 = P2(1:length(demod_signal30)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,3,5),plot(f,P1),xlabel('f, Hz'),ylabel('|X_r(j\omega)|'),grid on;

%plot demod_signal40
subplot(2,3,3), plot(demod_signal40), grid on;
xlabel('time (s)'), ylabel('x_r(t)'), title('Demod (without LPF)')

f = (0:length(demod_signal40)/2)*fs/length(demod_signal40);
XV = fft(demod_signal40);
P2 = abs(XV/length(demod_signal40));
P1 = P2(1:length(demod_signal40)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,3,6),plot(f,P1),xlabel('f, Hz'),ylabel('|X_r(j\omega)|'),grid on;


%% Lowpass Filter
R = 1000;
C = 2*10^(-6);
T = 0.002;
f = 44100;
delta_t = 1/f;

filter_low = filter(delta_t/T,[(R*C)/delta_t 1-((R*C)/delta_t)], demod_signal20);

figure();
plot(filter_low);

%% Matched Filter Receiver

%generate z(t)
xn_tilda = zeros(1,length(bits1));
p_neg = p_t(end:-1:1);
z_t = conv(demod_signal20,p_neg);

%generate estimate xn_tilda
for j=1:length(xn_tilda)
   if z_t((j-1)*(Ts/dt)+1+(length(p_t)+1)/2)<=0
       xn_tilda(j)=-1;
   else
       xn_tilda(j)=1;
   end
end

%space out xn_tilda
xn_spaced = zeros(1,N*(Ts/dt));
for i=1:length(xn_tilda)
    xn_spaced((i-1)*(Ts/dt)+1)=xn_tilda(i);
end
%shift xn_spaced
xn_spaced_new = zeros(1,length(z_t));
xn_spaced_new((2*length(p_t)+1)/2:(length(xn1)+(2*length(p_t)+1)/2)-1) = xn_spaced;

%calculate error rate
incorrect_count = 0;
for j=1:length(xn_tilda)
    if bits1(j)~=xn_tilda(j)
        incorrect_count = incorrect_count+1;
    end
end
error = incorrect_count/length(bits1);

t_new = (0:length(z_t)-1)*dt; %define new time vector

figure();
stem(bits1);
hold on;
stem(xn_tilda);
legend('bits1','xn tilda');

figure();
plot(t_new,z_t);
hold on
stem(t_new,xn_spaced_new);

%% Recover text message
for i = 1:length(xn_tilda)
    if xn_tilda(i)==-1
        xn_tilda(i)=0;
    end
end
messageOut = char(bin2dec(num2str(reshape(xn_tilda,7,[])')))';

%% Noise Levels and Error Rates

error_matched = zeros(1,11);
index = 1;
for i = 0:0.1:1
    sigma = i;
    n = sigma*randn(1,length(y));
    r = y(1:length(t))+n(1:length(t)); 
    xn_tilda_matched = ones(1,length(t));
    p_neg = p(end:-1:1);
    z = conv(r,p_neg);
    for j=1:length(z)
        if z(j)<=0
            xn_tilda_matched(j)=-1;
        end
    end
    %calculate error rate for matched filter
    incorrect_count_matched = 0;
    for j=1:length(x)
        if x(j)~=xn_tilda_matched(j)
            incorrect_count_matched = incorrect_count_matched+1;
        end
    end
    error_matched(index) = incorrect_count_matched/length(y);  
    index = index+1;
end

figure();
plot(sigma_vector,error_matched);
title('Sigma Value vs. Error Rate for Matched Filter');
xlabel('Sigma Value');
ylabel('Error Rate');


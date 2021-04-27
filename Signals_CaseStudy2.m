% Case Study 2
% Leandre Pestcoe and Julianne Wegmann
% ESE 351: Signals and Systems
% Created on: 4/20/21, Last Edited on: 4/26/21

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
xlabel('Time[s]'), ylabel('p(t)'), title('Truncated Signal p(t)')

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
xlabel('Time[s]'), ylabel('y(t)cos(wc*t)'), title('Modulated Signal with Three Messages')

f = (0:length(mod_signal)/2)*fs/length(mod_signal);
XV = fft(mod_signal);
P2 = abs(XV/length(mod_signal));
P1 = P2(1:length(mod_signal)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,1,2),plot(f,P1);
xlabel('Frequency[Hz]')
ylabel('|Y(j\omega)|')

%plot three signals separately
figure
sgtitle('Individual Modulated Signals');
subplot(2,3,1),plot(mod_signal20),title('Message 1'),grid on;
f = (0:length(mod_signal20)/2)*fs/length(mod_signal20);
XV = fft(mod_signal20);
P2 = abs(XV/length(mod_signal20));
P1 = P2(1:length(mod_signal20)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,3,4),plot(f,P1),title('Modulated at 20Hz'),xlabel('Frequency[Hz]'),ylabel('|Y(j\omega)|'),grid on;

subplot(2,3,2),plot(mod_signal30),title('Message 2'),grid on;
f = (0:length(mod_signal30)/2)*fs/length(mod_signal30);
XV = fft(mod_signal30);
P2 = abs(XV/length(mod_signal30));
P1 = P2(1:length(mod_signal30)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,3,5),plot(f,P1),title('Modulated at 30Hz'),xlabel('Frequency[Hz]'),ylabel('|Y(j\omega)|'),grid on;

subplot(2,3,3),plot(mod_signal40),title('Message 3'),grid on;
f = (0:length(mod_signal40)/2)*fs/length(mod_signal40);
XV = fft(mod_signal40);
P2 = abs(XV/length(mod_signal40));
P1 = P2(1:length(mod_signal40)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,3,6),plot(f,P1),title('Modulated at 40Hz'),xlabel('Frequency[Hz]'),ylabel('|Y(j\omega)|'),grid on;


%% Noisy Recieved Signal r(t)
sigma = 0.4;
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
sgtitle('Individual Demodulated Signals');
subplot(2,3,1), plot(demod_signal20), grid on;
xlabel('Time[s]'), ylabel('x_r(t)'), title('Message 1')

f = (0:length(demod_signal20)/2)*fs/length(demod_signal20);
XV = fft(demod_signal20);
P2 = abs(XV/length(demod_signal20));
P1 = P2(1:length(demod_signal20)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,3,4),plot(f,P1),title('Demodulated at 20 Hz'),xlabel('Frequency[Hz]'),ylabel('|X_r(j\omega)|'),grid on;

%plot demod_signal30
subplot(2,3,2), plot(demod_signal30), grid on;
xlabel('Time[s]'), ylabel('x_r(t)'), title('Message 2')

f = (0:length(demod_signal30)/2)*fs/length(demod_signal30);
XV = fft(demod_signal30);
P2 = abs(XV/length(demod_signal30));
P1 = P2(1:length(demod_signal30)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,3,5),plot(f,P1),title('Demodulated at 30 Hz'),xlabel('Frequency[Hz]'),ylabel('|X_r(j\omega)|'),grid on;

%plot demod_signal40
subplot(2,3,3), plot(demod_signal40), grid on;
xlabel('Time[s]'), ylabel('x_r(t)'), title('Message 3')

f = (0:length(demod_signal40)/2)*fs/length(demod_signal40);
XV = fft(demod_signal40);
P2 = abs(XV/length(demod_signal40));
P1 = P2(1:length(demod_signal40)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,3,6),plot(f,P1),title('Demodulated at 40 Hz'),xlabel('Frequency[Hz]'),ylabel('|X_r(j\omega)|'),grid on;


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
xn_tilda1 = zeros(1,length(bits1));
xn_tilda2 = zeros(1,length(bits2));
xn_tilda3 = zeros(1,length(bits3));
p_neg = p_t(end:-1:1);
z_t1 = conv(demod_signal20,p_neg);
z_t2 = conv(demod_signal30,p_neg);
z_t3 = conv(demod_signal40,p_neg);

%generate estimate xn_tilda
for j=1:length(xn_tilda1)
   if z_t1((j-1)*(Ts/dt)+1+(length(p_t)+1))<=0
       xn_tilda1(j)=-1;
   else
       xn_tilda1(j)=1;
   end
end
for j=1:length(xn_tilda2)
   if z_t2((j-1)*(Ts/dt)+1+(length(p_t)+1))<=0
       xn_tilda2(j)=-1;
   else
       xn_tilda2(j)=1;
   end
end
for j=1:length(xn_tilda3)
   if z_t3((j-1)*(Ts/dt)+1+(length(p_t)+1))<=0
       xn_tilda3(j)=-1;
   else
       xn_tilda3(j)=1;
   end
end

%space out xn_tilda
xn_spaced1 = zeros(1,N*(Ts/dt));
for i=1:length(xn_tilda1)
    xn_spaced1((i-1)*(Ts/dt)+1)=xn_tilda1(i);
end
xn_spaced2 = zeros(1,N*(Ts/dt));
for i=1:length(xn_tilda2)
    xn_spaced2((i-1)*(Ts/dt)+1)=xn_tilda2(i);
end
xn_spaced3 = zeros(1,N*(Ts/dt));
for i=1:length(xn_tilda3)
    xn_spaced3((i-1)*(Ts/dt)+1)=xn_tilda3(i);
end

%shift xn_spaced
xn_spaced_new1 = zeros(1,length(z_t1));
xn_spaced_new1((2*length(p_t)+1)/2:(length(xn_spaced1)+(2*length(p_t)+1)/2)-1) = xn_spaced1;
xn_spaced_new2 = zeros(1,length(z_t2));
xn_spaced_new2((2*length(p_t)+1)/2:(length(xn_spaced2)+(2*length(p_t)+1)/2)-1) = xn_spaced2;
xn_spaced_new3 = zeros(1,length(z_t3));
xn_spaced_new3((2*length(p_t)+1)/2:(length(xn_spaced3)+(2*length(p_t)+1)/2)-1) = xn_spaced3;

%calculate error rate
incorrect_count1 = 0;
for j=1:length(xn_tilda1)
    if bits1(j)~=xn_tilda1(j)
        incorrect_count1 = incorrect_count1+1;
    end
end
error1 = incorrect_count1/length(bits1);
incorrect_count2 = 0;
for j=1:length(xn_tilda2)
    if bits2(j)~=xn_tilda2(j)
        incorrect_count2 = incorrect_count2+1;
    end
end
error2 = incorrect_count2/length(bits2);
incorrect_count3 = 0;
for j=1:length(xn_tilda3)
    if bits3(j)~=xn_tilda3(j)
        incorrect_count3 = incorrect_count3+1;
    end
end
error3 = incorrect_count3/length(bits3);

t_new = (0:length(z_t1)-1)*dt; %define new time vector

figure();
sgtitle('Original Bit Messages vs. Estimates');
subplot(3,1,1),stem(bits1);
hold on;
stem(xn_tilda1);
title('Message 1'),legend('bits','xn tilda');
subplot(3,1,2),stem(bits2);
hold on;
stem(xn_tilda2);
title('Message 2'),legend('bits','xn tilda');
subplot(3,1,3),stem(bits3);
hold on;
stem(xn_tilda3);
title('Message 3'),legend('bits','xn tilda');

figure();
sgtitle('Recovered Individual Signals');
subplot(3,1,1),plot(t_new,z_t1),title('Message 1'),ylabel('z_1(t)'),xlabel('Time[s]');
hold on
stem(t_new,xn_spaced_new1);
subplot(3,1,2),plot(t_new,z_t2),title('Message 2'),ylabel('z_2(t)'),xlabel('Time[s]');
hold on
stem(t_new,xn_spaced_new2);
subplot(3,1,3),plot(t_new,z_t3),title('Message 3'),ylabel('z_3(t)'),xlabel('Time[s]');
hold on
stem(t_new,xn_spaced_new3);

%% Recover Text Message
for i = 1:length(xn_tilda1)
    if xn_tilda1(i)==-1
        xn_tilda1(i)=0;
    end
end
messageOut1 = char(bin2dec(num2str(reshape(xn_tilda1,7,[])')))';

for i = 1:length(xn_tilda2)
    if xn_tilda2(i)==-1
        xn_tilda2(i)=0;
    end
end
messageOut2 = char(bin2dec(num2str(reshape(xn_tilda2,7,[])')))';

for i = 1:length(xn_tilda3)
    if xn_tilda3(i)==-1
        xn_tilda3(i)=0;
    end
end
messageOut3 = char(bin2dec(num2str(reshape(xn_tilda3,7,[])')))';

disp(messageOut1);
disp(messageOut2);
disp(messageOut3);

%% Noise Levels and Error Rates

%random bit message
error1 = [0,0,0,0,0,0,0,0,.1,.1,.1];
error2 = [0,0,0,0,0,0,0,0,0,0,0];
error3 = [0,0,0,0,0,0,0,0,.1,.1,.1];

sigma = (0:0.1:1);

figure,plot(sigma,error1);
hold on
plot(sigma,error2);
hold on
plot(sigma,error3);
legend('1','2','3');


%%


error = zeros(1,N+1);
index = 1;
for i = 0:0.1:1
    sigma = i;
    n = sigma*randn(1,length(mod_signal));
    r = y(1:length(t))+n(1:length(t)); 
    xn_tilda_matched = ones(1,length(t));
    p_neg = p(end:-1:1);
    z = conv(r,p_neg);
    for j=1:length(xn_tilda1)
        if z_t1((j-1)*(Ts/dt)+1+(length(p_t)+1))<=0
            xn_tilda1(j)=-1;
        else
            xn_tilda1(j)=1;
        end
    end
    %calculate error rate for matched filter
    incorrect_count = 0;
    for j=1:length(bits1)
        if bits1(j)~=xn_tilda1(j)
            incorrect_count = incorrect_count+1;
        end
    end
    error(index) = incorrect_count/length(y);  
    index = index+1;
end

figure();
plot(sigma_vector,error_matched);
title('Sigma Value vs. Error Rate for Matched Filter');
xlabel('Sigma Value');
ylabel('Error Rate');


clc
%clear
close all;

%% SIMULATION PARAMETERS
UPSAMPLING = 10;       % Numbers of times that the clk period is bigger than the carrier_period
BIT_PERIOD = 50;       % Numbers of times that the bit period is bigger than the clk period
DATA_BITS = 50;        % Numbers of bits of the message
MOD = 1;               % MOD = 0 => sin carrier , MOD = 1 => square carrier 
%% BASE TIME SIGNALS

carrier_freq = 10e3;
carrier_period = 1/carrier_freq;

clk_freq = carrier_freq * UPSAMPLING;
clk_period = 1/clk_freq;

N = DATA_BITS * BIT_PERIOD;    

t = 0:clk_period:(N-1)*clk_period;



%% SIGNALS

% Digital signal from GNU RADIO
data_signal = randi([0 1],1,DATA_BITS);

% RF carrier for analog OOK
carrier_analog = sin(carrier_freq*2*pi*t);

% RF carrier from PLL2
carrier_digital = (square(carrier_freq*2*pi*t)+1)/2;


    if   (MOD==0)
         carrier = carrier_analog;
    else (MOD==1);
         carrier = carrier_digital;
    end
    
    
message=[];

% Loop to replicate bits from original signal
for n=1:length(data_signal);
    if data_signal(n)==0;
        ext=zeros(1,BIT_PERIOD);  
    else data_signal(n)==1;
        ext=ones(1,BIT_PERIOD);     
    end

     
    message=[message ext];
end


%% MODULATION
ook=message.*carrier;
ook_filtered = filter(BPF, ook); 


%% MODULATION TIME PLOTS 

figure('Name','Time Domain Modulation Stage','NumberTitle','off')

tm(1)=subplot(4,1,1);
plot(t,message,'LineWidth',1.25);grid on;
title('Binary Signal');
xlabel('t_{[seconds]}','Color','r');
axis([0 (max(t)) -2.5 2.5]);

tm(2)=subplot(4,1,2);
plot(t,carrier,'LineWidth',1.25);grid on;
title('RF Carrier');
xlabel('t_{[seconds]}','Color','r');
axis([0 max(t) -2 2]);


tm(3)= subplot(4,1,3);
plot(t,ook,'LineWidth',1.25);grid on;
title('OOK modulation');
xlabel('t_{[seconds]}','Color','r');
axis([0 max(t) -1 2]);

tm(4)= subplot(4,1,4);
plot(t,ook_filtered,'LineWidth',1.25);grid on;
title('OOK modulation Filtered Antenna');
xlabel('t_{[seconds]}','Color','r');
axis([0 (max(t)) -2.5 2.5]);

linkaxes(tm, 'x');



%% FRECUENCY MODULATION
spectral_resolution = clk_freq / N;
f = (0:N-1)*spectral_resolution;

data_spec = abs(fft(message)/100);
mod_spec  = abs(fft(carrier)/100);
ook_spec  = abs(fft(ook)/100);

%% FRECUENCY MODULATION PLOTS

figure('Name','Frecuency Domain Modulation Stage','NumberTitle','off')

fm(1) = subplot(3,1,1);
plot(f(1:length(data_spec)/2),data_spec(1:(length(data_spec)/2)),'LineWidth',1.25);grid on;
title('Binary Signal Spectrum');
xlabel('f_{[hertz]}','Color','r');
%axis([0 length(f)/2 -1 10]);


fm(2) = subplot(3,1,2);
plot(f(1:length(mod_spec)/2),mod_spec(1:(length(mod_spec)/2)),'LineWidth',1.25);grid on;
title('RF Carrier Spectrum');
xlabel('f_{[hertz]}','Color','r');
%axis([0 length(f)/2 -1 10]);

fm(3) = subplot(3,1,3);
plot(f(1:length(ook_spec)/2),ook_spec(1:(length(ook_spec)/2)),'LineWidth',1.25);grid on;
title('OOK modulation Spectrum');
xlabel('f_{[hertz]}','Color','r');
%axis([0 length(f)/2 -1 10]);

linkaxes(fm, 'x'); 


%% DEMODULATION TIME PLOTS
 
demod = carrier_analog.*ook_filtered;
message_demod = filter(LPF, demod);


figure('Name','Time Domain Demodulation Stage','NumberTitle','off')
 
td(1)= subplot(3,1,1);
plot(t,ook_filtered,'LineWidth',1.25);grid on;
title('OOK modulation Filtered Antenna');
xlabel('t_{[seconds]}','Color','r');
axis([0 (max(t)) -2.5 2.5]);

td(2)= subplot(3,1,2);
plot(t,carrier_analog,'LineWidth',1.25);grid on;
title('Analog Carrier');
xlabel('t_{[seconds]}','Color','r');
axis([0 (max(t)) -2.5 2.5]);

td(3)=subplot(3,1,3);
plot(t,message_demod,'LineWidth',1.25);grid on;
title('Message demodulated');
xlabel('t_{[seconds]}','Color','r');
axis([0 (max(t)) -1 2]);

linkaxes(td, 'x'); 


%% FRECUENCY DEMODULATION

 ook_filtered_spec = abs(fft(ook_filtered));
 demod_spec = abs(fft(demod));
 message_demod_spec = abs(fft(message_demod));
 

 %% FRECUENCY DEMODULATION PLOTS

figure('Name','Frecuency Domain Modulation Stage','NumberTitle','off')

fd(1)= subplot(3,1,1);
plot(f(1:length(ook_filtered_spec)/2),ook_filtered_spec(1:(length(ook_filtered_spec)/2)),'LineWidth',1.25);grid on;
title('OOK Received Spectrum');
xlabel('f_{[hertz]}','Color','r');
%axis([0 (max(t)) -2.5 2.5]);

fd(2)= subplot(3,1,2);
plot(f(1:length(demod_spec)/2),demod_spec(1:(length(demod_spec)/2)),'LineWidth',1.25);grid on;
title('Message demodulated Spectrum');
xlabel('f_{[hertz]}','Color','r');
%axis([0 (max(t)) -2.5 2.5]);


fd(3)=subplot(3,1,3);
plot(f(1:length(message_demod_spec)/2),message_demod_spec(1:(length(message_demod_spec)/2)),'LineWidth',1.25);grid on;
title('Message demodulated & filtered Spectrum');
xlabel('f_{[hertz]}','Color','r');
%axis([0 (max(t)) -2.5 2.5]);

linkaxes(fd, 'x'); 




clc; ; close all;

%% SIMULATION PARAMETERS
DATA_BITS = 6;
PWM_STEPS = 2^DATA_BITS;
UPSAMPLING = 100;
DATA_CYCLES = 2;

carrier_freq = 10e6;
carrier_period = 1/carrier_freq;

clk_freq = carrier_freq * UPSAMPLING;
clk_period = 1/clk_freq;

% TEST 01 -> pwm_freq is (clk_freq * 2^DATA_BITS)
pwm_period = carrier_period * PWM_STEPS;
pwm_freq = 1/pwm_period;

% % TEST 02 -> pwm_freq is data_samp_freq
% pwm_freq = 50e3;
% pwm_period = 1/pwm_freq;

data_freq = 1e3;
data_period = 1/data_freq;
data_amp = (2^DATA_BITS)/2;

N = DATA_CYCLES * (data_period/clk_period);


%% SIGNALS
% time vector
t = (0:N-1)*clk_period;

% RF carrier from PLL2
carrier = (square(carrier_freq*2*pi*t) + 1)/2;

% data from GNURadio
data_an = data_amp * (sin(data_freq*2*pi*t) + 1);
data_quant = max(data_an*2) / (2^DATA_BITS -1);
data_dig = uint8(round(data_an / data_quant));

% PWM from data_dig
% TEST 01 -> pwm_freq is clk_freq * 2^DATA_BITS

pwm = uint8(zeros(1,N));

for i = 1:PWM_STEPS:N
    for j = 0:PWM_STEPS-1
        if j < data_dig(i)
            pwm(i+j) = 1;
        else
            pwm(i+j) = 0;
        end
    end
end

% % PWM from data_dig
% % TEST 02 -> pwm_freq is data_samp_freq
% pwm = uint8(zeros(1,N));
% 
% for i = 1 : pwm_period*N/(DATA_CYCLES*data_period) : N
%     aux = data_dig(i);
%     for j = 0 : pwm_period*N/(DATA_CYCLES*data_period) - 1
%         data_dig(i+j) = aux;
%         if j < aux
%             pwm(i+j) = 1;
%         else
%             pwm(i+j) = 0;
%         end
%     end
% end


%% PLOT
figure(1);
subplot(221);   plot(t, carrier);       title("carrier");
subplot(222);   plot(t, data_an);       title("data - analog");
subplot(223);   plot(t, data_dig);      title("data - digital");
subplot(224);   plot(t, pwm);           title("PWM");


%% FREQUENCY
res_freq = clk_freq/N;
f = (0:N-1)*res_freq;

carrier_spec = abs(fft(carrier))/N;
data_an_spec = abs(fft(data_an))/N;
data_dig_spec = abs(fft(double(data_dig)))/N;
pwm_spec = abs(fft(double(pwm)))/N;


%% PLOT
figure(2);
subplot(221);   plot(f(1:N/2), mag2db(carrier_spec(1:N/2)));    title("carrier");
subplot(222);   plot(f(1:N/2), mag2db(data_an_spec(1:N/2)));    title("data - analog");
subplot(223);   plot(f(1:N/2), (data_dig_spec(1:N/2)));         title("data - digital");
subplot(224);   plot(f(1:N/2), (pwm_spec(1:N/2)));              title("PWM");



%% SIGNAL SIGNAL OUTPUT
am_dig = double(pwm) .* carrier;


%% FILTER
am_filtered = filter(BPF_NUM, 1, am_dig); 


%% FREQUENCY
am_filtered_spec = abs(fft(am_filtered))/N;


%% PLOT
% figure(3);
% subplot(111);   plot(f(1:N/2), mag2db(am_filtered_spec(1:N/2)));

figure(3);
subplot(211);   plot(f, am_filtered);
    title("am filtered");
    
subplot(212);   plot(f(.5*carrier_freq/res_freq : 1.5*carrier_freq/res_freq), ...
        mag2db(am_filtered_spec(.5*carrier_freq/res_freq : 1.5*carrier_freq/res_freq)));
    title("am filtered - spec");
% subplot(212);   plot(f, mag2db(am_filtered_spec));
%     title("am filtered - spec");


%% DEMODULATION
am_demod = amdemod(am_filtered, carrier_freq, clk_freq);
am_demod_spec = mag2db(abs(fft(am_demod))/N);


%% PLOT
figure(4);
subplot(311);   plot(t, am_demod);
    title("am demod");
subplot(312);   plot(f(1:N/2), am_demod_spec(1:N/2));
    title("am demod - full spec");
subplot(313);   plot(f(1:4*data_freq/res_freq), am_demod_spec(1:4*data_freq/res_freq));
    title("am demod - base band");




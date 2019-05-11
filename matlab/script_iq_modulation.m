clc; ; close all;

%% SIMULATION PARAMETERS
SYMBOL_RATE = 5;
UPSAMPLING = 100;
DATA_CYCLES = 20;

carrier_freq = 10;
carrier_period = 1/carrier_freq;

% for matlab use only, sampling frequency of all simulation
samp_freq = carrier_freq * UPSAMPLING;
samp_period = 1/samp_freq;

symbol_period = carrier_period * SYMBOL_RATE;
symbol_freq = 1/symbol_period;

N = round(DATA_CYCLES * (symbol_period/samp_period));


%% SIGNALS
% time vector
t = (0:N-1)*samp_period;

% compose a signal: cos , -cos , sin , -sin
symb1 = cos(carrier_freq*2*pi*t(1:N/4));
symb2 = -cos(carrier_freq*2*pi*t(1:N/4));
symb3 = sin(carrier_freq*2*pi*t(1:N/4));
symb4 = -sin(carrier_freq*2*pi*t(1:N/4));
qpsk = [symb1 symb2 symb3 symb4];


%% PLOT
figure(1);
subplot(121);   plot(t, qpsk);           title("qpsk");


%% FREQUENCY
res_freq = samp_freq/N;
f = (0:N-1)*res_freq;

qpsk_spec = abs(fft(double(qpsk)))/N;


%% PLOT
subplot(122);   plot(f(1:N/2), mag2db(qpsk_spec(1:N/2)));   title("qpsk spectrum");
  

%% DEMODULATION
osc = cos(carrier_freq*2*pi*t);
i_demod = qpsk .* osc;
osc = sin(carrier_freq*2*pi*t);
q_demod = qpsk .* osc;

i_filtered = filter(Num, 1, i_demod); 
q_filtered = filter(Num, 1, q_demod); 
%i_filtered = mean(i_demod);
%q_filtered = mean(q_demod);

% qpsk_demod = ceil(i_filtered*2)/2 + 1i * ceil(q_filtered*2)/2;
% qpsk_demod = round(i_filtered, 1) + 1i * round(q_filtered, 1);
qpsk_demod = (i_filtered) + 1i * (q_filtered);
qpsk_message = pskdemod(qpsk_demod, 4);


%% PLOT
figure(3);
subplot(121);   plot(t, i_filtered);        title("i in time");
subplot(122);   plot(t, q_filtered);        title("q in time");

scatterplot(qpsk_demod);	title("qpsk demod");
scatterplot(qpsk_message);	title("qpsk message");

figure(4);
subplot(121);   plot(qpsk_message);         title("qpsk in time");
subplot(122);   histogram(qpsk_message);    title("histogram");





clc; ; close all;

%% SIMULATION PARAMETERS
DATA_BITS = 1;
SYMBOLS = 2^DATA_BITS;
SYMBOL_RATE = 10;
UPSAMPLING = 100;
DATA_CYCLES = 10;

clk_freq = 10e6;
clk_period = 1/clk_freq;

% for matlab use only, sampling frequency of all simulation
samp_freq = clk_freq * UPSAMPLING;
samp_period = 1/samp_freq;

symbol_period = clk_period * SYMBOL_RATE;
symbol_freq = 1/symbol_period;

N = round(DATA_CYCLES * (symbol_period/samp_period));


%% SIGNALS
% time vector
t = (0:N-1)*samp_period;

% CLK system from PLL1
clk = (square(clk_freq*2*pi*t) + 1)/2;

% data from GNURadio
data = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2];
data_index = 1;

% QPSK from data
qpsk = double(zeros(1,N));

%%%% GENERATES QPSK -> 4 PSK
for i = 1: SYMBOL_RATE*UPSAMPLING :N
    switch data(data_index)
        % phase 0 grades
        case 1
            for j = 0:SYMBOL_RATE*UPSAMPLING-1
                switch mod(j,4)
                    case 0
                        qpsk(i+j) = 1;
                    case 1
                        qpsk(i+j) = 1;
                    case 2
                        qpsk(i+j) = 0;
                    case 3
                        qpsk(i+j) = 0;
                end
            end
        % phase 180 grades
        case 2
            for j = 0:SYMBOL_RATE*UPSAMPLING-1
                switch mod(j,4)
                    case 0
                        qpsk(i+j) = 0;
                    case 1
                        qpsk(i+j) = 0;
                    case 2
                        qpsk(i+j) = 1;
                    case 3
                        qpsk(i+j) = 1;
                end
            end
        % phase 90 grades
        case 3
            for j = 0:SYMBOL_RATE*UPSAMPLING-1
                switch mod(j,4)
                    case 0
                        qpsk(i+j) = 0;
                    case 1
                        qpsk(i+j) = 1;
                    case 2
                        qpsk(i+j) = 1;
                    case 3
                        qpsk(i+j) = 0;
                end
            end
        % phase -90 grades
        case 4
            for j = 0:SYMBOL_RATE*UPSAMPLING-1
                switch mod(j,4)
                    case 0
                        qpsk(i+j) = 1;
                    case 1
                        qpsk(i+j) = 0;
                    case 2
                        qpsk(i+j) = 0;
                    case 3
                        qpsk(i+j) = 1;
                end
            end
    end
    data_index = data_index + 1;
end


%% PLOT
figure(1);
subplot(311);   plot(t, clk);           title("clk");
subplot(312);   plot(data);             title("data");
subplot(313);   plot(t, qpsk);          title("qpsk");


%% FREQUENCY
res_freq = samp_freq/N;
f = (0:N-1)*res_freq;

clk_spec = abs(fft(clk))/N;
qpsk_spec = abs(fft(double(qpsk)))/N;


%% PLOT
figure(2);
subplot(211);   plot(f(1:N/2), mag2db(clk_spec(1:N/2)));    title("clk");
subplot(212);   plot(f(1:N/2), mag2db(qpsk_spec(1:N/2)));   title("qpsk");


%% FILTER
qpsk_filtered = filter(BPF_NUM, 1, qpsk); 


%% FREQUENCY
qpsk_filtered_spec = abs(fft(qpsk_filtered))/N;


%% PLOT
figure(3);
subplot(211);   plot(t, qpsk_filtered);
    title("qpsk filtered");
subplot(212);   plot(f(1:N/2), mag2db(qpsk_filtered_spec(1:N/2)));
    title("qpsk filtered - spect");
    

%% DEMODULATION
qpsk_demod = pskdemod(qpsk_filtered, 4, 0);
qpsk_demod_spec = mag2db(abs(fft(qpsk_demod))/N);


%% PLOT
scatterplot(qpsk_demod);	title("qpsk demod");




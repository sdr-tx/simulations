clc; ; close all;

%% SIMULATION PARAMETERS
DATA_BITS = 1;
SYMBOLS = 2^DATA_BITS;
SYMBOL_RATE = 10;
UPSAMPLING = 100;
DATA_CYCLES = 4;

clk_freq = 100e3;
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
data = [1, 2, 3, 4];
data_index = 1;

% QPSK from data
qpsk_digital = double(zeros(1,N));

%%%% GENERATES QPSK -> 4 PSK
for i = 1: SYMBOL_RATE*UPSAMPLING :N
    switch data(data_index)
        % phase 0 grades
        case 1
            for j = 0:SYMBOL_RATE*UPSAMPLING-1
                switch mod(j,16)
                    case 0
                        qpsk_digital(i+j) = 1;
                    case 1
                        qpsk_digital(i+j) = 1;
                    case 2
                        qpsk_digital(i+j) = 1;
                    case 3
                        qpsk_digital(i+j) = 1;
                    case 4
                        qpsk_digital(i+j) = 1;
                    case 5
                        qpsk_digital(i+j) = 1;
                    case 6
                        qpsk_digital(i+j) = 1;
                    case 7
                        qpsk_digital(i+j) = 1;
                    case 8
                        qpsk_digital(i+j) = 0;
                    case 9
                        qpsk_digital(i+j) = 0;
                    case 10
                        qpsk_digital(i+j) = 0;
                    case 11
                        qpsk_digital(i+j) = 0;
                    case 12
                        qpsk_digital(i+j) = 0;
                    case 13
                        qpsk_digital(i+j) = 0;
                    case 14
                        qpsk_digital(i+j) = 0;
                    case 15
                        qpsk_digital(i+j) = 0;
                end
            end
        % phase 180 grades
        case 2
            for j = 0:SYMBOL_RATE*UPSAMPLING-1
                switch mod(j,16)
                    case 0
                        qpsk_digital(i+j) = 0;
                    case 1
                        qpsk_digital(i+j) = 0;
                    case 2
                        qpsk_digital(i+j) = 0;
                    case 3
                        qpsk_digital(i+j) = 0;
                    case 4
                        qpsk_digital(i+j) = 0;
                    case 5
                        qpsk_digital(i+j) = 0;
                    case 6
                        qpsk_digital(i+j) = 0;
                    case 7
                        qpsk_digital(i+j) = 0;
                    case 8
                        qpsk_digital(i+j) = 1;
                    case 9
                        qpsk_digital(i+j) = 1;
                    case 10
                        qpsk_digital(i+j) = 1;
                    case 11
                        qpsk_digital(i+j) = 1;
                    case 12
                        qpsk_digital(i+j) = 1;
                    case 13
                        qpsk_digital(i+j) = 1;
                    case 14
                        qpsk_digital(i+j) = 1;
                    case 15
                        qpsk_digital(i+j) = 1;
                end
            end
        % phase 90 grades
        case 3
            for j = 0:SYMBOL_RATE*UPSAMPLING-1
                switch mod(j,16)
                    case 0
                        qpsk_digital(i+j) = 0;
                    case 1
                        qpsk_digital(i+j) = 0;
                    case 2
                        qpsk_digital(i+j) = 0;
                    case 3
                        qpsk_digital(i+j) = 0;
                    case 4
                        qpsk_digital(i+j) = 1;
                    case 5
                        qpsk_digital(i+j) = 1;
                    case 6
                        qpsk_digital(i+j) = 1;
                    case 7
                        qpsk_digital(i+j) = 1;
                    case 8
                        qpsk_digital(i+j) = 1;
                    case 9
                        qpsk_digital(i+j) = 1;
                    case 10
                        qpsk_digital(i+j) = 1;
                    case 11
                        qpsk_digital(i+j) = 1;
                    case 12
                        qpsk_digital(i+j) = 0;
                    case 13
                        qpsk_digital(i+j) = 0;
                    case 14
                        qpsk_digital(i+j) = 0;
                    case 15
                        qpsk_digital(i+j) = 0;
                end
            end
        % phase -90 grades
        case 4
            for j = 0:SYMBOL_RATE*UPSAMPLING-1
                switch mod(j,16)
                    case 0
                        qpsk_digital(i+j) = 1;
                    case 1
                        qpsk_digital(i+j) = 1;
                    case 2
                        qpsk_digital(i+j) = 1;
                    case 3
                        qpsk_digital(i+j) = 1;
                    case 4
                        qpsk_digital(i+j) = 0;
                    case 5
                        qpsk_digital(i+j) = 0;
                    case 6
                        qpsk_digital(i+j) = 0;
                    case 7
                        qpsk_digital(i+j) = 0;
                    case 8
                        qpsk_digital(i+j) = 0;
                    case 9
                        qpsk_digital(i+j) = 0;
                    case 10
                        qpsk_digital(i+j) = 0;
                    case 11
                        qpsk_digital(i+j) = 0;
                    case 12
                        qpsk_digital(i+j) = 1;
                    case 13
                        qpsk_digital(i+j) = 1;
                    case 14
                        qpsk_digital(i+j) = 1;
                    case 15
                        qpsk_digital(i+j) = 1;
                end
            end
    end
    data_index = data_index + 1;
end


%% PLOT
figure(1);
subplot(311);   plot(t, clk);           title("clk");
subplot(312);   plot(data);             title("data");
subplot(313);   plot(t, qpsk_digital);  title("qpsk");


%% FREQUENCY
res_freq = samp_freq/N;
f = (0:N-1)*res_freq;

clk_spec = abs(fft(clk))/N;
qpsk_digital_spec = abs(fft(double(qpsk_digital)))/N;


%% PLOT
figure(2);
subplot(211);   plot(f(1:N/2), mag2db(clk_spec(1:N/2)));            title("clk");
subplot(212);   plot(f(1:N/2), (qpsk_digital_spec(1:N/2)));         title("qpsk digital specter");


%% FILTER
qpsk = filter(BPF_NUM, 1, qpsk_digital); 


%% FREQUENCY
qpsk_filtered_spec = abs(fft(qpsk))/N;


%% PLOT
figure(3);
subplot(211);   plot(t, qpsk);
    title("qpsk filtered");
subplot(212);   plot(f(1:N/2), mag2db(qpsk_filtered_spec(1:N/2)));
    title("qpsk filtered - spect");
    

%% DEMODULATION
carrier_freq = 6.25e5;
osc = cos(carrier_freq*2*pi*t);
i_demod = qpsk .* osc;
osc = sin(carrier_freq*2*pi*t);
q_demod = qpsk .* osc;

i_filtered = filter(Num, 1, i_demod); 
q_filtered = filter(Num, 1, q_demod); 

qpsk_demod = (i_filtered) + 1i * (q_filtered);
qpsk_message = pskdemod(qpsk_demod, 4);




%% PLOT
figure(8);
% subplot(121);   plot(t, i_filtered);        title("i in time");
% subplot(122);   plot(t, q_filtered);        title("q in time");
plot(t, i_filtered, t, q_filtered);         title("i-q in time");

figure(7);
subplot(121);   histogram(i_filtered);      title("i histogram");
subplot(122);   histogram(q_filtered);      title("q histogram");

% figure(5);
subplot(121);   scatterplot(qpsk_demod);	title("qpsk demod");
subplot(122);   scatterplot(qpsk_message);	title("qpsk message");

figure(4);
subplot(121);   plot(qpsk_message);         title("message in time");
subplot(122);   histogram(qpsk_message);    title("message histogram");



clc; ; close all;

%% SIMULATION PARAMETERS
UPSAMPLING = 100;
PWM_RESOLUTION = 20;
SYMBOL_RATE = 3;

SYMBOLS = 4;
DATA_CYCLES = 10;

% DATA_BITS = 1;
% SYMBOLS = 2^DATA_BITS;

%% SIGNALS FREQUENCY
clk_freq = 3000e3;
clk_period = 1/clk_freq;

% for matlab use only, sampling frequency of all simulation
samp_freq = clk_freq * UPSAMPLING;
samp_period = 1/samp_freq;

pwm_period = clk_period * PWM_RESOLUTION;
pwm_freq = 1/pwm_period;

symbol_period = pwm_period * SYMBOL_RATE;
symbol_freq = 1/symbol_period;

N = round(DATA_CYCLES * SYMBOLS * (symbol_period/samp_period));


%% SIGNALS
% time vector
t = (0:N-1)*samp_period;

% data from GNURadio
data = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4  ...
        1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4  ...
        1, 1, 1, 1, 2, 2, 2, 2];
data_index = 1;

% QPSK from data
qpsk_digital = zeros(1, PWM_RESOLUTION);

% (0.25+0.25j) {d[0]: 0.332210726289157, d[1]: 0.555963784362220, d[2]: 0.611825489348623}
sym1 = [7, 11, 12];
% (-0.25+0.25j) {d[0]: 0.444036215637780, d[1]: 0.667789273710843, d[2]: 0.388174510651377}
sym2 = [9, 13, 8];
% (-0.25-0.25j) {d[0]: 0.686956813972096, d[1]: 0.461135716596282, d[2]: 0.351907469431621}
sym3 = [14, 9, 7];
% (0.25-0.25j) {d[0]: 0.538864283403718, d[1]: 0.313043186027904, d[2]: 0.648092530568379}
sym4 = [11, 6, 13];

symbol1 = function_QamCreateSymbol( PWM_RESOLUTION*UPSAMPLING, sym1(1)*UPSAMPLING, sym1(2)*UPSAMPLING, sym1(3)*UPSAMPLING );
symbol2 = function_QamCreateSymbol( PWM_RESOLUTION*UPSAMPLING, sym2(1)*UPSAMPLING, sym2(2)*UPSAMPLING, sym2(3)*UPSAMPLING );
symbol3 = function_QamCreateSymbol( PWM_RESOLUTION*UPSAMPLING, sym3(1)*UPSAMPLING, sym3(2)*UPSAMPLING, sym3(3)*UPSAMPLING );
symbol4 = function_QamCreateSymbol( PWM_RESOLUTION*UPSAMPLING, sym4(1)*UPSAMPLING, sym4(2)*UPSAMPLING, sym4(3)*UPSAMPLING );

%%%% GENERATES QPSK -> 4 PSK
for i = 1: SYMBOL_RATE*PWM_RESOLUTION*UPSAMPLING :N
    if i == 1
        switch data(data_index)
            case 1
                qpsk_digital = symbol1;
            case 2
                qpsk_digital = symbol2;
            case 3
                qpsk_digital = symbol3;
            case 4
                qpsk_digital = symbol4;
        end
    else
        switch data(data_index)
            case 1
                qpsk_digital = [qpsk_digital symbol1];
            case 2
                qpsk_digital = [qpsk_digital symbol2];
            case 3
                qpsk_digital = [qpsk_digital symbol3];
            case 4
                qpsk_digital = [qpsk_digital symbol4];
        end
    end
    
    data_index = data_index + 1;
end


%% PLOT
figure(1);
    subplot(311);   plot(data);             
        title("data");
    subplot(312);   plot(t, qpsk_digital);  
        title("qpsk digital");


%% FREQUENCY
res_freq = samp_freq/N;
f = (0:N-1)*res_freq;

qpsk_digital_spec = abs(fft(double(qpsk_digital)))/N;


%% PLOT
    subplot(313);   plot(f(1:N/20), (qpsk_digital_spec(1:N/20)));
        title("qpsk digital specter");


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
carrier_freq = pwm_freq;
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
scatterplot(qpsk_demod);	title("qpsk demod");
scatterplot(qpsk_message);	title("qpsk message");

figure(4);
subplot(121);   plot(qpsk_message);         title("message in time");
subplot(122);   histogram(qpsk_message);    title("message histogram");



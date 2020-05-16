clc; ; close all; 

%% SIMULATION PARAMETERS
DATA_BITS = 1;
SYMBOLS = 2^DATA_BITS;
SYMBOL_RATE = 10;
UPSAMPLING = 100;
DATA_CYCLES = 10;

carrier_freq = 10e6;
carrier_period = 1/carrier_freq;

clk_freq = carrier_freq * UPSAMPLING;
clk_period = 1/clk_freq;

symbol_period = carrier_period * SYMBOL_RATE;
symbol_freq = 1/symbol_period;

N = round(DATA_CYCLES * (symbol_period/clk_period));


%% SIGNALS
% time vector
t = (0:N-1)*clk_period;

% RF carrier from PLL2
carrier = (square(carrier_freq*2*pi*t) + 1)/2;

% data from GNURadio
% data = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2];
data = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2];
data_index = 1;

% BPSK from data
bpsk = double(zeros(1,N));

%%%% GENERATES 4PSK

% for i = 1: SYMBOL_RATE*UPSAMPLING :N
%     switch data(data_index)
%         case 1
%             for j = 0:SYMBOL_RATE*UPSAMPLING-1
%                 if j < SYMBOL_RATE*UPSAMPLING/2
%                     bpsk(i+j) = carrier(i+j);
%                 else
%                     bpsk(i+j) = 0;
%                 end
%             end
%         case 2
%             for j = 0:SYMBOL_RATE*UPSAMPLING-1
%                 if j < SYMBOL_RATE*UPSAMPLING/2
%                     if carrier(i+j) == 1
%                         bpsk(i+j) = 0;
%                     else
%                         bpsk(i+j) = 1;
%                     end
%                 else
%                     bpsk(i+j) = 0;
%                 end
%             end
%         case 3
%             for j = 0:SYMBOL_RATE*UPSAMPLING-1
%                 if j >= SYMBOL_RATE*UPSAMPLING/2
%                     bpsk(i+j) = carrier(i+j);
%                 else
%                     bpsk(i+j) = 0;
%                 end
%             end
%         case 4
%             for j = 0:SYMBOL_RATE*UPSAMPLING-1
%                 if j >= SYMBOL_RATE*UPSAMPLING/2
%                     if carrier(i+j) == 1
%                         bpsk(i+j) = 0;
%                     else
%                         bpsk(i+j) = 1;
%                     end
%                 else
%                     bpsk(i+j) = 0;
%                 end
%             end
%     end
%     data_index = data_index + 1;
% end


%%%% GENERATES 2PSK -> BPSK

for i = 1: SYMBOL_RATE*UPSAMPLING :N
    for j = 0:SYMBOL_RATE*UPSAMPLING-1
        switch data(data_index)
            case 1
                bpsk(i+j) = carrier(i+j);
            case 2
                for j = 0:SYMBOL_RATE*UPSAMPLING-1
                    if carrier(i+j) == 1
                        bpsk(i+j) = 0;
                    else
                        bpsk(i+j) = 1;
                    end
                end
        end
    end
    data_index = data_index + 1;
end



%% PLOT
figure(1);
subplot(311);   plot(t, carrier);       title('carrier');
subplot(312);   plot(data);             title('data');
subplot(313);   plot(t, bpsk);          title('bpsk');


%% FREQUENCY
res_freq = clk_freq/N;
f = (0:N-1)*res_freq;

carrier_spec = abs(fft(carrier))/N;
bpsk_spec = abs(fft(double(bpsk)))/N;


%% PLOT
figure(2);
subplot(211);   plot(f(1:N/2), mag2db(carrier_spec(1:N/2)));    title('carrier');
subplot(212);   plot(f(1:N/2), mag2db(bpsk_spec(1:N/2)));       title('bpsk');


%% FILTER
bpsk_filtered = filter(BPF_NUM, 1, bpsk); 


%% FREQUENCY
bpsk_filtered_spec = abs(fft(bpsk_filtered))/N;


%% PLOT
figure(3);
subplot(211);   plot(t, bpsk_filtered);
    title('bpsk filtered');
subplot(212);   plot(f, bpsk_filtered_spec);
    title('bpsk filtered - spect');
    

%% DEMODULATION
bpsk_demod = pskdemod(bpsk_filtered, 2);
bpsk_demod_spec = mag2db(abs(fft(bpsk_demod))/N);


%% PLOT
scatterplot(bpsk_demod);	title('bpsk demod');




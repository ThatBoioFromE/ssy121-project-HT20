close all
clear

%% Input parameters
N = 3*12;    % Number of bits to transmit
Rb = 440;   % bit rate [bit/sec]
fs = 4400; % Sampling frequency [Sample/s]
Ts = 1/fs; % Sample time [s/Sample]
fc = 1000;  % Carrier frequency [Hz]
span = 6;   % Pulse-width in symbol times.

rolloff = 0.4; % Roll-off factor (ie. alpha) for Raised Cosine pulse.

%% Transmitter
% Constellation
%constellation = [(1 + 1i), (1 - 1i), (-1 -1i), (-1 + 1i)]/sqrt(2);% Constellation 1 - QPSK/4-QAM
constellation = [sqrt(2), (1 + 1i), sqrt(2)*i, (1 - 1i), -sqrt(2), (-1 -1i), -sqrt(2)*i, (-1 + 1i)]/sqrt(2);  % 8-PAM
grid on;
voronoi(real(constellation), imag(constellation))
title('Constellation plot')
xlabel("Phase"); ylabel("Quadrature");


% Preamble (including initial delay where nothing is tx'd)
delay_time = randi([1,5], 1);       % Generate a random delay between 1 and 20 samples long.
delay = zeros(1, delay_time);
preamble = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]; % length-13 Barker code
preamble_symbs = preamble.*(1+1i)./sqrt(2); % Map preamble to symbols, so that they are 180deg apart! MUST BE MODIFIED IF USING DIFFERENT CONSTELLATION.

% Data source
b = randsrc(1,N,[0 1]); % Message bits

% Bits to Messages
M = length(constellation);      % M unique messages/symbols.
bpsymb = log2(M);       % Bits per symbol
Rsymb = Rb/bpsymb;
Tsymb = 1/Rsymb;
sps   = fs/Rsymb; % Samples per symbol. [Sample/s / Symbol/s = Sample/Symbol] 
m = buffer(b, bpsymb)'; % Group information bits into bits per symbol
m_idx = bi2de(m, 'left-msb')'+1; % Bits to symbol index, msb: the Most Significant Bit

% Message to symbol
message_symbs = constellation(m_idx);    % Map messages to constellation symbols using indices
%symbs = [delay preamble_symbs symbs]; % Prepend variable delay and preamble to message.
symbs = [delay preamble_symbs message_symbs]; % ONLY TRANSMIT PREAMBLE FOR NOW!


% Generate basic pulse
[pulse t] = rtrcpuls(rolloff, Tsymb, fs, span); % Generate basis pulse.
figure
plot(t, pulse); grid on; title('Basic pulse');

% Pulse-modulation
sps = fix(sps); % Convert explicitly to integer !!UGLY HACK!!
x_upsample = upsample(symbs, sps); % Upsample to fs.
s = conv(pulse, x_upsample); % Pulse shaping the symbol-train by convolving with basis pulse.
t_s = (0:length(s)-1).*Ts; % Signal time-axis ie. x[n]*n*Ts.

figure
subplot(221) 
plot(real(s))
grid on; title('Real pulse-modulated baseband');

subplot(223)
stem(real(x_upsample))
grid on; title('Real symbol baseband');

subplot(222) 
plot(imag(s))
grid on; title('Imag pulse-modulated baseband');

subplot(224)
stem(imag(x_upsample))
grid on; title('Imag symbol baseband');

I = real(s);
Q = imag(s);

tx_sig = I.*cos(2*pi*fc*(t_s)) - Q.*sin(2*pi*fc*(t_s));
figure
plot(t_s, tx_sig); 
title('Carrier signal')
grid on;


%% Channel (AWGN)
SNRdB = -5;
s_noisy = awgn(tx_sig, SNRdB, 'measured');
figure; plot(real(s_noisy))
title('Carrier signal with Noise')

%% Receiver
rx_I =  s_noisy.*cos(2*pi*fc*t_s);
rx_Q = -s_noisy.*sin(2*pi*fc*t_s);

figure
subplot(121)
plot(t_s, rx_I)
title('RX-I, mixed')
subplot(122)
plot(t_s, rx_Q)
title('RX-Q, mixed')

rx_I_filtered = lowpass(rx_I, fc/2, fs, 'ImpulseResponse','iir','Steepness',0.85);
rx_Q_filtered = lowpass(rx_Q, fc/2, fs, 'ImpulseResponse','iir','Steepness',0.85);

figure
subplot(121)
plot(t_s, rx_I_filtered)
title('RX-I, filtered')
subplot(122)
plot(t_s, rx_Q_filtered)
title('RX-Q, filtered')

% Matched filter
MF = fliplr(conj(pulse));        %create matched filter impulse response
MF_I_output = filter(MF, 1, rx_I_filtered); 
MF_I_output = MF_I_output(length(MF):end); %remove transient
MF_Q_output = filter(MF, 1, rx_Q_filtered); 
MF_Q_output = MF_Q_output(length(MF):end); %remove transient

figure
subplot(121)
plot(MF_I_output)
hold on
grid on; title('RX-I, matched-filtered')
subplot(122)
plot(MF_Q_output)
grid on; title('RX-Q, matched-filtered')

% Synchroniser (TODO)
rx_symbs = [MF_I_output(1:sps:end)', MF_Q_output(1:sps:end)']; % This assumes PERFECT synchronization.

% Plot sample points on baseband signal.
figure
subplot(121)
plot(MF_I_output)
hold on
stem(1:sps:length(MF_I_output), MF_I_output(1:sps:end))
grid on; title('RX-I, sampled')
subplot(122)
plot(MF_Q_output)
hold on
stem(1:sps:length(MF_Q_output), MF_Q_output(1:sps:end))
grid on; title('RX-Q, sampled')
hold off

% Unmapping to symbols
samples = rx_symbs(:,1) + 1i.*rx_symbs(:,2);
figure
plot(samples, 'o')
axis square; grid on
xlim([-1 1]);
ylim([-1 1]);

% Euclidean distances from each sample to each possible constellation point
% (INEFFICIENT!)
% d = zeros(length(rx_symbs), length(constellation))
% for idx = 1:length(rx_symbs)
%     for jdx = 1:length(constellation)
%         d(idx, jdx) = norm(rx_symbs(idx) - constellation(jdx));
%     end
% end 

d = abs(repmat(samples, 1, length(constellation)) - repmat(constellation, length(samples), 1)).^2;

[nearest_val, nearest_idx] = min(d, [], 2);  % Get index of value with minimum distance to a constellation point. Traverse alone column.

rx_mapped = constellation(nearest_idx); % Choose corresponding symbol.

% Preamble detection
corr = conv(rx_mapped, fliplr(preamble)); % Compare preamble to received symbol vector. Peak at end of preamble.
[~,peak_idx] = max(corr);

figure
plot( abs(corr) )
grid on; title('Preamble detection');
hold on
plot(peak_idx, abs(corr(peak_idx)), 'or')

% Symbols to messages
m_hat = rx_mapped(peak_idx+1:end); % Get symbols after preamble -> actual message.
SER = symerr(message_symbs, m_hat) 

% Messages to bits
b_hat_buffer = de2bi(nearest_idx(peak_idx+1:end)-1, bpsymb, 'left-msb')';
b_hat = b_hat_buffer(:)'; %write as a vector
BER = biterr(b, b_hat) %count of bit errors

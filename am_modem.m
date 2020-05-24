clear; close all; clc;

pkg load signal;

f_sample = 64e3; % sample frequency in samples/s
f_carrier = 1e3; % carrier frequency in Hz
t_symbol = 5e-3; % symbol duration in seconds
m = 1;           % modulation index
N_preamble = 3;  % preable length
snr = 10;        % desired SNR in dB

data = [1 0 1 0 0 1 0 0 1 1 1 0];
tx_data = [ones(1, N_preamble) data]; % Prepend preamble
N_sym = length(tx_data); % Number of symbols to be transmitted

sps = t_symbol * f_sample; % samples per symbol
N = N_sym * sps;           % Total number of required samples

%% Tx
% Generate baseband
M = []; % baseband vector
for i=1:N_sym
    M = [M tx_data(i)*ones(1, sps)];
end

% Scale and shift baseband according to the modulation index
M = M*m + (1-m);

% Generate modulated carrier
n = 0:N-1; % sample indices
rf = cos(2*pi*f_carrier/f_sample*n); % unmodulated carrier
Y = M .* rf;

% Add noise so snr of the Y_noisy is as desired
noise = rand(size(Y));
noise_var = var(Y)/10^(snr/10);
scaled_noise = noise * sqrt(noise_var);

Y_noisy = Y + scaled_noise;

%% Rx
% Simple envelope detector
y_abs = abs(Y_noisy);

taps = fir1(sps, 1/(f_sample*t_symbol), 'low');
bb = filter(taps, 0.5, y_abs);

% Normalize
bbn = bb - min(bb(:));
bbn = bbn ./ max(bbn(:));

% Convert to binary
rx_bin = round(bbn);

% add trailing samples for a single symbol
rx_bin = [rx_bin rx_bin(end)*ones(1, sps-1)];

% Index of first non-zero value
ix = find(rx_bin)(1);

% Undersample the normalized baseband starting from
% the middle of first symbol (assuming no dopler shift)
rx_data = rx_bin(ix+floor(sps/2):sps:end);

%% Bit error ratio
bit_error_ratio = sum(abs(tx_data-rx_data))/length(tx_data)

%% Plots
figure(1)
plot(n, Y_noisy)
title('Modulated carrier');

figure(2)
plot(n, bbn);
title('Normalized demodulated baseband');

figure(3)
stem(tx_data);
hold on;
stem(rx_data);
title('Transmitted and received binary data');
legend('tx', 'rx');

clear; close all; clc;

f_sample = 64e3; % sample frequency in samples/s
f_carrier = 1e3; % carrier frequency in Hz
t_symbol = 5e-3; % symbol duration in seconds
m = 1;           % modulation index

tx_data = [1 0 1 0 0 1 0 0 1 1 1 0];
N_sym = length(tx_data); % Number of symbols to be transmitted

sps = t_symbol * f_sample; % samples per symbol
N = N_sym * sps;           % Total number of required samples

% Generate baseband
M = []; % baseband vector
for i=1:N_sym
    M = [M tx_data(i)*ones(1, sps)];
end

% Scale and shift baseband according to the modulation index
M = M*m + (1-m);

% Generate modulated carrier
n = 0:N-1; % sample indices
Y = M .* cos(2*pi*f_carrier/f_sample*n);

%% Plots
figure(1)
plot(n, Y)
title('Modulated carrier');

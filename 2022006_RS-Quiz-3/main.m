% Aarya Gupta R.NO. -> 2022006
%---------------------------------------------------------%

% Parameters
KEY = 6; 
K = KEY * 1e15; % Chirp rate (Hz/s)
% speed of increasing signal frequency -> Key
B = 3e9; % Bandwidth of the transmitted signal (Hz)
PRI = 520e-6; % Pulse repetition interval (s)
% radar pulses ke beech ka time gap
R1 = KEY * 10; % Target 1 range (m)
% distance of target 1 from radar
R2 = KEY * 20; % Target 2 range (m)
c = 3e8; % Speed of light (m/s)

% Chirp duration T = B/K (chirp signal ki time duration for which it is active)
T = B / K; % 5e-8 seconds (50 ns)

% Sampling frequency (10x bandwidth) 
% (just for clarity, taking it 10 times of bandwith)
fs = 10 * B; % 3e10 Hz
dt = 1 / fs; % 3.333e-11 seconds

% Time vectors
t_tx = -T/2 : dt : T/2 - dt; % Transmitted signal time
t_rx = 0 : dt : 1e-6; % Received signal time (up to 1 µs)

% Transmitted signal
s_tx = cos(pi * K * t_tx.^2);
s_tx(abs(t_tx) > T/2) = 0; % Apply rectangular window

% Received signal (two delayed chirps)
% target se reflected chirps...
tau1 = 2 * R1 / c; % Delay for target 1 (4e-7 s) {time b/w radar and target 1}
tau2 = 2 * R2 / c; % Delay for target 2 (8e-7 s)

% Generate received signal
t1 = t_rx - tau1; % Adjusting time vector for 1st target
s_rx1 = cos(pi * K * t1.^2); % Generated 1st target ke reflected signal
s_rx1(abs(t1) > T/2) = 0; % Rectangular window lagake signal bound kiya

t2 = t_rx - tau2; 
s_rx2 = cos(pi * K * t2.^2); 
s_rx2(abs(t2) > T/2) = 0; 

s_rx = s_rx1 + s_rx2; % Combining the received signals

% Plot transmitted and received signals
figure;
subplot(2,1,1);
plot(t_tx, s_tx, 'Color', 'b');
title('Transmitted Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t_rx, s_rx, 'Color', 'g');
title('Received Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Matched filter output (cross-correlation)
[corr_output, lags] = xcorr(s_rx, s_tx);
time_lags = lags / fs; % Convert lags to seconds
range = (c * time_lags) / 2; % Convert to range (m)

% Find peaks
[pks, locs] = findpeaks(abs(corr_output), 'MinPeakHeight', 0.5*max(abs(corr_output)));

% Plot matched filter output
figure;
plot(range, abs(corr_output), 'Color', 'r');
xlabel('Range (m)');
ylabel('Amplitude');
title('Matched Filter Output');
grid on;
hold on;
plot(range(locs), pks, 'ro');
hold off;

% Display detected ranges
detected_ranges = round(range(locs), 2);
disp('Detected ranges (m):');
disp(detected_ranges);
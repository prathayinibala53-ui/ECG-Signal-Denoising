%% --- 1. Generate Synthetic ECG Signals ---
fs = 200;               % Sampling frequency (Hz)
t = 0:1/fs:10;          % 10-second signals

ECG3 = sin(2*pi*1.2*t) + 0.05*randn(1,length(t));
ECG4 = sin(2*pi*1.1*t) + 0.05*randn(1,length(t));
ECG5 = sin(2*pi*1.3*t) + 0.05*randn(1,length(t));
ECG6 = sin(2*pi*1.25*t) + 0.05*randn(1,length(t));

ECGsignals = {ECG3, ECG4, ECG5, ECG6};
ECGnames = {'ECG3','ECG4','ECG5','ECG6'};

%% --- 2. Define Filters for Pan-Tompkins ---
% Low-pass filter
lp_num = [1 zeros(1,5) -2 zeros(1,5) 1]; lp_den = [1 -2 1];

% High-pass filter
hp_num = [-1/32 zeros(1,15) 1 -1 zeros(1,14) 1/32]; hp_den = [1 -1];

% Derivative filter
deriv_num = (1/8)*[1 2 0 -2 -1]; deriv_den = 1;

% Moving Window Integrator
N = 30; mw_num = (1/N)*ones(1,N); mw_den = 1;

%% --- 3. Process Each ECG Signal ---
figure;
sgtitle('Simplified Pan-Tompkins Stages for Synthetic ECGs');

% Initialize arrays for metrics
total_beats_array = zeros(1,4);
BPM_array = zeros(1,4);
mean_RR_array = zeros(1,4);
SD_RR_array = zeros(1,4);

for k = 1:length(ECGsignals)
    x = ECGsignals{k};
    t_signal = (0:length(x)-1)/fs;

    % Band-pass filter
    x_lp = filter(lp_num, lp_den, x);
    x_bp = filter(hp_num, hp_den, x_lp);

    % Derivative
    x_deriv = filter(deriv_num, deriv_den, x_bp);

    %  Squaring
    x_sq = x_deriv.^2;

    % Moving Window Integration
    x_mwi = filter(mw_num, mw_den, x_sq);

    % Thresholding to detect peaks
    refractory = round(0.2*fs);             % 200 ms
    threshold = 0.3 * max(x_mwi);           % Simple threshold
    [pks, locs] = findpeaks(x_mwi,'MinPeakHeight',threshold,'MinPeakDistance',refractory);
    tp = t_signal(locs);

    % --- 4. Plot stages ---
    subplot(4,5,(k-1)*5+1); plot(t_signal,x); title('Original'); grid on;
    subplot(4,5,(k-1)*5+2); plot(t_signal,x_bp); title('Band-pass'); grid on;
    subplot(4,5,(k-1)*5+3); plot(t_signal,x_deriv); title('Derivative'); grid on;
    subplot(4,5,(k-1)*5+4); plot(t_signal,x_sq); title('Squared'); grid on;
    subplot(4,5,(k-1)*5+5); plot(t_signal,x_mwi); hold on; stem(tp,pks,'filled'); title('MW Integrator + Peaks'); grid on; hold off;

    % --- 5. Compute simple heart metrics ---
    total_beats = length(tp);
    total_time_minutes = t_signal(end)/60;
    BPM = total_beats/total_time_minutes;
    RR_intervals = diff(tp);
    mean_RR = mean(RR_intervals)*1000; % ms
    SD_RR = std(RR_intervals)*1000;    % ms

    total_beats_array(k) = total_beats;
    BPM_array(k) = BPM;
    mean_RR_array(k) = mean_RR;
    SD_RR_array(k) = SD_RR;
end

%% --- 6. Display Results Table ---
result = table(ECGnames', total_beats_array', BPM_array', mean_RR_array', SD_RR_array', ...
    'VariableNames', {'Signal','# Beats','Heart Rate','Mean RR Interval','SD RR Interval'});
disp(result);

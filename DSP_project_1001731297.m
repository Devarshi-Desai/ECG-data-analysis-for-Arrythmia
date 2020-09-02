% DEVARSHI DESAI
% UTA ID: 1001731297
% My last 5 digits of UTA ID is 31297 so 3+1+2+9+7 = 22 hence I chose 122
% as per the guidelines.
%% Question 1
load('122m.mat');
% Loading the data and calculating value from raw to mV
% We need to divide the value by Gain in info file and that is 200.
val = (val - 0)/200;
% We need to subtract the base from the val but since it is 0 we dont have
% to do it.


% The sampling frequency for the data when collected is 360 so we need to
% change the data as such
fs = 360;
t = (0:length(val)-1)/fs;
% Plotting the waveform with altered value
subplot(2,1,1);
plot(t,val); grid on
xlabel('Time (sec)')
ylabel('Amplitude (mV)')
title('Initial ECG plot')
subplot(2,1,2);
%  Total duration of data is 30 mins and that data consists of 650000
%  signals so to get 10 minutes data we divide it by 1800 and multiply by
%  600.
plot(t(1:216667),val(1:216667)); grid on
xlabel('Time (sec)')
ylabel('Amplitude (mV)')
title('ECG plot for 10 mins')
ecg_10min = val(1:216667);

%% Question 2

% Thus populating ecg_10min with first 10 mins of data.

% Iterating through ecg_10min to check for peaks. This can alternatively be
% done by findpeaks function but doing raw calculation gives us more
% control and we can enhance the algorithm if needed later.
beat_count = 0;

for k = 2 : length(ecg_10min)-1
    if (ecg_10min(k) > ecg_10min(k-1)) && (ecg_10min(k) > ecg_10min(k+1)) && (ecg_10min(k) > 0.7)
        % Registering the suitable peak here by incrementing the beat_count.
        beat_count = beat_count + 1;
    end
end

% Displaying the final beat count
disp('Total No of Heartbeats in 10 minutes :');
disp(beat_count);

% Calculating the Beats per minute
N = length(ecg_10min);
dur_sec = N/fs;
dur_min = dur_sec/60;
bpm = beat_count/dur_min;
disp('BPM is :');
disp(bpm);

%% RR interval
interval1 = 1/fs;
limit = interval1 : interval1: length(ecg_10min)/fs;
ten_signal = find(limit==600);
seconds = t(1:ten_signal);
data = ecg_10min(1:ten_signal);
[R,Rt] = findpeaks(data,'MinPeakHeight',0.7);
% Finding the peaks and locations and then finding the absolute rr
% distribution for plotting the histogram.
Rtx = Rt/360;
rrint = abs(diff(Rtx)); 
figure;
stem (Rtx(1:667),rrint); grid on
xlabel('Time(Sec)')
ylabel('R Interval(Sec)')
title('R-R Interval as a function of time')
legend('R Interval')
peak = seconds(Rt);


% Taking the peaks and locations and finding the mean of the 
% distribution for calculating the RR interval time in seconds.
rrint1=mean(diff(Rt));

disp(['RR Interval in secs =', num2str(rrint1/360)]);

%% Interpolating RR intervals
t = linspace(peak(1,1),peak(1,end),length(R)); % equally spaced time interval
d1 = R;
t1 = peak;
d2 = interp1(t1,d1,t);

newinterval = diff(t); % Interpolated interval
%% Histogram calculation of interpolated RR interval.
% Plotting the histogram with the found data.
figure;
hist(newinterval); % Finds histogram
xlabel('interval in seconds')
ylabel('frequency')
title('Histogram of RR interval')

%% Finding Autocorrelation

[acor,lag] = xcorr(rrint,rrint); % Autocorrelation of old interval
[acor1,lag1] = xcorr(newinterval,100); % Autocorrelation of new interval

figure;
plot(lag,acor,'r');
xlabel('Samples lag');
ylabel('Autocorrelation');
title('Autocorrelation of old interval');

figure;
stem(lag,acor,'red');
xlabel('Samples lag');
ylabel('Autocorrelation');
title('Autocorrelation of old interval');

figure;
plot(lag1,acor1,'black');
xlabel('Samples lag');
ylabel('Autocorrelation');
title('Autocorrelation of new interval');
xlim([-100 100]);

figure;
stem(lag1,acor1,'black');
xlabel('Samples lag');
ylabel('Autocorrelation');
title('Autocorrelation of new interval');

%% Power density spectrum
% old interval
figure;
pwelch(rrint,[],[],[],360);

% new interval
figure;
pwelch(newinterval,[],[],[],360);
%% Designing FIR filter (Method 1)
nyquist = fs/2;
cutoff = 50/nyquist; %defining cutoff as 50Hz
f = fir1(100,cutoff);
freqz(f,1);
figure;
new_filtered = filter(f,1,ecg_10min);
subplot(2,1,1);
plot(ecg_10min);
xlabel('time(s)');
ylabel('Amplitude');
title('ECG signal for 10 minutes');
subplot(2,1,2);
plot(new_filtered);
xlabel('time(s)');
ylabel('Amplitude');
title('Filtered 10 minutes signal');

%% Method 2
fc = 50;
tb = 10;
N = (3.3*fs)/tb;
N = 118;
N_mid = N/2;
fc1 = (fc+(tb/2))/fs;
a = ['Corner Frequency (cycles/sample) = ',num2str(fc)];
disp(a)
h_0 = 2*fc;
w_0 = 1;
h(1)=h_0*w_0;
for k=2:N+1
    hd(k)= (2*fc*sin(2*pi*fc*(k-1)))/(2*pi*fc*(k-1));
    wd(k)=0.54+(0.46*cos((2*pi*(k-1))/N)); h(k)=hd(k)*wd(k);
end
figure;
subplot(2,1,1);
plot(ecg_10min);
xlabel('Time(s)');
ylabel('Amplitude');
title('ECG signal for 10 minutes');
subplot(2,1,2);
new_filtered1 = conv(h,ecg_10min);
plot(new_filtered1);
xlabel('Time(s)');
ylabel('Amplitude');
title('Filtered ECG signal');
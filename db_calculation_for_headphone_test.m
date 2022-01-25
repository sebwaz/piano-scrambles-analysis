% Generate reference tone
fs = 48000;
f = 200;
t = linspace(0, 1 - 1/fs, fs);
x1 = sin(2 * pi * f * t);

% Power is sum of squares divided by duration (arbitrary time unit)
P_0 = sum(x1 .^ 2) / 1000;

% Want signal to be -6 dB 
P = P_0 * 10 ^(-6 / 10);

% Figure out the corresponding gain to apply to x1 to get this P
A = sqrt(P / P_0)

% Sanity check
P = sum((A * x1) .^2) / 1000;
10 * log10(P / P_0)
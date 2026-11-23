%% ROBO301: Robotic Sensor Signal Processing
% Scenario: a vibration sensor on a mobile robot; amplitude is in mm.
% All three sinusoidal components are useful sensor motion.
% The unwanted disturbance is additive white Gaussian noise.

clear;
close all;
clc;

%% SETUP: Parameters to change

% Signal and sampling

% 2*cos(4*pi*t)               Amplitude 2         Frequency 2
% 0.8*sin(12*pi*t)            Amplitude 0.8       Frequency 6
% 0.4*cos(20*pi*t)            Amplitude 0.4       Frequency 10

A        = [2, 0.8, 0.4];     % Component amplitudes (mm)
f        = [2, 6, 10];        % Frequencies (Hz): cosine, sine, cosine

Fs       = 500;               % Sampling frequency (Hz)
duration = 10;                % Record length (s)
FsAlias  = 15;                % Low sampling rate for aliasing (Hz)

% Noise
noiseStd   = 1;               % Noise standard deviation (mm)
randomSeed = 5;               % Repeat the same random experiment

% Filter candidates
windowLengths = [3, 5, 9, 15, 25];           % Moving-average windows
alphaValues   = [0.3, 0.5, 0.7, 0.85, 0.95]; % Recursive smoothing values

% Evaluation settings
testFreq           = [f, 30, 60, 100, 180]; % Include useful and higher frequencies
minimumGain        = 0.85;                % Preserve at least 85% of useful amplitude
transientTolerance = 1e-6;                % Threshold for settling
testDuration       = 8;                   % Test length after settling (s)
saveResults        = false;               % Export figures and CSV tables if true

%% SETUP: Check parameter choices
% assert checks whether a condition is true
assert(all(f > 0 & mod(f,1) == 0), 'Use positive integer frequencies.');
assert(numel(A) == 3 && numel(f) == 3, 'The signal uses three components.');
assert(numel(unique(f)) == 3, 'Use three distinct component frequencies.');

% We want to prevent aliasing ( happens when we sample too slowly) so a
% high frequency signal doesnt appear as lower frequency, here 20 Hz is our
% Nyquist rate, with 500 Hz as sampling.
assert(Fs > 2*max(f), 'Fs must exceed twice the highest signal frequency.');
assert(FsAlias < 2*max(f), 'Choose FsAlias below the Nyquist rate.');
assert(all(testFreq > 0 & testFreq < Fs/2), 'Test frequencies must be below Fs/2.');
assert(all(windowLengths >= 1 & mod(windowLengths,1) == 0), 'M must be a positive integer.');
assert(all(alphaValues > 0 & alphaValues < 1), 'Require 0 < alpha < 1.');

%% STAGE 1: Signal model and time-domain properties
% x(t) = 2*cos(4*pi*t) + 0.8*sin(12*pi*t) + 0.4*cos(20*pi*t)

% The greatest common divisor of 2,6,10 is 2 Hz, T0 = 1/f0 = 0.5 seconds 
% So the complete signal repeats every 0.5 seconds

% f0 = gcd(2,6,10) = 2 Hz; T0 = 0.5 seconds.
f0           = gcd(gcd(f(1), f(2)), f(3));
T0           = 1/f0;

% Real, continuous-time, periodic, deterministic; neither even nor odd.
% It is a power signal: finite average power and infinite total energy.
% Its DC/mean is zero. Sum(A) is an amplitude bound, not its exact peak.

w0           = 2*pi*f0;
harmonics    = f/f0; % [1,3,5] are harmonics for 2,6,10 repsectively


signalPower  = sum(A.^2)/2;     % Distinct harmonics are orthogonal

% Here we create the signal , @(t) defines a function whose input is time
% Each term creates one vibration using its ampltidue and frequency
sensorSignal = ...
    @(t) A(1)*cos(2*pi*f(1)*t) + ...
    A(2)*sin(2*pi*f(2)*t)+ ...
    A(3)*cos(2*pi*f(3)*t);

% When we call for example sensorSignal(0), it calculates the signal at t=0
% x(0) = 2cos(0) + 0.8sin(0) + 0.4cos(0) = 2.4 mm, so our clear signal plot
% would start at 2.4 mm


%linspace creates equally spaced values, here T0=0.5, so tFine will contain
% 10,001 time values from 0 to 1 seconds
tFine = linspace(0, 2*T0, 10001); % Dense approximation of continuous time
%contains the signal's value at each of those times
xFine = sensorSignal(tFine);
% We use many points so the plotted curve looks smooth

fprintf('Fundamental frequency = %.2f Hz; period = %.3f s\n', f0, T0);
fprintf('Component amplitudes = %.2f, %.2f, %.2f mm\n', A);
fprintf('Amplitude bound = %.2f mm; dense-grid peak |x| = %.4f mm\n', ...
    sum(A), max(abs(xFine)));
fprintf('Mean = 0 mm; average squared amplitude = %.4f mm^2\n', signalPower);

figure('Name','A: Clean sensor signal');
plot(tFine, xFine, 'b', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Vibration displacement (mm)');
title('Clean periodic sensor signal: two fundamental periods');

%% STAGE 2: Sampling and aliasing
fmax   = max(f);
Ts     = 1/Fs;

t      = (0:round(duration*Fs)-1)/Fs; % Equally spaced sampling times
x      = sensorSignal(t); % Evalutes the clean signal at those times

tAlias = 0:1/FsAlias:2*T0;
xAlias = sensorSignal(tAlias);
% Fold frequencies into [0, FsAlias/2].
aliasFreq     = abs(mod(f + FsAlias/2, FsAlias) - FsAlias/2);
samplingTable = table(f(:), aliasFreq(:), ...
    'VariableNames', {'Original_Hz','Apparent_Hz_at_low_Fs'});
fprintf('\nSampling: fmax = %.1f Hz; Nyquist rate = %.1f samples/s\n', fmax, 2*fmax);
fprintf('Chosen Fs = %.1f Hz; Ts = %.6f s\n', Fs, Ts);
disp(samplingTable);

figure('Name','C: Sampling and aliasing');

subplot(2,1,1);
plot(tFine, xFine, 'k');
hold on;
show = t <= 2*T0;
stem(t(show), x(show), 'b.', 'MarkerSize', 4);
grid on;
xlabel('Time (s)');
ylabel('Displacement (mm)');
title(sprintf('Valid sampling: Fs = %g Hz > 2 fmax', Fs));
legend('Clean continuous-time model','Samples');

subplot(2,1,2);
plot(tFine, xFine, 'k');
hold on;
stem(tAlias, xAlias, 'r', 'filled');
grid on;
xlabel('Time (s)');
ylabel('Displacement (mm)');
title(sprintf('Undersampling: Fs = %g Hz (see sampling table for aliases)', FsAlias));
legend('Clean continuous-time model','Undersampled values');

%% STAGE 3: Reproducible additive white Gaussian noise
% Now we need to add noise to simulate an imperfect sensor

% rng sets the random generator so every run produces the same noise,
% making the results reproducible
rng(randomSeed, 'twister');

% Creates one Gaussian random noise value for each sampling time.
%  noiseStd = 1 would set its std to 1 mm
% If we increase it the fluctations become larger
noise     = noiseStd*randn(size(t));

% Adds the noise to the clean signal, sample by sample
% xm[n] = x[n] + v[n], for example if a clean sample is 2.4 mm and its
% noise value is -0.3 mm , the measured sample becomes 2.1 mm
xMeasured = x + noise;

figure('Name','B: Noise and noisy measurement');

subplot(3,1,1);
plot(t, x, 'k');
grid on;
xlabel('Time (s)');
ylabel('Displacement (mm)');
title('Clean sampled signal');

subplot(3,1,2);
plot(t, noise, 'r');
grid on;
xlabel('Time (s)');
ylabel('Displacement (mm)');
title('Generated AWGN');

subplot(3,1,3);
plot(t, xMeasured, 'b');
grid on;
xlabel('Time (s)');
ylabel('Displacement (mm)');
title('Measured signal = clean signal + noise');

%% STAGE 4: Analytical and numerical CTFS

% x(t) = sum_k Ck*exp(j*k*w0*t)
% A*cos(m*w0*t): C_m = C_-m = A/2.
% A*sin(m*w0*t): C_m = A/(2j), C_-m = -A/(2j).
% Defaults: C(+/-1)=1, C(3)=-0.4j, C(-3)=+0.4j,
% C(+/-5)=0.2; all other coefficients (including C0) are zero.

% contains the harmonic numbers where we store the Fourier series coefficents
% since harmonics = [1,3,5] this gives k = [-7, -6, -5, ..., 0, ..., 5, 6, 7];
k           = -(max(harmonics)+2):(max(harmonics)+2);

% Each harmonic corresponds to freuqnecy kf0. For example k=3 means 3(2) = 6Hz
cAnalytical = complex(zeros(size(k)));

% creeates storage fo these coefficients
% size(k) gets the dimension of k
% zeroes() creates an equally sized array of zeroes
% complex() prepares it to hold complex numbers, such as 0.4j

% Assign the coefficients for the first cosine
% harmonics(1) = 1 and A(1) = 2, k selects the position corresponding to
% the harmonic
cAnalytical(k == harmonics(1)) = A(1)/2;    % Positive harmonic cAnalytical(k == 1) = 1
cAnalytical(k == -harmonics(1)) = A(1)/2;   % Negative harmonic cAnalytical(k == -1) = 1

% ^ Doesnt have imginary component because it is cos


% li is matlab's imaginary unit, both magnitudes are 0.4 but with oppositephases
% Second component: sine at harmonic 3
cAnalytical(k == harmonics(2)) = cAnalytical(k == harmonics(2)) + A(2)/(2*1i);
cAnalytical(k == -harmonics(2)) = cAnalytical(k == -harmonics(2)) - A(2)/(2*1i);

% Third component: cosine at harmonic 5
cAnalytical(k == harmonics(3)) = cAnalytical(k == harmonics(3)) + A(3)/2;
cAnalytical(k == -harmonics(3)) = cAnalytical(k == -harmonics(3)) + A(3)/2;

% Prepare time values and storage
tPeriod    = linspace(0, T0, 10001); %contains 10,001 points over one period, from 0 to 0.5 s.
cNumerical = complex(zeros(size(k))); % stores the coefficients MATLAB calculates.

%This repeats the integral for every harmonic:
for index = 1:numel(k)
    cNumerical(index) = trapz(tPeriod, ...
        sensorSignal(tPeriod).*exp(-1i*k(index)*w0*tPeriod))/T0;
end
%   sensorSignal(tPeriod) gives x(t)
%   exp(...) gives e^-jk\omega0t
%   .* multiplies their values point by point.
%   trapz approximates the integral using trapezoids.
%   /T0 divides by the period.

ctfsTable = table(k(:), (k(:)*f0), cAnalytical(:), cNumerical(:), ...
    abs(cAnalytical(:)-cNumerical(:)), 'VariableNames', ...
    {'Harmonic_k','Frequency_Hz','C_analytical','C_numerical','Absolute_error'});
disp(ctfsTable);
assert(max(abs(cAnalytical-cNumerical)) < 1e-8, 'CTFS verification failed.');


phaseC = angle(cAnalytical); % angle calculates each coefficient's phase in radians.
phaseC(abs(cAnalytical) < 1e-10) = NaN; % Zero coefficients have undefined phase

figure('Name','D: CTFS spectra');

subplot(2,1,1);

% Draws the magnitude spectrum, plots the harmonic numbers on horizontal
% and coefficient magnitues on vertical
stem(k, abs(cAnalytical), 'filled');

grid on;
xlabel('Harmonic index k');
ylabel('|C_k| (mm)');
title('CTFS magnitude spectrum');

subplot(2,1,2);

stem(k, phaseC, 'filled'); % plots phase versus harmonic number.

grid on;
xlabel('Harmonic index k');
ylabel('Phase (rad)');
title('CTFS phase: only nonzero coefficients');

%% STAGE 5: Compare filter parameters using the SAME noisy record
% MATLAB filter(b,a,input) uses zero initial conditions by default.
% MA: y[n] = (1/M)*sum_{r=0}^{M-1} xMeasured[n-r]
% Recursive: y[n] = alpha*y[n-1] + (1-alpha)*xMeasured[n]
% Use one common evaluation interval, even during parameter investigation.
transientSamples = max(max(windowLengths)-1, ...
    ceil(log(transientTolerance)/log(max(alphaValues))));
assert(numel(t)-transientSamples >= Fs*T0, 'Increase duration: insufficient steady-state data.');
steadySamples    = (transientSamples+1):numel(t);
errorMA          = zeros(size(windowLengths));
errorRecursive   = zeros(size(alphaValues));
minGainMA        = zeros(size(windowLengths));
minGainRecursive = zeros(size(alphaValues));

for index = 1:numel(windowLengths)
    trialM = windowLengths(index);
    yTrial = filter(ones(1,trialM)/trialM, 1, xMeasured);
    errorMA(index) = std(yTrial(steadySamples)-x(steadySamples));
    gains = zeros(size(f));
    for frequencyIndex = 1:numel(f)
        gains(frequencyIndex) = abs(mean(exp(-1i*2*pi*f(frequencyIndex)/Fs*(0:trialM-1))));
    end
    minGainMA(index) = min(gains);
end

for index = 1:numel(alphaValues)
    trialAlpha = alphaValues(index);
    yTrial     = filter(1-trialAlpha, [1,-trialAlpha], xMeasured);
    errorRecursive(index) = std(yTrial(steadySamples)-x(steadySamples));
    gains = abs((1-trialAlpha)./(1-trialAlpha*exp(-1i*2*pi*f/Fs)));
    minGainRecursive(index) = min(gains);
end

maParameterTable = table(windowLengths(:), errorMA(:), minGainMA(:), ...
    'VariableNames', {'M','Error_std_mm','Minimum_useful_gain'});
recursiveParameterTable = table(alphaValues(:), errorRecursive(:), minGainRecursive(:), ...
    'VariableNames', {'Alpha','Error_std_mm','Minimum_useful_gain'});
disp(maParameterTable);
disp(recursiveParameterTable);

%% STAGE 5 (continued): Choose and apply the filters
% Choose the lowest error while preserving the required useful amplitude.
scoreMA = errorMA;
scoreMA(minGainMA < minimumGain) = Inf;
scoreRecursive = errorRecursive;
scoreRecursive(minGainRecursive < minimumGain) = Inf;
assert(any(isfinite(scoreMA)) && any(isfinite(scoreRecursive)), ...
    'No feasible candidate: add smaller M/alpha values or revise the gain target.');
[~, bestMA] = min(scoreMA);
[~, bestRecursive] = min(scoreRecursive);
M          = windowLengths(bestMA);
alpha      = alphaValues(bestRecursive);
bMA        = ones(1,M)/M;
yMA        = filter(bMA, 1, xMeasured);
yRecursive = filter(1-alpha, [1,-alpha], xMeasured);
fprintf('\nSelected M = %g; alpha = %.2f; common discarded samples = %d\n', ...
    M, alpha, transientSamples);

figure('Name','Parameter investigation');

subplot(2,1,1);
plot(windowLengths, errorMA, '-o');
grid on;
xlabel('Window length M (samples)');
ylabel('Error std (mm)');
title('Moving-average parameter investigation');

subplot(2,1,2);
plot(alphaValues, errorRecursive, '-o');
grid on;
xlabel('Smoothing parameter alpha');
ylabel('Error std (mm)');
title('Recursive parameter investigation');

%% STAGE 6: LTI, impulse response, and measured amplitude ratios
% Both fixed-parameter filters are causal, stable LTI systems.
% y[n] = xMeasured[n] * h[n], where * denotes convolution.
% hMA[n] = 1/M for 0<=n<M, and 0 otherwise.
% hRecursive[n] = (1-alpha)*alpha^n*u[n].
% HMA(e^jw) = (1/M)*sum_{r=0}^{M-1} exp(-j*w*r)
% HRecursive(e^jw) = (1-alpha)/(1-alpha*exp(-j*w)).
nImpulse   = 0:max(M, ceil(log(1e-3)/log(alpha)));
hMA        = double(nImpulse < M)/M;
hRecursive = (1-alpha)*alpha.^nImpulse;

figure('Name','E: Impulse responses');

subplot(2,1,1);
stem(nImpulse, hMA, 'filled');
grid on;
xlabel('Sample n');
ylabel('h_{MA}[n] (dimensionless)');
title(sprintf('Moving average, M = %g',M));

subplot(2,1,2);
stem(nImpulse, hRecursive, 'filled');
grid on;
xlabel('Sample n');
ylabel('h_R[n] (dimensionless)');
title(sprintf('Recursive smoothing, alpha = %.2f',alpha));

%% STAGE 6 (continued): Test one frequency at a time
% Use a unit-amplitude sinusoid with no added noise.
% Measure amplitude by fitting a sine and cosine to the settled output:
%     output = c*cos(2*pi*f*t) + s*sin(2*pi*f*t) + DC
%     amplitude = sqrt(c^2 + s^2)
% This accounts for phase delay and avoids missing peaks between samples.
tTest           = (0:(transientSamples+round(testDuration*Fs)-1))/Fs;
testKeep        = (transientSamples+1):numel(tTest);
ratioMA         = zeros(size(testFreq));
ratioRecursive  = ratioMA;
theoryMA        = ratioMA;
theoryRecursive = ratioMA;

for index = 1:numel(testFreq)
    % Generate and filter one test signal.
    w             = 2*pi*testFreq(index)/Fs;
    testInput     = cos(2*pi*testFreq(index)*tTest);
    testMA        = filter(bMA, 1, testInput);
    testRecursive = filter(1-alpha, [1,-alpha], testInput);

    % Keep settled samples and fit the input and output waveforms.
    timeFit       = tTest(testKeep).';
    fitMatrix     = [cos(2*pi*testFreq(index)*timeFit), ...
        sin(2*pi*testFreq(index)*timeFit), ones(numel(timeFit),1)];
    inputFit       = fitMatrix\testInput(testKeep).';
    maFit          = fitMatrix\testMA(testKeep).';
    recursiveFit   = fitMatrix\testRecursive(testKeep).';
    inputAmplitude = hypot(inputFit(1), inputFit(2));

    % Compare measured amplitude ratios with the theoretical responses.
    ratioMA(index) = hypot(maFit(1),maFit(2))/inputAmplitude;
    ratioRecursive(index) = hypot(recursiveFit(1),recursiveFit(2))/inputAmplitude;
    theoryMA(index) = abs(sum(bMA.*exp(-1i*w*(0:M-1))));
    theoryRecursive(index) = abs((1-alpha)/(1-alpha*exp(-1i*w)));
end

responseTable = table(testFreq(:), ratioMA(:), theoryMA(:), ...
    ratioRecursive(:), theoryRecursive(:), 'VariableNames', ...
    {'Frequency_Hz','MA_measured','MA_theory','Recursive_measured','Recursive_theory'});
disp(responseTable);
assert(max(abs(ratioMA-theoryMA)) < 1e-5, 'MA amplitude check failed.');
assert(max(abs(ratioRecursive-theoryRecursive)) < 1e-5, 'Recursive amplitude check failed.');

figure('Name','F: Measured amplitude ratios');
plot(testFreq, ratioMA, '-o', testFreq, ratioRecursive, '-s', 'LineWidth',1.2);
grid on;
xlabel('Test frequency (Hz)');
ylabel('Output / input amplitude (dimensionless)');
title('Separate noise-free sinusoidal tests after transients');
legend('Moving average','Recursive smoothing');
% MA has a sinc-like response with zeros and sidelobes; attenuation need
% not increase monotonically with frequency. Recursive magnitude decreases.

%% STAGE 7: Time-domain and numerical performance comparison
% Errors use the SAME sample indices and original time reference.
% No delay compensation: phase delay and attenuation count as distortion.
errors = [xMeasured(steadySamples)-x(steadySamples); ...
    yMA(steadySamples)-x(steadySamples); yRecursive(steadySamples)-x(steadySamples)];
errorStd = std(errors, 0, 2);
% Filter the known generated noise separately to distinguish pure noise
% suppression from total error (which also contains useful-signal distortion).
noiseMA        = filter(bMA, 1, noise);
noiseRecursive = filter(1-alpha, [1,-alpha], noise);
outputNoiseStd = [std(noise(steadySamples)); std(noiseMA(steadySamples)); ...
    std(noiseRecursive(steadySamples))];
noiseReductionPct = 100*(1-outputNoiseStd/outputNoiseStd(1));
Method            = {'Noisy input'; 'Moving average'; 'Recursive smoothing'};

comparisonTable = table(Method, errorStd, outputNoiseStd, noiseReductionPct, ...
    [1; min(ratioMA(1:3)); min(ratioRecursive(1:3))], ...
    'VariableNames', {'Method','Error_std_mm','Noise_only_std_mm', ...
    'Noise_std_reduction_percent','Minimum_useful_gain'});
disp(comparisonTable);

figure('Name','G: Sensor output comparison');

subplot(2,1,1);
plot(t, xMeasured, 'Color',[0.8,0.8,0.8]);
hold on;
plot(t, x, 'k', t, yMA, 'b', t, yRecursive, 'r');
grid on;
xlabel('Time (s)');
ylabel('Displacement (mm)');
title('Complete record, including startup');
legend('Noisy input','Clean reference','Moving average','Recursive smoothing');

subplot(2,1,2);
plot(t, x, 'k', t, yMA, 'b', t, yRecursive, 'r', 'LineWidth',1.1);
grid on;
xlim([t(transientSamples+1), min(t(transientSamples+1)+2*T0,t(end))]);
xlabel('Time (s)');
ylabel('Displacement (mm)');
title('Steady-state detail: smoothing, attenuation, and delay');
legend('Clean reference','Moving average','Recursive smoothing');

figure('Name','I: Error comparison');
bar(errorStd);
grid on;
set(gca,'XTick',1:3,'XTickLabel',Method);
ylabel('Error standard deviation (mm)');
title('Common interval, initial transients excluded');

%% STAGE 8: Evidence-based engineering recommendation
% Both selected candidates pass the useful-frequency amplitude target.
% Use raw error std to rank these candidates; inspect the plots and response
% table as well. This is an experiment-specific decision, not a universal one.
if errorStd(2) <= errorStd(3)
    preferredFilter = 'Moving average';
else
    preferredFilter = 'Recursive smoothing';
end
fprintf('\nPreferred filter for this record and tested candidates: %s\n', preferredFilter);
fprintf('Both selected filters meet the %.0f%% useful-amplitude target.\n', 100*minimumGain);
fprintf('MA error std = %.4f mm; recursive error std = %.4f mm.\n', errorStd(2),errorStd(3));
fprintf('Use the response table and time-domain plots to discuss distortion.\n');
fprintf('MA: finite memory, linear phase, delay = %.1f samples (%.4f s).\n', ...
    (M-1)/2, (M-1)/(2*Fs));
fprintf('Recursive: one previous output, constant work per sample, frequency-dependent delay.\n');
fprintf('The built-in MA uses an FIR of length M; a running-sum MA can also use constant work.\n');
fprintf('Stronger smoothing removes more noise but can distort useful motion.\n');
fprintf('Selection uses one seeded record; other noise records can change the ranking.\n');

%% EXPORT: Save figures and tables (optional)
% Set saveResults = true in SETUP to enable this section.
if saveResults
    scriptFolder = fileparts(mfilename('fullpath'));
    outputFolder = fullfile(scriptFolder, 'project_results');
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    writetable(samplingTable, fullfile(outputFolder,'sampling.csv'));
    writetable(ctfsTable, fullfile(outputFolder,'ctfs.csv'));
    writetable(maParameterTable, fullfile(outputFolder,'ma_parameters.csv'));
    writetable(recursiveParameterTable, fullfile(outputFolder,'recursive_parameters.csv'));
    writetable(responseTable, fullfile(outputFolder,'amplitude_ratios.csv'));
    writetable(comparisonTable, fullfile(outputFolder,'filter_comparison.csv'));
    figures = findall(0,'Type','figure');
    for index = 1:numel(figures)
        saveas(figures(index), fullfile(outputFolder,sprintf('Figure_%02d.png',figures(index).Number)));
    end
    fprintf('Figures and CSV tables saved in: %s\n', outputFolder);
end

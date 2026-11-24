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

%------------------------------------------------------------------------------------

% Filter candidates
windowLengths = [3, 5, 9, 15, 25];           % Moving-average windows
alphaValues   = [0.3, 0.5, 0.7, 0.85, 0.95]; % Recursive smoothing values

%------------------------------------------------------------------------------------

% Evaluation settings
testFreq           = [f, 30, 60, 100, 180]; % Include useful and higher frequencies
minimumGain        = 0.85;                % Our design target: preserve at least 85% of useful amplitude
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
harmonics    = f/f0; % [1,3,5] are harmonics for 2,6,10 respectively


signalPower  = sum(A.^2)/2;     % Distinct harmonics are orthogonal

% Here we create the signal , @(t) defines a function whose input is time
% Each term creates one vibration using its amplitude and frequency
sensorSignal = ...
    @(t) A(1)*cos(2*pi*f(1)*t) + ...
    A(2)*sin(2*pi*f(2)*t)+ ...
    A(3)*cos(2*pi*f(3)*t);

% When we call for example sensorSignal(0), it calculates the signal at t=0
% x(0) = 2cos(0) + 0.8sin(0) + 0.4cos(0) = 2.4 mm, so our clean signal plot
% would start at 2.4 mm


%linspace creates equally spaced values, here T0=0.5, so tFine will contain
% 10,001 time values from 0 to 1 second
tFine = linspace(0, 2*T0, 10001); % Dense approximation of continuous time
%contains the signal's value at each of those times
xFine = sensorSignal(tFine);
% We use many points so the plotted curve looks smooth


%------------------------------------------------------------------------------------

% Outputs

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
x      = sensorSignal(t); % Evaluates the clean signal at those times

tAlias = 0:1/FsAlias:2*T0;
xAlias = sensorSignal(tAlias);
% Fold frequencies into [0, FsAlias/2].
aliasFreq     = abs(mod(f + FsAlias/2, FsAlias) - FsAlias/2);
samplingTable = table(f(:), aliasFreq(:), ...
    'VariableNames', {'Original_Hz','Apparent_Hz_at_low_Fs'});


%------------------------------------------------------------------------------------

% Figures & Output

fprintf('\nSampling: fmax = %.1f Hz; Nyquist rate = %.1f samples/s\n', fmax, 2*fmax);
fprintf('Chosen Fs = %.1f Hz; Ts = %.6f s\n', Fs, Ts);
disp(samplingTable);

figure('Name','C: Sampling and aliasing');

subplot(2,1,1);
plot(tFine, xFine, 'y');
hold on;
show = t <= 2*T0;
stem(t(show), x(show), 'b.', 'MarkerSize', 4);
grid on;
xlabel('Time (s)');
ylabel('Displacement (mm)');
title(sprintf('Valid sampling: Fs = %g Hz > 2 fmax', Fs));
legend('Clean continuous-time model','Samples');

subplot(2,1,2);
plot(tFine, xFine, 'y');
hold on;
stem(tAlias, xAlias, 'r', 'filled');
grid on;
xlabel('Time (s)');
ylabel('Displacement (mm)');
title(sprintf('Undersampling: Fs = %g Hz (see sampling table for aliases)', FsAlias));
legend('Clean continuous-time model','Undersampled values');

%------------------------------------------------------------------------------------

%% STAGE 3: Reproducible additive white Gaussian noise
% Now we need to add noise to simulate an imperfect sensor

% rng sets the random generator so every run produces the same noise,
% making the results reproducible
rng(randomSeed, 'twister');

% Creates one Gaussian random noise value for each sampling time.
% noiseStd = 1 sets the noise distribution standard deviation to 1 mm.
% The standard deviation of a finite generated record will differ slightly.
% If we increase it the fluctuations become larger
noise     = noiseStd*randn(size(t));

% Adds the noise to the clean signal, sample by sample
% xm[n] = x[n] + v[n], for example if a clean sample is 2.4 mm and its
% noise value is -0.3 mm , the measured sample becomes 2.1 mm
xMeasured = x + noise; % Noisy signal

%------------------------------------------------------------------------------------

% Figures

figure('Name','B: Noise and noisy measurement');

subplot(3,1,1);
plot(t, x, 'y');
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

% contains the harmonic numbers where we store the Fourier series coefficients
% since harmonics = [1,3,5] this gives k = [-7, -6, -5, ..., 0, ..., 5, 6, 7];
k           = -(max(harmonics)+2):(max(harmonics)+2);

% Each harmonic corresponds to frequency kf0. For example k=3 means 3(2) = 6Hz
cAnalytical = complex(zeros(size(k)));

% Creates storage for these coefficients.
% size(k) gets the dimension of k
% zeros() creates an equally sized array of zeros.
% complex() prepares it to hold complex numbers, such as 0.4j

% Assign the coefficients for the first cosine
% harmonics(1) = 1 and A(1) = 2, k selects the position corresponding to
% the harmonic
cAnalytical(k == harmonics(1)) = A(1)/2;    % Positive harmonic cAnalytical(k == 1) = 1
cAnalytical(k == -harmonics(1)) = A(1)/2;   % Negative harmonic cAnalytical(k == -1) = 1

% This cosine has no phase shift, so its coefficients are real.


% 1i is MATLAB's imaginary unit.
% The sine coefficients have magnitude 0.4 and opposite phases.
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

%------------------------------------------------------------------------------------

% Figures

figure('Name','D: CTFS spectra');

subplot(2,1,1);

% Draws the magnitude spectrum, plots the harmonic numbers on horizontal
% and coefficient magnitudes on vertical
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
% Moving Average: y[n] = (1/M)*sum_{r=0}^{M-1} xMeasured[n-r]

% Recursive: y[n] = alpha*y[n-1] + (1-alpha)*xMeasured[n]

% Use one common evaluation interval, even during parameter investigation.


% We need to compare different M and alpha values,  both filters start with
% zero initial conditions so their first outputs include a startup transient
% (temporary adjustment before normal filtering behavior settles)


% These lines choose which samples to evaluate

% Allows enough settling time for all tested filters
% Here it excludes the first 270 samples/ 0.54 seconds
transientSamples = max(max(windowLengths)-1, ...
    ceil(log(transientTolerance)/log(max(alphaValues))));

% Checks that enough data remains after skipping the startup transient
assert(numel(t)-transientSamples >= Fs*T0, 'Increase duration: insufficient steady-state data.');
% numel(t) = 5000 total samples; transientSamples = 270 samples to skip
% Fs*T0 = 250 samples (one complete signal period )

% selects the sample positions to use
steadySamples    = (transientSamples+1):numel(t); %271:5000


errorMA          = zeros(size(windowLengths));
errorRecursive   = zeros(size(alphaValues));
minGainMA        = zeros(size(windowLengths));
minGainRecursive = zeros(size(alphaValues));

% Tests each moving average window length  ([3,5,9,15,25])

for index = 1:numel(windowLengths)
    trialM = windowLengths(index); % picks the window length to test

    yTrial = filter(ones(1,trialM)/trialM, 1, xMeasured); % Store its filtered output.

    errorMA(index) = std(yTrial(steadySamples)-x(steadySamples)); % Error standard deviation over the settled samples.

    % Calculate the moving-average filter's theoretical amplitude gain at
    % each useful frequency.

    % Gain = output amplitude/input amplitude
    % Gain  1 means amplitude is fully preserved,
    %       0.95 -> 95% is preserved
    %       0.5 means only 50% remains
    gains = zeros(size(f));
    for frequencyIndex = 1:numel(f)
        gains(frequencyIndex) = abs(mean(exp(-1i*2*pi*f(frequencyIndex)/Fs*(0:trialM-1))));
    end
    minGainMA(index) = min(gains); % Store the lowest gain among the three useful frequencies.
    % gain tells us how much of a frequency component's amplitude survives filtering
end


% tests each recursive smoothing value
% Follow the same process as the moving-average loop.
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

% Selects the best eligible settings
% minimumGain = 0.85

scoreMA = errorMA;
scoreMA(minGainMA < minimumGain) = Inf;

scoreRecursive = errorRecursive;
scoreRecursive(minGainRecursive < minimumGain) = Inf;

assert(any(isfinite(scoreMA)) && any(isfinite(scoreRecursive)), ...
    'No feasible candidate: add smaller M/alpha values or revise the gain target.');

% Any candidate that preserves less than 85% gets a score of infinity, so it
% will not be selected as the smallest score.
% ~ ignores the score itself
[~, bestMA] = min(scoreMA);
[~, bestRecursive] = min(scoreRecursive);

% For the current settings and noise realization: M = 9, alpha = 0.7.
M          = windowLengths(bestMA);
alpha      = alphaValues(bestRecursive);

% Apply selected filters to the noisy signal
bMA        = ones(1,M)/M; % create M equal weights
yMA        = filter(bMA, 1, xMeasured); % Apply those weights to the noisy measurements.
% With the selected alpha = 0.7: 70% previous output, 30% current measurement.
yRecursive = filter(1-alpha, [1,-alpha], xMeasured);
% y[n] - alpha y[n-1]=(1-alpha)xm[n]
% For alpha = 0.7, combine 70% of the previous filtered output
% with 30% of the current noisy measurement.
% The higher alpha, the stronger smoothing, but with slower response


fprintf('\nSelected M = %g; alpha = %.2f; common discarded samples = %d\n', ...
    M, alpha, transientSamples);

%------------------------------------------------------------------------------------

% Figures
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

%------------------------------------------------------------------------------------

%% STAGE 6: LTI, impulse response, and measured amplitude ratios
% Both fixed-parameter filters are causal, stable LTI systems.
% y[n] = xMeasured[n] * h[n], where * denotes convolution.
% hMA[n] = 1/M for 0<=n<M, and 0 otherwise.
% hRecursive[n] = (1-alpha)*alpha^n*u[n].
% HMA(e^jw) = (1/M)*sum_{r=0}^{M-1} exp(-j*w*r)
% HRecursive(e^jw) = (1-alpha)/(1-alpha*exp(-j*w)).

% Check how filter respond to different frequencies
% We treat each filter as an LTI system:

% Linear: scaling or adding inputs scales or adds their outputs
% Time invariant: delaying the input delays the output by the same amount

% impulse = [1,0,0,0...]

nImpulse   = 0:max(M, ceil(log(1e-3)/log(alpha)));

hMA        = double(nImpulse < M)/M; % [1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 0, ...];
hRecursive = (1-alpha)*alpha.^nImpulse;

%------------------------------------------------------------------------------------
%Figures

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

% We measure input and output amplitudes after startup,
% then divide them to determine how strongly each frequency is preserved or attenuated.
for index = 1:numel(testFreq)


    % Generate and filter one test signal.
    w             = 2*pi*testFreq(index)/Fs;

    testInput     = cos(2*pi*testFreq(index)*tTest);

    testMA        = filter(bMA, 1, testInput);

    testRecursive = filter(1-alpha, [1,-alpha], testInput);

    % Keep settled samples and fit the input and output waveforms.
    timeFit       = tTest(testKeep).';  %testKeep selects samples after startup
    fitMatrix     = [cos(2*pi*testFreq(index)*timeFit), ...
        sin(2*pi*testFreq(index)*timeFit), ones(numel(timeFit),1)]; %contains a cosine, a sine, and a constant
    % we fit the settled output to
    % y(t) = ccos(2pi f t) + ssin (2pi f t ) + d
    % Using both sine and cosine lets us measure the amplitude even when the filter introduces phase delay.
    inputFit       = fitMatrix\testInput(testKeep).';
    maFit          = fitMatrix\testMA(testKeep).';
    recursiveFit   = fitMatrix\testRecursive(testKeep).';

     % Calculate the amplitude.
    inputAmplitude = hypot(inputFit(1), inputFit(2));

    % Compare measured amplitude ratios with the theoretical responses.

    ratioMA(index) = hypot(maFit(1),maFit(2))/inputAmplitude;


    ratioRecursive(index) = hypot(recursiveFit(1),recursiveFit(2))/inputAmplitude;

    % calculates the expected amplitude ratios from theory
    theoryMA(index) = abs(sum(bMA.*exp(-1i*w*(0:M-1))));
    theoryRecursive(index) = abs((1-alpha)/(1-alpha*exp(-1i*w)));


end

responseTable = table(testFreq(:), ratioMA(:), theoryMA(:), ...
    ratioRecursive(:), theoryRecursive(:), 'VariableNames', ...
    {'Frequency_Hz','MA_measured','MA_theory','Recursive_measured','Recursive_theory'});
disp(responseTable);

%This checks that the largest difference is below 0.00001. A similar check is performed for the recursive filter.
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

% Compare each signal with the clean reference x.
% Use the same settled samples for all three signals.
% Do not shift the outputs to remove delay: delay counts as distortion.
errors = [xMeasured(steadySamples)-x(steadySamples); ...
          yMA(steadySamples)-x(steadySamples); ...
          yRecursive(steadySamples)-x(steadySamples)];

% Calculate one error standard deviation for each row.
% This measures error variation; a constant bias is not detected by std.
% Also inspect the filtered outputs against the clean-reference plots.
% 0 = default sample standard deviation; 2 = calculate across columns.
errorStd = std(errors, 0, 2);


%% Measure noise reduction separately

% Filter only the known noise to measure how much noise survives.
% Unlike total error, this does not include distortion of useful motion.
noiseMA        = filter(bMA, 1, noise);
noiseRecursive = filter(1-alpha, [1, -alpha], noise);

% Store the noise standard deviation before and after filtering.
% Row order: original noise, moving average, recursive smoothing.
outputNoiseStd = [std(noise(steadySamples)); ...
                  std(noiseMA(steadySamples)); ...
                  std(noiseRecursive(steadySamples))];

% Percentage reduction relative to the original noise standard deviation.
% Example: reducing the standard deviation from 1 to 0.4 gives 60%.
noiseReductionPct = 100*(1-outputNoiseStd/outputNoiseStd(1));

% Labels for the table and bar chart, in the same row order.
Method = {'Noisy input'; 'Moving average'; 'Recursive smoothing'};


%% Display the numerical comparison table

% Include error, remaining noise, noise reduction, and useful-signal gain.
% The first three test frequencies are the useful frequencies: 2, 6, 10 Hz.
% min(...) selects the lowest gain among those three frequencies.
% The unfiltered input has gain 1 because no filter has been applied.
comparisonTable = table( ...
    Method, errorStd, outputNoiseStd, noiseReductionPct, ...
    [1; min(ratioMA(1:3)); min(ratioRecursive(1:3))], ...
    'VariableNames', {'Method', 'Error_std_mm', 'Noise_only_std_mm', ...
                     'Noise_std_reduction_percent', 'Minimum_useful_gain'});

% Print the table in the Command Window.
disp(comparisonTable);


%% Plot the signals over the full recording

figure('Name', 'G: Sensor output comparison');

% Create the upper plot in a layout with two rows and one column.
subplot(2, 1, 1);

% Draw the noisy measurement in light grey.
plot(t, xMeasured, 'Color', [0.8, 0.8, 0.8]);

% Keep the existing curve when adding the other signals.
hold on;

% Overlay the clean reference and both filtered outputs.
% k = black, b = blue, r = red.
plot(t, x, 'y', t, yMA, 'b', t, yRecursive, 'r');

grid on;
xlabel('Time (s)');
ylabel('Displacement (mm)');
title('Complete record, including startup');
legend('Noisy input', 'Clean reference', ...
       'Moving average', 'Recursive smoothing');


%% Zoom in on the settled outputs

% Create the lower plot.
subplot(2, 1, 2);

% Plot the reference and filtered signals with slightly thicker lines.
plot(t, x, 'y', t, yMA, 'b', t, yRecursive, 'r', 'LineWidth', 1.1);

grid on;

% Show two signal periods, starting after the excluded startup samples.
% min(..., t(end)) prevents the view from extending past the recording.
% xlim only changes the displayed interval; it does not change the data.
xlim([t(transientSamples+1), ...
      min(t(transientSamples+1)+2*T0, t(end))]);

xlabel('Time (s)');
ylabel('Displacement (mm)');
title('Steady-state detail: smoothing, attenuation, and delay');
legend('Clean reference', 'Moving average', 'Recursive smoothing');


%% Plot the error standard deviations

figure('Name', 'I: Error comparison');

% Draw one bar for each signal. Lower bars mean less error variation.
bar(errorStd);

grid on;

% Place the method names below the three bars.
set(gca, 'XTick', 1:3, 'XTickLabel', Method);

ylabel('Error standard deviation (mm)');
title('Common interval, initial transients excluded');


%% STAGE 8: Evidence-based engineering recommendation

% Both selected filters already meet our chosen design requirement to preserve
% at least 85% of the amplitude at every useful frequency.
%
% Now compare their error standard deviations:
% errorStd(1) = noisy input
% errorStd(2) = moving-average output
% errorStd(3) = recursive output

% Choose the filter with the smaller error standard deviation.
% If the values are equal, this rule selects the moving average.
if errorStd(2) <= errorStd(3)
    preferredFilter = 'Moving average';
else
    preferredFilter = 'Recursive smoothing';
end

% Display the recommendation for this signal and the tested settings.
% This does not mean the chosen filter is best for every application.
fprintf('\nPreferred filter for this record and tested candidates: %s\n', ...
    preferredFilter);

% Convert the amplitude-preservation target to a percentage.
fprintf('Both selected filters meet the %.0f%% useful-amplitude target.\n', ...
    100*minimumGain);

% Display the numerical evidence supporting the decision.
fprintf('MA error std = %.4f mm; recursive error std = %.4f mm.\n', ...
    errorStd(2), errorStd(3));

% Error standard deviation alone does not describe every filter effect.
% Also inspect amplitude preservation, visible distortion, and delay.
fprintf('Use the response table and time-domain plots to discuss distortion.\n');

% A causal moving average introduces a delay of (M-1)/2 samples.
% Divide by Fs to convert the delay from samples to seconds.
% For M = 9 and Fs = 500: delay = 4 samples = 0.008 seconds.
fprintf('MA: finite memory, linear phase, delay = %.1f samples (%.4f s).\n', ...
    (M-1)/2, (M-1)/(2*Fs));

% Recursive smoothing needs only the previous output as its stored state.
% Its phase delay depends on frequency.
fprintf(['Recursive: one previous output, constant work per sample, ' ...
         'frequency-dependent delay.\n']);

% The direct moving-average implementation uses M coefficients.
% A running-sum implementation can reduce the work needed per sample.
fprintf(['The built-in MA uses an FIR of length M; ' ...
         'a running-sum MA can also use constant work.\n']);

% State the main engineering trade-off.
fprintf(['Stronger smoothing removes more noise ' ...
         'but can distort useful motion.\n']);

% State the limitation: results are based on one fixed noise realization.
fprintf(['Selection uses one seeded record; ' ...
         'other noise records can change the ranking.\n']);

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

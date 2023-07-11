clear;
clc;

%{
fluxComparisons[i].xlsx layout

1. fluxVals | 2. fluxVals Interp. | 3. fluxVals Norm | 4. fluxVals Detrend 
5. diff Interp. | 6. diff Norm | 7. diff Detrend

%}

% cd("Flux Comparisons") % Enter directory with all data files

diffInterp = zeros(99, 10, 5);
diffNorm = zeros(99, 10, 5);
diffDetrend = zeros(99, 10, 5);


% Load % differences for each flux calculation
for i = 1:20
    rawFileName = append('Flux Comparisons/fluxComparisons', int2str(i), '.xlsx'); % Sets the file name to be imported
    diffInterp(:, :, i) = readmatrix(rawFileName, 'Sheet', 5); % Pulls data from 4th sheet
    % fprintf("Imported fluxComparisons " + i + " Interpolated Differences \n");

    diffNorm(:, :, i) = readmatrix(rawFileName, 'Sheet', 6); % Pulls data from 4th sheet
    % fprintf("Imported fluxComparisons " + i + " Normalized Differences \n");

    diffDetrend(:, :, i) = readmatrix(rawFileName, 'Sheet', 7); % Pulls data from 4th sheet
    % fprintf("Imported fluxComparisons " + i + " Detrended Differences \n");
end

%{

Column layout
1. Max | 2. 20% | 3. 30% | 4. 40% | 5. 50%
6. 60% | 7. 70% | 8. 80% | 9. 90% | 10. 100%

%}

cycles = (1:99).';

% Interpolation Flux Comparisons
interpAvgs = mean(diffInterp, 3); % calculate mean of each test at each cycle and flux percentage
interpMins = min(diffInterp, [], 3); % calculate min of each test at each cycle and flux percentage
interpMaxs = max(diffInterp, [], 3); % calculate max of each test at each cycle and flux percentage

% Normalized Flux Comparisons
normAvgs = mean(diffNorm, 3); % calculate mean of each test at each cycle and flux percentage
normMins = min(diffNorm, [], 3); % calculate min of each test at each cycle and flux percentage
normMaxs = max(diffNorm, [], 3); % calculate max of each test at each cycle and flux percentage

% Detrended Flux Comparisons
detrendAvgs = mean(diffDetrend, 3); % calculate mean of each test at each cycle and flux percentage
detrendMins = min(diffDetrend, [], 3); % calculate min of each test at each cycle and flux percentage
detrendMaxs = max(diffDetrend, [], 3); % calculate max of each test at each cycle and flux percentage

% TODO: Abstract and put into for loop for 20 40, 60, 80, 100
% probably 2:10:2

figure()
subplot(2, 2, 1)

col = 10;

plot(cycles, interpAvgs(:, col), 'b', 'LineWidth', 2.0); hold on;
fill([cycles; flipud(cycles)], [interpMins(:, col); flipud(interpMaxs(:, col))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel("Cycle")
ylabel("Percent Difference")
title("Interpolated Flux Difference")

subplot(2, 2, 2)
plot(cycles, normAvgs(:, col), 'b', 'LineWidth', 2.0); hold on;
fill([cycles; flipud(cycles)], [normMins(:, col); flipud(normMaxs(:, col))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel("Cycle")
ylabel("Percent Difference")
title("Normalized Flux Difference")

subplot(2, 2, 3)
plot(cycles, detrendAvgs(:, col), 'b', 'LineWidth', 2.0); hold on;
fill([cycles; flipud(cycles)], [detrendMins(:, col); flipud(detrendMaxs(:, 10))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel("Cycle")
ylabel("Percent Difference")
title("Detrended Flux Difference")



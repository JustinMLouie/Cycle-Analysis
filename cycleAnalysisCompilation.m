clear;
clc;

%{
fluxComparisons[i].xlsx layout

1. fluxVals | 2. fluxVals Interp. | 3. fluxVals Norm 
4. fluxVals Detrend | 5. fluxVals All | 6. diff Interp 
7. diff Norm | 8. diff Norm | 9. diff All

%}

% cd("Flux Comparisons") % Enter directory with all data files

diffInterp = zeros(99, 10, 20);
diffNorm = zeros(99, 10, 20);
diffDetrend = zeros(99, 10, 20);
diffAll = zeros(99, 10, 20);


% Load % differences for each flux calculation
for i = 1:20
    if i ~= 3 % Dataset 3 omitted due to missing data
        rawFileName = append('Flux Comparisons/fluxComparisons', int2str(i), '.xlsx'); % Sets the file name to be imported
        diffInterp(:, :, i) = readmatrix(rawFileName, 'Sheet', 10); % Pulls data from 6th sheet
        diffNorm(:, :, i) = readmatrix(rawFileName, 'Sheet', 11); % Pulls data from 7th sheet
        diffDetrend(:, :, i) = readmatrix(rawFileName, 'Sheet', 12); % Pulls data from 8th sheet
        diffAll(:, :, i) = readmatrix(rawFileName, 'Sheet', 13); % Pulls data from 9th sheet
    end
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

% All Flux Comparisons
allAvgs = mean(diffAll, 3); % calculate mean of each test at each cycle and flux percentage
allMins = min(diffAll, [], 3); % calculate min of each test at each cycle and flux percentage
allMaxs = max(diffAll, [], 3); % calculate max of each test at each cycle and flux percentage

for i = 2:2:10 % loop through to create graphs of the 20%, 40%, 60%, 80%, and 100% graphs
    
    percent = i * 10;

    figure()
    sgtitle(percent + "% Flux Avg Comparisons"); % Create title of overall Figure

    % Interpolated Diff Subplot
    subplot(2, 2, 1)
    plot(cycles, interpAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
    fill([cycles; flipud(cycles)], [interpMins(:, i); flipud(interpMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    ylim([-30 30])
    xlabel("Cycle");
    ylabel("Percent Difference");
    title("Interpolated Flux Difference");

    % Normalized Diff Subplot
    subplot(2, 2, 2)
    plot(cycles, normAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
    fill([cycles; flipud(cycles)], [normMins(:, i); flipud(normMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    ylim([-30 30])
    xlabel("Cycle");
    ylabel("Percent Difference");
    title("Normalized Flux Difference");

    % Detrended Diff Subplot
    subplot(2, 2, 3)
    plot(cycles, detrendAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
    fill([cycles; flipud(cycles)], [detrendMins(:, i); flipud(detrendMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    ylim([-30 30])
    xlabel("Cycle")
    ylabel("Percent Difference")
    title("Detrended Flux Difference")

    % All Diff Plot
    subplot(2, 2, 4)
    plot(cycles, allAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
    fill([cycles; flipud(cycles)], [allMins(:, i); flipud(allMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    ylim([-30 30])
    xlabel("Cycle")
    ylabel("Percent Difference")
    title("All Data Processing Flux Difference")

    fileName = append("Figures/fluxComparisons", int2str(percent), ".png"); % Creates file name for figure
    saveas(gcf,fileName); % Saves current figure 
end


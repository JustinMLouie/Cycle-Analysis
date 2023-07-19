
%% Setup

clear;
clc;

fluxVals = zeros(99, 10, 19);

diffInterp = zeros(99, 10, 19);
diffNorm = zeros(99, 10, 19);
diffDetrend = zeros(99, 10, 19);
diffAll = zeros(99, 10, 19);

diffNoFit = zeros(99, 10, 13);
diffInterpNoFit = zeros(99, 10, 13);
diffDetrendNoFit = zeros(99, 10, 13);
diffAllNoFit = zeros(99, 10, 13);

% Load % differences for each flux calculation
for i = 1:20
    rawFileName = append('Flux Comparisons No Fit/fluxComparisons', int2str(i), '.xlsx'); % Sets the file name to be imported

    % Biexponential comparisons: dataset 3 omitted due to missing data
    if i < 3 
        fluxVals(:, :, i) = readmatrix(rawFileName, 'Sheet', 1);
        diffInterp(:, :, i) = readmatrix(rawFileName, 'Sheet', 10); % Pulls data from 10th sheet
        diffNorm(:, :, i) = readmatrix(rawFileName, 'Sheet', 11); % Pulls data from 11th sheet
        diffDetrend(:, :, i) = readmatrix(rawFileName, 'Sheet', 12); % Pulls data from 12th sheet
        diffAll(:, :, i) = readmatrix(rawFileName, 'Sheet', 13); % Pulls data from 13th sheet

        % diffNoFit(:, :, i - 5) = readmatrix(rawFileName, 'Sheet', 14); % Pulls data from 14th sheet
        % diffInterpNoFit(:, :, i - 5) = readmatrix(rawFileName, 'Sheet', 15); % Pulls data from 15th sheet
        % diffDetrendNoFit(:, :, i - 5) = readmatrix(rawFileName, 'Sheet', 16); % Pulls data from 16th sheet
        % diffAllNoFit(:, :, i - 5) = readmatrix(rawFileName, 'Sheet', 17); % Pulls data from 17th sheet

    elseif i > 3
        fluxVals(:, :, i - 1) = readmatrix(rawFileName, 'Sheet', 1);
        diffInterp(:, :, i - 1) = readmatrix(rawFileName, 'Sheet', 10); % Pulls data from 10th sheet
        diffNorm(:, :, i - 1) = readmatrix(rawFileName, 'Sheet', 11); % Pulls data from 11th sheet
        diffDetrend(:, :, i - 1) = readmatrix(rawFileName, 'Sheet', 12); % Pulls data from 12th sheet
        diffAll(:, :, i - 1) = readmatrix(rawFileName, 'Sheet', 13); % Pulls data from 13th sheet
    end

    % Biexponential comparisons: 
    % Datasets 1-5 omitted due to missing data
    % Datasets 6, 10 omitted due to outliers

    if i > 6 
        diffNoFit(:, :, i - 6) = readmatrix(rawFileName, 'Sheet', 14); % Pulls data from 14th sheet
        diffInterpNoFit(:, :, i - 6) = readmatrix(rawFileName, 'Sheet', 15); % Pulls data from 15th sheet
        diffDetrendNoFit(:, :, i - 6) = readmatrix(rawFileName, 'Sheet', 16); % Pulls data from 16th sheet
        diffAllNoFit(:, :, i - 6) = readmatrix(rawFileName, 'Sheet', 17); % Pulls data from 17th sheet
    end

    if i > 10 % exclude fluxComparisons10
        diffNoFit(:, :, i - 7) = readmatrix(rawFileName, 'Sheet', 14); % Pulls data from 14th sheet
        diffInterpNoFit(:, :, i - 7) = readmatrix(rawFileName, 'Sheet', 15); % Pulls data from 15th sheet
        diffDetrendNoFit(:, :, i - 7) = readmatrix(rawFileName, 'Sheet', 16); % Pulls data from 16th sheet
        diffAllNoFit(:, :, i - 7) = readmatrix(rawFileName, 'Sheet', 17); % Pulls data from 17th sheet
    end
    
end

%{

Column layout
1. Max | 2. 20% | 3. 30% | 4. 40% | 5. 50%
6. 60% | 7. 70% | 8. 80% | 9. 90% | 10. 100%

%}

cycles = (1:99).';

%% Calculations

% Control Flux Vals
[controlAvgs, controlMins, controlMaxs] = calculateStats(fluxVals);

% Fitted Flux Comparisons 
[interpAvgs, interpMins, interpMaxs] = calculateStats(diffInterp);
[normAvgs, normMins, normMaxs] = calculateStats(diffNorm);
[detrendAvgs, detrendMins, detrendMaxs] = calculateStats(diffDetrend);
[allAvgs, allMins, allMaxs] = calculateStats(diffAll);

% Non-fitted Comparisons
[noFitAvgs, noFitMins, noFitMaxs] = calculateStats(diffNoFit);
[interpNoFitAvgs, interpNoFitMins, interpNoFitMaxs] = calculateStats(diffInterpNoFit);
[detrendNoFitAvgs, detrendNoFitMins, detrendNoFitMaxs] = calculateStats(diffDetrendNoFit);
[allNoFitAvgs, allNoFitMins, allNoFitMaxs] = calculateStats(diffAllNoFit);


%% Control Graphs
figure();
sgtitle("Control Flux Values"); % Create title of overall Figure

subplot(2, 3, 1)
col = 1;
percent = col * 10;
plot(cycles, controlAvgs(:, col), 'b', 'LineWidth', 2.0); hold on; 
fill([cycles; flipud(cycles)], [controlMins(:, col); flipud(controlMaxs(:, col))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
ylim([0 50]);
xlabel("Cycle");
ylabel("Flux (gCO2 / (m^2 * h))");
title("Max Flux")

for i = 1:5
    subplot(2, 3, i + 1)
    col = i * 2;
    percent = col * 10;
    plot(cycles, controlAvgs(:, col), 'b', 'LineWidth', 2.0); hold on; 
    fill([cycles; flipud(cycles)], [controlMins(:, col); flipud(controlMaxs(:, col))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
    ylim([0 50]);
    xlabel("Cycle");
    ylabel("Flux (gCO2 / (m^2 * h))");
    title(percent + "% Flux Avg")
end



%% Comparison Graphs

for i = 2:2:10 % loop through to create graphs of the 20%, 40%, 60%, 80%, and 100% graphs

    percent = i * 10;

    %%% Biexponentially Fitted Comparisons
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

    figure()
    sgtitle(percent + "% Flux Avg Comparisons: Non-Fitted Pressure Vals"); % Create title of overall Figure

    %%% Non Fitted Comparisons

    % NonFitted Diff Subplot
    subplot(2, 2, 1)
    plot(cycles, noFitAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
    fill([cycles; flipud(cycles)], [noFitMins(:, i); flipud(noFitMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    ylim([-100 100])
    xlabel("Cycle");
    ylabel("Percent Difference");
    title("Raw Transducer Flux Difference");

    % Nonfitted Interpolated Diff Subplot
    subplot(2, 2, 2)
    plot(cycles, interpNoFitAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
    fill([cycles; flipud(cycles)], [interpNoFitMins(:, i); flipud(interpNoFitMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    ylim([-100 100])
    xlabel("Cycle");
    ylabel("Percent Difference");
    title("Non-Fitted Interpolated Flux Difference");

    % Nonfitted Detrended Diff Subplot
    subplot(2, 2, 3)
    plot(cycles, detrendNoFitAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
    fill([cycles; flipud(cycles)], [detrendNoFitMins(:, i); flipud(detrendNoFitMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    ylim([-100 100])
    xlabel("Cycle")
    ylabel("Percent Difference")
    title("Non-Fitted Detrended Flux Difference")

    % Nonfitted All Diff Plot
    subplot(2, 2, 4)
    plot(cycles, allNoFitAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
    fill([cycles; flipud(cycles)], [allNoFitMins(:, i); flipud(allNoFitMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    ylim([-100 100])
    xlabel("Cycle")
    ylabel("Percent Difference")
    title("Non-Fitted All Data Processing Flux Difference")

    fileName = append("Figures/fluxComparisons", int2str(percent), "NoFit", ".png"); % Creates file name for figure
    saveas(gcf,fileName); % Saves current figure 

end

fprintf("Script completed");

%% Functions

function [avgs, mins, maxs] = calculateStats(diff)
    avgs = mean(diff, 3); % calculate mean of each test at each cycle and flux percentage
    mins = min(diff, [], 3); % calculate min of each test at each cycle and flux percentage
    maxs = max(diff, [], 3); % calculate max of each test at each cycle and flux percentage
end

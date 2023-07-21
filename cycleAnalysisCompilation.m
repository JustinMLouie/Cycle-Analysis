
%% Setup

clear;
clc;

fluxVals = zeros(99, 10, 19);
fluxValsDetrend = zeros(99, 10, 19);

diffInterp = zeros(99, 10, 19);
diffNorm = zeros(99, 10, 19);
diffDetrend = zeros(99, 10, 19);
diffAll = zeros(99, 10, 19);

diffNoFit = zeros(99, 10, 13);
diffInterpNoFit = zeros(99, 10, 13);
diffDetrendNoFit = zeros(99, 10, 13);
diffAllNoFit = zeros(99, 10, 13);

dPValsControl = zeros(99, 19);
dPValsNorm = zeros(99, 19);
dPValsDetrend = zeros(99, 19);
dPValsAll = zeros(99, 19);

diffDPNorm = zeros(99, 19);
diffDPDetrend = zeros(99, 19);
diffDPAll = zeros(99, 19);

ccVals = zeros(99, 2, 20);

% Load % differences for each flux calculation
for i = 1:20
    rawFileName = append('Flux Comparisons dP/fluxComparisons', int2str(i), '.xlsx'); % Sets the file name to be imported
    
    fprintf(rawFileName)

    ccVals(:, :, i) = readmatrix(rawFileName, 'Sheet', 25);

    % Biexponential comparisons: dataset 3 omitted due to missing data
    if i < 3 
        fluxVals(:, :, i) = readmatrix(rawFileName, 'Sheet', 1);
        fluxValsDetrend(:, :, i) = readmatrix(rawFileName, 'Sheet', 4);
        
        diffInterp(:, :, i) = readmatrix(rawFileName, 'Sheet', 10); % Pulls data from 10th sheet
        diffNorm(:, :, i) = readmatrix(rawFileName, 'Sheet', 11); % Pulls data from 11th sheet
        diffDetrend(:, :, i) = readmatrix(rawFileName, 'Sheet', 12); % Pulls data from 12th sheet
        diffAll(:, :, i) = readmatrix(rawFileName, 'Sheet', 13); % Pulls data from 13th sheet

        dPValsControl(:, i) = readmatrix(rawFileName, 'Sheet', 18);
        dPValsDetrend(:, i) = readmatrix(rawFileName, 'Sheet', 19);
        dPValsNorm(:, i) = readmatrix(rawFileName, 'Sheet', 20);
        dPValsAll(:, i) = readmatrix(rawFileName, 'Sheet', 21);

        diffDPNorm(:, i) = readmatrix(rawFileName, 'Sheet', 22);
        diffDPDetrend(:, i) = readmatrix(rawFileName, 'Sheet', 23);
        diffDPAll(:, i) = readmatrix(rawFileName, 'Sheet', 24);

    elseif i > 3
        fluxVals(:, :, i - 1) = readmatrix(rawFileName, 'Sheet', 1);
        fluxValsDetrend(:, :, i - 1) = readmatrix(rawFileName, 'Sheet', 4);

        diffInterp(:, :, i - 1) = readmatrix(rawFileName, 'Sheet', 10); % Pulls data from 10th sheet
        diffNorm(:, :, i - 1) = readmatrix(rawFileName, 'Sheet', 11); % Pulls data from 11th sheet
        diffDetrend(:, :, i - 1) = readmatrix(rawFileName, 'Sheet', 12); % Pulls data from 12th sheet
        diffAll(:, :, i - 1) = readmatrix(rawFileName, 'Sheet', 13); % Pulls data from 13th sheet

        dPValsControl(:, i  - 1) = readmatrix(rawFileName, 'Sheet', 18);
        dPValsDetrend(:, i  - 1) = readmatrix(rawFileName, 'Sheet', 19);
        dPValsNorm(:, i - 1) = readmatrix(rawFileName, 'Sheet', 20);
        dPValsAll(:, i - 1) = readmatrix(rawFileName, 'Sheet', 21);

        diffDPNorm(:, i - 1) = readmatrix(rawFileName, 'Sheet', 22);
        diffDPDetrend(:, i - 1) = readmatrix(rawFileName, 'Sheet', 23);
        diffDPAll(:, i - 1) = readmatrix(rawFileName, 'Sheet', 24);
    end

    % Non-Biexponential comparisons: 
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

% Flux Vals
[controlAvgs, controlMins, controlMaxs] = calculateStats(fluxVals, 3);
[detrendValAvgs, detrendValMins, detrendValMaxs] = calculateStats(fluxValsDetrend, 3);

% Fitted Flux Comparisons 
[interpAvgs, interpMins, interpMaxs] = calculateStats(diffInterp, 3);
[normAvgs, normMins, normMaxs] = calculateStats(diffNorm, 3);
[detrendAvgs, detrendMins, detrendMaxs] = calculateStats(diffDetrend, 3);
[allAvgs, allMins, allMaxs] = calculateStats(diffAll, 3);

% Non-fitted Comparisons
[noFitAvgs, noFitMins, noFitMaxs] = calculateStats(diffNoFit, 3);
[interpNoFitAvgs, interpNoFitMins, interpNoFitMaxs] = calculateStats(diffInterpNoFit, 3);
[detrendNoFitAvgs, detrendNoFitMins, detrendNoFitMaxs] = calculateStats(diffDetrendNoFit, 3);
[allNoFitAvgs, allNoFitMins, allNoFitMaxs] = calculateStats(diffAllNoFit, 3);

% dP Comparisons
[dPControlAvgs, dPControlMins, dPControlMaxs] = calculateStats(dPValsControl, 2);
[dPNormAvgs, dPNormMins, dPNormMaxs] = calculateStats(dPValsNorm, 2);
[dPDetrendAvgs, dPDetrendMins, dPDetrendMaxs] = calculateStats(dPValsDetrend, 2);
[dPAllAvgs, dPAllMins, dPAllMaxs] = calculateStats(dPValsAll, 2);

[diffNormAvgs, diffNormMins, diffNormMaxs] = calculateStats(diffDPNorm, 2);
[diffDetrendAvgs, diffDetrendMins, diffDetrendMaxs] = calculateStats(diffDPDetrend, 2);
[diffAllAvgs, diffAllMins, diffAllMaxs] = calculateStats(diffDPAll, 2);

% cc Comparisons
[ccAvgs, ccMins, ccMaxs] = calculateStats(ccVals, 3);

%% Flux Graphs
% figure();
% plot(cycles, controlAvgs(:, 8), 'b', 'LineWidth', 2.0); hold on; 
% fill([cycles; flipud(cycles)], [controlMins(:, 8); flipud(controlMaxs(:, 8))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
% ylim([0 50]);
% xlabel("Cycle");
% ylabel("Flux (gCO2 / (m^2 * h))");

% % Detrended Diff Subplot
% figure()
% plot(cycles, controlAvgs(:, 8), 'LineWidth', 2.0); hold on;
% plot(cycles, detrendValAvgs(:, 8), 'LineWidth', 2.0);
% xlim([1 20])
% xlabel("Cycle")
% ylabel("Flux (gCO2 / (m^2 * h))");
% legend("Control Flux", "Detrended Flux")
% hold off;

%% CC Graphs
figure();
plot(cycles, ccAvgs(:, 1), 'b', 'LineWidth', 2.0); hold on; 
fill([cycles; flipud(cycles)], [ccMins(:, 1); flipud(ccMaxs(:, 1))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
xlabel("Cycle");
ylabel("Charge Capacity (uAh)");


figure();
plot(cycles, ccAvgs(:, 2), 'b', 'LineWidth', 2.0); hold on; 
fill([cycles; flipud(cycles)], [ccMins(:, 2); flipud(ccMaxs(:, 2))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
xlabel("Cycle");
ylabel("Charge Capacity (Coulumbs)");


%% Pressure graphs

%%% Norm
figure();
subplot(2, 1, 1);
plot(cycles, dPControlAvgs,'b', 'LineWidth', 2.0); hold on;
fill([cycles; flipud(cycles)], [dPControlMins; flipud(dPControlMaxs)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
xlabel("Cycle"); 
ylabel("dP (psi)");
title("Control dP Values")

subplot(2, 1, 2)
plot(cycles, dPNormAvgs,'b', 'LineWidth', 2.0); hold on;
fill([cycles; flipud(cycles)], [dPNormMins; flipud(dPNormMaxs)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
xlabel("Cycle"); 
ylabel("dP (psi)");
title("Normalized dP Values")

figure()
plot(cycles, diffNormAvgs,'b', 'LineWidth', 2.0); hold on;
fill([cycles; flipud(cycles)], [diffNormMins; flipud(diffNormMaxs)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
xlabel("Cycle"); 
ylabel("Percent Difference")
title("Normalized dP Difference")


% %%% Detrend
% figure();
% subplot(2, 1, 1);
% plot(cycles, dPControlAvgs,'b', 'LineWidth', 2.0); hold on;
% fill([cycles; flipud(cycles)], [dPControlMins; flipud(dPControlMaxs)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
% xlabel("Cycle"); 
% ylabel("dP (psi)");
% title("Control dP Values")
% 
% subplot(2, 1, 2)
% plot(cycles, dPDetrendAvgs,'b', 'LineWidth', 2.0); hold on;
% fill([cycles; flipud(cycles)], [dPDetrendMins; flipud(dPDetrendMaxs)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
% xlabel("Cycle"); 
% ylabel("dP (psi)");
% title("Detrended dP Values")
% 
% figure()
% plot(cycles, diffDetrendAvgs,'b', 'LineWidth', 2.0); hold on;
% fill([cycles; flipud(cycles)], [diffDetrendMins; flipud(diffDetrendMaxs)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
% xlabel("Cycle"); 
% ylabel("Percent Difference")
% title("Detrended dP Difference")



%% Comparison Graphs
% 
% for i = 2:2:10 % loop through to create graphs of the 20%, 40%, 60%, 80%, and 100% graphs
% 
%     percent = i * 10;
% 
%     %%% Biexponentially Fitted Comparisons
%     figure() 
%     sgtitle(percent + "% Flux Avg Comparisons"); % Create title of overall Figure
% 
%     % Interpolated Diff Subplot 
%     subplot(2, 2, 1) 
%     plot(cycles, interpAvgs(:, i), 'b', 'LineWidth', 2.0); hold on; 
%     fill([cycles; flipud(cycles)], [interpMins(:, i); flipud(interpMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
%     ylim([-30 30])
%     xlabel("Cycle"); 
%     ylabel("Percent Difference"); 
%     title("Interpolated Flux Difference");
% 
%    % Normalized Diff Subplot
%     subplot(2, 2, 2)
%     plot(cycles, normAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
%     fill([cycles; flipud(cycles)], [normMins(:, i); flipud(normMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%     ylim([-30 30])
%     xlabel("Cycle");
%     ylabel("Percent Difference");
%     title("Normalized Flux Difference");
% 
%     % Detrended Diff Subplot
%     subplot(2, 2, 3)
%     plot(cycles, detrendAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
%     fill([cycles; flipud(cycles)], [detrendMins(:, i); flipud(detrendMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%     ylim([-30 30])
%     xlabel("Cycle")
%     ylabel("Percent Difference")
%     title("Detrended Flux Difference")
% 
%     % All Diff Plot
%     subplot(2, 2, 4)
%     plot(cycles, allAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
%     fill([cycles; flipud(cycles)], [allMins(:, i); flipud(allMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%     ylim([-30 30])
%     xlabel("Cycle")
%     ylabel("Percent Difference")
%     title("All Data Processing Flux Difference")
% 
%     fileName = append("Figures/fluxComparisons", int2str(percent), ".png"); % Creates file name for figure
%     saveas(gcf,fileName); % Saves current figure
% 
%     figure()
%     sgtitle(percent + "% Flux Avg Comparisons"); % Create title of overall Figure
% 
%     %%% Non Fitted Comparisons
% 
%     % NonFitted Diff Subplot
%     subplot(2, 2, 1)
%     plot(cycles, noFitAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
%     fill([cycles; flipud(cycles)], [noFitMins(:, i); flipud(noFitMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%     ylim([-30 30])
%     xlabel("Cycle");
%     ylabel("Percent Difference");
%     title("Raw Transducer Flux Difference");
% 
%     % Nonfitted Interpolated Diff Subplot
%     subplot(2, 2, 2)
%     plot(cycles, interpNoFitAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
%     fill([cycles; flipud(cycles)], [interpNoFitMins(:, i); flipud(interpNoFitMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%     ylim([-30 30])
%     xlabel("Cycle");
%     ylabel("Percent Difference");
%     title("Non-Fitted Interpolated Flux Difference");
% 
%     % Nonfitted Detrended Diff Subplot
%     subplot(2, 2, 3)
%     plot(cycles, detrendNoFitAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
%     fill([cycles; flipud(cycles)], [detrendNoFitMins(:, i); flipud(detrendNoFitMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%     ylim([-30 30])
%     xlabel("Cycle")
%     ylabel("Percent Difference")
%     title("Non-Fitted Detrended Flux Difference")
% 
%     % Nonfitted All Diff Plot
%     subplot(2, 2, 4)
%     plot(cycles, allNoFitAvgs(:, i), 'b', 'LineWidth', 2.0); hold on;
%     fill([cycles; flipud(cycles)], [allNoFitMins(:, i); flipud(allNoFitMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%     ylim([-30 30])
%     xlabel("Cycle")
%     ylabel("Percent Difference")
%     title("Non-Fitted All Data Processing Flux Difference")
% 
%     fileName = append("Figures/fluxComparisons", int2str(percent), "NoFit", ".png"); % Creates file name for figure
%     saveas(gcf,fileName); % Saves current figure 
% 
% end

fprintf("Script completed");

%% Functions

function [avgs, mins, maxs] = calculateStats(diff, dim)
    avgs = mean(diff, dim); % calculate mean of each test at each cycle 
    mins = min(diff, [], dim); % calculate min of each test at each cycle
    maxs = max(diff, [], dim); % calculate max of each test at each cycle
end
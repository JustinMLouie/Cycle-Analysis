
%% Setup

clear;
clc;

%%% High Performance Arrays 

numValidFiles = 14;
numCycles = 29;

fluxValsControl = zeros(numCycles, 10, numValidFiles);
fluxValsInterp = zeros(numCycles, 10, numValidFiles);
fluxValsNoFit = zeros(numCycles, 10, numValidFiles);
fluxValsInterpNoFit = zeros(numCycles, 10, numValidFiles);

% Read data from fluxComparison files

files = findDatasheets("Flux Comparisons High Performance", "fluxComparisons");
numFiles = length(files);

for i = 1:numFiles
    rawFileName = string(files(i));
    
    fprintf(rawFileName + "\n");

    % Pull values from each file's sheet into a 3D matrix
    fluxValsControl(:, :, i) = readmatrix(rawFileName, 'Sheet', 1);
    fluxValsInterp(:, :, i) = readmatrix(rawFileName, 'Sheet', 2);
    fluxValsNoFit(:, :, i) = readmatrix(rawFileName, 'Sheet', 6);
    fluxValsInterpNoFit(:, :, i) = readmatrix(rawFileName, 'Sheet', 7);
end 

%{

Column layout
1. Max | 2. 20% | 3. 30% | 4. 40% | 5. 50%
6. 60% | 7. 70% | 8. 80% | 9. 90% | 10. 100%

%}

%% Tabulate 

% Creates array of cycles
cycles = (1:numCycles).'; % High Performance Data

tabulateFluxes(fluxValsControl, numFiles, "Control Flux Values");
tabulateFluxes(fluxValsInterp, numFiles, "Interpolated Flux Values");
tabulateFluxes(fluxValsNoFit, numFiles, "No Fit Flux Values");
tabulateFluxes(fluxValsInterpNoFit, numFiles, "No Fit Interpolated Flux Values");


[controlAvgs, controlMins, controlMaxs] = calculateStats(fluxValsControl, 3);
[interpAvgs, interpMins, interpMaxs] = calculateStats(fluxValsControl, 3);
[noFitAvgs, noFitmins, noFitMaxs] = calculateStats(fluxValsControl, 3);
[interpNoFitAvgs, interpNoFitMins, interpNoFitMaxs] = calculateStats(fluxValsControl, 3);

fprintf("Script completed");

%% Comparison Graphs

for i = 2:2:10 % loop through to create graphs of the 20%, 40%, 60%, 80%, and 100% graphs

    percent = i * 10;

    %%% Flux Comparisons
    figure() 
    sgtitle(percent + "% Flux Avg Comparisons"); % Create title of overall Figure

    % Control: Fitted
    subplot(2, 2, 1) 
    plot(cycles, controlAvgs(:, i), 'b', 'LineWidth', 2.0); hold on; 
    fill([cycles; flipud(cycles)], [controlMins(:, i); flipud(controlMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
    ylim([0 100])
    xlabel("Cycle"); 
    ylabel("Flux (gCO2 / (m^2 * h))"); 
    title("Fitted (Control)");

    % Fitted + Interp
    subplot(2, 2, 2) 
    plot(cycles, interpAvgs(:, i), 'b', 'LineWidth', 2.0); hold on; 
    fill([cycles; flipud(cycles)], [interpMins(:, i); flipud(interpMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
    ylim([0 100])
    xlabel("Cycle"); 
    ylabel("Flux (gCO2 / (m^2 * h))"); 
    title("Fitted + Interpolated");

    % Not Fitted
    subplot(2, 2, 3) 
    plot(cycles, noFitAvgs(:, i), 'b', 'LineWidth', 2.0); hold on; 
    fill([cycles; flipud(cycles)], [noFitmins(:, i); flipud(noFitMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
    ylim([0 100])
    xlabel("Cycle"); 
    ylabel("Flux (gCO2 / (m^2 * h))"); 
    title("Non-Fitted");

    % Not Fitted + Interp
    subplot(2, 2, 4) 
    plot(cycles, interpNoFitAvgs(:, i), 'b', 'LineWidth', 2.0); hold on; 
    fill([cycles; flipud(cycles)], [interpNoFitMins(:, i); flipud(interpNoFitMaxs(:, i))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
    ylim([0 100])
    xlabel("Cycle"); 
    ylabel("Flux (gCO2 / (m^2 * h))"); 
    title("Non-Fitted + Interpolated");

    fileName = append("Figures/fluxComparisons", int2str(percent), ".png"); % Creates file name for figure
    saveas(gcf,fileName); % Saves current figure
 
end



%% Functions

function [avgs, mins, maxs] = calculateStats(diff, dim)
    avgs = mean(diff, dim); % calculate mean of each test at each cycle 
    mins = min(diff, [], dim); % calculate min of each test at each cycle
    maxs = max(diff, [], dim); % calculate max of each test at each cycle
end

function tabulateFluxes(fluxVals, numFiles, fileName)
    %{
    Create and save a file with the following format
    File (Data Processing Method) 
        --> Sheet (Avg % Flux) 
            --> Table (Cycle 17, 18 vals)
                       fluxComparison1 vals
                       fluxComparison2 vals
                       fluxComparison3 vals
                            ...
    %}
    excelFile = append('Flux Tables/', fileName, '.xlsx');

    percentage = ["Max", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "Full"];

    for i = 1:10 
        % Go through each percentage
        isolatedFluxVals = zeros(numFiles, 2);
        for j = 1:numFiles
            isolatedFluxVals(j, 1) = fluxVals(17, i, j); % Copies cycle 17
            isolatedFluxVals(j, 2) = fluxVals(18, i, j); % Copies Cycle 18
        end
        writematrix(isolatedFluxVals, excelFile, 'Sheet', percentage(i), 'Range', 'A1')
    end

    
end

function files = findDatasheets(folder, target_string)
    %{
    Searches through input folder to find files with target string
    Returns cell array of files
    %}

    % Initialize cell array
    files = {};

    % List all directories in parent folder
    subfolders = dir(folder);

    % Iterate over subfolders
    for i = 1:length(subfolders)
        if ~strcmp(subfolders(i).name, '.') && ~strcmp(subfolders(i).name, '..') && ~strcmp(subfolders(i).name, '100Co2')
            % if directory, enter directory to look for files (recursive)
            if subfolders(i).isdir()
                files = [files; findDatasheets(fullfile(folder, subfolders(i).name), target_string)];

            % If file, get file
            else
                if contains(subfolders(i).name, target_string)
                    files = [files; fullfile(folder, subfolders(i).name)];
                end
            end
        end
    end
end
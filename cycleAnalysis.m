clear;
clc;

% import matlab.io*

%% Loading Cell Volume File

% rawDataFiles = findDatasheets('.', "TemporaryTestName_Channel"); % used for poor performance cells
rawDataFiles = findDatasheets('.', "TemporaryTestName_Channel"); 
procDataFiles = findDatasheets('.', "Processed Data");

numFiles = length(rawDataFiles);

for n = 1:numFiles

    %% Set Up
    % Loading Raw Data File
    rawFileName = string(rawDataFiles(n));
    rawData = readmatrix(rawFileName, 'Sheet', 2); % Pulls data from 2nd sheet of rawData file

    % Loading Proc Data File 
    procFileName = string(procDataFiles(n));
    cellVolData = readmatrix(procFileName, 'Sheet', 1); % Pulls data from 2nd sheet of rawData file
  
    % Constants used in calculations
    cellVolume = cellVolData(8, 2) * 1e-3; % Processed Data cell volume, L
    uAhToCoulumbs = 0.0036;
    cellArea = 2e-4; % Elecrode area, m^2
    mwCO2 = 44; % Molecular weight of CO2, g/mol 
    F = 96485; % Faraday Constant, (s * A)/mol

    % Key Metrics calculated from Raw Data
    numCycles = max(rawData(:, 5));
    cycleNum = rawData(:, 5);
    mcVals = zeros(numCycles - 2, 9); % 20% 30% 40% 50% 60% 70 80 90 Full
    fluxVals = zeros(numCycles - 2, 10); % Max 20% 30% 40% 50% 60% 70 80 90 Full
    ccVals = zeros(numCycles - 2, 2); % uAh, Coulumbs
    feVals = zeros(numCycles - 2, 1);

    % Flux Values to Compare Against
    % Max 20% 30% 40% 50% 60% 70 80 90 Full
    fluxValsInterp = zeros(numCycles - 2, 10); % Fitted --> Interpolated Flux Values
    fluxValsNorm = zeros(numCycles - 2, 10); % Normalized --> Fitted --> Flux Values
    fluxValsDetrend = zeros(numCycles - 2, 10); % Detrended --> Fitted --> Flux Values
    fluxValsAll = zeros(numCycles - 2, 10); % Detrended --> Normalized --> Fitted --> Interpolated Values

    fluxValsNoFit = zeros(numCycles - 2, 10); % Unfitted --> Calc Flux Vals
    fluxValsInterpNoFit = zeros(numCycles - 2, 10); % Unfitted --> Interp Flux Vals
    fluxValsDetrendNoFit = zeros(numCycles - 2, 10); % Detrend --> Unfitted --> Calc Flux Vals
    fluxValsAllNoFit = zeros(numCycles - 2, 10); % Detrend --> Unfitted --> Interp Flux Vals

    % DP Vals
    dPValsControl = zeros(numCycles - 2, 1);
    dPValsDetrend = zeros(numCycles - 2, 1);
    dPValsNorm = zeros(numCycles - 2, 1);
    dPValsAll = zeros(numCycles - 2, 1);
    
    % Caclulate detrendData
    detrendPressure = detrend(rawData(:, 19), 10); % Detrend the pressure values
    detrendData = rawData; % Create a copy of rawData 
    detrendData(:, 19) = detrendPressure; % copy the detrended pressurevalues to detrendData

    % %% Detrend diagram
    % figure()
    % subplot(2, 1, 1)
    % plot(rawData(:, 3), rawData(:, 19), "LineWidth", 2.0)
    % ylabel("Pressure (psi)", "LineWidth", 2.0)
    % xlabel("Time (s)")
    % title("Transducer Pressures")
    % 
    % subplot(2, 1, 2)
    % plot(detrendData(:, 3), detrendData(:, 19), "LineWidth", 2.0)
    % ylabel("Pressure (psi)")
    % xlabel("Time (s)")
    % title("Detrended Pressures")

    %% Calculating Values

    fprintf("Starting calculations \n")

    % Loops through each cycle, skips the incomplete cycle 1
    % figure()
    for i = 2:numCycles - 1 % iterate through cycles

        % Extracts cycle data
        cycData = rawData(rawData(:, 5) == i, :); % Data in cycle i
        chargeRegionData = cycData(cycData(:, 7) > 0.00002, :);
        stepIndexVals = chargeRegionData(5, 6);
        chargeRegionData = cycData(cycData(:, 7) > 0, :); % Data in charging region, theoretical > 0, realistic > 0.0002
        chargeRegionData = chargeRegionData(chargeRegionData(:, 6) == stepIndexVals, :);

        timeVals = chargeRegionData(:, 4); % Step time values, s
        pressureVals = chargeRegionData(:, 19); % Auxilary Presure Values, psi
        chargeVals = chargeRegionData(:, 10) * 1e6; % Charge Capacity, uAh

        % Extract Detrend Cycle Data
        cycDataDetrend = detrendData(detrendData(:, 5) == i, :); % Data in cycle i
        chargeDataDetrend = cycDataDetrend(cycDataDetrend(:, 7) > 0.00002, :);
        stepIndexValsDetrend = chargeDataDetrend(5, 6);
        chargeDataDetrend = cycDataDetrend(cycDataDetrend(:, 7) > 0, :); % Data in charging region, theoretical > 0, realistic > 0.0002
        chargeDataDetrend = chargeDataDetrend(chargeDataDetrend(:, 6) == stepIndexValsDetrend, :);

        % Calculates time-pressure values of biexponential fit for pressure curve
        pFitTimePressure = createTimePressureCurves(timeVals, pressureVals);

        % Calculate Normalized Pressure Data
        [pFitTimePressureNorm, rescaleVal] = normalizeTimePressureVals(timeVals, pressureVals);
        dPValsNorm(i - 1) = rescaleVal;
        timeValsNorm = pFitTimePressureNorm(:, 1);
        pressureValsNorm = pFitTimePressureNorm(:, 2);
        
        % Calculate Detrended Data
        timeValsDetrend = chargeDataDetrend(:, 4);
        pressureValsDetrend = chargeDataDetrend(:, 19);
        pFitTimePressureDetrend = createTimePressureCurves(timeValsDetrend, pressureValsDetrend); % creates biexponential curve of detrended data
        dPValsDetrend(i - 1) = pFitTimePressureDetrend(1, 2) - pFitTimePressureDetrend(end, 2);

        % Calculate All Data Processing (Detrend --> Normalized --> Fitted --> Interpolated)
        [pFitTimePressureAll, rescaleValAll] = normalizeTimePressureVals(timeValsDetrend, pressureValsDetrend); % Normalize + Fit the detrended data
        dPValsAll(i - 1) = rescaleValAll;
        timeValsAll = pFitTimePressureAll(:, 1);
        pressureValsAll = pFitTimePressureAll(:, 2);

        % Create non-fitted time-pressure value arrays
        pressureTime = createTimePressureVals(timeVals, pressureVals);
        pressureTimeDetrend = createTimePressureVals(timeValsDetrend, pressureValsDetrend);

        % Calculates Max Flux Values
        fluxVals(i - 1, 1) = calculateFluxMax(timeVals, pressureVals, cellArea, cellVolume); % Biexponential
        fluxValsInterp(i - 1, 1) = calculateFluxMax(timeVals, pressureVals, cellArea, cellVolume); % Normalized
        fluxValsNorm(i - 1, 1) = calculateFluxMax(timeValsNorm, pressureValsNorm, cellArea, cellVolume) * rescaleVal; % Normalized
        fluxValsDetrend(i - 1, 1) = calculateFluxMax(timeVals, pressureValsDetrend, cellArea, cellVolume); % Detrended
        fluxValsAll(i - 1, 1) = calculateFluxMax(timeValsAll, pressureValsAll, cellArea, cellVolume) * rescaleValAll; % All data processing applied
        
        %%%%% Did not calculate max flux for non-fitted, cannot calc
        %%%%% without curve, no use doing so

        % Calculates minimum pressure drops for avg flux calcs
        minP = calculateMinP(pFitTimePressure);
        minPNorm = calculateMinP(pFitTimePressureNorm);
        minPDetrend = calculateMinP(pFitTimePressureDetrend);
        minPAll = calculateMinP(pFitTimePressureAll); 

        minPNoFit = calculateMinP(pressureTime);
        minPDetrendNoFit = calculateMinP(pressureTimeDetrend);

        % Calculates CC Vals
        ccuAh = max(chargeVals); % charge capacity, uAh
        ccCoulumbs = ccuAh * uAhToCoulumbs; % convert cc to coulumbs

        % Calculates FE Vals
        dPFit = pFitTimePressure(1, 2) - pFitTimePressure(end, 2);
        fe = 100 * (psiToGCO2(dPFit, cellVolume) / mwCO2) * F / (ccCoulumbs);

        % Assign values at end of loop
        ccVals(i - 1, 1) = ccuAh;
        ccVals(i - 1, 2) = ccCoulumbs;
        feVals(i - 1) = fe;
        dPValsControl(i - 1) = dPFit;

        % %%% Plot to demonstrate what a biexponential Fit is
        % figure()
        % scatter(timeVals(2:2:end), pressureVals(2:2:end), "LineWidth", 2.0); hold on;
        % plot(pFitTimePressure(:, 1), pFitTimePressure(:, 2), "LineWidth", 2.0);
        % xlabel("Time (s)")
        % ylabel("Pressure (psi)")
        % legend("Transducer Data", "Biexponentially Fitted Data")
        % 
        % %%% Plot to demonstrate non-fitted data
        % figure()
        % scatter(timeVals(2:4:end), pressureVals(2:4:end), "LineWidth", 2.0); hold on;
        % plot(timeVals, pressureVals, "LineWidth", 2.0);
        % xlabel("Time (s)")
        % ylabel("Pressure (psi)")
        % legend("Transducer Data", "Biexponentially Fitted Data")

        % subplot(3, 10, i)
        % plot(timeVals, pressureVals); hold on;
        % plot(pFitTimePressure(:, 1), pFitTimePressure(:, 2));
        % xlabel("Time");
        % ylabel("Pressure")
        % title("Cycle " + i)
        % 
        % if i == 15
        %     figure()
        %     plot(timeVals, pressureVals); hold on;
        %     plot(pFitTimePressure(:, 1), pFitTimePressure(:, 2));
        %     xlabel("Time");
        %     ylabel("Pressure")
        %     title("Cycle " + i)
        % end

        % % plots to demonstrate normalization
        % if i < 10
        %     figure()
        %     subplot(2, 1, 1)
        %     plot(timeVals, pressureVals, "LineWidth", 2.0)
        % 
        %     yline(pressureVals(1), '--k', "LineWidth", 1.0);
        %     yline(pressureVals(end), '--k', "LineWidth", 1.0);
        %     ylabel("Pressure (psi)")
        %     xlabel("Time (s)")
        %     subplot(2, 1, 2)
        %     plot(pFitTimePressureNorm(:, 1), pFitTimePressureNorm(:, 2), "LineWidth", 2.0)
        %     ylabel("Pressure (psi)")
        %     xlabel("Time (s)")
        % end
        
        % Plot to demonstrate 20% Cutoff calculation errors + interpolation
        % if i < 20
        %     figure()
        %     subplot(2, 1, 1);
        %     x = (1:100);
        %     y1 = (x - 50) .* (x - 50);
        %     y2 =  x * -50 + 2500;
        %     plot(x, y1, "LineWidth", 2.0); hold on;
        %     plot(x, y2, "LineWidth", 2.0);
        %     ylim([0 2500])
        %     xlim([0 60])
        %     xline(25, '--k', "LineWidth", 1.5);
        %     yline(1250, '--k', "LineWidth", 1.5);
        % 
        %     subplot(2, 1, 2);
        %     plot(pFitTimePressure(:, 1), pFitTimePressure(:, 2), 'LineWidth', 2.0); hold on;
        %     scatter(pFitTimePressure(:, 1), pFitTimePressure(:, 2), 'LineWidth', 1.5);
        %     xlim([0 10])
        %     interpTime = interp1(pFitTimePressure(:, 2), pFitTimePressure(:, 1), minP(1)); % End Time, Interpolated at min pressure
        %     xline(interpTime, '--k', 'LineWidth', 1.5);
        %     yline(minP(1), '--k', 'LineWidth', 1.5);
        %     xlabel("Time (s)")
        %     ylabel("Pressure (psi)")
        % end

        % if i == 20
        %     figure()
        %     plot(pressureTime(:, 1), pressureTime(:, 2), 'b', 'LineWidth', 1.0)
        %     hold on;
        %     yline(minPNoFit(7), '--b', 'LineWidth', 1.0);
        %     plot(pFitTimePressure(:, 1), pFitTimePressure(:, 2), 'r', 'LineWidth', 1.0);
        %     yline(minP(7), '--r', 'LineWidth', 1.0);
        %     ylabel("Pressure (psi)")
        %     xlabel("Time (s)")
        %     legend("Unfitted", "80% Unfitted Cutoff", "Fitted", "80% Fitted Cutoff")
        % 
        % end

        % if i == 20
        %     figure()
        %     yyaxis left
        %     scatter(timeVals, pressureVals, 'k', "LineWidth", 1.0)
        %     hold on;
        %     plot(pFitTimePressure(:, 1), pFitTimePressure(:, 2), '-b', 'LineWidth', 1.0); % Fitted curve
        %     ylabel("Pressure (psi)", 'Color', 'c')
        %     yline(minP(7), '--b', 'LineWidth', 1.0)
        %     yline(minP(6), '--b', 'LineWidth', 1.0)
        % 
        %     yyaxis right
        %     plot(timeValsNorm, pressureValsNorm, '-r', 'LineWidth', 1.0); % Normalized Curve
        %     yline(minPNorm(7), '--r', 'LineWidth', 1.0)
        %     yline(minPNorm(6), '--r', 'LineWidth', 1.0)
        % 
        %     % xlim([0 100])
        % 
        %     ax = gca; % Get current axes handle
        %     ax.YAxis(1).Color = 'b'; % Left y-axis color
        %     ax.YAxis(2).Color = 'r'; % Right y-axis color
        % 
        %     ylabel("Pressure (unitless)", 'Color', 'r')
        %     xlabel("Time (s)");
        %     legend("Transducer", "Control", "80% Control", "70% Control", "Normalized", "80% Normalized", "70% Normalized")
        %     title("rawdata" + n + " Cycle " + i);
        % 
        % end

        % Avg calculations at 20-100%
        minPData = [];
        for j = 1:9 % iterate through percentages

            % Calculates Mass Capture Values - not relevant to flux calcs
            minPData =  pFitTimePressure(pFitTimePressure(:, 2) >= minP(j), :);
            pStart = minPData(1, 1); % Start Time
            pEnd = minPData(end, 1); % End Time
            dP = pStart - pEnd; % Pressure drop
            mc = psiToGCO2(dP, cellVolume) / cellArea;

            % Flux Calculations
            flux = calculateFluxAvgs(pFitTimePressure, minP(j), cellVolume); % Biexponetial
            fluxInterp = interpolateFluxAvgs(pFitTimePressure, minP(j), cellVolume); % Biexponential --> Interpolated
            fluxNorm = calculateFluxAvgs(pFitTimePressureNorm, minPNorm(j), cellVolume); % Normalized --> Biexponential
            fluxDetrend = calculateFluxAvgs(pFitTimePressureDetrend, minPDetrend(j), cellVolume); % Detrend --> Biexponential
            fluxAll = interpolateFluxAvgs(pFitTimePressureAll, minPAll(j), cellVolume); % Detrend --> Normalized --> Biexponential --> Interpolated

            fluxNoFit = calculateFluxAvgs(pressureTime, minPNoFit(j), cellVolume); % No fit --> calc flux
            fluxInterpNoFit = interpolateFluxAvgs(pressureTime, minPNoFit(j), cellVolume); % No fit --> interp flux
            fluxDetrendNoFit = calculateFluxAvgs(pressureTimeDetrend, minPDetrendNoFit(j), cellVolume); % Detrend --> no fit --> calc flux
            fluxAllNoFit = interpolateFluxAvgs(pressureTimeDetrend, minPDetrendNoFit(j), cellVolume); % Detrend --> no fit --> interp flux

            % Rescale Normalized Flux vals
            fluxNorm = fluxNorm * rescaleVal; % Rescale flux vals
            fluxAll = fluxAll * rescaleValAll; % Rescale flux vals

            % Records final calculations
            mcVals(i - 1, j + 1) = mc;
            fluxVals(i - 1, j + 1) = flux;
            fluxValsInterp(i - 1, j + 1) = fluxInterp;
            fluxValsNorm(i - 1, j + 1) = fluxNorm;
            fluxValsDetrend(i - 1, j + 1) = fluxDetrend;
            fluxValsAll(i - 1, j + 1) = fluxAll;

            fluxValsNoFit(i - 1, j + 1) = fluxNoFit;
            fluxValsInterpNoFit(i - 1, j + 1) = fluxInterpNoFit;
            fluxValsDetrendNoFit(i - 1, j + 1) = fluxDetrendNoFit;
            fluxValsAllNoFit(i - 1, j + 1) = fluxAllNoFit;

        end
    end

    % Post Calculation Evaluations

    cycColumn = (1:size(fluxVals)).';
    columnTitles = {'Cycle', 'Max', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', 'Full'};

    % Round flux vals to 2 decimals
    fluxVals = round(fluxVals, 2);
    fluxValsInterp = round(fluxValsInterp, 2);
    fluxValsNorm = round(fluxValsNorm, 2);
    fluxValsDetrend = round(fluxValsDetrend, 2);
    fluxValsAll = round(fluxValsAll, 2);

    fluxValsNoFit = round(fluxValsNoFit, 2);
    fluxValsInterpNoFit = round(fluxValsInterpNoFit, 2);
    fluxValsDetrendNoFit = round(fluxValsDetrendNoFit, 2);
    fluxValsAllNoFit = round(fluxValsAllNoFit, 2);

    % Percent Diff Calculations
    diffInterp = round(100 * (fluxValsInterp - fluxVals) ./ fluxVals, 2);
    diffNorm = round(100 * (fluxValsNorm - fluxVals) ./ fluxVals, 2);
    diffDetrend = round(100 * (fluxValsDetrend - fluxVals) ./ fluxVals, 2);
    diffAll = round(100 * (fluxValsAll - fluxVals) ./ fluxVals, 2);

    diffNoFit = round(100 * (fluxValsNoFit - fluxVals) ./ fluxVals, 2);
    diffInterpNoFit = round(100 * (fluxValsInterpNoFit - fluxVals) ./ fluxVals, 2);
    diffDetrendNoFit = round(100 * (fluxValsDetrendNoFit - fluxVals) ./ fluxVals, 2);
    diffAllNoFit = round(100 * (fluxValsAllNoFit - fluxVals) ./ fluxVals, 2);
    
    diffDPNorm = round(100 * (dPValsNorm - dPValsControl) ./ dPValsControl, 2);
    diffDPDetrend = round(100 * (dPValsDetrend - dPValsControl) ./ dPValsControl, 2);
    diffDPAll = round(100 * (dPValsAll - dPValsControl) ./ dPValsControl, 2);

    fprintf("rawdata" + n + " ended" + "\n");

    %% Write data to Excel file
    if n < 10
        num = append("0", int2str(n));
    else
        num = int2str(n);
    end

    excelFile = append('fluxComparisons', num, '.xlsx');

    % Biexpontially fitted fluxVals 
    writematrix(fluxVals, excelFile, 'Sheet', 'Flux Vals (Control)', 'Range', 'A1'); % Sheet 1
    writematrix(fluxValsInterp, excelFile, 'Sheet', 'Flux Vals Interp', 'Range', 'A1'); % Sheet 2
    writematrix(fluxValsNorm, excelFile, 'Sheet', 'Flux Vals Norm', 'Range', 'A1'); % Sheet 3
    writematrix(fluxValsDetrend, excelFile, 'Sheet', 'Flux Vals Detrend', 'Range', 'A1'); % Sheet 4
    writematrix(fluxValsAll, excelFile, 'Sheet', 'Flux Vals All', 'Range', 'A1'); % Sheet 5

    % Non fitted fluxVals 
    writematrix(fluxValsNoFit, excelFile, 'Sheet', 'Flux Vals No Fit', 'Range', 'A1'); % Sheet 6
    writematrix(fluxValsInterpNoFit, excelFile, 'Sheet', 'Flux Vals Interp No Fit', 'Range', 'A1'); % Sheet 7
    writematrix(fluxValsDetrendNoFit, excelFile, 'Sheet', 'Flux Vals Detrend No Fit', 'Range', 'A1'); % Sheet 8
    writematrix(fluxValsAllNoFit, excelFile, 'Sheet', 'Flux Vals All No Fit', 'Range', 'A1'); % Sheet 9

    % Biexpontially fitted fluxVals diff 
    writematrix(diffInterp, excelFile, 'Sheet', 'Diff Interp', 'Range', 'A1'); % Sheet 10
    writematrix(diffNorm, excelFile, 'Sheet', 'Diff Norm', 'Range', 'A1'); % Sheet 11
    writematrix(diffDetrend, excelFile, 'Sheet', 'Diff Detrend', 'Range', 'A1'); % Sheet 12
    writematrix(diffAll, excelFile, 'Sheet', 'diff All', 'Range', 'A1'); % Sheet 13

    % Non fitted fluxVals diff
    writematrix(diffNoFit, excelFile, 'Sheet', 'Diff No Fit', 'Range', 'A1'); % Sheet 14
    writematrix(diffInterpNoFit, excelFile, 'Sheet', 'Diff Interp No Fit', 'Range', 'A1'); % Sheet 15
    writematrix(diffDetrendNoFit, excelFile, 'Sheet', 'Diff Detrend No Fit', 'Range', 'A1'); % Sheet 16
    writematrix(diffAllNoFit, excelFile, 'Sheet', 'Diff All No Fit', 'Range', 'A1'); % Sheet 17

    % Biexpontially fitted dP 
    writematrix(dPValsControl, excelFile, 'Sheet', 'dP Vals Control', 'Range', 'A1'); % Sheet 18
    writematrix(dPValsDetrend, excelFile, 'Sheet', 'dP Vals Detrend', 'Range', 'A1'); % Sheet 19
    writematrix(dPValsNorm, excelFile, 'Sheet', 'dP Vals Norm', 'Range', 'A1'); % Sheet 20
    writematrix(dPValsAll, excelFile, 'Sheet', 'dP Vals All', 'Range', 'A1'); % Sheet 21

    % Biexponentially fitted dP
    writematrix(diffDPNorm, excelFile, 'Sheet', 'Diff dP Norm', 'Range', 'A1'); % Sheet 22
    writematrix(diffDPDetrend, excelFile, 'Sheet', 'Diff dP Detrend', 'Range', 'A1'); % Sheet 23
    writematrix(diffDPAll, excelFile, 'Sheet', 'Diff dP All', 'Range', 'A1'); % Sheet 24

    writematrix(ccVals, excelFile, 'Sheet', 'cc vals', 'Range', 'A1') % Sheet 25

end

%% Extraneous Functions

function timePressureVals = createTimePressureCurves(timeVals, pressureVals)
    %{
    Calculate the biexpnential fitted curve, and stores values in a
    timePressure array
    %}
    % Create biexponential fit for pressure curve
    pFitCurve = fit(timeVals, pressureVals, 'exp2'); % Creates biexponential model of P vs t
    pFitEval = pFitCurve(timeVals); % Calculates pressures at given times on pFitCurve

    % Stores the values of the fitted exponential curve
    timePressureVals = createTimePressureVals(timeVals, pFitEval);
end

function timePressureVals = createTimePressureVals(timeVals, pressureVals)
    % Creates array of time and pressure
    timePressureVals = [];
    timePressureVals(:, 1) = timeVals; % set column 1 to time values
    timePressureVals(:, 2) = pressureVals; % sets column 2 to pressures values
end

function [normalizedTPVals, normalizingDP] = normalizeTimePressureVals(timeVals, pressureVals) 
    %{
    Calculate Normalized Pressure Values
    First Normalization: Normalized = (pressure - min) / dP
    %}
    normalizingValue = min(pressureVals); % Calculates baseline pressure
    normalizingDP = pressureVals(1) - mean(pressureVals(end - 4:end)); % Calculates range
    pressureValsNorm = (pressureVals - normalizingValue) / normalizingDP; % Normalizes pressures between 0 and 1
    timeValsNorm = timeVals - timeVals(1); % Sets time to start at 0
    normalizedTPVals = createTimePressureCurves(timeValsNorm, pressureValsNorm);
end

function minP = calculateMinP(timePressureVals)
    %{
    Calculates minimum pressure drops required for the following pressure calculations
    20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, Full
    %}
    percentages = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
    pStart = timePressureVals(1, 2); % Start Pressure
    pEnd = timePressureVals(end, 2); % End Pressure
    dPFit = pStart - pEnd; % Pressure drop
    minP = pStart - (dPFit * percentages); % Minimum pressure values for each avg flux calc
end

function gCO2 = psiToGCO2(psi, cellVolume)
    %{
    Converts psi to grams CO2
    Calculates using PV = NRT, then multiplies by molecular weight 
    %}
    psiToPa = 6894.76; % Pa psi-1

    R = 8.314 * 1e3; % Gas const ,  (L * Pa)/(mol * K)
    T = 298; %  erature, assume RT = 298K
    mwCO2 = 44; % Molecular weight of CO2, g/mol

    pascals = psi * psiToPa;
    molsCO2 = (pascals * cellVolume) / (R * T);
    gCO2 = molsCO2 * mwCO2;
end

function flux = pressureToFlux(dp, dt, cellVolume)
    %{
    Calculates flux values using change in pressure and time
    Flux = dP / (area * dt)
    Desired units: g CO2/ h
    %}
    secToHr = 1/3600;
    cellArea = 2e-4; % Elecrode area, m^2
    dtHours = dt * secToHr;
    dPGramsCO2 = psiToGCO2(dp, cellVolume);
    flux = dPGramsCO2 / (cellArea * dtHours); % calculate cycle flux, gCo2 / h
end

function fluxCalc = calculateFluxAvgs(timePressureVals, minP, cellVolume)
    %{
    Calculates flux by determining the pressure greater than minP
    Use time-pressure pairs to calculate flux
    %}
    % Constrains dataset to have minimum threshold pressure, subtract
    % 0.00001 to prevent rounding issues with flux
    minPData =  timePressureVals(timePressureVals(:, 2) >= (minP - 0.00001), :); 
    % minPData = timePressureVals;
    pStart = minPData(1, 2); % Start pressure
    pEnd = minPData(end, 2); % End Pressure
    tStart = minPData(1, 1); % Start Time
    tEnd = minPData(end, 1); % End Time
    dP = pStart - pEnd; % Pressure drop
    dt = tEnd - tStart; % Time elapsed
    fluxCalc = pressureToFlux(dP, dt, cellVolume); % Calculates flux
end

function fluxInterp = interpolateFluxAvgs(timePressureVals, minP, cellVolume)
    %{
    Calculates flux using interpolated points
    Utilize exact pressures 20 - 100% to interpolate time
    Utilize interpolated time-pressure pairs to calculate avg fluxes
    %}

    pEndInterp = 0;
    if (minP < timePressureVals(end, 2))
        pEndInterp = timePressureVals(end, 2);
    else
        pEndInterp = minP; % End Pressure
    end
    
    pStartInterp = timePressureVals(1, 2); % Start pressure
    tStartInterp = timePressureVals(1, 1); % Start time

    % Ensure P vs t is unique
    [~, uniqueIndices] = unique(timePressureVals(:, 2), 'stable');
    finalPressureVals = timePressureVals(uniqueIndices, :);

    tEndInterp = interp1(finalPressureVals(:, 2), finalPressureVals(:, 1), pEndInterp); % End Time, Interpolated at min pressure
    dPInterp = pStartInterp - pEndInterp; % Pressure drop
    dtInterp = tEndInterp - tStartInterp; % Time elapsed
    fluxInterp = pressureToFlux(dPInterp, dtInterp, cellVolume); % Calculates flux
end

function maxFlux = calculateFluxMax(timeVals, pressureVals, cellArea, cellVolume)
    %{
    Calculates Max flux
    Takes derivative of biexponential curve
    Converts flux to gCo2 / (m^2 * hr)
    Takes maximum derivative value
    %}
    secToHr = 1/3600;
    pFitCurve = fit(timeVals, pressureVals, 'exp2', 'StartPoint', [0,0,0,0]); % Creates biexponential model of P vs t
    pDeriv = differentiate(pFitCurve, timeVals); % Differentiates pFitCurve with datapoints at x = timeVals
    gDeriv = -1 * psiToGCO2(pDeriv, cellVolume) / cellArea / secToHr; % Convert pDeriv from psi to g CO2
    maxFlux = max(gDeriv);
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


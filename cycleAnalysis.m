clear;
clc;

% import matlab.io*

%% Loading Cell Volume File

cd("Raw Data") % Enter directory with all data files
volData = readmatrix("cellvolumes.xlsx", "Sheet", 1);

for n = 1:20
    % Loading Raw Data File
    rawFileName = append('rawdata', int2str(n), '.xlsx'); % Sets the file name to be imported
    rawData = readmatrix(rawFileName, 'Sheet', 2); % Pulls data from 2nd sheet of rawData file

    %% Set Up
    % Constants used in calculations
    cellVolume = volData(n, 2) * 1e-3; % Processed Data cell volume, L
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
    fluxValsInterp = zeros(numCycles - 2, 10); % Interpolated Flux Values
    fluxValsNorm = zeros(numCycles - 2, 10); % Normalized Flux Values
    fluxValsDetrend = zeros(numCycles - 2, 10); % Detrended Flux Values
    fluxValsAll = zeros(numCycles - 2, 10); % Detrended --> Normalized --> Fitted --> Interpolated Values

    % Caclulate detrendData
    detrendPressure = detrend(rawData(:, 19), 10); % Detrend the pressure values
    detrendData = rawData; % Create a copy of rawData 
    detrendData(:, 19) = detrendPressure; % copy the detrended pressurevalues to detrendData

    % Plots p vs time for entire data sets
    % figure()
    % subplot(2,1,1);
    % plot(rawData(:, 3), rawData(:, 19)); hold on;
    % xlabel("Test Time (s)")
    % ylabel("Pressure (gCO2)")
    % title("Transducer Pressure Values")
    % 
    % subplot(2,1,2);
    % plot(detrendData(:, 3), detrendData(:, 19)); hold on;
    % xlabel("Test Time (s)")
    % ylabel("Pressure (gCO2)")
    % title("Detrended Pressure Values")

    %% Calculating Values

    % fprintf("Starting calculations \n")

    % Loops through each cycle, skips the incomplete cycle 1
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
        pFitTimePressure = setTimePressureCurves(timeVals, pressureVals);

        % Calculate Normalized Pressure Data
        pFitTimePressureNorm = normalizeTimePressureVals(timeVals, pressureVals);
        rescaleVals = mean(pressureVals(1:5)) - mean(pressureVals(end - 5:end)); 
        timeValsNorm = pFitTimePressureNorm(:, 1);
        pressureValsNorm = pFitTimePressureNorm(:, 2);
        
        % Calculate Detrended Data
        timeValsDetrend = chargeDataDetrend(:, 4);
        pressureValsDetrend = chargeDataDetrend(:, 19);
        pFitTimePressureDetrend = setTimePressureCurves(timeValsDetrend, pressureValsDetrend); % creates biexponential curve of detrended data

        % Calculate All Data Processing (Detrend --> Normalized --> Fitted --> Interpolated)
        pFitTimePressureAll = normalizeTimePressureVals(timeValsDetrend, pressureValsDetrend); % Normalize + Fit the detrended data
        rescaleValsAll = mean(pressureValsDetrend(1:5)) - mean(pressureValsDetrend(end - 5:end)); % Create value to rescale detrended data
        timeValsAll = pFitTimePressureAll(:, 1);
        pressureValsAll = pFitTimePressureAll(:, 2);

        % Calculates Max Flux Values
        fluxVals(i - 1, 1) = calculateFluxMax(timeVals, pressureVals, cellArea, cellVolume); % Biexponential
        fluxValsInterp(i - 1, 1) = calculateFluxMax(timeVals, pressureVals, cellArea, cellVolume); % Normalized
        fluxValsNorm(i - 1, 1) = calculateFluxMax(timeValsNorm, pressureValsNorm, cellArea, cellVolume) * rescaleVals; % Normalized
        fluxValsDetrend(i - 1, 1) = calculateFluxMax(timeVals, pressureValsDetrend, cellArea, cellVolume); % Detrended
        fluxValsAll(i - 1, 1) = calculateFluxMax(timeValsAll, pressureValsAll, cellArea, cellVolume) * rescaleValsAll; % All data processing applied

        % Calculates minimum pressure drops for avg flux calcs
        minP = calculateMinP(pFitTimePressure);
        minPNorm = [0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0];
        % minPNorm = calculateMinP(pFitTimePressureNorm);
        minPDetrend = calculateMinP(pFitTimePressureDetrend);
        minPAll = [0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, pressureValsAll(end)]; 
        % minPAll = calculateMinP(pFitTimePressureAll);

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

        % Avg calculations at 20-100%
        minPData = [];
        for j = 1:9 % iterate through percentages

            % Graph Cycle 80
            % if (i == 80) 
            %     if (j == 1)
            %     figure();
            %     plot(pFitTimePressure(:, 1), pFitTimePressure(:, 2),  'b', 'LineWidth', 1.0); hold on;
            %     scatter(pFitTimePressure(:, 1), pFitTimePressure(:, 2),  'b', 'LineWidth', 2.0);
            %     xlim([0 10]);
            %     yline(minP(j), 'k', 'LineWidth', 2.0);
            %     xlabel("Step Time")
            %     ylabel("Pressure (gCO2) ")
            %     title("rawdata3 Cycle 80")
            %     legend("", "Transducer Pressure Measurements", "20% Flux");
            %     end
            % end

            % Calculates Mass Capture Values - not relevant to flux calcs
            minPData =  pFitTimePressure(pFitTimePressure(:, 2) > minP(j), :);
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

            % Rescale Normalized Flux vals
            fluxNorm = fluxNorm * rescaleVals; % Rescale flux vals
            fluxAll = fluxAll * rescaleValsAll; % Rescale flux vals

            % Records final calculations
            mcVals(i - 1, j + 1) = mc;
            fluxVals(i - 1, j + 1) = flux;
            fluxValsInterp(i - 1, j + 1) = fluxInterp;
            fluxValsNorm(i - 1, j + 1) = fluxNorm;
            fluxValsDetrend(i - 1, j + 1) = fluxDetrend;
            fluxValsAll(i - 1, j + 1) = fluxAll;

        end
    end

    %% Post Calculation Evaluations

    cycColumn = (1:size(fluxVals)).';
    columnTitles = {'Cycle', 'Max', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', 'Full'};

    % Round flux vals to 2 decimals
    fluxVals = round(fluxVals, 2);
    fluxValsInterp = round(fluxValsInterp, 2);
    fluxValsNorm = round(fluxValsNorm, 2);
    fluxValsDetrend = round(fluxValsDetrend, 2);
    fluxValsAll = round(fluxValsAll, 2);

    % Percent Diff Calculations
    diffInterp = round(100 * (fluxValsInterp - fluxVals) ./ fluxVals, 2);
    diffNorm = round(100 * (fluxValsNorm - fluxVals) ./ fluxVals, 2);
    diffDetrend = round(100 * (fluxValsDetrend - fluxVals) ./ fluxVals, 2);
    diffAll = round(100 * (fluxValsAll - fluxVals) ./ fluxVals, 2);

    diffInterpMax = max(abs(diffInterp));
    diffNormMax = max(abs(diffNorm));
    diffDetrendMax = max(abs(diffDetrend));
    diffAllMax = max(abs(diffAll));

    figure()
    plotName = append("rawdata", int2str(n));
    sgtitle(plotName); 
    subplot(2,1,1);
    plot(cycColumn, fluxVals(:, 10), 'b', 'LineWidth', 1.0); hold on;
    ylim([0 3.5])
    xlabel("Cycle")
    ylabel("Flux (gCO2 / (m^2 * h))")
    title("Flux Values")

    subplot(2,1,2)
    plot(cycColumn, fluxValsDetrend(:, 10), 'b', 'LineWidth', 1.0); hold on;
    ylim([0 3.5])
    xlabel("Cycle")
    ylabel("Flux (gCO2 / (m^2 * h))")
    title("Detrended Flux Values")

    %% Formatting Matrices
    % If columns and headers desired, uncomment section, change writematrix to writecell
    
    % % Add Cycle Column
    % fluxVals = [cycColumn, fluxVals];
    % fluxValsInterp = [cycColumn, fluxValsInterp];
    % fluxValsNorm = [cycColumn, fluxValsNorm];
    % fluxValsDetrend = [cycColumn, fluxValsDetrend];
    % fluxVals = [cycColumn, fluxValsAll];
    % diffInterp = [cycColumn, diffInterp];
    % diffNorm = [cycColumn, diffNorm];
    % diffDetrend = [cycColumn, diffDetrend];
    % diffAll = [cycColumn, diffAll];
    %
    % % Add column headers
    % fluxVals = [columnTitles; num2cell(fluxVals)];
    % fluxValsInterp = [columnTitles; num2cell(fluxValsInterp)];
    % fluxValsNorm = [columnTitles; num2cell(fluxValsNorm)];
    % fluxValsDetrend = [columnTitles; num2cell(fluxValsDetrend)];
    % fluxValsAll = [columnTitles; num2cell(fluxValsAll)];
    % diffInterp = [columnTitles; num2cell(diffInterp)];
    % diffNorm = [columnTitles; num2cell(diffNorm)];
    % diffDetrend = [columnTitles; num2cell(diffDetrend)];
    % diffAll = [columnTitles; num2cell(diffAll)];

    %% Write data to Excel file
    excelFile = append('fluxComparisons', int2str(n), '.xlsx');

    % Write each data array to a separate sheet
    % writematrix(fluxVals, excelFile, 'Sheet', 'Flux Vals', 'Range', 'A1');
    % writematrix(fluxValsInterp, excelFile, 'Sheet', 'Flux Vals Interp', 'Range', 'A1');
    % writematrix(fluxValsNorm, excelFile, 'Sheet', 'Flux Vals Norm', 'Range', 'A1');
    % writematrix(fluxValsDetrend, excelFile, 'Sheet', 'Flux Vals Detrend', 'Range', 'A1');
    % writematrix(diffDetrend, excelFile, 'Sheet', 'Flux Vals All', 'Range', 'A1');
    % writematrix(diffInterp, excelFile, 'Sheet', 'Diff Interp', 'Range', 'A1');
    % writematrix(diffNorm, excelFile, 'Sheet', 'Diff Norm', 'Range', 'A1');
    % writematrix(diffDetrend, excelFile, 'Sheet', 'Diff Detrend', 'Range', 'A1');
    % writematrix(diffDetrend, excelFile, 'Sheet', 'Diff All', 'Range', 'A1');
    % 
    fprintf("rawdata" + n + " ended" + "\n");
end

%% Extraneous Functions

function timePressureVals = setTimePressureCurves(timeVals, pressureVals)
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

function normalizedTPVals = normalizeTimePressureVals(timeVals, pressureVals) 
    %{
    Calculate Normalized Pressure Values
    First Normalization: Normalized = (pressure - min) / dP
    Second Normalization: 
        a' = a/(a + c)
        c' = c/(a + c) 
        P = a' * exp(b * t) + c' * exp(d*t)
    %}
    normalizingValue = min(pressureVals); % Calculates baseline pressure
    normalizingDP = range(pressureVals); % Calculates range
    pressureValsNorm = (pressureVals - normalizingValue) / normalizingDP; % Normalizes pressures between 0 and 1
    timeValsNorm = timeVals - timeVals(1); % Sets time to start at 0
    pFitCurveNorm = fit(timeValsNorm, pressureValsNorm, 'exp2'); % First Normalization
    % pFitEvalNorm = pFitCurveNorm(timeValsNorm); % TEMP: Only used to evaluate single normalization
    % normalizedTPVals = createTimePressureVals(timeValsNorm, pFitEvalNorm); % TEMP: Only used to evaluate single normalization

    curveCoeff = coeffvalues(pFitCurveNorm); % extract coeffs from biexponential fit
    aPrime = curveCoeff(1) / (curveCoeff(1) + curveCoeff(3)); % Normalizing Fitted value coeffs
    cPrime = curveCoeff(3) / (curveCoeff(1) + curveCoeff(3)); % Normalizing Fitted value coeffs
    eqtn = fittype('a * exp(b * x) + c * exp(d*x)'); % Defining biexponential eqtn
    pFitCurveNorm2 = cfit(eqtn, aPrime, curveCoeff(2), cPrime, curveCoeff(4));
    pFitCurveNormVals = pFitCurveNorm2(timeValsNorm); % Calculating Second Normalization p vals
    normalizedTPVals = createTimePressureVals(timeValsNorm, pFitCurveNormVals);
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
    T = 298; % Temperature, assume RT = 298K
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
    pEndInterp = minP; % End Pressure
    pStartInterp = timePressureVals(1, 2); % Start pressure
    tStartInterp = timePressureVals(1, 1); % Start time
    [finalPressureVals, ~, ~] =...
        unique(timePressureVals, 'rows', 'stable'); % Ensures P vs t is unique
    tEndInterp = interp1(finalPressureVals(:, 2), ...
        finalPressureVals(:, 1), pEndInterp); % End Time, Intepolated at min pressure
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
    pFitCurve = fit(timeVals, pressureVals, 'exp2'); % Creates biexponential model of P vs t
    pDeriv = differentiate(pFitCurve, timeVals); % Differentiates pFitCurve with datapoints at x = timeVals
    gDeriv = -1 * psiToGCO2(pDeriv, cellVolume) / cellArea / secToHr; % Convert pDeriv from psi to g CO2
    maxFlux = max(gDeriv);
end





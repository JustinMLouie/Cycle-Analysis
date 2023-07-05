clear;
clc;

%% Loading rawData file
[rawFile, rawPath] = uigetfile('.xlsx');
rawFileName = string(fullfile(rawPath, rawFile));

rawData = readmatrix(rawFileName, 'Sheet', 2); % Pulls data from 2nd sheet of rawData file


%% Setting Relevant Constants
% Constants used in calculations
cellVolume = 5.45884579 * 1e-3; % Processed Data cell volume, L
psiToPa = 6894.76; % Pa psi-1
uAhToCoulumbs = 0.0036;
cellArea = 2e-4; % Elecrode area, m^2
mwCO2 = 44.01; % Molecular weight of CO2, g/mol
R = 8.314 * 1e3; % Gas const ,  (L * Pa)/(mol * K)
T = 298; % Temperature, assume RT = 298K
F = 96485; % Faraday Constant, (s * A)/mol
secToHr = 1/3600;% 

%% Setting up 
numCycles = max(rawData(:, 5)); 
finalMCVals = zeros(numCycles - 2, 9); % 20% 30% 40% 50% 60% 70 80 90 Full
finalFluxVals = zeros(numCycles - 2, 10); % Max 20% 30% 40% 50% 60% 70 80 90 Full
finalCCVals = zeros(numCycles - 2, 2); % uAh, Coulumbs
finalFEVals = zeros(numCycles - 2, 1);

percentages = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];

%% Calculating Values

% Loops through each cycle, skips the incomplete cycle 1
for i = 2:numCycles - 1
    
    % Extracts cycle data
    cycData = rawData(rawData(:, 5) == i, :); % Data in cycle i
    chargeRegionData = cycData(cycData(:, 7) > 0, :); % Data in charging region, theoretical > 0, realistic > 0.0002
    stepIndexVals = chargeRegionData(5,6);
    chargeRegionData = chargeRegionData(chargeRegionData(:, 6) == stepIndexVals, :);

    timeVals = chargeRegionData(:, 4); % Test time values, s
    pressureVals = chargeRegionData(:, 19); % Auxilary Presure Values, psi
    chargeVals = chargeRegionData(:, 10) * 1e6; % Charge Capacity, uAh

    % Create biexponential fit for pressure curve
    pFitCurve = fit(timeVals, pressureVals, 'exp2'); % Creates biexponential model of P vs t
    pDeriv = differentiate(pFitCurve, timeVals); % Differentiates pFitCurve with datapoints at x = timeVals

    pFitEval = pFitCurve(timeVals); % Calculates pressures at given times on pFitCurve

    % Stores the values of the fitted exponential curve
    pFitTimePressure = [];
    pFitTimePressure(:, 1) = timeVals; % set column 1 to test time values
    pFitTimePressure(:, 2) = pFitEval;

    % Calculates minimum pressure drops required for the following pressure calculations
    % 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, Full
    minPressures = zeros(1, 9);
    pStart = pFitEval(1);
    pEnd = pFitEval(end);
    dPFit = pStart - pEnd; % Pressure drop

    minPressures = pStart - dPFit * percentages; % Minimum pressure values for each avg flux calc

    % Avg calculations at 20-100%
    minPData = [];
    for j = 1:9
        minPData =  pFitTimePressure(pFitTimePressure(:, 2) > minPressures(j), :); % Constrains dataset to have minimum threshold pressure
        % fprintf("Min Data for " + percentages(j) + " with minimum pressure " + minPressures(j));
        pStartMinPData = minPData(1, 2);
        pEndMinPData = minPData(end, 2);
        tStart = minPData(1, 1);
        tEnd = minPData(end, 1);
        dP = pStartMinPData - pEndMinPData; % pressure drop
        dt = tEnd - tStart;

         % Calculates dP (grams CO2)
        dPGramsCO2 = psiToGCO2(dP); % Convert dp from psi to g

        % Calculates Flux
        dtHours = dt * secToHr; % Convert dt from s to h
        flux = dPGramsCO2 / (cellArea * dtHours); % calculate cycle flux, gCo2 / h

        % Calculates Mass Capture Values
        mc = dPGramsCO2 / cellArea;
        
        % Records final calculations
        finalMCVals(i - 1, j + 1) = mc;
        finalFluxVals(i - 1, j + 1) = flux;    

    end

    % Calculates Max Flux
    gDeriv = -1 * psiToGCO2(pDeriv) / cellArea / secToHr; % Convert pDeriv from psi to g CO2
    finalFluxVals(i - 1, 1) = max(gDeriv);

    % Calculates CC Vals
    ccuAh = max(chargeVals); % charge capacity, uAh
    ccCoulumbs = ccuAh * uAhToCoulumbs; % convert cc to coulumbs
    
    % Calculates FE Vals
    fe = 100 * (psiToGCO2(dPFit) / mwCO2) * F / (ccCoulumbs);
    
    % Assign values at end of loop
    finalCCVals(i - 1, 1) = ccuAh;
    finalCCVals(i - 1, 2) = ccCoulumbs;
    finalFEVals(i - 1) = fe;

end

fprintf("Script ended" + "\n");

% Note: Difference in values between cycleAnalysis and
% hundredCycleParameterCalculatorSingleCase stems from using
% chargeRegionData vs cycData during calculations
% ChargeRegionData is a smaller subset of cycData, and excludes some
% values


% Converts psi to grams CO2
function gCO2 = psiToGCO2(psi) 
    psiToPa = 6894.76; % Pa psi-1
    cellVolume = 5.45884579 * 1e-3; % processed data cell volume (L)
    R = 8.314 * 1e3; % Gas const ,  (L * Pa)/(mol * K)
    T = 298; % Temperature, assume RT = 298K
    mwCO2 = 44.01; % Molecular weight of CO2, g/mol

    pascals = psi * psiToPa;
    molsCO2 = (pascals * cellVolume) / (R * T);
    gCO2 = molsCO2 * mwCO2;
end

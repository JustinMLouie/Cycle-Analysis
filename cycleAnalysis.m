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
mcVals = zeros(numCycles - 2, 9); % 20% 30% 40% 50% 60% 70 80 90 Full
fluxVals = zeros(numCycles - 2, 10); % Max 20% 30% 40% 50% 60% 70 80 90 Full
ccVals = zeros(numCycles - 2, 2); % uAh, Coulumbs
feVals = zeros(numCycles - 2, 1);

% Values to Compare Against
normalizedFlux = zeros(numCycles - 2, 10); % Max 20% 30% 40% 50% 60% 70 80 90 Full
interpolatedFlux = zeros(numCycles - 2, 10); % Max 20% 30% 40% 50% 60% 70 80 90 Full
movingMeanFlux = zeros(numCycles - 2, 10); % Max 20% 30% 40% 50% 60% 70 80 90 Full

percentages = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];

%% Calculating Values

% Loops through each cycle, skips the incomplete cycle 1
for i = 2:numCycles - 1
    
    % Extracts cycle data
    cycData = rawData(rawData(:, 5) == i, :); % Data in cycle i
    chargeRegionData = cycData(cycData(:, 7) > 0, :); % Data in charging region, theoretical > 0, realistic > 0.0002
    stepIndexVals = chargeRegionData(5,6);
    chargeRegionData = chargeRegionData(chargeRegionData(:, 6) == stepIndexVals, :);

    timeVals = chargeRegionData(:, 4); % Step time values, s
    pressureVals = chargeRegionData(:, 19); % Auxilary Presure Values, psi
    chargeVals = chargeRegionData(:, 10) * 1e6; % Charge Capacity, uAh

    % Create biexponential fit for pressure curve
    pFitCurve = fit(timeVals, pressureVals, 'exp2'); % Creates biexponential model of P vs t
    pFitEval = pFitCurve(timeVals); % Calculates pressures at given times on pFitCurve

    % Stores the values of the fitted exponential curve
    pFitTimePressure = [];
    pFitTimePressure(:, 1) = timeVals; % set column 1 to time values
    pFitTimePressure(:, 2) = pFitEval;

    % Calculates minimum pressure drops required for the following pressure calculations
    % 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, Full
    minP = zeros(1, 9);
    pStart = pFitEval(1);
    pEnd = pFitEval(end);
    dPFit = pStart - pEnd; % Pressure drop

    minP = pStart - dPFit * percentages; % Minimum pressure values for each avg flux calc

    % Avg calculations at 20-100%
    minPData = [];
    for j = 1:9
        minPData =  pFitTimePressure(pFitTimePressure(:, 2) > minP(j), :); % Constrains dataset to have minimum threshold pressure
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
        mcVals(i - 1, j + 1) = mc;
        fluxVals(i - 1, j + 1) = flux;    
    end

    % Calculates Max Flux
    pDeriv = differentiate(pFitCurve, timeVals); % Differentiates pFitCurve with datapoints at x = timeVals
    gDeriv = -1 * psiToGCO2(pDeriv) / cellArea / secToHr; % Convert pDeriv from psi to g CO2
    fluxVals(i - 1, 1) = max(gDeriv);

    % Calculates CC Vals
    ccuAh = max(chargeVals); % charge capacity, uAh
    ccCoulumbs = ccuAh * uAhToCoulumbs; % convert cc to coulumbs
    
    % Calculates FE Vals
    fe = 100 * (psiToGCO2(dPFit) / mwCO2) * F / (ccCoulumbs);
    
    % Assign values at end of loop
    ccVals(i - 1, 1) = ccuAh;
    ccVals(i - 1, 2) = ccCoulumbs;
    feVals(i - 1) = fe;

    % Calculate Normalized Flux Values 
    pNormalized = normalize(pressureVals); % normalizes raw pressure values
    pNFitCurve = fit(timeVals, pNormalized, 'exp2'); % Creates biexponential model of normalized P vs t
    pNFitVals = pNFitCurve(timeVals); % Calculates normalized pressure at each time
    pNormalizedTimePressure = []; 
    pNormalizedTimePressure(:, 1) = timeVals; % set col 1 to time values
    pNormalizedTimePressure(:, 2) = pNFitVals; % set col 2 to normalized p vals
 
    % Calculates minimum pressure drops required for the following pressure calculations
    minNP = zeros(1, 9);
    pNStart = pNFitVals(1);
    pNEnd = pNFitVals(end);
    dPNFit = pNStart - pNEnd; % Pressure drop
    minNP = pNStart - dPNFit * percentages; % Pressures at x%

    normalizedFluxData = [];
    for j = 1:9
        normalizedFluxData =  pNormalizedTimePressure(pNormalizedTimePressure(:, 2) > minNP(j), :); % Constrains dataset to have minimum threshold pressure
        pStartNoramlized = normalizedFluxData(1, 2);
        pEndNormalized = normalizedFluxData(end, 2);
        tStart = normalizedFluxData(1, 1);
        tEnd = normalizedFluxData(end, 1);
        dP = pStartMinPData - pEndMinPData; % pressure drop
        dt = tEnd - tStart;

         % Calculates dP (grams CO2)
        dPGramsCO2 = psiToGCO2(dP); % Convert dp from psi to g

        % Calculates Flux
        dtHours = dt * secToHr; % Convert dt from s to h
        flux = dPGramsCO2 / (cellArea * dtHours); % calculate cycle flux, gCo2 / h

        % Records final calculations
        normalizedFlux(i - 1, j + 1) = flux; 
    end

    % Calculates Max Flux Normalized
    pNDeriv = differentiate(pNFitCurve, timeVals); % Differentiates pFitCurve with datapoints at x = timeVals
    gNDeriv = -1 * psiToGCO2(pNDeriv) / cellArea / secToHr; % Convert pDeriv from psi to g CO2
    normalizedFlux(i - 1, 1) = max(gNDeriv);

end

fprintf("Script ended" + "\n");

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


%{

    TODO:
    1. Create a normalized version of flux - calculate and compare values
    2. Create a interpolated version of flux - calculate and compare values
    3. Create a moving average version of flux - calculate and compare
    values

%}


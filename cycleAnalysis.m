clear;
clc;

%% Loading rawData file
% fprintf("Add Raw file \n")
[rawFile, rawPath] = uigetfile('.xlsx');
rawFileName = string(fullfile(rawPath, rawFile));
rawData = readmatrix(rawFileName, 'Sheet', 2); % Pulls data from 2nd sheet of rawData file

[procFile, procPath] = uigetfile('.xlsx');
procFileName = string(fullfile(procPath, procFile));
procData = readmatrix(procFileName, 'Sheet', 2); % Pulls data from 2nd sheet of rawData file

% fprintf("File loaded \n")

%% Set Up
% Constants used in calculations
cellVolume = procData(8, 2) * 1e-3; % Processed Data cell volume, L
uAhToCoulumbs = 0.0036;
cellArea = 2e-4; % Elecrode area, m^2
mwCO2 = 44; % Molecular weight of CO2, g/mol
F = 96485; % Faraday Constant, (s * A)/mol

% Key Metrics calculated from Raw Data
numCycles = max(rawData(:, 5)); 
mcVals = zeros(numCycles - 2, 9); % 20% 30% 40% 50% 60% 70 80 90 Full
fluxVals = zeros(numCycles - 2, 10); % Max 20% 30% 40% 50% 60% 70 80 90 Full
ccVals = zeros(numCycles - 2, 2); % uAh, Coulumbs
feVals = zeros(numCycles - 2, 1);

% Flux Values to Compare Against
% Max 20% 30% 40% 50% 60% 70 80 90 Full
fluxValsInterp = zeros(numCycles - 2, 10); % Interpolated Flux Values
fluxValsNorm = zeros(numCycles - 2, 10); % Interpolated Flux Values
fluxValsDetrend = zeros(numCycles - 2, 10); % Interpolated Moving Mean Trailing Values

% Caclulate detrendData
detrendPressure = detrend(rawData(:, 19), 10);
detrendData = rawData; % Create a copy of rawData
detrendData(:, 19) = detrendPressure;

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
    normalizingValue = min(pressureVals); % Calculates baseline pressure 
    normalizingDP = range(pressureVals); % Calculates range
    pValsNorm = (pressureVals - normalizingValue) / normalizingDP; % Normalizes pressures between 0 and 1
    timeValsNorm = timeVals - timeVals(1); % Sets time to start at 0
    pFitCurveNorm = fit(timeValsNorm, pValsNorm, 'exp2');
    rescaleVals = mean(pressureVals(1:5)) - mean(pressureVals(end - 5:end));

    curveCoeff = coeffvalues(pFitCurveNorm); % extract coeffs from biexponential fit
    aPrime = curveCoeff(1) / (curveCoeff(1) + curveCoeff(3)); % Normalizing Fitted value coeffs
    cPrime = curveCoeff(3) / (curveCoeff(1) + curveCoeff(3)); % Normalizing Fitted value coeffs
    eqtn = fittype('a * exp(b * x) + c * exp(d*x)'); % Defining biexponential eqtn
    pFitCurveNorm2 = cfit(eqtn, aPrime, curveCoeff(2), cPrime, curveCoeff(4));
    pFitCurveNormVals = pFitCurveNorm2(timeValsNorm); % Calculating Second Normalization p vals
    pFitTimePressureNorm = createTimePressureVals(timeValsNorm, pFitCurveNormVals);

    % Calculate Detrended Dada
    timeValsDetrend = chargeDataDetrend(:, 4);
    pressureValsDetrend = chargeDataDetrend(:, 19);
    pFitTimePressureDetrend = setTimePressureCurves(timeValsDetrend, pressureValsDetrend); % creates biexponential curve of detrended data
    
    % Calculates Max Flux Values
    fluxVals(i - 1, 1) = calculateFluxMax(timeVals, pressureVals, cellArea); % Biexponential
    fluxValsInterp(i - 1, 1) = calculateFluxMax(timeVals, pressureVals, cellArea); % Normalized
    fluxValsNorm(i - 1, 1) = calculateFluxMax(timeValsNorm, pValsNorm, cellArea) * rescaleVals; % Normalized
    fluxValsDetrend(i - 1, 1) = calculateFluxMax(timeVals, pressureValsDetrend, cellArea); % Detrended  
    
    % Calculates minimum pressure drops for avg flux calcs
    minP = calculateMinP(pFitTimePressure);
    minPNorm = [0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0];
    minPDetrend = calculateMinP(pFitTimePressureDetrend);

    % Calculates CC Vals
    ccuAh = max(chargeVals); % charge capacity, uAh
    ccCoulumbs = ccuAh * uAhToCoulumbs; % convert cc to coulumbs
    
    % Calculates FE Vals
    dPFit = pFitTimePressure(1, 2) - pFitTimePressure(end, 2);
    fe = 100 * (psiToGCO2(dPFit) / mwCO2) * F / (ccCoulumbs);
    
    % Assign values at end of loop
    ccVals(i - 1, 1) = ccuAh;
    ccVals(i - 1, 2) = ccCoulumbs;
    feVals(i - 1) = fe;

    % Avg calculations at 20-100%
    minPData = [];
    for j = 1:9 % iterate through percentages
        
        % Calculates Mass Capture Values - not relevant to flux calcs
        minPData =  pFitTimePressure(pFitTimePressure(:, 2) > minP(j), :);
        pStart = minPData(1, 1); % Start Time
        pEnd = minPData(end, 1); % End Time
        dP = pStart - pEnd; % Pressure drop
        mc = psiToGCO2(dP) / cellArea;

        % Flux Calculations
        flux = calculateFluxAvgs(pFitTimePressure, minP(j)); % Biexponetial
        fluxInterp = interpolateFluxAvgs(pFitTimePressure, minP(j)); % Biexponential --> Interpolated 
        fluxNorm = calculateFluxAvgs(pFitTimePressureNorm, minPNorm(j)); % Normalized --> Biexponential
        fluxDetrend = calculateFluxAvgs(pFitTimePressureDetrend, minPDetrend(j)); % Detrend --> Biexponential 

        fluxNorm = fluxNorm * rescaleVals; % Rescale flux vals

        % Records final calculations
        mcVals(i - 1, j + 1) = mc;
        fluxVals(i - 1, j + 1) = flux;   
        fluxValsInterp(i - 1, j + 1) = fluxInterp;
        fluxValsNorm(i - 1, j + 1) = fluxNorm;
        fluxValsDetrend(i - 1, j + 1) = fluxDetrend;
      
    end
end

% Post Calculation Evaluations

% Round flux vals to 2 decimals
fluxVals = round(fluxVals, 2);
fluxValsInterp = round(fluxValsInterp, 2);
fluxValsNorm = round(fluxValsNorm, 2);
fluxValsDetrend = round(fluxValsDetrend, 2);

% Percent Diff Calculations
% neg diff means fluxVals is greater, fluxVals overshoots (diffInterp)
diffInterp = 100 * (fluxValsInterp - fluxVals) ./ fluxVals;
diffNorm = 100 * (fluxValsNorm - fluxVals) ./ fluxVals;
diffDetrend = 100 * (fluxValsDetrend - fluxVals) ./ fluxVals;

diffInterpMax = max(abs(diffInterp));
diffNormMax = max(abs(diffNorm));
diffDetrendMax = max(abs(diffDetrend));

% fprintf("Script ended" + "\n");

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

function gCO2 = psiToGCO2(psi) 
    %{
    Converts psi to grams CO2
    Calculates using PV = NRT, then multiplies by molecular weight 
    %}
    psiToPa = 6894.76; % Pa psi-1
    cellVolume = 5.45884579 * 1e-3; % processed data cell volume (L)
    R = 8.314 * 1e3; % Gas const ,  (L * Pa)/(mol * K)
    T = 298; % Temperature, assume RT = 298K
    mwCO2 = 44; % Molecular weight of CO2, g/mol

    pascals = psi * psiToPa;
    molsCO2 = (pascals * cellVolume) / (R * T);
    gCO2 = molsCO2 * mwCO2;
end

function flux = pressureToFlux(dp, dt) 
    %{
    Calculates flux values using change in pressure and time
    Flux = dP / (area * dt)
    Desired units: g CO2/ h
    %}
    secToHr = 1/3600;
    cellArea = 2e-4; % Elecrode area, m^2
    dtHours = dt * secToHr; 
    dPGramsCO2 = psiToGCO2(dp);
    flux = dPGramsCO2 / (cellArea * dtHours); % calculate cycle flux, gCo2 / h
end

function fluxCalc = calculateFluxAvgs(timePressureVals, minP)
    %{
    Calculates flux by determining the pressure greater than minP
    Use time-pressure pairs to calculate flux
    %}
    minPData =  timePressureVals(timePressureVals(:, 2) >= minP, :); % Constrains dataset to have minimum threshold pressure
    pStart = minPData(1, 2); % Start pressure
    pEnd = minPData(end, 2); % End Pressure
    tStart = minPData(1, 1); % Start Time
    tEnd = minPData(end, 1); % End Time
    dP = pStart - pEnd; % Pressure drop
    dt = tEnd - tStart; % Time elapsed
    fluxCalc = pressureToFlux(dP, dt); % Calculates flux
end

function fluxInterp = interpolateFluxAvgs(timePressureVals, minP)
    %{
    Calculates flux using interpolated points
    Utilize exact pressures 20 - 100% to interpolate time
    Utilize interpolated time-pressure pairs to calculate avg fluxes
    %}
    pEndInterp = minP; % End Pressure
    tEndInterp = interp1(timePressureVals(:, 2), timePressureVals(:, 1), pEndInterp); % End Time, Intepolated at min pressure
    pStartInterp = timePressureVals(1, 2); % Start pressure
    tStartInterp = timePressureVals(1, 1); % Start time
    dPInterp = pStartInterp - pEndInterp; % Pressure drop
    dtInterp = tEndInterp - tStartInterp; % Time elapsed
    fluxInterp = pressureToFlux(dPInterp, dtInterp); % Calculates flux
end

function maxFlux = calculateFluxMax(timeVals, pressureVals, cellArea)
    %{
    Calculates Max flux
    Takes derivative of biexponential curve
    Converts flux to gCo2 / (m^2 * hr)
    Takes maximum derivative value
    %}
    secToHr = 1/3600;
    pFitCurve = fit(timeVals, pressureVals, 'exp2'); % Creates biexponential model of P vs t
    pDeriv = differentiate(pFitCurve, timeVals); % Differentiates pFitCurve with datapoints at x = timeVals
    gDeriv = -1 * psiToGCO2(pDeriv) / cellArea / secToHr; % Convert pDeriv from psi to g CO2
    maxFlux = max(gDeriv);
end





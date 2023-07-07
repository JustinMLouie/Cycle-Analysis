clear;
clc;

%% Loading rawData file
fprintf("Add Raw file" + "\n")
[rawFile, rawPath] = uigetfile('.xlsx');
rawFileName = string(fullfile(rawPath, rawFile));
rawData = readmatrix(rawFileName, 'Sheet', 2); % Pulls data from 2nd sheet of rawData file

fprintf("Add Proc file" + "\n")

[procFile, procPath] = uigetfile('.xlsx');
procFileName = string(fullfile(procPath, procFile));
procData = readmatrix(procFileName, 'Sheet', 2); % Pulls data from 2nd sheet of processedData file

%% Set Up
% Constants used in calculations
cellVolume = 5.45884579 * 1e-3; % Processed Data cell volume, L
psiToPa = 6894.76; % Pa psi-1
uAhToCoulumbs = 0.0036;
cellArea = 2e-4; % Elecrode area, m^2
mwCO2 = 44; % Molecular weight of CO2, g/mol
R = 8.314 * 1e3; % Gas const ,  (L * Pa)/(mol * K)
T = 298; % Temperature, assume RT = 298K
F = 96485; % Faraday Constant, (s * A)/mol
secToHr = 1/3600; % 

% Key Metrics calculated from Raw Data
numCycles = max(rawData(:, 5)); 
mcVals = zeros(numCycles - 2, 9); % 20% 30% 40% 50% 60% 70 80 90 Full
rawFluxVals = zeros(numCycles - 2, 10); % Max 20% 30% 40% 50% 60% 70 80 90 Full
cc_Vals = zeros(numCycles - 2, 2); % uAh, Coulumbs
fe_Vals = zeros(numCycles - 2, 1);

% Flux Values to Compare Against
% Max 20% 30% 40% 50% 60% 70 80 90 Full
procFluxVals = procData(:, 2:end);
% n_fluxVals = zeros(numCycles - 2, 10); 
% i_fluxVals = zeros(numCycles - 2, 10); 
% mm_fluxVals = zeros(numCycles - 2, 10);
% mmt_fluxVals = zeros(numCycles - 2, 10);

% n_dpVals = zeros(numCycles - 2, 10);
% mm_dpVals = zeros(numCycles - 2, 10);
% mmt_dpVals = zeros(numCycles - 2, 10);

percentages = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];

%% Calculating Values

% Loops through each cycle, skips the incomplete cycle 1
for i = 2:numCycles - 1
    
    % Extracts cycle data
    cycData = rawData(rawData(:, 5) == i, :); % Data in cycle i
    chargeRegionData = cycData(cycData(:, 7) > 0.00002, :);
    stepIndexVals = chargeRegionData(5,6);
    chargeRegionData = cycData(cycData(:, 7) > 0, :); % Data in charging region, theoretical > 0, realistic > 0.0002
    chargeRegionData = chargeRegionData(chargeRegionData(:, 6) == stepIndexVals, :);

    timeVals = chargeRegionData(:, 4); % Step time values, s
    pressureVals = chargeRegionData(:, 19); % Auxilary Presure Values, psi
    chargeVals = chargeRegionData(:, 10) * 1e6; % Charge Capacity, uAh

    % Create biexponential fit for pressure curve
    pFitCurve = fit(timeVals, pressureVals, 'exp2'); % Creates biexponential model of P vs t
    pFitEval = pFitCurve(timeVals); % Calculates pressures at given times on pFitCurve

    % Calculates Max Flux
    pDeriv = differentiate(pFitCurve, timeVals); % Differentiates pFitCurve with datapoints at x = timeVals
    gDeriv = -1 * psiToGCO2(pDeriv) / cellArea / secToHr; % Convert pDeriv from psi to g CO2
    rawFluxVals(i - 1, 1) = max(gDeriv);


    % Stores the values of the fitted exponential curve
    pFitTimePressure = [];
    pFitTimePressure(:, 1) = timeVals; % set column 1 to time values
    pFitTimePressure(:, 2) = pFitEval;

    % Calculates minimum pressure drops required for the following pressure calculations
    % 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, Full
    pStart = pFitEval(1);
    pEnd = pFitEval(end);
    dPFit = pStart - pEnd; % Pressure drop

    minP = pStart - dPFit * percentages; % Minimum pressure values for each avg flux calc

    % Calculates CC Vals
    ccuAh = max(chargeVals); % charge capacity, uAh
    ccCoulumbs = ccuAh * uAhToCoulumbs; % convert cc to coulumbs
    
    % Calculates FE Vals
    fe = 100 * (psiToGCO2(dPFit) / mwCO2) * F / (ccCoulumbs);
    
    % Assign values at end of loop
    cc_Vals(i - 1, 1) = ccuAh;
    cc_Vals(i - 1, 2) = ccCoulumbs;
    fe_Vals(i - 1) = fe;

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
        rawFluxVals(i - 1, j + 1) = flux;    
    end

    rawFluxVals = round(rawFluxVals, 2);

    diffProc = 100 * (rawFluxVals - procFluxVals) ./ procFluxVals;
    % diffProc =  round(diffProc);

    %{

    % % Calculate Normalized Flux Values 
    % n_vals = normalize(pressureVals); % normalizes raw pressure values
    % n_pFitCurve = fit(timeVals, n_vals, 'exp2'); % Creates biexponential model of normalized P vs t
    % n_pFitVals = n_pFitCurve(timeVals); % Calculates normalized pressure at each time
    % n_TimePressure = []; 
    % n_TimePressure(:, 1) = timeVals; % set col 1 to time values
    % n_TimePressure(:, 2) = n_pFitVals; % set col 2 to normalized p vals
    % 
    % % Calculates minimum pressure drops required for the following pressure calculations
    % n_pStart = n_pFitVals(1);
    % n_pEnd = n_pFitVals(end);
    % n_dPFit = n_pStart - n_pEnd; % Pressure drop
    % n_minP = n_pStart - n_dPFit * percentages; % Pressures at x%
    % 
    % n_pData = []; % change this back later
    % for j = 1:9
    %     n_pData =  n_TimePressure(n_TimePressure(:, 2) > n_minP(j), :); % Constrains dataset to have minimum threshold pressure
    %     n_pStart = n_pData(1, 2);
    %     n_pEnd = n_pData(end, 2);
    %     n_tStart = n_pData(1, 1);
    %     n_tEnd = n_pData(end, 1);
    %     n_dP = n_pStart - n_pEnd; % pressure drop
    %     n_dt = n_tEnd - n_tStart;
    % 
    %     %{
    %     How should I normalize the difference in pressure?  
    %     If I leave it as is, the n_dP is too big, causing it to overshoot
    %     flux vals
    %     %}
    % 
    %      % Calculates dP (grams CO2)
    %     n_dPGramsCO2 = psiToGCO2(n_dP); % Convert dp from psi to g
    % 
    %     % Calculates Flux
    %     n_dtHours = n_dt * secToHr; % Convert dt from s to h
    %     n_flux = n_dPGramsCO2 / (cellArea * n_dtHours); % calculate cycle flux, gCo2 / h
    % 
    %     % Records final calculations
    %     n_dpVals(i - 1, j + 1) = n_dP;
    %     n_fluxVals(i - 1, j + 1) = n_flux; 
    % end
    % 
    % % Calculates Normalized Max Flux 
    % n_pDeriv = differentiate(n_pFitCurve, timeVals); % Differentiates pFitCurve with datapoints at x = timeVals
    % n_gDeriv = -1 * psiToGCO2(n_pDeriv) / cellArea / secToHr; % Convert pDeriv from psi to g CO2
    % n_fluxVals(i - 1, 1) = max(n_gDeriv);
    % 
    % % Calculate Interpolated Flux Values
    % % pStart = pFitEval(1);
    % 
    % %{
    % Theory: 
    % 1. Utilize pFitEval to obtain pStart (same as standard version)
    % 2. Calculate minP values (same as standard version)
    % 3. Find inverse of pFitCurve
    % 4. Plug in minP values to get the time start and finishes
    % 5. Calculate dt and dp
    % 6. Convert dp from psi to g Co2
    % 7. Calculate flux
    % %}
    % 
    % % Inverse equation formula = a*exp(b*x) + c*exp(d*x)
    % % pressure80 = minP(7);
    % % pFitInv = approximateInverseFit(pFitCurve, pressure80);
    % 
    % 
    % 
    % % Calculate Moving Mean (Centered) Flux Values
    % mm_Vals = movmean(pressureVals, 3); % Applies moving mean to presssureVals
    % mm_FitCurve = fit(timeVals, mm_Vals, 'exp2'); % Creates biexponential model of mmPvals
    % mm_FitVals = mm_FitCurve(timeVals); % Calculates mm pressure at each time
    % mm_TimePressure = []; 
    % mm_TimePressure(:, 1) = timeVals; % set col 1 to time values
    % mm_TimePressure(:, 2) = mm_FitVals; % set col 2 to normalized p vals
    % 
    %  % Calculates minimum pressure drops required for the following pressure calculations
    % mm_Start = mm_FitVals(1);
    % mm_End = mm_FitVals(end);
    % mm_dPFit = mm_Start - mm_End; % Pressure drop
    % mm_minP = mm_Start - mm_dPFit * percentages; % Pressures at x%
    % 
    % 
    % mm_PData = [];
    % for j = 1:9
    %     mm_PData =  mm_TimePressure(mm_TimePressure(:, 2) > mm_minP(j), :); % Constrains dataset to have minimum threshold pressure
    %     mm_PStart = mm_PData(1, 2);
    %     mm_PEnd = mm_PData(end, 2);
    %     mm_tStart = mm_PData(1, 1);
    %     mm_tEnd = mm_PData(end, 1);
    %     mm_dP = mm_PStart - mm_PEnd; % pressure drop
    %     mm_dt = mm_tEnd - mm_tStart;
    % 
    %      % Calculates dP (grams CO2)
    %     mm_dPGramsCO2 = psiToGCO2(mm_dP); % Convert dp from psi to g
    % 
    %     % Calculates Flux
    %     mm_dtHours = mm_dt * secToHr; % Convert dt from s to h
    %     flux = mm_dPGramsCO2 / (cellArea * mm_dtHours); % calculate cycle flux, gCo2 / h
    % 
    %     % Records final calculations
    %     mm_dpVals(i - 1, j + 1) = mm_dP;
    %     mm_fluxVals(i - 1, j + 1) = flux; 
    % 
    % end
    % 
    % % Calculates Moving Mean (Centered) Max Flux 
    % mm_pDeriv = differentiate(mm_FitCurve, timeVals); % Differentiates pFitCurve with datapoints at x = timeVals
    % mm_gDeriv = -1 * psiToGCO2(mm_pDeriv) / cellArea / secToHr; % Convert pDeriv from psi to g CO2
    % mm_fluxVals(i - 1, 1) = max(mm_gDeriv);
    % 
    % 
    % % Calculate Moving Mean (Trailing) Flux Values
    % mmt_Vals = movmean(pressureVals, [2 0]); % Applies moving mean to presssureVals
    % mmt_FitCurve = fit(timeVals, mmt_Vals, 'exp2'); % Creates biexponential model of mmPvals
    % mmt_FitVals = mm_FitCurve(timeVals); % Calculates mm pressure at each time
    % mmt_TimePressure = []; 
    % mmt_TimePressure(:, 1) = timeVals; % set col 1 to time values
    % mmt_TimePressure(:, 2) = mmt_FitVals; % set col 2 to normalized p vals
    % 
    %  % Calculates minimum pressure drops required for the following pressure calculations
    % mmt_Start = mmt_FitVals(1);
    % mmt_End = mmt_FitVals(end);
    % mmt_dPFit = mmt_Start - mmt_End; % Pressure drop
    % mmt_minP = mmt_Start - mmt_dPFit * percentages; % Pressures at x%
    % 
    % 
    % mm_PData = [];
    % for j = 1:9
    %     mmt_PData =  mmt_TimePressure(mmt_TimePressure(:, 2) > mmt_minP(j), :); % Constrains dataset to have minimum threshold pressure
    %     mmt_PStart = mmt_PData(1, 2);
    %     mmt_PEnd = mmt_PData(end, 2);
    %     mmt_tStart = mmt_PData(1, 1);
    %     mmt_tEnd = mmt_PData(end, 1);
    %     mmt_dP = mmt_PStart - mmt_PEnd; % pressure drop
    %     mmt_dt = mmt_tEnd - mmt_tStart;
    % 
    %      % Calculates dP (grams CO2)
    %     mmt_dPGramsCO2 = psiToGCO2(mmt_dP); % Convert dp from psi to g
    % 
    %     % Calculates Flux
    %     mmt_dtHours = mmt_dt * secToHr; % Convert dt from s to h
    %     mmt_flux = mmt_dPGramsCO2 / (cellArea * mmt_dtHours); % calculate cycle flux, gCo2 / h
    % 
    %     % Records final calculations
    %     mmt_dpVals(i - 1, j + 1) = mmt_dP;
    %     mmt_fluxVals(i - 1, j + 1) = mmt_flux; 
    % 
    % end
    % 
    % % Calculates Moving Mean (Trailing) Max Flux
    % mmt_pDeriv = differentiate(mmt_FitCurve, timeVals); % Differentiates pFitCurve with datapoints at x = timeVals
    % mmt_gDeriv = -1 * psiToGCO2(mmt_pDeriv) / cellArea / secToHr; % Convert pDeriv from psi to g CO2
    % mmt_fluxVals(i - 1, 1) = max(mmt_gDeriv);


    
    %}

end


fprintf("Script ended" + "\n");

%% Extraneous Functions

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


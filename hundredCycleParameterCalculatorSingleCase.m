function [fluxData, chargeCapacity, dischargeCapacity, ...
        coulombicEfficiency, massCaptureData, massReleaseData, ...
        timescalesCapture, timescalesRelease, timescalesCharge, ...
        timescalesDischarge, averageChargeRates, averageDischargeRates, ...
        faradaicEfficiency, fluxChargeVoltageDep,generalChargeVoltageDep] ...
        = hundredCycleParameterCalculatorSingleCase(dataFolder, internalVolume)

    %clear
    %clc

    %dataFolder = uigetdir('G://Shared drives/Data/Arbin data/_Sealed Cell Data/_TwoConcentrationFlux/_ProcessedData/');
    date = string(datetime("today", "Format", "uuuu-MM-dd"));
    plotDataFolder = fullfile(dataFolder, date);
    cd(dataFolder);
    % will this automaticlaly create the folder within the directory the script
    % is being run?
    mkdir(plotDataFolder, "Degradation Performance Plots");
    mkdir(plotDataFolder, "Cycle Performance Plots");
    mkdir(plotDataFolder, "Processed Datasheet");

    cyclePerformancePlotsLocation = fullfile(plotDataFolder, "Cycle Performance Plots/");
    degradationPerformancePlotsLocation = fullfile(plotDataFolder, "Degradation Performance Plots/");

    mkdir(cyclePerformancePlotsLocation, "Flux Plots");
    mkdir(cyclePerformancePlotsLocation, "Equivalent Current Plots");
    mkdir(cyclePerformancePlotsLocation, "Capacity Curve Plots");

    fluxPlotsLocation = fullfile(cyclePerformancePlotsLocation, "Flux Plots");
    equivCurrPlotsLocation = fullfile(cyclePerformancePlotsLocation, "Equivalent Current Plots");
    capacityCurvePlotsLocation = fullfile(cyclePerformancePlotsLocation, "Capacity Curve Plots");

    s = dir('*.xlsx'); % get the data files info
    scell = struct2cell(s); % converts s to cell for easier extraction
    rawData = [];

    for i = 1:size(s,1)
        % Read the raw data from sheet 2 of the xlsx files;
        fprintf("Reading file %d...\n", i);
        % For each xlsx, read data
        data = xlsread(scell{1,i},2);
        spreadsheetHeight = height(data);
        newData = []; % reset new data for each xlsx
        newData = data(2:spreadsheetHeight,:);
        rawData = [rawData; newData]; % append the data from each xlsx file to the end of the data variable
    end

    %% Initialize constants for future calculations
    psiToPascal = 6894.76;
    pascalToMol = double(internalVolume) / (298 * 8.314) / 1000000;
    molToGramsCO2 = 44;
    cmSquareToMSquare = 1/10000;
    electrodeSurfaceAreaCMSquare = 2;
    perSecondToHour = 3600;
    uAhToCoulombs = 0.0036;
    faradayConstant = 96485;

    %% Calculate key parameters by cycle
    %{
Select data

Input file structure, numbered columns:

1. Data Point | 2. Date_Time (MM/DD/YYYY HH:MM:SS) | 3. Test_Time (s) | 4. Step_Time (s)| 5. Cycle_Index 

6. Step_Index | 7. Current (A) | 8. Voltage (V) | 9. Power (W) | 10. Charge Capacity (Ah)

11. Discharge Capacity (Ah)| 12. Charge_Energy (Wh)| 13. Discharge_Energy (Wh)| 14. ACR (Ohm)| 15. dV/dt

16. Internal Resistance (Ohm) | 17. dQ/dV |18. dV/dQ | 19. Aux_Pressure_1 (PSI) | 20. Aux_Pressure/dt_1 (PSI)|

    %}

    % Initialize matrices for parfor loop
    numberCycles = max(rawData(:,5));
    cycles = zeros([numberCycles - 2]);
    capturePressureDropsPascal = zeros([(numberCycles - 2) 9]);
    timescalesCapture = zeros([(numberCycles - 2) 9]);
    massCaptureData = zeros([(numberCycles - 2) 9]);
    fluxDataFractions = zeros([(numberCycles - 2) 9]);
    chargeCapacityMicroAmpHour = zeros([(numberCycles - 2) 1]);
    chargeCapacityCoulombs = zeros([(numberCycles - 2) 1]);
    dischargeCapacityMicroAmpHour = zeros([(numberCycles - 2) 1]);
    dischargeCapacityCoulombs = zeros([(numberCycles - 2) 1]);
    coulombicEfficiency = zeros([(numberCycles - 2) 1]);
    averageChargeRates = zeros([(numberCycles - 2) 9]);
    averageDischargeRates = zeros([(numberCycles - 2) 9]);
    faradaicEfficiency = zeros([(numberCycles - 2) 1]);
    massReleaseData = zeros([(numberCycles - 2) 9]);
    chargeVoltage = zeros([(numberCycles - 2) 1]);
    dischargeVoltage = zeros([(numberCycles - 2) 1]);
    chargingEndTime = zeros([(numberCycles - 2) 1]);
    dischargingEndTime = zeros([(numberCycles - 2) 1]);
    timescalesRelease = zeros([(numberCycles - 2) 1]);
    timescalesCharge = zeros([(numberCycles - 2) 1]);
    timescalesDischarge = zeros([(numberCycles - 2) 1]);

    % Execute loop for each cycle
    parfor i = 2:max(rawData(:,5) - 1)

        cycles(i - 1, 1) = i - 1;
        cycleData = rawData(rawData(:,5) == i, :);

        chargingData = cycleData(cycleData(:, 7) > 0.00002, :);
        chargeVoltage(i - 1)  = round(median(chargingData(:, 8)), 1);
        stepIndexCharging = chargingData(5, 6);
        chargingData = cycleData(cycleData(:, 7) > 0, :);
        chargingData = chargingData(chargingData(:, 6) == stepIndexCharging, :);

        dischargingData = cycleData(cycleData(:, 7) < -0.00002, :);
        dischargeVoltage(i - 1, 1) = round(median(dischargingData(:, 8)), 1);

        % Collect time values for time-dependent analysis
        chargingStartTime(i - 1, 1) = min(cycleData(:, 3));
        chargingEndTime(i - 1, 1) = max(chargingData(:, 3));

        %% Calculate the capture, timescale, and flux at various percent capacities for each cycle.
        %{

           Method of calculation:
            - Select cycle data, identifying index in row
            - Select cycle data for the capture phase, identifying any positive
                current > 10 uA to avoid sensitivity errors 
            - Shift data for the capture cycle to set the initial time to 0,
                required for exponential fits
            - Generate a two-term (biexponential) fit to the data. Use of a
                one-term fit is generally ineffective, as a constant offset is
                required and the Matlab exp1 fit type does not include a constant
                parameter.
            - Create a variable and store function values over the timeframe of
                the cycle.
            - Evaluate the total pressure drop and the time taken to get to
                different percents of the total drop.
            -Convert pressure in psi to grams of CO2, referencing recorded volume
                measurements for each device.
        %}

        pressureFitCurve = fit(chargingData(:, 4), chargingData(:, 19), 'exp2');
        pressureFitCurveDerivative = differentiate(pressureFitCurve, chargingData(:, 4));
        %pressureFitDerivativeEval = pressureFitCurveDerivative(chargingData(:, 4));

        gramFitDerivativeEval = -1 * pressureFitCurveDerivative * psiToPascal *  ...
            pascalToMol * molToGramsCO2 * perSecondToHour / ...
            (electrodeSurfaceAreaCMSquare * cmSquareToMSquare);
        maxFlux(i - 1, 1) = max(gramFitDerivativeEval);

        % plot(chargingData(:,4), chargingData(:,19));
        % hold on
        % plot(pressureFitCurve);
        pressureFitEval = pressureFitCurve(chargingData(:, 4));

        pressureFitTimePressure = [];

        pressureFitTimePressure(:, 1) = chargingData(:, 4);
        pressureFitTimePressure(:, 2) = pressureFitEval;

        startPressure =  pressureFitEval(1);
        endPressure = pressureFitEval(end);

        pressureCaptured = startPressure - pressureFitEval;
        gramsPerM2Captured = pressureCaptured * psiToPascal * pascalToMol ...
            * molToGramsCO2 / (electrodeSurfaceAreaCMSquare * cmSquareToMSquare);

        % CODE ClEAN UP: is this if statement necessary, since
        % captureFailure is never used
        if endPressure > startPressure
            captureFailure = 1;
        else
            captureFailure = 0;
        end

        totalPressureDrop = startPressure - endPressure;
        thresholdPressure = zeros([1 9]);
        for k = 2:10
            fractionCalc = k/10;
            index = k - 1;
            threshold = startPressure - totalPressureDrop * fractionCalc;

            thresholdPressure(index) = threshold;
        end

        for j = 1:9

            fractionCycleData = [];
            fractionStartPressure = [];
            fractionEndPressure = [];
            fractionStartTime = [];
            fractionEndTime = [];
            try
                fractionCycleData = pressureFitTimePressure(pressureFitTimePressure(:, 2) > thresholdPressure(j), :);
                fractionStartPressure = fractionCycleData(1, 2);
                fractionEndPressure = fractionCycleData(end, 2);
                fractionStartTime = fractionCycleData(1, 1);
                fractionEndTime = fractionCycleData(end, 1);

                fractionPressureDrop = fractionStartPressure - fractionEndPressure;
                fractionTimeDifference = fractionEndTime - fractionStartTime;

                capturePressureDropsPascal(i - 1, j) = fractionPressureDrop * 6894.76;
                timescalesCapture(i - 1, j) = fractionTimeDifference;
                massCaptureData(i - 1, j) = fractionPressureDrop*psiToPascal * pascalToMol ...
                    * molToGramsCO2 / (electrodeSurfaceAreaCMSquare * cmSquareToMSquare);

                fluxDataFractions(i - 1, j) = fractionPressureDrop * psiToPascal * pascalToMol ...
                    * molToGramsCO2 / (electrodeSurfaceAreaCMSquare * cmSquareToMSquare) ...
                    / (fractionTimeDifference) * perSecondToHour;
            catch
                continue
            end
        end

        %% Collect the charge capacity and discharge capacity for each cycle. Also calculate average charging rates.
        %{

    Method of calculation:
    - Charge, discharge capacities provided by Arbin. They are the
    integrated current over time.
    - Average charging rate is found by determining time taken to get to a
    given percent charge and dividing the charge by the time.
    
    Input file structure, numbered columns:

    1. Data Point | 2. Date_Time (MM/DD/YYYY HH:MM:SS) | 3. Test_Time (s) | 4. Step_Time (s)| 5. Cycle_Index 

    6. Step_Index | 7. Current (A) | 8. Voltage (V) | 9. Power (W) | 10. Charge Capacity (Ah)

    11. Discharge Capacity (Ah)| 12. Charge_Energy (Wh)| 13. Discharge_Energy (Wh)| 14. ACR (Ohm)| 15. dV/dt

    16. Internal Resistance (Ohm) | 17. dQ/dV |18. dV/dQ | 19. Aux_Pressure_1 (PSI) | 20. Aux_Pressure/dt_1 (PSI)|
        %}

        chargeCapacityMicroAmpHour(i - 1, 1) = max(cycleData(:, 10)) * 1000000;
        chargeCapacityCoulombs(i - 1, 1) = max(cycleData(:, 10)) * 1000000 * uAhToCoulombs;

        dischargeCapacityMicroAmpHour(i - 1, 1) = max(cycleData(:, 11)) * 1000000;
        dischargeCapacityCoulombs(i - 1, 1) = max(cycleData(:, 11)) * 1000000 * uAhToCoulombs;

        coulombicEfficiency(i - 1, 1) = 100 * dischargeCapacityCoulombs(i - 1,1) ...
            / chargeCapacityCoulombs(i - 1, 1);
        faradaicEfficiency(i - 1, 1) = 100 * totalPressureDrop * psiToPascal * ...
            pascalToMol * faradayConstant / (chargeCapacityCoulombs(i - 1, 1));

        % Get average charging rates
        thresholdCharge = zeros([1 9]);
        for k = 2:10
            fractionCalc = k / 10;
            index = k - 1;
            thresholdC = chargeCapacityMicroAmpHour(i - 1, 1) * fractionCalc;

            thresholdCharge(index) = thresholdC;
        end
        thresholdDischarge = zeros([1 9]);
        for k = 2:10
            fractionCalc = k / 10;
            index = k - 1;
            thresholdD = dischargeCapacityMicroAmpHour(i - 1, 1) * fractionCalc;

            thresholdDischarge(index) = thresholdD;
        end

        for k = 1:9
            fraction = (k + 1) / 10;
            chargingDataFraction = cycleData(cycleData(:, 10) < thresholdCharge(k), :);
            averageChargeRates(i - 1, k) = max(chargingDataFraction(:, 10)) / ...
                max(chargingDataFraction(:, 4)) * 1000000 * uAhToCoulombs * 1000;

            dischargingDataFraction = cycleData(cycleData(:, 11) < thresholdDischarge(k), :);
            dischargingDataFraction = dischargingDataFraction(dischargingDataFraction(:, 7) < 0, :);
            averageDischargeRates(i-1,k) = max(dischargingDataFraction(:, 11)) ... 
                / max(dischargingDataFraction(:, 4)) * 1000000 * uAhToCoulombs * 1000;
        end

        %% Collect data on timescales of release. Currently lower priority, but better to generate a backlog.
        %{



        %}

        %% Generate plots of key per cycle information, such as actual pressure drop, flux, and model parameters


        if mod(i, 10) == 0
            % Plot pressure, model, and flux on graph
            fileName = fluxPlotsLocation + "\Cycle " + string(i) + "FluxPlot.png";
            figure1 = figure("Visible","off");
            yyaxis left
            plot(chargingData(:, 4),chargingData(:, 19), "Marker","o", "Color","r", "MarkerSize", 1);
            ylabel("Pressure (psi)");
            xlabel("Cycle Time (s)");
            hold on;
            plot(chargingData(:, 4), pressureFitEval);
            yyaxis right
            ylabel("Flux (g/m2-hr)");
            plot(chargingData(:, 4), gramFitDerivativeEval);
            legend("Actual Pressure (psi)", "Model Pressure (psi)", "Flux (g/m2-hr)", "Location", "eastoutside");
            saveas(figure1, fileName, "png")
            close(figure1);

            % Plot equivalent current versus actual current
            fileNameEquiv = equivCurrPlotsLocation + "\Cycle " + string(i) + "EquivCurrentPlot.png";
            equivCurrent = -1 * pressureFitCurveDerivative * psiToPascal * pascalToMol * faradayConstant * 1000;
            chargingMilliAmp = chargingData(:,7) * 1000;
            figure2 = figure("Visible", "off");
            yyaxis left
            plot(chargingMilliAmp, equivCurrent);
            ylabel("Equivalent Current (mA)")
            xlabel("Actual Current (mA)")
            ylim([0 max(chargingMilliAmp)]);
            xlim([0 max(chargingMilliAmp)]);
            hold on
            plot(chargingMilliAmp, chargingMilliAmp, "Color", "magenta");
            saveas(figure2, fileNameEquiv, "png");
            close(figure2);

            % Plot average flux and time versus capacity
            fileNameCapacity  = capacityCurvePlotsLocation + "\Cycle " + string(i) + "CapacityTimePlot.png";
            averageFlux = (gramsPerM2Captured./chargingData(:, 4)) * perSecondToHour;
            figure3 = figure("Visible", "off");
            yyaxis left
            plot(gramsPerM2Captured, averageFlux);
            ylabel("Average Flux (g/m2-hr)")
            xlabel("Capacity (g/m2)")
            yyaxis right
            plot(gramsPerM2Captured,chargingData(:, 4));
            ylabel("Cycle Time (s)")
            legend("Average Flux (g/m2-hr)", "Cycle Time (s)", "Location", "eastoutside");
            saveas(figure3, fileNameCapacity, "png");
            close(figure3)

        end
    end

    % Correct for discharging time
    dischargeLength = zeros([height(chargingStartTime) 1]);
    dischargeLengthCorrection = zeros([height(chargingStartTime) 1]);
    for i = 2:length(chargingStartTime)
        dischargeLength(i) = chargingStartTime(i) - chargingEndTime(i - 1);
        dischargeLengthCorrection(i) = sum(dischargeLength);
    end

    chargingStartTimeCorrected = chargingStartTime - dischargeLengthCorrection;

    % Round matrices
    maxFlux = round(maxFlux, 2);
    fluxDataFractions = round(fluxDataFractions, 2);
    chargeCapacityMicroAmpHour = round(chargeCapacityMicroAmpHour, 4);
    chargeCapacityCoulombs = round(chargeCapacityCoulombs, 4);
    dischargeCapacityMicroAmpHour = round(dischargeCapacityMicroAmpHour, 4);
    dischargeCapacityCoulombs = round(dischargeCapacityCoulombs, 4);
    coulombicEfficiency = round(coulombicEfficiency, 2);
    faradaicEfficiency = round(faradaicEfficiency, 2);
    massCaptureData = round(massCaptureData, 4);

    % Add cycle column
    fluxData = [cycles maxFlux fluxDataFractions];
    chargeCapacity = [cycles chargeCapacityMicroAmpHour chargeCapacityCoulombs];
    dischargeCapacity = [cycles dischargeCapacityMicroAmpHour dischargeCapacityCoulombs];
    coulombicEfficiency = [cycles coulombicEfficiency];
    faradaicEfficiency = [cycles faradaicEfficiency];
    averageChargeRates = [cycles averageChargeRates];
    averageDischargeRates = [cycles averageDischargeRates];
    chargingStartTime = [cycles chargingStartTime];
    massCaptureData = [cycles massCaptureData];
    sqrtChargingTime = sqrt(chargingStartTime(:, 2));
    sqrtChargingTimeCorrected = sqrt(chargingStartTimeCorrected);
    effectiveMicroMol = chargeCapacity(:, 3) / 96485 * 1000000;

    fluxChargeVoltageDep = [cycles chargeVoltage maxFlux fluxDataFractions];
    generalChargeVoltageDep = [cycles chargeVoltage dischargeVoltage maxFlux ...
        chargeCapacityCoulombs dischargeCapacityCoulombs coulombicEfficiency(:, 2) faradaicEfficiency(:, 2)];


    %% Generate plots of key over cycles information, demonstrating decay etc.

    % Flux
    for i = 2:10
        fluxTitle = degradationPerformancePlotsLocation + "\" + string(i * 10)+ " Percent Average Flux vs Cycle.png";
        yLabeler = string(i * 10) + "% Average Flux (g/m2-hr)";
        fluxCycle = figure("Visible","off");
        plot(fluxData(:,1), fluxData(:, i + 1), "Marker", "o", "Color", "b");
        ylabel(yLabeler);
        ylim([0 (max(fluxData(:,i + 1) + 10))])
        xlabel("Cycle");
        saveas(fluxCycle, fluxTitle, "png");
        close(fluxCycle);
    end

    % Charge Capacity
    chargeTitle = degradationPerformancePlotsLocation + "\" + "Charge Capacity vs Cycle.png";
    chargeCycle = figure("Visible","off");
    plot(chargeCapacity(:,1), chargeCapacity(:, 3))
    ylabel("Charge Capacity (C)")
    xlabel("Cycle")
    saveas(chargeCycle, chargeTitle, "png")
    close(chargeCycle);

    % Charge Capacity as f(t)
    chargeTimeTitle = degradationPerformancePlotsLocation + "\" + "Charge Capacity vs Total Time.png";
    chargeTimePlot = figure("Visible","off");
    plot(chargingStartTime(:,2), chargeCapacity(:, 3))
    ylabel("Charge Capacity (C)")
    xlabel("Time (s)")
    saveas(chargeTimePlot, chargeTimeTitle, "png")
    close(chargeTimePlot);

    % Charge Capacity as f(tcorrected)
    chargeTimeCTitle = degradationPerformancePlotsLocation + "\" + "Charge Capacity vs Charged Time.png";
    chargeTimeCPlot = figure("Visible","off");
    plot(chargingStartTimeCorrected, chargeCapacity(:, 3))
    ylabel("Charge Capacity (C)")
    xlabel("Time (s in charged state)")
    saveas(chargeTimeCPlot, chargeTimeCTitle, "png")
    close(chargeTimeCPlot);

    % Discharge Capacity
    dischargeTitle = degradationPerformancePlotsLocation + "\" + "Discharge Capacity vs Cycle.png";
    dischargeCycle = figure("Visible","off");
    plot(dischargeCapacity(:, 1), dischargeCapacity(:, 3))
    ylabel("Discharge Capacity (C)")
    xlabel("Cycle")
    saveas(dischargeCycle, dischargeTitle, "png")
    close(dischargeCycle)

    % Coulombic Efficiency
    coulombTitle = degradationPerformancePlotsLocation + "\" + "Coulombic Efficiency vs Cycle.png";
    coulombCycle = figure("Visible","off");
    plot(coulombicEfficiency(:, 1), coulombicEfficiency(:, 2));
    ylim([(min(coulombicEfficiency(:, 2)) - 10) 100]);
    ylabel("Coulombic Efficiency (%)", FontName="Arial")
    xlabel("Cycle")
    saveas(coulombCycle, coulombTitle, "png");
    close(coulombCycle);

    % Faradaic Efficiency
    faradTitle = degradationPerformancePlotsLocation + "\" + "Faradaic Efficiency vs Cycle.png";
    faradCycle = figure("Visible","off");
    plot(faradaicEfficiency(:,1), faradaicEfficiency(:, 2));
    ylabel("Faradaic Efficiency (%)")
    xlabel("Cycle")
    saveas(faradCycle, faradTitle, "png");
    close(faradCycle);

    % Flux vs Voltage
    chargeTitle = degradationPerformancePlotsLocation + "\" + "50% Average Flux vs Voltage.png";
    chargeVoltageFig = figure("Visible", "off");
    yyaxis left
    plot(chargeVoltage, fluxData(:, 9), "Marker", "o");
    ylabel("80% Average Flux (g/m2-hr)");
    xlabel("Voltage (V)");
    xlim([(min(chargeVoltage) - 0.2) (max(chargeVoltage) + 0.2)])
    saveas(chargeVoltageFig, chargeTitle, "png");
    close(chargeVoltageFig)

    % Charge Capacity vs Voltage
    chargeTitle = degradationPerformancePlotsLocation + "\" + "Charge Capacity vs Voltage.png";
    chargeVoltageFig = figure("Visible","off");
    yyaxis left
    plot(chargeVoltage, chargeCapacity(:, 3), "Marker", "o");
    ylabel("Charge Capacity (C)");
    xlabel("Voltage (V)");
    xlim([(min(chargeVoltage)-0.2) (max(chargeVoltage) + 0.2)])
    saveas(chargeVoltageFig, chargeTitle, "png");
    close(chargeVoltageFig)
end














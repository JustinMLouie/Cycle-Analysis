clear;
clc;
%% File Loading

% Import arbin file 
% folder = 'G:\Shared drives\Data\Arbin data\_Sealed Cell Data\_Cycling\_ProcessedData';
% [file, path] = uigetfile('.xlsx', 'Select Folder', folder);
[rawFile, rawPath] = uigetfile('.xlsx');
raw_filename = string(fullfile(rawPath, rawFile));

[proc_file, processed_path] = uigetfile('.xlsx');
proc_filename = string(fullfile(processed_path, proc_file));

fprintf("Justin's version")

% Import raw data, test parameters
raw = readmatrix(raw_filename, 'Sheet', 2);
fprintf("\n")
fprintf(raw_filename)
params = f_comment(raw_filename); 

% Import processed data, used to compare final values and obtain v0
proc_sum = readmatrix(proc_filename, 'Sheet', 1); % processed summary workbook sheet
proc_flux = readmatrix(proc_filename, 'Sheet', 2); % processed flux data
proc_mc = readmatrix(proc_filename, 'Sheet', 3); % processed mass capture
proc_cc = readmatrix(proc_filename, 'Sheet', 4); % processed charge capacity
proc_dc = readmatrix(proc_filename, 'Sheet', 5); % processed discharge capacity
proc_ce = readmatrix(proc_filename, 'Sheet', 6); % processed coulumbic efficiency 
proc_fe = readmatrix(proc_filename, 'Sheet', 7); % processed faradaic capacity

% Assign raw data to matrix columns
ex = 2:height(raw); % Create matrix height
data = raw(ex, [3, 7, 10, 16, 19]); % Pulls data from raw at selected columns
cyc = raw(ex, 5); % cycle number

% Detrend data
% FIGURE 1
% TODO: add labels to graph
figure();
plot(data(:, 5), LineWidth=1); % check for breakpoints
bp = [];
data(:, 6) = detrend(data(:, 5), 10, bp, Continuous=false); % detrends the data
plot(data(:, 6)); % check detrended data for kinks

% Column key
% [1,    2,       3,      4,          5       ]
% [time, current, charge, resistance, pressure]

%% Calculate Cell Data
% Parameters
psi_pa = 6894.76; % Pa psi-1
F = 96485; % Faraday constant, C mol-1 
V0 = proc_sum(8,2); % cell volume, read from the processed data file
A = 2e-4; % elecrode area, m^2
MW_CO2 = 44.01; % molecular weight of CO2, g/mol
MW_Q = [590, 100]; % molecular weight of quinone polymers, g/mol
% for future usage, should MW_Q be abstracted?
R_T = 8.314 * 298; % gas const. scaled by 298K, J/mol
psi_gCO2 = 1e-6 * psi_pa * MW_CO2 * V0 / R_T; % gCO2 psi-1 (given RT, cell volume)

% Initialize pressure data arrays
cycles = cyc(end) - 2; % ignore the class cycle due to truncation
params_p.dp = zeros(cycles,1);
params_p.dm = zeros(cycles,1);        
params_p.flux_max = zeros(cycles, 1);
params_p.flux_avg = zeros(cycles, 6);

% Initialize charge data arrays
params_c.c_range = zeros(cycles, 1);
params_c.c_abs = zeros(cycles, 1);
params_c.j_max = zeros(cycles, 1);
params_c.c_fe = zeros(cycles, 1);
params_c.c_util = zeros(cycles, 1); % not used, can be removed

% Initialize model data arrays
params_m.beta = zeros(cycles, 4);
params_m.r2 = zeros(cycles, 1);

% Calculate hypothetical charge capacity for utilization
p_quinone = 0.103; % wt% quinone in electrode 
c_hyp = 1e-3 * 2 * F * p_quinone * params.wt_w / (MW_Q(1)); % hypothetical charge, C

% Loop through cycles
for i = 1:cycles
    % assign columns to explicit fariables for readability
    [t, j, c, ~, p] = data_ex(i + 1, cyc, data); 
    p_detrend = data((cyc == i + 1), 6);
    
    % Isolate capture region
    ind_cap = j > 0; 

    t_cap = t(ind_cap); % capture window, s
    p_cap = p_detrend(ind_cap); % capture pressure profile, psi
    c_cap = c(ind_cap); % charge capacity, uAh
    
    % Perform window averaging
    % TODO: Once means are more aligned, update movmean to compare 2-3
    % numbers
    p_avg = movmean(p_cap, 1); 
    c_avg = movmean(c_cap, 1);
    dp = range(p_avg); % total pressure drop, psi
    dc = range(c_avg); % charge capacity, uAh

    % Calculate array of 25, 50, 75, 80, 100% avg fluxes
    p_norm = (p_avg - min(p_avg)) ./ dp; % normalized pressure profile
    p_drops = 0.01 * [20, 40, 60, 80, 90]; % unitless percentages
    t_avg = f_findt(t_cap, p_norm, p_drops); % returns all average times for target drops
    temp = [t_avg, t_cap(end)];
    temp2 = t_cap(end);
    temp3 = [p_drops, 1];
    flux_avg = dp * [p_drops, 1] ./ ([t_avg, t_cap(end)]); % psi/s
    flux_avg = 3600 * psi_gCO2 * flux_avg / A; % convert flux to g/(m^2 * h)

    % Pressure "interpolation" using 20deg polynomial for max flux calculation
    p_fit = polyfit(t_cap, p_norm, 20);
    p_der = polyder(p_fit); % derivative
    flux_max = -min(polyval(p_der, t_cap(t_cap < t_avg(4)))); % max flux under 80% capture time, psi/s
    flux_max = 3600 * psi_gCO2 * dp * flux_max / A; % convert flux to g/(m^2 * h)

    try
        % Performs pressure fitting using 2-regime exponential
        % beta0 = [0.95, 0.05, mean(t_avg(2:3)), 50 * t_avg(2)]; % return
        % model params, residuals
        % r2 = 1 - sum(resid.^2) / sum((p_norm - mean(p_norm)).^2);
        [beta, gof] = fit(t_cap, p_norm, 'exp2');

        if beta.a < beta.c
            beta = [beta.c, -beta.d^-1, beta.a, -beta.b^-1];
        else
            beta = [beta.a, -beta.b^-1, beta.c, -beta.d^-1];
        end

        r2 = gof.adjrsquare;

    catch
        % set params to NaNs in case of error
        beta = NaN(1,4);
        r2 = NaN;
    end

    % Return related pressure quantitites of interest
    params_p.dp(i) = dp; % psi
    params_p.dm(i) = psi_gCO2 * dp / A; % g CO2
    params_p.flux_max(i) = flux_max; % g/(m^2 * h)
    params_p.flux_avg(i, :) = flux_avg;

    % Return related model parameters
    params_m.beta(i, :) = beta;
    params_m.r2(i) = r2;

    % Return related charge parameters
    params_c.c_range(i) = 3600 * dc; % C
    params_c.c_abs(i) = 3600 * c_avg(end); % C
    params_c.j_max(i) = 1e3 * max(j); % mA
    params_c.fe(i) = 100 * psi_gCO2 * F * dp / MW_CO2 / (3600 * c_avg(end)); % percentage
    params_c.util(i) = 3600 * c_avg(end) / c_hyp;

end

%% Plot Relevant Quanitities
% Cycle array
x = 2:cycles + 1;

% FIGURE 2: Plot flux_max
figure(); 
subplot(2, 1, 1);
plot(x, params_p.flux_max, LineWidth=1.5);
ylabel('Max Flux g/(m^2 * h)')

% Plot avg fluxes
subplot(2, 1, 2);

plot(x, params_p.flux_avg(:, 1), LineWidth=1.5); hold on 
plot(x, params_p.flux_avg(:, 2), LineWidth=1.5);
plot(x, params_p.flux_avg(:, 3), LineWidth=1.5);
plot(x, params_p.flux_avg(:, 4), LineWidth=1.5);
legend('20% Flux', '40% Flux', '60% Flux', '80% Flux')
xlabel('Cycle Number')
ylabel('Avg Flux g/(m^2 * h)')

% FIGURE 3: Plot charge capacity
figure(); 
subplot(2, 1, 1);
plot(x, params_c.c_abs, 'k', LineWidth=1.5);
ylabel('Charge Capacity C')

% Plot FE
subplot(2, 1, 2);
plot(x, params_c.fe, 'k', LineWidth=1.5);
xlabel('Cycle Number')
ylabel('Faradaic Efficiency')

%% Processed Data Comparison
proc_flux_max = proc_flux(:, 2); % max flux matrix
proc_flux20 = proc_flux(:, 3); % 20% flux matrix
proc_flux40 = proc_flux(:, 5); % 40% flux matrix
proc_flux60 = proc_flux(:, 7); % 60% flux matrix
proc_flux80 = proc_flux(:, 9); % 80% flux matrix

err_flux_max = percent_error(params_p.flux_max, proc_flux_max);
err_flux_20 = percent_error(params_p.flux_avg(:, 1), proc_flux20);
err_flux_40 = percent_error(params_p.flux_avg(:, 2), proc_flux40);
err_flux_60 = percent_error(params_p.flux_avg(:, 3), proc_flux60);
err_flux_80 = percent_error(params_p.flux_avg(:, 4), proc_flux80);

err_mc = percent_error(params_p.dm, proc_mc(:, 10));

err_fe = percent_error(params_c.fe', proc_fe(:, 2));

% TODO: fix charge cap calc
% err_cc = percent_error(c_abs , proc_charge_cap);

fprintf("Script completed")

% Calculates relative percent error
function [error] = percent_error(analyzedData, processedData) 
    error = 100 * (analyzedData - processedData) ./ processedData;
end

% Extracts cycle data into individual variables (one-liner for readability)
function [t_cyc, j, c, r, p] = data_ex(i, cyc, data)
    % i = cycle number
    % cyc = cycle array
    % data = data matrix

    % Performs data extraction
    data_cyc = data((cyc == i), :);

    % Assigns columns to explicit vars for readability
    t = data_cyc(2:end, 1); % time
    j = data_cyc(2:end, 2); % current
    c = data_cyc(2:end, 3); % charge
    r = data_cyc(2:end, 4); % resistance
    p = data_cyc(2:end, 5); % pressure
    % p_detrend = data_cyc(2:end, 6);

    % Zero-out time
    t_cyc = t - t(1);

end








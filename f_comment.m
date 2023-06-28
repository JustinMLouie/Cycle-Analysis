%% Extract electrode info from comments, adapted from "Import Data" module
function params = f_comment(file)
    % sets up the Import Options and imports the data
    opts = spreadsheetImportOptions("NumVariables", 5);
    
    % Specify sheet and range
    opts.Sheet = "Global_Info";
    opts.DataRange = "A5:E5";

    % Specify column names and types
    opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "SerialNumber"];
    opts.SelectedVariableNames = "SerialNumber";
    opts.VariableTypes = ["char", "char", "char", "char", "string"];

    % Specify variable properties
    opts = setvaropts(opts, ["Var1", "Var2", "Var3", "SerialNumber"], "WhitespaceRule", "preserve");
    %opts = setvaropts(opts, opts.VariablesNames, "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Var1", "Var2", "Var3", "SerialNumber"], "EmptyFieldRule", "auto");
    %opts = setvaropts(opts, opts.VariablesNames, "EmptyFieldRule", "auto");

    % Read sheet
    raw = readtable(file, opts, "UseExcel", false);
    raw = table2array(raw); 
    % could combine these two lines?
    % raw = table2array(readtable(file, opts, "UseExcel", false));

    % Split lines, extract relevant parameters
    raw = splitlines(raw);
    params.working = string(upper(sscanf(raw(7), 'WE: %s'))); % working ID
    params.working = string(upper(sscanf(raw(4), 'CE: %s'))); % working ID
    params.gas = upper(raw(1)); % gas ID
    params.il = upper(raw(2)); % IL ID
    params.wt_w = sscanf(raw(8), 'WT: %f'); % working weight, mg
    params.wt_c = sscanf(raw(5), 'WT: %f'); % counter weight, mg
    params.wt_il = string(sscanf(raw(10), 'IL = %s')); % IL amount 
    % question, what unit is IL amount in?
    params.op = string(sscanf(raw(11), 'Operator: %s')); % operator

    % assertion for empty parameters
    assert(all(~structfun(@isempty, params)), 'Not all parameters covered.')

end
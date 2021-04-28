%% Script to demonstrate using OSCILOS_brass via function calls

clear all
addpath('C:\Users\amacl\OneDrive - Imperial College London\UROP 2020\OSCILOS\V5_OSCILOS_brass\');

% Temperature values ranging from -20°C to 30°C
T = (-20:5:30) + 273.2;

% For each Temperature value
for i = 1:length(T)
    fprintf('Calculating T = %.1f K ...', T(i))
    
    EigVals{i} = OSCILOS_brass('T1', T(i), ...
        ... % Name of output file
        'eig_filename', "EigVals_T_" + num2str(T(i)) + ".txt", ...
        ... % Suppress other outputs
        'no_popups', true, ...
        'cl_out', false, ...
        'log_out', false ...
        );
    
    fprintf(' done\n');
end

% Extract 4th eigenvalue of each and convert to frequency in Hz
FreqsMode4 = imag( cellfun( @(x) x(4), EigVals) )./2./pi;

% Plot f_4 vs T
plot(T - 273.2, FreqsMode4, 'g', 'LineWidth', 2);
xlabel('T (°C)'); ylabel('f_4 (Hz)');

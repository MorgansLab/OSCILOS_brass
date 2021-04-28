%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PARSE CONFIGURATION FILE
%
% This subroutine loads the 'Config.txt' file and initialises the config
% options
% 
% Last update : 09/09/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CFG] = Config_subfc(pthstr,varargin)

	tic;
    
    OSCILOS_BRASS_VERSION = '5'; % string or char vector
    
    CFG.pthstr = pthstr;

    %% Option name definitions - complete list of config options

    Plot_opts = {
        'geom_plot',    'lgc';  % Draw Geometry plot
        'mf_plot',      'lgc';  % Draw Mean Flow plot
        'bc_plot',      'lgc';  % Draw Boundary Conditions plot
        'efp_plot',     'lgc';  % Draw Equivalent Fundamental Pitch plot
        'contour_plot', 'lgc';  % Draw eigenspace contour map
        'plot_modes',   'int';  % Number of mode shapes to draw
        'save_figs',    'lgc';  % Save .fig files - only plots which are drawn are saved
        'save_pdfs',    'lgc';  % Save .pdf files - only plots which are drawn are saved
        'small_plots',  'lgc';  % Reduce plot size to accommodate smaller screens e.g. laptop screens
        };

    Mean_flow_opts = {
        'p1',           'dbl';  % Mean pressure at inlet
        'T1',           'dbl';  % Mean Temperature at inlet
        'choice_M1_u1', 'int';  % Choice between Mach number (value: 1) or mean flow velocity (value: 2) at the inlet
        'M1_u1',        'dbl';  % Value of M1 or u1 as specified by choice_M1_u1
        'choice_gamma', 'int';  % Choice between a constant adiabatic index (value: 1) or changing with temperature (value: 2)
        };

    Inlet_BC_opts = {
        'inlet_type',   'int';  % Inlet Boundary Condition type ID - see docs
        'inlet_param1', 'dbl';  % Inlet BC parameter 1 - see docs
        'inlet_param2', 'dbl';  % Inlet BC parameter 2 - see docs
        'inlet_param3', 'dbl';  % Inlet BC parameter 3 - see docs
        };

    Outlet_BC_opts = {
        'outlet_type',  'int';  % Outlet Boundary Condition type ID - see docs
        'outlet_param1','dbl';  % Outlet BC parameter 1 - see docs
        'outlet_param2','dbl';  % Outlet BC parameter 2 - see docs
        'outlet_param3','dbl';  % Outlet BC parameter 3 - see docs
        };

    Scan_range_opts = {
        'min_freq',     'dbl';  % Minimum frequency, in Hertz              
        'max_freq',     'dbl';  % Maximum frequency, in Hertz          
        'number_freq',  'int';  % Number of frequencies
        'min_gr',       'dbl';  % Minimum growth rate, in Hertz
        'max_gr',       'dbl';  % Minimum growth rate, in Hertz
        'number_gr'     'int';  % Number of growth rates
        };

    Geom_opts = {
        'geom_filename','str';  % Geometry filename, defaults to Geometry.txt if unspecified
        };

    Output_opts = {
        'run_name',     'str';  % Name of output subfolder - leave blank to write directly to output folder (no overwrite warning if blank)
        'warn_ovwrt',   'lgc';  % Warning dialog asking user if overwriting existing output is OK
        'cl_out' ,      'lgc';  % Output progress updates to command line
        'log_out',      'lgc';  % Output progress updates to logfile
        'wait_bars',    'lgc';  % Display waitbars during and messagebox on finishing calculation
        'fin_msg',      'lgc';  % Display calculation finished message
        'no_popups',    'lgc';  % Disable all popups (useful for external function calls)
        'eig_filename', 'str';  % Eigenvalue list filename, defaults to Eigenvalues.txt if unspecified
        'log_filename', 'str';  % Log filename, defaults to Logfile.txt if unspecified
        };
    
    Conf_opts = {
        'config_filename','str';% Configuration filename - may only be set using Name Value argument pair
        };

    % The Opts cell array contains config options as parsed from the config file
    % It has dimensions 2 x n where n is the no. of config options listed above
    %   Opts{1,:} contains option names
    %   Opts(2,:} contains expected option types
    %   Opts(3,:} contains option values (pre-validation)
    Opts = [Plot_opts', Mean_flow_opts', Inlet_BC_opts', Outlet_BC_opts', Scan_range_opts', Geom_opts', Output_opts', Conf_opts'];
    Opts(3,:) = cell(1,size(Opts,2));
    
    %% Read values from file into Opts array
    
    for i = 1:2:nargin-1 % Check arguments for alternative config_filename
        if strcmpi(varargin{i},'config_filename')
            Opts(3,strcmpi(Opts(1,:), 'config_filename')) = varargin(i+1);
            break
        end
    end
    CFG.config_filename = fullfile(CFG.pthstr,'Inputs', valOpt(Opts,'config_filename','Config.txt'));
    fid = fopen(CFG.config_filename);
    if fid == -1
        error("Config file "+CFG.config_filename+" not found");
    end        

    while ~feof(fid) % for each line in file
        C_cell = textscan(fid,"%s %[^%\t\r\n ] %*[^\r\n]",1);
        % unknown parameter names and empty values are ignored
        if isempty(C_cell{1}) || isempty(C_cell{2}); continue; end
        Opts(3,strcmpi(Opts(1,:), C_cell{1}{1})) = C_cell{2};
    end

    fclose(fid);
    
    %% Override values from config file with Name Value argument pairs
    
    for i = 1:2:nargin-1
        optloc = strcmpi(Opts(1,:), varargin{i});
        if any(optloc)
            Opts(3,optloc) = varargin(i+1);
        else
            warning("Argument '%s' is not a valid parameter name - this will be ignored, along with the following argument '%s'",varargin{i},string(varargin{i+1}));
        end
    end

    %% Initialise output options

    CFG.OUT.CL_OUT          = valOpt(Opts, 'cl_out', true);         % Output progress updates to command line
    CFG.OUT.LOG_OUT         = valOpt(Opts, 'log_out', true);        % Output progress updates to Logfile.txt
    CFG.OUT.RUN_NAME        = valOpt(Opts, 'run_name', '');         % Name of folder to save results to
    CFG.OUT.WARN_OVWRT      = valOpt(Opts, 'warn_ovwrt', true);     % Warning dialog if run_name directory already exists
    CFG.OUT.WAIT_BARS       = valOpt(Opts, 'wait_bars', true);      % Display waitbars
    CFG.OUT.FIN_MSG         = valOpt(Opts, 'fin_msg', true);        % Display message when calculation finished
    CFG.OUT.NO_POPUPS       = valOpt(Opts, 'no_popups', false);     % Disable all popups and plots
    %% Start CL output

    if CFG.OUT.CL_OUT
        fprintf("OSCILOS_lite brass version %s\n",OSCILOS_BRASS_VERSION)
        fprintf("\tLoading config parameters... ");
    end

    %% Configure output directory
    
    if CFG.OUT.NO_POPUPS; CFG.OUT.WARN_OVWRT = false; end
    if CFG.OUT.WARN_OVWRT && (exist(fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME,'Initialization',CFG.config_filename), 'file') || ...
            exist(fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME,'Results','Eigenvalues.txt'), 'file'))
        f = questdlg("Overwrite existing data in " + fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME) + "?",...
            'Overwrite existing directory?','Yes','No','No');
        if strcmp(f,'No')
            error("Aborted - Directory " + fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME) + " already exists.");
        end
    end
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    if ~isempty(CFG.OUT.RUN_NAME)
        mkdir(fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME));
    end
    mkdir(fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME),'Initialization');
    mkdir(fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME),'Results');

    CFG.OUT.log_filename = fullfile(CFG.pthstr,'Outputs', CFG.OUT.RUN_NAME, 'Results', valOpt(Opts, 'log_filename', 'Logfile.txt'));
    CFG.OUT.eig_filename = fullfile(CFG.pthstr,'Outputs', CFG.OUT.RUN_NAME, 'Results', valOpt(Opts, 'eig_filename', 'Eigenvalues.txt'));
    CFG.GEOM.geom_filename = fullfile(CFG.pthstr,'Inputs',valOpt(Opts, 'geom_filename', 'Geometry.txt'));

    copyfile(CFG.config_filename, fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME,'Initialization'));
    try
        copyfile(CFG.GEOM.geom_filename, fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME,'Initialization'));
    catch
        error("Geometry file "+CFG.GEOM.geom_filename+" not found");
    end

    %% First line of logfile output

    if CFG.OUT.LOG_OUT
        write_log(strcat("OSCILOS_brass version ",OSCILOS_BRASS_VERSION), CFG.OUT.log_filename, 'w');
    end

    %% Plot options

    CFG.PLT.GEOM_PLOT       = valOpt(Opts, 'geom_plot', true);    % Enables Geometry plot
    CFG.PLT.MF_PLOT         = valOpt(Opts, 'mf_plot', true);      % Enables Mean Flow plot
    CFG.PLT.BC_PLOT         = valOpt(Opts, 'bc_plot', true);      % Enables Boundary Condition plot
    CFG.PLT.EIGMAP_PLOT     = valOpt(Opts, 'contour_plot', true); % Enables Eigenvalue contour plot
    CFG.PLT.EFP_PLOT        = valOpt(Opts, 'efp_plot', true);     % Enables Equivalent Fundamental Pitch plot
    CFG.PLT.PLOT_MODES      = valOpt(Opts, 'plot_modes', 5);      % Number of modes to plot
    CFG.PLT.SAVE_FIGS       = valOpt(Opts, 'save_figs', true);    % Specify whether to save .fig files (speeds up output)   
    CFG.PLT.SAVE_PDFS       = valOpt(Opts, 'save_pdfs', true);    % Specify whether to save PDF copies (speeds up output)        
    CFG.PLT.SMALL_PLOTS     = valOpt(Opts, 'small_plots', false); % Small figures for smaller screens
    
    if CFG.OUT.NO_POPUPS % Disable all plots and wait bars
        CFG.PLT.GEOM_PLOT   = false;
        CFG.PLT.MF_PLOT     = false;
        CFG.PLT.BC_PLOT     = false;
        CFG.PLT.EIGMAP_PLOT = false;
        CFG.PLT.EFP_PLOT    = false;
        CFG.PLT.PLOT_MODES  = 0;
        CFG.OUT.WAIT_BARS   = false;
        CFG.OUT.FIN_MSG     = false;
    end
    
    %% Mean flow options

    CFG.MF.P1          = valOpt(Opts, 'p1');           % Mean pressure at the inlet        
    CFG.MF.T1          = valOpt(Opts, 'T1');           % Mean temperature at the inlet           
    CFG.MF.index_M1_u1 = valOpt(Opts, 'choice_M1_u1', 1); % Choice between Mach number (value: 1) or mean flow velocity (value: 2) at the inlet
    CFG.MF.M1_u1       = valOpt(Opts, 'M1_u1');        % Mach number or mean flow velocity at the inlet
    CFG.MF.index_gamma = valOpt(Opts, 'choice_gamma', 1); % Choice between a constant adiabatic index (value: 1) or changing with temperature (value: 2)

    %% Boundary Conditions

    % Inlet
    CFG.BC.Inlet_type    = valOpt(Opts, 'inlet_type');     % Inlet type              
    CFG.BC.Inlet_param1  = valOpt(Opts, 'inlet_param1', 0);   % Amplitude of the reflection coefficient         
    CFG.BC.Inlet_param2  = valOpt(Opts, 'inlet_param2', 0);   % Time lag of the reflection coefficient (ms)
    CFG.BC.Inlet_param3  = valOpt(Opts, 'inlet_param3', 0);   % Additional parameter

    % Outlet
    CFG.BC.Outlet_type    = valOpt(Opts, 'outlet_type');   % Outlet type              
    CFG.BC.Outlet_param1  = valOpt(Opts, 'outlet_param1', 0); % Amplitude of the reflection coefficient         
    CFG.BC.Outlet_param2  = valOpt(Opts, 'outlet_param2', 0); % Time lag of the reflection coefficient (ms)
    CFG.BC.Outlet_param3  = valOpt(Opts, 'outlet_param3', 0); % Additional parameter

    %% Scan range options

    CFG.SR.FreqMin = valOpt(Opts, 'min_freq', 0); % Minimum frequency, in Hertz              
    CFG.SR.FreqMax = valOpt(Opts, 'max_freq', 1000);  % Maximum frequency, in Hertz 
    CFG.SR.FreqNum = valOpt(Opts, 'number_freq', 20);  % Number of frequencies

    CFG.SR.GRMin   = valOpt(Opts, 'min_gr', -200); % Minimum growth rate, in Hertz              
    CFG.SR.GRMax   = valOpt(Opts, 'max_gr', 200);  % Maximum growth rate, in Hertz 
    CFG.SR.GRNum   = valOpt(Opts, 'number_gr', 10);  % Number of growth rates

    if CFG.OUT.CL_OUT; fprintf("done\n"); end
end

%% Function to validate options

function [r] = valOpt(Optcell, Optname, default)
    try % try to find option with Optname in Optcell
        ty = Optcell{2,strcmpi(Optcell(1,:), Optname)};
        st = Optcell{3,strcmpi(Optcell(1,:), Optname)};
    catch % if cannot find Optname in Optcell, throw wobbly
        error("Optname "+Optname+" is not recognised - check definitions in Config_subfc");
    end
    st = char(string(st)); % Handles arguments passed as nonstrings
    if isempty(st) % if option is missing
        if nargin < 3 % if no default given, option is mandatory
            error("Required option "+Optname+" not found while parsing config file");
        else
            r = default; % set missing options to their defaults
        end
    else
        switch ty
            case 'int' % integer
                if isempty(str2double(st)) || isnan(str2double(st)) || (str2double(st)-floor(str2double(st))) > 0
                    warning("Value of option "+Optname+" found in config file ("+st+") cannot be parsed as integer"+newline...
                        +"Default value ("+num2str(default)+") will be used");
                    r = default;
                else
                    r = floor(str2double(st));
                end
            case 'dbl' % double
                if isempty(str2double(st)) || isnan(str2double(st))
                    warning("Value of option "+Optname+" found in config file ("+st+") cannot be parsed as double"+newline...
                        +"Default value ("+num2str(default)+") will be used");
                    r = default;
                else
                    r = str2double(st);
                end
            case 'lgc' % logical
                t = strcmpi('true',st) || strcmpi('1',st) || strcmpi('yes',st) || strcmpi('on',st);
                f = strcmpi('false',st) || strcmpi('0',st) || strcmpi('no',st) || strcmpi('off',st);
                if ~t && ~f
                    warning("Value of option "+Optname+" found in config file ("+st+") cannot be parsed as logical"+newline...
                        +"Default value ("+num2str(default)+") will be used");
                    r = (default == true);
                else
                    r = t;
                end
            case 'str' % string
                r = st;
        end
    end
end

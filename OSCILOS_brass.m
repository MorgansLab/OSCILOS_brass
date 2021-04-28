%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OSCILOS_brass version 5
%
% This program is an extension of the cold-flow portion of the 'lite'
% version of OSCILOS_long. It does not have a graphical interface.
%
% It is designed to model the acoustic behaviour of brass wind instruments
% like the trumpet and trombone, using the OSCILOS solver.
%
% Please run this main script after the input files have been updated. More
% information about the required inputs is included in the User Guide.
%
% OSCILOS_brass(Name, Value, ...) may also be called as a function from any
% other script, providing the folder containing OSCILOS_brass.m is added to
% the path. Arguments are optional and to be given in Name, Value pairs to
% provide additional configuration parameters which override those in the
% config file.
%
% Last update : 09/09/2020
%
% Authors: R. Gaudron and A. MacLaren
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EIG] = OSCILOS_brass(varargin)
    %% INITIALISATION
    
    if mod(nargin,2) == 1 && ~isempty(varargin{1}) % Throw error if odd number of arguments passed
        error("Odd number of arguments (%d) given - even number required for Name, Value pairs",nargin)
    end
    
    % Add directories for subfunctions to path
    filepath = mfilename('fullpath');
    [pthstr,~,~] = fileparts(filepath);
    addpath(fullfile(pthstr, 'SubFunctions'))
    addpath(fullfile(pthstr, 'SubFunctions', 'SubsubFunctions'))
    

    %% PARAMETER CONFIGURATION
    %
    % This subroutine parses the 'Config.txt' file and initialises the mean
    % flow, boundary condition, scan range and plot config variables, which
    % are stored in the CFG data structure.

    CFG = Config_subfc(pthstr,varargin{:});
    
    %% GEOMETRY
    %
    % This subroutine loads the geometry from the file 'Geometry.txt',
    % plots it and saves it to the .pdf and .fig format. The discretised
    % geometry is returned to the GEOM data structure.

    GEOM = Geometry_subfc(CFG);
    
    %% MEAN FLOW
    %
    % This subroutine solves the mean flow conservation equations and
    % returns the distribution of the mean flow variables throughout the
    % domain to the MF data structure. The mean velocity and mean
    % temperature are then plotted and saved in .pdf and .fig format.

    MF = Mean_flow_subfc(CFG, GEOM);

    %% BOUNDARY CONDITIONS
    %
    % This subroutine initialises the acoustic boundary conditions at the
    % inlet and outlet of the geometry and constructs the function handles
    % stored by the BC structure. The acoustic reflection coefficients at
    % the inlet and outlet are plotted and saved in .pdf and .fig format.

    BC = BC_subfc(CFG, GEOM, MF);

    %% FREQUENCY-DOMAIN SOLVER
    %
    % This subroutine finds the eigenvalues of the system, and the corresponding
    % eigenmodes. The eigenvalues are saved as an output file called
    % 'Eigenvalues.txt' (configurable by the eig_filename config option).
    % A map of the eigenvalues in the range of interest is then plotted and
    % saved to the .pdf and .fig format. The number of mode shapes specified
    % by plot_modes, represented as the modulus of the acoustic pressure and
    % velocity, are also plotted and saved to the .pdf and .fig format.

    EIG = Solver_subfc(CFG, GEOM, MF, BC);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utility to import csv file for bore profile and output geometry file 
% for OSCILOS_brass
%
% Input file should be of format:
%   [Axial position], [Radial dimension]
%
% Units of the data can be converted to the required [m] axial and radial
% using the parameters LengthConv and RadConv respectively
%
% The Geometry file output is of the format required by OSCILOS_brass for
% Geometry.txt, namely tab-separated with the following column headers:
%   x[m]	r[m]	SectionIndex	TubeIndex	TubeSplit
%
% A nonzero TubeIndex parameter is used to instruct OSCILOS to interpolate
% between data points - for TubeIndex = 1, the interpolation is linear. 
% A maximum axial distance between output data points may be specified by
% MaxAxGap, and a maximum adjacent section radius ratio (unitary) by
% MaxFracStep. Where the input data contains gaps exceeding these, the 
% number of interpolation points required to smooth the profile to the
% required tolerance(s) is calculated and written to TubeSplit parameter in
% the output file
% 
% May be run as-is with no arguments using default values and UI filepicker
% or externally as a function with up to 7 arguments in the order specified
%
% Author: Alexander MacLaren
%
% Last modified: 06/08/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Geom_file_utility(inflnm, outflnm, MaxFracStep, MinFracStep, MaxAxGap, LengthConv, RadConv)

    %% Parameters
    
    if nargin<1
        [fin,fip] = uigetfile({'*.csv','Comma-Separated Values (*.csv)'; '*.txt','Text Files (*.txt)'; '*.*', 'All Files (*.*)'}, 'Select Data File');
        inflnm = strcat(fip,fin);
    end
    
    if nargin<2
        [fon, fop] = uiputfile('./Inputs/Geometry.txt');
        outflnm = strcat(fop,fon);
    end
    
    % Maximum radius ratio between adjacent sections e.g. 0.1 for 10% increase/decrease (zero to disable, overrides MinFracStep if smaller) refines discretisation
    if nargin < 3; MaxFracStep = 0.02; end
    
    % Minimum radius ratio between adjacent sections e.g. 0.01 for 1% increase/decrease (zero to disable) coarsens discretisation
    if nargin < 4; MinFracStep = 0.01; end
    
    % Maximum axial section length in m (zero to disable)
    if nargin < 5; MaxAxGap = 0; end
    
    % Unit conversion to axial length in m, output value / input value e.g. 0.001 for mm input
    if nargin < 6; LengthConv = 0.001; end
    
    % Unit conversion to radius in m, output value / input value e.g. 0.0005 for mm diameter input
    if nargin < 7; RadConv = 0.0005; end
    
    %% Read csv file
    A = readmatrix(inflnm);

    %% Convert dimensions
    A(:,1) = A(:,1) .* LengthConv;
    A(:,2) = A(:,2) .* RadConv;
    
    %% Sort by x-coordinate
    [A(:,1), xord] = sort(A(:,1));
    A(:,2) = A(xord,2);

    %% Add SectionIndex and TubeIndex information
    A(:,3:5)=zeros(size(A,1),3);

    if MinFracStep > 0 && MinFracStep < MaxFracStep
        B = [];
        j = 1;
        for i=2:size(A,1)
            if A(i,2)/A(j,2) > MinFracStep + 1 || A(j,2)/A(i,2) > MinFracStep + 1
                B = [B,j];
                j=i;
            end
        end
        B = [B,size(A,1)];
        A = A(B,:); % Remove unnecessary elements
    end

    for i=1:size(A,1)-1
        if MaxAxGap > 0
            if A(i+1,1)-A(i,1) > MaxAxGap
                A(i,4) = 1; % Set tube to interpolated cone
                A(i,5) = ceil( (A(i+1,1)-A(i,1)) / MaxAxGap); % Num sections required
            end
        end
        if MaxFracStep > 0
            if A(i+1,2)/A(i,2) > MaxFracStep + 1 || A(i,2)/A(i+1,2) > MaxFracStep + 1
                A(i,4) = 1; % Set tube to interpolated cone
                A(i,5) = ceil( abs(A(i+1,2)-A(i,2)) / (A(i,2) .* MaxFracStep)); % Num sections required
            end
        end
    end

    %% Write geometry file
    H = ["x[m]","r[m]","SectionIndex","TubeIndex","TubeSplit"];
    writematrix(H,outflnm,'Delimiter','\t');
    writematrix(A,outflnm,'Delimiter','\t','WriteMode','append'); % MATLAB R2020a
    %dlmwrite(outflnm,A,'-append','Delimiter','\t'); % MATLAB R2019a

end


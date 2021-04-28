%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Equivalent Fundamental Pitch calculator subsubfunction for OSCILOS_brass
%
%
% Author: Alexander MacLaren
%
% Last Modified: 06/08/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate EFP of vector of frequencies based on 4th harmonic or ref. freq
% given in optional argument

function [efp, F] = EFP_calc(f, F)
    if nargin < 2
        if length(f)<4
            error("efp vector less than 4 in length and no reference frequency provided (4th peak used by default)")
        end
        F = f(4) / 4; % Use 4th frequency peak as reference if none provided
    end
    for i = 1:length(f)
        efp(i) = 1200/log(2) .* log( f(i) ./ (i*F) ); % EFP in cents
    end
end
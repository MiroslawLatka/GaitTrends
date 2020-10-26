function [H] = variogram(x,p)
% Computes  madogram  estimator  of  the  scaling  exponents.   
% Implementation  based  on  the algorithm described by
% Gneiting, T., Ševèíková, H., & Percival, D. B. (2012). Estimators 
% of fractal dimension: Assessing the roughness of time series 
% and spatial data. Statistical Science, 247-277.
% 
% Inputs:
% x: scalar time series (vector)
% p: variogram's order
%
% 
% Outputs:
% H: scaling exponent estimate

% GaitTrends: 
% Authors: Klaudia Kozlowska (Klaudia.Kozlowska@pwr.edu.pl)
%          Miroslaw Latka    (Miroslaw.Latka@pwr.edu.pl)
% URL: https://github.com/mlatka/GaitTrends.git
%
% Copyright (C) 2020  Klaudia Kozlowska and Miroslaw Latka
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% <http://www.gnu.org/licenses/>.

% =========================================================================
% Last update: July 21, 2020
% =========================================================================

% Citing the GaitTrends:
% https://doi.org/10.1371/journal.pcbi.1007180


% =========================================================================


n = length(x);


V1 = 0;
V2 = 0;

for i = 1 : n-1
    V1 = V1 + (abs(x(i+1)-x(i)))^p;
end

for i = 1 : n-2
    V2 = V2 + (abs(x(i+2)-x(i)))^p;
end

V1 = 0.5*V1/(n-1);
V2 = 0.5*V2/(n-2);

H = (1./p)* (log(V2)-log(V1))/log(2);

end


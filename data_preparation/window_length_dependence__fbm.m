% =========================================================================
%  This script generates an ensemble of 500 fractional Brownian motion (fbm) 
%  random walks with a chosen Hurst exponent. The length of time series 
%  was set to 260. Each trajectory was divided into non-overlapping windows
%  of length k from k=40 to k=260 with step 20. For each window, the madogram
%  estimator  and  detrended fluctuation analysis of order n=1 to n=3  were 
%  used to compute the scaling exponents. The boxplots of scaling exponents 
%  for all four methods are plotted as a function of window length.
%  For a computer with 2 cores it may take up to 15 minutes to complete the
%  calculations! The execution time decreases with the number of workers in
%  the MATLAB parallel pool.
%
% Before running the script, please set the value of Hurst exponent (alpha).
% By default the variable saveResults is set to true  so that the scaling
% exponents are saved to a MAT-file.
% Please ensure that you have added WFDB Toolbox folder (in libs/ folder) 
% to the MATLAB search path. dfa function from WFDB library is used to 
% perform detrended fluctuation analysis.
% =========================================================================
%
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
% https://doi.org/10.1101/677948

% =========================================================================


clc, clear, close all

alpha = 0.75; %Hurst exponent
saveResults = true;

outputDir = '../data/window_length_dependence/';
window_lengths = 40:20:260;
alpha_matrix = [];

parfor i = 1 : 500
    
    % generate fbm time series
    series = diff(wfbm(alpha,260));
    
    for w = 1 : length(window_lengths)
        
        n = floor(length(series)/window_lengths(w));
        chunks = reshape(series(1:n*window_lengths(w)), [], n);

        sc = size(chunks);
        
        % calculate exponents for chunks
        for c = 1 : sc(2)
             
             Hrow = [];
             fragment = chunks(:,c);
             
            if(length(fragment) > 1)
                 
                % DFA1
                [ln,lf] = dfa(fragment,1);
                fit = polyfit(ln,lf,1);
                H1 = fit(1);
                % DFA2
                [ln,lf] = dfa(fragment,2);
                fit = polyfit(ln,lf,1);
                H2 = fit(1);
                % DFA3
                [ln,lf] = dfa(fragment,3);
                fit = polyfit(ln,lf,1);
                H3 = fit(1);  
                % MAD
                H_MD = variogram(cumsum(fragment),1);
                 
                row =  [H1 H2 H3 H_MD window_lengths(w) i];
                alpha_matrix = [alpha_matrix; row];
  
            end % end if
             
         end % end for c (chunks)
        
        
    end % end for w (windows)
    
    
end % end for i (subjects)

% save results (optional)
if(saveResults)
    file = strcat(outputDir,'fbm_',num2str(alpha),'_alphas_windowed.mat');
    save(file,'alpha_matrix');
    disp(strcat('Data saved to: ',file));
end

% visualize results
labels = {'{\alpha}^(^1^)','{\alpha}^(^2^)','{\alpha}^(^3^)', ...
    '{\alpha}^(^M^D^)'};

window_lengths = 40:20:260;

for col = 1:4 % 1 - alpha1, 2 - alpha2, 3 - alpha3

    dataBoxplot = [];
    groups = [];
    
    for w = 1 : length(window_lengths)

        wl = window_lengths(w);
        ind = find(alpha_matrix(:,5) == wl);
        selectedWindow = alpha_matrix(ind,col);
        dataBoxplot = [dataBoxplot; selectedWindow];
        groups = [groups; wl*ones(size(selectedWindow))];

    end
    
    figure;
    boxplot(dataBoxplot,groups);
    title(labels{col});
    xlabel('window length');
    ylabel('scaling exponent');

end



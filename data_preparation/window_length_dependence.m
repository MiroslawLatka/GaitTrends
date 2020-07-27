% =========================================================================
% This script divides first 260 samples of ST and SL experimental time series
% from ../data/mat_data folder into non-overlapping windows of length k 
%(from k = 40 to k = 260 with step 20) and then for each such window 
% calculates the scaling exponents using DFA1-3 and madogram algorithms.
% The corresponding exponents are stored in the first four columns of 
% the matrix. The window length and the trial’s id may be found 
% in columns 5 and 6, respectively. Set fileName variable to choose 
% the desired parameter (SL  or ST) and the  speed. E.g. fileName = Ln_SPD1.mat
% corresponds to stride length at preferred walking speed while Tn_SPD1.mat}
% is the stride time series for the same speed. The output MAT files are saved
% in ../data/window_length_dependence folder. For a given scaling exponents,
% the script generates the boxplots of its values for all window lengths.
% The script uses dfa function from WFDB library. Please ensure 
% that you have added WFDB Toolbox folder (in libs/ folder) 
% to the MATLAB search path. dfa  function from WFDB library is used to
% perform detrended fluctuation analysis.
% =========================================================================

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

% 1 - 100, 2 - 110, 3 - 90, 4 - 120, 5 - 80 [%PWS]
fileName = 'Ln_SPD1.mat';
saveResults = true;

addpath('../data/mat_data');
addpath('../utils/');

data = load(fileName);
s = size(data.residualsAll);
alpha_matrix = [];
window_lengths = 40:20:260;

tic
parfor i = 1 : s(2)
    
    % gait parameter time series
    series = data.seriesAll{i};
    
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
                H_MD = variogram(cumsum(fragment-mean(fragment)),1);
                 
                row =  [H1 H2 H3 H_MD window_lengths(w) i];
                alpha_matrix = [alpha_matrix; row];
  
             end % end if
             
        end % end for c (chunks)
        
        
    end % end for w (windows)
    
    
end % end for i (subjects)
toc

% save results (optional)
if(saveResults)
    outputDir = '../data/window_length_dependence/';
    file = strcat(outputDir,'dingwell_alphas_windowed_',fileName);
    save(file,'alpha_matrix');
    disp(strcat('Data saved to: ',file));
end

% visualize results
labels = {'{\alpha}^(^1^)','{\alpha}^(^2^)','{\alpha}^(^3^)', ...
    '{\alpha}^(^M^D^)'};

window_lengths = 40:20:260;

for col = 1 : 4 % 1 - alpha1, 2 - alpha2, 3 - alpha3

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

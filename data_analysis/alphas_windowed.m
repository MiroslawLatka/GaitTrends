% =========================================================================
% This script visualizes the dependence of scaling exponents of Dingwell’s 
% gait time series  and  fractional Brownian motion on data window length. 
% The boxplots are displayed for the DFA and madogram scaling exponents. 
% Change fileName variable to analyze either experimental data or fractional
% Brownian motion time series. Both types of time series are stored in
% ../data/window_length_dependence folder. For example, with
% fileName='dingwell_alphas_windowed_Ln_SPD1.mat' the calculations are 
% performed for SL at 100%PWS whereas for fileName='fbm_0.4_alphas_windowed.mat'
% for fractional Brownian motion with the Hurst exponent equal to 0.4.
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
% SPD: 1 - 100, 2 - 110, 3 - 90, 4 - 120, 5 - 80 [%PWS]
fileName = 'dingwell_alphas_windowed_Ln_SPD1.mat';
% fileName = 'fbm_0.4_alphas_windowed.mat';
% fileName = 'fbm_0.75_alphas_windowed.mat';

addpath('../data/window_length_dependence');
labels = {'{\alpha}^(^1^)','{\alpha}^(^2^)','{\alpha}^(^3^)', ...
    '{\alpha}^(^M^D^)'};
dataBoxplot = [];
groups = [];
% load data
data = load(fileName);

for col = 1 : 4 % 1 - alpha1, 2 - alpha2, 3 - alpha3, 4 - alphaMD

	window_lengths = 40:20:260;

	for w = 1 : length(window_lengths)
		
		wl = window_lengths(w);
		ind = find(data.alpha_matrix(:,5) == wl);
		selectedWindow = data.alpha_matrix(ind,col);
		dataBoxplot = [dataBoxplot; selectedWindow];
		groups = [groups; wl*ones(size(selectedWindow))];

	end

	figure;
	boxplot(dataBoxplot,groups);
	title(labels{col});
	xlabel('window length');
	ylabel('scaling exponent');

end
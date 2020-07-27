% =========================================================================
% This script  uses DFA1-3 and madogram algorithms to calculate  the scaling
% indices  of Dingwell's experimental time series  and those of the corresponding 
% MARS residuals. The input files must be generated first using 
% ../data_preparation/prepare_data.m script. For a given treadmill speed, 
% the script also determines whether the scaling is anti-persistent (i.e. 
% the median or mean of scaling exponents is statistically smaller than 0.5). 
% Set attributeNumber variable  to 1, 2, and 3 to  select  SL, ST,
% and SS, respectively.
%
% Please ensure that you have added WFDB Toolbox folder 
% (in libs/ folder) to the MATLAB search path. dfa function from WFDB library
% is used to perform detrended fluctuation analysis.
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
attributeNumber = 1; % 1 - SL, 2 - ST, 3 - SS


addpath('../data/mat_data/');
addpath('../utils/');

switch attributeNumber
    case 1
        fileNamesCell = {'Ln_SPD1.mat','Ln_SPD2.mat',...
        'Ln_SPD3.mat','Ln_SPD4.mat',...
        'Ln_SPD5.mat'};
    case 2
        fileNamesCell = {'Tn_SPD1.mat','Tn_SPD2.mat',...
        'Tn_SPD3.mat','Tn_SPD4.mat',...
        'Tn_SPD5.mat'};
    case 3
        fileNamesCell = {'Sn_SPD1.mat','Sn_SPD2.mat',...
        'Sn_SPD3.mat','Sn_SPD4.mat',...
        'Sn_SPD5.mat'};
    otherwise
        error('Error. attributeNumber must be an integer between 1 and 3.')
end

s = size(fileNamesCell);
alphas_cell = {};

for i = 1 : s(2)
    
	% load data
    data = load(fileNamesCell{i});
	trialSize = size(data.residualsAll);
	alphas_matrix = [];

	for j = 1 : trialSize(2)
		
		% get series and residuals
        series = data.seriesAll{j};
		residuals = data.residualsAll{j};

		row = [];

		% calculate exponents for original series
		% DFA1
		[ln,lf] = dfa(series,1);
		fit = polyfit(ln,lf,1);
		H1o = fit(1);
		% DFA2
		[ln,lf] = dfa(series,2);
		fit = polyfit(ln,lf,1);
		H2o = fit(1);
		% DFA3
		[ln,lf] = dfa(series,3);
		fit = polyfit(ln,lf,1);
		H3o = fit(1);
		% madogram
		H_MDo = variogram(cumsum(series-mean(series)),1);
		 
		% calculate exponents for residuals series
		% DFA1
		[ln,lf] = dfa(residuals,1);
		fit = polyfit(ln,lf,1);
		H1d = fit(1);
		% DFA2
		[ln,lf] = dfa(residuals,2);
		fit = polyfit(ln,lf,1);
		H2d = fit(1);
		% DFA3
		[ln,lf] = dfa(residuals,3);
		fit = polyfit(ln,lf,1);
		H3d = fit(1);
		% madogram
		H_MDd = variogram(cumsum(residuals-mean(residuals)),1);

		row =  [H1o H2o H3o H_MDo H1d H2d H3d H_MDd];
			
		alphas_matrix = [alphas_matrix; row]; 

	end % end trial loop

	alphas_cell{end+1} = alphas_matrix;

end % end speed loop

% checking if exponents mean/median < 0.5
clc
% L denotes detrended series
alphas_labels = {'alpha1','alpha2','alpha3','alphaMD',...
        'alphaL1','alphaL2','alphaL3','alphaLMD'};
    
for SPD = 1:5

    % alphas matrix
    alphas_mat = alphas_cell{SPD};
    disp(strcat('SPD',num2str(SPD)))

    d = size(alphas_mat);
    pVec = [];

    for k = 1 : d(2)

        v = alphas_mat(:,k);

        % Shapiro-Wilk test for normality
        [hs, ps] = swtest(v);

        if hs == 0 % t-test
            [h,p] = ttest(v,0.5,'tail','left');
        else % Wilcoxon signed rank test
            [p,h] = signrank(v,0.5,'tail','left');
        end

        disp(strcat(alphas_labels{k},': H=',num2str(h),', p=',num2str(p)))

        pVec = [pVec; p];

    end

end

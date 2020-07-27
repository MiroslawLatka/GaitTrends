% =========================================================================
% This script  uses DFA1-3 and madogram algorithms to calculate  the scaling
% indices  of surrogates of Dingwell's experimental time series (SL,ST,SS)  
% and those of the corresponding  MARS residuals.  The input files must be 
% generated first using ../data_preparation/prepare_surrogates.m script.
% For a given treadmill speed the scripts also determines whether the scaling
% is anti-persistent (i.e.  the median or mean of scaling exponents is 
% statistically smaller than 0.5). Before running the script,  please set
% the type of surrogates (1 - independent, 2 - cross-correlated, 3 - shuffled). 
% The boxplots of madogram scaling exponents of gait surrogates are displayed.
%
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

surrogates_type = 1; % 1 - independent, 2 - cross-correlated, 3 - shuffled

addpath('../utils/');

if(surrogates_type < 1 || surrogates_type > 3)
    error('Error. surrogates_type must be an integer between 1 and 3.')
end

alphasST_cell = {};
alphasSL_cell = {};
alphasSS_cell = {};

for i = 1 : 5 % 1-5 SPD
    
	switch surrogates_type
		case 1
			dir = '../data/surrogates/independent/'; 
			data = load(strcat(dir,'independent_surrogates_SPD',num2str(i),'.mat'));
			ver = 'independent';
		case 2
			dir = '../data/surrogates/cross_correlated/';
			data = load(strcat(dir,'cross_correlated_surrogates_SPD',num2str(i),'.mat'));
			ver = 'cross-correlated';
		case 3
			dir = '../data/surrogates/shuffled/';
			data = load(strcat(dir,'shuffled_surrogates_SPD',num2str(i),'.mat'));
			ver = 'shuffled';
		otherwise
			error('Error. Attribute must be value between 1 and 3.')
	end

	trialSize = size(data.data_surrogatesSL.seriesAll);
	alphas_ST = [];
	alphas_SL = [];
	alphas_SS = [];

	for j = 1 : trialSize(2)
		
		% load SL, ST, SS
        sl = data.data_surrogatesSL.seriesAll{j};
		st = data.data_surrogatesST.seriesAll{j};
		ss = sl./st;
			
		% calculate scaling exponent using madogram
		H_MDsl = variogram(cumsum(sl-mean(sl)),1);
		H_MDst = variogram(cumsum(st-mean(st)),1);
		H_MDss = variogram(cumsum(ss-mean(ss)),1);


		alphas_SL = [alphas_SL; H_MDsl];
		alphas_ST = [alphas_ST; H_MDst];
		alphas_SS = [alphas_SS; H_MDss];

	end % end trial loop

alphasST_cell{end+1} = alphas_ST;
alphasSL_cell{end+1} = alphas_SL;
alphasSS_cell{end+1} = alphas_SS;

end % end speed loop

% visualize results

% SL madogram
boxplots_for_all_speeds(alphasSL_cell, strcat(ver,' surrogates'),...
    '{\alpha}^(^M^D^ for SL');


% ST madogram
boxplots_for_all_speeds(alphasST_cell, strcat(ver,' surrogates'),...
    '{\alpha}^(^M^D^ for ST');

% SS madogram
boxplots_for_all_speeds(alphasSS_cell, strcat(ver,' surrogates'),...
    '{\alpha}^(^M^D^ for SS');


% checking if exponents mean/median < 0.5
clc
disp('madogram exponents:')
check_alphas(alphasSL_cell, 'SL')
check_alphas(alphasST_cell, 'ST')
check_alphas(alphasSS_cell, 'SS')

function check_alphas(data_cell, param_name)
% Determines whether scaling exponents mean or median is lower than 0.5.
% 
% Inputs:
% data_cell:  cell with five vectors (for each treadmill speed)
% param_name: parameter name ('SL'/'ST'/'SS')
%
% Outputs:    none

    for SPD = 1:5

        % alphas matrix
        alphas_vec = data_cell{SPD};

        % Shapiro-Wilk test for normality
        [hs, ps] = swtest(alphas_vec);

        if hs == 0 % t-test
            [h,p] = ttest(alphas_vec,0.5,'tail','left');
        else % Wilcoxon signed rank test
            [p,h] = signrank(alphas_vec,0.5,'tail','left');
        end

        disp(strcat('SPD',num2str(SPD),':',param_name,' H=',num2str(h),', p=',num2str(p)))

    end

end




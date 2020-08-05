% =========================================================================
% Two-sample Kolmogorov-Smirnov test is used to compare the
% distributions of trend durations and slopes in experimental data and
% their phase-randomized surrogates. The script uses MAT-files located in 
% ../data/surrogates/surrogates_trend_distributions/ folder. These input
% files must be created first by running ../data_preparation/multiple_surrogates.m
% Before running the script,  please choose:
%       1) type of surrogates (cross-correlated = true, independent = false)
%       2) whether you want to analyze trend durations (param=1) or trend
%          slopes (param=2).
% The value of repetitions variable  may not be greater than 
% the number of surrogates created by the multiple_surrogates.m script.
% The number of cases in which distributions are not statistically
% different is printed in the MATLAB's command window.
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
clc, clear, close all

cross_correlated = false;
param = 1; % 1 - trend durations, 2 - trend slopes

% make sure that the value of repetitions variable  is not greater than 
% the number of surrogates created by the multiple_surrogates.m script.
repetitions = 50;


addpath('../data/surrogates/surrogates_trend_distributions/');
addpath('../data/trend_stats/');

if(~(param == 1 || param == 2))
    error('Error. param must be 1 or 2.')
end

all_trend_statsSL = cell(repetitions,0);
all_trend_statsST = cell(repetitions,0);

for x = 1 : repetitions
    all_trend_statsSL{x} = [];
    all_trend_statsST{x} = [];
end

if(cross_correlated)  
        fileNamesCell = {'distr_cross_correlated_surrogates_SPD1.mat',...
        'distr_cross_correlated_surrogates_SPD2.mat',...
        'distr_cross_correlated_surrogates_SPD3.mat',...
        'distr_cross_correlated_surrogates_SPD4.mat',...
        'distr_cross_correlated_surrogates_SPD5.mat'};  
else
    fileNamesCell = {'distr_independent_surrogates_SPD1.mat',...
        'distr_independent_surrogates_SPD2.mat',...
        'distr_independent_surrogates_SPD3.mat',...
        'distr_independent_surrogates_SPD4.mat',...
        'distr_independent_surrogates_SPD5.mat'};    
end

% aggregate surrogate data
s = size(fileNamesCell);
for i = 1 : s(2)  
    
    data = load(fileNamesCell{i}); 

	for r = 1 : repetitions	
        
        if(param == 1)
            all_trend_statsSL{r} = [all_trend_statsSL{r};  ...
                data.trend_durations_statsSL{r}];
            all_trend_statsST{r} = [all_trend_statsST{r}; ...
                data.trend_durations_statsST{r}]; 
        else
            all_trend_statsSL{r} = [all_trend_statsSL{r};  ...
                data.trend_slopes_statsSL{r}];
            all_trend_statsST{r} = [all_trend_statsST{r}; ...
                data.trend_slopes_statsST{r}]; 
        end
    end
	
end

% load experimental data
SL_trends_Org = load('SL_trends.mat');
ST_trends_Org = load('ST_trends.mat');

orgSL = [];
orgST = [];

% aggregate experimental data
for l = 1 : length(SL_trends_Org.trend_durations_cell)
    if(param == 1) 
        orgSL = [orgSL; SL_trends_Org.trend_durations_cell{l}];
        orgST = [orgST; ST_trends_Org.trend_durations_cell{l}];
    else
        orgSL = [orgSL; SL_trends_Org.trend_slopes_cell{l}];
        orgST = [orgST; ST_trends_Org.trend_slopes_cell{l}];
    end
end

countSL = 0;
countST = 0;
meansSL = zeros(1,repetitions);
meansST = zeros(1,repetitions);

% compare distributions
for r = 1 : repetitions
    
	% H0 - same distributions
	[hSL,pSL,~] = kstest2(orgSL, all_trend_statsSL{r});
	
	if(hSL == 0)
        countSL = countSL + 1;	
    end
    meansSL(r) = mean(all_trend_statsSL{r});
	
	% H0 - same distributions
	[hST,pST,~] = kstest2(orgST, all_trend_statsST{r});
	
	if(hST == 0)
		countST = countST + 1;	
    end
    meansST(r) = mean(all_trend_statsST{r});
	
    % visualize results   
    figure;
    subplot(1,2,1)
    histogram(orgSL,30,'Normalization','pdf');
    title(strcat('experimental SL mean:',{' '},num2str(mean(orgSL))));
    subplot(1,2,2)
    histogram(all_trend_statsSL{r},30,'Normalization','pdf');
    title(strcat('SL surrogate mean:',{' '},num2str(mean(all_trend_statsSL{r}))));
    
    figure;
    subplot(1,2,1)
    histogram(orgST,30,'Normalization','pdf');
    title(strcat('experimental ST mean:',{' '},num2str(mean(orgST))));
    subplot(1,2,2)
    histogram(all_trend_statsST{r},30,'Normalization','pdf');
    title(strcat('ST surrogate mean:',{' '},num2str(mean(all_trend_statsST{r}))));


end

disp('same distribution:')
disp(strcat('SL:',{' '},num2str(countSL),'/',num2str(repetitions),...
    {' '},'mean',{' '},num2str(mean(meansSL)),'+/-',...
    num2str(std(meansSL)),{' '},'range',{' '},num2str(min(meansSL)),'-',...
    num2str(max(meansSL)))) 
disp(strcat('ST:',{' '},num2str(countST),'/',num2str(repetitions),...
    {' '},'mean',{' '},num2str(mean(meansST)),'+/-',...
    num2str(std(meansST)),{' '},'range',{' '},num2str(min(meansST)),'-',...
    num2str(max(meansST)))) 

% =========================================================================
% This script loads surrogate data created with ../data_preparation/prepare_surrogates.m
% script and calculates Pearson's correlation coefficients between
% ST and SL surrogates as well as their MARS trends and the corresponding 
% MARS residuals. For trends, the boxplots of Pearson's correlation coefficients 
% are plotted for each treadmill speed. Before running the script, 
% please set type of surrogates (cross-correlated=true, independent = false).
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

cross_correlated = true;

addpath('../data/surrogates/');
corr_cell = {};

for i = 1 : 5
    
	corrMatrix = [];

    % load data for given surrogates
	if(cross_correlated)
		dir = '../data/surrogates/cross_correlated/';
		data = load(strcat(dir,'cross_correlated_surrogates_SPD',num2str(i),'.mat'));
		ver = 'cross-correlated';
    else
        dir = '../data/surrogates/independent/'; 
		data = load(strcat(dir,'independent_surrogates_SPD',num2str(i),'.mat'));
		ver = 'independent';
	end

	% loading SL data
	SLdata = data.data_surrogatesSL;
	% loading ST data
	STdata = data.data_surrogatesST;

	s = size(SLdata.residualsAll);

	for j = 1 : s(2)
		
		corrRow = [];

		% original series
		orgST = STdata.seriesAll{j};
		orgSL = SLdata.seriesAll{j};
		
        % correlation coefficient and p-value for original series
		[rhoOrg, pOrg] = corr(orgST,orgSL,'type','Pearson');
			
		% MARS trends
		trST = STdata.trendsAll{j};
		trSL = SLdata.trendsAll{j};
		
        % correlation coefficient and p-value for MARS trend
		[rhoTr, pTr] = corr(trST,trSL,'type','Pearson');
		 
		% MARS residuals
		resST = STdata.residualsAll{j};
		resSL = SLdata.residualsAll{j};
		
         % correlation coefficient and p-value for MARS residuals
		[rhoRes, pRes] = corr(resST,resSL,'type','Pearson');  
		 
		corrRow =  [rhoOrg pOrg rhoTr pTr rhoRes pRes];     
		corrMatrix = [corrMatrix; corrRow]; 
		
	end % end trial loop

	corr_cell{end+1} = corrMatrix;

end % end speed loop

% visualize results for trends
ind = [5, 3, 1, 2, 4];
dat = [];
group = [];

for l = 1 : length(ind)
    dat = [dat; corr_cell{ind(l)}(:,3)];
    group = [group; l*ones(size(corr_cell{ind(l)}(:,3)))];
end
figure;
boxplot(dat,group);
xlabel('treadmill speed [%PWS]');
ylabel('trend correlation coefficient');
title(strcat(ver,' surrogates'));
set(gca,'XTickLabel',{'80','90','100','110','120'});


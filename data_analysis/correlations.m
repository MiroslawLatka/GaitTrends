% =========================================================================
% This script calculates Pearson's correlation coefficient between ST and
% SL for: experimental time series, their piecewise linear MARS trends, 
% and their MARS residuals. First, the input files must be generated using
% ../data_preparation/prepare_data.mscript. For trends, the boxplots 
% of Pearson's correlation coefficients are plotted for each treadmill speed.
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

addpath('../data/mat_data/');

corr_cell = {};

for i = 1 : 5
    
	% load SL data
	SLdata = load(strcat('Ln_SPD',num2str(i),'.mat'));
	% load ST data
	STdata = load(strcat('Tn_SPD',num2str(i),'.mat'));
		
	trialSize = size(SLdata.residualsAll);
	corrMatrix = [];

	for j = 1 : trialSize(2)
		
		corrRow = [];

		% original series
		orgST = STdata.seriesAll{j};
		orgSL = SLdata.seriesAll{j};
		
        % correlation coefficient and p-value for original series
		[rhoOrg, pOrg] = corr(orgST,orgSL,'type','Pearson');
			
		% MARS trens
		trST = STdata.trendsAll{j};
		trSL = SLdata.trendsAll{j};
		
        % correlation coefficient and p-value for MARS trends
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
set(gca,'XTickLabel',{'80','90','100','110','120'});


% =========================================================================
% This script loads gait data from Dingwell’s MAT-files and calculates 
% shuffled surrogate data. Then it finds 
% piecewise linear MARS trends in surrogate time series.
% Before running the script,  please set:
%       1) speed (SPD)
%       2) random number generator seeds.
% By default the variable generateFigures is set to true 
% so that the surrogate time series for ST and SL and their trends 
% are plotted for all subjects. 
% Please ensure that you have added ARESLab folder (in libs/ folder) 
% to the MATLAB search path. 
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

SPD = 1; % 1 - 100, 2 - 110, 3 - 90, 4 - 120, 5 - 80 [%PWS]
seed=456
rng(seed)
generateFigures  = true;

if(SPD < 1 || SPD > 5)
    error('Error. SPD must be an integer between 1 and 5.')
end

addpath('../data/mat_data/');
% load SL data
SLdata = load(strcat('Ln_SPD',num2str(SPD),'.mat'));
% load ST data
STdata = load(strcat('Tn_SPD',num2str(SPD),'.mat'));

s = size(SLdata.residualsAll);
surrogatesSL = {};
surrogatesST = {};
surrogatesSL_trends = {};
surrogatesST_trends = {};
timestampsSL = {};
timestampsST = {};
surrogatesSL_residualsAll = {};
surrogatesST_residualsAll = {};
modelsSL = {};
modelsST = {};
infosSL = {};
infosST = {};
knotIndicesST = {};
knotIndicesSL = {};
valuesAtKnotSL = {};
valuesAtKnotST = {};

% MARS params
params = aresparams2('useMinSpan',-1,'useEndSpan',-1,'cubic',false,...
    'c',2,'threshold',1e-3,'maxFuncs',50);

for i = 1 : s(2)
		
		% raw series
		orgST = STdata.seriesAll{i};
		orgSL = SLdata.seriesAll{i};
		T = STdata.timestampsAll{i};
		
		% shuffled surrogates
        surST =  orgST(randperm(length(orgST)));
		surSL =  orgSL(randperm(length(orgSL)));
		
		surrogatesST{end+1} = surST;
		surrogatesSL{end+1} = surSL;
        
		% MARS for SL       
		X = T(1:length(surSL));
		Y = surSL; 
        % build MARS model using ARESLab
		model = aresbuild(X,Y,params);
        % calculate trends
		predSL = arespredict(model,X);
        % calculate MARS residuals
		detrendedSL = Y - predSL;
		
		surrogatesSL_trends{end+1} = predSL;
        surrogatesSL_residualsAll{end+1} = detrendedSL;
		modelsSL{end+1} = model;
		infosSL{end+1} = SLdata.infosAll{i};
		timestampsSL{end+1} = {X};
		
		% find MARS knots
        knotSites = sort(cell2mat(model.knotsites));
		knotIndices = zeros(1,length(knotSites));
		
		for k = 1 : length(knotSites)
			for j = 1 : length(X)
				if (knotSites(k) == X(j))
					knotIndices(k) = j;
					break;
				 end
			end

		end
		
		valuesAtKnot = [];
		
		if(knotSites > 0)
			knotIndices = [1 knotIndices];
			knotSites = [1; knotSites];
			valuesAtKnot = predSL(knotIndices);
		end
		
		knotIndicesSL{end+1} = knotIndices;
		valuesAtKnotSL{end+1} = valuesAtKnot;
		 
        % MARS for ST 
		X = T(1:length(surST));
		Y = surST; 
        % build MARS model using ARESLab
		model = aresbuild(X,Y,params);
        % calculate trends
		predST = arespredict(model,X);
		% calculate residuals (noise)
		detrendedST = Y - predST;
		
        surrogatesST_trends{end+1} = predST;
		surrogatesST_residualsAll{end+1} = detrendedST;
		modelsST{end+1} = model;
		infosST{end+1} = STdata.infosAll{i};
		timestampsST{end+1} = {X};

		% find MARS knots
        knotSites = sort(cell2mat(model.knotsites));
		knotIndices = zeros(1,length(knotSites));
		
		for k = 1 : length(knotSites)
			for j = 1 : length(X)
				if (knotSites(k) == X(j))
					knotIndices(k) = j;
					break;
				 end
			end

		end
		
		valuesAtKnot = [];
		
		if(knotSites > 0)
			knotIndices = [1 knotIndices];
			knotSites = [1; knotSites];
			valuesAtKnot = predST(knotIndices);
		end
		
		knotIndicesST{end+1} = knotIndices;
		valuesAtKnotST{end+1} = valuesAtKnot;
		
		if(generateFigures)
			figure;
			plot(surSL,'b--'); hold on;
			plot(surST,'r--');
			plot(predSL,'b');
			plot(predST,'r'); hold off;
			title('surrogate data');  		
		end
  
end

% prepare data for saving
data_surrogatesSL.trendsAll = surrogatesSL_trends;
data_surrogatesSL.seriesAll = surrogatesSL;
data_surrogatesSL.timestampsAll = timestampsSL;
data_surrogatesSL.residualsAll = surrogatesSL_residualsAll;
data_surrogatesSL.modelsAll = modelsSL;
data_surrogatesSL.infosAll = infosSL;
data_surrogatesSL.knotIndicesAll = knotIndicesSL;
data_surrogatesSL.valuesAtKnotAll = valuesAtKnotSL;

data_surrogatesST.trendsAll = surrogatesST_trends;
data_surrogatesST.seriesAll = surrogatesST;
data_surrogatesST.timestampsAll = timestampsST;
data_surrogatesST.residualsAll = surrogatesST_residualsAll;
data_surrogatesST.modelsAll = modelsST;
data_surrogatesST.infosAll = infosST;
data_surrogatesST.knotIndicesAll = knotIndicesST;
data_surrogatesST.valuesAtKnotAll = valuesAtKnotST;

% save data
file = strcat('../data/surrogates/shuffled/shuffled_surrogates','_SPD',num2str(SPD),'.mat');
save(file,'data_surrogatesSL','data_surrogatesST');
disp(strcat('Data saved to: ', file));

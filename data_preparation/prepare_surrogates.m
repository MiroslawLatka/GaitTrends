% =========================================================================
% This script loads gait data from Dingwell’s MAT-files and calculates 
% phase-randomized surrogates using TISEAN library. Then it finds 
% piecewise linear MARS trends in surrogate time series.
% Before running the script,  please set:
%       1) speed (SPD)
%       2) type of surrogates (cross-correlated true, independent = false)
%       3) random number generator seeds.
% By default the variable generateFigures is set to true 
% so that the surrogate time series for ST and SL and their trends 
% are plotted for all subjects. 

% This script uses surrogates.exe file from TISEAN library.
% Two text auxiliary files: in.txt and out.txt are used 
% as input and output to this function.
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
cross_correlated = true;
seed = [567 789]; % random number generator seeds
generateFigures = true;

if(SPD < 1 || SPD > 5)
    error('Error. SPD must be an integer between 1 and 5.')
end

currentPath = pwd;
addpath('../libs/');
addpath('../data/mat_data/');
% load SL data
SLdata = load(strcat('SL_SPD',num2str(SPD),'.mat'));
% load ST data
STdata = load(strcat('ST_SPD',num2str(SPD),'.mat'));

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

    % original series
    orgST = STdata.seriesAll{i};
    orgSL = SLdata.seriesAll{i};
    T = STdata.timestampsAll{i};
    
    if(cross_correlated) % correlated surrogates
		str = 'cross_correlated_surrogates';
		
		tot = [orgSL  orgST]; 
        inFile = strcat(pwd, '\in.txt');
        save(inFile,'tot','-ascii');
        cmd = strcat('surrogates.exe -o',{' '},pwd, ...
            '\out.txt -m 2 -I',{' '},num2str(seed(1)),{' '},pwd,'\in.txt');
        cd('../libs/');
        system(cmd{1});
        cd(currentPath)
		surrogatesTot = load('out.txt');
 
        surSL = surrogatesTot(:,1);
		surST = surrogatesTot(:,2);
		
		surrogatesSL{end+1} = surSL;
		surrogatesST{end+1} = surST;
		
    else % independent surrogates
		 
		str = 'independent_surrogates';       
        
		inFile = strcat(pwd, '\in.txt');
        save(inFile,'orgSL','-ascii');
        cmd = strcat('surrogates.exe -o',{' '},pwd, ...
            '\out.txt -I',{' '},num2str(seed(1)),{' '},pwd,'\in.txt');
        cd('../libs/');
        system(cmd{1});
        cd(currentPath)
		surSL = load('out.txt'); 
		surrogatesSL{end+1} = surSL;
		
		inFile = strcat(pwd, '\in.txt');
        save(inFile,'orgST','-ascii');
        cmd = strcat('surrogates.exe -o',{' '},pwd, ...
            '\out.txt -I',{' '},num2str(seed(2)),{' '},pwd,'\in.txt');
        cd('../libs/');
        system(cmd{1});
        cd(currentPath)
		surST = load('out.txt');
		surrogatesST{end+1} = surST;		 
    end
    
    % perform calculations for SL surrogates
     plotTitle = strcat('Y',num2str(i),{'SPD'}, num2str(SPD));
     outputData = perform_MARS(T(1:length(surSL)),surSL,...
        generateFigures,plotTitle,'SL surrogate',params);
    
     surrogatesSL_residualsAll{end+1} = outputData.residuals;
     surrogatesSL_trends{end+1} = outputData.trends;
     timestampsSL{end+1} = outputData.timestamps;
     modelsSL{end+1} = outputData.model;
     infosSL{end+1} = outputData.info;    
     knotIndicesSL{end+1} = outputData.knotIndices;
     valuesAtKnotSL{end+1} = outputData.valuesAtKnot;
     
    % perform calculations for ST surrogates
     plotTitle = strcat('Y',num2str(i),{'SPD'}, num2str(SPD));
     outputData = perform_MARS(T(1:length(surST)),surST,...
        generateFigures,plotTitle,'ST surrogate',params);
    
     surrogatesST_residualsAll{end+1} = outputData.residuals;
     surrogatesST_trends{end+1} = outputData.trends;
     timestampsST{end+1} = outputData.timestamps;
     modelsST{end+1} = outputData.model;
     infosST{end+1} = outputData.info;    
     knotIndicesST{end+1} = outputData.knotIndices;
     valuesAtKnotST{end+1} = outputData.valuesAtKnot;
    
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
if(cross_correlated)
	file = strcat('../data/surrogates/cross_correlated/',str,'_SPD',num2str(SPD),'.mat');
else
	file = strcat('../data/surrogates/independent/',str,'_SPD',num2str(SPD),'.mat');
end

save(file,'data_surrogatesSL','data_surrogatesST');
  
disp(strcat('Data saved to: ',file));

% =========================================================================
% This script loads gait data from Dingwell�s MAT-files and calculates 
% multiple phase-randomized surrogates using TISEAN library. Then it finds 
% piecewise linear MARS trends in surrogate time series. For a given 
% treadmill speed, the normalized trend durations and normalized trend slopes
% are saved to
% '../data/surrogates/surrogates_trend_distributions' folder.
% Before running the script,  please set:
%       1) speed (SPD)
%       2) number of surrogates for each experimental time series (repetitions)
%       3) type of surrogates (cross-correlated true, independent = false)
%
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


clc, clear, close all

SPD = 1; % 1 - 100, 2 - 110, 3 - 90, 4 - 120, 5 - 80 [%PWS]
repetitions = 20; % number of surrogates for each experimental time series
cross_correlated = true;

currentPath = pwd;
addpath('../libs/');
addpath('../utils/');
addpath('../data/mat_data/');

trend_durations_statsSL = cell(repetitions,0);
trend_durations_statsST = cell(repetitions,0);
trend_slopes_statsSL = cell(repetitions,0);
trend_slopes_statsST = cell(repetitions,0);

for x = 1 : repetitions
    trend_durations_statsSL{x} = [];
    trend_durations_statsST{x} = [];
    trend_slopes_statsSL{x} = [];
    trend_slopes_statsST{x} = [];
end


% MARS params
params = aresparams2('useMinSpan',-1,'useEndSpan',-1,'cubic',false,...
    'c',2,'threshold',1e-3,'maxFuncs',50);
    
% load SL data
SLdata = load(strcat('SL_SPD',num2str(SPD),'.mat'));
% load ST data
STdata = load(strcat('ST_SPD',num2str(SPD),'.mat'));
    
s = size(SLdata.residualsAll);

for r = 1 : repetitions
        
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
    
    
    for i = 1 : s(2)

        orgST = STdata.seriesAll{i};
        orgSL = SLdata.seriesAll{i};
        T = STdata.timestampsAll{i};
        
  
        seeds = randi([1 30000],2,1);
         
        if(cross_correlated) % correlated surrogates
       
            str = 'cross_correlated_surrogates';

            tot = [orgSL  orgST]; 
            inFile = strcat(pwd, '\in.txt');
            save(inFile,'tot','-ascii');
            cmd = strcat('surrogates.exe -o',{' '},pwd, ...
                        '\out.txt -m 2 -I',{' '},num2str(seeds(1)),{' '},...
                        pwd,'\in.txt');
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
                '\out.txt -I',{' '},num2str(seeds(1)),{' '},pwd,'\in.txt');
            cd('../libs/');
            system(cmd{1});
            cd(currentPath)
            surSL = load('out.txt'); 
            surrogatesSL{end+1} = surSL;

            inFile = strcat(pwd, '\in.txt');
            save(inFile,'orgST','-ascii');
            cmd = strcat('surrogates.exe -o',{' '},pwd, ...
                '\out.txt -I',{' '},num2str(seeds(2)),{' '},pwd,'\in.txt');
            cd('../libs/');
            system(cmd{1});
            cd(currentPath)
            surST = load('out.txt');
            surrogatesST{end+1} = surST;		
         
        end % end if
            
        % SL surrogates 
         outputData = perform_MARS(T(1:length(surSL)),...
            surSL,false,[],[],params);
        

        surrogatesSL_residualsAll{end+1} = outputData.residuals;
        surrogatesSL_trends{end+1} = outputData.trends;
        timestampsSL{end+1} = outputData.timestamps;
        modelsSL{end+1} = outputData.model;
        infosSL{end+1} = outputData.info;    
        knotIndicesSL{end+1} = outputData.knotIndices;
        valuesAtKnotSL{end+1} = outputData.valuesAtKnot;
        
        data_surrogatesSL.trendsAll = surrogatesSL_trends;
        data_surrogatesSL.seriesAll = surrogatesSL;
        data_surrogatesSL.timestampsAll = timestampsSL;
        data_surrogatesSL.residualsAll = surrogatesSL_residualsAll;
        data_surrogatesSL.modelsAll = modelsSL;
        data_surrogatesSL.infosAll = infosSL;
        data_surrogatesSL.knotIndicesAll = knotIndicesSL;
        data_surrogatesSL.valuesAtKnotAll = valuesAtKnotSL;

    
        % ST surrogates
        outputData = perform_MARS(T(1:length(surSL)),...
            surST,false,[],[],params);

        surrogatesST_residualsAll{end+1} = outputData.residuals;
        surrogatesST_trends{end+1} = outputData.trends;
        timestampsST{end+1} = outputData.timestamps;
        modelsST{end+1} = outputData.model;
        infosST{end+1} = outputData.info;    
        knotIndicesST{end+1} = outputData.knotIndices;
        valuesAtKnotST{end+1} = outputData.valuesAtKnot;

        data_surrogatesST.trendsAll = surrogatesST_trends;
        data_surrogatesST.seriesAll = surrogatesST;
        data_surrogatesST.timestampsAll = timestampsST;
        data_surrogatesST.residualsAll = surrogatesST_residualsAll;
        data_surrogatesST.modelsAll = modelsST;
        data_surrogatesST.infosAll = infosST;
        data_surrogatesST.knotIndicesAll = knotIndicesST;
        data_surrogatesST.valuesAtKnotAll = valuesAtKnotST;
    
    
         % trend stats
        [trend_durationsSL, trend_slopesSL, ~, ~] = ...
            calculate_trend_stats(data_surrogatesSL);

        [trend_durationsST, trend_slopesST, ~, ~] = ...
            calculate_trend_stats(data_surrogatesST);
     
     end % end trial loop
     
     trend_durations_statsSL{r} = [trend_durations_statsSL{r}; ...
         trend_durationsSL];
     trend_durations_statsST{r} = [trend_durations_statsST{r}; ...
         trend_durationsST]; 
     
     trend_slopes_statsSL{r} = [trend_slopes_statsSL{r}; ...
         trend_slopesSL];
     trend_slopes_statsST{r} = [trend_slopes_statsST{r}; ...
         trend_slopesST];

     
end % end repetitions loop    

% save data 
file = strcat('../data/surrogates/surrogates_trend_distributions/distr_',...
              str,'_SPD',num2str(SPD),'.mat');

save(file,'trend_durations_statsSL','trend_durations_statsST',...
    'trend_slopes_statsSL','trend_slopes_statsST');
  
disp(strcat('Data saved to: ',file));

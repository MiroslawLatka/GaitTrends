% =========================================================================
% This script loads Dingwell’s MAT-file and calculates piecewise linear MARS
% trends.  The calculations are performed for a given treadmill speed (SPD)
% and gait parameter (attributeNumber). The output file contains:
% original series, MARS series, MARS residuals, time stamps,
% MARS models, MARS knot indices, and values at MARS knots.

% Before running the script,  please set: speed (SPD) and attribute (SL/ST/SS).
% By default, the variable generateFigures is set to true so that the time
% series and their MARS trends are plotted for all subjects. 
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
attributeNumber = 1; % 1 - SL, 2 - ST, 3 - SS
generateFigures = true;


dataDir = '../data/original_Dingwell_data/';
addpath('../utils/');
fileList = dir(dataDir);
fileList = fileList(~[fileList.isdir]);

switch attributeNumber
    case 1
        param = 'Ln';
    case 2
        param = 'Tn';
    case 3
        param = 'Sn';
    otherwise
        error('Error. attributeNumber must an integer between 1 and 3.')
end

% MARS params
params = aresparams2('useMinSpan',-1,'useEndSpan',-1,'cubic',false,...
         'c',2,'threshold',1e-3,'maxFuncs',50);

residualsAll = {};
trendsAll  = {};
seriesAll  = {};
timestampsAll = {};
modelsAll  = {};
infosAll  = {};
knotIndicesAll = {};
valuesAtKnotAll = {};
        
for i = 1 : length(fileList)
    
    inputFileName = strcat(dataDir,fileList(i).name);
    load(inputFileName);

    % select matrices corresponding to the chosen treadmill speed (SPD)
    switch SPD
        case 1
            data1 = SPD1TR1;
            data2 = SPD1TR2;
        case 2
            data1 = SPD2TR1;
            data2 = SPD2TR2;
        case 3
            data1 = SPD3TR1;
            data2 = SPD3TR2;
        case 4
            data1 = SPD4TR1;
            data2 = SPD4TR2;
        case 5
            data1 = SPD5TR1;
            data2 = SPD5TR2;
        otherwise
            error('Error. SPD must be an integer between 1 and 5.')
    end
    
    % trial 1
    strTrial = '1';
    
    % perform calculations for trial 1
    if length(data1(:,1)) > 1

        % cumulate ST to obtain time stamps
        X = cumsum(data1(:,2))-data1(1,2);
        Y = data1(:,attributeNumber);

        % build MARS model using ARESLab
        model = aresbuild(X,Y,params);
        % calculate trends
        predY = arespredict(model,X);
        % calculate MARS residuals
        detrendedY = Y - predY;

        residualsAll{end+1} = detrendedY;
        trendsAll{end+1} = predY;
        seriesAll{end+1} = Y;
        timestampsAll{end+1} = X;
        modelsAll{end+1} = model;
        infosAll{end+1} = strcat(fileList(i).name,'/trial',strTrial);

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
            valuesAtKnot = predY(knotIndices);
        end

        knotIndicesAll{end+1} = knotIndices;
        valuesAtKnotAll{end+1} = valuesAtKnot;


        if(generateFigures)
			plotTitle = strcat('Y',num2str(i),{'SPD'}, num2str(SPD), ...
				{'trial'}, strTrial);
			fig = figure;
			plot(X,Y,'LineWidth',2); hold on;
			plot(X,predY,'r','LineWidth',2); 
			plot(X(knotIndices), predY(knotIndices),'go',...
                'MarkerSize',8,'LineWidth',2);
			title(plotTitle,'FontSize', 20);
			xlim([min(X) max(X)]);
			xlabel('time [s]','FontSize', 18);
			ylabel(param,'FontSize', 18);
			grid on;
			set(gca,'FontWeight','bold','FontSize', 13);
			set(gcf, 'PaperPositionMode', 'auto');
			hold off;    
        end

    end
        
    % trial 2
    strTrial = '2';
    
    % perform calculations for trial 2
	if length(data2(:,1)) > 1

        % cumulate ST to obtain time stamps
		X = cumsum(data2(:,2))-data2(1,2);
		Y = data2(:,attributeNumber);
		
		% build MARS model using ARESLab
        model = aresbuild(X,Y,params);
        % calculate MARS trends
        predY = arespredict(model,X);
        % calculate MARS residuals 
		detrendedY = Y - predY;
		
		residualsAll{end+1} = detrendedY;
		trendsAll{end+1} = predY;
		seriesAll{end+1} = Y;
		timestampsAll{end+1} = X;
		modelsAll{end+1} = model;
		infosAll{end+1} = strcat(fileList(i).name,'/trial',strTrial);
		
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
			valuesAtKnot = predY(knotIndices);
		end
		
		knotIndicesAll{end+1} = knotIndices;
		valuesAtKnotAll{end+1} = valuesAtKnot;
		
		if(generateFigures)
			plotTitle = strcat('Y',num2str(i),{'SPD'}, num2str(SPD), ...
				{'trial'}, strTrial);
			fig = figure;
			plot(X,Y,'LineWidth',2); hold on;
			plot(X,predY,'r','LineWidth',2); 
			plot(X(knotIndices), predY(knotIndices),'go',...
                'MarkerSize',8,'LineWidth',2);
			title(plotTitle,'FontSize', 20);
			xlim([min(X) max(X)]);
			xlabel('time [s]','FontSize', 18);
			ylabel(param,'FontSize', 18);
			grid on;
			set(gca,'FontWeight','bold','FontSize', 13);
			set(gcf, 'PaperPositionMode', 'auto');
			hold off;    
         end  

	end


end


% save experimental time series and calculated quantities

outputDir = '../data/mat_data/';
file = strcat(outputDir,param,'_SPD',num2str(SPD),'.mat');

save(file,'residualsAll','trendsAll','seriesAll','timestampsAll',...
      'modelsAll','infosAll','knotIndicesAll','valuesAtKnotAll');
  
disp(strcat('Data saved to: ',file));
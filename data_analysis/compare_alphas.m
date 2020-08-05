% =========================================================================
% For each time series at a given treadmill speed, the script 
% generates one of three types of surrogates: a cross-correlated (SL-ST)
% phase-randomized surrogate, an independently phase-randomized surrogate, 
% and shuffled surrogates. Then, it uses either 
% a paired t-test or Wilcoxon signed rank test to compare the madogram scaling 
% index of experimental and surrogate data. This procedure is repeated 
% using different realizations of surrogates.
% Before running the script,  please choose:
%       1) treadmill speed (SPD)
%       2) type of surrogates (surrogates_type)
%       3) number of surrogates for each experimental time series (repetitions). 
% The script prints in MATLAB's command window the number of realizations 
% for which the scaling exponent for the experimental and surrogate data 
% are not different. The boxplots for the experimental and surrogates are
% plotted.
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

SPD = 1; % 1 - 100, 2 - 110, 3 - 90, 4 - 120, 5 - 80 [%PWS]
repetitions = 10; % number of surrogates for each experimental time series
% 1 - independent, 2 - cross-correlated, 3 - shuffled
surrogates_type = 1;


currentPath = pwd;
addpath('../libs/');
addpath('../utils/');
addpath('../data/mat_data/');

if(surrogates_type < 1 || surrogates_type > 3)
    error('Error. surrogates_type must be an integer between 1 and 3.')
end

% load data
% load SL data
SLdata = load(strcat('SL_SPD',num2str(SPD),'.mat'));
% load ST data
STdata = load(strcat('ST_SPD',num2str(SPD),'.mat'));

trialSize = size(SLdata.residualsAll);
MD_OrgSL = zeros(1,trialSize(2));
MD_OrgST = zeros(1,trialSize(2));
MD_SurSL_cell = cell(repetitions,0);
MD_SurST_cell = cell(repetitions,0);

for x = 1 : repetitions
    MD_SurSL_cell{x} = [];
    MD_SurST_cell{x} = [];
end

for i = 1 : trialSize(2)
    
    % calculate madogram for experimental datal
    orgSL = SLdata.seriesAll{i};
    MD_OrgSL(i) = variogram(cumsum(orgSL-mean(orgSL)),1);
    orgST = STdata.seriesAll{i};
    MD_OrgST(i) = variogram(cumsum(orgST-mean(orgST)),1);
       
    % generate proper surrogates
    for r = 1 : repetitions
        
       seeds = randi([1 30000],2,1);
       
        switch surrogates_type
            case 1
                str = 'independent';
                
                inFile = strcat(pwd, '\in.txt');
                save(inFile,'orgSL','-ascii');
                cmd = strcat('surrogates.exe -o',{' '},pwd, ...
                    '\out.txt -I',{' '},num2str(seeds(1)),{' '},pwd,'\in.txt');
                cd('../libs/');
                system(cmd{1});
                cd(currentPath)
                surSL = load('out.txt'); 

                inFile = strcat(pwd, '\in.txt');
                save(inFile,'orgST','-ascii');
                cmd = strcat('surrogates.exe -o',{' '},pwd, ...
                    '\out.txt -I',{' '},num2str(seeds(2)),{' '},pwd,'\in.txt');
                cd('../libs/');
                system(cmd{1});
                cd(currentPath)
                surST = load('out.txt');
	
            case 2
                str = 'cross-correlated';

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


            case 3
                str = 'shuffled';
                
                % shuffled surrogates
                surST =  orgST(randperm(length(orgST)));
                surSL =  orgSL(randperm(length(orgSL)));


            otherwise
                error('Error. Attribute must be value between 1 and 3.')
        end % end switch
        
        % calculate madogram for experimental data 
         MD_SurSL_cell{r} = [MD_SurSL_cell{r}; ...
             variogram(cumsum(surSL-mean(surSL)),1);];
         MD_SurST_cell{r} = [MD_SurST_cell{r}; ...
             variogram(cumsum(surST-mean(surST)),1);];
         
    end % end repetitions loop
    
   
end

% compare data
meanMD_SurSL = zeros(1,repetitions);
meanMD_SurST = zeros(1,repetitions);
countSL = 0;
countST = 0;

for r = 1 : repetitions
    
    meanMD_SurSL(r) = mean(MD_SurSL_cell{r});
    meanMD_SurST(r) = mean(MD_SurST_cell{r});
    
    % Shapiro-Wilk test for normality
    [hOSL, pOSl] = swtest(MD_OrgSL);
    [hSSL, pSSL] = swtest(MD_SurSL_cell{r});
    [hOST, pOST] = swtest(MD_OrgST);
    [hSST, pSST] = swtest(MD_SurST_cell{r});
    
     if hOSL + hSSL == 0 % t-test
        [hSL,pSL] = ttest(MD_SurSL_cell{r},MD_OrgSL','tail','both');
     else % Wilcoxon signed rank test
        [pSL,hSL] = signrank(MD_SurSL_cell{r},MD_OrgSL','tail','both');
     end
     
     if hSL == 0
         countSL = countSL + 1;
     end
     
     if hOST + hSST == 0 % t-test
        [hST,pST] = ttest(MD_SurST_cell{r},MD_OrgST','tail','both');
     else % Wilcoxon signed rank test
        [pST,hST] = signrank(MD_SurST_cell{r},MD_OrgST','tail','both');
     end
     
     if hST == 0
         countST = countST + 1;
     end
        
end

% visualize data
datSL = MD_OrgSL;
datST = MD_OrgST;
group = 1*ones(length(MD_OrgSL),1);
labels = {'exp'};

for l = 1 : repetitions
    datSL = [datSL  MD_SurSL_cell{l}'];
    datST = [datST  MD_SurST_cell{l}'];
    group = [group; (l+1)*ones(length(MD_SurSL_cell{r}),1)];
    labels{end+1} = strcat('s',num2str(l));
end

figure;
boxplot(datSL,group);
xlabel('data type');
ylabel('{\alpha}^(^M^D^)');
title('SL');
set(gca,'XTickLabel',labels,'FontWeight','bold','FontSize', 13);
grid on;
set(gcf, 'PaperPositionMode', 'auto');
set(findobj(gca,'type','line'),'linew',1.5);
hold off; 

figure;
boxplot(datST,group);
xlabel('data type');
ylabel('{\alpha}^(^M^D^)');
title('ST');
set(gca,'XTickLabel',labels,'FontWeight','bold','FontSize', 13);
grid on;
set(gcf, 'PaperPositionMode', 'auto');
set(findobj(gca,'type','line'),'linew',1.5);
hold off; 


disp('same mean/median:')
disp(strcat('SL:',{' '},num2str(countSL),'/',num2str(repetitions)))
disp(strcat('ST:',{' '},num2str(countST),'/',num2str(repetitions)))
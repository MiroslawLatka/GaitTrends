% =========================================================================
% This script calculates:
% 1) coefficient of variation (COV) for gait parameters (SL/ST/SS),
% 2) trend speed and its COV,
% 3) trend speed control.
% The boxplots of these quantities are plotted for each treadmill speed.
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

COV_trend_speed_cell = {};
COV_SL = {};
COV_ST = {};
COV_SS = {};
speed_control = {};

for i = 1 : 5 % 1-5 SPD
    
	% load SL data
	SLdata = load(strcat('Ln_SPD',num2str(i),'.mat'));
	% load SL data
	STdata = load(strcat('Tn_SPD',num2str(i),'.mat'));
	% load SL data
	SSdata = load(strcat('Sn_SPD',num2str(i),'.mat'));

	trialSize = size(SLdata.residualsAll);
	COV_trend_speed_vec = [];
	COV_SL_vec = [];
	COV_ST_vec = [];
	COV_SS_vec = [];
	speed_control_vec = [];

	for j = 1 : trialSize(2)
		
		tr_sl = SLdata.trendsAll{j};
		tr_st = STdata.trendsAll{j};
		sl = SLdata.seriesAll{j};
		st = STdata.seriesAll{j};
		ss = SSdata.seriesAll{j};
		x = 1:length(tr_sl);
        % calculate trend speed
		trendSpeed = tr_sl./tr_st;


		% calculate COVs
		COV_SL_vec =[COV_SL_vec; 100*std(sl)/mean(sl)];
		COV_ST_vec =[COV_ST_vec; 100*std(st)/mean(st)];
		COV_SS_vec =[COV_SS_vec; 100*std(ss)/mean(ss)];
		COV_trend_speed_vec = [COV_trend_speed_vec; ...
			100*std(trendSpeed)/mean(trendSpeed)];

		% calculate speed control
        speed_control_vec = [speed_control_vec; sum((trendSpeed-mean(ss)).^2)/...
			sum((ss-mean(ss)).^2)];


	end % end trial loop

	COV_trend_speed_cell{end+1} = COV_trend_speed_vec;
	COV_SL{end+1} = COV_SL_vec;
	COV_ST{end+1} = COV_ST_vec;
	COV_SS{end+1} = COV_SS_vec;
	speed_control{end+1} = speed_control_vec;


end % end speed loop

% visualize data

% speed control
boxplots_for_all_speeds(speed_control, 'experimental data',...
    'trend speed control');

% COV trend speed
boxplots_for_all_speeds(COV_trend_speed_cell, 'experimental data',...
    'COV trend speed [%]');

% COV SL
boxplots_for_all_speeds(COV_SL, 'experimental data',...
    'COV SL [%]');

% COV ST
boxplots_for_all_speeds(COV_ST, 'experimental data',...
    'COV ST [%]');

% COV SS
boxplots_for_all_speeds(COV_SS, 'experimental data',...
    'COV SS [%]');

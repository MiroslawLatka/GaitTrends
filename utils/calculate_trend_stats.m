function [trend_durations, trend_slopes, long_trend_durations, long_trend_slopes] = ...
         calculate_trend_stats(data, threshold)
% calculate the normalized duration and slope of MARS trends in gait time
% series (ST/SL/SS).
%
% Inputs:
% data: the collection of cells returned by one of the following scripts:
%       ../data_preparation/prepare_data.m,
%       ../data_preparation/prepare_surrogates.m,
%       ../data_preparation/prepare_shuffled_surrogates.m
% threshold (optional): if this optional argument has been passed 
%       the function also returns the normalized durations and slopes of all
%       MARS trends longer than the threshold. 
%
% Outputs:
% trend_durations: normalized durations of MARS trends (vector)
% trend_slopes: normalized slopes of MARS trends (vector)
% long_trend_durations:  normalized durations of MARS trends longer than
%                        threshold (vector)
% long_trend_slopes: normalized slopes of MARS trends longer than
%                    threshold (vector)

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

trialSize = size(data.residualsAll);
trend_durations = [];
trend_slopes = [];
long_trend_durations = [];
long_trend_slopes = [];
trendDurationsLong = [];
trendSlopesLong = [];

for j = 1 : trialSize(2)
    
    % MARS knot indices
    kia = data.knotIndicesAll{j};

    % check if there are any knots
    if ~isempty(kia)

        vaka = data.valuesAtKnotAll{j};

        % remove duplicate MARS knots (this may happen ocassionally) 
        [~, ind] = unique(kia);

        knotInd = kia(ind);

        if(iscell(data.timestampsAll{j}))
            timeStamps = cell2mat(data.timestampsAll{j});
        else
            timeStamps = data.timestampsAll{j};  
        end

        series = data.seriesAll{j};

        knotSites = timeStamps(knotInd);

        % recover stride times from time stamps
        STSeries = diff(timeStamps);

        % calculate trend durations
        trendDurations = diff(knotSites)/mean(STSeries);

        valuesAtKnot = vaka(ind);

        trendSlopes = [];

        % calculate trend slopes
        for i = 1 : length(trendDurations)
            trendSlopes = [trendSlopes; 
                (valuesAtKnot(i+1)-valuesAtKnot(i))/...
                (mean(data.seriesAll{j})*trendDurations(i))];    
        end

       % find MARS trends longer than threshold
       if nargin > 1
           
           locs = [];
           for i = 1 : length(trendDurations)    
                if(trendDurations(i) > threshold && abs(trendSlopes(i)) < 0.001)
                    locs = [locs; i];      
                end    
           end

           if (length(locs) > 0)

               trendDurationsLong = trendDurations(locs);
               trendSlopesLong = trendSlopes(locs);

               for l = 1 : length(locs)  

                   timeFragment = timeStamps(knotInd(locs(l)):knotInd(locs(l)+1));
                   seriesFragment = series(knotInd(locs(l)):knotInd(locs(l)+1));

               end

           end
           
          long_trend_durations = [long_trend_durations; trendDurationsLong];
          long_trend_slopes = [long_trend_slopes; trendSlopesLong];
           
       else % no threshold set
           
           long_trend_durations = [];
           long_trend_slopes = [];
           
       end
       
       	trend_durations = [trend_durations; trendDurations];
		trend_slopes = [trend_slopes; trendSlopes];
      
        
    else % no knots

        trendDurations = []; 
        trendSlopes = [];
        valuesAtKnot = [];
        trendDurationsLong = [];
        trendSlopesLong = [];  

    end
 
  
end % end trial loop


end


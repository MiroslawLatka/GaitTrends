function boxplots_for_all_speeds(cell_data,title_str, ylabel_str)
% Makes boxplots of five vectors stored in a cell. Used to visualize
% the quantities calculated from Dingwell's experimental data. 
% 
% Inputs:
% cell_data:  cell with five vectors
% title_str:  plot's title
% ylabel_str: y-axis label
%
% Outputs:    none     
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
% Last update: July 21, 2020
% =========================================================================

% Citing the GaitTrends:
% https://doi.org/10.1101/677948

ind = [5, 3, 1, 2, 4]; % 5 - 80, 3 - 90, 1 - 100, 2 - 110, 4 - 120 [%PWS]
dat = [];
group = [];
speeds = {'80','90','100','110','120'};

for l = 1 : length(ind)
    dat = [dat; cell_data{ind(l)}];
    group = [group; l*ones(size(cell_data{ind(l)}))];
end

figure;
boxplot(dat,group);
xlabel('treadmill speed [%PWS]');
ylabel(ylabel_str);
title(title_str);
set(gca,'XTickLabel',speeds,'FontWeight','bold','FontSize', 13);
grid on;
set(gcf, 'PaperPositionMode', 'auto');
set(findobj(gca,'type','line'),'linew',1.5);
hold off; 

end


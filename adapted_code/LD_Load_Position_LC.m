function [POS,POS_maze,Pos_labels] = LD_Load_Position_LC(pos_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the position data.
% - returns in microseconds as that's what it is saved as.
% Cowen 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    pos_file = 'PositionStruct.mat';
end
%%
GP = LD_Globals_LC;
load(pos_file);
POS = [Range(pos.All.X)*100 Data(pos.All.X) Data(pos.All.Y)];
% Compute speed as the 4th column
S = Speed_from_xy([POS(:,1)/1e6 POS(:,2:3)*GP.cm_per_pixel]); %make sure output is in cm/sec

POS = [POS S];

%sFreq = 1e6/median(diff(POS(:,1)))
if nargout > 1
    POS_maze = cell(2,1);
    for ii = 1:2
        theta = Data(pos.Maze(ii).Theta);
        d_theta = movmedian(diff([theta(1);theta]),30);
        IX = d_theta < 0;
        new_theta = theta;
        new_theta(IX) = 720 - new_theta(IX);
        %         heading = movmedian(diff(Data(pos.Maze(ii).Theta) ),60);
        POS_maze{ii} = [Range(pos.Maze(ii).X)*100 Data(pos.Maze(ii).X) ...
            Data(pos.Maze(ii).Y) Data(pos.Maze(ii).Theta) Data(pos.Maze(ii).Radius) new_theta];
    end
    %sFreq = 1e6/median(diff(POS_maze{1}(:,1)))
    Pos_labels = {'t_usec' 'x' 'y' 'theta' 'radius' 'new_theta'};
end

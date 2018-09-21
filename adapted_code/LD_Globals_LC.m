function GP = LD_Globals_LC()
% Contains global data used by all functions for the LD analysis.
[GP.Root_dir, GP.SessionList_dir, GP.Data_dir, GP.Processed_data_dir ] = LD_get_directories_LC;
GP.Bad_LFP_Rats = [8700]; % bad lfp on this rat.
% NOTE: in original Session list files - it listed that CSC13 was the EMG
% for 8957. This is not true. It was CSC14. I fixed the master session
% lists to correct for this.
GP.Young_remapping_rats = [8419 8820];
GP.Old_remapping_rats = [8417 8886 8981]; % 8700 is a bad rat
GP.Young_non_remapping_rats = [8570 8645 8646 8957];
GP.Old_non_remapping_rats = [8564 8778];
GP.Young  = [8419 8570 8645 8646 8820 8957];
GP.Old = [8417 8564 8700 8778 8886 8981];
GP.Sessions_with_bad_EMG = [8646 15];
GP.Watermaze =[];
GP.Color.Old = [.7 0 .7];
GP.Color.Young = [.1 .7 .1];
GP.Color.Shuffle = [0.6 0.6 0.6];
GP.Maze_Diameter_cm = 85;
GP.Maze_Circumference_cm = 2*pi*(GP.Maze_Diameter_cm /2);
GP.Maze_Diameter_pixels = 392;
GP.cm_per_pixel = GP.Maze_Diameter_cm/GP.Maze_Diameter_pixels; % conversion from camera pixels to cm.
GP.cm_per_degree = GP.Maze_Circumference_cm/360; % ?? CM per degree. 
GP.Tracking_Sample_Rate_Hz = 30;
% How many degrees makes up 20 cm on the track.
GP.Control_Zone_Size_cm = 22.5; % Note this is 30 degrees. How convenient.
GP.Control_Zone_Size_deg = GP.Control_Zone_Size_cm/GP.cm_per_degree; % Note this is 30 degrees. How convenient.
GP.Lin_Pos_Bounds = [7 353];
GP.RewardZoneLinPos = [36 315];
GP.StartZoneLinPos = [80 270];
GP.StimZoneLinPos = [140 256]; % eyeballed from TI data.
GP.Control_Zone_Pos(1) = GP.StimZoneLinPos(1) + diff( GP.StimZoneLinPos)/2 - GP.Control_Zone_Size_deg/2 ;
GP.Control_Zone_Pos(2) = GP.Control_Zone_Pos(1) +  GP.Control_Zone_Size_deg;


% 
% [Rat	Group	Age	21	22	23	24	Average	MeanModulationDepth
% 8419	Y	10mo	17.2	28.8	2.1	6.7	13.7	86.42987246
% 8570	Y	10mo	5.1	7	38.4	11.4	15.475	50.62691357
% 8645	Y	9mo	0.6	331.7	0.8	0.2	83.325	73.51039195
% 8646	Y	9mo	0.8	60.6	63.4	-1.4	30.85	82.85545479
% 8820	Y	8-10mo?	400.4	330.7	310.6	168.1	302.45	60.48531732
% 8957	Y	8mo	0.2	23.4	24.8	1.5	12.475	56.88262714
% 								
% 								
% 								
% 								
% 								MeanModulationDepth
% 8417	O	24mo	408.8	446.2	539.1	344.5	434.65	77.31658965
% 8564	O	24mo	1.4	49.3	153.8	110.4	78.725	60.23067907
% 8700	O	25mo	403.1	16.9	451.5	420.1	322.9	61.93649096
% 8778	O	25mo	450.2	472.2	0.4	488.4	352.8	51.79506744
% 8886	O	24-25mo?	60.7	103.3	372	149.1	171.275	56.03799494
% 8981	O	25mo	447	431.5	44.6	502.3	356.35	47.18550324
% ]
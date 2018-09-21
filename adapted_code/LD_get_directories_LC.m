function [ROOTDIR, SES_DIR, DATA_DIR, PROC_DIR] = LD_get_directories_LC()
%
% JP: You are not allowed to edit this file IF it is in the cowen_code
% subdirectory. Modify your own copy please.
%
cur_dir = pwd;
[p] = fileparts(which('LD_Session_Iterator'));  
cd (p)
cd ..
ROOTDIR = pwd;

%
cd .\SessionList_files\Updated_SessionLists\
SES_DIR = pwd;
cd(cur_dir)

% Determine the directory where the data lives.
DATA_DIR = [];
if exist('C:\Cowen\Data\Lesley_Eyeblink\Lipa-HDD\Eyeblink_Cut-Data','dir')
    DATA_DIR = 'C:\Cowen\Data\Lesley_Eyeblink\Lipa-HDD\Eyeblink_Cut-Data';
elseif exist('F:\Lesley_Data\Eyeblink_Cut-Data','dir')
    DATA_DIR = 'F:\Lesley_Data\Eyeblink_Cut-Data';
elseif exist('E:\Cowen\Data\Lesley_Eyeblink\Lipa-HDD\Eyeblink_Cut-Data','dir')
    DATA_DIR = 'E:\Cowen\Data\Lesley_Eyeblink\Lipa-HDD\Eyeblink_Cut-Data';
elseif exist('G:\Lesley_Data\Eyeblink_Cut-Data','dir')
    DATA_DIR = 'G:\Lesley_Data\Eyeblink_Cut-Data';
elseif exist('E:\Lipa-HDD\Eyeblink_Cut-Data\Lipa-HDD\Eyeblink_Cut-Data','dir')
    DATA_DIR = 'E:\Lipa-HDD\Eyeblink_Cut-Data\Lipa-HDD\Eyeblink_Cut-Data';
elseif exist('C:\Users\dtgray\Documents\MATLAB\Lesley_Data\Processed_data','dir')
    DATA_DIR = 'C:\Users\dtgray\Documents\MATLAB\Lesley_Data\Processed_data';
elseif exist('C:\Cowen\Data\Lesley_Eyeblink\Lipa-HDD\Eyeblink_Cut-Data','dir')
    DATA_DIR = 'C:\Cowen\Data\Lesley_Eyeblink\Lipa-HDD\Eyeblink_Cut-Data';
elseif exist('E:\Lindsey\Lesley_Data\Eyeblink_Cut-Data','dir')
    DATA_DIR = 'E:\Lindsey\Lesley_Data\Eyeblink_Cut-Data';
end
% ANA_DIR = fullfile(Analysis_dir(),'Lesley_data');
% 
% 
% if exist('C:\Users\Cowen.Stephen\Dropbox\Foldershare\Analysis_Results_Dropbox\Lesley_Data','dir')
%     ANA_DIR = 'C:\Users\Cowen.Stephen\Dropbox\Foldershare\Analysis_Results_Dropbox\Lesley_Data';
% elseif exist('C:\Cowen\Analysis_Results','dir')
%     ANA_DIR = 'C:\Cowen\Analysis_Results';
% elseif exist('E:\JP\Lesley Data (original)\Lipa-HDD\Eyeblink_Cut-Data\New_processed_data\Summary_Files','dir')
%     ANA_DIR = 'E:\JP\Lesley Data (original)\Lipa-HDD\Eyeblink_Cut-Data\New_processed_data\Summary_Files';
% else 
%     error('Could not find your output results directory')
% end

PROC_DIR = fullfile(DATA_DIR,'Processed_Data');
if ~exist(PROC_DIR,'dir')
    msgbox(['You should create ' PROC_DIR ' and fill it with processed data'])
end
% 
% PLACEFIELD_DIR = fullfile(DATA_DIR,'\Data1Volume8-Aging_Eye_Blink\AnalyzedDataArchive\LesleyPlaceField_Results\');
% if ~exist(PLACEFIELD_DIR,'dir')
%     msgbox(['No ' PLACEFIELD_DIR ' on your computer.'])
% end
% 
% PLACEFIELD_DIR = fullfile(DATA_DIR,'\Data1Volume8-Aging_Eye_Blink\AnalyzedDataArchive\LesleyPlaceField_Results\');
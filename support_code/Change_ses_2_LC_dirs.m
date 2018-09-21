function [newses]=Change_ses_2_LC_dirs(ses)

%Session Info gives you a structure called ses that has
% 
%           positionStruct: 'F:\Lesley_Data\Eyeblink_Cut-Data\8417_cut\Jie\Sessions_1-2done\PositionStruct.mat'
%        ThetaEEG_FileName: 'F:\Lesley_Data\Eyeblink_Cut-Data\8417_cut\Jie\Sessions_1-2done\CSC01.ncs'
%           BestRippleFile: 'F:\Lesley_Data\Eyeblink_Cut-Data\8417_cut\Jie\Sessions_1-2done\CSC11.ncs'
%     SecondBestRippleFile: 'F:\Lesley_Data\Eyeblink_Cut-Data\8417_cut\Jie\Sessions_1-2done\CSC08.ncs'
%                  EMGfile: 'F:\Lesley_Data\Eyeblink_Cut-Data\8417_cut\Jie\Sessions_1-2done\CSC13.ncs'
%                Eventfile: 'F:\Lesley_Data\Eyeblink_Cut-Data\8417_cut\Jie\Sessions_1-2done\Events.Nev'
%      BoundingBoxFileName: 'BoundingBox_Room_311A.mat'

fields_all=fieldnames(ses);
for ifield=1:numel(fields_all)
    if ~isequal(ses.(fields_all{ifield}),ses.other)
        newses.(fields_all{ifield})=ses.(fields_all{ifield});
    end
end

answ=ses.tfiledir;
peice=strsplit(answ,'\');
peice{1}='E:\Lindsey';
newfile=strjoin(peice,'\');
newses.tfiledir=newfile;

fields=fieldnames(ses.other);
for ifield=1:numel(fields)
   currentf= sprintf('ses.other.%s', fields{ifield});
   if ~isequal(fields{ifield},'BoundingBoxFileName')
   path=eval(currentf);
   bits=strsplit(path,'\');
   bits{1}='E:\Lindsey';
   newpath=strjoin(bits,'\');
   newses.other.(fields{ifield})=newpath;
   end
   newses.other.BoundingBoxFileName=ses.other.BoundingBoxFileName;
end
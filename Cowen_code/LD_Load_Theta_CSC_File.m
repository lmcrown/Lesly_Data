function [LFP,head, artifact_times_usec, BIX] = LD_Load_Theta_CSC_File(theta_fname, time_interval, downsample_fq, Stim_times_usec, artifact_thresh)

if ~exist(theta_fname,'file')
    if exist([theta_fname '.gz'],'file')
        gunzip([theta_fname '.gz'])
        disp(['uncompressed ' theta_fname])
        delete([theta_fname '.gz'])
    elseif exist([theta_fname '.zip'],'file')
        unzip([theta_fname '.zip'])
        disp(['uncompressed ' theta_fname])
        delete([theta_fname '.zip'])
    else
        error('could not find file')
        theta_fname
    end
end

[LFP,~,sFreq] = nlx2matCSC_Matrix(theta_fname);
head = Read_nlx_header(theta_fname);
% get rid of slow DC shifts. (Should we reall do this?) I did play with
% this and movmedian. they are all similar for different widths ranging
% from 2 to 20 sec.
 LFP(:,2) = LFP(:,2) -  movmean(LFP(:,2),round(sFreq*4));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert to uVolts.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LFP(:,2) = LFP(:,2) * head.ADBitVolts * 1e6;
head.LFP_units = 'uVolts';
head.sFreq = downsample_fq;
% Restrict.
if nargin > 1
    LFP = Restrict(LFP,time_interval);
end
BIX = [];
artifact_times_usec = [];
if nargin > 2
    % do additional artifact removal and resampling...
    [BIX,artifact_times_usec] = LD_Clean_LFP(LFP,Stim_times_usec,artifact_thresh,downsample_fq);
    perc_bad = sum(BIX)/length(BIX);
    fprintf('BAD percent: %2.2f\n',perc_bad*100)
    if perc_bad > .2
        error('Too much bad data')
    end
    if 0
        figure
        stem(BIX)
        plot(LFPep(:,2))
    end
    % Zero out the artifacts.
    LFP(BIX,2) = 0;  % could do nans
    % Reduce to a reasonable sampling rate. Since we may be interested in
    % ripples, I think 500 would be the highest we should go. Go lower if you
    % just want up to 100 Hz (say 250).
    
    [y,t] = resample(LFP(:,2),LFP(:,1)/1e6,downsample_fq);
    t = t*1e6;
    GIX = t>LFP(1,1);
    t = t(GIX); y = y(GIX);
    %     figure
    %     plot(LFPep(:,1),LFPep(:,2),t,y)
    LFP = zeros(length(y),2);
    LFP(:,1) = t;
    LFP(:,2) = y;
    
end
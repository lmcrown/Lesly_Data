function [ OUT ] = Q2_phase_precession(ses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Q1_theta_gamma_coupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLOT_IT = true;
if nargin == 0
    load('Session_Info.mat')
else
    save('Session_Info.mat','ses')
end
GP = LD_Globals;
EPOCHS_NAMES = {'Maze1' 'Maze2'};
EPOCHS_TO_LOAD = [2 5];
downsample_fq = 400;
NAR = 120;
% run_speed_thresh = 5; % minium running speed to be considered running.
speed_intervals = [.5 2; 10 20; 25 35; 40 50];
artifact_thresh = 1000; % This may not work on some sets for reasons I do not understand.

fqs = 1:.2:downsample_fq/2;
low_fqs = 2:.5:11;
high_fqs = 12:.5:downsample_fq/2;
% Create a theta filter, but center it on the peak of the theta band.

N = 8;      % Order of the filter.
bpFilt = designfilt('bandpassiir','FilterOrder',N, ...
    'HalfPowerFrequency1',min(high_fqs),'HalfPowerFrequency2',max(high_fqs)-10, ...
    'SampleRate',downsample_fq,'DesignMethod' ,'butter'); % can't go to nyquist - has to be a little lower.

%%
try
    TR = load(fullfile(GP.Processed_data_dir,[ses.animal '_' ses.name '_TrialInfo.mat']));
catch
    disp('Failed to load trialinfo file. Trying to create again.')
    LD_Create_trial_info_file
    TR = load(fullfile(GP.Processed_data_dir,[ses.animal '_' ses.name '_TrialInfo.mat']));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load position and event data and epoch.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('epochs.mat') % NOTE: these are in 0.1msec units.
[EVT, STIMTIMES] = LD_Load_Events; % these are the starts of the stims.
% All_Stim_times_usec = unique([EVT.ts_usec])';
STIM_TIMES.Maze1 = STIMTIMES{2};
STIM_TIMES.Maze2 = STIMTIMES{5};
[POS,POS_Lin] = LD_Load_Position;
POS_Lin = [POS_Lin{1}(:,[1 4]);POS_Lin{2}(:,[1 4])];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE - that POS_Lin already breaks thing down into heading.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sl = interp1(POS(:,1),POS(:,4),POS_Lin(:,1));
POS_Lin = [POS_Lin sl];
POS_Lin = sortrows(POS_Lin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT.ses = ses;
OUT.epochs = epochs;
OUT.fqs = fqs;
OUT.aborted = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the theta electrode. Debating wither to use the fissure, but
% it's often not good. But if we only do it once in a while and use CA1,
% then wer are not consistent.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_fname = ses.other.ThetaEEG_FileName{1};
if ses.other.ThetaEEG_ThetaDelta(1) < 1 % not sure if this is the best threshold.
    OUT.aborted = true;
    return
end
% if ses.other.
% if strfind(lower(theta_fname),'csc15')
%     theta_fname = ses.other.ThetaEEG_FileName{2};
% end
[~,nm] = fileparts(theta_fname);
OUT.theta_ch_num = str2double(nm(end-1:end));
[LFP,head,artifact_times_usec] = LD_Load_Theta_CSC_File(theta_fname,[0 inf], downsample_fq, [STIM_TIMES.Maze1; STIM_TIMES.Maze2],artifact_thresh);
not_artifact_times_usec(:,1) = [LFP(1,1); artifact_times_usec(:,2)];
not_artifact_times_usec(:,2) = [artifact_times_usec(:,1); LFP(end,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the spike data. ELIMINATE INTERNEURONS! This just keeps things
% easier to interpret. This should be analyzed in isolation later on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SS = Load_tfiles(ses.tfilefullnames);
ClustSum = load(fullfile(GP.Processed_data_dir,[ses.animal '_' ses.name '_ClustSummarys.mat']));
S = struct([]); % Create a new structure that has the data relevant for this analysis.
for iS = 1:length(ClustSum.S)
    S(iS).TS = ClustSum.ClustSummary(iS).t;
    S(iS).Type_ID = ClustSum.ClustSummary(iS).Type_ID;
    S(iS).Type_label = ClustSum.ClustSummary(iS).Type_label;
    S(iS).RateRest1 = ClustSum.ClustSummary(iS).rate_Rest1;
end
IX_P = [S.Type_ID]==2;
% IX_IN = [S.Type_ID]==1;

S_PYR = S(IX_P);
TS = {S_PYR.TS};
OUT.S_PYR_Original = S_PYR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 60 Hz is not typically that bad. The notch does not get rid of harmonics.
%  notch_filt = Notch_filter_cowen(downsample_fq);
%  LFP(:,2) = filtfilt(notch_filt,LFP(:,2));
%  notch_filt = Notch_filter_cowen(downsample_fq,119,121);
%  LFP(:,2) = filtfilt(notch_filt,LFP(:,2));
% remember to zero out stim times.
%%
for iEpoch = 1:length(EPOCHS_TO_LOAD)
    ep_name = EPOCHS_NAMES{iEpoch};
    cur_epoch = EPOCHS_TO_LOAD(iEpoch);
    OUT.epoch(iEpoch) = cur_epoch;
    %     OUT.RF.(ep_name) = RF{cur_epoch};
    EPOCH_TIMES_US = epochs.(ep_name)*100;
    % Restrict position and LFP to this epoch
    LFPep = Restrict(LFP,EPOCH_TIMES_US);
    % find the peak theta frequency for running and use this as the center
    % for the filter...
    
    sFreq = 1e6/median(diff(LFPep(:,1)));
    POSep = Restrict(POS,EPOCH_TIMES_US);
    GIX = POSep(:,4) >= 10 & POSep(:,4) < 50;
    %     pwelch(LFPep(GIX,2),sFreq*2,sFreq,[],sFreq)
    th_fq = 2:.2:20;
    [pxx,f]=pburg(LFPep(GIX,2),51,th_fq,sFreq);
    [pxpow,ix] = findpeaks(pxx);
    theta_fq = th_fq(ix(1));
    
    N = 8;      % Order of the filter.
    bpFilt = designfilt('bandpassiir','FilterOrder',N, ...
        'HalfPowerFrequency1',theta_fq-2,'HalfPowerFrequency2',theta_fq+2, ...
        'SampleRate',sFreq,'DesignMethod' ,'butter'); % can't go to nyquist - has to be a little lower.
    %     freqz(bpFilt)
    
    LFPep(:,3) = filtfilt(bpFilt,LFPep(:,2));
    BIX = LFPep(:,2) == 0;
    LFPep(BIX,3) = 0;
    %     LFPep(:,4)'
    LFPep(:,4) = envelope(LFPep(:,3));
    LFPep(BIX,4) = 0;
    LFPep(:,5) = angle(hilbert(LFPep(:,3)));
    LFPep(BIX,5) = nan;
    GIX = POSep(:,4) >= 10 & POSep(:,4) < 50;
    
    intervals = find_intervals(GIX, .5);
    dif_sec = (intervals(:,2) - intervals(:,1))/30; % convert to time.
    intervals = intervals(dif_sec > 0.5,:);
    st = POSep(intervals(:,1),1);
    ed = POSep(intervals(:,2),1);
    LFPrun = Restrict(LFPep,[st ed]);
    POSrun = Restrict(POSep,[st ed]);
    SPrun = Restrict(TS,[st ed]/100);
    % Goal: spkphase and spkpos
    pos_pha = [];
    % for kicks
    for iCell = 1:length(SPrun)
        SPrun{iCell} = Burst_detector(SPrun{iCell},50*10,2);
        nspk(iCell) = length( SPrun{iCell} );
    end
    pos_pha = [];
    for iCell = 1:length(SPrun)
        if nspk(iCell) > 5
        pos_pha{iCell}(:,1) = SPrun{iCell}*100;
        pos_pha{iCell}(:,2) = interp1(POS_Lin(:,1),POS_Lin(:,2),SPrun{iCell}*100);
        pos_pha{iCell}(:,3) = interp1(LFPep(:,1),LFPep(:,5),SPrun{iCell}*100);
        else
            pos_pha{iCell} = [];
        end
    end
    % Restrict by direction
    xpos = 10:10:340;
    for iDirect = 1:2
        IX = TR.TI.(ep_name).TRL.TrlTypeAndZone == iDirect;
        new_pos_pha = Restrict(pos_pha, TR.TI.(ep_name).TRL.TrlStartEnd(IX,:));
        
        for iCell = 1:length(SPrun)
            if Rows(new_pos_pha{iCell}) < 10
                continue
            end
            figure
            subplot(4,1,1:2)
            plot( new_pos_pha{iCell}(:,2), new_pos_pha{iCell}(:,3),'.k')
            
            hold on
            plot( new_pos_pha{iCell}(:,2), new_pos_pha{iCell}(:,3)+ 2*pi,'.k')
            set(gca,'XLim',xpos([1 end]))
            xlabel('pos')
            ylabel('phase')
            title(sprintf('Direc %d',iDirect))
            subplot(4,1,3)
            histogram(new_pos_pha{iCell}(:,2),xpos);
            set(gca,'XLim',xpos([1 end]))

            xlabel('pos')
            subplot(4,1,4)
            histogram([new_pos_pha{iCell}(:,3); new_pos_pha{iCell}(:,3) + 2*pi],40);
            xlabel('phase')
            
            %          axis tight
        end
    end
    
    %      yup2 = envelope_cowen(LFPep(:,3));
    % Just include running times.
    for iInterval = 1:Rows(speed_intervals)
        GIX = POSep(:,4) >= speed_intervals(iInterval,1) & POSep(:,4) < speed_intervals(iInterval,2);
        GIX = convn(GIX,ones(20,1),'same')>0;
        if sum(GIX) < 200
            continue
        end
        intervals = find_intervals(GIX, .5);
        % get rid of really short intervals.
        dif_sec = (intervals(:,2) - intervals(:,1))/30; % convert to time.
        intervals = intervals(dif_sec > 0.5,:);
        st = POSep(intervals(:,1),1);
        ed = POSep(intervals(:,2),1);
        LFPrun = Restrict(LFPep,[st ed]);
        POSrun = Restrict(POSep,[st ed]);
        if 0
            figure
            plotyy(LFPrun(:,1),LFPrun(:,2),POSrun(:,1),POSrun(:,4))
            plot_markers_simple(Stim_times_usec)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the PSds of running times.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [lpsd_welch] = pwelch(LFPrun(:,2), 4*sFreq,0,fqs,sFreq);
        lpsd_welch = 10*log10(lpsd_welch);
        [lpsd] = pyulear(LFPep(:,2),NAR,fqs,sFreq);
        lpsd = 10*log10(lpsd);
        OUT.PSD(iEpoch,iInterval,:) = lpsd;
        OUT.PSDwelch(iEpoch,iInterval,:) = lpsd_welch;
        
        if 0
            % test some parameters for the yulear...
            params = 34:4:250;
            %     x = arrayfun(@(a) pyulear(LFPep(:,2),a,fqs,downsample_fq), params)
            for ii = 1:length(params)
                pyy(ii,:) = pyulear(LFPrun(:,2),params(ii),fqs,sFreq);
            end
            figure
            imagesc(fqs,params,log(pyy))
            axis xy
        end
        
        %         lspsd = sgolayfilt(lpsd,3,9); % not really necessary with the pyulear as it's smooth already.
        [pkpow,pks_ix] = findpeaks(lspsd);
        [pow_sort,six] = sort(pkpow, 'descend' );
        pks_ix = pks_ix(six);
        ix = 1:min([length(pks_ix) 2]);
        
        OUT.peak_fqs(iEpoch,iInterval,1:2) = nan;
        OUT.peak_pow(iEpoch,iInterval,1:2) = nan;
        
        OUT.peak_fqs(iEpoch,iInterval,ix) = fqs(pks_ix(ix));
        OUT.peak_pow(iEpoch,iInterval,ix) = pow_sort(ix);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CFC!. We might have to do this for each interval and then do a
        % weighted mean of these values - at least for the dupre method.
        O = SPEC_cross_fq_coupling_comod_dupre2017(LFPrun(:,2),sFreq,low_fqs,CFC_method);
        % Find the peak in the slow range and use this to align a filtered trace
        % of the high range.
        OUT.CM(:,:,iEpoch,iInterval) = O.CM;
        
        if PLOT_IT
            figure(1)
            clf
            subplot(1,2,1)
            plot(fqs,lpsd,fqs,lspsd,fqs,lpsd_welch)
            hold on
            plot_markers_simple(fqs(pks_ix(ix)));
            set(gca, 'XScale', 'log')
            legend('raw','smooth','welch')
            subplot(1,2,2)
            imagesc(O.low_fq_range,O.high_fq_range, O.CM)
            xlabel('Hz');ylabel('Hz');colorbar;colormap(viridis)
            axis xy
            drawnow
        end
        
        % Find the peak theta frequency for low and high gamma.
        LGIX = O.high_fq_range > OSC.low_gamma_colgin(1) & O.high_fq_range < OSC.low_gamma_colgin(2);
        HGIX = O.high_fq_range > OSC.high_gamma_colgin(1) & O.high_fq_range < OSC.high_gamma_colgin(2);
        
        lgamma_peak_pac_pow = max(max(O.CM(LGIX,:)));
        hgamma_peak_pac_pow = max(max(O.CM(HGIX,:)));
        
        [~, ix] = max(mean(O.CM(LGIX,:)));
        lf_for_lgamma = O.low_fq_range(ix);
        
        [~, ix] = max(mean(O.CM(HGIX,:)));
        lf_for_hgamma = O.low_fq_range(ix);
        
        lg_fqs = O.high_fq_range(LGIX);
        hg_fqs = O.high_fq_range(HGIX);
        
        [~, ix] = max(mean(O.CM(LGIX,:),2));
        hf_for_lgamma = lg_fqs(ix);
        
        [~, ix] = max(mean(O.CM(HGIX,:),2));
        hf_for_hgamma = hg_fqs(ix);
        
        OUT.lgamma_peak_pac_pow(iEpoch,iInterval) = lgamma_peak_pac_pow;
        OUT.hgamma_peak_pac_pow(iEpoch,iInterval) = hgamma_peak_pac_pow;
        OUT.lf_for_lgamma(iEpoch,iInterval) = lf_for_lgamma;
        OUT.lf_for_hgamma(iEpoch,iInterval) = lf_for_hgamma;
        OUT.hf_for_lgamma(iEpoch,iInterval) = hf_for_lgamma;
        OUT.hf_for_hgamma(iEpoch,iInterval) = hf_for_hgamma;
        
        % Now that we know the pacs for the different bands, refilter the
        % low LFP to those targets and compute the peaks of this filtered
        % trace. Use these peaks to align the filtered high and low gamma
        % to see what these bands are doing. Be sure to store running speed
        % so that the relationship between this and running speed can be
        % determined.
        bp = designfilt('bandpassiir','FilterOrder',12, ...
            'HalfPowerFrequency1',lf_for_lgamma - 1.5,'HalfPowerFrequency2',lf_for_lgamma+1.5, ...
            'SampleRate',sFreq,'DesignMethod' ,'butter'); % can't go to nyquist - has to be a little lower.
        %         freqz(bp,[],sFreq)
        LFPrun(:,3) = filtfilt(bp, LFPrun(:,2));
        ZIX = LFPrun(:,2)==0;
        LFPrun(ZIX,3) = 0;
        bp = designfilt('bandpassiir','FilterOrder',12, ...
            'HalfPowerFrequency1',hf_for_lgamma-10,'HalfPowerFrequency2',hf_for_lgamma + 10, ...
            'SampleRate',sFreq,'DesignMethod' ,'butter'); % can't go to nyquist - has to be a little lower.
        LFPrun(:,4) = filtfilt(bp, LFPrun(:,2));
        LFPrun(ZIX,4) = 0;
        
        bp = designfilt('bandpassiir','FilterOrder',12, ...
            'HalfPowerFrequency1',hf_for_hgamma-10,'HalfPowerFrequency2',hf_for_hgamma + 10, ...
            'SampleRate',sFreq,'DesignMethod' ,'butter'); % can't go to nyquist - has to be a little lower.
        
        LFPrun(:,5) = filtfilt(bp, LFPrun(:,2));
        LFPrun(ZIX,5) = 0;
        
        [pk,ix] = findpeaks(LFPrun(:,3));
        ix = ix(pk > median(pk));
        
        if 0
            figure
            plot(LFPrun(:,2))
            hold on
            plot(LFPrun(:,3))
            plot(LFPrun(:,4))
            plot(LFPrun(:,5))
            
            figure
            [M, ixw, x] = PETH_EEG_simple(LFPrun(:,[1 3]), LFPrun(ix,1), sFreq/2, sFreq/2, sFreq, true);
        end
        
        [M, ixw, x] = PETH_EEG_simple(LFPrun(:,[1 3]), LFPrun(ix,1), sFreq/2, sFreq/2, sFreq, false);
        OUT.ETA_x_sec = x;
        OUT.mean_ETA(:,iEpoch,iInterval) = mean(M);
        [M, ixw, x] = PETH_EEG_simple([LFPrun(:,1) abs(hilbert(LFPrun(:,4)))], LFPrun(ix,1), sFreq/2, sFreq/2, sFreq, false);
        OUT.mean_ETA_lgam(:,iEpoch,iInterval) = mean(abs(M));
        [M, ixw, x] = PETH_EEG_simple([LFPrun(:,1) abs(hilbert(LFPrun(:,5)))], LFPrun(ix,1), sFreq/2, sFreq/2, sFreq, false);
        OUT.mean_ETA_hgam(:,iEpoch,iInterval) = mean(abs(M));
        
        % Determine phase alignment to theta.
        [M, ixw, x] = PETH_EEG_simple([LFPrun(:,1) angle(hilbert(LFPrun(:,4)))], LFPrun(ix,1), sFreq/2, sFreq/2, sFreq, false);
        OUT.mean_ETA_ph_lgam(:,iEpoch,iInterval) = nanmean(M);
        [M, ixw, x] = PETH_EEG_simple([LFPrun(:,1) angle(hilbert(LFPrun(:,5)))], LFPrun(ix,1), sFreq/2, sFreq/2, sFreq, false);
        OUT.mean_ETA_ph_hgam(:,iEpoch,iInterval) = nanmean(M);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT.PSD_fqs = fqs;
OUT.low_fq_range = O.low_fq_range;
OUT.high_fq_range = O.high_fq_range;
OUT.method = O.method;
OUT.theta_fname = theta_fname;
OUT.CFC_method = CFC_method;
OUT.speed_intervals = speed_intervals;


function OUT = LD_Create_trial_info_file_LC(ses)
% Dissects the behavior to identify all trials when the rat was or was not
% zapped and also indicate the direction of travel.
%
% TODO: Needs to have the start and end time of each trial and for each
% trial, heading for trial (directo of travel) whether the trial zapped, where zapped, and when zapped.
%
%
% NEED THE CONTROL ZONE INFO AS WELL.
% Cowen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ==0
    load Session_Info.mat
else
    save('Session_Info','ses')
end
PLOT_IT = false;
GP = LD_Globals_LC;
OUT.Cutoff_point_for_stim_zones_theta = 200;
OUT.mfilename = mfilename;
OUT.pwd = pwd;
OUT.ses = ses;
OUT.abort_session = false;
end_of_trial_points = [16 342];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stimulation times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[EVT, STIMTIMES, ALLTIMES ]= LD_Load_Events;
OUT.Stim_starts.Maze1 = STIMTIMES{2}(:,1);
OUT.Stim_starts.Maze2 = STIMTIMES{5}(:,1);
%
% OUT.All_Stim_times_usec.Maze2 = EVT(2).ts_usec;
%
% [OUT.Stim_starts.Maze1, durations, n_per_burst, burst_id] = Burst_detector(OUT.All_Stim_times_usec.Maze1, 1e6, 2);
% [OUT.Stim_starts.Maze2, durations, n_per_burst, burst_id] = Burst_detector(OUT.All_Stim_times_usec.Maze2, 1e6, 2);

% OUT.Stim_starts.Maze1 = unique(OUT.Stim_starts.Maze1,'rows');
% OUT.Stim_starts.Maze2 = unique(OUT.Stim_starts.Maze2,'rows');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Epoch times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EP = load('epochs.mat'); % NOTE: these are in 0.1msec units.
OUT.All_Stim_times_usec.Maze1 = Restrict(ALLTIMES,EP.epochs.Maze1*100);
OUT.All_Stim_times_usec.Maze2 = Restrict(ALLTIMES,EP.epochs.Maze2*100);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load position and event data and epoch.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, POS_maze,Pos_labels] = LD_Load_Position_LC;
epoch_names = {'Maze1' 'Maze2' };
for iEpoch = 1:2
    
    ename = epoch_names{iEpoch};
    EP_st_ed_usec(1) = EP.epochs.(ename)(1)*100 + 5e6; % add some time for slop in getting things started.
    EP_st_ed_usec(2) = EP.epochs.(ename)(2)*100 - 10e6; % add some time for slop in getting things finished.
    
    
    T_stim_usec = Restrict(unique(OUT.Stim_starts.(ename)),EP_st_ed_usec);
    if isempty(T_stim_usec)
        OUT.abort_session = true;
        save(fullfile(GP.Processed_data_dir,[ ses.animal '_' ses.name '_TrialInfo']),'-struct','OUT')
        msgbox(sprintf('Aborted in %s \n Due to screwed up stim times',pwd))
        return
    end
    [SF] = ScatterFields_cowen(POS_maze{iEpoch},T_stim_usec);
    SF = unique(SF,'rows');
    
    IXlt = SF(:,4) < OUT.Cutoff_point_for_stim_zones_theta;
    IXgt = SF(:,4) > OUT.Cutoff_point_for_stim_zones_theta;
    OUT.Stimzone_theta.(ename) = zeros(2,3);
    OUT.Stimzone_theta.(ename)(1,:) = [mean(SF(IXlt,4))  min(SF(IXlt,4)) max(SF(IXlt,4))];
    OUT.Stimzone_theta.(ename)(2,:) = [mean(SF(IXgt,4))  min(SF(IXgt,4)) max(SF(IXgt,4))];
    df = OUT.Stimzone_theta.(ename)(2,1) - OUT.Stimzone_theta.(ename)(1,1);
    cz = OUT.Stimzone_theta.(ename)(1,1) + df/2;
    OUT.Stimzone_theta.(ename)(3,:) = [cz  cz - 15 cz + 15]; % need to remember that the 3rd is the control zone.
    
    v = sortrows(POS_maze{iEpoch}(:,[4 2 3]),1);
    ixc = find(v(:,1) > cz,1,'first');
    ixl = find(v(:,1) > cz-GP.Control_Zone_Size_deg/2,1,'first');
    ixu = find(v(:,1) > (cz+GP.Control_Zone_Size_deg/2),1,'first');
    
    OUT.Stimzone_x.(ename) = zeros(2,3);
    OUT.Stimzone_x.(ename)(1,:) = [mean(SF(IXlt,2))  min(SF(IXlt,2)) max(SF(IXlt,2))];
    OUT.Stimzone_x.(ename)(2,:) = [mean(SF(IXgt,2))  min(SF(IXgt,2)) max(SF(IXgt,2))];
    OUT.Stimzone_x.(ename)(3,:) = [v(ixc,2)  v(ixl,2) v(ixu,2) ]; % need to remember that the 3rd is the control zone.
    
    OUT.Stimzone_y.(ename) = zeros(2,3);
    OUT.Stimzone_y.(ename)(1,:) = [mean(SF(IXlt,3))  min(SF(IXlt,3)) max(SF(IXlt,3))];
    OUT.Stimzone_y.(ename)(2,:) = [mean(SF(IXgt,3))  min(SF(IXgt,3)) max(SF(IXgt,3))];
    OUT.Stimzone_y.(ename)(3,:) = [v(ixc,3)  v(ixl,3) v(ixu,3)]; % need to remember that the 3rd is the control zone.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % What is a trial? A trial is a journey in one or the other direction.
    % On that journey you need to know 1) direction of travel, 2) whether
    % he got zapped or if he did not, if he was still entering a stim zone.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find all points in which there was a position before and after Stimzone_theta
    OUT.TI.(ename).ColumnLabels = {Pos_labels 'zone' 'heading' 'stimtime'};
    OUT.TI.(ename).ZoneCrossTimeUsecAndPos = [];
    
%     d = diff(POS_maze{iEpoch}(:,4))
    smpos = sgolayfilt( movmedian( POS_maze{iEpoch}(:,4),50),3,231);
    smpos = POS_maze{iEpoch}(:,4);
    bix = smpos > 130 & smpos < 220;
%     smpos(bix) = ;
    all_ix = 1:length(smpos);
    smpos(bix) = interp1(all_ix(~bix),smpos(~bix),all_ix(bix),'linear');
    
%     smpos = movmean( movmedian( POS_maze{iEpoch}(:,4),450),400);

    for iZone = 1:3
        THETAn = rad2deg(circ_dist(deg2rad(smpos),deg2rad(OUT.Stimzone_theta.(ename)(iZone,1))));
%         THETAn = smpos - 185;
        %         THETAns = POS_maze{iEpoch}(:,4) - OUT.Stimzone_theta.(ename)(iZone,1); % subtract the
        [~,~,ZeroUpIdx, ZeroDownIdx] = Find_peaks_troughs_zeros(THETAn); % Zero crossings.
        diff(ZeroUpIdx)
        % zeroup are the times when the rat enters the field from the
        % left, where the zerodown are when he comes through this field from
        % the right.
        %         figure;plot(THETAn);hold on;plot_markers_simple(ZeroUpIdx);plot_markers_simple(ZeroDownIdx,[],[],'r');
        M1 = [ POS_maze{iEpoch}(ZeroUpIdx,:)   ones(length(ZeroUpIdx),1)*iZone  ones(length(ZeroUpIdx),1)     ];
        M2 = [ POS_maze{iEpoch}(ZeroDownIdx,:) ones(length(ZeroDownIdx),1)*iZone  -1*ones(length(ZeroDownIdx),1) ];
        OUT.TI.(ename).ZoneCrossTimeUsecAndPos = [OUT.TI.(ename).ZoneCrossTimeUsecAndPos;M1;M2];
    end
    OUT.TI.(ename).ZoneCrossTimeUsecAndPos = sortrows(OUT.TI.(ename).ZoneCrossTimeUsecAndPos);
    % Some datasets are just fucked up....
    
    %%THIS DOWN HERE I WILL WANT TO KEEP IN BUT IT REQUIRES THAT SES BE A
    %%STRUCTURE AND I"M NOT SURE WHERE THAT GETS PUT IN, ITS THE INPUT
    %%VARIABLE TO THE FUNCTION BUT WHERE DOES IT COME FROM
%     if strcmpi(ses.animal,'8417') && strcmpi(ses.name,'Day21') && strcmpi(ename,'Maze1')
%         %Shit
%         GIX =   OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) > 3.7343e+09 & OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) < 5.0081e+09 ;
%         OUT.TI.(ename).ZoneCrossTimeUsecAndPos  = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(GIX,:);
%     end
%     if strcmpi(ses.animal,'8417') && strcmpi(ses.name,'Day21') && strcmpi(ename,'Maze2')
%         %Shit
%         GIX =   OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) < 1.0186e+10 ;
%         OUT.TI.(ename).ZoneCrossTimeUsecAndPos  = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(GIX,:);
%     end
%     if strcmpi(ses.animal,'8419') && strcmpi(ses.name,'Day12') && strcmpi(ename,'Maze1')
%         %Shit
%         GIX =   OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) >  4.4802e+09  & OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) < 5.1825e+09 ;
%         OUT.TI.(ename).ZoneCrossTimeUsecAndPos  = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(GIX,:);
%     end
%     
%     if strcmpi(ses.animal,'8419') && strcmpi(ses.name,'Day20') && strcmpi(ename,'Maze1')
%         %Shit
%         GIX =   OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) >  3.3199e+09  & OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) < 4.002e+09 ;
%         OUT.TI.(ename).ZoneCrossTimeUsecAndPos  = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(GIX,:);
%     end
%     if strcmpi(ses.animal,'8419') && strcmpi(ses.name,'Day20') && strcmpi(ename,'Maze2')
%         %Shit
%         GIX =   OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) >  8.0218e+09   & OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) < 8.6331e+09  ;
%         OUT.TI.(ename).ZoneCrossTimeUsecAndPos  = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(GIX,:);
%     end
%     
%     if strcmpi(ses.animal,'8564') && strcmpi(ses.name,'Day01') && strcmpi(ename,'Maze2')
%         %Shit
%         GIX =   OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) >   9.9935e+09    & OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) < 1.1358e+10  ;
%         OUT.TI.(ename).ZoneCrossTimeUsecAndPos  = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(GIX,:);
%     end
%     
%     if strcmpi(ses.animal,'8886') && strcmpi(ses.name,'Day26') && strcmpi(ename,'Maze1')
%         %Shit
%         GIX =   OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) >   4.2455e+09    & OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) < 5.8243e+09  ;
%         OUT.TI.(ename).ZoneCrossTimeUsecAndPos  = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(GIX,:);
%     end
%     if strcmpi(ses.animal,'8981') && strcmpi(ses.name,'Day09') && strcmpi(ename,'Maze2')
%         %Shit
%         BIX =   OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) >   1.2983e+10    & OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1) < 1.2987e+10 ;
%         OUT.TI.(ename).ZoneCrossTimeUsecAndPos(BIX,:) = [];
%     end

    
%     figure(1);clf
%     plot( OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1), OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,4))
%     hold on
%     plot(POS_maze{iEpoch}(:,1),POS_maze{iEpoch}(:,4))
    % Deterimine each run's start and end time and it's direction....
    %     t_head = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,[1 7]);
    %     t_theta = POS_maze{iEpoch}(:,[1 4]);
    %     badix = diff(t_head(2:end,2));
    %     %     t_start_end(count,:)
    %
    %     d = diff(diff(t_head(:,2)));
    %     badix = find(abs(d) > 2) + 2;
    %     t_head(badix,:) = [];
    %     OUT.TI.(ename).ZoneCrossTimeUsecAndPos(badix,:) = [];
    %     d = diff([t_head(1,2)*10; t_head(:,2)]);
    %     st_ix = find(d > 1);
    %     ed_ix = find(d < 0);
    %     if ed_ix(1) < st_ix(1)
    %         tmp = st_ix;
    %         st_ix = ed_ix;
    %         ed_ix = tmp;
    %     end
    %
    %
    
    
    
    
    % Now that we have determined the times when the rat passed through
    % each stimzone, determine if a stim actually happenened there. How?
    % Look for the closest stim time. If it happened within 3s of being
    % around the stim zone, then call it a stim trial. Also record the stim
    % time.
    for iR = 1:length(T_stim_usec)
        m = abs(T_stim_usec(iR) - OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,1));
        [min_time, ix] = min(m);
        %         min_time/1e6
        %         figure(1)
        %         hist( min_time/1e6,40)
        %         pause
        %
        if min_time/1e6 < 3
            OUT.TI.(ename).ZoneCrossTimeUsecAndPos(ix,8) = 1;
            OUT.TI.(ename).ZoneCrossTimeUsecAndPos(ix,9) = T_stim_usec(iR);
        else
            OUT.TI.(ename).ZoneCrossTimeUsecAndPos(iR,8) = 0;
            OUT.TI.(ename).ZoneCrossTimeUsecAndPos(iR,9) = nan;
        end
    end
    %
    IX =  OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,7)==1  & ...
        OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,8)==1  & ...
        OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,9)==0;
    % but this does not get at heading. Nonstim should only be in the same
    % heading.
    OUT.NonStimtimes_matched.(ename).Zone1 =  OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX,1);
    
    IX =  OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,7)==2  & ...
        OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,8)==-1  & ...
        OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,9)==0;
    
    OUT.NonStimtimes_matched.(ename).Zone2 =  OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX,1);
    
    
    OUT.TI.(ename).ZoneCrossTimeUsecAndPos = sortrows(OUT.TI.(ename).ZoneCrossTimeUsecAndPos);
    Z_Head_Stim = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,6:8);
    U = unique(Z_Head_Stim,'rows');
    U = U(U(:,3)==1,:);
    IX1 = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,6) == U(1,1) & OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,7) == U(1,2);
    IX2 = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,6) == U(2,1) & OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,7) == U(2,2);
    
    OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,10) = zeros(Rows(OUT.TI.(ename).ZoneCrossTimeUsecAndPos),1)*nan;
    OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX1,10)  = 1;
    OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX2,10)  = 1;
    % Find trial starts and ends
    z = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,7);
    IX1 = z(1:end-2,:) == 1 & z(2:end-1,:) == 3 & z(3:end,:) == 2;
    IX2 = z(1:end-2,:) == 2 & z(2:end-1,:) == 3 & z(3:end,:) == 1;
    TR123_1st_zcross = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX1,1);
    TR321_1st_zcross = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX2,1);
    TR123_last_zcross = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(find(IX1)+2,1);
    TR321_last_zcross = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(find(IX2)+2,1);
    % now move this forward.
    T1 = [TR123_1st_zcross TR123_last_zcross ones(size(TR123_1st_zcross))];
    T2 = [TR321_1st_zcross TR321_last_zcross 2*ones(size(TR321_1st_zcross))];
    TRSTED = sortrows([T1;T2]);
    TRSTED = [[1:Rows(TRSTED)]' TRSTED];
    % For each trial, determine if the rat was stimmed and ID the
    % trajectory (movement direction.
    if iEpoch == 1
        epnum = 2;
    elseif iEpoch == 2
        epnum = 5;
    end
    
    all_nonstim_usec{1} = OUT.NonStimtimes_matched.(ename).Zone1;
    all_nonstim_usec{2} = OUT.NonStimtimes_matched.(ename).Zone2;
    
    MULTISTIM = false(size(TRSTED(:,1)));
    MULTINONSTIM = false(size(TRSTED(:,1)));
    for iR = 1:Rows(TRSTED)
        % Determine if stimmed and store stim time.
        ix = find(STIMTIMES{epnum}(:,1)> (TRSTED(iR,2)-1e6) & STIMTIMES{epnum}(:,1)< (TRSTED(iR,3) + 1e6));
        %         tr_dur_sec = (TRSTED(iR,3)-TRSTED(iR,2))/1e6;
        if isempty(ix)
            TRSTED(iR,5:6) = [nan nan];
        else
            if length(ix) > 1
                MULTISTIM(iR) = true;
            end
            TRSTED(iR,5:6) = STIMTIMES{epnum}(ix(1),1:2);
        end
        % Determine if stimmed and store stim time.
        ZONE = TRSTED(iR,4);
        
        ix = find(all_nonstim_usec{ZONE}> (TRSTED(iR,2)-1e6) & all_nonstim_usec{ZONE} < (TRSTED(iR,3) + 1e6));
        %         tr_dur_sec = (TRSTED(iR,3)-TRSTED(iR,2))/1e6;
        if isempty(ix)
            TRSTED(iR,7:8) = [nan nan];
        else
            if length(ix) > 1
                MULTINONSTIM(iR) = true;
            end
            TRSTED(iR,7:8) = [all_nonstim_usec{ZONE}(ix(1)) all_nonstim_usec{ZONE}(ix(1)) + 90*1000];
        end
        
    end
    OUT.TI.(ename).TRL.TrlNum = TRSTED(:,1);
    OUT.TI.(ename).TRL.TrlTypeAndZone = TRSTED(:,4);
    OUT.TI.(ename).TRL.TimeEnterZone1or3 = TRSTED(:,2:3);
    OUT.TI.(ename).TRL.Stimmed = ~isnan(TRSTED(:,5));
    OUT.TI.(ename).TRL.StimStEd = TRSTED(:,5:6);
    OUT.TI.(ename).TRL.NonStimStEd = TRSTED(:,7:8);
    OUT.TI.(ename).TRL.TrDurationSec =(TRSTED(:,3)-TRSTED(:,2))/1e6;
    OUT.TI.(ename).TRL.MultipleStimOnTrial = MULTISTIM; % something probably went wrong.
    OUT.TI.(ename).TRL.BadTrial = MULTISTIM | OUT.TI.(ename).TRL.TrDurationSec > 40 | OUT.TI.(ename).TRL.TrDurationSec < .5; % something probably went wrong.
    %%%%%%%%%%%%%%%%%%%%%
    PM = POS_maze{iEpoch}(:,[1 4]);
    
    OUT.TI.(ename).TRL.StimLoc = nan(size(OUT.TI.(ename).TRL.StimStEd(:,1)));
    IX = ~isnan(OUT.TI.(ename).TRL.StimStEd(:,1));
    T = OUT.TI.(ename).TRL.StimStEd(IX,1);
    ix = binsearch_vector(PM(:,1),T);
    OUT.TI.(ename).TRL.StimLoc(IX) = PM(ix,2);
    
    IX = ~isnan(OUT.TI.(ename).TRL.NonStimStEd(:,1));
    T = OUT.TI.(ename).TRL.NonStimStEd(IX,1);
    ix = binsearch_vector(PM(:,1),T);
    OUT.TI.(ename).TRL.StimLoc(IX) = PM(ix,2);
    
    
    % Get a good estimate of the true start and end time.
    for ii = 1:length(OUT.TI.(ename).TRL.TrlNum)
        st = OUT.TI.(ename).TRL.TimeEnterZone1or3(ii);
        ix = binsearch(PM(:,1),st);
        if OUT.TI.(ename).TRL.TrlTypeAndZone(ii) == 2
            % Go back in time and find the first crossing to 300.
            lastcross = find(PM(1:ix,2) > GP.RewardZoneLinPos(2),1,'last');
            
            lastcross2 = find(PM(ix:end,2) < GP.RewardZoneLinPos(1),1,'first');
            lastcross2 = lastcross2 + ix-1;
        elseif OUT.TI.(ename).TRL.TrlTypeAndZone(ii) == 1
            lastcross = find(PM(1:ix,2) < GP.RewardZoneLinPos(1),1,'last');
            lastcross2 = find(PM(ix:end,2) > GP.RewardZoneLinPos(2),1,'first');
            lastcross2 = lastcross2 + ix-1;
        else
            error('wtf')
        end
        if isempty(lastcross)
            lastcross = 1;
        end
        if isempty(lastcross2)
            lastcross2 = Rows(PM);
        end
        OUT.TI.(ename).TRL.TrlStartEnd(ii,1) = PM(lastcross,1);
        OUT.TI.(ename).TRL.TrlStartEnd(ii,2) = PM(lastcross2,1);
        
    end
    IX = diff(OUT.TI.(ename).TRL.TrlStartEnd(:,1))<=0;
    if any(IX)
        if find(IX) < 4
            OUT.TI.(ename).TRL.BadTrial(find(IX)) = true;
        else
            OUT.TI.(ename).TRL.BadTrial(find(IX)+1) = true;
        end
        disp('asdfp')
    end
    
     GIX = ~OUT.TI.(ename).TRL.BadTrial;
     
     % CLean shit up and get rid of all bad trials from all recs.
     f = fieldnames(OUT.TI.(ename).TRL);
     for iF = 1:length(f)
         OUT.TI.(ename).TRL.(f{iF}) = OUT.TI.(ename).TRL.(f{iF})(GIX,:);
     end
     
    IX = diff(OUT.TI.(ename).TRL.TrlStartEnd(~ OUT.TI.(ename).TRL.BadTrial,1))<=0;

    % For each trial, le'ts identify the start and end of the food zone
    % period - but let's guistimate the first trial (1s) as he just starts running
    % then.
%     food_times = [OUT.TI.(ename).TRL.TrlStartEnd(1:end-1,2)  OUT.TI.(ename).TRL.TrlStartEnd(2:end,1)]
%     food_times = [OUT.TI.(ename).TRL.TrlStartEnd(1,1)-1e6 OUT.TI.(ename).TRL.TrlStartEnd(1,1);food_times];
%     long_ix = (food_times(:,2) - food_times(:,1))/1e6 > 10;
%     if any(long_ix)
%         food_times(long_ix,2) = food_times(long_ix,1)+5e6;
%         error('aahsahf too long')
%     end
%     OUT.TI.(ename).TRL.EatingFoodStartEnd = food_times;
    
    % Compute rates per trial for each food zone.
    OUT.TI.(ename).TRL.EatAtStartOfRun(:,2) = OUT.TI.(ename).TRL.TrlStartEnd(:,1);
    OUT.TI.(ename).TRL.EatAtStartOfRun(1,1) = OUT.TI.(ename).TRL.TrlStartEnd(1,1)-1e6;
    OUT.TI.(ename).TRL.EatAtStartOfRun(2:end,1) = OUT.TI.(ename).TRL.TrlStartEnd(1:(end-1),2);
    
    OUT.TI.(ename).TRL.EatAtEndOfRun(:,1) = OUT.TI.(ename).TRL.TrlStartEnd(:,2);
    OUT.TI.(ename).TRL.EatAtEndOfRun(:,2) = [OUT.TI.(ename).TRL.TrlStartEnd(2:end,1); OUT.TI.(ename).TRL.TrlStartEnd(end,2) + 1e6];
   
    if any(diff(OUT.TI.(ename).TRL.TrlStartEnd(:,1))<0)
        error('negative times TrlStartEnd')
    end
    
    if any(diff(OUT.TI.(ename).TRL.EatAtStartOfRun(:,1))<0)
        error('negative times EatAtStartOfRun')
    end
 
    if any(OUT.TI.(ename).TRL.EatAtEndOfRun(:,2) - OUT.TI.(ename).TRL.EatAtEndOfRun(:,1)<0)
        error('negative EatAtEndOfRun INTERVAL')
    end
 
    
    
    %     OUT.TI.(ename).TRL.TrialStEd = TRSTED(:,2:3);
    
    if PLOT_IT
        figure(iEpoch)
        clf
        subplot(2,2,1)
        plot(POS_maze{iEpoch}(:,2),POS_maze{iEpoch}(:,3),'c',SF(:,2),SF(:,3),'ro')
        a = axis; hold on;
        plot(  OUT.Stimzone_x.(ename)(3,1),OUT.Stimzone_y.(ename)(3,1),'r*')
        plot(  OUT.Stimzone_x.(ename)(3,2),OUT.Stimzone_y.(ename)(3,2),'rx')
        plot(  OUT.Stimzone_x.(ename)(3,3),OUT.Stimzone_y.(ename)(3,3),'rx')
        
        
        subplot(2,2,2)
        plot(POS_maze{iEpoch}(:,4),POS_maze{iEpoch}(:,5),'c',SF(:,4),SF(:,5),'rx')
        hold on
        IX1 = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,7) == 1;
        IX2 = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,7) == -1;
        plot(OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX1,4),OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX1,5),'>g');
        plot(OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX2,4),OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX2,5),'<m');
        IX1 = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,8) == 1;
        plot(OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX1,4),OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX1,5),'r.');
        IX1 = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,7) == 1 & OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,6) == 3;
        IX2 = OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,7) == -1 & OUT.TI.(ename).ZoneCrossTimeUsecAndPos(:,6) == 3;
        plot(OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX1,4),OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX1,5),'>k');
        plot(OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX2,4),OUT.TI.(ename).ZoneCrossTimeUsecAndPos(IX2,5),'<r');
        
        
        a = axis; hold on;
        plot(  OUT.Stimzone_theta.(ename)(1,[1 1]),a(3:4),'r:')
        plot(  OUT.Stimzone_theta.(ename)(2,[1 1]),a(3:4),'g:')
        plot(  OUT.Stimzone_theta.(ename)(3,[1 1]),a(3:4),'m:')
        plot(  OUT.Stimzone_theta.(ename)(3,[2 2]),a(3:4),'m')
        plot(  OUT.Stimzone_theta.(ename)(3,[3 3]),a(3:4),'m')
        subplot(2,2,3:4)
        % Plot trial start and end times...
        plot(POS_maze{iEpoch}(:,1)/60e6,POS_maze{iEpoch}(:,4),'k');
        axis tight
        hold on
        %         plot_markers_simple([TR123_1st_zcross TR123_last_zcross]/60e6)
        %         plot_markers_simple([TR321_1st_zcross TR321_last_zcross]/60e6,[],[],[.3 .9 .0])
        %         plot_markers_simple([TR321_1st_zcross TR321_last_zcross]/60e6,[],[],[.3 .9 .0])
        IX1 = OUT.TI.(ename).TRL.TrlTypeAndZone == 1;
        IX2 = OUT.TI.(ename).TRL.TrlTypeAndZone == 2;
        R1 = Restrict([POS_maze{iEpoch}(:,1)/60e6,POS_maze{iEpoch}(:,4)],OUT.TI.(ename).TRL.TrlStartEnd(IX1,:)/60e6);
        R2 = Restrict([POS_maze{iEpoch}(:,1)/60e6,POS_maze{iEpoch}(:,4)],OUT.TI.(ename).TRL.TrlStartEnd(IX2,:)/60e6);
        plot(R1(:,1),R1(:,2),'g.')
        plot(R2(:,1),R2(:,2),'r.')
        plot(OUT.TI.(ename).TRL.EatAtStartOfRun(IX1,1)/60e6, ones(size(OUT.TI.(ename).TRL.EatAtStartOfRun(IX1,1))),'g>')
        plot(OUT.TI.(ename).TRL.EatAtStartOfRun(IX1,2)/60e6, ones(size(OUT.TI.(ename).TRL.EatAtStartOfRun(IX1,2))),'r<')
        plot(OUT.TI.(ename).TRL.EatAtStartOfRun(IX2,1)/60e6, 300*ones(size(OUT.TI.(ename).TRL.EatAtStartOfRun(IX2,1))),'m>')
        plot(OUT.TI.(ename).TRL.EatAtStartOfRun(IX2,2)/60e6, 300*ones(size(OUT.TI.(ename).TRL.EatAtStartOfRun(IX2,2))),'c<')
        %         plot_markers_simple(OUT.TI.(ename).TRL.TrlStartEnd(IX2,:)/60e6,[],[],[.9 .0 .1])
%         saveas(gcf,fullfile('C:\Temp\',sprintf('Rat%sSes%sEp%dZONES.png',ses.animal,ses.name,iEpoch)))
        saveas(gcf,sprintf('Rat%sSes%sEp%dZONES.png',ses.animal,ses.name,iEpoch)) % save with the data - good to have when playing with individual sets.
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(fullfile(GP.Processed_data_dir,[ ses.animal '_' ses.name '_TrialInfo']),'-struct','OUT')


% function Q1_theta_gamma_coupling_Ana
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clearvars
mfile = 'Q1_theta_gamma_coupling_Ana';
GP = KT_Globals;
PLOT_IT = true;
ses_to_ana = 'Q1_theta_gamma_coupling';
adir = fullfile(GP.Analysis_dir,'Lesley_Data',ses_to_ana);
d = dir(fullfile(adir,'Dset*.mat'));
load(fullfile(adir,d(1).name),'Dset'); % load the first dataset to be sure things work.
%%
ses_cnt = 1;
CM = [];
group = [];
day = [];
animal = [];
peak_fq = [];
peak_pow = [];
lf_for_lgamma = [];
ses = {};
for iF = 1:length(d)
    load(fullfile(adir,d(iF).name),'Dset');
    if Dset.aborted
        continue
    end
    %
    [p,n] = fileparts(Dset.theta_fname);
    ses{ses_cnt} = [num2str(iF) ',' Dset.ses.animal ',' Dset.ses.name ',' Dset.ses.group ',' n];
    fprintf('%s\n',ses{ses_cnt})
    %
    group{ses_cnt} = Dset.ses.group;
    day(ses_cnt) = str2double( Dset.ses.name(4:end) );
    animal(ses_cnt) = str2double( Dset.ses.animal );
    % x y epoch param.    OUT.peak_fqs(iEpoch,iInterval,ix)
    %        OUT.CM(:,:,iEpoch,iInterval) = O.CM;
    % OUT.hf_for_hgamma(iEpoch,iInterval)
    for iInterval = 1:Rows(Dset.speed_intervals)
        if size(Dset.CM,4) < iInterval
            CM{iInterval}(:,:,ses_cnt) = nan;
            peak_fq{iInterval}(ses_cnt,:)  = nan;
            peak_pow{iInterval}(ses_cnt,:) = nan;
            lf_for_lgamma{iInterval}(ses_cnt) = nan;
            continue
        end
        CM{iInterval}(:,:,ses_cnt) = squeeze(nanmean(Dset.CM(:,:,:,iInterval),3));
        peak_fq{iInterval}(ses_cnt,:) = mean(Dset.peak_fqs(:,iInterval,:),1);
        peak_pow{iInterval}(ses_cnt,:) = mean(Dset.peak_pow(:,iInterval,:),1);
        lf_for_lgamma{iInterval}(ses_cnt) = mean(Dset.lf_for_lgamma(:,iInterval),1);
    end
    ses_cnt = ses_cnt + 1;
end
group = categorical(group);
%%
% func = @(x) mean(x,3);
% let's plot everything
spd_labs = {'stop' 'slow' 'med' 'fast'};
for iSes = 1:length(group)
    clf
    h = [];
    for iI = 1:Rows(Dset.speed_intervals)
        h(iI) = subplot(1,Rows(Dset.speed_intervals),iI);
        %         imagesc(Dset.low_fq_range,Dset.high_fq_range, log(CM{iI}(:,:,iSes)))
        imagesc(Dset.low_fq_range,Dset.high_fq_range, CM{iI}(:,:,iSes))
        axis xy
        colorbar; colormap(viridis)
        xlabel('Hz');ylabel('Hz')
        hold on
        plot(peak_fq{iI}(iSes,1),Dset.high_fq_range(1),'w*')
        plot(lf_for_lgamma{iI}(iSes),Dset.high_fq_range(1),'m^')
        if iI == 1
            title([ses{iSes} ' ' spd_labs{iI}])
        else
            title([spd_labs{iI}])
        end
    end
    equalize_color_axes(h);
    pause
end
%%
YIX = group == 'yng';
OIX = group == 'old';

figure
h = [];
INTERVAL = 2;
h(1) = subplot(2,2,1);
imagesc(Dset.low_fq_range,Dset.high_fq_range, (nanmean(CM{INTERVAL}(:,:,YIX),3)))
title('CFC Young')
axis xy
colorbar
xlabel('Hz');ylabel('Hz')

h(2) = subplot(2,2,2);
imagesc(Dset.low_fq_range,Dset.high_fq_range, (nanmean(CM{INTERVAL}(:,:,OIX),3)))
axis xy
colorbar
xlabel('Hz');ylabel('Hz')
title('CFC Old')
equalize_color_axes(h);
INTERVAL = 3;

h(3) = subplot(2,2,3);
imagesc(Dset.low_fq_range,Dset.high_fq_range, (nanmean(CM{INTERVAL}(:,:,YIX),3)))
title('CFC Young')
axis xy
colorbar
xlabel('Hz');ylabel('Hz')

h(4) = subplot(2,2,4);
imagesc(Dset.low_fq_range,Dset.high_fq_range, (nanmean(CM{INTERVAL}(:,:,OIX),3)))
axis xy
colorbar
xlabel('Hz');ylabel('Hz')
title('CFC Old')
equalize_color_axes(h);

%%%%%%%%%%%%%%%
%% What was the mean fq for each group?
%%%%%%%%%%%%%%%
INTERVAL = 3;
peak_fq1(peak_fq{INTERVAL} > 15 | peak_fq{INTERVAL} < 4) = nan; % wacko values of peak 'theta';
pkfq_grp = grpstats(peak_fq{INTERVAL},group);
unique(group)
%%%%%%%%%%%%%%%
%% Plot the peak frequencies for CFC for each session
%%%%%%%%%%%%%%%
for ii = 1:length(peak_fq)
    peak_fq{ii}(peak_fq{ii}(:,1) > 20,1) = nan;
end
[~,psw] = swtest(peak_fq{2}(group == 'yng',1));
[~,psw_DOUBLECHECK] = swtest(randn(100,1));
[~,psw_DOUBLECHECK2] = swtest(rand(100,1));
figure
h = [];
h(1) = subplot(2,2,1)
boxplot(peak_fq{2}(:,1),group)
ylabel('peak Hz');

h(2) = subplot(2,2,2)
boxplot(peak_fq{3}(:,1),group)
ylabel('peak Hz');

equalize_y_axes(h);

subplot(2,2,3)
histogram_cowen({peak_fq{2}(YIX,1) peak_fq{2}(OIX,1)},.25)
legend('yng','old')
subplot(2,2,4)
histogram_cowen({peak_fq{3}(YIX,1) peak_fq{3}(OIX,1)},.25)
legend('yng','old')
%%%%%%%%%%%%%%%%
%% Do the same but for each animal. Do we use parametric or non-parametric measures?
%%%%%%%%%%%%%%%%
pkfq_animal = grpstats(peak_fq1,animal,'nanmean');

[an,ix] = unique(animal);
grp = group(ix);
[~,psw] = swtest(pkfq_animal(grp == 'yng')); % young pass
[~,psw] = swtest(pkfq_animal(grp == 'old')); % old don't so we should probably default to non-parametric.

figure
boxplot(pkfq_animal,grp)
ylabel('peak Hz');
IXO = grp == 'old';
IXY = grp == 'yng';
% [~,p] = ttest2(pkfq_animal(IXO),pkfq_animal(IXY),'Vartype','unequal');
[p] = ranksum(pkfq_animal(IXO),pkfq_animal(IXY));
title(sprintf('p = %f',p))




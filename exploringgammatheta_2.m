GP = LD_Globals_LC();
load('Session_Info.mat')

fs=1000; %pick this - for how you want the fs to read out from ReadCR_cowen
%%
%want the theta file
theta_fname = ses.other.ThetaEEG_FileName; %getting rid of curly brackets {1}, for some reason i dont need it?

[LFP, sFreq, saturation_times, head] = ReadCR_cowen(theta_fname,[],[],fs,1);
%the actual sample rate is 1890.4, but can input something different so
%I'll make it come out in 1000
%%
%plot with turned into hours
figure
plot(LFP(:,1)/3600e6,LFP(:,2))

%Notes that Peter hasa couple Gamma analysis functions....
% % %Peter's "Gamma Analysis" - can choose Morlet or Hilbert
% % [lfp,phs,env,pow, fh1, fh2] = GammaEpisodeAnalysis_Morlet_v1p1(LFP, sFreq, tit, par, mazeNum)
% % [lfp,phs,env, fh1, fh2] = GammaAnalysis_v0p1(LFP,sFreq)
%%
LowGamma = designfilt('bandpassiir','FilterOrder',12, ...
         'HalfPowerFrequency1',25, 'HalfPowerFrequency2',50, ...
         'SampleRate',sFreq,'designmethod', 'butter');
                    as_lg=filtfilt(LowGamma,LFP(:,2));
                    lowpwr=abs(hilbert(as_lg));

HighGamma = designfilt('bandpassiir','FilterOrder',12, ...
         'HalfPowerFrequency1',65, 'HalfPowerFrequency2',100, ...
         'SampleRate',sFreq,'designmethod', 'butter');
                    as_hg= filtfilt(HighGamma,LFP(:,2));
                    highpwr=abs(hilbert(as_hg)); 
                    
Theta = designfilt('bandpassiir','FilterOrder',12, ...
         'HalfPowerFrequency1',6, 'HalfPowerFrequency2',10, ...
         'SampleRate',sFreq,'designmethod', 'butter');
                    as_theta= filtfilt(Theta,LFP(:,2));
                    thetapwr=abs(hilbert(as_theta));                    
                    thetaphase=angle(hilbert(as_theta));
%%
figure;
plot(LFP(:,1),lowpwr)
hold on;
plot(LFP(:,1),thetaphase)

phasevals=linspace(min(thetaphase),max(thetaphase),360);
% phase_rounded = interp1(phasevals,phasevals,thetaphase,'nearest');
% 
% phasepow=[phase_rounded lowpwr];
n_hist_bins = 30;

phase_edges=linspace(min(thetaphase),max(thetaphase),n_hist_bins+1);
amp_by_phases_low=zeros(1,n_hist_bins);
amp_by_phases_high=zeros(1,n_hist_bins);

for i=1:n_hist_bins-1
    amp_by_phases_low(i) = mean(lowpwr(thetaphase>phase_edges(i) & thetaphase<phase_edges(i+1)));
    amp_by_phases_high(i) = mean(highpwr(thetaphase>phase_edges(i) & thetaphase<phase_edges(i+1)));
end

figure
bar(phase_edges(1:end-1),amp_by_phases_low,'histc');
hold on
bar(phase_edges(1:end-1),amp_by_phases_high,'histc');

    
%from peters
% Find peak indexes higher than 2 std
% phs.lg_pk_ind = find((pk.lg(:,2) > 2*std(env.lg)));
% phs.hg_pk_ind = find((pk.hg(:,2) > 2*std(env.hg)));
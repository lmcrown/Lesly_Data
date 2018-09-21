GP = LD_Globals_LC();
load('Session_Info.mat')
% Theta_file=ses.other.ThetaEEG_FileName;


fs=1000;
%%
[LFP, sFreq, saturation_times, head] = ReadCR_cowen('CSC08.ncs',[],[],fs,1);
%the actual sample rate is 1890.4, but can input something different so
%I'll make it come out in 1000
%%
%plot with turned into hours
figure
plot(LFP(:,1)/3600e6,LFP(:,2))

%%

[pxx,f] = pwelch(LFP(:,2),500,300,500,fs);
plot(f,10*log10(pxx))

xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

% do I need to notch this? It looks like I am seeing 60Hz harmonics
%%
%slow gamma according to Zheng,Colgin: 25-50 Hz      high gamma 65-100 
LowGamma = designfilt('bandpassiir','FilterOrder',12, ...
         'HalfPowerFrequency1',25, 'HalfPowerFrequency2',50, ...
         'SampleRate',fs,'designmethod', 'butter');
                    as= filtfilt(LowGamma,LFP(:,2));
                    lowpwr=envelope_cowen(abs(as));

HighGamma = designfilt('bandpassiir','FilterOrder',12, ...
         'HalfPowerFrequency1',65, 'HalfPowerFrequency2',100, ...
         'SampleRate',fs,'designmethod', 'butter');
                    as= filtfilt(HighGamma,LFP(:,2));
                    highpwr=envelope_cowen(abs(as)); %this vs hilbert...
                    
Theta = designfilt('bandpassiir','FilterOrder',12, ...
         'HalfPowerFrequency1',6, 'HalfPowerFrequency2',10, ...
         'SampleRate',fs,'designmethod', 'butter');
                    as= filtfilt(Theta,LFP(:,2));
                    thetapwr=abs(hilbert(as));                    
                    
                    figure;
                    plot(LFP(:,1),highpwr)
                    hold on
                    plot(LFP(:,1),lowpwr)
                       hold on
                    plot(LFP(:,1),thetapwr)
                    
%%
%epoch times are in tenths of milliseconds. This means that to get them to microseconds you will need to multiply them by 100

                    
                    
%%
MazeoneIX= LFP(:,1)>ses.epochs.Maze1(1)*100 & LFP(:,1)< ses.epochs.Maze1(2)*100;
data=LFP(MazeoneIX,:);
frex=1:150;
 [pow]=real(SPEC_cwt_cowen(data(1:length(data)/200),fs,frex,32,0));
 figure;contourf(1:Cols(pow),frex,pow)
 
%%
GP.cm_per_pixel;
POS = load('vtm1.pvd');
% I think the last 2 cols are always crap.
POS = POS(:,1:3);
% This file will also be quite useful...
load PositionStruct.mat % loads into 'pos' lowercase.
POS = [ Range(pos.All.X) Data(pos.All.X) Data(pos.All.Y)];
SPEED=Speed_from_xy(POS,19,GP.cm_per_pixel);
figure;
plot(POS(:,1),SPEED)
ylabel('Speed (cm/s')
%thing is there seems to be a whole FindVelocity.m that peter did
[vs1, vs2, t1, t2] = FindVelocity(pos, GP.Maze_Diameter_pixels); %must be done with the lower case pos because its a structure
velocity=zeros(length(POS)); %set to zero and then fill in with velocity values
velocity

%WILL NEEED to to make logicals for sleep and wake based on epochs that work for velocity, because of the tenths of milliseconds things 
MazeoneIX= POS(:,1)>ses.epochs.Maze1(1) & POS(:,1)< ses.epochs.Maze1(2); %if epoch time and POS time are the same....
%test by adding vertical lines on the POS data
x1=epochs.Maze1(1);
line([x1 x1], get(gca, 'ylim'));
figure;
plot(vs1) %vs1 has positive and negative values to indicate direction, thats cool- you get diff ones for maze1 and 2
% OK SO POS DATA WORKS WITH EPOCHS ITS ONLY LFP THAT IS WEIRD
%this is the one that seems to give you the linearized position- from
%Stephen's demo
figure;
plot(POS(1:3:end,1)/1e6/60,POS(1:3:end,2))
% Load the epoch data:
% these data are very handy - all of the intervals are here for every rest
% and maze running epoch.
load('epochs.mat');
epochs


figure
subplot(2,1,1)
plot(LFP(:,1)/3600e6,LFP(:,2))
subplot(2,1,2)

%%
% Find Phase of theta
phs.theta = angle(hilbert(Data(lfp.theta)));
%%

% x = 0:.1:4*pi;
% % d = sin(x);
% x= new_pos_pha{iCell}(:,2)
% 
% % d_noise = d + (rand(size(x))-.5)/2;
% % d_noise2 = d + (rand(size(x))-.5);
% figure
% %x has to be phase bins,
% % y has to be counts in bins- so need to bin the data
% y=new_pos_pha{iCell}(:,3) 
%% %now you have the shape of the histogram with phase at the bottom
%  h=histogram([new_pos_pha{iCell}(:,3); new_pos_pha{iCell}(:,3) + 2*pi],40);

bins=-pi:0.25:pi
h=histogram([new_pos_pha{iCell}(:,3)],bins);

x=h.BinEdges(1:end-1)
y=h.Values

figure;
plot(x,y)
% perhaps I need to dramatically smooth this to get it into just 2 peaks
% num=9
% w = hamming(num)
% m=conv(y,w(1:(end-(num-1))))

% figure;
% plot(x,m)

%ehh nan, probably better off trying to specify that i want the model to
%have 1 cycle
yu = max(y);
yl = min(y);
yr = (yu-yl);                               % Range of 'y'
yz = y-yu+(yr/2);
zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
% per = 2*mean(diff(zx)); % Estimate period   
%maybe i can make it have a period of around 3.5
per3=pi*2

ym = mean(y);                               % Estimate offset
% fit = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);    % Function to fit
fit = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);    % Function to fit
 % Function to fit

fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
s = fminsearch(fcn, [yr;  per3;  -1;  ym]);


% s = fminsearch(fcn, [yr;  per;  -1;  ym]);

%%
xp = linspace(min(x),max(x));
figure(1)
clf
plot(x,y,'b',  xp,fit(s,xp), 'r')
grid
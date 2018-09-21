function [BIX,sted_usec] = LD_Clean_LFP(LFP,Stim_times_usec, thresh, smooth_win)
% get rid of stim artifact and other artifact. LFP - 2 cols (time in usec,
% lfp. stim times in usec, thresh for artifact in units of lfp.
% return index of bad times to blank out or avoid.
BIX = false(Rows(LFP),1);
for iS = 1:Rows(Stim_times_usec)
    IX = LFP(:,1) >= Stim_times_usec(iS,1) & LFP(:,1) <= Stim_times_usec(iS,2);
    BIX(IX) = true;
end
moreBIX = abs(LFP(:,2))>thresh | [0;abs(diff( LFP(:,2)))] > thresh - 200;
moreBIX = moreBIX(:);
moreBIX = convn(moreBIX,ones(smooth_win,1),'same');
moreBIX = moreBIX>0;

BIX = moreBIX | BIX;

if nargout > 1
    ix = find_intervals( BIX, .5);
    sted_usec = [LFP(ix(:,1)) LFP(ix(:,2))];
end
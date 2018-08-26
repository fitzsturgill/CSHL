% dev fixing missing pupil video files by alignment of time stamps

fcell = {wtf(:).files};
tcell = {wtf(:).trials};
nSessions = length(fcell);
maxLength = max(max(cellfun(@length, fcell)), max(cellfun(@length, tcell)));

padLength = 20; % how many NaNs to add on to the end of maxLength
files = NaN(maxLength + padLength, nSessions); % to hold differences between successive timestamps for pupil files
trials = NaN(maxLength + padLength, nSessions); % ditto for trials

for counter = 1:nSessions
    tsf = length(fcell{counter});
    tst = length(tcell{counter});
    files(1:tsf,counter) = fcell{counter}(:);
    trials(1:tst,counter) = tcell{counter}(:);
end


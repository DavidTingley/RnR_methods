function [data counts] = process_replayData(spkmat, ts, keep, binSize)
% function that 'pulls out' a start and stop time around each
% ripple/population event. This is an entirely arbitrary heuristic, I have
% attempted to match the behavior of this function to some of the most
% common heuristics used by the field. 
%
% INPUTS
%
% ts - timestamp of peak of event (in seconds)


start = round((round(ts(1) * 1000)-50) ./ (spkmat.dt*1000));
stop = round((round(ts(2) * 1000)+50) ./ (spkmat.dt*1000));


for spk = 1:size(spkmat.data,2)
data(:,spk) = (spkmat.data(start:stop,spk)')';  
counts(:,spk) = (spkmat.data(start:stop,spk)')';  
end

% cut 0 and 1 spk count bins from the start/end
while sum(counts(1,keep)) < 2 & size(counts,1) > 1
    data = data(2:end,:);
    counts = counts(2:end,:);
    start = start + 1;
end
while sum(counts(end,keep)) < 2 & size(counts,1) > 1
    data = data(1:end-1,:);
    counts = counts(1:end-1,:);
    stop = stop-1;
end
data = data ./ binSize;
% this script attempts to calculate false positive rates for the rank order 
% correlation method given a number typical numbers of neurons and events
% recorded
%
% Copyright (C) 2019 Adrien Peyrache & David Tingley.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


seqRange = 3:10; % the typical number of unique neurons that spike in a given ripple or candidate replay event
numEvents = 100; % typical number of candidate events to examine, within a single session
seqLeng = nan(length(seqRange),1);
seqLengMax = nan(length(seqRange),1);
numIterations = 100;

for seqLen = seqRange %% iterate through sequences of different length
    co = nan(numIterations,numEvents);
    pval = nan(numIterations,numEvents);
    
    parfor iter = 1:numIterations
        seq1=1:seqLen; % this is our 'target' sequence (i.e. behavioral template)
        for rip = 1:numEvents
            seq2 = seq1(randperm(seqLen));
            [co(iter,rip),pval(iter,rip)] = corr(seq1',seq2');
            
            %% simulate some 'within-event' shuffling
            for sh = 1:100
                [co_shuf(sh,iter,rip),pval_shuf(sh,iter,rip)] = corr(seq1',seq2(randperm(seqLen))');
            end
        end
    end
    % across event
    seqLeng(seqLen) = mean(sum(pval<.05,2)); % average FP rate
    seqLengMax(seqLen) = max(sum(pval<.05,2)); % the worst case FP rate
    % within event
    seqLeng_within(seqLen) = mean(mean(sum(pval_shuf<.05,3))); % average FP rate
    seqLengMax_within(seqLen) = max(mean(sum(pval_shuf<.05,3))); % the worst case FP rate
end

subplot(2,2,1)
histogram(sum(pval<.05,2))

subplot(2,2,2)
plot(seqLeng)
hold on
plot(seqLengMax)




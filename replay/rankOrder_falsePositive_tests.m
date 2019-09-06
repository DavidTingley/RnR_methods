% this script attempts to calculate false positive rates for the rank order 
% correlation method given a number typical numbers of neurons and events
% recorded
%
% david tingley 2019

seqRange = 3:10; % the typical number of unique neurons that spike in a given ripple or candidate replay event
numEvents = 100; % typical number of candidate events to examine, within a single session

for seqLen = seqRange
    for iter = 1:1000
        seq1=1:seqLen;
        for rip = 1:numEvents
            seq2 = seq1(randperm(seqLen));
            [co(iter,rip),pval(iter,rip)] = corr(seq1',seq2');
        end
    end
    seqLeng(seqLen) = mean(sum(pval<.05,2)); % average FP rate
    seqLengmax(seqLen) = max(sum(pval<.05,2)); % the worst case FP rate
end

subplot(2,2,1)
histogram(sum(pval<.05,2))

subplot(2,2,2)
plot(seqLeng)
hold on
plot(seqLengmax)
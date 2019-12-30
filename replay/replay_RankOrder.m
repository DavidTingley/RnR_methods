function [replayScores] = replay_RankOrder(spikes,ripples,template,include)
% Estimates the rank order correlation method for replay quantification
%
% USAGE
%
% [] = compareReplayMethods(spikes,ripples,template,include)
%
% INPUTS
%
% spikes  - buzcode cellinfo file. A structure that only requires three fields: 
%              -spikes.times, an Nx1 cell array of timestamps in seconds for each neuron
%              -spikes.UID, Nx1 vector of unique identifiers for each neuron in session
%              -spikes.spindices, a sorted vector of [spiketime UID], for all spikes. 
%                               useful input to some functions and plotting rasters
% ripples - buzcode ripples events file. A structure that only requires two fields:
%             -ripples.timestamps, an Nx2 matrix of start/stop times
%             -ripples.peaks, an Nx1 vector of peak timestamps for each event
% template -NxD matrix of N cells and D positions, average firing ratemaps
%           for each cell
% include - indices (1:N) of cells (i.e.place cells) to keep
%
% OUTPUTS 
% 
% replayScore.rankOrd - rank order correlations w/ template for each event
%
% HELP
%
% See bz_GetSpikes from the buzcode repo for help with the spikes/ripples
%   data structures
%
% Copyright (C) 2019 Adrien Peyrache & David Tingley.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

nCellsPerEvt = 4; % the minimum number of cells firing in an event to examine it
binSize = .01;
overlap = 1;
spkmat = bz_SpktToSpkmat(spikes.times,'overlap',overlap,'binSize',binSize * overlap); % converts spike timestamps into matrix format

if ~isempty(include)
    keep = intersect(include,find(sum(template')>0)); % uncomment to keep all HPC cells that fired
else
    keep = find(sum(template')>0); % just remove cells w/ zeros spikes in template
end

for event = 1:size(ripples.timestamps,1)
    
    ts = ripples.timestamps(event,:); % gets start/stop timestamp of ripple event
    [data counts] = process_replayData(spkmat, ts, keep, binSize); % processing func to get out the 'event' using common FR heuristics
    
    
    idx = intersect(find(sum(data)>0),keep); % find cells that fired and are in the 'include' array (i.e. active place cells)
    if length(idx) >= nCellsPerEvt
        [~,~,ord_template] = sort_cells(template(idx,:));
        [~,ord_firstSpk] = sortrows(data(:,idx)','descend');
        [rankOrd(event),pvals(event)] = corr(ord_template,ord_firstSpk,'rows','complete');
        [rankOrd_shuf(event),pvals_shuf(event)] = corr(ord_template(randperm(length(idx))),ord_firstSpk,'rows','complete');
    else
        rankOrd(event) = nan;
        pvals(event) = nan;
        rankOrd_shuf(event) = nan;
        pvals_shuf(event) = nan;
    end       
    
    % extra data to send up
    nCells(event) = sum(sum(counts(:,keep))>0);
    nSpks(event) = sum(sum(counts(:,keep)));
end

% create data struct to return
replayScores.rankOrd = rankOrd;
replayScores.rankOrd_shuf = rankOrd_shuf;
replayScores.nCells = nCells;
replayScores.nSpks = nSpks;


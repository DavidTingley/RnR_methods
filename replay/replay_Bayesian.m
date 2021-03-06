function [replayScores] = replay_Bayesian(spikes,ripples,template,include)
% Estimates the Bayesian method for replay quantification
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
% bayesRadon - integral under the line of best fit, using the Radon
%              transform (Davidson 2009)
% bayesLinearWeighted - linear weighted correlation of the posterior
%              probability matrix (Grosmark 2016)
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

nBinsThresh = 5; % the minimum number of bins, per event, to analyze event (binSize * nBinsThresh = total event duration)
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
    
    if size(data,1) >= nBinsThresh
        [Pr, prMax] = placeBayes(data(:,keep), template(keep,:), spkmat.dt); % generate posterior probability matrix using template and event FRs
        Pr(isnan(Pr)) = 0; % bad form... but some events have 0 spks in a particular time bin, doing this rather than filtering those events out
        [bayesLinearWeighted(event),outID] = makeBayesWeightedCorr1(Pr,ones(size(Pr,1),1)); % linear weight correlation method of quantification (Grosmark/Buzsaki 2016)
        [slope_hpc(event),bayesRadon(event)] = Pr2Radon(Pr'); % Radon transform method of quantification (Davidson/Frank 2009)   
        
        %% let's add some shuffling
        % cell ID
        [Pr, prMax] = placeBayes(data(:,keep), bz_shuffleCellID(template(keep,:)), spkmat.dt); % generate posterior probability matrix using template and event FRs
        Pr(isnan(Pr)) = 0; % bad form... but some events have 0 spks in a particular time bin, doing this rather than filtering those events out
        [bayesLinearWeighted_cellID_shuf(event),outID] = makeBayesWeightedCorr1(Pr,ones(size(Pr,1),1)); % linear weight correlation method of quantification (Grosmark/Buzsaki 2016)
        [slope_hpc_cellID_shuf(event),bayesRadon_cellID_shuf(event)] = Pr2Radon(Pr'); % Radon transform method of quantification (Davidson/Frank 2009) 
        % circular
        [Pr, prMax] = placeBayes(data(:,keep), bz_shuffleCircular(template(keep,:)), spkmat.dt); % generate posterior probability matrix using template and event FRs
        Pr(isnan(Pr)) = 0; % bad form... but some events have 0 spks in a particular time bin, doing this rather than filtering those events out
        [bayesLinearWeighted_circular_shuf(event),outID] = makeBayesWeightedCorr1(Pr,ones(size(Pr,1),1)); % linear weight correlation method of quantification (Grosmark/Buzsaki 2016)
        [slope_hpc_circular_shuf(event),bayesRadon_circular_shuf(event)] = Pr2Radon(Pr'); % Radon transform method of quantification (Davidson/Frank 2009)  
        
    else
        bayesLinearWeighted(event) = nan;
        slope_hpc(event) = nan;
        bayesRadon(event) = nan;
        
        bayesLinearWeighted_cellID_shuf(event) = nan;
        slope_hpc_cellID_shuf(event) = nan;
        bayesRadon_cellID_shuf(event) = nan;
        bayesLinearWeighted_circular_shuf(event) = nan;
        slope_hpc_circular_shuf(event) = nan;
        bayesRadon_circular_shuf(event) = nan;
    end
    % extra data to send up
    nCells(event) = sum(sum(counts(:,keep))>0);
    nSpks(event) = sum(sum(counts(:,keep)));
end

% create data struct to return
replayScores.bayesLinearWeighted = bayesLinearWeighted;
replayScores.bayesRadon = bayesRadon;

replayScores.bayesLinearWeighted_cellID_shuf = bayesLinearWeighted_cellID_shuf;
replayScores.bayesRadon_cellID_shuf = bayesRadon_cellID_shuf;
replayScores.bayesLinearWeighted_circular_shuf = bayesLinearWeighted_circular_shuf;
replayScores.bayesRadon_circular_shuf = bayesRadon_circular_shuf;


replayScores.nCells = nCells;
replayScores.nSpks = nSpks;


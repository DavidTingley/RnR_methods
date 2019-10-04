% this script runs the rank-order correlation and Radon integral replay methods
% on simulated data where replay events are systematically degraded, either
% by removed 'real' spikes, or adding spurious 'noise' spikes.
%
% Copyright (C) 2019 Adrien Peyrache & David Tingley.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


%% initial parameters to play with 
numCells = 100; % # of place fields for simulations
radonBinSize = 15; % 15 ms bins for the Bayesian method
numIterations = 2; % number of sims to run for each parameter set, 3-10 is usually a good range (smaller #s run faster)
inFieldMeanRate = 15; % average FR for distribution of fields, measured in Hz;

% pick one of the following 
fieldRateDistro = ones(numCells,1) .* inFieldMeanRate; % everyone is equal
% fieldRateDistro = ones(numCells,1) .* normrnd(inFieldMeanRate,inFieldMeanRate./3,numCells,1); % gaussian distributed in-field rates
% fieldRateDistro = ones(numCells,1) .* lognrnd(log(inFieldMeanRate),log(inFieldMeanRate)./3,numCells,1); % lognorm distributed in-field rates
fieldRateDistro(fieldRateDistro>120) = 120; % cap for lognorm distro


%% first, let's make some place fields.
offsets_rate = {round(sigmoid(1:100,50,.09) .* 100)... % reward/end over-representation
          1:100 ...                                    % linearly spaced PFs
          randi(100,1,100)};                           % randomly spaced PFs

for o = 1:length(offsets_rate)
    for cell = 1:numCells
        rateMaps{o}(cell,:) = minmax_norm(Smooth([zeros(1,offsets_rate{o}(cell)) 1 zeros(1,100-offsets_rate{o}(cell))]',5)) .* fieldRateDistro(cell); % this multiplier makes the in-field FR's about ~15Hz for each cell 
    end
end

%% ok, now let's make some ripple examples..
offsets_rip = {round([(numCells/2-logspace(2,0.8,50)./2)+3 (logspace(0.8,2,50)./2)+numCells/2-3])... % logarithmicaly spaced spikes to replicate population FR during ripples
               1:100}; % linearly spaced spikes within each ripple

%% we're going to iterate thru adding noise, and subtracting signal           
noise_rankOrd = nan(100,100,numIterations);
noise_integral = nan(100,100,numIterations);
noise_reactICA = nan(100,100,numIterations);
noise_reactPCA = nan(100,100,numIterations);

for nSub = 1:100 % subtract N 'real' spikes
    for nAdd = 1:100 % add N 'noise' spikes
        for o = 2        %% this can be varied [1-3] if you would like a different 'in-ripple' firing rate pattern (see lines 27-28)
            for oo = 2   %% this can be varied [1-3] if you would like a different PF tiling of space (see lines 12-21)
                for cell =1:numCells
                   rippleEvent{o}(cell,:) = ([zeros(1,offsets_rip{o}(cell)) 1 zeros(1,100-offsets_rip{o}(cell))]);
                end
                for iter = 1:numIterations 
                    spks = find(rippleEvent{o}==1);
                    rip = rippleEvent{o}; 
                    r = randperm(100);
                    rip(spks(r(1:nSub))) = 0;
                    r = randperm(length(rip(:)));
                    rip(r(1:nAdd)) = 1;
                    keep = find(sum(rip')>0); % used for rank order corr

                    %% discretize for radon integral here
                    for c = 1:size(rip,1)
                       rip_smooth(c,:) = rebin(rip(c,:),round(size(rip,2)/radonBinSize)); % 15 ms bins default (change on line 9)
                    end
                    % radon transform
                    [Pr prMax] = placeBayes(rip_smooth',rateMaps{oo},radonBinSize/1000);
                    [slope integral{o,oo}(iter)] = Pr2Radon(Pr');

%                     shuf = bz_shuffleCircular(rateMaps{oo});
%                     [Pr_shuf prMax] = placeBayes(rip_smooth',shuf,radonBinSize/1000);
%                     [slope_shuffle integral_shuffle{o,oo}(iter)] = Pr2Radon(Pr_shuf');

                    %% rank-order correlations
                    [a b ord] = sort_cells(rateMaps{oo}(keep,:));
%                     [a b ord_shuf] = sort_cells(shuf(keep,:));

                    [a b ord2] = sort_cells(rip(keep,:));

                    rankOrder{o,oo}(iter) = corr(ord,ord2);
%                     rankOrder_shuf{o,oo}(iter) = corr(ord_shuf,ord2);

                    %% reactivation analyses
                    [R,phi] = ReactStrength(rateMaps{oo}',[ rip_smooth]','method','pca');
                    reactPCA{o,oo}(iter) = mean(R(:,1));
                    [R,phi] = ReactStrength(rateMaps{oo}',[ rip_smooth]','method','ica');
                    reactICA{o,oo}(iter) = mean(R(:,1));
                end
                noise_rankOrd(nSub,nAdd,:) = (rankOrder{o,oo});
                noise_integral(nSub,nAdd,:) = (integral{o,oo});
                noise_reactICA(nSub,nAdd,:) = reactICA{o,oo};
                noise_reactPCA(nSub,nAdd,:) = reactPCA{o,oo};
            end
        end
        
        
        figure(1) % reactivation stuff
        subplot(3,2,1)
        imagesc(rip)
        title('ripple example')
        
        subplot(3,2,2)
        imagesc(rateMaps{oo})
        title('PF ratemaps')
        
        subplot(3,2,3)
        plot(ord,ord2,'.k') % visualize what 'example' events look like as the simulations run
        xlabel('PF order')
        ylabel('ripple order')
        title(['removed: ' num2str(nSub) ', added: ' num2str(nAdd) ', rank order: ' num2str(rankOrder{o,oo}(end))])
        
        subplot(3,2,4)
        Pr2Radon(Pr',1); % visualize what 'example' events look like as the simulations run
        title(['removed: ' num2str(nSub) ', added: ' num2str(nAdd) ', radon integral: ' num2str(integral{o,oo}(end))])
        ylabel('decoded position')
        xlabel(['timebins (' num2str(radonBinSize) ' ms)'])
        
        subplot(3,2,5)
        imagesc(squeeze(mean(noise_rankOrd,3)));
        title('rank order')
        xlabel('# added "noise" spks')
        ylabel('# removed "real" spks')
        
        subplot(3,2,6)
        imagesc(squeeze(mean(noise_integral,3)));
        title('radon integral')
        xlabel('# added "noise" spks')
        ylabel('# removed "real" spks')
        
        figure(2) % reactivation stuff
        
        subplot(2,2,1)
        imagesc(squeeze(mean(noise_reactPCA,3)))
        title('react strength PCA')
        
        subplot(2,2,2)
        imagesc(squeeze(mean(noise_reactICA,3)))
        title('react strength ICA')
        
        subplot(2,2,3)
        imagesc(squeeze(mean(noise_reactPCA,3))-squeeze(mean(noise_reactICA,3)))
        title('PCA minus ICA')
        
        pause(.001)
    end
end





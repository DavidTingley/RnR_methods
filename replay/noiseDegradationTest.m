% this script runs the rank-order correlation and Radon integral replay methods
% on simulated data where replay events are systematically degraded, either
% by removed 'real' spikes, or adding spurious 'noise' spikes.
%
% david tingley 2019

%% initial parameters to play with 
numCells = 100; % # of place fields for simulations
radonBinSize = 15; % 15 ms bins for the Bayesian method
numIterations = 2; % number of sims to run for each parameter set, 3-10 is usually a good range (smaller #s run faster)
inFieldMeanRate = 15; % average FR for distribution of fields, measured in Hz;

% pick one of the following 
% fieldRateDistro = ones(numCells,1) .* inFieldMeanRate; % everyone is equal
fieldRateDistro = ones(numCells,1) .* normrnd(inFieldMeanRate,inFieldMeanRate./3,numCells,1); % gaussian distributed in-field rates
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

for nSub = 1:100 % subtract one 'real' spike
    for nAdd = 1:100 % add one 'noise' spike
        for o = 2 %1:length(offsets_rip)            %% this can be varied if you would like a different 'in-ripple' firing rate pattern (see lines 27-28)
            for oo = 2 %1:length(offsets_rate)      %% this can be varied if you would like a different PF tiling of space (see lines 12-21)
                for cell =1:numCells
                   rippleEvent{o}(cell,:) = ([zeros(1,offsets_rip{o}(cell)) 1 zeros(1,100-offsets_rip{o}(cell))]);
                end
                for iter = 1:numIterations % seems to be enough for averaging, without being too slow
                    spks = find(rippleEvent{o}==1);
                    rip = rippleEvent{o}; 
                    r = randperm(100);
                    rip(spks(r(1:nSub))) = 0;
                    r = randperm(length(rip(:)));
                    rip(r(1:nAdd)) = 1;
                    keep = find(sum(rip')>0);

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
                end
                noise_rankOrd(nSub,nAdd,:) = (rankOrder{o,oo});
                noise_integral(nSub,nAdd,:) = (integral{o,oo});
            end
        end
        
        subplot(2,2,1)
        plot(ord,ord2,'.k') % visualize what 'example' events look like as the simulations run
        xlabel('PF order')
        ylabel('ripple order')
        title(['removed: ' num2str(nSub) ', added: ' num2str(nAdd) ', rank order: ' num2str(rankOrder{o,oo}(end))])
        
        subplot(2,2,2)
        Pr2Radon(Pr',1); % visualize what 'example' events look like as the simulations run
        title(['removed: ' num2str(nSub) ', added: ' num2str(nAdd) ', radon integral: ' num2str(integral{o,oo}(end))])
        ylabel('decoded position')
        xlabel(['timebins (' num2str(radonBinSize) ' ms)'])
        
        subplot(2,2,3)
        imagesc(squeeze(mean(noise_rankOrd,3)));
        title('rank order')
        xlabel('# added "noise" spks')
        ylabel('# removed "real" spks')
        
        subplot(2,2,4)
        imagesc(squeeze(mean(noise_integral,3)));
        title('radon integral')
        xlabel('# added "noise" spks')
        ylabel('# removed "real" spks')
        pause(.001)
    end
end





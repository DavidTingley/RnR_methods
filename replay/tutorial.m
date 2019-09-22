% a general script that creates place fields and 'replay' events, then shows 
% the rank order correlation and radon integral for these events, compared to shuffled data.
%
% david tingley 2019


%% initial parameters
numCells = 100; % # of place fields for simulations



%% first, let's make some place fields.
offsets_rate = {round(sigmoid(1:100,50,.09) .* 100)... % Dupret style, reward/end over-representation
          1:100 ...                                    % linearly spaced PFs
          randi(100,1,100)};                           % randomly spaced PFs

%A/ I've pre-allocated rateMaps. Sorry, can't let it happen :D
rateMaps = cell(length(offsets_rate),1);

for o = 1:length(offsets_rate)
    for neuron = 1:numCells
        rateMaps{o}(neuron,:) = Smooth([zeros(1,offsets_rate{o}(neuron)) 1 zeros(1,100-offsets_rate{o}(neuron))]',5); 
    end
end


%% ok, now let's make some ripple examples..
offsets_rip = {round([(numCells/2-logspace(2,0.8,50)./2)+3 (logspace(0.8,2,50)./2)+numCells/2-3])... 
               1:100};
       
for o =1:length(offsets_rip)
    for oo = 1:length(offsets_rate)
        
    for neuron =1:numCells
       rippleEvent{o}(neuron,:) = ([zeros(1,offsets_rip{o}(neuron)) 1 zeros(1,100-offsets_rip{o}(neuron))]);
    end
    for iter = 1:100
        spks = find(rippleEvent{o}==1);
        rip = rippleEvent{o}; 
        r = randperm(100);
        rip(spks(r(1:80))) = 0;
%         imagesc(rip)
        
        % radon transform
        [Pr, prMax] = placeBayes(rip',rateMaps{oo},1); %A/ prMax is not used
        [slope, integral{o,oo}(iter)] = Pr2Radon(Pr);
        
        shuf        = bz_shuffleCircular(rateMaps{oo});
        [Pr, prMax] = placeBayes(rip',shuf,1); %A/ prMax is not used
        [slope_shuffle, integral_shuffle{o,oo}(iter)] = Pr2Radon(Pr); %A/ slope_shuffle is not used
        
        % rank-order correlations
        [~, ~, ord]      = sort_cells(rateMaps{oo});
        [~, ~, ord_shuf] = sort_cells(shuf);
        [~, ~, ord2]     = sort_cells(rippleEvent{o});
        
        rankOrder{o,oo}(iter) = corr(ord,ord2);
        rankOrder_shuf{o,oo}(iter) = corr(ord_shuf,ord2);
    end
    end
end

%% Plotting results

conditions = length(offsets_rate)*length(offsets_rip);
cond = 1;
for o = 1:length(offsets_rip)
    for oo = 1:length(offsets_rate)
        subplot(conditions,4,cond*4-3)
        imagesc(rateMaps{oo})
        xlabel('position')
        ylabel('neuron #')
        
        subplot(conditions,4,cond*4-2)
        imagesc(rippleEvent{o})
        title('ripple template')
        xlabel('time bin')
        
        subplot(conditions,4,cond*4-1)
        histogram(integral{o,oo},[0:.001:.05])
        hold on
        histogram(integral_shuffle{o,oo},[0:.001:.05])    
        title('radon integral')
        
        subplot(conditions,4,cond*4)
        histogram(rankOrder{o,oo}) %A/ don't understand this part. it's all the same values in my hands.
        hold on
        histogram(rankOrder_shuf{o,oo})
        title('rank order correlation')
        
        cond = 1+cond;
    end
end


% subplot(4,2,1)
% imagesc(rateMaps)
% subplot(4,2,2)
% plot(mean(rateMaps))





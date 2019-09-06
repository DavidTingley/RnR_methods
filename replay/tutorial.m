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

for o = 1:length(offsets_rate)
    for cell = 1:numCells
        rateMaps{o}(cell,:) = Smooth([zeros(1,offsets_rate{o}(cell)) 1 zeros(1,100-offsets_rate{o}(cell))]',5); 
    end
end




%% ok, now let's make some ripple examples..
offsets_rip = {round([(numCells/2-logspace(2,0.8,50)./2)+3 (logspace(0.8,2,50)./2)+numCells/2-3])... 
               1:100};
       
for o =1:length(offsets_rip)
    for oo = 1:length(offsets_rate)
        
    for cell =1:numCells
       rippleEvent{o}(cell,:) = ([zeros(1,offsets_rip{o}(cell)) 1 zeros(1,100-offsets_rip{o}(cell))]);
    end
    for iter = 1:100
        spks = find(rippleEvent{o}==1);
        rip = rippleEvent{o}; 
        r = randperm(100);
        rip(spks(r(1:80))) = 0;
        imagesc(rip)
        
        % radon transform
        [Pr prMax] = placeBayes(rip',rateMaps{oo},1);
        [slope integral{o,oo}(iter)] = Pr2Radon(Pr);
        
        shuf = bz_shuffleCircular(rateMaps{oo});
        [Pr prMax] = placeBayes(rip',shuf,1);
        [slope_shuffle integral_shuffle{o,oo}(iter)] = Pr2Radon(Pr);
        
        % rank-order correlations
        [a b ord] = sort_cells(rateMaps{oo});
        [a b ord_shuf] = sort_cells(shuf);
        
        [a b ord2] = sort_cells(rippleEvent{o});
        
        rankOrder{o,oo}(iter) = corr(ord,ord2);
        rankOrder_shuf{o,oo}(iter) = corr(ord_shuf,ord2);
    end
    end
end

% plotting
conditions = length(offsets_rate)*length(offsets_rip);
cond = 1;
for o = 1:length(offsets_rip)
    for oo = 1:length(offsets_rate)
        subplot(conditions,4,cond*4-3)
        imagesc(rateMaps{oo})
        xlabel('position')
        ylabel('cell #')
        
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
        histogram(rankOrder{o,oo})
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





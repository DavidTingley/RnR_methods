% This script uses a decoder to estimate the 'position' when two place
% cells fire a single spike each at the same moment in time. It assumes
% these two neurons have different in-field firing rates. It
% demonstrates that the choice of bin size can determine the decoded
% position.

% Copyright (C) 2019 Adrien Peyrache & David Tingley.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


%% first, let's make some place fields.
offsets_rate = {round(sigmoid(1:100,50,.09) .* 100)... % reward/end over-representation
          1:100 ...                                    % linearly spaced PFs
          randi(100,1,100)}; 
      
for o = 1:length(offsets_rate)
    for cell = 1:100 
        rateMaps{o}(cell,:) = minmax_norm(Smooth([zeros(1,offsets_rate{o}(cell)) 1 zeros(1,100-offsets_rate{o}(cell))]',5)); % this multiplier makes the in-field FR's about ~15Hz for each cell 
    end
end


for rate =1:50 % 1 to 50 Hz
    r = rateMaps{2}([20 80],:);
    r(1,:) = rate*minmax_norm(r(1,:));
    r(2,:) = 50*minmax_norm(r(2,:));
    r(1,1:10)=0;
    r(1,31:end)=0;
    r(2,1:70)=0;
    r(2,91:end)=0;
    cl = poisson_naive_bayes_CL;
    cl = train(cl, r, 1:101);
    c=1; 
    for i=.001:.001:1
        [predicted_labels decision_values] = test(cl, [1 1]'/i);
        p(c,:)= decision_values;
        c=1+c;
    end
    for i=1:1000
        pz(i,:)=zscore(p(i,:));
    end
    for i=1:1000
        [a(i) b(i)] = max(pz(i,:));
    end
    bb(rate,:)= b;
    
    if rate == 25 % from paper, see what happens w/ 25Hz cell
        subplot(2,2,1)
        imagesc(r)
        title('place fields')
        xlabel('position')
        ylabel('cell #')

        subplot(2,2,2)
        imagesc([zeros(50,2); [1 1]; zeros(50,2)]')
        title('ripple event')
        xlabel('time (ms')
        ylabel('cell #')
        
        subplot(2,1,2)
        plot(bb(25,:))
        ylabel('decoded position')
        xlabel('bin size (ms)')
    end
end






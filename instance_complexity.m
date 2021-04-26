function instance_complexity(data)
    
    if nargin < 1; data = load_data; end
    
    cond = [data.cond];
    for c = 1:2
        subs = find(cond==c-1);
        i = 0;
        for s = subs
            for k = 2:6
                blocks = unique(data(s).learningblock(data(s).ns==k));
                for b = 1:length(blocks)
                    i = i + 1;
                    ix = data(s).learningblock == blocks(b);
                    p{k-1}(i,1) = length(unique(data(s).corchoice(ix)))/k;
                    r{k-1}(i,1) = nanmean(data(s).reward(ix));
                end
            end
        end
        
        for j = 2:4
            [R,P] = corr(p{j},r{j},'type','spearman')
        end
        clear p r
    end
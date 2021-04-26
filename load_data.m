function data = load_data
    
    load Collins_JN2014.mat
    
    T = {'ID' 'learningblock' 'ns' 'trial' 'state' 'image' 'folder' 'iter' 'corchoice' 'action' 'key' 'cor' 'reward' 'rt' 'cond' 'pcor' 'delay'};
    S = unique(expe_data(:,1));
    for s = 1:length(S)
        ix = expe_data(:,1)==S(s);
        for j = 1:length(T)
            data(s).(T{j}) = expe_data(ix,j);
        end
        data(s).ID = data(s).ID(1);
        data(s).cond = data(s).cond(1);
        data(s).C = 3;
        data(s).N = length(data(s).learningblock);
    end
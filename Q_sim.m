function simdata = Q_sim(x,data)
    
    % Simulation function for Q-learning agent.
    %
    % USAGE: simdata = Q_sim(x,data)
    
    if length(x)==2
        agent.C = [];
        agent.lrate_beta = 0;
        agent.lrate_p = 0;
        agent.beta = x(1);
        agent.lrate_Q = x(2);
    elseif length(x)==3
        agent.C = [];
        agent.lrate_beta = 0;
        agent.beta = x(1);
        agent.lrate_Q = x(2);
        agent.lrate_p = x(3);
    elseif length(x)==4
        agent.lrate_p = 0;
        agent.beta = x(1);
        agent.lrate_Q = x(2);
        agent.C = x(3);
        agent.lrate_beta = x(4);
    elseif length(x)==5
        agent.beta = x(1);
        agent.lrate_Q = x(2);
        agent.C = x(3);
        agent.lrate_beta = x(4);
        agent.lrate_p = x(5);
    end
    
    simdata = data;
    B = unique(data.learningblock);
    for b = 1:length(B)
        ix = find(data.learningblock==B(b));
        state = data.state(ix);
        action = zeros(size(state));
        corchoice = data.corchoice(ix);     % correct choice on each trial
        setsize = length(unique(state));    % number of distinct stimuli
        Q = zeros(setsize,3);           % policy xeters
        beta = agent.beta;
        p = ones(1,3)/3;                    % marginal action probabilities
        for t = 1:length(state)
            s = state(t);
            d = beta*Q(s,:) + log(p);
            logpolicy = d - logsumexp(d);
            policy = exp(logpolicy);    % softmax policy
            a = fastrandsample(policy); % action
            if a == corchoice(t)
                r = 1;                  % reward
            else
                r = 0;
            end
            cost = logpolicy(a) - log(p(a));    % policy complexity cost
            rpe = r - Q(s,a);   % reward prediction error
            action(t) = a;
            %R = mutual_information(state(1:t),action(1:t),0.1);
            if agent.lrate_beta > 0
                beta = beta + agent.lrate_beta*(agent.C-cost)*(Q(s,a)-(Q(s,:)*policy'));
                %beta = beta + agent.lrate_beta*(agent.C-R)*(Q(s,a)-(Q(s,:)*policy'));
                beta = max(min(beta,50),0);
            end
            if agent.lrate_p > 0
                p = p + agent.lrate_p*(policy - p); p = p./sum(p);                                  % marginal update
            end
            Q(s,a) = Q(s,a) + agent.lrate_Q*rpe;    % policy parameter update
            
            simdata.action(ix(t)) = a;
            simdata.reward(ix(t)) = r;
            simdata.beta(ix(t)) = beta;
            simdata.cost(ix(t)) = cost;
        end
    end
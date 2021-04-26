function simdata = actor_critic(agent,data)
    
    % Simulate actor-critic agent.
    %
    % USAGE: simdata = actor_critic(agent,data)
    
    if ~isfield(agent,'beta')
        agent.beta = agent.beta0;
    end
    
    simdata = data;
    B = unique(data.learningblock);
    for b = 1:length(B)
        ix = find(data.learningblock==B(b));
        state = data.state(ix);
        corchoice = data.corchoice(ix);     % correct choice on each trial
        setsize = length(unique(state));    % number of distinct stimuli
        theta = zeros(setsize,3);           % policy parameters
        V = zeros(setsize,1);               % state values
        beta = agent.beta;
        p = ones(1,3)/3;                    % marginal action probabilities
        for t = 1:length(state)
            s = state(t);
            d = beta*theta(s,:) + log(p);
            logpolicy = d - logsumexp(d);
            policy = exp(logpolicy);    % softmax policy
            a = fastrandsample(policy); % action
            if a == corchoice(t)
                r = 1;                  % reward
            else
                r = 0;
            end
            cost = logpolicy(a) - log(p(a));    % policy complexity cost
            rpe = beta*r - cost - V(s);   % reward prediction error
            g = rpe*beta*(1 - policy(a)); % policy gradient
            V(s) = V(s) + agent.lrate_V*rpe;                                                    % state value update
            if agent.lrate_beta > 0
                beta = beta + agent.lrate_beta*(agent.C-cost)*(theta(s,a)-(theta(s,:)*policy'));
                beta = max(min(beta,50),0);
            end
            if agent.lrate_p > 0
                p = p + agent.lrate_p*(policy - p); p = p./sum(p);                                  % marginal update
            end
            theta(s,a) = theta(s,a) + agent.lrate_theta*g/t;    % policy parameter update
            simdata.action(ix(t)) = a;
            simdata.reward(ix(t)) = r;
            simdata.expreward(ix(t)) = policy(corchoice(t));
            simdata.beta(ix(t)) = beta;
        end
    end
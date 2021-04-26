function [lik,latents] = actor_critic0_lik(x,data)
    
    % Likelihood function for actor-critic agent with no cost.
    %
    % USAGE: [lik,latents] = actor_critic0_lik(x,data)
    
    agent.beta = x(1);
    agent.lrate_theta = x(2);
    agent.lrate_p = x(3);
    agent.lrate_V = x(4);
    
    B = unique(data.learningblock);
    lik = 0;
    for b = 1:length(B)
        ix = data.learningblock==B(b);
        reward = data.reward(ix);
        action = data.action(ix);
        state = data.state(ix);
        setsize = length(unique(state));    % number of distinct stimuli
        theta = zeros(setsize,3);           % policy xeters
        V = zeros(setsize,1);               % state values
        beta = agent.beta;
        p = ones(1,3)/3;                    % marginal action probabilities
        if nargout > 1
            ii = find(ix);
        end
        for t = 1:length(state)
            s = state(t); a = action(t); r = reward(t);
            if a > 0
                d = beta*theta(s,:) + log(p);
                logpolicy = d - logsumexp(d);
                policy = exp(logpolicy);    % softmax policy
                lik = lik + logpolicy(a);
                rpe = beta*r - V(s);   % reward prediction error
                g = rpe*beta*(1 - policy(a)); % policy gradient
                V(s) = V(s) + agent.lrate_V*rpe;                                                    % state value update
                if agent.lrate_p > 0
                    p = p + agent.lrate_p*(policy - p); p = p./sum(p);                                  % marginal update
                end
                theta(s,a) = theta(s,a) + agent.lrate_theta*g/t;    % policy parameter update
                if nargout > 1
                    latents.rpe(ii(t)) = rpe;
                end
            end
        end
    end
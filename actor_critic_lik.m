function [lik,latents] = actor_critic_lik(x,data)
    
    % Likelihood function for actor-critic agent.
    %
    % USAGE: [lik,latents] = actor_critic_lik(x,data)
    
    if length(x)==3
        agent.C = [];
        agent.lrate_beta = 0;
        agent.lrate_p = 0;
        agent.beta0 = x(1);
        agent.lrate_theta = x(2);
        agent.lrate_V = x(3);
    elseif length(x)==4
        agent.C = [];
        agent.lrate_beta = 0;
        agent.beta0 = x(1);
        agent.lrate_theta = x(2);
        agent.lrate_V = x(3);
        agent.lrate_p = x(4);
    elseif length(x)==5
        agent.lrate_p = 0;
        agent.C = x(1);
        agent.lrate_theta = x(2);
        agent.lrate_V = x(3);
        agent.beta0 = x(4);
        agent.lrate_beta = x(5);
    elseif length(x)==6
        agent.C = x(1);
        agent.lrate_theta = x(2);
        agent.lrate_V = x(3);
        agent.beta0 = x(4);
        agent.lrate_beta = x(5);
        agent.lrate_p = x(6);
    end
    
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
        beta = agent.beta0;
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
                if nargout > 1
                    latents.rpe(ii(t)) = rpe;
                end
            end
        end
    end
function [lik,latents] = Q_lik(x,data)
    
    % Likelihood function for Q-learning agent.
    %
    % USAGE: [lik,latents] = Q_lik(x,data)
    
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
    
    B = unique(data.learningblock);
    lik = 0;
    for b = 1:length(B)
        ix = data.learningblock==B(b);
        reward = data.reward(ix);
        action = data.action(ix);
        state = data.state(ix);
        setsize = length(unique(state));    % number of distinct stimuli
        Q = zeros(setsize,3);           % policy xeters
        beta = agent.beta;
        p = ones(1,3)/3;                    % marginal action probabilities
        if nargout > 1
            ii = find(ix);
        end
        for t = 1:length(state)
            s = state(t); a = action(t); r = reward(t);
            if a > 0
                d = beta*Q(s,:) + log(p);
                logpolicy = d - logsumexp(d);
                policy = exp(logpolicy);    % softmax policy
                lik = lik + logpolicy(a);
                cost = logpolicy(a) - log(p(a));    % policy complexity cost
                %R = mutual_information(state(1:t),action(1:t),0.1);
                rpe = r - Q(s,a);   % reward prediction error
                if agent.lrate_beta > 0
                    beta = beta + agent.lrate_beta*(agent.C-cost)*(Q(s,a)-(Q(s,:)*policy'));
                    %beta = beta + agent.lrate_beta*(agent.C-R)*(Q(s,a)-(Q(s,:)*policy'));
                    beta = max(min(beta,50),0);
                end
                if agent.lrate_p > 0
                    p = p + agent.lrate_p*(policy - p); p = p./sum(p);                                  % marginal update
                end
                Q(s,a) = Q(s,a) + agent.lrate_Q*rpe;    % policy parameter update
                if nargout > 1
                    latents.rpe(ii(t)) = rpe;
                end
            end
        end
    end
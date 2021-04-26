function [simdata, simresults] = sim_collins14(data,results)
    
    % Simulate Collins et al. (2014) experiment using fitted actor-critic agent.
    
    rng(1); % set random seed for reproducibility
    
    if nargin < 1
        data = load_data;
    end
    
    if nargin < 2
        load model_fits;
        results = results(1);
    end
    
    for s = 1:length(data)
        
        agent.lrate_beta = 0;
        agent.lrate_p = 0;
        agent.C = [];
        for k = 1:length(results.param)
            agent.(results.param(k).name) = results.x(s,k);
        end
        
        simdata(s) = actor_critic(agent,data(s));
    end
    
    simresults = analyze_collins14(simdata);
    
    plot_figures('fig2',simresults,simdata);
    plot_figures('fig3',simresults,simdata);
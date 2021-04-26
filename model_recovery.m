function model_recovery
    
    load simresults_m3
    [results,bms_results] = fit_models(simdata);
    d = load('model_fits');
    [r,p] = corr(d.results(3).x,results(3).x)
    save model_fits_sim_m3 results bms_results
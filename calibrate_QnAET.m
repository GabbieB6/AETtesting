function [par_opt, of_cal, stopflag, caloutput]...
               = calibrate_QnAET(mmodel, obs, cal_idx, par_ini, optim_opts)
     
    % mmodel is a populated MARRMoT_model (i.e. it should contain
    %   mmodel.input_climate, .delta_t, .S0 and .solver_opts)
    % obs is a struct with obs.Q and obs.AET being arrays of the same
    %   length as mmodel.input_climte
    % cal_idx is a struct with obs.Q and obs.AET being either logical
    %   arrays of the same length of obs.Q and obs.AET or a numerical array
    %   of indices
    % par_ini is an array of initial parameters, same as the
    %   mmodel.calibrate method
    % optim_opts are the options to be passed to my_cmaes, same as the
    %   mmodel.calibrate method
    
    % this is a helper function that needs to be optimised, notice how it 
    % is only a function of the parameters, because that is the variable 
    % that needs to be varied every time the optimiser steps
    function fitness = fitness_fun(par)
        % run the model with the given parameter set
        output = mmodel.get_output([],[],par);
        
        % get the simulated Q and AET
        Qsim = output.Q;
        AETsim = output.Ea;
        
        % calculate the fitness, this will be returned by the helper funcn,
        % the fitness is inversed because my_cmaes is a minimiser, and we
        % need to maximise the fitness
        fitness = -1* OF_KGE_QnAET(obs.Q,  Qsim,  cal_idx.Q, ...
                                   obs.AET,AETsim,obs.IntYM,cal_idx.AET);
    end
    
    % here is where we perform the actual optimisation, using my_cmaes
     [par_opt,...                                              % optimal parameter set at the end of the optimisation
     of_cal,...                                                % value of the objective function at par_opt
     stopflag,...                                              % flag indicating reason the algorithm stopped
     caloutput] = ...                                          % output, see fminsearch for detail
               my_cmaes(...                                    % run my_cmaes
                     @fitness_fun,...                          % function to optimise is the fitness function
                     par_ini,...                               % initial parameter set
                     optim_opts);                              % optimiser options
 
 % re-inverse the final fitness (i.e. the calibration OF)
 of_cal = -1 * of_cal;

end
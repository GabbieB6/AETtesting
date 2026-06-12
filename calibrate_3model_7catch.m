%==========================================================================
% Main calibration workflow for:
%
% Burns et al. (2025)
% "A systematic evaluation of 15 actual evapotranspiration formulations
% within conceptual hydrological models"
%
% Author: Gabrielle Burns
% Institution: University of Melbourne
% Date: February 2024
%
% Description:
% This script implements the calibration and evaluation workflow used in
% Burns et al. (2026). The MARRMoT framework is used to evaluate 15
% evapotranspiration (ET) formulations within three conceptual hydrological
% models:
%
%   - GR4J
%   - SIMHYD
%   - VIC
%
% Models are calibrated against both observed streamflow and flux tower
% actual evapotranspiration (AET) using a combined multi-objective
% calibration framework.
%
% Workflow:
%   1. Data preparation
%   2. Model setup
%   3. Solver configuration
%   4. Calibration setup
%   5. Multi-objective calibration
%   6. Model evaluation
%   7. Output generation and visualisation 
%
% Requirements:
%   - MARRMoT toolbox
%   - Catchment forcing and evaluation datasets
%   - Supporting calibration functions included in this repository
%
%==========================================================================


% Study catchments:
% Wombat Forest (407221)
% Whroo (405229)
% Litchfield (G8150180)
% Dry River (G8140011)
% Robson Creek (111007A)
% Gingin (617003)
% Tumbarumba (401009)

clear all; close all;
%% add paths 
% addpath(genpath('../2_MARRMoT'))
addpath(genpath('../2_MARRMoT'))
addpath(genpath('data'))
addpath(genpath('calibration_codes'))

catchCodes = ["407221" "405229" "G8150180" "G8140011" "111007A" "617003" "401009"];
% [Wombat Forest, Whroo,  Litchfield, Dry River, Robson Creek, Gingin, Tumbarumba]

for iCatch = 1:size(catchCodes,2)
    
    %% 1. Prepare data
    % get this catchment 
    
    thisCatch = catchCodes(iCatch);

    % Load the data
    raw_data_table = readtable(append(thisCatch,"_ALLdaily.csv"));
    
    % chop dataset to exclude first years without streamflow
    data_table = raw_data_table(find(~isnan(raw_data_table.obsQ),1,'first'):end,:);

    % find the indexs for obs AET data
    ETStart = find(~isnan(data_table.FluxAET_adjust),1,"first");
    ETEnd   = find(~isnan(data_table.FluxAET_adjust),1,"last");

    % Create a climatology data input structure. 
    % NOTE: the names of all structure fields are hard-coded in each model
    % file. These should not be changed.
    input_climatology.precip   = data_table.precip;                             % Daily data: P rate  [mm/d]
    input_climatology.temp     = mean([data_table.tmax, data_table.tmin],2);    % Daily data: mean T  [degree C]
    input_climatology.pet      = data_table.pet;                                % Daily data: Ep rate [mm/d]
    input_climatology.delta_t  = 1;                                             % time step size of the inputs: 1 [d]
    
    % Create an integer vector to identify days in a given month year
    % (This is for the ET calibration teqnique, which utilises both a monthly
    % and annual time step)
    Years = data_table.date.Year; UniqueYears = unique(Years); NumYears = size(UniqueYears,1);
    Months = data_table.date.Month; NumMonths = size(Months,1); idx_yearmonth = nan(NumMonths,1);
    k = 0;
    
    for iYear = 1:NumYears
        ThisYear = UniqueYears(iYear);
        for iMonth = 1:12
            Ind = find(Years==ThisYear & Months ==iMonth);
            if isempty(Ind)
                continue;
            end
            k = k+1;
            idx_yearmonth(Ind) = k;               
       end  
    end
    
    % add the integer vector into the data table
    data_table = addvars(data_table,idx_yearmonth);
    
    % Extract observed streamflow & ET
    Q_obs = data_table.obsQ;
    aet_obs = data_table.FluxAET_adjust;
    
    %% 3. Define the solver settings  
    % Create a solver settings data input structure. 
    % NOTE: the names of all structure fields are hard-coded in each model
    % file. These should not be changed.
    input_solver_opts.resnorm_tolerance = 0.1;                                 % Root-finding convergence tolerance;
    % users have reported differences in simulation accuracy (KGE scores) during calibration between Matlab and Octave for a given tolerance.
    % In certain cases, Octave seems to require tigther tolerances to obtain the same KGE scores as Matlab does.
    input_solver_opts.resnorm_maxiter   = 6;                                   % Maximum number of re-runs
    
    %% 2. Define the model settings and create the model object
    % Model is GR4J_ET - which includes a flag to swap out the ET equation (GB)
    models = ["m_07_gr4j_4p_2s_ET" "m_18_simhyd_7p_3s_ET" "m_22_vic_10p_3s_ET"];
    model_names = ["GR4J" "Simhyd" "Vic"];
    
    for iModel = 1:3
        model     = models(iModel);          % Name of the model function (these can be found in Supporting Material 2)
        model_name = model_names(iModel);
        m = feval(model);
        
        % Time periods for calibration.
        % Five-year warm-up period used to minimise sensitivity to initial
        % storage conditions.
        n = length(Q_obs);
        warmup = 365*5;
        cal_idx.Q = warmup:n;
        cal_idx.AET = ETStart:ETEnd; % Restrict AET calibration to period with available flux tower data.
           
        obs.Q = Q_obs; obs.AET = aet_obs;
        obs.IntYM = data_table.idx_yearmonth;  % integer values assign for each day for a given month and year
        
        output_table = [];
        
        for iEqn = [1 2 3 4 6 7 8 11 13 16 19 20 21 22 23]
            % ET formulations selected for evaluation in Burns et al. (2025).
            % Formulations requiring additional state variables not available within
            % the tested model structures were excluded.

            model     = models(iModel);                                          
            m = feval(model);
            ThisETeqn = append('evap_', string(iEqn));
        
            [RunReq, OutFilePath] = IsRunReq(ThisETeqn, model_name, thisCatch);
	        
	        % if yes:
	        if RunReq
                disp(append('Running ', thisCatch,' - ',model_name, ' model, Cal QnAET, ' , ThisETeqn))
                m.evap_eqn = ThisETeqn;
                

                % Add ET-specific calibration parameters.
                % Additional parameters are appended to the base model parameter set and
                % subsequently transformed within get_evap_eqn.m to the formulation-
                % specific parameter ranges required by each ET equation.
                if (iEqn==2)||(iEqn==3)||(iEqn==6)||(iEqn==8)||(iEqn==16)||(iEqn==17)
                        m.numParams = m.numParams+1;
                        m.parRanges = [m.parRanges;      % model base params
                                            0.05   , 1];    % p1 [-]
                elseif iEqn==23
                        m.numParams = m.numParams+1;
                        m.parRanges = [m.parRanges;      % model base params
                                            0.05   , 1];    % p1 [-]
                elseif (iEqn==4)||(iEqn==13)||(iEqn==14)||(iEqn==19)||(iEqn==20)||(iEqn==21)||(iEqn==22)
                        m.numParams = m.numParams+2;
                        m.parRanges = [m.parRanges;      % model base params
                                           0.05   , 1;      % p1 [-]
                                          0.05   , 1];      % p2 [-]
                end

        
                parRanges = m.parRanges;                                                   % Parameter ranges
                numParams = m.numParams;                                                   % Number of parameters
                numStores = m.numStores;                                                   % Number of stores
                input_s0  = zeros(numStores,1);                                            % Initial storages (see note in paragraph 5 on model warm-up)
                        
                %% 4. Define calibration settings
                % Settings for 'my_cmaes'
                % the opts struct is made up of two fields, the names are hardcoded, so
                % they cannot be changed:
                %    .sigma0:     initial value of sigma
                %    .cmaes_opts: struct of options for cmaes, see cmaes documentation
                %                 or type cmaes to see list of options and default values
                
                % starting sigma
                optim_opts.insigma = .3*(parRanges(:,2) - parRanges(:,1));                 % starting sigma (this is default, could have left it blank)
                
                % other options
                optim_opts.LBounds  = parRanges(:,1);                                      % lower bounds of parameters
                optim_opts.UBounds  = parRanges(:,2);                                      % upper bounds of parameters
                optim_opts.PopSize  = 4 + floor(3*log(numParams));                         % population size (default)
                optim_opts.TolX       = 1e-4 * min(optim_opts.insigma);                    % stopping criterion on changes to parameters 
                optim_opts.TolFun     = 1e-4;                                              % stopping criterion on changes to fitness function
                optim_opts.TolHistFun = 1e-4;                                              % stopping criterion on changes to fitness function
                optim_opts.SaveFilename      = 'wf_ex_4_cmaesvars.mat';                    % output file of cmaes variables
                optim_opts.LogFilenamePrefix = 'wf_ex_4_';                                 % prefix for cmaes log-files
                % note that saving to ".mat" file of CMA-ES output is disabled for Octave
                % (lines 1795 and 1839 of my_cmaes.m) due to a bug that otherwise crashes the calibration.
                
                % Other useful options
                    % change to true to run in parallel on a pool of CPUs (e.g. on a cluster)
                optim_opts.EvalParallel = false;                                           
                
                % debugging options
                %optim_opts.MaxIter = 1;                                                  % just do 5 iterations, to check if it works
                %optim_opts.Seed =  12345;                                                % for reproducibility
                
                % initial parameter set
                par_ini = mean(parRanges,2);                                               % same as default value
                
                % Choose the objective function
                of_name      = 'OF_KGE_QnAET';                                             % Custom objective function combining streamflow and AET performance metrics.
                optim_fun    = 'my_cmaes';
                weights      = [1,1,1];                                                    % Weights for the three KGE components
                
                 
                %% 5. Calibrate the model
                % MARRMoT model objects have a "calibrate" method that takes uses a chosen
                % optimisation algorithm and objective function to optimise the parameter
                % set. See MARRMoT_model class for details.
                
                % first set up the model
                m.input_climate = input_climatology;
                %m.delta_t       = input_climatology.delta_t;                              % unnecessary if input_climate already contains .delta_t
                m.solver_opts   = input_solver_opts;
                m.S0            = input_s0;
            
            
                
            
                [par_opt,...                                                               % optimal parameter set
                    of_cal,...                                                             % value of objective function at par_opt
                    stopflag,...                                                           % flag indicating reason the algorithm stopped
                    output] = ...                                                          % other info about parametrisation
                              calibrate_QnAET(m,...                                        % a populated MARRMoT_model (Contain: .input_climate, .delta_t, .S0 & .solver_opts)
                                          obs, ...
                                          cal_idx, ...
                                          par_ini, ...
                                          optim_opts);
                
                % of_cal = Objective function:
                    %   Objective 1: Streamflow performance
                    %   Objective 2: Flux tower AET performance
                    %   Combined using equal weighting (0.5 / 0.5)
                
                %% 6. Evaluate the calibrated parameters on unseen data
                % Run the model with calibrated parameters, get only the streamflow
                [fluxOut, fluxInt, storeInt] = m.get_output([],[],par_opt);
                Q_sim = fluxOut.Q;
                AET_sim = fluxOut.Ea;
    
                if iModel == 1
                    SM_sim = storeInt.S1';   % soil_moisture store
                elseif iModel == 2 || iModel == 3
                    SM_sim = storeInt.S2';   % soil_moisture store
                end
                            
                % compute AET KGE function
                of_cal_aet = of_KGE_sqrt(aet_obs,AET_sim,obs.IntYM,cal_idx.AET);
                of_cal_q = of_bias_penalised_log(Q_obs ,Q_sim, cal_idx.Q, 'of_mean_hilo_root5_KGE');
            
                % BIAS 
                mean_Qobs = mean(Q_obs(warmup:n)); mean_Qsim = mean(Q_sim(warmup:n)); 
                mean_AETobs = mean(aet_obs(ETStart:ETEnd)); mean_AETsim = mean(AET_sim(ETStart:ETEnd));
                Q_bias = (1-((mean_Qobs - mean_Qsim) / mean_Qobs)) * 100;
                AET_bias = (1-((mean_AETobs - mean_AETsim) / mean_AETobs)) * 100;
                
                %% 8. save all the outputs
                outputs = [];
                % save the objective function value:
                outputs.(m.evap_eqn).of_cal_q   = of_cal_q;      %   - for streamflow (of_cal_q)
                outputs.(m.evap_eqn).of_cal_aet = of_cal_aet;    %   - for ET (of_cal_aet)
                outputs.(m.evap_eqn).of_cal     = of_cal;        %   - for ET & Q (of_cal_aet)
                % save the bias
                outputs.(m.evap_eqn).bias_q   = Q_bias;      %   - for streamflow (of_cal_q)
                outputs.(m.evap_eqn).bias_aet = AET_bias;    %   - for ET (of_cal_aet)
                % save the optimum parameters (par_opt)
                outputs.(m.evap_eqn).par_opt = par_opt;
        
        
                % save the simulated timeseries 
                sim_all = data_table(:,1);
                sim_all{:,'AET_sim'} = AET_sim;
                sim_all{:,'Q_sim'} = Q_sim;
                sim_all{:,'SM_sim'} = SM_sim;
                outputs.(m.evap_eqn).sim_timeseries = sim_all;
                    
                save(append('outputs/',thisCatch,'/',model_name,'_CalQnAET/',m.evap_eqn,'_outputs.mat'),'outputs')
        
        
                %% 7. Visualise the results
                              
                % Prepare a time vector
                t = data_table.date;
                
                % plot figure of streamflow (obs vs simulated)
                figure;
                plot(t, Q_obs)
                hold on
                plot(t, Q_sim)
                legend('Q_{obs}','Q_{sim}');
                title({'Streamflow calibration results', append('Cal Q & AET, ' ,m.evap_eqn)})
                ylabel('Streamflow [mm/d]')
                xlabel('Time [d]')
                txt_calQ  = sprintf('Q KGE = %.2f ',of_cal_q);
                txt_QBias = sprintf('Q Bias = %.2f ',Q_bias);
                text(t(warmup+365*5),40,txt_calQ,'fontsize',16,'HorizontalAlignment', 'left');
                text(t(warmup+365*5),37,txt_QBias,'fontsize',16,'HorizontalAlignment', 'left');
                set(gca,'fontsize',16);
                xlim([t(1), t(end)])
                datetick;
                ylim([0,45])
                set(gca,'TickLength',[0.005,0.005])
                savefig(append('outputs/',thisCatch,'/',model_name,'_CalQnAET/',m.evap_eqn,'_Q.fig'))
                
                % plot figure of ET (obs vs simulated)
                figure;
                plot(t, aet_obs)
                hold on
                plot(t, AET_sim)
                legend('AET_{obs}','AET_{sim}');
                title({'Actual ET calibration results', append('Cal Q & AET, ' ,m.evap_eqn)})
                ylabel('ET [mm/d]')
                xlabel('Time [d]')
                txt_calAET  = sprintf('AET KGE = %.2f ',of_cal_aet);
                txt_AETBias = sprintf('AET Bias = %.2f ',AET_bias);
                text(t(ETStart+365*2),7,txt_calAET,'fontsize',16,'HorizontalAlignment', 'left');
                text(t(ETStart+365*2),6.5,txt_AETBias,'fontsize',16,'HorizontalAlignment', 'left');
                set(gca,'fontsize',16);
                xlim([t(ETStart), t(ETEnd)])
                datetick;
                ylim([0,8])    
                set(gca,'TickLength',[0.005,0.005])
                
                savefig(append('outputs/',thisCatch,'/',model_name,'_CalQnAET/',m.evap_eqn,'_ET.fig'))

            end
        
        end
    end
    
end

%% functions

function [RunReq, OutFilePath] = IsRunReq(ThisETeqn, model_name, thisCatch)

    % created July 2020 by Keirnan Fowler, University of Melbourne
    % Utility function developed by Keirnan Fowler (University of Melbourne)
    % to manage batch model runs and prevent duplicate executions.
    
    % define file path
    OutFilePath   = 'outputs/'+thisCatch+'/'+model_name+'_CalQnAET/'+ ThisETeqn + '_outputs.mat';
    OutFigurePath = 'outputs/'+thisCatch+'/'+model_name+'_CalQnAET/'+ ThisETeqn + '_ET.fig';

    % infer whether this run is required based on existence of output file.
    FileExist = exist(OutFilePath, 'file'); 
    FigureExist = exist(OutFigurePath, 'file'); 
    if FileExist == 0 
        % it doesn't exist, and thus hasn't been run  
        CreateDummyFile = true; RunReq = true; 
    elseif FigureExist ~= 0 
        % the figure does exist - so no need to check dummy file the run
        % has successfully completed
        CreateDummyFile = false; RunReq = false; 
    else
        % File exists.  Get its size.
        s = dir(OutFilePath); FileSize_KB = s.bytes/1000; 
        
        if FileSize_KB > 1
            % we infer from its size that the run has successfully completed.
            CreateDummyFile = false; RunReq = false; 
        else
            % we infer that it is just the dummy file created by a
            % earlier instantiation of this code.  Get its age.
            FileAge = now() - s.datenum;
            
            if FileAge < (70/24) % (hours want/ 24)
                % it was created recently, so we infer that the run is in progress.  
                % note, both optimisers create a new one of these every
                % generation, so we are assuming a generation of model evaluations 
                % (here 200) can occur in an hour or less.  
                CreateDummyFile = false; RunReq = false; 
            else
                % it is old.  We infer that the run crashed.  For the 
                % Fortran-based models, this is common due to instabilities 
                % in the Matlab-Fortran interface. 
                CreateDummyFile = true; % create a dummy file with a newer timestamp.
                RunReq = true; 
            end 
        end
    end
    
    if CreateDummyFile
        FileID = fopen(OutFilePath, 'w+'); fprintf(FileID, '%s\r\n', 'In progress'); fclose(FileID);
    end
    
end
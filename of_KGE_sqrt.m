function [val,val_daily,val_monthly,c,idx,w] = of_KGE_sqrt(obs, sim,int_YM,idx, w)
    % of_KGE Calculates Kling-Gupta Efficiency of simulated streamflow (Gupta
    % et al, 2009). Ignores time steps with negative flow values.
    
    % Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
    % This file is part of the Modular Assessment of Rainfall-Runoff Models
    % Toolbox (MARRMoT).
    % MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
    % WARRANTY. See <https://www.gnu.org/licenses/> for details.
    
    % In:
    % obs       - time series of observations       [nx1]
    % sim       - time series of simulations        [nx1]
    % idx       - optional vector of indices to use for calculation, can be
    %               logical vector [nx1] or numeric vector [mx1], with m <= n
    % w         - optional weights of components    [3x1]
    %
    % Out:
    % val       - objective function value          [1x1]
    % c         - components [r,alpha,beta]         [3x1]
    % idx       - indices used for the calculation
    % w         - weights    [wr,wa,wb]             [3x1]
    %
    % Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009). 
    % Decomposition of the mean squared error and NSE performance criteria: 
    % Implications for improving hydrological modelling. Journal of Hydrology, 
    % 377(1–2), 80–91. https://doi.org/10.1016/j.jhydrol.2009.08.003
    
    %% Check inputs and select timesteps
    if nargin < 3
        error('Not enugh input arguments')    
    end
    
    if nargin < 4; idx = []; end
    [sim, obs, idx] = check_and_select(sim, obs, idx);
    
    %% Set weights
    w_default = [1,1,1];          % default weights
    
    % update defaults weights if needed  
    if nargin < 5 || isempty(w)
        w = w_default;
    else
        if ~min(size(w)) == 1 || ~max(size(w)) == 3                            % check weights variable for size
            error('Weights should be a 3x1 or 1x3 vector.')                    % or throw error        
        end
    end   
    
    %% calculate components
    obs_sqrt = obs.^.5; 
    sim_sqrt = sim.^.5;
    
    
    c(1) = corr(obs_sqrt,sim_sqrt);                                             % r: linear correlation
    c(2) = std(sim_sqrt)/std(obs_sqrt);                                         % alpha: ratio of standard deviations
    c(3) = mean(sim_sqrt)/mean(obs_sqrt);                                       % beta: bias 
    
    %% calculate value
    val_daily = 1-sqrt((w(1)*(c(1)-1))^2 + (w(2)*(c(2)-1))^2 + (w(3)*(c(3)-1))^2);    % weighted KGE
    
    %% 2) Calculate objective function value at monthly scale
    
    % 1) Calculate number of integers to aggregate daily values to monthly
     Int_yearmonth = int_YM(idx);
     
     % Find unique yearmonth integer values
     Unique_int = unique(Int_yearmonth); NumUniqueInt = size(Unique_int,1);
     obs_month = []; sim_month = [];
     for iUniqueInt = 1:NumUniqueInt
         Ind_thisyearmonth = [];
         ThisInt = Unique_int(iUniqueInt);
         
         % Find positions of the given integer to select observed and predicted values for the given month and year
         Ind_thisyearmonth = find(Int_yearmonth==ThisInt);
         
         % Monthly observed and simulated AET
         obs_month(iUniqueInt) = nansum(obs(Ind_thisyearmonth),'all');
         sim_month(iUniqueInt) = nansum(sim(Ind_thisyearmonth),'all');
         
     end
     
     %% calculate monthly compoenents
     obs_monthly_sqrt = obs_month.^.5;
     sim_monthly_sqrt = sim_month.^.5;
     
     %note: column input vectors are needed in e(1),e(2) calculation
     e(1) = corr(obs_monthly_sqrt',sim_monthly_sqrt');                         % r: linear correlation
     e(2) = std(sim_monthly_sqrt')/std(obs_monthly_sqrt');                     % alpha: ratio of standard deviations
     e(3) = mean(sim_monthly_sqrt')/mean(obs_monthly_sqrt');                   % beta: bias 
    
     %% calculate monthly value
     val_monthly = 1-sqrt((w(1)*(e(1)-1))^2 + (w(1)*(e(2)-1))^2 + (w(3)*(e(3)-1))^2);    % weighted KGE
    
     %% Mean of the daily and monthly objective function values
     val = (val_daily+val_monthly)./2;

 
end
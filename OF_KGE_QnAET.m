function val = OF_KGE_QnAET(obs, sim,idx_Q,obsAET,simAET,int_YM,idx_AET,varargin)

    % OF_KGE_QnAET Combined objective function for calibration against
    % streamflow and actual evapotranspiration (AET).
    % This objective function combines:
    %   1) A bias-penalised log-transformed KGE objective for streamflow.
    %   2) A daily/monthly square-root transformed KGE objective for AET
    %      (implemented in of_KGE_sqrt).
    % The final objective value is the equally weighted mean of the two
    % metrics:
    %   val = 0.5 * (val_1 + val_2)
    %   % Inputs:
    %   obs      - Observed streamflow
    %   sim      - Simulated streamflow
    %   idx_Q    - Index used for streamflow objective evaluation
    %   obsAET   - Observed AET
    %   simAET   - Simulated AET
    %   int_YM   - Integer year-month identifier used to aggregate daily AET
    %              to monthly values
    %   idx_AET  - Index used for AET objective evaluation
    % See of_bias_penalised_log / of_KGE_sqrt for details of the objective functions.

    val_1 = of_bias_penalised_log(obs ,sim,idx_Q,'of_mean_hilo_root5_KGE');
	val_2 = of_KGE_sqrt(obsAET,simAET,int_YM,idx_AET);               % int_YM - assigned interger value each day for a
                                                                            %     given year month to calculate monthly values   
    % mean value of the objective function
    val = 0.5*(val_1+val_2);
end
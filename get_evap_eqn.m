
function [out] = get_evap_eqn(evap_eqn, S, Ep, dt, p1, Smax, p2)

% get_evap_eqn Wrapper function for selecting evapotranspiration formulations
% within the MARRMoT framework.
%
% This function was developed to enable comparison of multiple
% evapotranspiration formulations within a common conceptual hydrological
% modelling structure. 
%
% Developed for:
% Burns, G., Fowler, K., Peel, M., and Stephens, C.: A systematic evaluation
% of 15 actual evapotranspiration formulations within conceptual hydrological
% models, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-3122, 2025. 
%
% Based on the MARRMoT framework:
% Trotter, L., Knoben, W. J. M., Fowler, K. J. A., Saft, M., and Peel,
% M. C.: Modular Assessment of Rainfall–Runoff Models Toolbox (MARRMoT) v2.1:
% an object-oriented implementation of 47 established hydrological models for
% improved speed and readability, Geosci. Model Dev., 15, 6359–6369, 
% https://doi.org/10.5194/gmd-15-6359-2022, 2022. 
%
% Inputs:
%   evap_eqn - Evapotranspiration formulation ('evap_1', etc.)
%   S        - Soil water storage [mm]
%   Ep       - Potential evapotranspiration [mm d^-1]
%   dt       - Model time step [d]
%   p1,p2    - Dimensionless scaling parameters (0–1)
%   Smax     - Maximum soil water storage [mm]
%
% Notes:
%   Additional calibration parameters (p1, p2) are sampled in a common 0–1 space and
%   transformed within this function to the parameter ranges required by
%   each evapotranspiration formulation. Some formulations additionally
%   require secondary storage terms (S2), which are fixed to zero in this
%   implementation to match the single-store model structure.


S1 = S;
S1max = Smax;
S2min = 0.01;



switch evap_eqn
    case 'evap_1', [out] = evap_1(S,Ep,dt);
    case 'evap_2', p1 = p1*20; % Scale p1 from [0,1] to [0,20] 
        [out] = evap_2(p1,S,Smax,Ep,dt);  
    case 'evap_3', [out] = evap_3(p1,S,Smax,Ep,dt);
    case 'evap_4', [out] = evap_4(Ep,p1,S,p2,Smax,dt);

    case 'evap_6', [out] = evap_6(1,p1,S,Smax,Ep,dt);
    case 'evap_7', [out] = evap_7(S,Smax,Ep,dt);
    case 'evap_8', p1 = p1*Smax; % Scale threshold storage parameter
        S2 = 0; % No secondary storage in this implementation
        [out] = evap_8(S1,S2,1,p1,Ep,dt); % Forest fraction fixed at 1

    case 'evap_11', [out] = evap_11(S,Smax,Ep);

    case 'evap_13', [out] = evap_13(p1,p2,Ep,S,dt);
    case 'evap_14', S2 = 0; [out] = evap_14(p1,p2,Ep,S1,S2,S2min,dt);  % S2 == 0
    case 'evap_15', S2 = 0; [out] = evap_15(Ep,S1,S1max,S2,S2min,dt);  % S2 == 0
    case 'evap_16', S2 = 0; [out] = evap_16(p1,S1,S2,S2min,Ep,dt);     % S2 == 0
    case 'evap_17', [out] = evap_17(p1,S,Ep);

    case 'evap_19', [out] = evap_19(p1,p2,S,Smax,Ep,dt);
    case 'evap_20', p1 = p1*20; p2 = p2*Smax; [out] = evap_20(p1,p2,S,Smax,Ep,dt);
    case 'evap_21', p1 = p1*Smax; [out] = evap_21(p1,p2,S,Ep,dt);  % p1 = (0-1)*(Smax-S2)
    case 'evap_22', p1 = p1*Smax; p2 = p2*p1; [out] = evap_22(p1,p2,S,Ep,dt);
    case 'evap_23', [out] = evap_23(1,p1,S,Smax,Ep,dt); % Forest fraction fixed at 1
end

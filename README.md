# MARRMoT Multi-Objective ET Calibration Framework

## Overview

This repository contains MATLAB code used in:

**Burns et al. (2025)**
*A systematic evaluation of 15 actual evapotranspiration formulations within conceptual hydrological models*

The code extends the **MARRMoT hydrological modelling framework** to evaluate multiple evapotranspiration (ET) formulations within conceptual rainfall–runoff models. Models are calibrated using a **combined objective function combining streamflow and actual evapotranspiration (AET)** performance.

---

## Key Features

* Multi-model support:

  * GR4J
  * SIMHYD
  * VIC
* Evaluation of multiple ET formulations (15+ variants)
* Multi-objective calibration using:

  * Streamflow (KGE with bias penalty)
  * Flux tower AET (KGE on transformed fluxes)
* Integration with MARRMoT model structures

---

## Repository Contents

### Core calibration and workflow scripts

* **calibrate_3model_7catch.m**
  Main driver script that runs calibration across multiple models and catchments.

* **calibrate_QnAET.m**
  Multi-objective calibration wrapper for streamflow + AET.

---

### Objective functions

* **OF_KGE_QnAET.m**
  Combined objective function used for calibration.
  Equal weighting of:

  * Streamflow performance
  * AET performance

* **of_KGE_sqrt.m**
  Kling–Gupta Efficiency (KGE) implementation applied to square-root transformed fluxes, including daily and monthly aggregation.

* **of_bias_penalised_log.m**
  Streamflow objective function based on log-transformed KGE with bias penalisation.

---

### Evapotranspiration formulation interface

* **get_evap_eqn.m**
  Wrapper function that selects and parameterises different ET equations.
  Maps normalised calibration parameters (0–1) to formulation-specific parameter ranges.

---

## Requirements

* MATLAB 
* MARRMoT toolbox 
* Catchment forcing data (precipitation, temperature, PET, observed streamflow, and AET)

---

## How to Run

1. Add MARRMoT and repository folders to MATLAB path:

   ```matlab
   addpath(genpath('2_MARRMoT'))
   addpath(genpath('calibration_codes'))
   ```

2. Prepare input data in the required format (see `calibrate_3model_7catch.m`).

3. Run the main calibration script:

   ```matlab
   calibrate_3model_7catch
   ```

4. Outputs will be saved by catchment, model, and ET formulation in the `outputs/` directory.

---

## Method Summary

Each model is calibrated using a combined objective function:

* **Objective 1:** Streamflow performance (KGE with bias penalty)
* **Objective 2:** Flux tower AET performance (KGE on transformed fluxes)

Final objective:

[
OF = 0.5 \cdot OF_{Q} + 0.5 \cdot OF_{AET}
]

---

## Notes

* ET formulations are selected via `get_evap_eqn.m`
* AET evaluation is restricted to periods with available flux tower observations
* A 5-year warm-up period is used to reduce sensitivity to initial conditions

---

## Attribution

Developed as part of:

Burns, G. et al. (2025)

Built on the **MARRMoT hydrological modelling framework**:
Knoben, W.J.M., et al. (2019) & Trotter et al., (2022)

Individual objective functions and ET formulations are based on prior work within MARRMoT and associated literature, as referenced in the corresponding function headers.

---

## License / Use

This code is provided for academic and research use. Please cite the associated publication when using or adapting this repository.

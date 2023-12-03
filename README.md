# Abnormal_Grain_Growth_MATLAB
Stored energy driven grain growth codes


## Parameters
Simulation and material parameters are set in ```/codes/Params.m``` <br />
Timestepping is set in grain_growth.m and ```initial_condition.m``` <br />

## Running AGG Simulation

In ```/codes/```, <br />
Run ```initial_condition.m``` <br />
Run ```rho_grain_id.m``` or a custom initialization for the stored energy field <br />
Run ```grain_growth.m``` <br />

Figures will be saved if saveFigs==true. They can also be created separately later using ```make_figures.m``` if saveFiles==true.

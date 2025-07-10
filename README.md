# aew_mse
Repository for all codes used to generate a moist static energy budget for African Easterly Waves. Python scripts created by A. Aiyyer &amp; A. Thornton.
Recommended order:
  1. era5_regrid.py: downloads and regrids era5 data to 0.5 deg grid spacing for specified subset/domain.
  2. 2D_flux_regrid_6h_updated.py: downloads, calculates, and regrids era5 flux data to 0.5 deg grid spacing, 6 hourly, for specified subset/domain.
  3. q_term.py: calculates and saves total heating in a single term for the MSE budget equation, reflects atmospheric gains and losses.
  4. static_energy.py: calculates and saves both DSE and MSE from era5 data.
  5. vert_integration_static_energy.py: vertically integrates and saves calculated DSE and MSE.
  6. jas_means.py: calculates and saves the JAS mean for each variable for specified range of years.
  7. MSE_Budget_Part1.py: calculates and saves all advection terms of the MSE budget equation.
  8. filter_era5.py: filters and saves each variable in MSE budget equation using 2-10 BP Lancozs filter; saved as yearly files for each variable.

How to use:
- Submit batch jobs for each variable separately over a long period of time, use .py scripts.
- To get a quick look at the process, use a smaller subset of data, and view detailed breakdown of code, use the .ipynb files. 

Post-processing
- in post processing folder, find vector projections for maintenance and propagation, showing budget closure

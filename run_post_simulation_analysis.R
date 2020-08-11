### THis is the master script to run post-simulation analysis
### for GDAY simulations
### Mingkai Jiang (m.jiang@westernsydney.edu.au)

#### call functions
#### clear wk space
rm(list=ls(all=TRUE))

#### Source functions and packages
source("analysis/prepare.R")


############################## GDAY simulation output ############################## 
### Plot historic validation
plot_historic_validation()

### Plot current future trajectory
plot_current_future_trajectory()

############################## Met data ############################## 
### Plot annual met data
#plot_annual_met_data()

### Plot daily met data over the transition periods
#plot_daily_met_data()


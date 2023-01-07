#
# Replication code
# Athey, Grabarz, Luca, and Wernerfelt (2023)
#

# This file loads our simulated dataset and then runs the main analysis of the
# paper on the simulated data. Since all the data are simulated, the results
# do not match the ones in the paper, but we hope this provides a sense of how
# the analysis was done.

rm(list=ls())

# Load relevant packages
library(data.table)
library(dplyr)
library(survey)
library(meta)


###############################################################################
# Load data
###############################################################################

# Specify your own path to the downloaded .rdata file here
load('/Users/nilsw/Downloads/COVID_replication.rdata')
options(warn=-1) 

# Two data tables are in COVID_replication.rdata:
# *rdt* is the simulated response-level data table
# *wdt* is the simulated data table for weights
# The variable names are self-explanatory, but see the main text if confused

###############################################################################
# Function definitions: Used to analyze each individual experiment
###############################################################################

# Analyze individual experiment
experiment_level_regression_hte <- function(data, pop_proportions) {
  
  # Create survey design
  svy.data.unweighted <- svydesign(ids=~1, data=data)
  
  # Get the marginal probabilities for the variables that we want to weight the data by. 
  pop_proportions$Freq <- pop_proportions$weight_vals*nrow(data)
  
  # Get the breakdown distribution
  breakdown_dist <- pop_proportions[,.(breakdown_key, Freq)]
  
  # Will use this in raking below
  pop_test_control <- data.table(
    test_group=c(1,0), 
    Freq=c(0.5*nrow(data), 0.5*nrow(data))
  )
  
  # Use the rake function in survey package to weight current data by pop values
  data.svy.rake <- rake(design = svy.data.unweighted,
                        sample.margins = list(~breakdown_key, ~test_group), 
                        population.margins = list(breakdown_dist, pop_test_control)
  )
  
  # Trim the weights
  data.svy.rake.trim <- trimWeights(data.svy.rake, lower=0.3, upper=3, strict=TRUE)
  
  # Put weights in a column
  data[, weights := weights(data.svy.rake)]
  data[, weights_trim := weights(data.svy.rake.trim)]
  
  # Model specification
  formula_dr = 'response_score ~ test_group+
    demeaned_age2534+demeaned_age3544+demeaned_age4554+demeaned_age5564+demeaned_age65p+
    test_group*demeaned_age2534+test_group*demeaned_age3544+test_group*demeaned_age4554+
    test_group*demeaned_age5564+test_group*demeaned_age65p+demeaned_male+test_group*demeaned_male+ 
    demeaned_ectr + test_group * demeaned_ectr+ demeaned_ecvr + test_group * demeaned_ecvr'
  
  # Fit model
  dr_mod <- svyglm(formula = formula_dr, design = data.svy.rake.trim)
  
  # Grab main treatment effect coefficient
  be_coef <- summary(dr_mod)$coefficients[2,1] # same
  
  # Grab treatment effect SE 
  be_coef_se <- summary(dr_mod)$coefficients[2,2] # same
  
  # Grab p-value
  be_p_val <- summary(dr_mod)$coefficients[2,4] # same
  
  output <- c(be_coef, be_coef_se,be_p_val)
  
  return(output)
}



# Function for looping analysis across all experiments
execute_experiment_level_regression_hte <- function(unique_experiments, pop_props, dataset) {
  
  # Prepare output: Create empty output data table
  output_dt <- dataset %>%
    dplyr::select(study_id, cell_id, experiment_id, question_short_name) %>%
    dplyr::filter(experiment_id %in% unique_experiments) %>%
    unique() %>%
    data.table()
  
  output_dt[,reweighted_lift_coef := 0.0]
  output_dt[,reweighted_lift_coef_se := 0.0]
  output_dt[,reweighted_lift_p_val := 0.0]
  
  
  # Now loop through experiments and execute           
  for (i in 1:length(unique_experiments)) {
    
    # Define experiment
    exp_i <- unique_experiments[i]
    
    # Filter data
    exp_data_i <- dataset[experiment_id == exp_i]
    
    # Filter weight data
    weight_data_i <- pop_props[cell_id == unique(exp_data_i$cell_id)]
    
    # Execute function to perform regression
    output <- experiment_level_regression_hte(exp_data_i, weight_data_i)
    
    output_dt[experiment_id==exp_i,5:7] <- as.list(output)
    
  }
  
  return(output_dt)
}

###############################################################################
# Run analysis of each experiment
###############################################################################

# Declare list of experiments
experiment_list <- unique(rdt$experiment_id) 

# Run analysis
experiment_level_output_hte <- execute_experiment_level_regression_hte(experiment_list, wdt, rdt)


###############################################################################
# Take the results of each experiment, run meta-analysis, and print out results
###############################################################################


# Overall and Outcome-specific Fixed Effects Results

# Overall meta-analysis
meta_overall <- metagen(TE = reweighted_lift_coef, 
                        seTE = reweighted_lift_coef_se,
                        studlab = experiment_id, 
                        data = experiment_level_output_hte)

# Print out results
results_print <- paste('  Overall Fixed Effects Result (TE, SE, p-value):\n', 
                       meta_overall$TE.fixed, 
                       meta_overall$seTE.fixed, 
                       meta_overall$pval.fixed, sep='  ')
cat(results_print)


# Outcome-specific: redo on each specific outcome question
category_meta <- function(data, category) {
  
  # Subset data to just this category
  this_data <- data %>%
    dplyr::filter(question_short_name==category)
  
  # Perform IVW using metagen
  meta <- metagen(TE = reweighted_lift_coef, 
                  seTE = reweighted_lift_coef_se,
                  studlab = experiment_id, 
                  data = this_data) 
  
  # Fixed effects results
  meta_treatment_effect <- meta$TE.fixed
  meta_se <- meta$seTE.fixed
  meta_p_val <- meta$pval.fixed
  print(paste0('Results for ', category))
  print(paste0('Treatment Effect: ',meta_treatment_effect))
  print(paste0('Treatment Effect SE: ',meta_se))
  print(paste0('Treatment Effect p-val: ',meta_p_val))
  print('')
}

# Use function to perform analyses pooled by category
category_list <- sort(unique(experiment_level_output_hte$question_short_name))

# Call function to run on each category
for (this_category in category_list) {
  category_meta(experiment_level_output_hte, this_category)
}


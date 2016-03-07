#####SPRTwoState package for analyzing and simulating SPR data
#  this R package is used to do the analysis under the general 
#  two model, especially for antigen and antibody interactions
#  
#  Version 0.1
#  For now, 1)it will only fit the unified two state model
#  2)It could also generalized the two state model for 
#    multiple conformations of antibodies
#  3) assume only antibodies (ligands) have multiple conformations, but
#     not antigens (analytes).
#
#### Developed By Feng @ BU 2016. All right reserved###

#Below are roxygen2 comments

##########Dependency###########
##  library(MASS)
####################
#this is where the files to be include if necessary
#' @include SPRTwoState.R

#'
#' @title SPR Analyzer for the Two State model
#'
#' @description An R package to characterize the Antibody Antigen 
#'		interactions under the Two State model by 
#'		surface plasmon resonance
#'
#' @details This package is defined to run analysis of SPR data 
#'		under the two state model to estimate mean equilibrium
#'		dissociation constant as well as the mean dissociation rate
#'		constant. It takes in the exported text input data. Output 
#'		is estimated by non-linear regression. The modelling can be
#'		generilized to multi-state conditions. 
#'		This package can also provide functions to run simulations
#'		for different models, such as Langmuir model, induced fit
#'		model, conformation selection model and two state model. 
#'
#' 		Please refer to the vignettes to see details.
#' 		
"_PACKAGE"

#' @editor Dr. Martin C. Arostegui \email{martin.arostegui@@whoi.edu}
#' calculate_attenuation_eTUFF is a modified version of calculate_attenuation from the 'trawllight' package.
#' The edits were to accommodate the different diffuse attenuation coefficient, K, calculation when using tag manufacturer raw light units.
#' Edit 1: This uses Eq. 3 from Teo et al. 2009 Aquatic Biology 5:195-207 to calculate K and includes a modification:
#'    In the Teo paper they always calculated K at 1-m intervals, but I edited the code here to automatically adjust the answer by the depth bin interval size.
#'    This ensures the K output is always standardized back to units of m^-1.
#' Edit 2: @param tag.manufacturer *New param: options are ("Wildlife Computers") and ("Lotek")
#'    Each manufacturer has a unique 'D' value in Eq. 3 from Teo et al., which is the number of raw units per order of magnitude change in light level (units/decade).
#'    Wildlife Computers 'D' value of 20 can be cited with Jaud et al. 2012 PLoS ONE 7:e47444 (see slope of equation in Fig. 1).
#'    Lotek 'D' value of 32 can be cited with Ekstrom 2004 Memoirs of the National Institute of Polar Research, Special Issue, 58:210–226
#' Edit 3: Increased the minimum number of binned light-depth points that the loess must have available to model a profile (at least 7)
#' Edit 4: Added a QC for the loess light-depth predictions that excludes profiles with any non-positive K values at 1-m intervals
#' Edit 5: Changed name of column "trans_llight" to "raw_llight" and eliminated associated light transformations
#' Edit 6: @param plot *New param: options are ("T") and ("F"). T yields plot of loess light prediction (and input data) and derived K values (at bin intervals)
#' Edit 7: Added x$quality filter at beginning of function to exclude input data not meeting bin continuity standards
#'         This function should directly take the output of filter_stepwise_eTUFF, where the quality value is assigned (1 is good, -999 is bad)
#' Edit 8: @param bin.gap *New param: allowable gap (in meters) for binned light data used to generate attenuation profiles (must be same value as used in preceding call to filter_stepwise_eTUFF)
#' Edit 9: Adds '$meta' entry to output list containing tag manufacturer, bin size (m), allowable bin gap (m), minimum allowable depth range of profile (m) from filter_stepwise_eTUFF. and mean absolute error of loess light profile
#' Edit 10: @param depth.subset *Optional param: vector of length 2 containing (min, max) depths of desired profile subset. Input depths (m) must correspond to $attenuation$depth in data frame 'output_dfs'

#' Please see GitHub page below and *cite Sean Rohan's work* if using my edit of his function:
#' https://github.com/sean-rohan-NOAA/trawllight
#' Rohan et al. 2020 NOAA Technical Memorandum NMFS-AFSC-403
#' Rohan et al. 2021 Progress in Oceanography 194:102554
#' 
#'  Diffuse attenuation coefficient of downwelling irradiance
#'
#' \code{calculate_attenuation} fits a loess model between depth and light using an AICc-based span selection model adapated from \code{fANCOVA::loess.as}, then estimates the first derivative of the resultant model slope to approximate the diffuse attenuation coefficient of downwelling irradiance (vertical attenuation coefficient).
#'
#' @param x Data frame containing depth and light for a single cast.
#' @param loess.criterion Criterion for choosing the most parsimonious model. Options are bias-corrected Akaike's Information Criterion ("aicc") or generalized cross-validation ("gcv").
#' @param loess.degree Degrees for loess model. Default = 1.
#' @param kz.binsize Depth interval for estimating instantaneous diffuse attenuation coefficient of downwelling irradiance. Default = 0.2.
#' @param min.range Minimum range of depths necessary for model fitting. Default = 10.
#' @param light.predict Logical indicating whether predicted values for light should be returned.
#' @param ... Additional arguments passed to loess fitting function
#' @return Returns a list containing three data frames: \code{attenuation} contains depth and fitted values of vertical diffuse attenuation coefficient, \code{loess.fit} contains the model summary statistics, and \code{fit_residuals} contains model fit residuals.
#' @references Hurvich, C.M., Simonoff, J.S., and Tsai, C.-L. 1998. Smoothing parameter selection in nonparametric regression using an improved Akaike information criterion. J. R. Stat. Soc. B 60(2): 271-293.
#' @references Xiao-Feng Wang (2010). fANCOVA: Nonparametric Analysis of Covariance. R package version 0.5-1. https://CRAN.R-project.org/package=fANCOVA
#' @author Sean Rohan \email{sean.rohan@@noaa.gov}
#' @export

calculate_attenuation_eTUFF.BGC <- function(x,
                                        light.col = "raw_llight",
                                        depth.col = "cdepth",
                                        loess.criterion = "aicc",
                                        loess.degree = 1,
                                        kz.binsize = 5,
                                        min.range = 50,
                                        bin.gap = NULL,
                                        light.predict = T,
                                        tag.manufacturer = NULL,
                                        depth.subset = NULL,
                                        plot=F,
                                        ...) {
  
  names(x)[which(names(x) == light.col)] <- "raw_llight"
  names(x)[which(names(x) == depth.col)] <- "cdepth"
  
  # Assign 'D' value according to tag manufacturer
  if (tag.manufacturer == "Wildlife Computers") {
    D <- 20
  } else if (tag.manufacturer == "Lotek") {
    D <- 32
  } else if (tag.manufacturer == "BGC-Argo") {
    D <- 10
  } else {
    warning("calculate_attenuation_eTUFF: Your tag manufacturer is not yet recognized by this function, which will need to be updated to include its D value.")
    return(NA)
  }
  
  # Remove profiles with only a small portion of the water column sampled
  if (any(x$quality == 1)) { # Do not fit if input profile data from filter_stepwise_eTUFF does not meet bin continuity standards
    if((max(x$cdepth) - min(x$cdepth)) >= min.range) { # Do not fit if depth range is < min.range
      if(length(unique(x$cdepth)) > 4) { # Cannot fit loess to fewer than five data points
        
        # Fit loess model
        N_depths <- seq(min(x$cdepth), max(x$cdepth), kz.binsize)
        profile_light_loess <- loess.as2(x = x$cdepth, y = x$raw_llight, criterion = loess.criterion, degree = loess.degree)#, ...
        
        # Adjust for overfitting caused by omitted data
        k <- 3
        while(profile_light_loess$s == Inf) {
          profile_light_loess <- loess.as2(x = x$cdepth, y = x$raw_llight, criterion = loess.criterion,
                                           degree = loess.degree,
                                           min.bins = k+1)
        }
        
        # Quality control K values of loess prediction at 1-m intervals (all must be positive to be approved)
        light_fit_check <- predict(profile_light_loess, newdata = seq(min(x$cdepth), max(x$cdepth), 1))
        
          kdz_check = (1 - 10^(diff(light_fit_check) / D)) / 1
        
        
        if (any(kdz_check <= 0)) {
          warning("calculate_attenuation_eTUFF: Did not build profile. Attenuation coefficients (K) from 1-m intervals of loess prediction include zero or negative values")
          return(NA)
        } else {
          light_fit <- predict(profile_light_loess, newdata = N_depths)
          
          # Output data
            output <- data.frame(depth =  N_depths[1:(length(N_depths)-1)] + kz.binsize / 2,
                                kdz = (1 - 10^(diff(light_fit) / D)) / kz.binsize) # The K calculation for tag manufacturer raw light units
          
          # Output residuals
          resids <- data.frame(residual = residuals(profile_light_loess),
                               raw_llight = x$raw_llight,
                               predicted_llight = predict(profile_light_loess),#
                               cdepth = x$cdepth)
          
          # Mean Absolute Error (mae) of loess light profile
          # If subsetting the profile, the 'mae' is of the entire profile from which the subset is taken
          # This ensures the 'mae' and subsequent profile QC filtering is consistent regardless of the subset taken
          mae <- sum(abs(resids$residual)) / length(resids$residual)
          
          if (light.predict) {
            output$predict_light <- predict(profile_light_loess, newdata = output$depth)
          }
          
          # Output metadata
          meta <- data.frame(tag.manufacturer = tag.manufacturer,
                             bin.size = kz.binsize,
                             bin.gap = bin.gap,
                             min.range = min.range,
                             mae = mae)
          
          if (plot) {
            plot(x=light_fit_check,y=seq(min(x$cdepth), max(x$cdepth), 1)*-1,type="l",xlim = c(50,205),ylim = c(-450,0))
            points(y = x$cdepth*-1, x = x$raw_llight)
            par(new=T)
            plot(y = output$depth*-1, x = output$kdz,col=2,lwd=2,type="l",axes=F,xlim = c(0,.1),ylim = c(-450,0))
            axis(side=3)
          }
          
          # Subset attenuation and light profiles by desired min and max depths (these values must correspond to the K profile depths)
          if (is.null(depth.subset)){
            
          } else if (depth.subset[1] < min(output$depth) | depth.subset[2] > max(output$depth)){
            warning("calculate_attenuation_eTUFF: Depth subset mismatch. Input depth bound(s) not within available profile depth range")
            return(NA)
          } else {
            output <- output[which(output$depth >= depth.subset[1] & output$depth <= depth.subset[2]),]
            resids <- resids[which(resids$cdepth >= (depth.subset[1] - kz.binsize/2) & resids$cdepth <= (depth.subset[2] + kz.binsize/2)),]
          }
          
          output_dfs <- list(attenuation = output, fit_residuals = resids, meta = meta)
          return(output_dfs)
        }
        
      } else {
        warning("calculate_attenuation_eTUFF: Cannot fit loess. Cast has fewer than seven Unique depth bins and/or does not meet bin continuity standards")
        return(NA)
      }
    } else {
      warning("calculate_attenuation_eTUFF: Did not fit loess. Total depth range of cast < min.range.")
      return(NA)
    }
  } else {
    warning("calculate_attenuation_eTUFF: Did not fit loess. Input binned data from filter_stepwise_eTUFF does not meet bin continuity standards (quality == -999).")
    return(NA)
  }
}
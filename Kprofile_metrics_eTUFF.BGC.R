#' @author Dr. Martin C. Arostegui \email{martin.arostegui@@whoi.edu}
#' \code{Kprofile_metrics_eTUFF} calculates numerous attenuation (K) profile metrics and aggregates them into a data frame
#'
#' @param profile.data List of the output from calculate_attenuation_eTUFF (includes sublists '$attenuation', '$fit_residuals', and '$meta')
#' @param tag.manufacturer Options are ("Wildlife Computers") and ("Lotek")
#'    Each manufacturer has a unique 'D' value in Eq. 3 from Teo et al., which is the number of raw units per order of magnitude change in light level (units/decade).
#'    Wildlife Computers 'D' value of 20 can be cited with Jaud et al. 2012 PLoS ONE 7:e47444 (see slope of equation in Fig. 1).
#'    Lotek 'D' value of 32 can be cited with Ekstrom 2004 Memoirs of the National Institute of Polar Research, Special Issue, 58:210–226
#' @return Data frame with all K profile metrics
#' @export

#' Useful equations for converting among raw/absolute light and changes in orders of magnitude vs proportion
#' Aabs = Babs * 10^n where 'Aabs' is light (absolute) at deeper depth, 'Babs' is light (absolute) at shallower depth, n is change in orders of magnitude
#' Rearranged: n = log10(Aabs/Babs)
#' Also: n = (Arel - Brel)/D where 'Arel' is light (raw relative) at deeper depth,'Brel' is light (raw relative) at shallower depth, D is the units (raw relative)/decade
#' LLlost = 1 - 10^n where Ktotal is the proportion of light (absolute) lost over the depth interval (depthB, depthA)
#' 
Kprofile_metrics_eTUFF.BGC <- function(profile.data) {
  
  # Assign 'D' value according to tag manufacturer
  tag.manufacturer <- profile.data$meta$tag.manufacturer
  
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
  
  kdz <- profile.data$attenuation$kdz # where k has units of m^-1
  depth <- profile.data$attenuation$depth # depths at which K were calculated
  predicted_llight <- profile.data$fit_residuals$predicted_llight # predicted light at depth limits of bin bounding kdz estimates
  
  ## Calculate Metrics 
  # Mean K
  Kmean <- mean(kdz)
  
  # Standard deviation of K
  Ksd <- sd(kdz)
  
  # Kurtosis of K
  Kkurt <- moments::kurtosis(kdz)
  
  # Skewness of K
  Kskew <- moments::skewness(kdz)
  
  # Vertically Integrated K (Estrada et al. 1993 Mar Ecol Prog Ser 92: 289-300)
  Ktotal <- sum((kdz[-1] + kdz[-length(kdz)]) * (depth[-1] - depth[-length(depth)]) / 2)
  
  # Light level lost: proportion of light (absolute) lost [i.e., 1 - LLlost = LLratio]
  LLlost <- 1-10^((predicted_llight[length(predicted_llight)] - predicted_llight[1])/D)
  
  # Light level ratio (LLdeepest/LLshallowest): proportion of light (absolute) remaining
  LLratio <- 1-LLlost  
  
  # Depth of max K
  Zkmax <- depth[which.max(kdz)]
  
  # Depth centroid of K (Estrada et al. 1993 Mar Ecol Prog Ser 92: 289-300)
  Zck <- sum((kdz[-1] + kdz[-length(kdz)]) * (depth[-1] - depth[-length(depth)]) * (depth[-1] + depth[-length(depth)]) / (4 * Ktotal))
  
  # Depth (geometric mean) of 80th+ percentile of K values (i.e., top 20%)
  Zgmk80 <- geometric.mean(depth[which(kdz>=quantile(kdz,0.8))])
  
  # Depth (geometric mean) of 20th- percentile of K values (i.e., bottom 20%)
  Zgmk20 <- geometric.mean(depth[which(kdz<=quantile(kdz,0.2))])
  
  # Vertical Stratification Index of K (Napp 1987 Oceanologica Acta 10: 329-337)
  Kbar <- (1 / (2 * (max(depth)-min(depth)))) * sum((kdz[-1] + kdz[-length(kdz)]) * (depth[-1] + depth[-length(depth)]))
  VSIk <- (100 / (max(depth)-min(depth))) * (sum(abs(kdz[-1] - kdz[-length(kdz)])) / Kbar)
  
  # Output
  Kmetrics <- data.frame(Kmean = Kmean,
                     Ksd = Ksd,
                     Kkurt = Kkurt,
                     Kskew = Kskew,
                     LLlost = LLlost,
                     LLratio = LLratio,
                     Ktotal = Ktotal,
                     Zkmax = Zkmax,
                     Zck = Zck,
                     Zgmk80 = Zgmk80,
                     Zgmk20 = Zgmk20,
                     VSIk = VSIk)
  return(Kmetrics)
}
  
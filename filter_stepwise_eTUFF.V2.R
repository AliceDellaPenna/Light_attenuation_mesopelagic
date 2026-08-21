#' @editor Dr. Martin C. Arostegui \email{martin.arostegui@@whoi.edu}
#' filter_stepwise_eTUFF is a modified version of filter_stepwise from the 'trawllight' package.
#' The only edits were to accommodate the different formatting needs of eTUFF data:
#' This included changing column name "trans_llight" to "raw_llight" and eliminating NOAA-specific formatting in relation to their Bering Sea surveys.
#' Edit 1: @param surface.gap *New param: makes allowable gap from surface to profile start (surface.gap) to be different from allowable gap within profile (bin.gap)
#'      If this is not specified, profile will only be filtered by gaps within profiles (no surface gap filtering)
#' Edit 2: @param max.cutoff *New param: only attempt profile generation down to this depth)
#'      
#' Please see GitHub page below and *cite Sean Rohan's work* if using my edit of his function:
#' https://github.com/sean-rohan-NOAA/trawllight
#' Rohan et al. 2020 NOAA Technical Memorandum NMFS-AFSC-403
#' Rohan et al. 2021 Progress in Oceanography 194:102554
#' 
#' Bin and filter light measurements
#'
#' \code{filter_stepwise} aggregates light measurements from a single cast into depth bins using a specificed function (e.g. 'median'), removes light measurements using a stepwise algorithm, and assigns data continuity grades based on a threshold criteria.
#'
#' @param cast.data Data frame containing light measurements and depth.
#' @param bin.size The size of the depth bin used for aggregation. Default = 2.
#' @param bin.gap The maximum size of data gap before a profile is considered to not meet continuity standards. Units are in units of depth, not the number of bins. Default = 6.
#'                ***FORMERLY bin.gap defined the allowable gap WITHIN a light profile (among depth bin midpoints) as well as FROM THE SURFACE (depth bin midpoint of shallowest bin)
#'                Now changed to use bin.gap for WITHIN and surface.gap FROM THE SURFACE
#' @param agg.fun Function used to aggregate light measurements for each depth bin.
#' @param ... Additional arguments passed to findInterval function for binning light measurements
#' @return Data frame with light by depth bin and continuity grade (1 = Good, -999 = Bad)
#' @author Sean Rohan \email{sean.rohan@@noaa.gov}
#' @export

filter_stepwise_eTUFF.V2 <- function(cast.data,
                                  light.col,
                                  depth.col,
                                  bin.size = 2,
                                  bin.gap = 6,
                                  surface.gap = NULL,
                                  agg.fun,
                                  max.cutoff = NULL,
                                  light.floor = NULL,
                                  filter = T, ...) {
  
  names(cast.data)[which(names(cast.data) == light.col)] <- "raw_llight"
  names(cast.data)[which(names(cast.data) == depth.col)] <- "cdepth"
  
  max.depth <- max(ceiling(cast.data$cdepth), na.rm = T)
  
  # Bin by depth with bins centered
  cast.data$cdepth <- findInterval(x = cast.data$cdepth, vec = seq(0, max.depth, bin.size), rightmost.closed = T, left.open = F) * bin.size - bin.size/2
  
  # Calculate binned light level using user-specified function
  light_at_depth <- aggregate(raw_llight~cdepth,
                              data = cast.data,
                              FUN = agg.fun)
  
  
  light_at_depth <- light_at_depth[order(light_at_depth$cdepth),]
  
  
  # Stepwise measurement removal loop
  if(filter) {
    p2 <- 1
    while(p2 < nrow(light_at_depth) ) {
      if(nrow(light_at_depth) >= (p2 + 1)) {
        if((light_at_depth$raw_llight[p2 + 1] > light_at_depth$raw_llight[p2])) {
          light_at_depth <- light_at_depth[-p2,]
          p2 <- 0 # Index back to start
        }
      }
      p2 <- p2 + 1
    }
    
    
    
    if (is.null(max.cutoff) == F) {
      light_at_depth <- subset(light_at_depth, cdepth <= max.cutoff)
      if (nrow(light_at_depth) == 0) {
        light_at_depth <- as.data.frame(matrix(nrow=2,ncol=3));colnames(light_at_depth) <- c("cdepth","raw_llight","quality")
        light_at_depth$cdepth[1:2] <- c(10,1000); light_at_depth$raw_llight[1:2] <- c(200,25); light_at_depth$quality[1:2] <- -999
      } 
    }
    
    ###
    if (is.null(light.floor) == F) {
      light_at_depth <- subset(light_at_depth, raw_llight > (light.floor))
      if (nrow(light_at_depth) == 0) {
        light_at_depth <- as.data.frame(matrix(nrow=2,ncol=3));colnames(light_at_depth) <- c("cdepth","raw_llight","quality")
        light_at_depth$cdepth[1:2] <- c(10,1000); light_at_depth$raw_llight[1:2] <- c(200,25); light_at_depth$quality[1:2] <- -999
      } 
    }
    ###
    
      # Assign data continuity codes. -999 indicates gap >= bin.gap and/or (optionally) surface.gap
      if (is.null(surface.gap)){
        if(max(diff(light_at_depth$cdepth)) <= bin.gap) {
          light_at_depth$quality <- 1
        } else {
          light_at_depth$quality <- -999
        }
      } else {
        if(max(diff(light_at_depth$cdepth)) <= bin.gap & min(light_at_depth$cdepth + 1) <= surface.gap) {
          light_at_depth$quality <- 1
        } else {
          light_at_depth$quality <- -999
        }
      }
    }
  
  names(light_at_depth)[which(names(light_at_depth) == "raw_llight")] <- light.col
  names(light_at_depth)[which(names(light_at_depth) == "cdepth")] <- depth.col
  
  return(light_at_depth)
}

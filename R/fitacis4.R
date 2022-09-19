#' Fit A-Ci curves with custom kinetics and temperature response of mesophyll
#' conductance
#'
#' @param data data frame with A/Ci curves. Requires net CO2 assimilation
#' (Anet/Photo/ALEAF in umol m-2 s-1), leaf temperature (Tleaf in Celsius),
#' intercellular CO2 concentration (Ci, in umol mol-1), PPFD incident on the
#' leaf (PPFD/PARi in umol m-2 s-1), atmospheric pressure (Patm/Press in kPa),
#' and (optional, set useRd = TRUE for this option) respiration (Rd, in
#' umol m-2 s-1).
#' @param group1 grouping variable 1, could be species, temperature, ID
#' @param group2 grouping variable 2
#' @param group3 grouping variable 3
#' @param gm25 mesophyll conductance at 25 Celsius in mol m-2 s-1 bar-1
#' @param Egm activation energy of mesophyll conductance in kJ mol-1
#' @param K25 Km in 21% O2 (aka Kcair, aka apparent Km in 21% O2) at 25 Celsius
#' in umol mol-1 (equivalent to ubar bar-1)
#' @param Ek activation energy of Kcair in kJ mol-1
#' @param Gstar25 photorespiratory CO2 compensation point at 25 Celsius in
#' umol mol-1 (equivalent to ubar bar-1)
#' @param Egamma activation energy of GammaStar in kJ mol-1
#' @param fitmethod Set to either "bilinear" or "default". Default option
#' in this package is "default". See ?fitaci in plantecophys for more
#' details.
#' @param fitTPU Should TPU limitations be fit? Set to TRUE/FALSE. See
#' ?fitaci in plantecophys for more details.
#' @param Tcorrect Should outputs be temperature corrected? Default here
#' is FALSE. See ?fitaci in plantecophys for more details.
#' @param useRd Should respiration be used? Default is FALSE. See ?fitaci
#' in plantecophys for more details.
#' @param citransition Pre-specify Ci transition point? Units in umol mol-1
#' (ubar bar-1) Default is FALSE. See ?fitaci in plantecophys for more details.
#' @param alphag Fraction of photorespiratory glycolate carbon that is not
#' returned to the chloroplast (von Caemmerer, 2000). If ACi curves show
#' high-CO2 decline, then this value should be > 0. See ?fitaci in plantecophys
#' for more details.
#' @param PPFD Light intensity? Can be retrieved from dataframe. Default is
#' NULL. Units are umol m-2 s-1. See ?fitaci in plantecophys for more details.
#' @param Tleaf Leaf temperature? Can be retrieved from dataframe. Default is
#' NULL. Units are Celsius. See ?fitaci in plantecophys for more details.
#' @param alpha Quantum yield of electron transport. Default is 0.24. Units are
#' umol electrons / umol incident photons. See  ?fitaci in plantecophys for
#' more details.
#' @param theta Curvature of the photosynthetic light response. Default is
#' 0.85. If light response has sharper transition, increase up to 1. If light
#' response has shallower curves, decrease towards 0. See ?fitaci in
#' plantecophys for more details.
#' @param varnames Variable names in your dataframe. ALEAF is net CO2
#' assimilation in umol m-2 s-1, Tleaf is leaf temperature in Celsius,
#' Ci is intercellular CO2 concentration in umol mol-1, PPFD is light
#' intensity in umol m-2 s-1, Rd is respiration rate in umol m-2 s-1,
#' and Press is atmospheric pressure in kPa. See ?fitaci in plantecophys
#' for more details.
#' @param ... Further arguments for plantecophys::fitaci(). See ?fitaci for
#' details.
#' @return fitacis2 allows gmeso, GammaStar, and Km to vary with Tleaf.
#' Output matches the fitacis function from plantecophys. Note that the
#' temperature response function of Km is derived from the temperature
#' responses of Ko and Kc in Bernacchi et al.2001, as is the GammaStar
#' temperature response defaults. The gm defaults are from Bernacchi et
#' al. 2002 fitted between 1and 35 Celsius. Also note that this ALWAYS
#' uses gm. To fit data on a "Ci-basis", set gm25 really high (e.g.
#' 10000 mol m-2 s-1 bar-1) and Egm to 0 kJ mol-1.
#' 
#' In some instances (e.g. very low stomatal conductance), fitacis2 will fail.
#' In these cases, the output for that curve will be "Failed", rather than an
#' A-Ci curve fit object.
#' 
#' REFERENCES
#'
#' Bernacchi CJ, Singsaas EL, Pimentel C, Portis AR, Long SP.
#' 2001. Improved temperature response functions for models of
#' rubisco-limited photosynthesis. Plant Cell Environment 24:253-259.
#'
#' Bernacchi CJ, Portis AR, Nakano H, von Caemmerer S, Long SP. 2002.
#' Temperature response of mesophyll conductance. Implications for the
#' determination of rubisco enzyme kinetics and for limitations to
#' photosynthesis in vivo. Plant Physiology 130:1992-1998.
#' von Caemmerer S. 2000. Biochemical models of leaf photosynthesis.
#' CSIRO Publishing, Collingwood.
#'
#' @importFrom tidyr unite
#' @importFrom plantecophys fitaci
#' @export
#' @examples \donttest{
#' #Read in data
#' data <- read.csv(system.file("extdata", "example_2.csv",
#' package = "plantecowrap"), stringsAsFactors = FALSE)
#' #Run ACi curve fitting
#' fits <- fitacis2(data, group1 = "Grouping",
#' varnames = list(ALEAF = "A",
#'                 Tleaf = "Tleaf",
#'                 Ci = "Ci",
#'                 PPFD = "PPFD",
#'                 Rd = "Rd",
#'                 Press = "Press"),
#' fitmethod = "bilinear", fitTPU = TRUE, Tcorrect = FALSE)
#' }
#'


##==============================================================================
## FUNCTION fitacis4
##
## Three options for mesophyll conductance temperature response
##
## [1] gm_method = "Quad"  : quadratic equation fitted to actual gm values
##
##     Input is the quadratic equation fitted to measured (actual) gm values.
##     The equation is entered with the three coefficients a, b, c in the
##     equation of form y = a + bx + cx^2
##
##     There is no gm25 input needed/used
##
## [2] gm_method = "Arr1"  :  Arrhenius equation
##
##     Input is gm at 25 C and 'activation energy'. The underlying temperature
##     response, characterized by the activation energy, is scaled to 1 at 25 C;
##     the temperature respose is then multiplied by gm25
##
## [3] gm_method = "Arr2"  :  'double' Arrhenius equation
##
##     Input is gm at 25 C. The underlying temperature response is characterized
##     by the parameters from Bernacchi et al. (2002) Plant Physiology 130:
##     1992-1998, which is scaled to 1 at 25 C. That temperature response is
##     then multiplied by gm25 (analogous to Arr1)
##
##     PLEASE NOTE: Rgas is truncated to four digits (0.008314) to match the
##     Bernacchi equation, BUT 'c' is written as 20.0098424
##     For Bernacchi et al. (2002), with 'c' = 20.0, the relative rate
##     at 25 C is 0.9902058 (NOT 1.0). For 'c' = 20.01 (as in Sharkey et al.
##     [2007] PCE 30: 1035-1040), the relative rate at 25 C is 1.0001576. For
##     'c' = 20.0098424 (herein), the relative rate at 25 C is 1.0000000 (or
##     with one more digit: 0.99999996)
##
## DEFAULT values for each method are in the native function, and for Arr2, if
## the intent is use Bernacchi et al.'s equation, leave them alone and only
## enter a new gm25
##
## For method [1] the quadratic equation should already be in the correct units
## For methods [2] and [3] the gm25 entered should be given with correct units;
## the temperature response is independent of units (e.g., bar v. Pa)
##
## Arguments for mesophyll conductance calculations using quadratic temperature
## response:
##
##   gm_coef0    : intercept in gmeso v. T relationship
##   gm_coef1    : coefficient on linear term in gmeso v. T relationship
##   gm_coef2    : coefficient on squared term in gmeso v. T relationship
##------------------------------------------------------------------------------

  fitacis4 = function (data,
                     group1 = "Plant",
                     group2 = NA,
                     group3 = NA,
#-------------------------------------------------------------------------------
# temperature response for gm from one of three functional forms

                     gm_method = "Arr2",      # can be "Quad", "Arr1", "Arr2"

                     gm_coef0  = -0.3403706,  # from Demi's cotton data
                     gm_coef1  =  0.0522454,  # from Demi's cotton data
                     gm_coef2  = -0.0006509,  # from Demi's cotton data

                     Egm       = 47.650,   # kJ mol-1          [existing value]
                     gm25      = 0.08701,  # mol m-2 s-1 bar-1 [existing value]

                     Arr2_c    = 20.0098424,  # Bernacchi ['corrected' by Jeff]
                     Arr2_Ha   = 49.6,        # Bernacchi
                     Arr2_Hd   = 437.4,       # Bernacchi
                     Arr2_S    = 1.4,         # Bernacchi

# end of temperature response for gm
#-------------------------------------------------------------------------------
                     K25          = 718.40,    #umol mol-1 (ubar bar-1)
                     Ek           = 65.50828,  #kJ mol-1
                     Gstar25      = 42.75,     #umol mol-1 (ubar bar-1)
                     Egamma       = 37.83,     #kJ mol-1
                     fitmethod    = "default",
                     fitTPU       = TRUE,
                     Tcorrect     = FALSE,
                     useRd        = FALSE,
                     citransition = NULL,   #umol mol-1
                     alphag       = 0,
                     PPFD         = NULL,   #umol m-2 s-1
                     Tleaf        = NULL,   #Celsius
                     alpha        = 0.24,   #umol CO2 / umol photons
                     theta        = 0.85,
                     varnames     = list(ALEAF = "Photo",  #umol m-2 s-1
                                         Tleaf = "Tleaf",  #Celsius
                                         Ci = "Ci",        #umol mol-1
                                         PPFD = "PARi",    #umol m-2 s-1
                                         Rd = "Rd",        #umol m-2 s-1
                                         Press = "Press"), #kPa
                     ...) {

# Assign group names and pressure

  data$group1 <- data[, group1]
  data$Press  <- data[, varnames$Press]

  if (!is.na(group2)) { data$group2 <- data[, group2] }
  if (!is.na(group3)) { data$group3 <- data[, group3] }

  if (!is.na(group2) & !is.na(group3)) {
  data <- unite(data, col = "group", c("group1", "group2", "group3"), sep = "_")
  } else {
    if (!is.na(group2) & is.na(group3)) {
      data <- unite(data, col = "group", c("group1", "group2"), sep = "_")
    } else {
      data$group <- data$group1
    }
  }

  data <- split(data, data$group)           # Split data by group
  fits <- as.list(1:length(data))           # Create empty list for curve fits

# Fit an A-Ci curve for each 'dataset'/experiment in 'data'

  for (i in 1:length(data)) {  

#-------------------------------------------------------------------------------
# Re-written and new formulations for Km, GammaStar, and gmeso [Jeff]

    Patm   = mean(data[[i]]$Press)          # measured air pressure
    Rgas   = 0.008314                       # gas constant (kJ / K.mole)
                                            # (truncated at 4 digits for
                                            # compatibility with Bernacchi)
    Abs_0C = -273.15                        # absolute zero in degree Celsius
    T25K   = 25 - Abs_0C                    # absolute temperature @ 25 deg C
    Tleaf  = mean(data[[i]]$Tleaf)          # measured Tleaf (degree Celsius)
    TleafK = Tleaf - Abs_0C                 # measured Tleaf (K)
    RgasTK = Rgas * TleafK                  # R.Tk

# Tfactor = 0 at 25 degrees Celsius (and e^0 = 1) for the Arrhenius temperature
# response. The 'activation energies' and values at 25 deg C are then used
# with Tfactor for Km, GammaStar, and gmeso temperaure responses

    Tfactor = (TleafK - T25K) / (T25K * RgasTK) # TleafK factor relative to 25 C

# Calculate Km and GammaStar from Tfactor (and gmeso below for gm_method 'Arr1')

    Km        = K25     * exp(Ek     * Tfactor)
    GammaStar = Gstar25 * exp(Egamma * Tfactor)

# Calculate gmeso using one of the three temperature response methods
#
#   gm_method = "Quad"           :: use Tleaf (deg C) (fast quadratic)
#   gm_method = "Arr1" or "Arr2" :: use TleafK Arrhenius approaches

    if        (gm_method == 'Quad') {        ## quadratic fitted to data

      gmeso = gm_coef0 + (gm_coef1 + gm_coef2 * Tleaf) * Tleaf

    } else if (gm_method == 'Arr1') {        ## Arrhenius temperature response

      gmeso = gm25 * exp(Egm * Tfactor)

    } else if (gm_method == 'Arr2') {        ## Bernacchi temperature response

      gmeso = gm25 * exp(Arr2_c - Arr2_Ha / RgasTK) /
                     (1 + exp((Arr2_S * TleafK - Arr2_Hd) / RgasTK))

    } else                          {        ## BOGUS gm_method value: ERROR

      print ('ERROR--Invalid gm_method')
      print ('gm_method entered is')
      print (gm_method)
    }

#-------------------------------------------------------------------------------
# Fit ACi curve by calling the 'fitaci' function with desired parameter values

    fits[[i]] <- tryCatch(fitaci(data[[i]], Patm = Patm, varnames = varnames,
                   fitmethod = fitmethod, Tcorrect = Tcorrect, fitTPU = fitTPU,
                   gmeso = gmeso, Km = Km, GammaStar = GammaStar, useRd = useRd,
                   citransition = citransition, alphag = alphag, PPFD = PPFD,
                   Tleaf = Tleaf, alpha = alpha, theta = theta,  ...),
                   error = function(e) paste("Failed") )

    names(fits)[i] <- data[[i]]$group[1]    # Assign names
  }                                         # End of this curve fit

  return(fits)                              # Return all curve fits in list
}                                           # End of function

##==============================================================================
#
# A more precise (and accurate) calculation of saturation vapor pressurve (SVP)
# than generally used in plant physiology modeling when vapor pressure deficit
# is of concern -- ASSUMING ATMOSPHERIC PRESSURE OF 100 kPa which will usually
# for most locations be incorrect.
#
# Saturation Vapor Pressure Calculation
#
# Saturation water vapor pressure (SVP) of moist air with respect to a plane
# surface of water is calculated according to equations (17), (21), and (22) in
# Alduchov & Eskridge (1996), who used pressure in hPa. Instead, we are using
# Pa here. The resulting equation for SVP is:
#
# SVP = 1.00071 × EXP(4.5E-8 P) × 610.94 × EXP(17.625 T / (243.04 + T))
#
# where P is atmospheric pressure (Pa) and T is air temperature (°C). We are
# assuming P = 100,000 Pa, so 1.00071 × EXP(4.5E-8 P) × 610.94 reduces to about
# 614.131 , giving
#
# SVP = 614.131 × EXP(17.625 T / (243.04 + T)).
#
# The assumption of P = 100 kPa would, on average, apply to an altitude of about
# 113 m (at 20–21°C).
#
# REFERENCE:  Alduchov OA, Eskridge RE (1996) Improved Magnus form approximation
# of saturation vapor pressure. Journal of Applied Meteorology 35: 601-609.
#
#
# SVP   : saturation vapor pressure, Pa
# T     : air temperature, degrees Celsius

  preciseSVP = function(T) { SVP = 614.131 * exp(17.625 * T / (243.04 + T)) }

##==============================================================================

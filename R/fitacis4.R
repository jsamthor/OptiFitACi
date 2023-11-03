#==============================================================================
# FUNCTION fitacis4
#
# The function 'fitacis4' reads a C3 leaf photosynthetic gas exchange dataset
# composed of CO2 exchange rates at different intercellular CO2 concentrations
# (A-Ci data) and prepares the data for input to the function 'fitaci' in the R
# package 'plantecophys,' authored by Remko Duursma, that fits functions to the
# A-Ci data to estimate photosynthetic parameters including V<c,max> and J<max>,
# the maximum rate of carboxylation (umol m-2 s-1) and light-saturated potential
# rate of photosynthetic electron transport (umol m-2 s-1), respectively [see
# https://cran.r-project.org/web/packages/plantecophys/index.html].
#
# 'fitacis4' is an extension of the function 'fitacis2' in the R package
# 'plantecowrap,' authored by Joseph Stinziano, Demi Gamble [Demi Sargent],
# Robert Sharwood, and Warren Conaty, that provides extended flexibility in the
# use of inputs to the function 'fitaci' in the R package 'plantecophys' [see
# https://cran.r-project.org/web/packages/plantecowrap/index.html
#
# One purpose of 'fitacis4' is to allow input of a value for leaf mesophyll
# conductance, g<m> (mol m-2 s-1 bar-1), and the response of g<m> to
# temperature should the input dataset contain measurements at multiple
# temperatures. The 'fitaci' default g<m> value is infinite, meaning that there
# would be no resistance to CO2 transport between intercellular spaces and
# the chloroplast stroma where CO2 assimilation biochemistry occurs in C3 
# photosynthesis. That default infinite value for g<m> can be replaced by a
# user-supplied value (variable 'gmeso') that is constant across temperatures.
# The function 'fitacis2' always uses a finite value of g<m> at 25 deg Celsius,
# with a default value from published data collected from tobacco leaves, though
# it can be set to any value, and includes a temperature response function of
# the form g<m> = g<m,25> * exp(TleafK - 298.15) / (298.15 * TleafK * 0.008314)
# where g<m,25> is the input value of g<m> at 25 deg Celsius, TleafK is leaf
# temperature (kelvin), and 0.008314 approximates the gas constant R (MJ mol-1
# K-1). The function 'fitacis2' also allows for temperature responses of The
# input variables Km, the Michaelis-Menten constant for carboxylation by
# rubisco (umol mol-1 or ppm), and gammastar, the photorespiratory CO2
# compensation point (umol mol-1).
#
# The function 'fitacis4' extends the treatment of the temperature response of
# g<m>, allowing for either the 'fitaci' default value of infinite g<m> at all
# temperatures, the function used in fitacis2 (but with a faster calculation),
# and two new functions, labeled methods "Quad" and "Arr2" and described below.
#
# The three (four) options for mesophyll conductance and its response to
# temperature in 'fitacis4' are [three g<m> treatments and a fourth option to
# ignore g<m>]:
#
# [1] gm_method = "Quad"  : quadratic equation fitted to actual gm values
#
#     Input is the quadratic equation fitted to measured (actual) gm values.
#     The equation is entered with the three coefficients a, b, c in the
#     equation of form g<m> = a + b.Tleaf + c.Tleaf^2 where Tleaf is leaf
#     temperature in deg Celsius
#
#     The variable names for coefficients a, b, and c used for mesophyll
#     conductance calculations using a quadratic temperature response are:
#
#     gm_coef0    : intercept in gmeso v. T relationship
#     gm_coef1    : coefficient on linear term in gmeso v. T relationship
#     gm_coef2    : coefficient on squared term in gmeso v. T relationship
#
#     There is no gm25 input needed/used for method = "Quad"
#
# [2] gm_method = "Arr1"  :  Arrhenius equation
#
#     Input is gm at 25 C (gm25) and 'activation energy' (Egm). The underlying
#     temperature response, characterized by the activation energy, is scaled to
#     1 at 25 C; the temperature respose is then multiplied by gm25. This is the
#     formulation in 'fitacis2'
#
# [3] gm_method = "Arr2"  :  Arrhenius equation with 'deactivation'
#
#     Input is gm at 25 C (gm25), an activation energy (Arr2_Ha), an entropy
#     term (Arr2_S), and a deactivation term (Arr2_Hd). The deault temperature
#     response is parameters are from Bernacchi et al. (2002) Plant Physiology 130:
#     1992-1998, which is scaled to 1 at 25 C. That temperature response is
#     then multiplied by gm25 (analogous to Arr1)
#
# [4] gm_method = "ignore"  :  infinite mesophyll conductance
#
#     This sets gmeso = NULL, which is the 'fitaci' default in the R package
#     'plantecophys'
#
# DEFAULT values for each method are in the native function, and for Arr2, if
# the intent is use Bernacchi et al.'s equation, leave them alone and only
# enter a new gm25
#
# For method [1] the quadratic equation should already be in the correct units
# For methods [2] and [3] the gm25 entered should be given with correct units
# (i.e., mol m-2 s-1 bar-1); the temperature response is independent of units
# (i.e., it is a multiplier/scaler relative to gm25)
#
# Calculations of the temperature responses of Km and gammastar in fitacis4
# are the same as in fitacis2 (but slightly faster)
#
# J.S. Amthor, November 2023 [update]
#------------------------------------------------------------------------------

  fitacis4 = function (data,
               group1 = "Plant",
               group2 = NA,
               group3 = NA,
#-------------------------------------------------------------------------------
# temperature response for gm from one of three functional forms

               gm_method = "Quad",      # can be "Quad", "Arr1", "Arr2", "ignore"

               gm_coef0  =  0.571,      # cotton data
               gm_coef1  = -0.009,      # cotton data
               gm_coef2  =  0.000355,   # cotton data

               gm25      = 0.53,        # mol m-2 s-1 bar-1 [cotton data]

               Egm       = 8.59,        # kJ mol-1 [cotton data for Arr1 method]
                       
               Arr2_c    = 20.0098,     # based on Bernacchi et al. (2002)
               Arr2_Ha   = 49.6,        # based on Bernacchi et al. (2002)
               Arr2_Hd   = 437.4,       # based on Bernacchi et al. (2002)
               Arr2_S    = 1.4,         # based on Bernacchi et al. (2002)

# end of temperature response for gm
#-------------------------------------------------------------------------------
               K25          = 541.47,    # umol mol-1 @ 25 deg C [cotton]
               Ek           = 53.03,     # kJ mol-1              [cotton]
               Gstar25      = 39.98,     # umol mol-1 @ 25 C     [cotton]
               Egamma       = 28.11,     # kJ mol-1              [cotton]
               fitmethod    = "default",
               fitTPU       = TRUE,
               alphag       = 0,         # default, fraction (0-1) of glycolate NOT
                                         # returned to the chloroplast in
                                         # photorespiration. alphag = 0 is the MOST
                                         # efficient option (flat line TPU limitation)
               Tcorrect     = FALSE,
               useRd        = FALSE,
               citransition = NULL,      # umol CO2 mol-1
               PPFD         = NULL,      # umol photons (PAR) m-2 s-1
               Tleaf        = NULL,      # degress Celsius
               alpha        = 0.24,      # umol CO2 / umol photons
               theta        = 0.85,      # [plantecophys/fitaci default]

               varnames     = list(ALEAF = "Photo",  # umol CO2 m-2 s-1
                                   Tleaf = "Tleaf",  # degrees Celsius
                                   Ci = "Ci",        # umol CO2 mol-1
                                   PPFD = "PARi",    # umol photon m-2 s-1
                                   Rd = "Rd",        # umol CO2 m-2 s-1
                                   Press = "Press"), # kPa
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
# Re-written and new formulations for Km, GammaStar, and gmeso

    Patm   = mean(data[[i]]$Press)          # measured air pressure
    Rgas   = 0.00831446                     # gas constant (kJ / K.mole)
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

    if        (gm_method == 'Quad') {        # quadratic fitted to data

      gmeso = gm_coef0 + (gm_coef1 + gm_coef2 * Tleaf) * Tleaf

    } else if (gm_method == 'Arr1') {        # Arrhenius temperature response

      gmeso = gm25 * exp(Egm * Tfactor)

    } else if (gm_method == 'Arr2') {        # Bernacchi temperature response

      gmeso = gm25 * exp(Arr2_c - Arr2_Ha / RgasTK) /
                     (1 + exp((Arr2_S * TleafK - Arr2_Hd) / RgasTK))
					 
	} else if (gm_method == 'ignore') {      # gmeso = NULL ; plantecophys default	

      gmeso = NULL	

    } else                          {        # BOGUS gm_method value: ERROR

      print ('ERROR--Invalid gm_method')
      print ('gm_method entered is')
      print (gm_method)
    }

#-------------------------------------------------------------------------------
# If we are using experimentally measured/estimated Rd (rather than letting the
# fitting algorithm estimate Rd) -- i.e., using useRd = TRUE -- then ALL the Rd
# values for an individual set of A-Ci measurements MUST BE IDENTICAL (e.g., all
# the data lines for a rep of a species/treatment at a single temperature must
# be the same). So if there are different Rd values for different measurements
# contributing to a single A-Ci fit, we use their average (we also use their
# average if they are all the same)

    if (useRd) data[[i]]$Rd = mean(data[[i]]$Rd)

#-------------------------------------------------------------------------------
# Fit A-Ci curve by calling the 'fitaci' function with desired input values

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

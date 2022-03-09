#' @title Example Data Frame for Quantity-intensity relationship
#' @description User is advised to prepare the data as suggested in the example to derive Quantity-Intensity (Q/I) relationship parameters of soil potassium (K) using linear and polynomial (second order) regression
#' @usage df
#' @format Data frame with Solution: soil ratio in the first column, Initial K concentration (mg/L) in the second column, Final or equilibrium K concentration (mg/L) in the third column, Final or equilibrium 'Ca+Mg' concentration (mol/L) in the fourth column. Write the following notations on the spreadsheet:
#'
#' Solution_to_Soil_Ratio - for Solution: soil ratio
#'
#' Initial_K - for Initial K concentration (mg/L)
#'
#' Final_K - for Final or equilibrium K concentration (mg/L)
#'
#' Final_Ca_and_Mg - for Final or equilibrium 'Ca+Mg' concentration (mol/L)
#'
#' @export
df = read.table(text = "Solution_to_Soil_Ratio	Initial_K	Final_K	Final_Ca_and_Mg
40	0	0.795333333	0.01
20	0	1.154333333	0.01
10	0	1.73	0.01
10	5	4.256666667	0.01
10	10	6.4	0.01
10	20	10.13666667	0.01
10	40	17.71333333	0.01
10	60	29.11333333	0.01
10	80	39.51666667	0.01
10	120	66.23333333	0.01", header = TRUE)

#' @title Quantity-intensity relationship derived through linear regression
#'
#' @description The quantity-intensity (Q/I) relationships of soil K, introduced by Beckett (1964), is implemented in this function using linear regression equation as used by some earlier workers (Zhang et al., 2011; Islam et al., 2017; Das et al., 2019; 2021).
#'
#' @importFrom  ggplot2 ggplot aes element_blank element_text
#' @importFrom stats lm predict
#'
#' @usage QIlin(Solution2Soil = Solution2Soil, CKi = CKi, CKf = CKf, CCaMg = CCaMg,
#' NH4OAC_K = NH4OAC_K)
#'
#' @param  Solution2Soil Ratio of solution volume to soil mass (mL/g or L/kg)
#' @param  CKi Initial K concentration (mg/L)
#' @param  CKf Final or equilibrium K concentration (mg/L)
#' @param  CCaMg Final or equilibrium 'Ca+Mg' concentration (mol/L)
#' @param  NH4OAC_K K extracted from soil by 1 N ammonium acetate (NH4OAc) of pH 7 (mg/kg)
#' @return AReK - Equilibrium activity ratio (unitless)
#' -deltaK0 - Non-specifically held K (cmolc/kg)
#' Ks - Specifically held K (cmolc/kg)
#' PBCK - Potential buffering capacity (cmolc/kg)
#' deltaG0 - The standard free energy of exchange (cal/mol)
#'
#'
#' @references Beckett, P.H.T., 1964. The immediate Q/I relations of labile potassium in the soil. European Journal of Soil Science 19, 9-23.
#'
#' Bilias, F., Barbayiannis, N., 2018. Contribution of non-exchangeable potassium on its quantity-intensity relationships under K-depleted soils. Archives of Agronomy and Soil Science 64, 1988-2004.
#'
#' Das, D., Dwivedi, B.S., Datta, S.P., Datta, S.C., Meena, M.C., Agarwal, B.K., Shahi, D.K., Singh, M., Chakraborty, D., Jaggi S., 2019. Potassium supplying capacity of a red soil from eastern India after forty-two years of continuous cropping and fertilization. Geoderma 341: 76-92.
#'
#' Das D, Dwivedi BS, Datta SP, Datta SC, Meena MC, Dwivedi AK, Singh M, Chanraborty D, Jaggi S (2021) Long-term differences in nutrient management under intensive cultivation alter potassium supplying ability of soils. Geoderma 393:114983.
#'
#' Das, D., Nayak, A.K., Thilagam, V.K., Chatterjee, D., Shahid, M., Tripathi, R., Mohanty, S., Kumar, A., Lal, B., Gautam, P., Panda, B.B., Biswas, S.S., 2018. Measuring potassium fractions is not sufficient to assess the long-term impact of fertilization and manuring on soil's potassium supplying capacity. Journal of Soils and Sediments 18, 1806-1820.
#'
#' Evangelou, V.P., Blevins, R.L., 1988. Effect of long-term tillage systems and nitrogen addition on potassium quantity-intensity relationships. Soil Sci. Soc. Am. J. 52, 1047-1054.
#'
#' Islam, A., Karim, A.J.M.S., Solaiman, A.R.M., Islam, M.S., Saleque, M.A., 2017. Eight-year long potassium fertilization effects on quantity/intensity relationship of soil potassium under double rice cropping. Soil Till. Res. 169, 99-117.
#'
#' Le Roux, J., Summer, M.E., 1968. Labile potassium in soils, I: Factors affecting the quantity-intensity (Q/I) parameters. Soil Science 106, 35-41.
#'
#' Sparks, D.L., Liebhardt, W.C., 1981. Effect of long-term lime and potassium application on quantity-intensity (Q/I) relationships in sandy soil. Soil Science Society of America Journal 45,786-790.
#'
#' Zhang, H., Xu, M., Zhu, P., Peng, C., 2011. Effect of 15-year-long fertilization on potassium quantity/intensity relationships in black soil in Northeastern China. Communications in Soil Science and Plant Analysis 42, 1289-1297.
#'
#' @details A number of parameters related to soil K availability can be obtained from the Q/I plot, e.g., equilibrium activity ratio (AReK), total labile K (KL), non-specifically held K (-deltaK0), specifically held K (Ks), potential buffering capacity (PBCK), and standard free energy of exchange (deltaG0). The equilibrium activity ratio (AReK) is defined as the activity ratio of K to Ca or 'Ca+Mg' when there is no net adsorption or desorption of K between soil solution and exchange phases. It is a measure of the intensity factor. Total labile K is the amount of K held on the soil solids which is capable of ion exchange reactions during the time period provided for equilibration between soil solution and soil solids. It is a measure of the quantity factor. Conventionally, the total labile K has been sub-divided into non-specifically held K, which is mainly bound to the planar sites; and specifically held K, which is mainly bound to the edge/wedge positions of 2:1 clay minerals (Sparks and Liebhardt, 1981). The potential buffering capacity (PBCK) is a measure of the ability of a soil to resist the changes in intensity factor after additions or losses of K from the system.
#' @export
#' @examples
#' with(data = df, QIlin(Solution2Soil = Solution_to_Soil_Ratio, CKi = Initial_K,
#' CKf = Final_K, CCaMg = Final_Ca_and_Mg, NH4OAC_K = 55))

QIlin <- function(Solution2Soil = Solution2Soil, CKi = CKi, CKf = CKf, CCaMg = CCaMg,
               NH4OAC_K = NH4OAC_K){
  #delta K calculation in cmol/kg
  deltak <- (CKi - CKf)*Solution2Soil/390

  #Conversion of final K concentration from mg/L to mol/L
  Ck_f_mol_L <- CKf/39000

  ARk <- (Ck_f_mol_L/CCaMg^0.5)

  #Concentration of Cl (mol/L)
  C_Cl <- Ck_f_mol_L + 2*CCaMg

  #Ionic strength calculation
  I = 0.5*(Ck_f_mol_L*1^2 + CCaMg*2^2 + C_Cl*1^2)

  #Activity coefficient calculation
  activity_coef <- exp((0.5*I^0.5)/(1 + I^0.5) - 0.086*C_Cl)

  #ARk calculation
  ARk <- (Ck_f_mol_L/CCaMg^0.5)*(activity_coef)
  dat <- cbind.data.frame(deltak, ARk)
  #Plotting using ggplot2
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    myplot <- function(mydf, xcol, ycol){
      ggplot2::ggplot(data = mydf, aes(x = {{xcol}}, y = {{ycol}})) +
        ggplot2::geom_point(color='red', alpha=0.3, size=2)+
        ggplot2::stat_smooth(method='lm', se = F, formula = y ~ x) +
        ggplot2::theme_bw() +
        ggplot2::xlab(expression(paste("AR"["K"]))) +
        ggplot2::ylab(expression(paste(Delta,"K", " cmol"["c"], " kg"^-1)))+
        ggplot2::theme(panel.grid = element_blank()) +
        ggplot2::theme(text=element_text(family = "serif", size=15),
                       axis.title.x = element_text(size = 18),
                       axis.title.y = element_text(size = 18),
                       axis.text.x = element_text(colour="black",face="bold"),
                       axis.text.y = element_text(colour="black",face="bold"))
    }
  } else {
    ## Please install the "ggplot2" package
  }

  my_plot <- myplot(dat, ARk, deltak)

  mod <- lm(deltak ~ ARk)
  PBC <- mod$coefficients[[2]]
  deltaK0 <- abs(mod$coefficients[[1]])
  AReK <- deltaK0/PBC

  #mg/kg to	cmolc/kg conversion
  NH4OAC_K_cmolc_kg <- NH4OAC_K/390

  Ks <- NH4OAC_K_cmolc_kg - deltaK0
  #deltaG0 = RTln(AReK)
  #T : 298 Kelvin
  #R: Universal gas constant
  deltaG0 <- 1.9858775*298*(log(AReK))

  return(list(`PBC` = PBC, `DeltaK0` = deltaK0,
           `AReK` = AReK, `Ks` = Ks, `deltaG0` = deltaG0,
           `Plot` = my_plot))
}

#' @title Quantity-intensity (Q/I) relationship of soil K derived through a second order polynomial i.e., quadratic equation
#'
#' @description A quadratic equation of the form "y = ax2 + bx + c" can be fitted to Q/I data to find out different Q/I parameters
#' @importFrom  ggplot2 ggplot aes element_blank element_text
#' @importFrom stats lm predict
#'
#' @usage QIPoly(Solution2Soil = Solution2Soil, CKi = CKi, CKf = CKf, CCaMg = CCaMg)
#'
#' @param  Solution2Soil Ratio of solution volume to soil mass (mL/g or L/kg)
#' @param  CKi Initial K concentration (mg/L)
#' @param  CKf Final or equilibrium K concentration (mg/L)
#' @param  CCaMg Final or equilibrium 'Ca+Mg' concentration (mol/L)
#' @return AReK - Equilibrium activity ratio (unitless)
#' Kl - Total labile K (cmolc/kg)
#' PBCK - Potential buffering capacity (cmolc/kg)
#' deltaG0 - The standard free energy of exchange (cal/mol)
#'
#'
#' @references Wang, J.J., Harrell, D.L., Bell, P.F., 2004. Potassium buffering characteristics of three soils low in exchangeable potassium. Soil Science Society of America Journal 68, 654-661.
#'
#' Wang, J.J., Scott, A.D., 2001. Effect of experimental relevance on potassium Q/I relationships and its implications for surface and subsurface soils. Communications in Soil Science and Plant Analysis 32, 2561-2575.
#'
#' @export
#' @examples
#' with(data = df, QIPoly(Solution2Soil = Solution_to_Soil_Ratio, CKi = Initial_K,
#' CKf = Final_K, CCaMg = Final_Ca_and_Mg))

QIPoly <- function(Solution2Soil = Solution2Soil, CKi = CKi, CKf = CKf, CCaMg = CCaMg){
  #delta K calculation in cmol/kg
  deltak <- (CKi - CKf)*Solution2Soil/390

  #Conversion of final K concentration from mg/L to mol/L
  Ck_f_mol_L <- CKf/39000

  ARk <- (Ck_f_mol_L/CCaMg^0.5)

  #Concentration of Cl (mol/L)
  C_Cl <- Ck_f_mol_L + 2*CCaMg

  #Ionic strength calculation
  I = 0.5*(Ck_f_mol_L*1^2 + CCaMg*2^2 + C_Cl*1^2)

  #Activity coefficient calculation
  activity_coef <- exp((0.5*I^0.5)/(1 + I^0.5) - 0.086*C_Cl)

  #ARk calculation
  ARk <- (Ck_f_mol_L/CCaMg^0.5)*(activity_coef)

  #Plot deltak vs. ARk
  # plot(y = deltak, x = ARk,
  #      main='Quantity Intensity relationship',
  #      xlab=expression('AR'["K"]),
  #      ylab = expression(paste(Delta,"K", " cmol"["c"], " kg"^-1)))
  dat <- cbind.data.frame(deltak, ARk)
  #Fitting quadratic model
  #mod <- lm(deltak ~ ARk + I(ARk^2))
  mod <- lm(deltak ~ poly(ARk, 2, raw = TRUE))
  dat$pred = predict(mod, newdata = dat)
  # lines(sort(ARk), fitted(mod)[order(ARk)], col='red', type='l')
  # with(dat, lines(x = ARk, y = pred, col='red'))

  a <- mod$coefficients[[3]]
  b <- mod$coefficients[[2]]
  c <- mod$coefficients[[1]]
  #Value of ARK at deltak = 0
  AReK <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
  #Slope of the curve at deltak = 0
  PBC <- a*2*AReK + b

  #Plotting using ggplot2
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    myplot <- function(mydf, xcol, ycol){
      ggplot2::ggplot(data = mydf, aes(x = {{xcol}}, y = {{ycol}})) +
        ggplot2::geom_point(color='red', alpha=0.3, size=2)+
        ggplot2::stat_smooth(method='lm', formula = y~poly(x,2), se = F) +
        ggplot2::theme_bw() +
        ggplot2::xlab(expression(paste("AR"["K"]))) +
        ggplot2::ylab(expression(paste(Delta,"K", " cmol"["c"], " kg"^-1)))+
        ggplot2::theme(panel.grid = element_blank()) +
        ggplot2::theme(text=element_text(family = "serif", size=15),
                       axis.title.x = element_text(size = 18),
                       axis.title.y = element_text(size = 18),
                       axis.text.x = element_text(colour="black",face="bold"),
                       axis.text.y = element_text(colour="black",face="bold"))
    }
  } else {
    ## Please install the "ggplot2" package
  }

  my_plot <- myplot(dat, ARk, deltak)

  Kl <- abs(mod$coefficients[[1]])

  #deltaG0 = RTln(AReK)
  #T : 298 Kelvin
  #R: Universal gas constant
  deltaG0 <- 1.9858775*298*(log(AReK))

  return(list('PBC' = PBC, 'AReK' = AReK, 'Kl' = Kl, 'deltaG0' = deltaG0,
           'Plot' = my_plot))
}

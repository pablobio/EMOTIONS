#' Performs the model fitting and the weight assignment based on different strategies for each individual ID
#'
#' @param data A data frame containing the daily milking records
#' @param ID_col The name of the column containing the unique IDs of the individuals
#' @param ID The individual ID that is being analyzed
#' @param trait The name of the column containing daily milking records
#' @param dim  The name of the column containing days in milk records
#' @param alpha A penalization factor, ranging from 0 to 1, for the estimation of the model`s weight
#' @param models A vector describing the models to be included in the analysis. In total, 47 models are included in EmsembleLacs. The default option is "All", which results in the inclusion of the 47 models. Alternatively, a vector containing any subset of the following models can be provided: "MMR","MME","brody23","brody24", "SCH","SCHL","PBE","wood","DHA", "CB","QP","CLD","PapBo1","PapBo2", "PapBo3", "PapBo4", "PapBo6", "GS1",  "GS2","LQ", "wil", "wilk", "wilycsml", "BC", "DJK","MG2", "MG4", "MG", "KHN", "AS", "FRP","PTmult","PTmod", "MonoG", "MonoGpw", "DiG", "DiGpw","legpol3", "legpol4", "legpolWil", "cubsplin3", "cubsplin4", "cubsplin5", "cubsplindef", "wilminkPop", "qntReg", "Legpol4Poppe"
#' @importFrom stats predict AIC var sd acf as.formula dgamma lm median nls
#' @importFrom orthopolynom polynomial.values
#' @importFrom orthopolynom legendre.polynomials
#' @importFrom orthopolynom scaleX
#' @importFrom stats nls.control
#' @importFrom quantreg rq
#' @importFrom minpack.lm nlsLM
#' @return A list containing the fitted models, the model`s weigths and ranks, and the predicted daily production obtained through the model ensemble
ModelsLac<-function(data,ID_col,ID,trait,dim, alpha,models){

  x<-data

  out_list <- list(
    converged_models = list(),
    models_weight = list(),
    production = list()
  )

  #(1) Michaelis-Menten
  MM = try(nls(as.formula(paste0(trait, " ~ (a", " * ", dim, ") / (b", " + ", dim, ")")), data = x, na.action = na.exclude, start = list(a = 19.8, b = -1.65), control = list(maxiter = 100,warnOnly = TRUE)), silent = TRUE)


  #(1a) Michaelis-Menten (Rook)
  MMR = try(nls(as.formula(paste0(trait, "~",1, "/(",1,"+ a", "/(b","+",dim,"))")), data = x, na.action = na.exclude, start = list(a = -10, b = -1.4), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  # #(1b) Michaelis-Menten + Exponential (Rook)
  MME = try(nls(as.formula(paste0(trait, "~","a*(",1,"/(",1,"+b","/(c","+", dim,")))", "*","exp(-d*",dim,")")), data = x, na.action = na.exclude, start = list(a = -0.06608, b = 317.49, c = -328.06, d = -0.027), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(2) Brody 1923
  brody23 = try(nls(as.formula(paste0(trait, "~","a*exp(-b*",dim,")")), data = x, na.action = na.exclude, start = list(a = 25.6, b = 0.0015), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  # #(3) Brody 1924
  brody24 = try(nls(as.formula(paste0(trait, "~","a*exp(-b*",dim,")-","a*exp(-c*",dim,")")), data = x, na.action = na.exclude, start = list(a = 26.127, b = 0.0017, c = 0.2575), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(4) Schumacher
  SCH = try(nls(as.formula(paste0(trait, "~","exp(a+b/",dim,")")), data = x, na.action = na.exclude, start = list(a = 1, b = 16), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(4a) Schumacher (Lopez)
  SCHL = try(nls(as.formula(paste0(trait, "~","a*exp(b*",dim,"/(",dim,"+",1,"))")), data = x, na.action = na.exclude, start = list(a = min(x[,trait]), b = -2.6), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(5) Parabolic Exponential (Sikka, Adediran)
  PBE = try(nls(as.formula(paste0(trait, "~","a*exp((b*",dim,")-(c*",dim,"^2","))")), data = x, na.action = na.exclude, start = list(a = 25, b = -8e-04, c = 2e-06), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(6) Wood
  wood = try(nls(as.formula(paste0(trait, "~","a*",dim,"^b*exp(-c*",dim,")")), data = x, na.action = na.exclude, start = list(a = median(x[,trait]), b = log(max(x[,trait]) / median(x[,trait])) / max(x[,dim]), c = 0.05), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(6a) Wood (Dhanoa)
  DHA = try(nls(as.formula(paste0(trait, "~","a*(",dim,"^(b*c))*exp(-c*",dim,")")), data = x, na.action = na.exclude, start = list(a = median(x[,trait]), b = log(max(x[,trait]) / median(x[,trait])) / max(x[,dim]), c = 0.05), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(6b) Wood (Cappio-Borlino)
  CB = try(nls(as.formula(paste0(trait, "~","a*(",dim,"^b)*exp(-c*",dim,")")), data = x, na.action = na.exclude, start = list(a = min(x[,trait]), b = 0.07, c = 0.002), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(7) Quadratic Polynomial (Dave)
  QP = try(nls(as.formula(paste0(trait, "~","a+b*",dim,"+c*",dim,"^2")), data = x, na.action = na.exclude, start = list(a = 25, b = -0.02, c = -1.5e-05), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(8) Cobby & Le Du
  CLD = try(nls(as.formula(paste0(trait, "~","a-b*",dim,"-a*exp(-c*",dim,")")), data = x, na.action = na.exclude, start = list(a = 25, b = 0.03, c = 0.28), control = list(maxiter = 100,minFactor = 1e-6,  warnOnly = TRUE)), silent = TRUE)


  #(9) Papajcsik and Bodero 1
  PapBo1 = try(nls(as.formula(paste0(trait, "~","a*",dim,"^b/cosh(c*",dim,")")), data = x, na.action = na.exclude, start = list(a = min(x[,trait]), b = -0.0098, c = 0.0033), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(10) Papajcsik and Bodero 2
  PapBo2 = try(nls(as.formula(paste0(trait, "~","a*(1-exp(-b*",dim,"))/cosh(c*",dim,")")), data = x, na.action = na.exclude, start = list(a = 23.896, b = 0.398, c = 0.0034), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(11) Papajcsik and Bodero 3
  PapBo3 = try(nls(as.formula(paste0(trait, "~","a*atan(b*", dim,")/cosh(c*",dim,")")), data = x, na.action = na.exclude, start = list(a = 15.27, b = 2.587, c = 0.00346), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(12) Papajcsik and Bodero 4
  PapBo4 = try(nls(as.formula(paste0(trait, "~","a*log(b*", dim, ")*exp(-c*",dim,")")), data = x, na.action = na.exclude, start = list(a = 1.763, b = 76690, c = 0.00219), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(13) Papajcsik and Bodero 6
  PapBo6 = try(nls(as.formula(paste0(trait, "~","a*atan(b*",dim,")*exp(-c*", dim, ")")), data = x, na.action = na.exclude, start = list(a = 17.25, b = 0.493, c = 0.0018), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(14) Mixed log model 1 (Guo & Swalve)
  GS1 = try(nls(as.formula(paste0(trait, "~","a+b*",dim,"^0.5+c*log(", dim,")")), data = x, na.action = na.exclude, start = list(a = 18.28, b = -1.58, c = 4.33), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(15) Mixed log model 3 (Guo & Swalve)
  GS2 = try(nls(as.formula(paste0(trait, "~","a+b*",dim,"^0.5+c*log(", dim,")+d*", dim, "^4")), data = x, na.action = na.exclude, start = list(a = 17.75, b = -1.72, c = 4.81, d = 7e-11), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  # #(16) Log-quadratic (Adediran)
  LQ = try(nls(as.formula(paste0(trait, "~","exp(a*(b-log(", dim, "))^2+c)")), data = x, na.action = na.exclude, start = list(a = min(x[,trait]), b = (max(x[,trait], na.rm = TRUE) - min(x[,trait], na.rm = TRUE)) / max(x[,dim], na.rm = TRUE), c = (max(x[,trait], na.rm = TRUE) - min(x[,trait], na.rm = TRUE)) / 2), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(17) Wilmink
  wil = try(nls(as.formula(paste0(trait, "~","a+b*exp(-k*",dim,")+c*",dim)), data = x, na.action = na.exclude, start = list(a = 25, b = -7, c = -0.03, k = 0.1), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(17a) Wilmink (Kheidirabadi)
  wilk = try(nls(as.formula(paste0(trait, "~","a+b*exp(-k*",dim,")+c*(",dim,"/100)")), data = x, na.action = na.exclude, start = list(a = 25, b = -7, c = -3, k = 0.1), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(17b) Wilmink (Laurenson & Strucken)
  wilycsml = try(nls(as.formula(paste0(trait, "~","a+(b-a)*(1-exp(-k*",dim,"))-c*",dim)), data = x, na.action = na.exclude, start = list(a = 20, b = 30, c = 0.005, k = 0.08), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(18) Bicompartemental (Ferguson & Boston)
  BC = try(nls(as.formula(paste0(trait, "~","a*exp(-b*", dim,")+c*exp(-d*",dim,")")), data = x, na.action = na.exclude, start = list(a = min(x[,trait], na.rm = TRUE), b =0.0017, c = -7.68, d = 0.08), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(19) Dijkstra
  DJK = try(nls(as.formula(paste0(trait, "~","a*exp(b*(1-exp(-c*",dim,"))/c-(d*",dim,"))" )), data = x, na.action = na.exclude, start = list(a = 19, b = 0.027, c = 0.08, d = 0.0017), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(20) Morant & Gnanasakthy (Pollott)
  MG2 = try(nls(as.formula(paste0(trait, "~","exp(a+b*((",dim,"-","150)/100)+c*((",dim,"-150)/100)^2+d/",dim,")")), data = x, na.action = na.exclude, start = list(a = 3, b = -0.18, c = 0.002, d = -1.4), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(21) Morant & Gnanasakthy (Vargas)
  MG4 = try(nls(as.formula(paste0(trait, "~","a*exp(b*((",dim,"-150)/100)/2+c/",dim,"-d*(1+((",dim,"-21.4)/100)/2)*((",dim,"-21.4)/100))")), data = x, na.action = na.exclude, start = list(a = 20.5, b = -0.38, c = -1.45, d = -0.005), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(22) Morant & Gnanasakthy (Adediran)
  MG = try(nls(as.formula(paste0(trait, "~","exp(a-b*((",dim,"-150)/100)+c*((",dim,"-150)/100)^2+d/",dim,")")), data = x, na.action = na.exclude, start = list(a = 3, b = 0.18, c = 0.002, d = -1.4), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(23) Khandekar (Guo & Swalve)
  KHN = try(nls(as.formula(paste0(trait, "~","a+b*",dim,"+c*",dim,"^2+d*",dim,"^3+f*log(",dim,")")), data = x, na.action = na.exclude, start = list(a = 15, b = -0.15, c = 0.00043, d = 5e-07, f = 4.05), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(24) Ali & Schaeffer
  AS = try(nls(as.formula(paste0(trait, "~","a+b*(",dim,"/340)+c*(",dim,"^2/340)+d*log(340/",dim,")+f*log(340/",dim,")^2")), data = x, na.action = na.exclude, start = list(a = 19, b = -5, c = 0.003, d = 5.3, f = -1.1), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(25) Fractional Polynomial (Elvira)
  FRP = try(nls(as.formula(paste0(trait, "~","a+b*",dim,"+c*log(",dim,")+d*",dim,"^0.5+f*",dim,"^2")), data = x, na.action = na.exclude, start = list(a = 15, b = -0.08, c = 8, d = -3, f = -4e-05), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(26) Pollott multiplicative reduced (Elvira)
  PTmult = try(nls(as.formula(paste0(trait, "~","a/(1+((1 - 0.999999)/0.999999) *exp(-0.1 *",dim,"))*(2-exp(-b*",dim,"))")), data = x, na.action = na.exclude, start = list(a = 25, b = -0.001), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(27) Pollott modified (Adediran)
  PTmod = try(nls(as.formula(paste0(trait, "~","a/((1+b*(exp(-c*",dim,")))) * (2 -exp(-d*",dim,"))")), data = x, na.action = na.exclude, start = list(a = min(x[,trait],na.rm = TRUE), b = (max(x[,trait], na.rm = TRUE) - min(x[,trait], na.rm = TRUE)) / max(x[,dim], na.rm = TRUE), c = 0.13, d = -0.001), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(28) Monophasic Grossman
  MonoG = try(nls(as.formula(paste0(trait, "~","a*b*(1-tanh(b*(",dim,"-c)))")), data = x, na.action = na.exclude, start = list(a = 7703, b = 0.002, c = 277), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(29) Monophasic Grossman power
  MonoGpw = try(nls(as.formula(paste0(trait, "~","a*(1 - tanh(b*(",dim,"^k-c))^2)")), data = x, na.action = na.exclude, start = list(a = 24.5707, b = 0.6292, c = 1.9977, k = 0.1978), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(30) Diphasic Grossman
  DiG = try(nls(as.formula(paste0(trait, "~","a*b*(1 - tanh(b*(",dim,"-c))^2) +d*f*(1 - tanh(f*(",dim,"-g))^2)")), data = x, na.action = na.exclude, start = list(a = 353.8, b = 0.016, c = 36.3, d = 7371, f = 0.003, g = 113.5), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(31) Diphasic Grossman power
  DiGpw = try(nls(as.formula(paste0(trait, "~","a*b*(1 - tanh(b*(",dim,"^k-c))^2) +d*f*(1 - tanh(f*(",dim,"-g))^2)")), data = x, na.action = na.exclude, start = list(a = 15.99, b = 0.3351, c = 4.693, d = 8687, f = 0.00226, g = 80.33, k = 0.435), control = list(maxiter = 100, warnOnly = TRUE)), silent = TRUE)


  #(32) Legendre Polynomial 3th (Jacobson)
  y <- 1:nrow(x)
  leg3 <- as.matrix(as.data.frame(polynomial.values(polynomials = legendre.polynomials(n = 3, normalized = TRUE), x = scaleX(y, u = -1, v = 1))))
  colnames(leg3) <- c("leg0", "leg1", "leg2", "leg3")
  leg3 <- leg3[, 2:ncol(leg3)]
  legpol3 = try(lm(as.formula(paste0(trait, "~","leg3")), data = x, na.action = na.exclude),
                silent = TRUE)


  #(33) Legendre Polynomial 4th (Jacobson)
  y <- 1:nrow(x)
  leg4 <- as.matrix(as.data.frame(polynomial.values(polynomials = legendre.polynomials(n = 4,normalized = TRUE), x = scaleX(y, u = -1, v = 1))))
  colnames(leg4) <- c("leg0", "leg1", "leg2", "leg3", "leg4")
  leg4 <- leg4[, 2:ncol(leg4)]
  legpol4 = try(lm(as.formula(paste0(trait, "~","leg4")), data = x, na.action = na.exclude),
                silent = TRUE)



  #(34) Legendre + Wilmink (Lindauer)
  y <- 1:nrow(x)
  leg4 <- as.matrix(as.data.frame(polynomial.values(polynomials = legendre.polynomials(n = 4, normalized = TRUE), x = scaleX(y, u = -1, v = 1))))
  colnames(leg4) <- c("leg0", "leg1", "leg2", "leg3", "leg4")
  leg4 <- leg4[, 2:ncol(leg4)]
  legpolWil = try(nls(as.formula(paste0(trait, "~","a*","leg4[,1]","+b*","leg4[,2]","+c*","leg4[,3]","+d*exp(-k*",dim,")")), data = x, na.action = na.exclude, start = list(a = -0.8, b = -0.6, c = 0.1, d = 25.7, k = 0.002)), silent = TRUE)


  #(35) Natural Cubic Spline (3 percentiles)
  cubsplin3 = try(lm(as.formula(paste0(trait, "~","ns(",dim, ",df = 4)")), data = x, na.action = na.exclude), silent = TRUE)

  #(36) Natural Cubic Spline (4 percentiles)
  cubsplin4 = try(lm(as.formula(paste0(trait, "~","ns(",dim, ",df = 5)")), data = x, na.action = na.exclude), silent = TRUE)

  #(37) Natural Cubic Spline (5 percentiles)
  cubsplin5 = try(lm(as.formula(paste0(trait, "~","ns(",dim, ",df = 6)")), data = x, na.action = na.exclude),silent = TRUE)

  #(38) Natural Cubic Spline (defined Harrell)
  cubsplindef = try(lm(as.formula(paste0(trait, "~","ns(",dim, ",knots = c(49, 78, 112, 157, 210))")), data = x, na.action = na.exclude), silent = TRUE)


  ###Wilmink Poppe

  start_values <- list(
    b0 = min(x[,trait], na.rm = TRUE),
    b1 = (max(x[,trait], na.rm = TRUE) - min(x[,trait], na.rm = TRUE)) / max(x[,dim], na.rm = TRUE),
    b2 = (max(x[,trait], na.rm = TRUE) - min(x[,trait], na.rm = TRUE)) / 2
  )

  wilminkPop <- try(nls(as.formula(paste0(trait, "~","b0+b1*",dim,"+b2*exp(-0.5*",dim,")")),
                        data = x,
                        start = start_values,
                        algorithm = "port",
                        control = nls.control(maxiter = 500)),silent=TRUE)

  ###Quantile regression

  qntReg <- try(rq(as.formula(paste0(trait, "~","poly(",dim,",4,raw=TRUE)")), tau = 0.7, data = x),silent=TRUE)


  ###Legendre 4

  # Crear las funciones de Legendre de cuarto grado
  legendre_poly <- function(x) {
    p0 <- rep(1, length(x))
    p1 <- x
    p2 <- 0.5 * (3 * x^2 - 1)
    p3 <- 0.5 * (5 * x^3 - 3 * x)
    p4 <- 0.125 * (35 * x^4 - 30 * x^2 + 3)
    return(data.frame(p0, p1, p2, p3, p4))
  }

  # Aplicar las funciones a tu columna DIM
  legendre_vars <- legendre_poly(x[,dim])

  # Definir la función para el modelo
  legendre_function <- function(b0, b1, b2, b3, b4, p0, p1, p2, p3, p4) {
    b0 * p0 + b1 * p1 + b2 * p2 + b3 * p3 + b4 * p4
  }


  # Combinar la tabla original con las funciones de Legendre
  d_legendre <- cbind(x, legendre_vars)

  # Ajustar el modelo con nlsLM
  Legpol4Poppe <- try(nlsLM(as.formula(paste0(trait, "~","legendre_function(b0, b1, b2, b3, b4, p0, p1, p2, p3, p4)")),
                            data = d_legendre,
                            start = list(b0 = 0.5, b1 = 0.5, b2 = 0.5, b3 = 0.5, b4 = 0.5)),silent=TRUE)


  all_objects <- c("MMR","MME","brody23","brody24",
                   "SCH","SCHL","PBE","wood","DHA",
                   "CB","QP","CLD","PapBo1","PapBo2",
                   "PapBo3", "PapBo4", "PapBo6", "GS1",
                   "GS2","LQ", "wil", "wilk", "wilycsml",
                   "BC", "DJK","MG2", "MG4", "MG", "KHN",
                   "AS", "FRP","PTmult","PTmod", "MonoG",
                   "MonoGpw", "DiG", "DiGpw","legpol3",
                   "legpol4", "legpolWil", "cubsplin3",
                   "cubsplin4", "cubsplin5", "cubsplindef",
                   "wilminkPop", "qntReg", "Legpol4Poppe")

  if("All"%in%models){
    all_objects<-all_objects
  }else{
    all_objects<-all_objects[which(all_objects%in%models)]
  }

  # Retrieve the objects by their names from the environment
  object_list <- mget(all_objects, ifnotfound = NA)

  # Filter for objects that are of class 'lm', 'nls', or 'rq'
  model_objects <- names(object_list)[sapply(object_list, function(obj) {
    class(obj) %in% c("lm", "nls", "rq")
  })]


  # Create a list of model objects
  model_list <- mget(model_objects)


  # Assuming `model_list` is your list of model objects
  converged_models <- Filter(function(model) {
    # Check if model is from lm, nls, or rq
    if (inherits(model, "lm")) {
      # lm models generally don't have convergence issues
      return(TRUE)
    } else if (inherits(model, "nls")) {
      # Check convergence for nls models
      return(model$convInfo$isConv)
    } else if (inherits(model, "rq")) {
      # Check convergence for rq models
      return(model$converged)
    } else {
      # If it's not one of these types, assume non-converged
      return(FALSE)
    }
  }, model_list)


  out_list$converged_models[[ID]] <- converged_models

  model.metrics<-GetLacModelsMetrics(converged_models, x, trait)


  ########
  ####Estimating delta for the metrics
  ########
  model.metrics$delta_AIC <- model.metrics$AIC - min(model.metrics$AIC)
  model.metrics$delta_BIC <- model.metrics$BIC - min(model.metrics$BIC)
  model.metrics$delta_RMSPE <- model.metrics$RMSPE - min(model.metrics$RMSPE)
  model.metrics$delta_MAE <- model.metrics$MAE - min(model.metrics$MAE)


  ##########
  ####Calculating weights
  #########

  # Calculate the weights for RMSPE and MAE using the exponential formula
  model.metrics$AIC_weight <- exp(-alpha * model.metrics$delta_AIC) / sum(exp(-alpha * model.metrics$delta_AIC))

  model.metrics$BIC_weight <- exp(-alpha * model.metrics$delta_BIC) / sum(exp(-alpha * model.metrics$delta_BIC))

  model.metrics$RMSPE_weight <- exp(-alpha * model.metrics$delta_RMSPE) / sum(exp(-alpha * model.metrics$delta_RMSPE))

  model.metrics$MAE_weight <- exp(-alpha * model.metrics$delta_MAE) / sum(exp(-alpha * model.metrics$delta_MAE))


  # Weight for variance
  model.metrics$Var_weight <- VarWeight(converged_models, x)

  model.metrics$CosSquared_weight<-CosSquaredWeight(converged_models, x)

  model.metrics$BAM_weight<-BMAweight_gamma(converged_models, x,trait)

  #########
  ##### Ordering and creating sequence
  #########
  model.metrics<-model.metrics[order(model.metrics$RMSPE, decreasing = F),]
  model.metrics$RMSPE_rank<-seq(1,nrow(model.metrics),1)

  model.metrics<-model.metrics[order(model.metrics$MAE, decreasing = F),]
  model.metrics$MAE_rank<-seq(1,nrow(model.metrics),1)

  model.metrics<-model.metrics[order(model.metrics$AIC, decreasing = F),]
  model.metrics$AIC_rank<-seq(1,nrow(model.metrics),1)

  model.metrics<-model.metrics[order(model.metrics$BIC, decreasing = F),]
  model.metrics$BIC_rank<-seq(1,nrow(model.metrics),1)

  model.metrics<-model.metrics[order(model.metrics$BAM_weight, decreasing = T),]
  model.metrics$BAM_rank<-seq(1,nrow(model.metrics),1)

  model.metrics<-model.metrics[order(model.metrics$Var_weight, decreasing = T),]
  model.metrics$Var_rank<-seq(1,nrow(model.metrics),1)

  model.metrics<-model.metrics[order(model.metrics$CosSquared_weight, decreasing = T),]
  model.metrics$CosSquared_rank<-seq(1,nrow(model.metrics),1)


  out_list$models_weight[[ID]] <- model.metrics


  weight.AIC.out<-matrix(NA,nrow=length(x[,dim]),ncol=length(names(converged_models)), dimnames = list(x[,dim], names(converged_models)))

  weight.BIC.out<-matrix(NA,nrow=length(x[,dim]),ncol=length(names(converged_models)), dimnames = list(x[,dim], names(converged_models)))

  weight.RMSPE.out<-matrix(NA,nrow=length(x[,dim]),ncol=length(names(converged_models)), dimnames = list(x[,dim], names(converged_models)))

  weight.MAE.out<-matrix(NA,nrow=length(x[,dim]),ncol=length(names(converged_models)), dimnames = list(x[,dim], names(converged_models)))

  weight.BAM.out<-matrix(NA,nrow=length(x[,dim]),ncol=length(names(converged_models)), dimnames = list(x[,dim], names(converged_models)))

  weight.Var.out<-matrix(NA,nrow=length(x[,dim]),ncol=length(names(converged_models)), dimnames = list(x[,dim], names(converged_models)))

  weight.CosSquared.out<-matrix(NA,nrow=length(x[,dim]),ncol=length(names(converged_models)), dimnames = list(x[,dim], names(converged_models)))

  SMA.out<-matrix(NA,nrow=length(x[,dim]),ncol=length(names(converged_models)), dimnames = list(x[,dim], names(converged_models)))


  for(model_name in names(converged_models)){
    model <- converged_models[[model_name]]

    # Get predictions for the model
    predictions <- predict(model)  # Replace `data` with your actual data

    # Get the corresponding weights for AIC, RMSPE, and MAE
    weight.AIC <- model.metrics$AIC_weight[model.metrics$Model == model_name]

    if(weight.AIC<0.01){
      weight.AIC<-0
    }

    # Get the corresponding weights for AIC, RMSPE, and MAE
    weight.BIC <- model.metrics$BIC_weight[model.metrics$Model == model_name]

    if(weight.BIC<0.01){
      weight.BIC<-0
    }


    rmspe_weight <- model.metrics$RMSPE_weight[model.metrics$Model == model_name]

    if(rmspe_weight<0.01){
      rmspe_weight<-0
    }

    mae_weight <- model.metrics$MAE_weight[model.metrics$Model == model_name]

    if(mae_weight<0.01){
      mae_weight<-0
    }

    bam_weight <- model.metrics$BAM_weight[model.metrics$Model == model_name]

    if(bam_weight<0.01){
      bam_weight<-0
    }

    var_weight<- model.metrics$Var_weight[model.metrics$Model == model_name]

    if(var_weight<0.01){
      var_weight<-0
    }


    CosSquared_weight<-model.metrics$CosSquared_weight[model.metrics$Model == model_name]

    if(CosSquared_weight<0.01){
      CosSquared_weight<-0
    }

    weight.AIC.out[,model_name]<-predictions*weight.AIC
    weight.BIC.out[,model_name]<-predictions*weight.BIC
    weight.RMSPE.out[,model_name]<-predictions*rmspe_weight
    weight.MAE.out[,model_name]<-predictions*mae_weight
    weight.BAM.out[,model_name]<-predictions*bam_weight
    weight.Var.out[,model_name]<-predictions*var_weight
    weight.CosSquared.out[,model_name]<-predictions*CosSquared_weight
    SMA.out[,model_name]<-predictions

  }


  weighted_sum<-matrix(NA,nrow=nrow(x),ncol=8)


  colnames(weighted_sum)<-c("weight_AIC","weight_BIC","Weight_RMSPE","weight_MAE","weight_BAM", "weight_Var","weight_CosSquared","SMA")

  weighted_sum[,"weight_AIC"]<-rowSums(weight.AIC.out)

  weighted_sum[,"weight_BIC"]<-rowSums(weight.BIC.out)

  weighted_sum[,"Weight_RMSPE"]<-rowSums(weight.RMSPE.out)

  weighted_sum[,"weight_MAE"]<-rowSums(weight.MAE.out)

  weighted_sum[,"weight_BAM"]<-rowSums(weight.BAM.out)

  weighted_sum[,"weight_Var"]<-rowSums(weight.Var.out)

  weighted_sum[,"weight_CosSquared"]<-rowSums(weight.CosSquared.out)

  weighted_sum[,"SMA"]<-rowMeans(SMA.out)

  weighted_sum<-cbind(data[,c(trait,dim)],weighted_sum)

  out_list$production[[ID]] <- weighted_sum

  return(out_list)

}

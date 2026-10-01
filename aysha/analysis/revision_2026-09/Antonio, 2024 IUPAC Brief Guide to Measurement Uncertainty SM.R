######################################################################

## Supplementary Material
## A brief guide to measurement uncertainty (IUPAC technical report)
## A. Possolo, D. B. Hibbert, J. Stohner, O. Bodnar, J. Meija

## REVISION DATE: 2023-Jan-04

## Computer codes written in the R programming language. The R
## environment for statistical modeling and data analysis is freely
## available for download from https://www.r-project.org/, with
## supporting documentation at https://cran.r-project.org/manuals.html

## The copious comments interspersed with the code do not aim to teach
## the R language, which the reader is assumed to be familiar with
## already.

## Once the required add-on packages will have been installed, this R
## code should execute if it is submitted via the R console, under
## "Emacs Speaks Statistics" (ESS) in emacs, or in the RStudio
## integrated development environment that is available for free
## download from https://posit.co/products/open-source/rstudio/

## References appear between square brackets, and their numbers are as
## listed under "References" in the main article. Any reference not in
## the main article is fully specified where it is made.

######################################################################
##
## EXAMPLE 3E (MOLAR MASS OF IRIDIUM -- LABORATORY EFFECTS MODEL)
##
######################################################################

## Chang et al. (1992)   : Ir191/Ir193 = 0.59399, u(Ir191/Ir193) = 0.00103
## Walczyk et al. (1993) : Ir191/Ir193 = 0.59418, u(Ir191/Ir193) = 0.00037
## Zhu et al. (2017)     : Ir191/Ir193 = 0.59290, u(Ir191/Ir193) = 0.00021

r.C = 0.59399; ur.C = 0.00103 ## mol/mol
r.W = 0.59418; ur.W = 0.00037
r.Z = 0.59290; ur.Z = 0.00021

## AME 2020 Atomic Mass Evaluation (Wang et al. 2021)
Ir191.x = 190960591.5 / 1e6 ## dalton
Ir191.u =         1.4 / 1e6
Ir193.x = 192962923.8 / 1e6
Ir193.u =         1.4 / 1e6

## Isotopic amount fractions corresponding to isotope ratio C,
## corresponding estimate of the molar mass of iridium, and evaluation
## of the associated uncertainty by application of the Monte Carlo method

Ir191.C.x = r.C/(1+r.C)
Ir193.C.x = 1-Ir191.C.x

K = 1e7
rB = rnorm(K, mean=r.C, sd=ur.C)
Ir191.C.xB = rB/(1+rB)
Ir193.C.xB = 1 - Ir191.C.xB
Ir.C.x = (Ir191.C.x*Ir191.x + Ir193.C.x*Ir193.x)
Ir191.xB = rnorm(K, mean=Ir191.x, sd=Ir191.u)
Ir193.xB = rnorm(K, mean=Ir193.x, sd=Ir193.u)
Ir.C.u = sd(Ir191.C.xB*Ir191.xB + Ir193.C.xB*Ir193.xB)

## Isotopic amount fractions corresponding to isotope ratio W,
## corresponding estimate of the molar mass of iridium, and evaluation
## of the associated uncertainty by application of the Monte Carlo method

Ir191.W.x = r.W/(1+r.W)
Ir193.W.x = 1-Ir191.W.x

K = 1e7
rB = rnorm(K, mean=r.W, sd=ur.W)
Ir191.W.xB = rB/(1+rB)
Ir193.W.xB = 1 - Ir191.W.xB
Ir.W.x = (Ir191.W.x*Ir191.x + Ir193.W.x*Ir193.x)
Ir191.xB = rnorm(K, mean=Ir191.x, sd=Ir191.u)
Ir193.xB = rnorm(K, mean=Ir193.x, sd=Ir193.u)
Ir.W.u = sd(Ir191.W.xB*Ir191.xB + Ir193.W.xB*Ir193.xB)

## Isotopic amount fractions corresponding to isotope ratio Z,
## corresponding estimate of the molar mass of iridium, and evaluation
## of the associated uncertainty by application of the Monte Carlo method

Ir191.Z.x = r.Z/(1+r.Z)
Ir193.Z.x = 1-Ir191.Z.x

K = 1e7
rB = rnorm(K, mean=r.Z, sd=ur.Z)
Ir191.Z.xB = rB/(1+rB)
Ir193.Z.xB = 1 - Ir191.Z.xB
Ir.Z.x = (Ir191.Z.x*Ir191.x + Ir193.Z.x*Ir193.x)
Ir191.xB = rnorm(K, mean=Ir191.x, sd=Ir191.u)
Ir193.xB = rnorm(K, mean=Ir193.x, sd=Ir193.u)
Ir.Z.u = sd(Ir191.Z.xB*Ir191.xB + Ir193.Z.xB*Ir193.xB)

## Estimates of the molar mass of iridium corresponding to the C, W,
## and Z determinations of the isotope ratio, and associated standard
## uncertainties

options(digits=8)
c(Ir.C.x, Ir.W.x, Ir.Z.x)
## 192.21677 192.21662 192.21763

round(c(Ir.C.u, Ir.W.u, Ir.Z.u), 5)
## 0.00081 0.00029 0.00017

## Naive evaluation of dark uncertainty: excess variability of
## measured values, above and beyond what the standard uncertainties
## associated with them suggest that variability should be

## Standard deviation of the C, W, and Z estimates of the molar mass
s = sd(c(Ir.C.x, Ir.W.x, Ir.Z.x))             ## 0.00054
## Geometric average of the associated standard uncertainties
g = exp(mean(log(c(Ir.C.u, Ir.W.u, Ir.Z.u)))) ## 0.00034
s/g ## 1.6

## Naive evaluation of tau (dark uncertainty)
sqrt(s^2 - g^2)             ## 0.00043 g/mol
sqrt(0.00054^2 - 0.00034^2) ## 0.00042

## Cochran's [20] Q test of mutual consistency of the measurement
## results and DerSimonian-Laird evaluation of dark uncertainty
require(metafor)
rma(yi=c(Ir.C.x, Ir.W.x, Ir.Z.x), sei=c(Ir.C.u, Ir.W.u, Ir.Z.u), method="DL")
## Test for Heterogeneity: Q(df = 2) = 9.61, p-val = 0.0082

## Since the p-value is very small, Cochran's Q test rejects the
## hypothesis of mutual consistency of the two measurement
## results. The larger the Q criterion (9.61 in this case), the
## stronger the evidence against the hypothesis of mutual consistency

## The p-value of the statistical test of the hypothesis of mutual
## consistency is the probability of observing a value of the test
## criterion (Q in this case) more extreme than the value that was
## observed, when the measurement results in fact are mutually
## consistent

## The estimate of dark uncertainty depends on the estimation method
rma(yi=c(Ir.C.x, Ir.W.x, Ir.Z.x), sei=c(Ir.C.u, Ir.W.u, Ir.Z.u), method="DL")
## tau (square root of estimated tau^2 value):      0.0006
rma(yi=c(Ir.C.x, Ir.W.x, Ir.Z.x), sei=c(Ir.C.u, Ir.W.u, Ir.Z.u), method="REML")
## tau (square root of estimated tau^2 value):      0.0007
rma(yi=c(Ir.C.x, Ir.W.x, Ir.Z.x), sei=c(Ir.C.u, Ir.W.u, Ir.Z.u), method="ML")
## tau (square root of estimated tau^2 value):      0.0004

## Compare with the naive estimate above:           0.00041

######################################################################
##
## EXAMPLE 4B (ISOTOPES OF SILICON)
##
######################################################################

## Valkiers et al. [27, Table 4] -- Mixture III

R2928.x = 0.0500657 ## mol/mol
R2928.u = 0.0000025

R3028.x = 0.0329035
R3028.u = 0.0000034

## Relations between isotope ratios and isotopic composition: x28,
## x29, and x30 denote amount fractions of the three stable isotopes
## of silicon
x28 =       1 / (1 + R2928.x + R3028.x)
x29 = R2928.x / (1 + R2928.x + R3028.x)
x30 = R3028.x / (1 + R2928.x + R3028.x)

## Amount fractions are non-negative and add to 1
any(c(x28, x29, x30) < 0) ## FALSE
sum(c(x28, x29, x30)) ## 1

J = x28^2 * rbind(c(-1, 1), c(1+R3028.x, -R2928.x), c(-R3028.x, 1+R2928.x))
print(J)
## -0.85264410  0.85264410
##  0.88069907 -0.04268822
## -0.02805498  0.89533232

## RHO is the correlation between the two ratios, R2928.x and R3028.x,
## which we assume to be 0
rho = 0

sigmaR = rbind(c(R2928.u^2, rho*R2928.u*R3028.u),
               c(rho*R2928.u*R3028.u, R3028.u^2))
## Covariance matrix of the amount fractions (x28, x29, x30)
sigmaX = J %*% sigmaR %*% t(J)

cbind(round(c(x28=x28, x29=x29, x30=x30), 7),
      round(sqrt(diag(sigmaX)), 7))
## x28 0.9233873 3.6e-06 ## mol/mol
## x29 0.0462300 2.2e-06
## x30 0.0303827 3.0e-06

######################################################################
##
## EXAMPLE 4D (MODELING TAIL HEAVINESS)
##
######################################################################

## Laplace
require(extraDistr)
## Standard deviation of Laplace distribution with scale sigma is
## sqrt(2)*sigma 
sd(rlaplace(1e6, mu=0, sigma=1/sqrt(2))) ## 1
2 * plaplace(-3, mu=0, sigma=1/sqrt(2))  ## 0.0143696

## Gaussian
2 * pnorm(-3, mean=0, sd=1) ## 0.002699796

## Student's t[3]
## Standard deviation of Student's t with nu degrees of freedom is
## sqrt(nu/(nu-2)) for nu > 2
sd(rt(1e6, df=3)/sqrt(3/(3-2))) ## 1
2 * pt(-3 * sqrt(3/(3-2)), df=3) ## 0.01384683

######################################################################
##
## EXAMPLE 4E (SiO2 IN LIMESTONE)
##
######################################################################

## Maximum likelihood estimation of the mean and variance of
## a Gaussian distribution, given a sample of size 2 from such
## distribution. The maximum likelihood estimates of mu
## and sigma2 can be computed analytically, but below we use the
## general approach, which involves numerical optimization of the
## likelihood function

w.SiO2 = c(0.1811, 0.1818) ## g/g

negLogLik = function (theta, w, tol=sqrt(.Machine$double.eps))
{
    mu = theta[1]; sigma2 = theta[2]
    if (sigma2 < tol) { return(Inf)
    } else { return(-1*sum(dnorm(w, mean=mu, sd=sqrt(sigma2), log=TRUE))) }
}

SiO2.optim = optim(par=c(mean(w.SiO2), var(w.SiO2)), fn=negLogLik,
                   method="Nelder-Mead", w=w.SiO2)
signif(SiO2.optim$par[1], 5) ## 0.18145 ## g/g
signif(SiO2.optim$par[2], 2) ## 1.2e-7  ## (g/g)^2

## The maximum likelihood estimate (MLE) of the mean is the average of the
## measured values, but the MLE of the standard deviation is different
## from the sample standard deviation
signif(mean(w.SiO2), 5)                    ## 0.18145 ## g/g
signif(var(w.SiO2), 2)                     ## 2.4e-7  ## (g/g)^2
signif(mean((w.SiO2 - mean(w.SiO2))^2), 2) ## 1.2e-7  ## (g/g)^2

## NOTE: The maximum likelihood estimate of sigma is the square root
## of the maximum likelihood estimate of sigma2

######################################################################
##
## EXAMPLE 4F (INTERLABORATORY COMPARISON)
##
######################################################################

## The following measurement results of the mass fraction of nickel in
## bovine liver differ from those in the Final Report of CCQM-K145
## [34] only in that the uncertainty reported by INMC has been rounded
## to 2 significant digits

Ni = read.table(header=TRUE, text="
   lab     w    uw
   JSI 1.930 0.110
  INMC 1.940 0.064
  GLHK 1.942 0.092
   LNE 1.958 0.075
  NIST 1.984 0.020
 KRISS 1.993 0.033
INACAL 2.010 0.060
  NMIA 2.020 0.050
   NIM 2.022 0.023
  NMIJ 2.050 0.020
  RISE 2.055 0.052
   NRC 2.070 0.050
   PTB 2.077 0.035
  LATU 2.080 0.059
   LGC 2.131 0.042
   UME 2.150 0.030
   HSA 2.180 0.080") ## mg/kg

## Maximum likelihood estimate of the mean mu and dark uncertainty tau

require(bbmle)

negLogLik2 = function (mu, tau, w, uw)
{
    return(-1*sum(dnorm(w, mean=mu, sd=sqrt(uw^2 + tau^2), log=TRUE)))
}

## tauNaive is a naive estimate of the "excess" dispersion of the
## measured values, by comparison with the reported uncertainties,
## which are summarized by their geometric average
tauNaive = sqrt(var(Ni$w) - exp(mean(2*log(Ni$uw))))
Ni.mle = mle2(negLogLik2, start=list(mu=mean(Ni$w), tau=tauNaive), 
              method="L-BFGS-B", lower=c(mu=min(Ni$w), tau=0),
              skip.hessian=FALSE, data=list(w=Ni$w, uw=Ni$uw))

summary(Ni.mle)
##     Estimate Std. Error
## mu  2.043    0.015      ## mg/kg
## tau 0.044    0.014

## The covariance matrix of the estimates of mu and tau is obtained
## using the large-sample approximation based on the Hessian of the
## log-likelihood function (requested by specifying
## "skip.hessian=FALSE" above), even if we have only 17 measurement
## results
Ni.mle@vcov
##                mu           tau
## mu   2.384880e-04 -7.760112e-06
## tau -7.760112e-06  2.027808e-04

## Assuming the the Gaussian model is appropriate for these measured
## values, the restricted maximum likelihood (REML) estimator is
## generally preferable to the MLE. 

require(metafor)
rma(yi=w, sei=uw, data=Ni, method="REML", level=95, test="knha")

## tau = 0.0463 ## mg/kg
## Test for Heterogeneity: Q(df = 16) = 40.1205, p-val = 0.0007
## estimate      se    ci.lb   ci.ub      
##   2.0428  0.0159   2.0091  2.0765 ## mg/kg

## The interval whose endpoints are ci.lb and ci.ub is a 95 % coverage
## interval for the true value of the mass fraction of nickel in the
## material, obtained using the Knapp-Hartung adjustment, requested by
## specifying 'test="knha"' above.

## Knapp, G., & Hartung, J. (2003). Improved tests for a random
## effects meta-regression with a single covariate. Statistics in
## Medicine, 22(17), 2693-2710. DOI 10.1002/sim.1482

## IntHout, J., Ioannidis, J. P. A., & Borm, G. F.  (2014).  The
## Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis
## is straightforward and considerably outperforms the standard
## DerSimonian-Laird method.  BMC Medical Research Methodology, 14,
## 25, DOI 10.1186/1471-2288-14-25

######################################################################
##
## EXAMPLE 4G (MOLNUPIRAVIR)
##
######################################################################

## By the end of the trial reported by Butler et al. (2023), 105 among
## 12529 patients who had been on molnupiravir plus usual care, and 98
## among 12525 patients in the usual care group, had died or had been
## hospitalized

M.events   =  48
M.patients = 709

P.events   =  68
P.patients = 699

M.p = M.events/M.patients ## 0.06770099
P.p = P.events/P.patients ## 0.09728183

OR = (M.p/(1-M.p)) / (P.p/(1-P.p)) ##  0.6738453
logOR = log(OR)                    ## -0.3947547

logOR.u = sqrt(1/(M.patients*(M.p*(1-M.p))) + 1/(P.patients*(P.p*(1-P.p))))
## 0.1965626
## Note that OR.u yields the same value as the classical approximation
## for the standard error of the log odds ratio
sqrt(1/M.events + 1/(M.patients-M.events) +
     1/P.events + 1/(P.patients-P.events)) ## 0.1965626

## Z-score
logOR / logOR.u ## -2.00829
pnorm(logOR / logOR.u, mean=0, sd=1) ## 0.02230626

## An alternative analysis of the results of the clinical trial is
## based on logistic regression, which is a particular kind of a
## generalized linear model.

## A patient's "status" describes whether the patient is taking
## molnupiravir (M) or a placebo (P), and the corresponding "outcome"
## is 1 if the patient has been hospitalized or died, or 0 otherwise

status = relevel(factor(c(rep("M", M.patients),
                          rep("P", P.patients))), ref="P")
outcome = c(rep(1, M.events), rep(0, M.patients-M.events),
            rep(1, P.events), rep(0, P.patients-P.events))
z = data.frame(status, outcome)

z.glm = glm(outcome ~ status, family=binomial, data=z)
summary(z.glm)
##             Estimate Std. Error
## (Intercept)  -2.2278     0.1276
## statusM      -0.3948     0.1966

## Log odds ratio       = -0.3948
## Standard uncertainty =  0.1966

## Odds of hospitalization or death among patients on placebo
oddsP = exp(-2.2278)
## Probability of hospitalization or death among patients on placebo
probP = oddsP/(1+oddsP) ## 0.0973

## Odds Ratio
oddsRatioMP = exp(-0.3948)

## Odds of hospitalization or death among patients on molnupiravir
oddsM = oddsRatioMP * oddsP
## Probability of hospitalization or death among patients on molnupiravir
probM = oddsM/(1+oddsM) ## 0.0677

######################################################################
##
## EXAMPLE 4H (CADMIUM CALIBRATION STANDARD)
##
######################################################################

## Example A1 of the EURACHEM/CITAC Guide CG-4: preparation of a
## calibration standard for the determination of cadmium using 
## atomic absorption spectroscopy

m.x = 100.28     ## mg
P.x = 0.9999     ## g/g
V.x = 100.0      ## mL

m.u   = 0.05     ## mg
P.u   = 0.000058 ## g/g
V.u   = 0.07     ## mL

K = 1e7

## Gaussian distribution for m truncated at 0 mg
require(truncnorm)
m = rtruncnorm(K, a=0, b=+Inf, mean=m.x, sd=m.u)

## Beta distribution for P
alpha = P.x * ( P.x * (1 - P.x)/P.u^2 - 1) 
beta = (1 - P.x) * alpha / P.x
curve(dbeta(x, alpha, beta), n=2^12, from=0.9995, to=1)
P = rbeta(K, shape1=alpha, shape2=beta)

## Triangular distribution for V
require(triangle)
a = V.x - V.u*sqrt(6) ## 99.82854
b = V.x + V.u*sqrt(6) ## 100.1715
c = (a + b)/2
curve(dtriangle(x, a=a, b=b, c=c), from=a, to=b)
V = rtriangle(n=K, a=a, b=b, c=c)

## Calculate concentration of cadmium
cCd = m*P/V

## Always examine the results graphically before summarizing them, to
## ensure that the chosen summaries are meaningful and fit for purpose

plot(density(cCd), bty="n", lwd=3, col="DodgerBlue", main="",
     xlab=expression(italic(c)(plain(Cd))~{}/{}~(plain(mg)/plain(mL))),
     ylab="Prob. Density")
curve(dnorm(x, mean=mean(cCd), sd=sd(cCd)), add=TRUE, n=2^12, lwd=1, col="Red")
## NOTE: Probability distribution of cCd deviates significantly from
## Gaussian shape

## Summarize results (mean, std. uncertainty, 95 % expanded uncertainty)
round(c(mean(cCd), sd(cCd)), 4) ## 1.0027 0.0009 ## mg/mL
round(diff(quantile(cCd, probs=c(0.025, 0.975)))/2, 4) ## 0.0017 ## mg/mL

######################################################################
##
## EXAMPLE 4I (ARGENTOMETRIC TITRATION)
##
######################################################################

V0.x = 100/1000  ## L     (100 mL of NaCl solution)
Ca.x = 100/1000  ## mol/L (100 mmol/L of AgNO3 -- Titrant)
Va.x = 10/1000   ## L     (10 mL of AgNO3 solution consumed)
Ck.x = 1.75/1000 ## mol/L (1.75 mmol/L of K2CrO4 -- Indicator)

## V0 (volume of NaCl solution ) measured with a Class A volumetric
## pipette whose measurement error is modeled using a symmetric
## triangular distribution of half-width 0.08 mL 

K = 1e7
require(triangle)
V0 = V0.x + rtriangle(K, a=-0.08/1000, b=+0.08/1000, c=0)

## Va (volume of AgNO3 consumed to reach endpoint) measured using an
## electronic burette with 0.01 mL resolution, whose measurement error is
## modeled using a symmetric triangular distribution of half-width 0.01 mL.
## Va is also affected by the error in ascertaining that the endpoint
## was reached, which is modeled using a Gaussian distribution
## centered at 0 mL and with standard deviation 0.025 mL
Va = Va.x + rtriangle(K, a=-0.01/1000, b=+0.01/1000, c=0) +
    rnorm(K, mean=0, sd=0.025/1000)

## The uncertainty associated with the amount concentration of AgNO3
## in the titrant is assumed to be negligible
Ca = Ca.x

C0 = Ca * Va/V0
C0.x = mean(C0)
C0.u = sd(C0)
c(C0.x, C0.u) ## 0.010000 0.000026 ## mol/L

## Solubility products of AgCl (Ksp1) and Ag2CrO4 (Ksp2)
Ksp1 = 10^(-9.75)
Ksp2 = 10^(-11.9)
## VE is the true value of Va (added volume of titrant) that takes
## into account these solubility products
f = function (VE, Ca, Ksp1, Ksp2, V0, C0, Ck)
{
    (Ksp1 / sqrt(Ksp2*(V0+VE)/(Ck*V0))) - sqrt(Ksp2*(V0+VE)/(Ck*V0)) -
        (C0*V0 - Ca*VE) / (V0 + VE)
}

VE = uniroot(f, interval=c(0, 0.05), tol=sqrt(.Machine$double.eps),
             C=Ca.x, Ksp1=Ksp1, Ksp2=Ksp2, V0=V0.x, C0=C0.x, Ck=Ck.x)$root
## 0.010024 ## L

## Naive estimate of amount concentration of NaCl
C0.x = Va.x * Ca.x / V0.x

## Estimate of amount concentration of NaCl that takes into account
## the solubility products of AgCl and of Ag2CrO4
C0.STAR = VE * Ca.x / V0.x

## Relative error of naive estimate of amount concentration of NaCl
(C0.STAR - C0.x) / C0.x ## 0.0024

######################################################################
##
## EXAMPLE 4J (MICHAELIS-MENTEN)
##
######################################################################

## The beta-galactosidase enzyme catalizes the hydrolization of
## o-nitrophenyl-beta-galactoside (ONPG) to produce galactose and
## o-nitrophenol at a rate, v, that depends on the amount
## concentration of ONPG, c, as described by the Michaelis-Menten
## equation [51, 52]

c = c(0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 8, 20, 30)             ## mmol/L
v = c(3, 5.2, 14.4, 30.3, 49, 86.2, 112.6, 136.2, 170, 177.7) ## umol/L/min

v.nls = nls(v ~ vMax*c/(KM + c), start=list(vMax=195, KM=3),
            nls.control(maxiter=5000, tol=1e-07, minFactor= 1/2^12),
            data=data.frame(c=c, v=v))

summary(v.nls)
##      Estimate Std. Error
## vMax 195.3556     3.2342  ## umol/L/min
## KM     3.2662     0.1837  ## mmol/L
## Residual standard error: 3.085 on 8 degrees of freedom ## umol/L/min

vMaxHAT = coef(v.nls)[1]
KMHAT = coef(v.nls)[2] 
sigmaHAT = summary(v.nls)$sigma ## 3.085365

## Monte Carlo evaluation of the uncertainties associated with the
## estimates of vMax and KM -- essentially reproducing, in an
## independent manner, the values listed under "Std. Error" above

n = length(c)
K = 10000
theta = array(rep(NA, K*2), dim=c(K,2))
for (k in 1:K)
{
    if (k %% 1000 == 0) {cat(k, "of", K, "\n")}
    vB = vMaxHAT*c/(KMHAT+c) + sigmaHAT*rt(n, df=8)/sqrt(8/(8-2))

    vB.nls = try(nls(v ~ vMax * c/(KM+c),
                     start=list(vMax=195, KM=3),
                     nls.control(maxiter=5000, tol=1e-07, minFactor= 1/2^12),
                     data=data.frame(c=c, v=vB)))
    if (inherits(vB.nls, "try-error")) {next
    } else {theta[k,] = coef(vB.nls)}
}

iNA = apply(theta, 1, function(x){is.na(x)})
if (sum(iNA) > 0) {theta = theta[!iNA,]}

## Always examine the results graphically before summarizing them, to
## ensure that the chosen summaries are meaningful and fit for purpose

par(mfrow=c(1,2))
plot(density(theta[,1], n=2^12, adjust=1.5), lwd=2, col="DodgerBlue",
     main="", bty="n",
     xlab=expression(italic(v)[plain(MAX)]~{}/{}~##
                         ((mu*plain(mol)/plain(L))/plain(min))),
     ylab="Prob. Density")
plot(density(theta[,2], n=2^12, adjust=1.5), lwd=2, col="Tomato", 
main="", bty="n",
     xlab=expression(italic(K)[plain(M)]~{}/{}~(plain(mmol)/plain(L))),
     ylab="Prob. Density")
par(mfrow=c(1,1))

apply(theta, 2, sd) ##  3.2  0.18  ## umol/L/min, mmol/L
## Compare with standard errors listed above following "summary(v.nls)"

######################################################################
##
## EXAMPLE 4K (GREEN SEAMOUNT BASALTS)
##
######################################################################

## Evaluating the uncertainty of the median of replicated
## determinations of the mass fraction of titanium dioxide in samples
## of volcanic glass from the Green Seamount, in the Pacific Ocean,
## about 435 km to the west of Puerto Vallarta, Mexico

w = c(15.2, 14.2, 15.1, 15, 12.4, 12.2, 12.7, 14.9, 14, 14, 14.4,
      20.1, 20.5, 11.2, 13.7, 14.1, 12.7, 12.8, 17.9, 11.7, 13.1,
      13.5) ## / (mg/g)

## Sample is not likely to be from a Gaussian distribution
shapiro.test(w) ## p-value = 0.003052
## Sample is consistent with a symmetric distribution
require(symmetry)
symmetry_test(w, stat="MGG", B=25000)$p.value ## 0.32

## Adopt measurement error model of Example 3D, but with Laplace
## measurement errors. The median is the maximum likelihood estimator
## of the mean of the Laplace distribution

median(w) ## 14.0 mg/g

## Evaluation of the uncertainty associated with the median via
## nonparametric bootstrap resampling requires a "large" set of
## replicates determinations -- here we have only 30 but, as we shall
## see, it still suffices to obtain an acceptably accurate evaluation

m = length(w)
muB = replicate(1e5, { median(sample(w, size=m, replace=TRUE)) })
round(sd(muB), 2)                       ## 0.41 mg/g
round(quantile(muB, c(0.025,0.975)), 1) ## 12.9, 14.7 mg/g

## To validate the bootstrap results, we also obtain another,
## independent evaluation of the uncertainty associated with the
## median, by computing a 95 % coverage interval for the true median
## via inversion of the sign test

## The sign test of the hypothesis that the median is a specified
## value, say 15 mg/g, is based on the number of positive signs of
## the differences {w[i] - 15 mg/g}.

## If 15 mg/g indeed were the true median, then each determination
## would be equally likely to be above or below 15 mg/g, and the
## number of positive difference would have a binomial distribution
## based on 22 trials with 50 % probability of "success" in each
## trial.

## Since 5 determinations are greater than 15 mg/g, the probability
## of observing these many or fewer positive differences is pbinom(5,
## size=22, prob=0.5) = 0.0085. However, we are testing whether there
## are either too few positive differences or too few negative
## differences. Therefore the p-value of the test is 2 * 0.0085 =
## 0.017.

## Now, suppose that we choose to reject the hypothesis of the true
## median being M when the p-value of the test is 0.05 or smaller. The
## interval of values of M that the test does not reject is a 95 %
## coverage interval for the true median: determining this interval is
## called "inverting the test," and can be done using R function
## "SIGN.test" defined in R package BSDA

require(BSDA)
round(SIGN.test(w, conf.level=0.95)$conf.int, 1)
##  12.8, 14.9 mg/g

## In many applications, the standard uncertainty of the median is
## evaluated as mad(w)*sqrt(pi/2)/sqrt(length(w)) = 0.46 mg/g in this
## case, which happens to be 0.46/0.41 = 1.12 times too large (but
## could just as well have been too small).

## This widely misused formula is valid only when n is large and the
## data are a sample from a Gaussian distribution. But if the data are
## a sample from a Gaussian distribution, then the best estimate is
## the average of the determinations, not the median, whose
## corresponding standard uncertainty is s/sqrt(n).

######################################################################
##
## EXAMPLE 4N (MOLAR MASS OF IRIDIUM -- MARKOV CHAIN MONTE CARLO)
##
######################################################################

## We recommend that, instead of implementing the statistical
## procedures mentioned in this Brief Guide themselves, users should
## rely on tried-and-true implementations available in trustworthy,
## preexisting software, which is available from many sources, for R
## and for most other environments for statistical modeling and data
## analysis. This recommendation is particularly emphatic for any and
## all procedures that involve Markov Chain Monte Carlo sampling

Ir.M  = c(C=192.2168, W=192.2166, Z=192.2176) ## g/mol
Ir.uM = c(C=  0.0008, W=  0.0003, Z=  0.0002) ## g/mol

require(MCMCpack)

require(extraDistr)
median(rhcauchy(1e6, sigma=0.0004)) ## 0.0003997044

## User defined function "lognum" computes the logarithm of the
## numerator of Bayes's Rule, which comprises the product of the prior
## probability density function of the parameters and of the
## likelihood function

lognum = function (theta, Ir.M, Ir.uM)
{
    M = theta[1]; tau = theta[2]
    ## Prior distribution for the molar mass of ytterbium
    s = dnorm(M, mean=192.217, sd=0.0015, log=TRUE)
    ## Prior distribution for the dark uncertainty, tau, has median
    ## 0.0004 g/mol
    s = s + dhcauchy(tau, sigma=0.0004, log=TRUE)
    ## Gaussian likelihood
    s = s + dnorm(Ir.M["C"], mean=M,
                  sd=sqrt(tau^2 + Ir.uM["C"]^2), log=TRUE)
    s = s + dnorm(Ir.M["W"], mean=M,
                  sd=sqrt(tau^2 + Ir.uM["W"]^2), log=TRUE)
    s = s + dnorm(Ir.M["Z"], mean=M,
                  sd=sqrt(tau^2 + Ir.uM["Z"]^2), log=TRUE)
    return(s)
}

Ir.mcmc =
    MCMCmetrop1R(fun=lognum,

                 ## Initial values where the MCMC chain will start,
                 ## drawn randomly from the corresponding prior
                 ## distributions 
                 theta.init=c(rnorm(1, mean=192.217, sd=0.0015),
                              rhcauchy(1, sigma=0.0004)),

                 ## The number of burnin steps should be about one
                 ## half of the number of actual MCMC sampling
                 ## steps. The "thin" input specifies that only 1 in
                 ## every 25 samples generated by MCMC will be kept --
                 ## to reduce correlations between consecutive values
                 ## of the parameters
                 burnin=250000, mcmc=250000, thin=25,

                 ## The parameter "tune" (which can be a scalar or a
                 ## vector with as many components as there are
                 ## parameters) should be chosen so as to achieve
                 ## Metropolis acceptance rate of about 0.3 -- the
                 ## function prints this acceptance rate upon
                 ## termination, and "tune" should be adjusted as
                 ## needed, and the run repeated, until a satisfactory
                 ## acceptance rate is achieved
                 tune=c(4, 1.5),

                 ## If there are box constraints on the parameters, as
                 ## there are here for tau, which must be
                 ## non-negative, then the optimization method needs
                 ## to be "L-BFGS-B"
                 optim.method="L-BFGS-B",
                 ## The following inputs define the lower and upper
                 ## bounds for the "legal" values of the parameters
                 optim.lower=c(-Inf, 0), optim.upper=c(Inf, Inf),
                 optim.control=list(fnscale=-1, trace=0, REPORT=10,
                                    maxit=5000, ndeps=c(0.5e-4,0.5e-5)), 

                 ## Values needed to compute the numerator of Bayes's Rule
                 Ir.M=Ir.M, Ir.uM=Ir.uM)

dimnames(Ir.mcmc)[[2]] = c("M", "tau")
dim(Ir.mcmc) ## 10000     2

summary(Ir.mcmc)

## Iterations = 250001:499976
## Thinning interval = 25 
## Number of chains = 1 
## Sample size per chain = 10000 
## 
## 1. Empirical mean and standard deviation for each variable,
##    plus standard error of the mean:
## 
##          Mean        SD  Naive SE Time-series SE
## M   1.922e+02 0.0004173 4.173e-06      4.173e-06
## tau 5.598e-04 0.0004138 4.138e-06      6.654e-06
## 
## 2. Quantiles for each variable:
## 
##          2.5%       25%       50%       75%     97.5%
## M   1.922e+02 1.922e+02 1.922e+02 1.922e+02 1.922e+02
## tau 6.064e-05 2.979e-04 4.638e-04 7.060e-04 1.652e-03

M.TILDE = Ir.mcmc[,"M"]
tau.TILDE = Ir.mcmc[,"tau"]

par(mfrow=c(1,2))
plot(density(M.TILDE, adj=1.5, from=192.2145, to=192.2195), main="",
     xlab=expression(italic(M)[plain(Ir)]~{}/{}~(plain(g)/plain(mol))),
     ylab="Posterior Prob. Density", lwd=2, col="DodgerBlue")
plot(density(tau.TILDE, adj=1.5, from=0, to=0.002), main="",
     xlab=expression(tau~{}/{}~(plain(g)/plain(mol))),
     ylab="Posterior Prob. Density", lwd=2, col="Tomato")
par(mfrow=c(1,1))

round(c(M=mean(M.TILDE), "u(M)"=sd(M.TILDE)), 4)
##        M     u(M) 
## 192.2171   0.0004  ## g/mol

round(median(tau.TILDE), 4) ## 0.0005 ## g/mol

######################################################################
##
## EXAMPLE 4O (NITRITES IN SEAWATER)
##
######################################################################

## Four spectrophotometric determinations of the mass fraction of
## nitrites in a seawater sample obtained under repeatability
## conditions of measurement using Griess’s method

w = c(0.1514, 0.1523, 0.1531, 0.1545) ## mg/kg

## Relative measurement uncertainty could be as low as 0.33 % or as
## high as 3 %, lying somewhere between these two values with 90 %
## probability 

## Calibration of a gamma prior distribution for sigma that captures
## this vague prior knowledge about its true value

q = c(0.0033, 0.03)*median(w)
p = c(0.05,   0.95)
g = function (theta, q, p)
{
    if (any(theta < 0)) { return(Inf)
    } else {
        alpha = theta[1]; beta = theta[2]
        return(sum((qgamma(p, shape=alpha, rate=beta) - q)^2))
    }
}
optim(par=c(2, 500), fn=g, gr=NULL, q=q, p=p, method="Nelder-Mead")$par
## Shape = 2.619, Rate = 1248 kg/mg

## Bayesian model implemented using the probabilistic programming
## language Stan [B. Carpenter et al., 2017, "Stan: A Probabilistic
## Programming Language." Journal of Statistical Software 76(1):
## 1-32. DOI 10.18637/jss.v076.i01]

require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write=TRUE, javascript=FALSE)

model =
 "data
    { 
      // The determinations are stored in a vector of length 4
      real w[4]; 
    }
  parameters
    { 
      // True mean mass fraction of nitrites, which must be non-negative
      real<lower=0> omega;

      // True standard deviation, which must be non-negative,
      // of measurement errors using Griess's method
      real<lower=0> sigma; 
    }
  model
    { 
      // Prior distribution chosen for the true mean mass fraction 
      // of nitrites in the sample is truncated Gaussian with mean 0 mg/kg
      // and standard deviation 1 mg/kg (both before truncation), 
      // truncated at 0 mg/kg because omega must be non-negative
      omega ~ normal(0, 1); 

      // Prior distribution chosen for std. deviation of measurement errors
      // that takes into account knowledge of the method's repeatability
      sigma ~ gamma(2.619, 1248);  

      // Likelihood function corresponding to the assumption that the
      // measured values are like a sample from a Gaussian distribution
      // whose mean is the true value of the mass fraction and whose
      // standard deviation is the true value of the repeatability
      w ~ normal(omega, sigma); 
    }"

fit = stan(model_code = model, data = list(w=w), 

           ## It is recommended to use the initial half of the steps
           ## (the total number of steps requested here is 500000)
           ## to warm-up and auto-tune the samplers
           warmup=250000, iter=500000, 

           ## Run four Markov Chain samplers independently of one
           ## another, in parallel, using four CPU cores (if the CPU
           ## being used does not have 4 cores, then "cores=4,"
           ## should be removed from the next line
           chains=4, cores=4,

           ## Keep only every 25th set of parameter values sampled, to
           ## reduce (serial) correlations between them
           thin=25)

print(fit, digits=4)

## The values of Rhat should be very close to 1: otherwise, the model
## should be fitted using a larger number of iterations. The values of
## n_eff are effective sample sizes that the following summaries are
## based on

##          mean      sd    2.5%     50%   97.5%  n_eff     Rhat
## omega  0.1528  0.0010  0.1509  0.1528  0.1548  40478  0.99997
## sigma  0.0018  0.0007  0.0008  0.0016  0.0036  39757  0.99999

fit.post = extract(fit)

## The (posterior) standard uncertainty associated with the Bayesian
## estimate of the true mass fraction is 1.45 times larger than the
## conventional Type A evaluation computed according to the GUM

sd(fit.post$omega) / (sd(w)/sqrt(length(w))) ## 1.45

######################################################################
##
## EXAMPLE 4P (MOLNUPIRAVIR -- SKEPTICAL RECONSIDERATION)
##
######################################################################

## Mean and standard deviation of the skeptical prior for the log odds ratio
muP = 0; sigmaP = 0.175
## Log odds ratio observed in the MOVe-OUT trial reported by Bernal et
## al. (2022), and associated standard uncertainty (from EXAMPLE 4G)
logOR = -0.395; logOR.u = 0.197

## Posterior mean and standard deviation
muQ = (logOR/logOR.u^2 + muP/sigmaP^2) / (1/logOR.u^2 + 1/sigmaP^2)
sigmaQ = sqrt(1/(1/logOR.u^2 + 1/sigmaP^2))

c(muQ, sigmaQ) ## -0.174  0.131
## Posterior probability that molnupiravir is no more efficacious than placebo
1-pnorm(0, mean=muQ, sd=sigmaQ) ## 0.091

######################################################################
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
######################################################################

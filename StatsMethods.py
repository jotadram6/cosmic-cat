import numpy as np
from scipy import optimize
import math

#Implementation following 1007.1727

def LikelihoodFunction(mu,theta,n,s,b,m,u,NoSys):
    """
    Functio that calculates the Likelihood for a given:
    mu -> Signal strength
    theta -> Vector of systematic uncertainties (theta_s,theta_b,b_tot) in general, only implemented as a single paramter controling b normalization
    n -> Vector (histogram) of n expected counts
    s -> Vector (histogram) of expected signal events, must be the same dimension as n
    b -> Vector (histogram) of expected background events, must be the same dimension as n
    m -> Vector (histogram) of m expected counts in a control region
    u -> Vector (histogram) of modelling events in the control region modeling systematic uncertainties may be functions
    """
    if len(n)!=len(s) or len(n)!=len(b):
        print("Size of n, s and b must be the same!")
        return 0
    if len(m)!=len(u):
        print("Size of m and u must be the same!")
        return 0
    L1=1.0; L2=1.0
    for i in range(len(n)):
        L1*=((((mu*s[i])+b[i])**(n[i]))/(math.factorial(math.floor(n[i]))))*(np.exp(-1*((mu*s[i])+b[i])))

    if not NoSys:
        for i in range(len(m)):
            L2*=((u[i](theta,b[i])**m[i])/(math.factorial(math.floor(m[i]))))*(np.exp(-1*u[i](theta,b[i])))
        return L1*L2
    else:
        return L1

def SimplifiedLikelihoodFunction(mu,n,s,b):
    """
    Functio that calculates the Likelihood for a given:
    mu -> Signal strength
    theta -> Vector of systematic uncertainties (theta_s,theta_b,b_tot) in general, only implemented as a single paramter controling b normali
zation
    n -> Vector (histogram) of n expected counts
    s -> Vector (histogram) of expected signal events, must be the same dimension as n
    b -> Vector (histogram) of expected background events, must be the same dimension as n
    """
    if len(n)!=len(s) or len(n)!=len(b):
        print("Size of n, s and b must be the same!")
        return 0
    L1=1.0
    for i in range(len(n)):
        L1*=((((mu*s[i])+b[i])**(n[i]))/(math.factorial(math.floor(n[i]))))*(np.exp(-1*((mu*s[i])+b[i])))
    return L1

def OrigSimplifiedLogL(mu,n,s,b):
    return math.log(SimplifiedLikelihoodFunction(mu,n,s,b))

################################################################
################################################################
################################################################

def SimplifiedLogL(mu,n,s,b):
    """
    Function that calculates the Likelihood for a given:
    mu -> Signal strength
    theta -> Vector of systematic uncertainties (theta_s,theta_b,b_tot) in general, only implemented as a single paramter controling b normali
zation
    n -> Vector (histogram) of n expected counts
    s -> Vector (histogram) of expected signal events, must be the same dimension as n
    b -> Vector (histogram) of expected background events, must be the same dimension as n
    """
    if len(n)!=len(s) or len(n)!=len(b):
        print("Size of n, s and b must be the same!")
        return 0
    L1=0.0
    for i in range(len(n)):
        L1+=math.log(((((mu*s[i])+b[i])**(n[i]))/(math.factorial(math.floor(n[i]))))*(np.exp(-1*((mu*s[i])+b[i]))))
    return L1

def SimplifiedProfileLogL(mu,muhat,n,s,b):
    return SimplifiedLogL(mu,n,s,b)-SimplifiedLogL(muhat,n,s,b)

def SimplifiedTestStatistics(mu,muhat,n,s,b):
    return -2*SimplifiedProfileLogL(mu,muhat,n,s,b)

def ShiftSignSimplifiedLogL(mu,n,s,b):
    #return -1*math.log(SimplifiedLikelihoodFunction(mu,n,s,b))
    return -1*SimplifiedLogL(mu,n,s,b)

def MaxSimpLogL(mu,n,s,b,minm,maxm):
    return optimize.minimize(ShiftSignSimplifiedLogL,[mu],args=(n,s,b),bounds=[(minm,maxm)])

###########################################################
###########################################################
###########################################################

#Based on section 3.6 from 1007.1727
def UpperLimitsTestStats(muhat,mu,sigma):
    if mu>muhat:
        return ((mu-muhat)**2)/sigma**2
    else:
        return 0

def UperLimitsSignificance(muhat,mu,sigma):
    return np.sqrt(UpperLimitsTestStats(muhat,mu,sigma))

def UppeLimitsMu(muhat,mu,sigma):
    return muhat+sigma*UperLimitsSignificance(muhat,mu,sigma)

def FindUpperLimit(muhat,sigma,alpha=0.05,Smin=0.0,Smax=None,Steps=None):
    """
    Upper limits on signal strength
    Based on the looking for mu thath gives the test statistics which is closer to alpha
    """
    if Smax==None: Smax=muhat*10
    if Steps==None: Steps=int((Smax-Smin)/0.1)

    UpMuLimit=[]
    UpMuScan=np.linspace(Smin,Smax,Steps)
    
    for i in UpMuScan:
        UpMuLimit.append(UpperLimitsTestStats(muhat,i,sigma))
    #print(UpMuLimit)
    
    UpperLimitI=0
    for i in range(len(UpMuScan)):
        if UpMuLimit[i]>alpha:
            UpperLimitI=i-1
            break
    #print(i-1,UpMuLimit[i-1],UpMuScan[UpperLimitI])
    
    return UppeLimitsMu(muhat,UpMuScan[UpperLimitI],sigma)


"""
Usage example:

import StatsMethods
n_3b=[50.,20.,10.]                                                                         
b_3b=[48.,17.,9.]                                                  
s_3b=[5.,4.,2.]
Mymu=1.5
StatsMethods.SimplifiedLikelihoodFunction(Mymu,n_3b,s_3b,b_3b)
StatsMethods.SimplifiedLogL(Mymu,n_3b,s_3b,b_3b)
Maximum=StatsMethods.MaxSimpLogL(Mymu,n_3b,s_3b,b_3b,0.01,None)
StatsMethods.FindUpperLimit(Maximum.x[0],2.5)
"""

###############################################################################################################
###############################################################################################################
###############################################################################################################

#Implementation following 1708.01007
#For unbinned data

def LogLikelihoodParamEstimResonance(x,pdf,M0,Gamma):
    """
    Likelihood function to obtain parameter values from a 
    pdf with two parameters, especially a resonance: Mass and wwidth
    PDF function must be of the form pdf(data,var1,var2)
    """
    LLT=0.0;
    for i in range(len(x)):
        LLT+=math.log(pdf(x[i],M0,Gamma))
    return LLT

#def ShiftSignLogLParamEstimResonance(x,pdf,M0,Gamma):
#    return -1*LogLikelihoodParamEstimResonance(x,pdf,M0,Gamma)

def ShiftSignLogLParamEstimResonance(R,x,pdf):
    return -1*LogLikelihoodParamEstimResonance(x,pdf,R[0],R[1])

def MaxLogLParamEstimResonance(x,pdf,M0,Gamma,M0min,M0max,Gammamin,Gammamax):
    return optimize.minimize(ShiftSignLogLParamEstimResonance,(M0,Gamma),args=(x,pdf),bounds=[(M0min,M0max),(Gammamin,Gammamax)])

"""
Usage example:

import StatsMethods
x=[10.,11.,9.,5.,15.]                                          
def GaussianPDF(xi,mu,sigma):
    return (1./(2*math.pi*sigma*sigma))*(math.exp((-1*((xi-mu)**2))/(2*sigma*sigma)))
StatsMethods.LogLikelihoodParamEstimResonance(x,GaussianPDF,10.,1.0)
StatsMethods.ShiftSignLogLParamEstimResonance(x,GaussianPDF,10.,1.0)
StatsMethods.MaxLogLParamEstimResonance(x,GaussianPDF,10.,1.0,0.0,20.0,1.0,5.0)
"""

##########################################################################################################
##########################################################################################################
##########################################################################################################

#Implementation of parameter estimation following Maximum Likelihood approach discussed on section 6.10
#of Statitical Data Analysis from Glen Cowan

def ParamsLogLikelihood(n_observed,b_expected,s_expected):
    """
    All three histograms must have same number of bins
    This function will calculate the log likelihood for parameter estimation
    describing n=s+b
    """
    if len(n_observed)!=len(b_expected) or len(n_observed)!=len(s_expected):
        print("All three histograms must have same number of bins")
        return

    logL=(-1*(sum(b_expected)+sum(s_expected)))
    for i in range(len(n_observed)):
        logL+=n_observed[i]*math.log(b_expected[i]+s_expected[i])

    return logL

def FindMaxParam(n_observed,b_expected,SignalSet):
    """
    This function will run ParamsLogLikelihood(n_observed,b_expected,s_expected) over all given signal histograms
    and return the one which maximizes the LogLikelihood
    """
    LogLVector=[]
    for i in range(len(SignalSet)):
        LogLVector.append(ParamsLogLikelihood(n_observed,b_expected,SignalSet[i]))

    return max(LogLVector), LogLVector.index(max(LogLVector)), SignalSet[LogLVector.index(max(LogLVector))]

######################################################################################################
######################################################################################################
######################################################################################################

#Implementation of chi2 method and p-value
#As discussed in chapter 7 of Statitical Data Analysis from Glen Cowan

from scipy.stats import chi2 as scipy_chi2

def Chi2pValue(chi2value,NDF):
    """
    chi2value: chi2 value obtained from the test statistics
    NDF: Number of degrees of freedom used to calculate chi2value
    """
    return scipy_chi2.sf(chi2value, NDF, loc=0, scale=1)

def BinnedChi2(Obs,Exp,Normalized=False):
    """
    Obs: 1-d array of observed events bins
    Exp: 1-d array of expected events bins
    If Normalized flag is set to True, the function returns chi2/ndf, 
    otherwise it only returns chi2.
    """
    Chi2=np.sum(((Obs-Exp)**2)/Exp)
    if Normalized:
        return Chi2/len(Exp)
    else:
        return Chi2


"""
Example:
b=[2,2,2,2,2,2,2]
n=[3,3,3,3,4,6,4]
s=[1,1,1,1,2,4,2]
In [10]: StatsMethods.ParamsLogLikelihood(n,b,s)
Out[10]: 9.024259168344772

s1=[0,1,0,1,2,1,0]
In [12]: StatsMethods.ParamsLogLikelihood(n,b,s1)
Out[12]: 6.65999671409633

Signals=random.randint(0,7,(20,7))
In [21]: StatsMethods.FindMaxParam(n,b,Signals)
Out[21]: (8.5769807206118909, 15, array([0, 1, 1, 2, 2, 3, 2]))

BinnedChi2(np.array(n),np.array(b)+0.2*np.array(s),Normalized=True)
Chi2pValue(BinnedChi2(np.array(n),np.array(b)),7)
"""
def MuChi2(mu,n,s,b):
    return BinnedChi2(n,(mu*s)+b)

def MinChi2(mu,n,s,b,minm,maxm):
    return optimize.minimize(MuChi2,[mu],args=(n,s,b),bounds=[(minm,maxm)])

def FindUpperLimitChi2(n,s,b,alpha=0.05,Mumin=0.0,Mumax=10.0,Steps=1000):
    for mui in np.linspace(Mumax,Mumin,Steps):
        chi2v=MuChi2(mui,n,s,b)
        pvalue=Chi2pValue(chi2v,len(n))
        #print(mui,pvalue)
        if pvalue>alpha:
            return mui-((Mumax-Mumin)/Steps), pvalue


######################################################################################################
######################################################################################################
######################################################################################################

#Computational approach without factorials

def CompSimpLogL(mu,n,s,b):
    """
    Function that calculates the Likelihood for a given:
    mu -> Signal strength
    theta -> Vector of systematic uncertainties (theta_s,theta_b,b_tot) in general, only implemented as a single paramter controling b normali
zation
    n -> Vector (histogram) of n expected counts
    s -> Vector (histogram) of expected signal events, must be the same dimension as n
    b -> Vector (histogram) of expected background events, must be the same dimension as n
    """
    if len(n)!=len(s) or len(n)!=len(b):
        print("Size of n, s and b must be the same!")
        return 0
    L1=0.0
    for i in range(len(n)):
        theory=(mu*s[i])+b[i]
        L1+=(n[i]*math.log(theory))-theory
    return L1

def ProfileSimpCompLogL(mu,muhat,n,s,b):
    return CompSimpLogL(mu,n,s,b)-CompSimpLogL(muhat,n,s,b)

def SimpCompTestStatistics(mu,muhat,n,s,b):
    return -2*ProfileSimpCompLogL(mu,muhat,n,s,b)

def ShiftSignCompSimpLogL(mu,n,s,b):
    return -1*CompSimpLogL(mu,n,s,b)

def MaxCompSimpLogL(mu,n,s,b,minm,maxm):
    return optimize.minimize(ShiftSignCompSimpLogL,[mu],args=(n,s,b),bounds=[(minm,maxm)])

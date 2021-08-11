# BestSurvCutPoints
### Finding the optimal number and locations of cutpoints for survival analysis.


The origional publication and the corresponding scripts were retrived based on Attribution 4.0 International (CC BY 4.0) from

Chang, C., Hsieh, M. K., Chang, W. Y., Chiang, A. J., & Chen, J. (2017). Determining the optimal number and location of cutoff points with application to data of cervical cancer. PloS one, 12(4), e0176231.(https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0176231)

The source code was deposted on https://osf.io/ef7na/ by the authors and a step-by-step user manual for using the script can be found in the URLs below, https://github.com/paytonyau/BestSurvCutPoints/blob/main/User_Manual_mod.pdf


The authors also deposited additional scripts using COX on (http://www.math.nsysu.edu.tw/~cchang/) with a simulated dataset and an example

This program is open for everyone to edit and I am preparing to wrap the functions into a R package, 

(1) `findcutnum` - calculate akaike information criterion (AIC) values for both survival and binary outcomes

(2) `findcut`- use contingency tables (X^2) approach to make a decision on how many cutoff numbers are needed before it can be employed

(3) `findcutCox` - estimating the optimal locations of cutpoints

(4) `findnumCox` -  estimating the optimal number of cutpoints



# 20180328 Thoughts on analyzing the EPF data

We've talked about analyzing the EPF data with treatment as an explantory factor and also with pH as an explanatory variable instead.

The other question is how to treat external-corrected-EPF. I chatted with Tarik about this, because I was worried about introducing 
spurious correlation (see Biostats notes from Correlation chapter) when we subtract the external pH from the EPF, and then analyze 
it with external pH as an explanatory variable. 

We can analyze external-corrected-EPF as a response variable with treatment as an explanatory variable,  
because the SS calculation is not affected by the correction.

We have to be careful when we analyze external-corrected-EPF as a response variable and external pH as an explanatory variable, 
because we have to worry about spurious correlation introduced by the correction. 
There are Monte-carlo methods for testing against the null hypothesis, which in this case would be the spurious correlation expected 
by random (instead of 0).
(see Biostats notes from Correlation chapter).

fastxtile: just like Stata's xtile, but faster
----------------------------------------------

fastxtile is a drop-in replacement for Stata's built-in command xtile. It has the same syntax and produces identical results, but runs substantially faster in large datasets. fastxtile also has a few added features.

### Installation ###

Install fastxtile in Stata from the SSC repository: `ssc install fastxtile`

After installing fastxtile, you can read the documentation by running `help fastxtile`.

### Technical Details ###

fastxtile has been optimized to be more computationally efficient than xtile. It avoids the creation of unnecessary temporary variables (tempvars) and unnecessary sorting of the dataset.  The difference in running time is substantial in large datasets.

fastxtile also supports computing the quantile boundaries using a random sample of the data. This further increases the speed, but generates approximate quantiles due to sampling error.
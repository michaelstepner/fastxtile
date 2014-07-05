/*
fastxtile, version 0.4, 26 June 2012.
Created by Michael Stepner, e-mail at michaelstepner@gmail.com

Drop in replacement for xtile. Runs significantly faster than xtile on large datasets.

- Supports computing the quantile boundaries using a random sample of observations
	- randvar(<varname>) specifies the variable used to define the random sample
		- defaults to the variable "rand" if it exists. Otherwise, if randcut() in (0,1) the sample will be defined by a uniform random variable.
	- randcut(<upper bound>) gives the upper bound on randvar
	- Quantile boundaries are computed using the observations for which randvar<randcut

- Supports either a varname or an ascending list of numbers in cutpoints(), whereas xtile only supports a varname
*/

capture program drop fastxtile
program define fastxtile
	version 11

	* Parse weights, if any
	_parsewt "aweight fweight pweight" `0' 
	local 0  "`s(newcmd)'" /* command minus weight statement */
	local wt "`s(weight)'"  /* contains [weight=exp] or nothing */

	* Extract parameters
	syntax newvarname=/exp [if] [in] [,Nquantiles(integer 2) Cutpoints(string) randvar(varname) randcut(real 1) ALTdef]

	* Mark observations which will be placed in quantiles
	marksample touse, novarlist
	markout `touse' `exp'


	if "`cutpoints'"=="" { /***** NQUANTILES *****/
		if `"`wt'"'!="" & "`altdef'"!="" {
			di as error "altdef option cannot be used with weights"
			exit 198
		}

		* Prepare randvar defaults
		if "`randvar'"=="" {
			capture confirm variable rand
			if _rc==0 {
				local randvar rand
			}
			else {
				if `randcut'<1 & `randcut'>0 {
					tempvar randvar
					gen `randvar'=uniform()
				}
				else if `randcut'!=1 {
					di as error "if you do not specify randvar() and there is no existing variable named 'rand', randcut() must in (0,1)"
					exit 198
				}			
			}
		}


		* Mark observations used to calculate quantile boundaries
		tempvar randsample
		if ("`randvar'"!="") {
			mark `randsample' `wt' if `touse' & `randvar'<`randcut'
		}
		else {
			mark `randsample' `wt' if `touse'
		}

		* Error checks
		qui count if `randsample'
		if (`nquantiles' > r(N) + 1) {
			if ("`randvar'"=="") di as error "nquantiles() must be less than or equal to number of nonmissing observations [`r(N)'] plus one"
			else di as error "nquantiles() must be less than or equal to number of sampled nonmissing observations [`r(N)'] plus one"
			exit 198
		}
		else if (`nquantiles' < 2) {
			di as error "nquantiles() must be greater than or equal to 2"
			exit 198
		}

		* Compute quantile boundaries
		_pctile `exp' if `randsample' `wt', nq(`nquantiles') `altdef'

		local prefix "r(r"
		local suffix ")"
	}
	else { /***** CUTPOINTS *****/
		if "`wt'"!="" | "`randvar'"!="" | "`ALTdef'"!="" | `randcut'!=1 | `nquantiles'!=2 {
			di as error "cutpoints() cannot be used with nquantiles(), altdef, weights, randvar() or randcut()"
			exit 198
		}

		tokenize `cutpoints'
		if (real("`1'")==.) { /* if the first token is a word */
			capture confirm variable `1'
			if ("`2'"!="") | (_rc!=0) {
				di as error "cutpoints() must either consist of one varname or an ascending list of numbers"
				exit 198
			}

			qui count if `1'<.
			if r(N)==0 {
				di in red "cutpoints() all missing"
				exit 2000
			}
			local nquantiles = r(N) + 1

			tempname cutpointvals
			qui tab `1', matrow(cutpointvals)

			local prefix "cutpointvals["
			local suffix ",1]"
		}
		else { /* if the first token is a number */
			tempname r1
			scalar define r1=`1'

			local i=2
			while "``i''"!="" {
				if (real("``i''")==.) {
					di as error "cutpoints() must either consist of one varname or an ascending list of numbers"
					exit 198
				}
				local j=`i'-1
				if (``i''<=``j'') {
					di as error "the specified cutpoints must be strictly ascending from left to right"
					exit 198
				}

				tempname r`i'
				scalar define r`i'=``i''

				local i=`i'+1
			}
			local nquantiles=`i'
			local prefix "r"
			local suffix ""
		}
	}

	local nqm1=`nquantiles'-1

	* Create quantile variable
	qui gen `varlist'=1 if `touse' &  `exp'<=`prefix'1`suffix'
	
	if `nquantiles'>2 {
		tempname i j
		forvalues i = 2/`nqm1' {
			local j=`i'-1
			qui replace `varlist'=`i' if `touse' & `exp'<=`prefix'`i'`suffix' & `exp'>`prefix'`j'`suffix'
		}			
	}
	
	qui replace `varlist'=`nquantiles' if `touse' & `exp'>`prefix'`nqm1'`suffix'

	label var `varlist' "`nquantiles' quantiles of `exp'"

end


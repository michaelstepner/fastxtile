*! version 1.0  4nov2012  Michael Stepner, michaelstepner@gmail.com

program define fastxtile
	version 11

	* Parse weights, if any
	_parsewt "aweight fweight pweight" `0' 
	local 0  "`s(newcmd)'" /* command minus weight statement */
	local wt "`s(weight)'"  /* contains [weight=exp] or nothing */

	* Extract parameters
	syntax newvarname=/exp [if] [in] [,Nquantiles(integer 2) Cutpoints(varname) ALTdef ///
		CUTValues(numlist ascending) randvar(varname) randcut(real 1) version]
		
	if ("`version'"=="version") di "fastxtile, version 1.0, 4 November 2012."

	* Mark observations which will be placed in quantiles
	marksample touse, novarlist
	markout `touse' `exp'


	if "`cutpoints'"=="" & "`cutvalues'"=="" { /***** NQUANTILES *****/
		if `"`wt'"'!="" & "`altdef'"!="" {
			di as error "altdef option cannot be used with weights"
			exit
		}

		* Check if need to gen a temporary uniform random var
		if `randcut'<1 & `randcut'>0 & "`randvar'"=="" {
			tempvar randvar
			gen `randvar'=uniform()
		}
		* Randcut sanity check
		else if `randcut'!=1 {
			di as error "if randcut() is specified without randvar(), a uniform r.v. will be generated and randcut() must be in (0,1)"
			exit
		}

		* Mark observations used to calculate quantile boundaries
		if ("`randvar'"!="") {
			tempvar randsample
			mark `randsample' `wt' if `touse' & `randvar'<=`randcut'
		}
		else {
			local randsample `touse'
		}

		* Error checks
		qui count if `randsample'
		if (`nquantiles' > r(N) + 1) {
			if ("`randvar'"=="") di as error "nquantiles() must be less than or equal to number of nonmissing observations [`r(N)'] plus one"
			else di as error "nquantiles() must be less than or equal to number of sampled nonmissing observations [`r(N)'] plus one"
			exit
		}
		else if (`nquantiles' < 2) {
			di as error "nquantiles() must be greater than or equal to 2"
			exit
		}

		* Compute quantile boundaries
		_pctile `exp' if `randsample' `wt', nq(`nquantiles') `altdef'

		local prefix "r(r"
		local suffix ")"
	}
	else if "`cutpoints'"!="" { /***** CUTPOINTS *****/
	
		* Parameter checks
		if "`cutvalues'"!="" {
			di as error "cannot specify both cutpoints() and cutvalues()"
			exit
		}		
		if "`wt'"!="" | "`randvar'"!="" | "`ALTdef'"!="" | `randcut'!=1 | `nquantiles'!=2 {
			di as error "cutpoints() cannot be used with nquantiles(), altdef, randvar() or randcut() or weights"
			exit
		}


		qui count if `cutpoints'<.
		if r(N)==0 {
			di as error "cutpoints() all missing"
			exit 2000
		}
		local nquantiles = r(N) + 1

		tempname cutpointvals
		qui tab `1', matrow(cutpointvals)

		local prefix "cutpointvals["
		local suffix ",1]"
	}
	else { /***** cutvalues *****/
		if "`wt'"!="" | "`randvar'"!="" | "`ALTdef'"!="" | `randcut'!=1 | `nquantiles'!=2 {
			di as error "cutvalues() cannot be used with nquantiles(), altdef, randvar(), randcut() or weights"
			exit
		}
		
		* parse numlist
		numlist "`cutvalues'"
		local nquantiles=`: word count `r(numlist)''
		tokenize `r(numlist)'
		
		* store in scalars
		forvalues i=1/`nquantiles' {
			tempname r`i'
			scalar define r`i'=``i''
		}
		
		local prefix "r"
		local suffix ""
	
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


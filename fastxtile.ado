*! version 1.11  2sep2013  Michael Stepner, stepner@mit.edu

program define fastxtile, rclass
	version 11

	* Parse weights, if any
	_parsewt "aweight fweight pweight" `0' 
	local 0  "`s(newcmd)'" /* command minus weight statement */
	local wt "`s(weight)'"  /* contains [weight=exp] or nothing */

	* Extract parameters
	syntax newvarname=/exp [if] [in] [,Nquantiles(integer 2) Cutpoints(varname numeric) ALTdef ///
		CUTValues(numlist ascending) randvar(varname numeric) randcut(real 1) randn(integer -1)]

	* Mark observations which will be placed in quantiles
	marksample touse, novarlist
	markout `touse' `exp'
	qui count if `touse'
	local popsize=r(N)

	if "`cutpoints'"=="" & "`cutvalues'"=="" { /***** NQUANTILES *****/
		if `"`wt'"'!="" & "`altdef'"!="" {
			di as error "altdef option cannot be used with weights"
			exit 198
		}
		
		if `randn'!=-1 {
			if `randcut'!=1 {
				di as error "cannot specify both randcut() and randn()"
				exit 198
			}
			else if `randn'<1 {
				di as error "randn() must be a positive integer"
				exit 198
			}
			else if `randn'>`popsize' {
				di as text "randn() is larger than the population. using the full population."
				local randvar=""
			}
			else {
				local randcut=`randn'/`popsize'
				
				if "`randvar'"!="" {
					qui sum `randvar', meanonly
					if r(min)<0 | r(max)>1 {
						di as error "with randn(), the randvar specified must be in [0,1] and ought to be uniformly distributed"
						exit 198
					}
				}
			}
		}

		* Check if need to gen a temporary uniform random var
		if "`randvar'"=="" {
			if (`randcut'<1 & `randcut'>0) { 
				tempvar randvar
				gen `randvar'=runiform()
			}
			* randcut sanity check
			else if `randcut'!=1 {
				di as error "if randcut() is specified without randvar(), a uniform r.v. will be generated and randcut() must be in (0,1)"
				exit 198
			}
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
		local samplesize=r(N)
		if (`nquantiles' > r(N) + 1) {
			if ("`randvar'"=="") di as error "nquantiles() must be less than or equal to the number of observations [`r(N)'] plus one"
			else di as error "nquantiles() must be less than or equal to the number of sampled observations [`r(N)'] plus one"
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
	else if "`cutpoints'"!="" { /***** CUTPOINTS *****/
	
		* Parameter checks
		if "`cutvalues'"!="" {
			di as error "cannot specify both cutpoints() and cutvalues()"
			exit 198
		}		
		if "`wt'"!="" | "`randvar'"!="" | "`ALTdef'"!="" | `randcut'!=1 | `nquantiles'!=2 | `randn'!=-1 {
			di as error "cutpoints() cannot be used with nquantiles(), altdef, randvar(), randcut(), randn() or weights"
			exit 198
		}

		tempname cutvals
		qui tab `cutpoints', matrow(`cutvals')
		
		if r(N)==0 {
			di as error "cutpoints() all missing"
			exit 2000
		}
		else {
			local nquantiles = r(N) + 1

			local prefix "`cutvals'["
			local suffix ",1]"
		}
	}
	else { /***** CUTVALUES *****/
		if "`wt'"!="" | "`randvar'"!="" | "`ALTdef'"!="" | `randcut'!=1 | `nquantiles'!=2 | `randn'!=-1 {
			di as error "cutvalues() cannot be used with nquantiles(), altdef, randvar(), randcut(), randn() or weights"
			exit 198
		}
		
		* parse numlist
		numlist "`cutvalues'"
		local nqm1=wordcount(`"`r(numlist)'"')
		local nquantiles=`nqm1'+1

		tokenize `"`r(numlist)'"'
		
		* store in matrix
		tempname cutvals
		matrix `cutvals'=J(`nqm1',1,.)
		forvalues i=1/`nqm1' {
			matrix `cutvals'[`i',1]=``i''
		}
		
		local prefix "`cutvals'["
		local suffix ",1]"
	
	}

	local nqm1=`nquantiles'-1

	* Create quantile variable
	qui gen `varlist'=1 if `touse' &  `exp'<=`prefix'1`suffix'
	
	if `nquantiles'>2 {
		forvalues i = 2/`nqm1' {
			qui replace `varlist'=`i' if `touse' & `exp'<=`prefix'`i'`suffix' & `exp'>`prefix'`=`i'-1'`suffix'
		}			
	}
	
	qui replace `varlist'=`nquantiles' if `touse' & `exp'>`prefix'`nqm1'`suffix'

	label var `varlist' "`nquantiles' quantiles of `exp'"
	
	* Return values
	if ("`samplesize'"!="") return scalar n = `samplesize'
	else return scalar n = .
	return scalar N = `popsize'
	forvalues i=`nqm1'(-1)1 {
		return scalar r`i' = `prefix'`i'`suffix'
	}

end


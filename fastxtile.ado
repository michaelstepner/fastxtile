*! version 1.22  4jul2014  Michael Stepner, stepner@mit.edu

/* CC0 license information:
To the extent possible under law, the author has dedicated all copyright and related and neighboring rights
to this software to the public domain worldwide. This software is distributed without any warranty.

This code is licensed under the CC0 1.0 Universal license.  The full legal text as well as a
human-readable summary can be accessed at http://creativecommons.org/publicdomain/zero/1.0/
*/

* Why did I include a formal license? Jeff Atwood gives good reasons: http://www.codinghorror.com/blog/2007/04/pick-a-license-any-license.html


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

		* Store quantile boundaries in list
		local maxlist 248
		forvalues i=1/`=`nquantiles'-1' {
			local cutvallist`=ceil(`i'/`maxlist')' `cutvallist`=ceil(`i'/`maxlist')'' r(r`i')
		}
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
		
		if r(r)==0 {
			di as error "cutpoints() all missing"
			exit 2000
		}
		else {
			local nquantiles = r(r) + 1
			
			* Store quantile boundaries in list
			local maxlist 245
			forvalues i=1/`r(r)' {
				local cutvallist`=ceil(`i'/`maxlist')' `cutvallist`=ceil(`i'/`maxlist')'' `cutvals'[`i',1]
			}
		}
	}
	else { /***** CUTVALUES *****/
		if "`wt'"!="" | "`randvar'"!="" | "`ALTdef'"!="" | `randcut'!=1 | `nquantiles'!=2 | `randn'!=-1 {
			di as error "cutvalues() cannot be used with nquantiles(), altdef, randvar(), randcut(), randn() or weights"
			exit 198
		}
		
		* parse numlist
		numlist "`cutvalues'"
		local maxlist=-1
		local cutvallist1 `"`r(numlist)'"'
		local nquantiles=wordcount(`"`r(numlist)'"')+1
	}

	* Pick data type for quantile variable
	if (`nquantiles'<=maxbyte()) local qtype byte
	else if (`nquantiles'<=maxint()) local qtype int
	else local qtype long

	* Create quantile variable
	local cutvalcommalist : subinstr local cutvallist1 " " ",", all
	qui gen `qtype' `varlist'=1+irecode(`exp',`cutvalcommalist') if `touse'
	
	forvalues i=2/`=ceil((`nquantiles'-1)/`maxlist')' {
		local cutvalcommalist : subinstr local cutvallist`i' " " ",", all
		qui replace `varlist'=1 + `maxlist'*(`i'-1) + irecode(`exp',`cutvalcommalist') if `varlist'==1 + `maxlist'*(`i'-1)
	}

	label var `varlist' "`nquantiles' quantiles of `exp'"
	
	* Return values
	if ("`samplesize'"!="") return scalar n = `samplesize'
	else return scalar n = .
	
	return scalar N = `popsize'
	
	local c=`nquantiles'-1
	forvalues j=`=max(ceil((`nquantiles'-1)/`maxlist'),1)'(-1)1 {
		tokenize `"`cutvallist`j''"' /*XX*/
		forvalues i=`: word count `cutvallist`j'''(-1)1 {
			return scalar r`c' = ``i''
			local --c
		}
	}

end


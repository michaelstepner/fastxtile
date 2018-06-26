*! version 2.0.0beta2  26jun2018  Michael Stepner, stepner@mit.edu

/* CC0 license information:
To the extent possible under law, the author has dedicated all copyright and related and neighboring rights
to this software to the public domain worldwide. This software is distributed without any warranty.

This code is licensed under the CC0 1.0 Universal license.  The full legal text as well as a
human-readable summary can be accessed at http://creativecommons.org/publicdomain/zero/1.0/
*/

* Why did I include a formal license? Jeff Atwood gives good reasons: http://www.codinghorror.com/blog/2007/04/pick-a-license-any-license.html


program define fastxtile, rclass byable(recall, noheader)
	version 11

	* Parse weights, if any
	_parsewt "aweight fweight pweight" `0' 
	local 0  "`s(newcmd)'" /* command minus weight statement */
	local wt "`s(weight)'"  /* contains [weight=exp] or nothing */

	* Extract parameters
	syntax name=/exp [if] [in] [,Nquantiles(integer 2) Cutpoints(varname numeric) ALTdef ///
		CUTValues(numlist ascending) randvar(varname numeric) randcut(real 1) randn(integer -1)]
	local qvar `namelist'
	
	* Check that the variable to be created does not already exist
	if (_byindex()==1) confirm new variable `qvar'
	
	* If there's a by-var, use dots
	if (_by()==1 & _byindex()==1) noisily _dots 0, title(Looping over by-groups:)
	
	* Mark observations which will be placed in quantiles
	tempvar touse
	mark `touse' `if' `in' `wt'
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
			tempvar sample
			mark `sample' `wt' if `touse' & `randvar'<=`randcut'
			
			qui count if `sample'
			local samplesize=r(N)
		}
		else {
			local sample `touse'
			local samplesize `popsize'
		}

		* Error checks
		if (`nquantiles' < 2) {
			di as error "nquantiles() must be greater than or equal to 2"
			exit 198
		}
		else if (`nquantiles' > `samplesize' + 1) {
			if (_by()==0) {
				if ("`randvar'"=="") di as error "nquantiles() must be less than or equal to the number of observations [`r(N)'] plus one"
				else di as error "nquantiles() must be less than or equal to the number of sampled observations [`r(N)'] plus one"
				exit 198
			}
			else {
				noisily _dots `=_byindex()' 1
				exit
			}
		}

		* Compute quantile boundaries
		_pctile `exp' if `sample' `wt', nq(`nquantiles') `altdef'

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

		* Find quantile boundaries from cutpoints var
		mata: process_cutp_var("`cutpoints'")

		* Store quantile boundaries in list
		if r(nq)==1 {
			di as error "cutpoints() all missing"
			exit 2000
		}
		else {
			local nquantiles = r(nq)
			
			local maxlist 248
			forvalues i=1/`=`nquantiles'-1' {
				local cutvallist`=ceil(`i'/`maxlist')' `cutvallist`=ceil(`i'/`maxlist')'' r(r`i')
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
		local nquantiles : word count `r(numlist)'
		local ++nquantiles
	}

	* Pick data type for quantile variable
	if (`nquantiles'<=maxbyte()) local qtype byte
	else if (`nquantiles'<=maxint()) local qtype int
	else local qtype long

	* Create quantile variable
	local cutvalcommalist : subinstr local cutvallist1 " " ",", all
	
	cap confirm new variable `qvar'
	if (_rc==0) qui gen `qtype' `qvar'=1+irecode(`exp',`cutvalcommalist') if `touse'
	else qui replace `qvar'=1+irecode(`exp',`cutvalcommalist') if `touse'
	
	forvalues i=2/`=ceil((`nquantiles'-1)/`maxlist')' {
		local cutvalcommalist : subinstr local cutvallist`i' " " ",", all
		qui replace `qvar'=1 + `maxlist'*(`i'-1) + irecode(`exp',`cutvalcommalist') if `touse' & `qvar'==1 + `maxlist'*(`i'-1)
	}

	if (_by()==1) nois _dots `=_byindex()' 0
	
	if _bylastcall() {
	
		* Label var
		if (_by()==1) label var `qvar' "`nquantiles' quantiles of `exp'; by `_byvars'"
		else label var `qvar' "`nquantiles' quantiles of `exp'"
		
		* Return values
		if ("`samplesize'"!="") return scalar n = `samplesize'
		else return scalar n = .
		
		return scalar N = `popsize'
		
		local c=`nquantiles'-1
		forvalues j=`=max(ceil((`nquantiles'-1)/`maxlist'),1)'(-1)1 {
			tokenize `"`cutvallist`j''"'
			forvalues i=`: word count `cutvallist`j'''(-1)1 {
				return scalar r`c' = ``i''
				local --c
			}
		}
		
	}
	
end


version 11
set matastrict on

mata:

void process_cutp_var(string scalar var) {

	// Load and sort cutpoints	
	real colvector cutp
	cutp=sort(st_data(.,st_varindex(var)),1)
	
	// Return them to Stata
	stata("clear results")
	real scalar ind
	ind=1
	while (cutp[ind]!=.) {
		st_numscalar("r(r"+strofreal(ind)+")",cutp[ind])
		ind=ind+1
	}
	st_numscalar("r(nq)",ind)
	
}

end


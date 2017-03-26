
**********************
***** Initialize *****
**********************
clear all
cap log close _all
discard
set seed 8943159


*******************************
***** Define test program *****
*******************************

cap program drop test_fastxtile
program define test_fastxtile

	version 12.1

	syntax, [Length(integer 1000) Width(integer 1) integer reps(integer 1) BYgroups(integer 1) csvwrite(string) Nquantiles(integer 2) verbose *]

	* Prep dataset
	clear
	qui set obs `length'
	
	if ("`integer'"=="integer") {
		gen int rand=ceil(runiform()*500)  // random integers between {1,...,500}
	}
	else {
		local integer float
		gen float rand=rnormal()
	}
	
	if (`bygroups'>1) {
		gen int byvar=ceil(runiform()*`bygroups')
		sort byvar
	}
	
	if (`width'>1) {
		forvalues i=2/`width' {
			gen float w`i'=rnormal()
		}
	}
	
	di "Testing fastxtile: length `length', width `width', `integer', `reps' reps, `nquantiles' quantiles, `bygroups' by-groups"
		
	
	* Perform testing
	
	timer clear
	
	tempvar fxt xt
	
	forvalues rep=1/`reps' {
	
		local coinflip=round(runiform())
		
		if (`bygroups'>1) {
			timer on 1
			by byvar: fastxtile `fxt'=rand, nq(`nquantiles') `options'
			timer off 1
			
			timer on 2
			forvalues b=1/`bygroups' {
				xtile xt`b'=rand if byvar==`b', nq(`nquantiles') `options'

				if (`b'==1) local qlist xt`b'
				else local qlist `qlist',xt`b'
			}
			gen `xt'=min(`qlist')
			forvalues b=1/`bygroups' {
				drop xt`b'
			}
			timer off 2	
		}
		else if (`coinflip'==0) {
			timer on 1
			fastxtile `fxt'=rand, nq(`nquantiles') `options'
			timer off 1
			
			timer on 2
			xtile `xt'=rand, nq(`nquantiles')  `options'
			timer off 2
		}
		else {
			timer on 2
			xtile `xt'=rand, nq(`nquantiles')  `options'
			timer off 2
			
			timer on 1
			fastxtile `fxt'=rand, nq(`nquantiles')  `options'
			timer off 1		
		}
				
		* Check that output is identical
		assert `fxt'==`xt'
		drop `fxt' `xt'
	}
	
	* Display timings
	qui timer list
	if ("`verbose'"=="verbose") {
		di "fastxtile: " (r(t1)/r(nt1)) " seconds"
		di "xtile: " (r(t2)/r(nt2)) " seconds"
	}
	
	
	* Output timings
	if "`csvwrite'"!="" {
		file write `csvwrite' (`length') "," (`width') ",`integer'," (`reps') "," (`nquantiles') "," (`bygroups') ","
		file write `csvwrite' (r(t2)/r(nt2)) "," (r(t1)/r(nt1)) "," ( (r(t2)/r(nt2)) / (r(t1)/r(nt1)) ) _n
	}

end


**********************
*** Accuracy tests ***
**********************

*** Non-unique quantile (example data taken from pctile PDF manual)

clear
set obs 11
gen bp=98 in 1
replace bp = 100 in 2
replace bp = 104 in 3
replace bp = 110 in 4
replace bp = 120 in 5
replace bp = 120 in 6
replace bp = 120 in 7
replace bp = 120 in 8
replace bp = 125 in 9
replace bp = 130 in 10
replace bp = 132 in 11


* Test manually: nq(5)
xtile quint_xt = bp, nq(5)
fastxtile quint_fxt = bp, nq(5)

assert quint_xt==quint_fxt

* Test manually: pctile, nq(5)
pctile pct = bp, nq(5) genp(percent)
xtile quint2_xt = bp, cutp(pct)
fastxtile quint2_fxt = bp, cutp(pct)

assert quint_xt==quint2_xt
assert quint2_xt==quint2_fxt

list, sepby(quint_xt)

* What about a cutpoint outside the range of the data?
replace pct=-1 in 5
xtile quint3_xt = bp, cutp(pct)
fastxtile quint3_fxt = bp, cutp(pct)
assert quint3_xt==quint3_fxt

list, sepby(quint3_xt)


*** By-able

clear
set obs 4
gen lvl=_n*10
gen grp=_n>2

bys grp: fastxtile fxt=lvl, nq(2)

assert fxt[1]==1
assert fxt[2]==2
assert fxt[3]==1
assert fxt[4]==2


*******************
*** Speed tests ***
*******************

*** cutpoints vs cutvalues

clear
set obs 1000000
gen rand=rnormal()
gen cutp=.
replace cutp=-0.5 in 1
replace cutp=0 in 2
replace cutp=0.5 in 3


timer clear

timer on 1
fastxtile fxt_cutp=rand, cutpoints(cutp)
timer off 1

timer on 2
fastxtile fxt_cutv=rand, cutv(-0.5 0 0.5)
timer off 2

assert fxt_cutp==fxt_cutv

timer list


*** Random sampling

clear
set obs 10000000
gen rand=rnormal()

timer clear
local i 0

forvalues r=0.1(0.1)1 {
	local ++i
	
	timer on `i'
	fastxtile fxt_rand`i'=rand, nq(20) randcut(`r')
	timer off `i'
	
	return list
}


timer list


*** Automatic speed tests

* Create date and time stamp
local date = date(c(current_date), "DMY")
local time = subinstr(c(current_time),":","",.)
local dd : di %02.0f day(`date')
local mm : di %02.0f month(`date')
local yyyy = year(`date')

local timestamp "`yyyy'-`mm'-`dd'-`time'"

* Open output files
tempname outcsv outlog

log using "testresults/fastxtile_tests_`timestamp'.txt", replace text name(`outlog')

file open `outcsv' using "testresults/fastxtile_tests_`timestamp'.csv", write text replace
file write `outcsv' "length,width,datatype,reps,nquantiles,by-groups,xtile secs,fastxtile secs,faster factor" _n

* Print system characteristics to log
which fastxtile

di "Stata version `c(stata_version)'; `c(born_date)'"
di "`c(bit)' bit"

if c(MP)==1 di "Stata MP"
else if c(SE)==1 di "Stata SE"
else di "Stata `c(flavor)'"

di "Processors: `c(processors)'"

di "Operating system: `c(os)' `c(osdtl)'. `c(machine_type)'"

di "Running at: `timestamp'"

* Perform tests
foreach b in 1 2 20 200 {
	foreach l in 1e3 1e5 1e7 {
		foreach w in 1 10 100 {
			foreach nq in 2 10 100 {
				foreach dtype in integer "" {
				
					if (`l'<=1e3) local reps 200
					else if (`l'<=1e5) local reps 20
					else local reps 2
				
					test_fastxtile, l(`l') w(`w') nq(`nq') bygroups(`b') `dtype' reps(`reps') csvwrite(`outcsv')
					
				}
			}
		}
	}
}

* Close outputs
file close `outcsv'
log close `outlog'

*! Eugene Lee, Myoung-jae Lee, and Goeun Lee, "Extending Policies to Control Groups or Past Periods in Nonlinear Difference in Differences", Applied Economics, 2025, forthcoming.
*! version 0.0.1 04aug2026

set more off
clear all
*cd ""  // change directory in which ohie2021.txt or bumi2021.txt exist

/* when using ohie2021.txt */
import delimited id wage y s q d age agesq medicost female white msa diabetes ///
	copd chf educoll employed famincat using ohie2021.txt, clear

global y "y"
global s "s"
global q "q"
global d "d"
global w "age agesq medicost female white msa diabetes copd chf educoll employed famincat"
* */

/* when using bumi2021.txt *
import delimited year y s q d male single china india ed2 ed3 ed4 ind2 ind3 ind4 ///
	using bumi2021.txt, clear

global y "y"
global s "s"
global q "q"
global d "d"
global w "male single china india ed2 ed3 ed4 ind2 ind3 ind4"
* */

//---------- begin: programs ----------
// syntax: polextndd depvar s q d covariates [if] [in]
capture program drop polextndd
program define polextndd
	version 17.0
	syntax varlist [if] [in]
	marksample touse
	
	local yv "`1'"
	local sv "`2'"
	local qv "`3'"
	macro shift 4
	local wv "`*'"
	
	// probit
	qui probit `varlist' if `touse', vce(opg)
	tempname bhat
	tempvar touse1 xb cdf ss0
	mat `bhat' = e(b)
	qui gen `touse1' = e(sample)
	qui predict `xb' if `touse1', xb
	qui predict `cdf' if `touse1', pr
	qui gen `ss0' = ((`yv'-`cdf')*normalden(`xb'))/(`cdf'*(1-`cdf')) if `touse1'
	
	// mata
	mata: polextndd_est("`sv'", "`qv'", "`wv'", "`bhat'", "`ss0'", "`touse1'")
	
	// display
	local namelist mte mce1 mce2 mce3 mce4 mce5
	local txtmte  = "Actual, Treat Q = 1 at t = 1 (tv):"
	local txtmce1 = "    Treat everybody at t = 0 (tv):"
	local txtmce2 = "    Treat everybody at t = 1 (tv):"
	local txtmce3 = "        Treat Q = 1 at t = 0 (tv):"
	local txtmce4 = "        Treat Q = 0 at t = 0 (tv):"
	local txtmce5 = "        Treat Q = 0 at t = 1 (tv):"
	di ""
	foreach v of local namelist {
		local k = "`v'_tv"
		di as txt "`txt`v''" _skip(2) as res %5.4fc `v' ///
			_skip(2) as txt "( " as res %3.2fc `k' as txt " )"
	}
end

version 17.0
mata:
mata clear
mata set matastrict on

void polextndd_est(string scalar sv, string scalar qv, string scalar wv, ///
	string scalar bhat, string scalar ss0, string scalar touse1) {
	
	real scalar nn, n2, bhs, bhq, bhd, kk, mte01, mte01tv, mte00, mce00, mce00tv, ///
		mte11, mte11tv, mce_1, mce_1tv, mte10, mce10, mce10tv, mce_0, mce_0tv
	real rowvector bh, mhe01, mhe00, mhe11, mhe10
	real colvector used, onev, zero, q1, q0, s1, s0, ss1, wb, f00, f01, f10, f11n, ///
		te01, te00, te11, te10
	real matrix ww, ss, inf_p, i01a, i01b, i00a, i00b, i11a, i11b, i10a, i10b
	
	// variables
	used = st_data(.,touse1)
	nn = sum(used)
	n2 = nn^2
	onev = J(nn,1,1)
	zero = J(nn,1,0)
	q1 = select(st_data(.,qv),used:==1)
	q0 = onev:-q1
	s1 = select(st_data(.,sv),used:==1)
	s0 = onev:-s1
	ww = (select(st_data(.,wv),used:==1),onev)
	
	// influence function
	ss1 = select(st_data(.,ss0),used:==1)
	ss = ss1:*(s1,q1,s1:*q1,ww)
	inf_p = invsym(cross(ss,ss)/nn)*ss'
	
	// estimates
	bh = st_matrix(bhat)
	bhs = bh[1,1]
	bhq = bh[1,2]
	bhd = bh[1,3]
	kk = cols(bh)
	
	// fitted values
	wb = ww*bh[1,4..kk]'
	f00 = wb
	f01 = wb:+bhs
	f10 = wb:+bhq
	f11n = f01:+bhq
	
	// treat Q = 0 at t = 1
	te01 = gettei(q0,s1,f01:+bhd,f01)
	mte01 = mean(te01)
	i01a = (onev,zero,onev,ww)
	i01b = (onev,zero,zero,ww)
	mhe01 = getmhe(q0,s1,f01:+bhd,f01,i01a,i01b)
	mte01tv = mtetv(te01,mte01,mhe01,inf_p,n2)
	
	// treat Q = 0 at t = 0
	te00 = gettei(q0,s0,f00:+bhd,f00)
	mte00 = mean(te00)
	mce00 = mte00:+mte01
	i00a = (zero,zero,onev,ww)
	i00b = (zero,zero,zero,ww)
	mhe00 = getmhe(q0,s0,f00:+bhd,f00,i00a,i00b)
	mce00tv = mtetv(te00:+te01,mce00,mhe00:+mhe01,inf_p,n2)
	
	// treat Q = 1 at t = 1 (actual)
	te11 = gettei(q1,s1,f11n:+bhd,f11n)
	mte11 = mean(te11)
	i11a = (onev,onev,onev,ww)
	i11b = (onev,onev,zero,ww)
	mhe11 = getmhe(q1,s1,f11n:+bhd,f11n,i11a,i11b)
	mte11tv = mtetv(te11,mte11,mhe11,inf_p,n2)
	
	// treat everybody at t = 1
	mce_1 = mte01:+mte11
	mce_1tv = mtetv(te01:+te11,mce_1,mhe01:+mhe11,inf_p,n2)
	
	// treat Q = 1 at t = 0
	te10 = gettei(q1,s0,f10:+bhd,f10)
	mte10 = mean(te10)
	mce10 = mte10:+mte11
	i10a = (zero,onev,onev,ww)
	i10b = (zero,onev,zero,ww)
	mhe10 = getmhe(q1,s0,f10:+bhd,f10,i10a,i10b)
	mce10tv = mtetv(te10:+te11,mce10,mhe10:+mhe11,inf_p,n2)
	
	// treat everybody at t = 0
	mce_0 = mce00:+mce10
	mce_0tv = mtetv(te00:+te01:+te10:+te11,mce_0,mhe00:+mhe01:+mhe10:+mhe11,inf_p,n2)
	
	// mata to stata
	st_numscalar("mte",mte11)
	st_numscalar("mte_tv",mte11tv)
	st_numscalar("mce1",mce_0)
	st_numscalar("mce1_tv",mce_0tv)
	st_numscalar("mce2",mce_1)
	st_numscalar("mce2_tv",mce_1tv)
	st_numscalar("mce3",mce10)
	st_numscalar("mce3_tv",mce10tv)
	st_numscalar("mce4",mce00)
	st_numscalar("mce4_tv",mce00tv)
	st_numscalar("mce5",mte01)
	st_numscalar("mce5_tv",mte01tv)
}

real colvector gettei(real colvector qi, real colvector si, ///
	real colvector fia, real colvector fib) {
	return(qi:*si:*(normal(fia):-normal(fib)))
}

real rowvector getmhe(real colvector qi, real colvector si, ///
	real colvector fia, real colvector fib, real matrix i01a, real matrix i01b) {
	return(mean(qi:*si:*(normalden(fia):*i01a:-normalden(fib):*i01b)))
}

real scalar mtetv(real colvector tei, real scalar mte0, real rowvector hte0, ///
	real matrix inf_p0, real scalar n0) {
	real colvector inf0
	real scalar var0
	inf0 = (hte0*inf_p0)':+tei:-mte0
	var0 = cross(inf0,inf0):/n0
	return(mte0:/sqrt(var0))
}

mata set matastrict off
end
//---------- end: programs ------------

// probit
probit ${y} ${s} ${q} ${d} ${w}, vce(opg)

// estimation results
polextndd ${y} ${s} ${q} ${d} ${w}

set more on

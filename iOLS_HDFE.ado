* 15/12/2021 : corrected "cross" in S.E. inversion to increase speed. Note: this required deleting the diagonalization step.
* 15/12/2021 : corrected iteration logical specification
* 16/12/2021 : corrected absorb from varlist to string, as in ppmlhdfe
* 22/12/2021 : coded with matrix multiplication instead of pre-canned program
* 22/12/2021 : added convergence control (limit and maximum)
* 04/01/2022 : added constant + checks for convergence
cap program drop iOLS_HDFE
program define iOLS_HDFE, eclass 
syntax varlist [if] [in] [aweight pweight fweight iweight] [, DELta(real 1) LIMit(real 1e-8)  MAXimum(real 10000) ABSorb(string)  Robust CLuster(string)]        

//	syntax [anything] [if] [in] [aweight pweight fweight iweight] [, DELta(real 1)  ABSorb(string) LIMit(real 0.00001) MAXimum(real 1000) Robust CLuster(varlist numeric)]
	marksample touse
	markout `touse'  `cluster', s     

	preserve
	quietly keep if `touse'
	if  "`robust'" !="" {
		local opt1  = "`robust' "
	}
	if "`cluster'" !="" {
		local opt2 = "vce(cluster `cluster') "
	}
    if "`absorb'" !="" {
		local opt3 = "absorb(`absorb') "
	}
	local option = "`opt1'`opt2'"	
	local list_var `varlist'
	* get depvar and indepvar
	gettoken depvar list_var : list_var
	//gettoken indepvar list_var : list_var, p("(")
	gettoken _rhs list_var : list_var, p("(")
foreach var of varlist  `_rhs' {
quietly drop if missing(`var')	
}
*** check seperation : code from "ppml"
 tempvar zeros                            																						// Creates regressand for first step
 quietly: gen `zeros'=1 	
foreach var of varlist `absorb'{
tempvar group 
quietly egen `group' = group(`var')
tempvar max_group 
bys `group' : egen `max_group' = max(`depvar') if `touse'
tempvar min_group 
bys `group' : egen `min_group' = min(`depvar') if `touse'
quietly: replace `zeros' = 0 if `min_group' == 0  & `max_group' == 0 & `touse'
cap drop `group'
}
 tempvar logy                            																						// Creates regressand 
 quietly: gen `logy'=log(`depvar') if (`touse')&(`depvar'>0)
 quietly: reghdfe `logy' `_rhs'    if (`touse')&(`depvar'>0)	, absorb(`absorb')
// quietly: replace `touse' = 0 if (e(sample)==0) & (`touse')&(`depvar'>0)
// Initialize observations selector
 local _drop ""                                                                                                // List of regressors to exclude
 local indepvar ""     
    foreach x of varlist `_rhs' {   
      if (_se[`x']==0) {                                                                                       // Try to include regressors dropped in the
          qui summarize `x' if (`depvar'>0)&(`touse'), meanonly
          local _mean=r(mean)
          qui summarize `x' if (`depvar'==0)&(`touse')
          if (r(min)<`_mean')&(r(max)>`_mean'){                                            // Include regressor if conditions met and
              local indepvar "`indepvar' `x'"                                                                        // strict is off 
          }
          else{
              qui su `x' if `touse', d                                                                         // Otherwise, drop regressor
              local _mad=r(p50)
              qui inspect  `x'  if `touse'                                                                         
              qui replace `zeros'=0 if (`x'!=`_mad')&(r(N_unique)==2)&(`touse')                                // Mark observations to drop
              local _drop "`_drop' `x'"
          }
      }                                                                                                        // End of LOOP 1.1 
      if _se[`x']>0 {                                                                                          // Include safe regressors: LOOP 1.2
      local indepvar "`indepvar' `x'" 
      }                                                                                                        // End LOOP 1.2
    }   
 qui su `touse' if `touse', mean                                                                               // Summarize touse to obtain N
 local _enne=r(sum)                                                                                            // Save N
 qui replace `touse'=0 if (`zeros'==0)&("`keep'"=="")&(`depvar'==0)&(`touse')                                  // Drop observations with perfect fit
 di                                                                                                              // if keep is off
 local k_excluded : word count `_drop'                                                                      // Number of variables causing perfect fit
 di in green "Number of non-absorbed regressors excluded to ensure that the estimates exist: `k_excluded'" 
 if ("`_drop'" != "") di "Excluded regressors: `_drop'"                                                        // List dropped variables if any
 qui su `touse' if `touse', mean
 local _enne = `_enne' - r(sum)                                                                                // Number of observations dropped
 di in green "Number of observations excluded: `_enne'" 
 local _enne =  r(sum)
quietly keep if `touse'	
	** drop collinear variables
    _rmcoll `indepvar', forcedrop 
	local var_list `r(varlist)'
		*** FWL Theorem Application
	cap drop M0_*
	cap drop Y0_*
	cap drop xb_hat*
	tempvar new_sample
	quietly hdfe `r(varlist)' [`weight'] , absorb(`absorb') generate(M0_) sample(`new_sample') 
quietly:	keep if `new_sample'
	tempvar y_tild  
	quietly gen `y_tild' = log(`depvar' + 1)
	cap drop `new_sample'
	quietly	hdfe `y_tild' , absorb(`absorb') generate(Y0_) sample(`new_sample') 
quietly:	keep if `new_sample'
	mata : X=.
	mata : PX=.
	mata : y_tilde =.
	mata : Py_tilde =.
	mata : y =.
	mata : st_view(X,.,"`var_list'")
	mata : st_view(PX,.,"M0_*")
	mata : st_view(y_tilde,.,"`y_tild'")
	mata : st_view(Py_tilde,.,"Y0_")
	mata : st_view(y,.,"`depvar'")	
	* prepare  future inversions 
	mata : invPXPX = invsym(cross(PX,PX))
	mata : beta_initial = invPXPX*cross(PX,Py_tilde)
	mata : beta_t_1 = beta_initial // needed to initialize
	mata : beta_t_2 = beta_initial // needed to initialize
	mata : q_hat_m0 = 0
	local k = 1
	local eps = 1000	
	mata: q_hat = J(`maximum', 1, .)
	*** Iterations iOLS
	_dots 0
	while ((`k' < `maximum') & (`eps' > `limit' )) {
	mata: xb_hat_M = PX*beta_initial 
	mata: xb_hat_N = X*beta_initial
	mata: fe = y_tilde - Py_tilde + xb_hat_M - xb_hat_N
	mata: xb_hat = xb_hat_N + fe
	* Update d'un nouveau y_tild et regression avec le nouvel y_tild
	mata: alpha = log(mean(y:*exp(-xb_hat)))
	mata: y_tilde = log(y + `delta'*exp(xb_hat :+ alpha )) :-mean(log(y + `delta'*exp(xb_hat :+ alpha)) -xb_hat :- alpha  )
	cap drop `y_tild' 
	quietly mata: st_addvar("double", "`y_tild'")
	mata: st_store(.,"`y_tild'",y_tilde)
	cap drop Y0_
    quietly hdfe `y_tild' [`weight'] , absorb(`absorb') generate(Y0_)
	mata : st_view(Py_tilde,.,"Y0_")
	* OLS
	mata: beta_new = invPXPX*cross(PX,Py_tilde)
	mata: criteria = mean(abs(beta_initial - beta_new):^(2))
	mata: st_numscalar("eps", criteria)
	mata: st_local("eps", strofreal(criteria))
		* safeguard for convergence.
	if `k'==`maximum'{
		  di "There has been no convergence so far: increase the number of iterations."  
	}
	if `k'>4{
	mata: q_hat[`k',1] = mean(log( abs(beta_new-beta_initial):/abs(beta_initial-beta_t_2)):/log(abs(beta_initial-beta_t_2):/abs(beta_t_2-beta_t_3)))	
	mata: check_3 = abs(mean(q_hat)-1)
		if mod(`k'-4,50)==0{
    mata: q_hat_m =  mm_median(q_hat[((`k'-49)..`k'),.] ,1)
	mata: check_1 = abs(q_hat_m - q_hat_m0)
	mata: check_2 = abs(q_hat_m-1)
	mata: st_numscalar("check_1", check_1)
	mata: st_local("check_1", strofreal(check_1))
	mata: st_numscalar("check_2", check_2)
	mata: st_local("check_2", strofreal(check_2))
	mata: st_numscalar("check_3", check_3)
	mata: st_local("check_3", strofreal(check_3))
	mata: q_hat_m0 = q_hat_m
		if ((`check_1'<1e-4)&(`check_2'>1e-2)) {
di "delta is too small to achieve convergence -- update to larger value"
	local k = `maximum'
		}
		if ((`check_3'>0.5) & (`k'>500)) {
	local k = `maximum'
di "q_hat too far from 1"
		}
					  }
	}
	if `k'>2 { // keep in memory the previous beta_hat for q_hat 
	mata:   beta_t_3 = beta_t_2
	mata:   beta_t_2 = beta_initial
	}
	mata: beta_initial = beta_new
	local k = `k'+1
	_dots `k' 0
	}
	*** Calcul de la bonne matrice de variance-covariance

 *** Calcul de la matrice de variance-covariance
 	foreach var in `indepvar' {     // rename variables for last ols
	quietly	rename `var' TEMP_`var'
	quietly	rename M0_`var' `var'
	}	
cap _crcslbl Y0_ `depvar'
	quietly reg Y0_ `indepvar' if `touse' [`weight'`exp'], `option'  noconstant
	* Calcul du "bon" residu
	quietly predict xb_hat, xb 
	quietly gen ui = `depvar'*exp(-xb_hat)
	quietly gen weight = ui/(`delta'+ ui)
	mata : weight= st_data(.,"weight")
*** rename variables
	foreach var in `indepvar' {      // rename variables back
	quietly	rename `var' M0_`var'
	quietly	rename TEMP_`var' `var'
	}
* Calcul de Sigma_0, de I-W, et de Sigma_tild
	matrix beta_final = e(b) // 	mata: st_matrix("beta_final", beta_new)
	matrix Sigma = e(V)
	mata : Sigma_hat = st_matrix("Sigma")
	mata : Sigma_0 = (cross(PX,PX):/rows(PX))*Sigma_hat*(cross(PX,PX):/rows(PX)) // recover original HAC 
	mata : invXpIWX = invsym(cross(PX, weight, PX):/rows(X)) 
	mata : Sigma_tild = invXpIWX*Sigma_0*invXpIWX*rows(X)	
	mata: Sigma_tild = (Sigma_tild+Sigma_tild'):/2
    mata: st_matrix("Sigma_tild", Sigma_tild) // used in practice
	*** Stocker les resultats dans une matrice
	local names : colnames e(b)
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
   ereturn post beta_final Sigma_tild , obs(`=e(N)') depname(`depvar') esample(`touse')  dof(`=e(df_r)')
   restore 
   
ereturn scalar delta = `delta'
ereturn scalar eps =   `eps'
ereturn scalar niter =  `k'
ereturn local cmd "iOLS_HDFE"
ereturn local vcetype `option'
di in gr _col(55) "Number of obs = " in ye %8.0f e(N)
ereturn display

end

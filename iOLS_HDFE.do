program define iOLS_HDFE, eclass 
	syntax [anything] [if] [in] [aweight pweight fweight iweight] [, DELta(real 1)  ABSorb(varlist) Robust CLuster(varlist numeric)]
	marksample touse
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
	local list_var `anything'
	* get depvar and indepvar
	gettoken depvar list_var : list_var
	gettoken indepvar list_var : list_var, p("(")
    * get endogenous variables and instruments
	gettoken endog list_var : list_var, bind
	gettoken endog endog : endog, p("(")
    gettoken endog instr_temp : endog , p("=")
    gettoken equalsign instr_temp : instr_temp , p("=")
	gettoken instr instr_temp : instr_temp, p(")")
	
	*di `"`depvar'"'
	*di `"`indepvar'"'
	*di `"`endog'"'
	*di `"`instr'"'
	
	*** FWL Theorem Application
	cap drop M0_*
	cap drop fe
	cap drop Y0_*
	cap drop xb_hat*
	quietly hdfe `indepvar'  [`weight'] , absorb(`absorb') generate(M0_)
		tempname DF_ADAPT
		scalar `DF_ADAPT' = e(df_a) //- e(N_hdfe)
	*** Initialisation de la boucle
	tempvar fe
	quietly gen fe = 0
	tempvar y_tild  
	quietly gen `y_tild' = log(`depvar' + `delta'*exp(fe))
	quietly	hdfe `y_tild' [`weight'] , absorb(`absorb') generate(Y0_) 
	quietly reg Y0_ M0_* if `touse', `option'  noconstant  
	matrix beta_new = e(b)
	local k = 0
	local eps = 1000	
	quietly predict xb_hat,xb
			*** ItÃ©rations iOLS
	_dots 0
	while (`k' < 1000 & `eps' > 1e-100){	
		matrix beta_initial = beta_new
		* calcul de constante
		tempvar temp1
	quietly	gen `temp1' = `depvar'*exp(-xb_hat)
	quietly sum `temp1' if e(sample)
	tempname phi_hat
		scalar `phi_hat' = log(`r(mean)')
		cap drop `temp1'
		* Calcul de c_hat
		tempvar temp2
		quietly	gen `temp2' = log(`depvar' + `delta'*exp(`phi_hat' + xb_hat)) - xb_hat - `phi_hat'
		quietly sum `temp2' if e(sample)
		tempname c_hat
		scalar `c_hat' = `r(mean)'
		* Update d'un nouveau y_tild et regression avec le nouvel y_tild
		quietly replace `y_tild' = log(`depvar' + `delta'*exp(xb_hat)) - `c_hat'
		cap drop Y0_
	    quietly hdfe `y_tild' [`weight'] , absorb(`absorb') generate(Y0_) 
	quietly	reg Y0_ M0_* if `touse' [`weight'`exp'], `option'   noconstant
		matrix beta_new = e(b)
		* Update fixed effects	
	cap drop xb_hat_M xb_hat_N
	quietly predict xb_hat_M, xb		// predict with demeaned
		foreach var in `indepvar' { // rename variables for prediction
	quietly	rename M0_`var' TEMP`var'
	quietly	rename `var' M0_`var'
	}
	quietly predict xb_hat_N, xb // predict without demeaned
	foreach var in `indepvar' { // rename variables back to normal
	quietly	rename M0_`var' `var'
	quietly	rename TEMP`var' M0_`var'
	}
	quietly	replace fe = `y_tild'- Y0_ + xb_hat_M - xb_hat_N
	quietly	replace xb_hat = xb_hat_N + fe
		* DiffÃ©rence entre les anciens betas et les nouveaux betas
		matrix diff = beta_initial - beta_new
		mata : st_matrix("abs_diff", abs(st_matrix("diff")))
		mata : st_matrix("abs_diff2", st_matrix("abs_diff"):*st_matrix("abs_diff"))
		mata : st_matrix("criteria", rowsum(st_matrix("abs_diff2"))/cols(st_matrix("abs_diff2")))
		local eps = criteria[1,1]
		local k = `k'+1
		_dots `k' 0	
	}
 *** Calcul de la matrice de variance-covariance
	preserve
	keep if e(sample)	
	* matrice de beta
	matrix beta_final = e(b)
	* Calcul de Sigma_0, de I-W et de Sigma_tild
	local N_DF = e(df_r) -`DF_ADAPT' 
	matrix Sigma = e(V)
	* Calcul du "bon" residu
    quietly sum `y_tild' [`weight'`exp'] if e(sample)
	tempname nobs
	scalar  `nobs' = e(N)
	tempvar ui
    gen `ui' = `depvar'*exp(- xb_hat - `phi_hat')
	tempvar cste
	gen `cste' = 1
	tempvar ui_bis
	quietly gen `ui_bis' = 1 - `delta'/(`delta' + `ui')
	local var_list M0_* //`cste'
	mata : X=.
	mata : IW=.
	mata : st_view(X,.,"`var_list'")
	mata : st_view(IW,.,"`ui_bis'")
	mata : IW = diag(IW)
	mata : Sigma_hat = st_matrix("Sigma")
	mata : Sigma_0 = (X'*X)*Sigma_hat*(X'*X)
	mata : invXpIWX = invsym(X'*IW*X)
	mata : Sigma_tild = invXpIWX*Sigma_0*invXpIWX
    mata : st_matrix("Sigma_tild", Sigma_tild)
	matrix Sigma_tild = Sigma_tild*((`e(df_r)')/(`N_DF')) // adjust DOF
	*matrix list Sigma_tild
	*di `e(df_r)'
	di `N_DF'
	*** Stocker les resultats dans une matrice
	local names : colnames e(b)
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
    ereturn post beta_final Sigma_tild , obs(`=r(N)') depname(`depvar') esample(`touse')  dof(`N_DF')
	restore 

ereturn scalar delta = `delta'
ereturn scalar eps =   `eps'
ereturn scalar niter =  `k'
ereturn local cmd "iOLS_HDFE"
ereturn local vcetype `option'
ereturn display
* drop 
	cap drop M0_* 
	cap drop fe
	cap drop Y0_*
	cap drop xb_hat*
end
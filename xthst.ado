*! xthst
*! Version 1.0 - 13.11.2019
*! Tore Bersvendsen (University of Agder) tore.bersvendsen@uia.no
*! Jan Ditzen (Heriot-Watt University) j.ditzen@hw.ac.uk www.jan.ditzen.net 

capture program drop xthst
program define xthst, rclass sortpreserve
	syntax varlist(min=2 ts) [if] , [partial(varlist ts) NOCONStant ar hac bw(integer -999) WHITEning kernel(string) CRosssectional(string) NOOUTput ]
	
	version 14
	
	qui{
		if "`whitening'" != "" & "`hac'" == "" {
			local hac hac
		}
		
		if "`bw'" != "-999" & "`hac'" == "" {
			local hac hac
		}
		
		if "`hac'" != "" & "`bw'" == "-999" {
			local bw = 0
		}
		
		if "`hac'" == "hac" & "`ar'" == "ar" {
			noi disp as error "Option hac and ar" ,_c
			error 184
			exit
		}
		
		if "`kernel'" == "" {
			local kernel "bartlett"
		}
		
		tempvar touse	
		marksample touse
		
		qui xtset
		local idvar "`r(panelvar)'"
		local tvar "`r(timevar)'"	
		sort `idvar' `tvar'
		
		*** Create cross-sectional averages if requested
		if "`crosssectional'" != "" {
			local 0 `crosssectional'
			syntax varlist(ts) , [cr_lags(numlist)]
			local crosssectional `varlist'
			tempname csa
				if "`cr_lags'" == "" {
					local cr_lags = 0
				}
				xtdcce2_csa `crosssectional' , idvar(`idvar') tvar(`tvar') cr_lags(`cr_lags') touse(`touse') csa(`csa')			
				local csa `r(varlist)'	
				local cross_structure "`r(cross_structure)'"
				markout `touse' `csa'	
		}

		
		
		*** check for partial vars
		if "`partial'" != "" {
			** make sure partialled out vars do not appear on rhs
			local varlist: list varlist - partial
			
			tsrevar `partial'
			local partial `r(varlist)'
		}

		*** check for time series variables and generate tempvar
		tsrevar `varlist'
		tokenize `r(varlist)'
			
		local lhs `1'
		macro shift
		local rhs `*'
		
		if "`noconstant'" == "" {
			tempvar const
			gen double `const' = 1
			local partial `partial' `const'
		}
		
		*** start mata program here
		tempname delta delta_st delta_adj 
		if "`hac'" == "" {
			mata st_matrix("`delta'",deltatest("`lhs'","`rhs'","`partial' `csa'","`idvar' `tvar'","`touse'",`=("`ar'"=="ar")'))
		}
		else {
			mata st_matrix("`delta'",deltatesthac("`lhs'","`rhs'","`partial' `csa'","`idvar' `tvar'","`touse'",`bw',`=("`whitening'"=="whitening")',"`kernel'"))
			local bw = `delta'[3,1]	
		}

	}
	*** Output	
	scalar `delta_adj' = `delta'[2,1]
	scalar `delta_st' = `delta'[1,1]	
	
	if "`nooutput'" == "" {
		noi disp as text "Test for slope homogeneity"
		if "`hac'" == "" {
			noi disp as text "(Pesaran, Yamagata. 2008. Journal of Econometrics)"
			
		}
		else {
			noi disp as text "(Blomquist, Westerlund. 2013. Economic Letters)"
		}
		
		
		
		noi disp "H0: slope coefficients are homogenous"
		di as text "{hline 37}"
		noi disp as result _col(10) "Delta" _col(25) "p-value"
		noi disp as result  _col(7) %9.3f `delta_st' _col(23) %9.3f 2*(1-normal(abs(`delta_st')))
		noi disp as result  _col(2) "adj." _col(7) %9.3f `delta_adj'  _col(23) %9.3f 2*(1-normal(abs(`delta_adj')))
		di as text "{hline 37}"
		if "`hac'" != "" {
			if "`kernel'" == "qs" { 
				local kernel "quadratic spectral (QS)"
			}
			noi disp as txt "HAC Kernel: `kernel' with bandwith " `bw'
		}
		if "`partial'" != "" {
			local partial = subinstr("`partial'","`const'","constant",.)
			noi disp "Variables partialled out: `partial'"
		}
		if "`crosssectional'" != "" {
			if wordcount("`cr_lags'") > 1 {
				local crosssectional_output "`cross_structure'"
			}
			else {
				local crosssectional_output "`crosssectional'"
			}
			display  as text "Cross Sectional Averaged Variables: `crosssectional_output'"
		}
	}
	*** Return
	return clear
	matrix `delta' = (`delta_st' \ `delta_adj')
	matrix rownames `delta' = Delta Delta_adjusted
	matrix colnames `delta' = TestStat.
	return matrix delta = `delta'
	
	tempname delta_p
	matrix `delta_p' = 2*(1-normal(abs(`delta_st'))) \ 2*(1-normal(abs(`delta_adj')))
	matrix rownames `delta_p' = Delta Delta_adjusted
	matrix colnames `delta_p' = p-Value
	
	return matrix delta_p = `delta_p'

	
	if "`hac'" != "" {
		return scalar bw = `bw'
		return local kernel "`kernel'"
	}
	if "`partial'" != "" {
		return local partial "`partial'"
	}
	if "`crosssectional_output'" != "" {
		return local crosssectional "`crosssectional_output'"
	}
	
	 
end
/*
Steps
1. partial out
2. calculate fixed effect estimator
3. calculate sigma2i, beta2i, gives beta2wfe
4. calcualte s_tilde
5. calculate delta


*/

mata:
	function deltatest ( string scalar lhsname,		/// lhs variable
							string scalar rhsname,		/// rhs variables
							string scalar rhspartialname,	/// variables to be partialled out
							string scalar idtname, 		/// id and t variables
							string scalar tousename,	/// touse variable			
							real scalar ar)				/// 1 if ar, 0 if not ar
							
	{
		real matrix Y
		real matrix X
		real matrix idt
		real scalar Nuniq
		real scalar N_g

		/// load data
		Y = st_data(.,lhsname,tousename)
		X = st_data(.,rhsname,tousename)
		
		idt = st_data(.,idtname,tousename)

		
		K1 = 0
		Z = .

		if (rhspartialname[1,1]:!= " ") {
			Z = st_data(.,rhspartialname,tousename)
			K1 = cols(Z)		
		}
		
		Nuniq = uniqrows(idt[.,1])
		N_g = rows(Nuniq)
		K = cols(X)		
		Kpartial = 0
		/// set it as panel
		index = panelsetup(idt[.,1],1)
		
		/// 1. Partialling out
		if (Z[1,1] != .) {
			i = 1
			Kpartial = cols(Z)
			while (i<=N_g) {
				starti = index[i,1]
				endi = index[i,2]
				
				Yi = Y[(starti..endi),.]
				Xi = X[(starti..endi),.]
				Zi = Z[(starti..endi),.]
				
				/// partialling out
				tmp_zz = quadcross(Zi,Zi)
				tmp_zz1 = invsym(tmp_zz)
				
				Y[(starti..endi),.] = Yi - Zi * tmp_zz1*quadcross(Zi,Yi)
				X[(starti..endi),.] = Xi - Zi * tmp_zz1*quadcross(Zi,Xi)
				
				i++
			}

		}
		
		if (ar==1) {
			Kpartial = 0
		}
		
		//// 2 Fe estimates
		tmp_xx = quadcross(X,X)
		tmp_xy = quadcross(X,Y)
		tmp_xx1 = invsym(tmp_xx)
		b_fe = tmp_xx1 * tmp_xy
		resid = Y - X * b_fe
		
		/// 3 calcualte sigma2i, beta2i, gives beta2wfe
		sigma2 = J(N_g,1,.)
		beta2i = J(N_g,K,.)
		
		beta2wfe_up = 0
		beta2wfe_low = J(1,K,0)
		Tavg = 0
		i = 1
		while (i<=N_g) {
			starti = index[i,1]
			endi = index[i,2]

			Yi = Y[(starti..endi),.]
			Xi = X[(starti..endi),.]
			residi = resid[(starti..endi),.]
			
			Ti = rows(Xi)
			Ki = cols(Xi)
	
			sigma2[i] =  residi'residi :/ (Ti - Kpartial - 1)
			
			tmp_xx = quadcross(Xi,Xi)
			tmp_xx1 = invsym(tmp_xx)
			tmp_xy = quadcross(Xi,Yi)
						
			beta2i[i,.] = (tmp_xx1*tmp_xy)'
			beta2wfe_up = beta2wfe_up :+ tmp_xy :/ sigma2[i]
			beta2wfe_low = beta2wfe_low :+  tmp_xx :/sigma2[i]
			
			Tavg = Tavg + Ti
			
			i++
		}
		
		Tavg = Tavg / N_g
		
		beta2wfe_low = invsym(beta2wfe_low)
		
		beta2wfe = beta2wfe_low * beta2wfe_up 
		
		/// 4. calcualte s_tilde
		S_tilde = 0
		i = 1

		while (i <= N_g) {
			starti = index[i,1]
			endi = index[i,2]
			Xi = X[(starti..endi),.]
			beta_i = beta2i[i,.]'
			
			tmp_xx = quadcross(Xi,Xi) :/ sigma2[i]
			
			S_tilde = S_tilde + (beta_i - beta2wfe)' * tmp_xx * (beta_i - beta2wfe)			
			i++
		}

		delta = sqrt(N_g) * (S_tilde/N_g - K) / sqrt(2*K)

		var = 2 * K * (Tavg-K-Kpartial-1)/ (Tavg-Kpartial+1)
		
		delta_adj = sqrt(N_g)*(((S_tilde/N_g)-K)/sqrt(var))
		
		return(delta\delta_adj)
		
		
		
	}
end	


mata:
	function deltatesthac ( string scalar lhsname,		/// lhs variable
							string scalar rhsname,		/// rhs variables
							string scalar rhspartialname,	/// variables to be partialled out
							string scalar idtname, 		/// id and t variables
							string scalar tousename, /// touse variable	
							real scalar bandwith,	///
							real scalar whitening, ///
							string scalar kernel ///
							)			
														
	{

		real matrix Y
		real matrix X
		real matrix idt
		real scalar Nuniq
		real scalar N_g

		/// load data
		Y = st_data(.,lhsname,tousename)
		X = st_data(.,rhsname,tousename)
		
		idt = st_data(.,idtname,tousename)
		
		
		K1 = 0
		Z = .
		if (rhspartialname[1,1]:!= " " ) {
			Z = st_data(.,rhspartialname,tousename)
			K1 = cols(Z)		
		}
		
		Nuniq = uniqrows(idt[.,1])
		N_g = rows(Nuniq)
		K = cols(X)		
		Kpartial = 0
		/// set it as panel
		index = panelsetup(idt[.,1],1)
		Xbar = J(rows(X),cols(X),.)

		/// 1. Partialling out
		if (Z[1,1] != .) {
			i = 1
			Kpartial = cols(Z)
			while (i<=N_g) {
				starti = index[i,1]
				endi = index[i,2]
				
				Yi = Y[(starti..endi),.]
				Xi = X[(starti..endi),.]
				Zi = Z[(starti..endi),.]
				
				/// partialling out
				tmp_zz = quadcross(Zi,Zi)
				tmp_zz1 = invsym(tmp_zz)
				
				Y[(starti..endi),.] = Yi - Zi * tmp_zz1*quadcross(Zi,Yi)
				X[(starti..endi),.] = Xi - Zi * tmp_zz1*quadcross(Zi,Xi)

				Xbar[(starti..endi),.] = J(rows(Xi),1,mean(Xi))	
				
				i++
			}

		}		

		//// 2. Fe estimates
		tmp_xx = quadcross(X,X)
		tmp_xy = quadcross(X,Y)
		tmp_xx1 = invsym(tmp_xx)
		b_fe = tmp_xx1 * tmp_xy
		eps = Y - X * b_fe
		
		uhat = (X - Xbar) :* eps
		//Gamma = (N_g*T,K,.)

		id2 = Nuniq#J(K,1,1)
		index2 = panelsetup(id2[.,1],1)

		// here V is inverse of V!!
		V = J(N_g*K,K,.)
		
		i=1
		while (i<=N_g) {
			starti = index[i,1]
			endi = index[i,2]
			
			start2i = index2[i,1]
			end2i = index2[i,2]
			
			uhati = uhat[(starti..endi),.] 
			
			Ti= rows(uhati)
			
			
			if (whitening == 1) {
				uhatiy =uhati[1..rows(uhati)-1,.]
				uhatix =uhati[2..rows(uhati),.]
				
				tmp_uu = quadcross(uhatix,uhatix)
				tmp_uu1 = invsym(tmp_uu)
				tmp_uxy= quadcross(uhatix,uhatiy)
				
				A = tmp_uu1 * tmp_uxy
				
				uhati = uhatiy - uhatix*A

				Ti = rows(uhati)
				
			}
			
			
			Gammaj = 1/Ti * ((uhati[1..Ti,.]') * (uhati[1..Ti,.]))
			Vi =  (Gammaj + Gammaj')

			if (bandwith == 0 ) {
				bandwith = floor( 4 * (Ti:/100)^(2/9))		
			}
			
			
			j=1
			while (j<=bandwith) {	
				Gammaj = 1/Ti * ((uhati[j+1..Ti,.]') * (uhati[1..Ti-j,.]))

				/// use bartlett kernel 
				if (kernel == "bartlett") {
					kxi = 1-j/bandwith
				}
				else if (kernel == "qs") {
					kxi = j/bandwith
					kxi = 25 / (12 * pi()^2 * kxi^2) * (sin(6*pi()*kxi / 5) / (6*pi() * kxi/5)  - cos(6 * pi() * kxi/5))
				}
				else if (kernel == "truncated") {
					kxi = (j/bandwith <= 1)
				}

				Vi = Vi + kxi * (Gammaj + Gammaj')
				j++
			}	
			
			if (whitening == 1) {
				Vi = invsym(I(K)-A)* Vi * invsym(I(K)-A)'
			}	

			V[start2i..end2i,.] = invsym(Vi)

			i++
		}

		
		
		///estimate beta
		beta_low = J(K,K,0)
		beta_up = J(K,1,0)
		
		i = 1
		while (i<=N_g) {
			starti = index[i,1]
			endi = index[i,2]
		
			Yi = Y[(starti..endi),.]
			Xi = X[(starti..endi),.]
		
			start2i = index2[i,1]
			end2i = index2[i,2]
			
			Vi1 = V[start2i..end2i,.]
			
			Ti = rows(Xi)
			
			Qi = quadcross(Xi,Xi) / Ti
			QiY = quadcross(Xi,Yi)
			beta_low = beta_low + Ti * Qi * Vi1 * Qi
			beta_up = beta_up + Qi * Vi1 * QiY
			i++
		}
		
		beta = invsym(beta_low) * beta_up
		
		/// estimate S_HAC
		
		S_HAC = 0
		Tavg = 0
		i = 1
		while (i<=N_g) {
			starti = index[i,1]
			endi = index[i,2]
		
			Yi = Y[(starti..endi),.]
			Xi = X[(starti..endi),.]
		
			start2i = index2[i,1]
			end2i = index2[i,2]
			
			Vi1 = V[start2i..end2i,.]
			
			Ti = rows(Xi)
			
			tmp_xx = quadcross(Xi,Xi)
			Qi = tmp_xx / Ti
			tmp_xx1 = invsym(tmp_xx)
			tmp_xy = quadcross(Xi,Yi)
			
			betai = tmp_xx1 * tmp_xy

			S_HAC = S_HAC+ Ti * (betai - beta)' * Qi * Vi1 * Qi * (betai - beta)
			
			Tavg = Tavg + Ti
			i++
		}

		Tavg = Tavg / N_g
		
		delta_hac = sqrt(N_g) * (S_HAC / N_g - K) / sqrt(2*K)
		
		var = 2 * K * (Tavg-K-Kpartial-1)/ (Tavg-Kpartial+1)
		
		delta_adj = sqrt(N_g)*(((S_HAC/N_g)-K)/sqrt(var))
		
		return(delta_hac\delta_adj\bandwith)
	}
end

/* Program from xtdcce2 to calculate CSA; creates csa and returns list with tempvars */ 
capture program drop xtdcce2_csa
program define xtdcce2_csa, rclass
        syntax varlist(ts) , idvar(varlist) tvar(varlist) cr_lags(numlist) touse(varlist) csa(string) 
                tsrevar `varlist'
                local varlist `r(varlist)'
                foreach var in `varlist' {
                                local ii `=strtoname("`var'")'
                                tempvar `ii'
                                by `tvar' (`idvar'), sort: egen ``ii'' = mean(`var') if `touse'                         
                                local clist `clist' ``ii''
                        }
                        if "`cr_lags'" == "" {
                                local cr_lags = 0
                        }
                        local i = 1
                        local lagidef = 0
                        foreach var in `clist' {
                                local lagi = word("`cr_lags'",`i')
                                if "`lagi'" == "" {
                                        local lagi = `lagidef'
                                }
                                else {
                                        local lagidef = `lagi'                                  
                                }
                                sort `idvar' `tvar'
                                tsrevar L(0/`lagi').`var'
                                
                                local cross_structure "`cross_structure' `=word("`varlist'",`i')'(`lagi')"
                                local clistfull `clistfull' `r(varlist)'
                                local i = `i' + 1
                        }
                        local i = 1
                        foreach var in `clistfull' {
                                rename `var' `csa'_`i'
                                local clistn `clistn' `csa'_`i'
                                local i = `i' + 1
                        }
                        
                return local varlist "`clistn'"
                return local cross_structure "`cross_structure'"
end


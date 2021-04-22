{smcl}
{* *! version 1.0 22march2021}{...}
{vieweralsosee "[R] poisson" "help poisson"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "ppml" "help ppml"}{...}
{vieweralsosee "ppmlhdfe" "help ppmlhdfe"}{...}
{viewerjumpto "Syntax" "iOLS_HDFE##syntax"}{...}
{viewerjumpto "Description" "iOLS_HDFE##description"}{...}
{viewerjumpto "Citation" "iOLS_HDFE##citation"}{...}
{viewerjumpto "Authors" "iOLS_HDFE##contact"}{...}
{viewerjumpto "Examples" "iOLS_HDFE##examples"}{...}
{viewerjumpto "Description" "iOLS_HDFE##Testing"}{...}
{viewerjumpto "Stored results" "iOLS_HDFE##results"}{...}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:iOLS_HDFE} {hline 2}} Iterated Ordinary Least Squares (iOLS) with High Dimensional Fixed Effects (HDFE) with delta {p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:iOLS_HDFE}
{depvar} [{indepvars}][{absorb}]
{ifin} {it:{weight}} {cmd:,} [{help iOLS_HDFE##options:options}] {p_end}

{marker opt_summary}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab: Standard Errors: Classical/Robust/Clustered}
{synopt:{opt vce}{cmd:(}{help iOLS_HDFE##opt_vce:vcetype}{cmd:)}}{it:vcetype}
may be classical (assuming homoskedasticity), {opt r:obust}, or {opt cl:uster} (allowing two- and multi-way clustering){p_end}
{syntab: Delta}
{synopt:{opt delta}{cmd:(}{help iOLS_HDFE##delta:delta}{cmd:)}}{it:delta} is any strictly positive constant. {p_end}
{syntab: Delta}
{synopt:{opt absorb}{cmd:(}{help iOLS_HDFE##absorb:absorb}{cmd:)}}{it:absorb} include binary and categorical variables. You can generate an interaction between two categories with
"egen group_fe = group(var1 var2)". {p_end}



{marker description}{...}
{title:Description}

{pstd}{cmd:iOLS_HDFE} iterated Ordinary Least Squares with High Dimensional Fixed Effects and delta, as described by {browse "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3444996":Bellego, Benatia, and Pape (2021)}.

{pstd} This package:

{pmore} 1. relies on Stata's OLS reg procedure for estimation. {p_end}

{pmore} 2. assumes the iOLS exogeneity condition with delta E(X'log(delta+U)) = constant. {p_end}

{pmore} 3. allows for high dimensional fixed effects. {p_end}


{title:Background}

{pstd}
{cmd: iOLS_HDFE} estimates iOLS_delta, a solution to the problem of the log of zero.  This method relies on running the "regress" function iteratively whilst simultaneously partialling out
the fixed effects. This provides the reader with the final OLS estimates and allows the use the post-estimation commands available under regress (using M_FE * Y_tilde = log(Y + delta*exp(xb))) as a 
dependent variable with M_FE being the residual making matrix associated with the fixed effect variables. 

{synoptset 22}{...}
{synopthdr: variables}
{synoptline}
{synopt:{it:depvar}} Dependent variable{p_end}
{synopt:{it:indepvars}} List of explanatory variables {p_end}
{synoptline}
{p2colreset}{...}


{marker caveats}{...}
{title:Caveats}

{pstd} Convergence is decided based on coefficients (sum of squared coefficients < 1e-15) and not on the modulus of the contraction mapping.

{pstd} The {help test} postestimation commands are available after {cmd:iOLS_HDFE}.  This command yields 'x_1b' using "predict xb, xb"  and the aggregate fixed effect 'fe'. To obtain y_hat, you will need to also run "gen y_hat = exp(xb+fe)".

{marker contact}{...}
{title:Authors}

{pstd} Christophe Bellego, David Benatia, Louis Pape {break}
CREST - ENSAE - HEC Montréal - Ecole Polytechnique {break}
Contact: {browse "mailto:louis.pape@polytechnique.edu":louis.pape@polytechnique.edu} {p_end}

{marker citation}{...}
{title:Citation}

{pstd}
Bellégo Christophe, Benatia David, and Pape Louis-Daniel, Dealing with Logs and Zeros in Regression Models (2019).
Série des Documents de Travail n° 2019-13.
Available at SSRN: https://ssrn.com/abstract=3444996 

{marker examples}{...}
{title:Examples}


{pstd} We use data on households' trips away from home, as used in {browse "https://www.stata.com/manuals/rivpoisson.pdf":ivpoisson manual}.
to study the effect of cost of transportation (tcost). 
{p_end}
{hline}
{phang2}{cmd:. clear all}{p_end}
{phang2}{cmd:. webuse trip }{p_end}
{phang2}{cmd:. gen outside = trips>0 }{p_end}

{phang2}{cmd:. iOLS_HDFE trips cbd ptn tcost, delta(1) robust absorb(worker weekend) }{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:iOLS_HDFE} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:e(N)}} number of observations{p_end}
{synopt:{cmd:e(sample)}} marks the sample used for estimation {p_end}
{synopt:{cmd:e(eps)}} sum of the absolute differences between the parameters from the last two iterations of iOLS {p_end}
{synopt:{cmd:e(k)}} number of iterations of iOLS{p_end}

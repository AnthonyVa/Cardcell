/****************************************************************************
 * Moreno2011Cell.cpp
 * Based on code from Moreno et al. 2011
 *
 *  Created on: Sep 30, 2014
 *      Author: Anthony Varghese
 ***************************************************************************/



#include <cell/Cells.h>
#include <cell/Moreno2011Cell.h>
#include <string>



#include <gsl/gsl_linalg.h>

#ifndef HAVE_INLINE /* enable compilation of inline get and set functions for GSL vector and matrix */
#define HAVE_INLINE
#endif



/**
 * Cell constructor
 */
Moreno2011Cell::Moreno2011Cell(){
	allocTemps();
	SetInitialConditions();
}

/**
 * Cell destructor
 */
Moreno2011Cell::~Moreno2011Cell(){
	freeTemps();
}

/**
 * Cell copy constructor
 */
Moreno2011Cell::Moreno2011Cell(const Moreno2011Cell& source){
	allocTemps();
	Copy( source );
}

/**
 * Assignment operator
 */
Moreno2011Cell& Moreno2011Cell::operator=(const Moreno2011Cell& rhs){
	if (this != &rhs){
		freeTemps();
		allocTemps();
		Copy( rhs );
	}
	return *this;
}


/*
 * allocTemps - allocate heap memory for Cell object
 */
void Moreno2011Cell::allocTemps(){
	Yn			  = gsl_vector_calloc( nNaChStates );
	Ynph		  = gsl_vector_calloc( nNaChStates );
	Ynp1		  = gsl_vector_calloc( nNaChStates );
	Ynp1i		  = gsl_vector_calloc( nNaChStates );
	DiagonalAn    = gsl_vector_calloc( nNaChStates );
	DiagonalAnInv = gsl_vector_calloc( nNaChStates );
	OffDiagonalAn = gsl_matrix_calloc( nNaChStates, nNaChStates );
}

/*
 * freeTemps - deallocate heap memory for Cell object
 */
void Moreno2011Cell::freeTemps(){
	if ( Yn				!= nullptr) gsl_vector_free( Yn   );
	if ( Ynph			!= nullptr) gsl_vector_free( Ynph );
	if ( Ynp1			!= nullptr) gsl_vector_free( Ynp1 );
	if ( Ynp1i			!= nullptr) gsl_vector_free( Ynp1i);
	if ( DiagonalAn 	!= nullptr) gsl_vector_free( DiagonalAn );
	if ( DiagonalAnInv	!= nullptr) gsl_vector_free( DiagonalAnInv );
	if ( OffDiagonalAn	!= nullptr) gsl_matrix_free( OffDiagonalAn );
}


void Moreno2011Cell::Copy(const Moreno2011Cell& in){
	if (this == &in)
		return;

	Cell_type = in.Cell_type;

	// Copy all the states, auxes and rates over
	for (unsigned int i=0; i<nstates; i++) {
		state[i]	= in.state[i];
		rates[i]	= in.rates[i];
	}
	for (unsigned int i=0; i<naux; i++) {
		currents[i] = in.currents[i];
	}
	// When the data member "" is set to true,
	// it avoids resetting Na channel states in
	// the WT_SCN5A_Initial_Conditions method
	initializedExternally = true;

	peak_slope	= in.peak_slope;
	V_min		= in.V_min;
	t_min		= in.t_min;
	V_thr		= in.V_thr;
	t_thr		= in.t_thr;
	V_max		= in.V_max;
	t_max		= in.t_max;
	t_EAD		= in.t_EAD;
	V_EAD		= in.V_EAD;
	t_EAD2		= in.t_EAD2;
	V_EAD2		= in.V_EAD2;
	V_90		= in.V_90;
	t_90		= in.t_90;
	t_90_old	= in.t_90_old;
	dV_old		= in.dV_old;
	flag2		= in.flag2;
	flag_EAD	= in.flag_EAD;
	flag_EAD2	= in.flag_EAD2;

	MarkovIterations	 = in.MarkovIterations;
	MarkovIterationError = in.MarkovIterationError;

	I_total_old	= in.I_total_old;
	dI_total	= in.dI_total;
	dI_total_old= in.dI_total_old;
	t_I_total	= in.t_I_total;
	I_total_pt	= in.I_total_pt;
	V_I_total	= in.V_I_total;
	flagI_total	= in.flagI_total;
	t_90_old	= in.t_90_old;
	flag_EAD	= in.flag_EAD;
	flag_EAD2	= in.flag_EAD2;
	CV			= in.CV;
	CV_EAD		= in.CV_EAD;

	// memory allocated on heap
	// disabled because there is no need to copy these temporary variables right now.
	//copyTemps( in );
}


/*
 * copyTemps - copy contents of heap memory for Cell object
 */
void Moreno2011Cell::copyTemps(const Moreno2011Cell& source){
	gsl_vector_memcpy ( Yn  , source.Yn   );
	gsl_vector_memcpy (Ynph , source.Ynph );
	gsl_vector_memcpy (Ynp1 , source.Ynp1 );
	gsl_vector_memcpy (Ynp1i, source.Ynp1i);
	gsl_vector_memcpy (DiagonalAn , source.DiagonalAn );
	gsl_vector_memcpy (DiagonalAnInv , source.DiagonalAnInv );
	gsl_matrix_memcpy (OffDiagonalAn , source.OffDiagonalAn );
}

/*
 * setParams - set cell parameters
 *
 */
void Moreno2011Cell::setParams(const double* par){
	G_Na = par[0];
	GK1	 = par[1];

}

void Moreno2011Cell::NaChannelParameters::Initialize(cellmodel::Cells* env ){
	const Blockers b  = env->blocker;
	const double drug = env->drug;

	// Default
	drug_charged = 0;
	drug_neutral = 0;

	a4	= 0*a2;
	b4	= 0*a3;
	a5	= 0*a2;
	b5	= 0*a3;

	a44	= 0  *a2;
	b44	= 0  *a3;
	a55	= 0;
	b55	= 0;

	a_44= 0  *a2;
	b_44= 0  *a3;
	a_55= 0;
	b_55= 0;

	b13c= 0;
	b22	= 0;
	a_33= 0;
	b13n= 0;
	b_22= 0;
	bx2	= 0;

	kd0	= 0;
	kd_open= 0;
	// charged drug
	koff = 0;
	kcoff= 0;
	// charged drug
	kon	 = 0;
	kcon = 0;
	// neutral drug
	k_on	= 0;
	k_off	= 0;
	ki_on	= 0;
	ki_off	= 0;
	kc_on	= 0;
	kc_off	= 0;
	if (b == Lidocaine) {
		pKa		= 7.6;
		portion = 1/(1+ pow(10, (pH-pKa)) );
		diffusion = 500;

		drug_charged = drug * portion;
		drug_neutral = drug * (1 - portion);
		drug_distance = -0.7;
		kd0	= 318e-6;

		// charged drug
		kon	 = drug_charged*diffusion;
		kcon = kon;

		// neutral drug
		k_on	= drug_neutral*diffusion;
		k_off	= 400e-6 *diffusion;
		ki_on	= k_on/2;
		ki_off	= 3.4e-6 *diffusion;
		kc_on	= k_on/2;
		kc_off	= 900e-6 *diffusion;
	} else if (b == Flecainide) {
		pKa		= 9.3;
		portion = 1/(1+ pow(10, (pH-pKa)) );
		diffusion=5500;

		drug_charged = drug*portion;
		drug_neutral = drug*(1-portion);
		drug_distance = -0.7;
		kd0	= 11.2e-6;

		// charged drug
		kon	 = drug_charged*diffusion;
		kcon = kon;

		// neutral drug
		k_on	= drug_neutral*diffusion;
		k_off	= 400e-6 *diffusion;
		ki_on	= k_on/2;
		ki_off	= 5.4e-6 *diffusion;
		kc_on	= k_on/2;
		kc_off	= 800e-6 *diffusion;
	}
}


// Rate Constants **********************************************************
// WT Fits Reduced Model (no IM1, IM2)
void Moreno2011Cell::NaChannelParameters::ComputeVDepTerms(const double V, cellmodel::Cells* env ){
	const Blockers b  = env->blocker;
	const double drug = env->drug;
	const double Tfactor = 1.0/(pow(Q10, (37.0-(T-273))/10.0));

	// Default
	a11 = Tfactor*8.5539/(7.4392e-2*exp(-V/17.0)+ 2.0373e-1*exp(-V/150));
	a12 = Tfactor*8.5539/(7.4392e-2*exp(-V/15.0)+ 2.0373e-1*exp(-V/150));
	a13 = Tfactor*8.5539/(7.4392e-2*exp(-V/12.0)+ 2.0373e-1*exp(-V/150));
	b11 = Tfactor*7.5215e-2*exp(-V/20.3);
	b12 = Tfactor*2.7574*exp(-(V-5)/20.3);
	b13 = Tfactor*4.7755e-1*exp(-(V-10)/20.3);

	a2	= Tfactor*(13.370*exp(V/43.749));
	a3  = Tfactor*5.1458e-6*exp(-V/8.2471);
	b3  = Tfactor*6.1205*exp((V)/13.542);
	b2	= a13*a2*a3/(b13*b3);

	a4	= 0*a2;
	b4	= 0*a3;
	a5	= 0*a2;
	b5	= 0*a3;

	ax	= 3.4229e-2*a2;
	bx	= 1.7898e-2*a3;

	ax1	= 0;
	bx1	= 0;
	a13c= 0;
	a22	= 0;
	b33	= 0;
	a33	= 0;
	a55	= 0;
	b55	= 0;

	ax2	= 0;
	a13n= 0;
	a_22= 0;
	b_33= 0;
	a_55= 0;
	b_55= 0;


	if (b == Lidocaine) {
		ax1	= 6.3992e-07 * ax;
		bx1	= 1.3511e+00 * bx;
		a13c= 5.6974e-03 * a13;
		a22	= 6.7067e-06 * a2;
		b33	= 1.9698e-05 * b3;
		a33	= 3.2976e+00 * a3;

		a44	= 0 * a2;
		b44	= 0 * a3;

		ax2	= 1.3110e-01 * ax;
		a13n= 8.4559e+01 * a13;
		a_22= 1.7084e-05 * a2;
		b_33= 4.8477e+00 * b3;

		a_44= 0 * a2;
		b_44= 0 * a3;
	} else if (b == Flecainide) {
		ax1 = 5.7839e-05 * ax;
		bx1 = 1.6689e-08 * bx;
		a13c= 3.6324e-03 * a13;
		a22 = 1.4847e+03 * a2;
		b33 = 1.7352e-06 * b3;
		a33 = 6.7505e-05 * a3;

		a44 = 2.4135e+00 * a2;
		b44 = 4.9001e-02 * a3;

		ax2 = 2.6126e-01 * ax;
		a13n= 2.6452e+00 * a13;
		a_22= 4.2385e+01 * a2;
		b_33= 2.1181e+00 * b3;

		a_44= 1.0326e-03 * a2;
		b_44= 2.1378e-02 * a3;
	}

	kd_open= kd0*exp( (drug_distance*V*F) /(R*T));

	// charged drug
	koff	= kd_open*diffusion;
	kcoff	= koff;

	b13c = (drug ==0 || drug_charged ==0 ) ?  0 : (b13*kcon*koff*a13c)/(kon*kcoff*a13);
	b22  = (b13c ==0) ? 0 : (a13c*a22*a33)/(b13c*b33);

	a_33 = (drug ==0 || drug_neutral ==0 ) ? 0 : (ki_off*a3*kc_on*b_33)/(ki_on*kc_off*b3);
	b13n = (drug ==0 || drug_neutral ==0 ) ? 0 : (b13*kc_on*a13n*k_off)/(kc_off*a13*k_on);
	b_22 = (b13n==0) ? 0 : (a_33*a13n*a_22)/(b_33*b13n);
	bx2  = (drug ==0 || drug_neutral ==0)  ? 0 : (bx*k_on*ax2*ki_off)/(ax*ki_on*k_off);
}

/*
 * Get a pointer "back" up to the Cells environment so
 * that each cell knows the environment parameter values
 */
void Moreno2011Cell::setEnvironment( cellmodel::Cells* env ){
	cellenv = env;
	NaChParams.Initialize( cellenv );
}



void Moreno2011Cell::Calculate_All(const double time){
	//Updating current calculations for each cell
	CalculateCurrents(time);

	//Updating ionic concentrations
	Calculate_Na_in();
	Calculate_K_in ();
	Calculate_Ca_sr();
	Calculate_Ca_ss();
	Calculate_Ca_in();

	//Calculates t_min, V_min, t_max, V_max etc. etc.; Update calculates axial currents and updates the voltage
	Calculate_Points( time );
}



/*
 * Resets the cell characteristics -- it is the reponsibility of the
 * calling method to call this at the right time - at the end of an old
 * stimulus cycle so that the data will be reset for new cycle.
 */
void Moreno2011Cell::Calculate_Reset(double t) {
	if (cellenv == nullptr )
		return;

	if (t >= cellenv->waitTime ) {
		t_90_old = t_90;
		flag2		= false;
		flag_EAD	= false;
		flag_EAD2	= false;
		peak_slope	= 0;
		V_min	= Vmin_default;
		t_min	= 0;
		V_thr	= Vthr_default;
		t_thr	= 0;
		V_max	= Vthr_default;
		t_max	= t;
		t_EAD	= t;
		t_EAD2	= t;
		t_90	= t;
		V_EAD	= Vthr_default;
		V_EAD2	= Vthr_default;
		V_90	= Vthr_default;
		dV_old	= 0;
		flagI_total = false;
	}
}

/*
 * First-order Rush-Larsen
 */
void Moreno2011Cell::computeRushLarsenStepState(const double t, const double dt ) {
	semi_impl_I_Na(t, dt);	//Fast sodium current
	rush_I_Na_L(dt);

	calcI_K1();
	calcI_Kp();
	calcI_Na_Ca();
	calcI_Na_K();
	calcI_p_Ca();
	calcI_Ca_b();
	calcI_Na_b();


	rush_I_Ca_L(dt);	//L type calcium channel Ca contribution
	rush_I_Kr(dt);		//Rapid rectifier
	rush_I_Ks(dt);		//Slow rectifier
	rush_I_to(dt);		//Transient outward current

	calcI_tr ();
	calcI_leak();
	calcI_up();
	rush_I_rel(dt);		//Release from SR

	//Updating ionic concentrations
	rush_Na_in(dt);		//Dynamic [Na] in Myoplasm
	rush_K_in(dt);		//Dynamic [K] in Myoplasm
	rush_Ca_in(dt);		//Dynamic [Ca] in Myoplasm
	rush_Ca_sr(dt);		//Dynamic [Ca] in JSR
	rush_Ca_ss(dt);		//Dynamic [Ca] in NSR

	calcI_total(t);

	//Calculates t_min, V_min, t_max, V_max etc. etc.; Update calculates axial currents and updates the voltage
	Calculate_Points( t );
}

/*
 * Second-order Rush-Larsen - nor sure if this works
 *
 */
void Moreno2011Cell::computeRushLarsen2StepState(const double t, const double dt, const int step, const double Vhalf ) {
	double timestep = dt;
	if (step ==1){
		// Save the state before the Rush-Larsen half-step
		oldm	= m;		oldh	= h;
		oldj	= j;
		oldd	= d;		oldf	= f;
		oldf2	= f2;		oldf_Ca	= f_Ca;
		oldxr1	= xr1;		oldxr2	= xr2;
		oldxs	= xs;
		oldr	= r;		olds	= s;
		oldR_bar = R_bar;
		oldmL    = mL;		oldhL	= hL;
		oldCa_in = Ca_in;	oldCa_sr= Ca_sr;
		oldCa_ss = Ca_ss;
		oldNa_in = Na_in,	oldK_in	= K_in;

		exponential_I_Na1(t, dt);	//Fast sodium current
		// The first step should be a half-step and the state should not be changed.
		timestep = dt/2;
	} else {
		// Restore the states for the second (full) step
		m	= oldm;		h	= oldh;
		j	= oldj;
		d	= oldd;		f	= oldf;
		f2	= oldf2;	f_Ca= oldf_Ca;
		xr1	= oldxr1;	xr2	= oldxr2;
		xs	= oldxs;
		r	= oldr;		s	= olds;
		R_bar = oldR_bar;
		mL    = oldmL;		hL		= oldhL;
		Ca_in = oldCa_in;	Ca_sr	= oldCa_sr;
		Ca_ss = oldCa_ss;
		Na_in = oldNa_in,	K_in	= oldK_in;
		exponential_I_Na2(t, dt, Vhalf);	//Fast sodium current
	}

	rush_I_Na_L( timestep );
	calcI_K1();
	calcI_Kp();
	calcI_Na_Ca();
	calcI_Na_K();
	calcI_p_Ca();
	calcI_Ca_b();
	calcI_Na_b();

	rush_I_Ca_L( timestep );	//L type calcium channel Ca contribution
	rush_I_Kr  ( timestep );		//Rapid rectifier
	rush_I_Ks  ( timestep );		//Slow rectifier
	rush_I_to  ( timestep );		//Transient outward current
	calcI_tr ();
	calcI_leak();
	calcI_up();
	rush_I_rel( timestep );	//Release from SR
	//Updating ionic concentrations
	rush_Na_in( timestep );	//Dynamic [Na] in Myoplasm
	rush_K_in ( timestep );		//Dynamic [K] in Myoplasm
	rush_Ca_in( timestep );	//Dynamic [Ca] in Myoplasm
	rush_Ca_sr( timestep );	//Dynamic [Ca] in JSR
	rush_Ca_ss( timestep );	//Dynamic [Ca] in NSR
	calcI_total(t);

	if (step == 2){
		//Calculates t_min, V_min, t_max, V_max etc. etc.; Update calculates axial currents and updates the voltage
		Calculate_Points( t );
	}
}


void Moreno2011Cell::setI_stimulus( double stim ) {
	I_stim = stim;
}

void Moreno2011Cell::calcI_total(double t ) {

	I_total_old = I_total;
	dI_total_old = dI_total;

	I_Na_ion_total = I_Na	+ I_Na_L + I_Na_b	+ 3 * I_Na_K + 3 * I_Na_Ca;
	I_K_ion_total  = I_Kr	+ I_Ks	 + I_K1		+ I_Kp		 - 2 * I_Na_K	+ I_to;
	I_Ca_ion_total = I_Ca_L	+ I_p_Ca + I_Ca_b	- 2 * I_Na_Ca;

	I_total = I_stim + I_Na_ion_total + I_K_ion_total + I_Ca_ion_total; // used to compute dV/dt

	/*
	 * dI_total is used to detect peak inward current:
	 * i.e. to detect local minimum Itotal - where dI_total_old is < 0 and dI_total is > 0
	 */
	dI_total = I_total - I_total_old;
}

void Moreno2011Cell::calcI_Na(double t) {
	if (t==0)
		WT_SCN5A_Initial_Conditions();
	else
		if (cellenv->blocker == Lidocaine)
			WT_SCN5A_Lidocaine();
		else // if cellenv->blocker is Flecainide || None
			WT_SCN5A_Flecainide();
}

/*
 * Na current in heart failure cells
 */
void Moreno2011Cell::calcI_Na_L() {
	double G_Na_L = (cellenv->HeartFailure) ? 5.0 * 0.0065 : 0;
	const double tau_h_L = 600;

	const double E_Na = RTF * log(Na_out / Na_in);
	const double amL = 0.32 * (V + 47.13) / (1.0 - exp(-0.1 * (V + 47.13)));
	const double bmL = 0.08 * exp(-V / 11.0);
	const double hLss = 1.0 / (1 + exp((V + 91.0) / 6.1));

	const double tau_m_L = 1. / (amL + bmL);
	const double mLss = amL / (amL + bmL);

	dmL = (mLss - mL)/tau_m_L;
	dhL = (hLss - hL)/tau_h_L;

	I_Na_L = G_Na_L * mL * mL * mL * hL * (V - E_Na);
}

void Moreno2011Cell::calcI_Na_Ca() {
	const double k_Na_Ca = (cellenv->HeartFailure) ? 1.65 * 1000 : 1000; // pA/pF
	const double Km_Na = 87.5; // Na_in half-saturation concentration of NaCa exhanger (mM)
	const double KmNa3 = Km_Na * Km_Na * Km_Na; // KmNa^3 (mM^3)
	const double Km_Ca = 1.38; // Ca_in half-saturation concentration of NaCa exhanger (mM)
	const double k_Na_Ca_sat = 0.1;		// Saturation factor
	const double gamma = 0.35;// Position of energy barrier controlling voltage dependance of inaca
	const double alpha = 2.5;			// Factor enhancing outward nature of I_Na_Ca
	const double Nao3 = Na_out * Na_out * Na_out;
	const double Nai3 = Na_in  * Na_in  * Na_in;
	const double VmRTF = V / RTF;
	I_Na_Ca = k_Na_Ca
			* (1. / (KmNa3 + Nao3))
			* (1. / (Km_Ca + Ca_out))
			* (1. / (1 + k_Na_Ca_sat * exp((gamma - 1) * VmRTF )))
			* (exp(gamma * VmRTF ) * Nai3 * Ca_out - exp((gamma - 1) * VmRTF ) * Nao3 * Ca_in * alpha);
}

void Moreno2011Cell::calcI_Ca_L() {
	//const double G_Ca_L = 3.980E-5;		// cm/(ms*uF)
	const double G_Ca_L = (cellenv->HeartFailure) ? 2.5 * 3.980E-5 : 3.980E-5;

	const double d_inf = 1. / (1. + exp((-8 - V) / 7.5));
	const double a_d = 1.4 / (1. + exp((-35 - V) / 13)) + 0.25;
	const double b_d = 1.4 / (1. + exp((V + 5) / 5));
	const double c_d = 1. / (1. + exp((50 - V) / 20));
	const double tau_d = a_d * b_d + c_d;

	// d d/dt
	dd = (d_inf - d)/tau_d;

	const double f_inf = 1. / (1. + exp((V + 20) / 7));
	const double a_f = 1102.5 * exp(-(V + 27) * (V + 27) / 225);
	const double b_f = 200. / (1 + exp((13 - V) / 10.));
	const double c_f = (180. / (1 + exp((V + 30) / 10))) + 20;
	const double tau_f = a_f + b_f + c_f;

	// df/dt
	df = (f_inf - f)/tau_f;

	const double f2_inf = 0.67 / (1. + exp((V + 35) / 7)) + 0.33;
	const double a_f2 = 600 * exp(-(V + 25) * (V + 25) / 170);
	const double b_f2 = 31 / (1. + exp((25 - V) / 10));
	const double c_f2 = 16 / (1. + exp((V + 30) / 10));
	const double tau_f2 = a_f2 + b_f2 + c_f2;

	// df2/dt
	df2 = (f2_inf - f2)/tau_f2;

	const double f_Ca_inf = 0.6 / (1 + (Ca_ss / 0.05) * (Ca_ss / 0.05)) + 0.4;
	const double tau_f_Ca = 80. / (1 + (Ca_ss / 0.05) * (Ca_ss / 0.05)) + 2.;

	// dfCa/dt
	df_Ca = (f_Ca_inf - f_Ca)/tau_f_Ca;

	const double Vm15RT2F = z_Ca * (V - 15)/RTF;
	// Need to apply L'Hopital's rule if V is close to 15
	I_Ca_L = G_Ca_L * F * d * f * f2 * f_Ca * 2 * Vm15RT2F
			* (0.25 * exp( Vm15RT2F ) * Ca_ss - Ca_out)
			/ (exp( Vm15RT2F ) - 1.);
}

void Moreno2011Cell::calcI_Kr() {
	const double G_Kr = 0.153;		// nS/pF
	const double E_Kr = RTF * log(K_out / K_in);

	const double a_xr1 = 450. / (1. + exp((-45. - V) / 10.));
	const double b_xr1 = 6. / (1. + exp((V - (-30.)) / 11.5));
	const double xr1_ss = 1. / (1. + exp((-26. - V) / 7.));
	const double tau_xr1 = a_xr1 * b_xr1;

	// dxr1/dt
	dxr1 = (xr1_ss - xr1)/tau_xr1;

	const double a_xr2 = 3. / (1. + exp((-60. - V) / 20.));
	const double b_xr2 = 1.12 / (1. + exp((V - 60.) / 20.));
	const double xr2_ss = 1. / (1. + exp((V - (-88.)) / 24.));
	const double tau_xr2 = a_xr2 * b_xr2;

	// dxr2/dt
	dxr2 = (xr2_ss - xr2)/tau_xr2;

	I_Kr = G_Kr * sqrt(K_out / 5.4) * xr1 * xr2 * (V - E_Kr);
}

void Moreno2011Cell::calcI_Ks() {
	double G_Ks = 0.392;
	switch (Cell_type) {
	case 1:	G_Ks = 0.392;	break;	// Cell type == 1 for endo
	case 2:	G_Ks = 0.098;	break;	// Cell type == 2 for M
	case 3:	G_Ks = 0.392;	break;	// Cell type == 3 for epi
	default: ;
	}

	const double PR_NaK = 0.03;
	const double E_Ks = RTF * log((K_out + PR_NaK * Na_out) / (K_in + PR_NaK * Na_in));

	const double xs_ss = 1. / (1. + exp((-5. - V) / 14.));
	const double ax_s = (1400. / (sqrt(1. + exp((5. - V) / 6))));
	const double bx_s = (1. / (1. + exp((V - 35.) / 15.)));
	const double tau_xs = (ax_s * bx_s) + 80;

	// dxs/dt
	dxs = (xs_ss - xs)/tau_xs;
	I_Ks = G_Ks * xs * xs * (V - E_Ks);
}

void Moreno2011Cell::calcI_K1() {
	//const double G_K1 = 5.405;			// nS/pF
	const double G_K1 = (cellenv->HeartFailure) ? 0.75 * GK1 : GK1;

	const double E_K1 = RTF * log(K_out / K_in);
	const double VmEK1 = V - E_K1;
	const double a_K1 = 0.1 / (1. + exp(0.06 * (VmEK1 - 200)));
	const double b_K1 = (3. * exp(0.0002 * (VmEK1 + 100)) + exp(0.1 * (VmEK1 - 10)))
			/ (1. + exp(-0.5 * VmEK1));

	const double K1_s = a_K1 / (a_K1 + b_K1);

	I_K1 = G_K1 * K1_s * (V - E_K1);
}

void Moreno2011Cell::calcI_to() {
	const double endoGtomax  = 0.073;
	const double McellGtomax = 0.294;
	const double epiGtomax   = 0.294;
	double G_to;
	// default values for s_inf and tau_s:
	double s_inf = 1. / (1. + exp((V + 20) / 5.));
	double tau_s = 85. * exp(-(V + 45.) * (V + 45.) / 320.)
			+ 5. / (1. + exp((V - 20.) / 5.)) + 3.;

	switch (Cell_type) {
	case 1:				// Cell type = (1) for endocardial
		G_to = (cellenv->HeartFailure) ? 0.64 * endoGtomax : endoGtomax;		// nS/pF
		s_inf = 1. / (1. + exp((V + 28) / 5.));
		tau_s = 1000. * exp(-(V + 67) * (V + 67) / 1000.) + 8.;
		break;
	case 2:				// Cell type = (2) for M cells
		G_to = (cellenv->HeartFailure) ? 0.64 * McellGtomax : McellGtomax;	// nS/pF
		// default values for s_inf and tau_s
		break;
	case 3:				// Cell type = (3) for epicardial
	default:
		G_to = (cellenv->HeartFailure) ? 0.64 * epiGtomax : epiGtomax;		// nS/pF
	}

	const double E_to = ((R * T) / (z_K * F)) * log(K_out / K_in);

	const double r_inf = 1. / (1. + exp((20 - V) / 6.));
	const double tau_r = 9.5 * exp(-(V + 40.) * (V + 40.) / 1800.) + 0.8;

	// dr/dt
	dr = (r_inf - r)/tau_r;
	ds = (s_inf - s)/tau_s;

	I_to = G_to * r * s * (V - E_to);
}

void Moreno2011Cell::calcI_p_Ca() {
	const double G_p_Ca = 0.1238;		// nS/pF
	const double Km_p_Ca = 0.0005;		// Half saturation constant (mM)

	I_p_Ca = G_p_Ca * (Ca_in / (Km_p_Ca + Ca_in));
}

void Moreno2011Cell::calcI_Kp() {
	const double G_Kp = 0.0146;			// nS/pF
	const double E_Kp = RTF * log(K_out / K_in);

	I_Kp = G_Kp * (1. / (1. + exp((25 - V) / 5.98))) * (V - E_Kp);
}

void Moreno2011Cell::calcI_Ca_b() {
	const double G_Ca_b = 0.000592;
	const double E_Ca_b = (RTF/z_Ca) * log(Ca_out / Ca_in);
	I_Ca_b = G_Ca_b * (V - E_Ca_b);
}

void Moreno2011Cell::calcI_Na_b() {
	const double G_Na_b = 0.000290;		// nS/pF

	double E_Na_b = ((R * T) / (z_Na * F)) * log(Na_out / Na_in);
	I_Na_b = G_Na_b * (V - E_Na_b);
}

void Moreno2011Cell::calcI_Na_K() {
	//const double I_Na_K_bar = 2.724;	// Maximal I_Na_K (pA/pF)
	const double I_Na_K_bar = (cellenv->HeartFailure) ? 1.0 * 2.724 : 2.724;
	const double Km_Na_in = 40;			// Na_in half saturation constant
	const double Km_K_out = 1;			// K_out half saturation constant

	I_Na_K = I_Na_K_bar	* (K_out * Na_in / ((K_out + Km_K_out) * (Na_in + Km_Na_in)))
			/ (1. + 0.1245 * exp( -0.1 * V /RTF ) + 0.0365 * exp(-V / RTF) );
}

void Moreno2011Cell::calcI_rel() {
	double k_Ca_sr, k1_rel, k2_rel;
	const double G_rel   = 0.102;	// Max. rate constant of Ca release from JSR due to overload (mM/ms)
	const double k1_rel_ = 0.15;// R to O and RI to I I_rel transition rate (mM^2/ms)
	const double k2_rel_ = 0.045;// O to I and R to RI I_rel transition rate (mM^2/ms)
	const double k3_rel  = 0.060;// O to R and I to RI I_rel transition rate (ms^-1)
	const double k4_rel  = 0.005;// I to O and RI to I I_rel transition rate (ms^-1)
	const double EC_sr  = 1.5;		// Ca_sr half-saturation constant of k_Ca_sr
	const double max_sr = 2.5;			// Max value of k_Ca_sr (dimensionless)
	const double min_sr = 1;			// Min value of k_Ca_sr (dimensionless)

	k_Ca_sr = max_sr - ((max_sr - min_sr) / (1 + (EC_sr / Ca_sr) * (EC_sr / Ca_sr)));

	k1_rel = k1_rel_ / k_Ca_sr;
	k2_rel = k2_rel_ * k_Ca_sr;

	dR_bar = -k2_rel * Ca_ss * R_bar + k4_rel * (1 - R_bar);
	const double OO = (k1_rel * Ca_ss * Ca_ss * R_bar) / (k3_rel + k1_rel * Ca_ss * Ca_ss);

	I_rel = G_rel * OO * (Ca_sr - Ca_ss);
}

void Moreno2011Cell::calcI_up() {
	//const double G_up = 0.006375;		// mM/ms
	const double G_up = (cellenv->HeartFailure) ? 0.64 * 0.006375 :  0.006375;
	const double Km_up = 0.00025;			//mM
	const double Cai2 = Ca_in * Ca_in;
	I_up = G_up * Cai2 / (Cai2 + (Km_up * Km_up) );
}

void Moreno2011Cell::calcI_leak() {
	//const double G_leak = 0.00036;		// mM/ms
	const double G_leak = (cellenv->HeartFailure) ? 0.00025 : 0.00036;
	I_leak = G_leak * (Ca_sr - Ca_in);
}

void Moreno2011Cell::calcI_tr() {
	const double G_tr = 0.0038;		// mM/ms

	I_tr = G_tr * (Ca_ss - Ca_in);
}

/************************************** Dynamic Ion concentrations ************/

//New K concentration in myoplasm
void Moreno2011Cell::Calculate_K_in() {
	dK_in = -((I_K_ion_total + I_stim) * CAP) / (V_cyto * z_K * F);
}

//New Na concentration in myoplasm
void Moreno2011Cell::Calculate_Na_in() {
	dNa_in = -I_Na_ion_total * CAP / (V_cyto * z_Na * F);
}

/*
 * New version - from cellml.org - see:
 * https://models.cellml.org/exposure/a7179d94365ff0c9c0e6eb7c6a787d3d/ten_tusscher_model_2006_IK1Ko_epi_units.cellml/view
 */
void Moreno2011Cell::Calculate_Ca_in() {
	const double K_bufc = 0.001; // Ca_in half-saturation constant for cytoplasmic buffer (mM)
	const double Bufc   = 0.2;   //Total cytoplasmic buffer concentration (mM)
	const double cai_p_kbufc = Ca_in + K_bufc;
	const double denominator = 1.0 + Bufc * K_bufc /( cai_p_kbufc * cai_p_kbufc );
	const double ca_i_bufc = 1.0/denominator;

	dCa_in = (-(I_Ca_b + I_p_Ca - 2 * I_Na_Ca) * CAP / (V_cyto * z_Ca * F))
						+ (V_sr / V_cyto) * (I_leak - I_up) + I_tr;
	dCa_in = ca_i_bufc * dCa_in;
}

/*
 * New version - from cellml.org - see:
 * https://models.cellml.org/exposure/a7179d94365ff0c9c0e6eb7c6a787d3d/ten_tusscher_model_2006_IK1Ko_epi_units.cellml/view
 */
void Moreno2011Cell::Calculate_Ca_sr() {
	const double Bufsr   = 10;	//Total SR buffer concentration (mM)
	const double K_bufsr = 0.3;	// Ca_sr half-saturation constant for SR buffer (mM)
	const double casr_p_kbufsr = Ca_sr + K_bufsr;
	const double denominator = 1.0 + Bufsr * K_bufsr /( casr_p_kbufsr * casr_p_kbufsr );
	const double ca_sr_bufsr = 1.0/denominator;

	dCa_sr = I_up - I_rel - I_leak;
	dCa_sr = ca_sr_bufsr * dCa_sr;
}

/*
 * New version - from cellml.org - see:
 * https://models.cellml.org/exposure/a7179d94365ff0c9c0e6eb7c6a787d3d/ten_tusscher_model_2006_IK1Ko_epi_units.cellml/view
 */
void Moreno2011Cell::Calculate_Ca_ss() {
	const double Bufss	 = 0.4;		//Total SS buffer concentration (mM)
	const double K_bufss = 0.00025;	// Ca_ss half-saturation constant for SR buffer (mM)
	const double cass_p_kbufss = Ca_ss + K_bufss;
	const double denominator = 1.0 + Bufss * K_bufss /( cass_p_kbufss * cass_p_kbufss );
	const double ca_ss_bufss = 1.0/denominator;

	dCa_ss = ((-I_Ca_L * CAP) / (V_ss * z_Ca * F) + (V_sr / V_ss) * I_rel - (V_cyto / V_ss) * I_tr);
	dCa_ss = ca_ss_bufss * dCa_ss;
}

void Moreno2011Cell::Calculate_Points(double t) {
	if (t <= cellenv->waitTime )
			return; // nothing to report yet

	const double cycle_time = fmod ( t, cellenv->BCL );

	//Minimum voltage catcher
	if (cycle_time < 0.01) {
		t_min = t;
		V_min = V;
	}

	//Threshold voltage point catcher needed for t_min
	if ( dVdt > peak_slope) {
		peak_slope = dVdt;
		V_thr = V;
		t_thr = t;
	}

	//Maximum voltage point catcher
	if ( dV_old > 0 && dVdt < 0 && V > -5 ) {
		V_max = V;
		t_max = t;
	}

	if (cycle_time > 200) {
		//I_total minimum point
		if (dI_total_old < 0 && dI_total > 0 && !flagI_total) {
			flagI_total = true;
			t_I_total = t;
			I_total_pt = I_total;
			V_I_total = V;
		}

		//Local minima
		if ( dV_old < 0 && dVdt > 0 && V > -40 && !flag_EAD ) {
			flag_EAD = true;
			V_EAD = V;
			t_EAD = t;
		}

		//Local maxima
		if (t > t_EAD && dV_old > 0 && dVdt < 0 && !flag_EAD2 ) {
			flag_EAD2 = true;
			V_EAD2 = V;
			t_EAD2 = t;
		}

		//V_90 voltage point catcher, cycle number and APD_90 calculation
		V_90 = V_max - .90 * (V_max - V_min);

		//Diastolic Interval
		DI = t_thr - t_90_old;

		//Height (or length) of the EAD
		L_EAD = V_EAD2 - V_EAD;
	}

	if ( fabs(V - V_90) < 0.015 && dVdt < 0 && cycle_time > 30
			&& t > t_max && !flag2) {
		flag2 = true;
		t_90 = t;
		APD_90 = t_90 - t_thr;
	}

	//Put in two conditions here for the check to see if there was an AP initiated.

}



void Moreno2011Cell::SetInitialConditions() {
	Cell_type = 3;
	V		= -86.2;
	dVdt	= 0;
	Na_in	= 7.67;
	K_in	= 138.3;
	Ca_in	= 0.00007;
	Ca_sr	= 1.3;
	Ca_ss	= 0.00007;
	m	= 0;
	h	= 0.75;
	j	= 0.75;

	mL	= 0.00111859;
	hL	= 0.339310414;

	d	= 0;
	f	= 1;
	f2	= 1;
	f_Ca = 1;
	r	= 0;
	s	= 1;
	xr1 = 0;
	xr2 = 1;
	xs	= 0;
	R_bar = 1;
	WT_SCN5A_Initial_Conditions();

	for (unsigned int i=0; i<naux; i++)
		currents[i] = 0;

	Calculate_Reset( 0 );
}


double Moreno2011Cell::CalculateCurrents(double t) {
	calcI_Na(t);
	calcI_Ca_L();
	calcI_Kr();
	calcI_Ks();
	calcI_K1();
	calcI_Kp();
	calcI_to();
	calcI_Na_Ca();
	calcI_Na_K();
	calcI_p_Ca();
	calcI_Ca_b();
	calcI_Na_b();

	calcI_tr ();
	calcI_leak();
	calcI_up();
	calcI_rel();

	calcI_total(t);
	return I_total;
}

/*
 * For finite difference of dVdt equation
 */
double Moreno2011Cell::CalculateCurrents(double t, const double Vx) {
	const double oldV = V;
	V = Vx;
	calcI_Na(t);
	calcI_Ca_L();
	calcI_Kr();
	calcI_Ks();
	calcI_K1();
	calcI_Kp();
	calcI_to();
	calcI_Na_Ca();
	calcI_Na_K();
	calcI_p_Ca();
	calcI_Ca_b();
	calcI_Na_b();

	calcI_tr ();
	calcI_leak();
	calcI_up();
	calcI_rel();

	calcI_total(t);
	// Reset V
	V = oldV;
	return I_total;
}
/******************************************************************************
 WT Markov Model with Zheng's new drug binding states
 January 14, 2009
 ******************************************************************************/
void Moreno2011Cell::WT_SCN5A_Initial_Conditions() {

	if (initializedExternally)
		return;

	/* Initial Conditions for Sodium channel states:  Set all the states to 0 except drug-free C3  */
	for (unsigned int i=nNaChFirstState; i<=nNaChLastState; i++)
		state[i] = 0;
	C3 = 1; // all channels in closed state C3

	// The above statement does not work! if the references are messed up.
	// Intermittent error: although the state array is in the heap,
	//   for some reason, within the same executable, the C3 variable moves
	//   into the stack!
	// Fix:
	//state[5] = 1;
}

void Moreno2011Cell::WT_SCN5A_Lidocaine(){
	checkNaChannelStates();

	const double pKa  = 7.6;
	const double portion = 1/(1+ pow(10, (pH-pKa)) );
	const double diffusion = 500;

	const double drug_charged = cellenv->drug * portion;
	const double drug_neutral = cellenv->drug * (1 - portion);
	const double drug_distance = -0.7;

	const double E_Na = RTF*log(Na_out/Na_in);

	//Rate Constants **********************************************************
	//WT Fits Reduced Model (no IM1, IM2)
	const double Tfactor = 1.0/(pow(Q10, (37.0-(T-273))/10.0));

	const double a11 = Tfactor*8.5539/(7.4392e-2*exp(-V/17.0)+ 2.0373e-1*exp(-V/150));
	const double a12 = Tfactor*8.5539/(7.4392e-2*exp(-V/15.0)+ 2.0373e-1*exp(-V/150));
	const double a13 = Tfactor*8.5539/(7.4392e-2*exp(-V/12.0)+ 2.0373e-1*exp(-V/150));
	const double b11 = Tfactor*7.5215e-2*exp(-V/20.3);
	const double b12 = Tfactor*2.7574*exp(-(V-5)/20.3);
	const double b13 = Tfactor*4.7755e-1*exp(-(V-10)/20.3);

	const double a3  = Tfactor*5.1458e-6*exp(-V/8.2471);
	const double b3  = Tfactor*6.1205*exp((V)/13.542);

	const double a2= Tfactor*(13.370*exp(V/43.749));
	const double b2= ((a13*a2*a3)/(b13*b3));

	const double a4 = 0*a2;
	const double b4 = 0*a3;
	const double a5 = 0*a2;
	const double b5 = 0*a3;

	const double ax	= 3.4229e-2*a2;
	const double bx	= 1.7898e-2*a3;


	const double ax1	= 6.3992e-07  *ax;
	const double bx1	= 1.3511e+00  *bx;
	const double a13c	= 5.6974e-03  *a13;
	const double a22	= 6.7067e-06  *a2;
	const double b33	=  1.9698e-05 *b3;
	const double a33	= 3.2976e+00  *a3;

	const double a44	= 0  *a2;
	const double b44	= 0  *a3;
	const double a55	= 0;
	const double b55	= 0;

	const double ax2	= 1.3110e-01  *ax;
	const double a13n	= 8.4559e+01  *a13;
	const double a_22	= 1.7084e-05  *a2;
	const double b_33	= 4.8477e+00  *b3;

	const double a_44	= 0  *a2;
	const double b_44	= 0  *a3;
	const double a_55	= 0;
	const double b_55	= 0;


	const double kd0	= 318 * 1e-6;
	const double kd_open= kd0*exp( drug_distance*V / RTF );

	// charged drug
	const double kon	= drug_charged*diffusion;
	const double koff	= kd_open*diffusion;
	const double kcoff	= koff;
	const double kcon	= kon;

	const double b13c	= (cellenv->drug ==0 || drug_charged ==0 ) ? 0 : (b13*kcon*koff*a13c)/(kon*kcoff*a13);
	const double b22	= (b13c ==0) ? 0 : (a13c*a22*a33)/(b13c*b33);

	// neutral drug
	const double k_on	= drug_neutral*diffusion;
	const double k_off	= 400* 1e-6 *diffusion;
	const double ki_on	= k_on/2;
	const double ki_off	= 3.4* 1e-6 *diffusion;
	const double kc_on	= k_on/2;
	const double kc_off	= 900* 1e-6 *diffusion;

	const double a_33	= (cellenv->drug ==0 || drug_neutral ==0 ) ? 0 : (ki_off*a3*kc_on*b_33)/(ki_on*kc_off*b3);
	const double b13n	= (cellenv->drug ==0 || drug_neutral ==0 ) ? 0 : (b13*kc_on*a13n*k_off)/(kc_off*a13*k_on);
	const double b_22	= (b13n ==0) ? 0 : (a_33*a13n*a_22)/(b_33*b13n);
	const double bx2	= (cellenv->drug ==0 || drug_neutral ==0) ? 0 : (bx*k_on*ax2*ki_off)/(ax*ki_on*k_off);

	const double coef_O = (b13 + a2  + kon + k_on + ax);
	const double coef_C1 = (b12 + b3 + a13  + kcon + kc_on);
	const double coef_C2 = (b11 + b3 + a12  + kcon + kc_on);
	const double coef_C3 = (b3 + a11  + kcon + kc_on);
	const double coef_IC3 = (a11 + a3 + ki_on);
	const double coef_IC2 = (b11 + a3 + a12 + ki_on);
	const double coef_IF = (b12 + b2 + a3 + a4 + ki_on);
	const double coef_IM1 = (b4 + a5);
	const double coef_IM2 = b5;
	const double coef_OS = (bx + ki_on);

	const double coef_DO = (koff + b13c + a22 + ax1 );
	const double coef_DC1 = (kcoff + b12 + b33 + a13c );
	const double coef_DC2 = (kcoff + b11 + b33 + a12 );
	const double coef_DC3 = (kcoff+ b33 + a11 );
	const double coef_DOS = bx1;
	const double coef_DIC3 = (a11 + a33);
	const double coef_DIC2 = (a33 + b11 + a12);
	const double coef_DIF = (a33 + b12 + a44 + b22);
	const double coef_DIM1 = ( b44 + a55 );
	const double coef_DIM2 = b55 ;

	const double coef_D_O = (k_off + b13n + a_22 + ax2 );
	const double coef_D_C1 = (kc_off + b12 + b_33 + a13n );
	const double coef_D_C2 = (kc_off + b11 + b_33 + a12 );
	const double coef_D_C3 = (kc_off + b_33 + a11 );
	const double coef_D_OS = (bx2 + ki_off);
	const double coef_D_IC3 = (a_33 + a11 + ki_off);
	const double coef_D_IC2 = (a_33 + b11 + a12 + ki_off);
	const double coef_D_IF = (a_33 + a_44 + b_22 + b12 + ki_off);
	const double coef_D_IM1 = (b_44 + a_55);
	const double coef_D_IM2 = b_55;

	//Drug Free States
	dO   =  a13 * C1  + b2  * IF  + koff * DO  + k_off * D_O + bx * OS  - O * coef_O;
	dC1  =  a12 * C2  + a3  * IF  + b13  * O   + kcoff * DC1 + kc_off * D_C1 - C1 * coef_C1;
	dC2  =  a11 * C3  + a3  * IC2 + b12  * C1  + kcoff * DC2 + kc_off * D_C2 - C2 * coef_C2;
	dC3  =  a3  * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 - C3 * coef_C3;
	dIC3 =  b3  * C3  + b11 * IC2 + ki_off * D_IC3 - IC3 * coef_IC3;
	dIC2 =  a11 * IC3 + b3  * C2  + b12  * IF  + ki_off * D_IC2 - IC2 * coef_IC2;
	dIF  =  a12 * IC2 + b3  * C1  + b4   * IM1 + a2 * O + ki_off * D_IF - IF * coef_IF;
	dIM1 =  a4  * IF  + b5  * IM2 - IM1  * coef_IM1;
	dIM2 =  a5  * IM1 - IM2 * coef_IM2;
	dOS  =  ax  * O   + ki_off * D_OS - OS * coef_OS;

	//Charged Drug Bound States
	dDO   =  kon * O + a13c * DC1 + bx1 * DOS + b22 * DIF  - DO * coef_DO;
	dDC1  =  kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO  - DC1 * coef_DC1;
	dDC2  =  kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1  - DC2 * coef_DC2;
	dDC3  =  kcon * C3 + b11 * DC2 + a33 * DIC3  - DC3 * coef_DC3;
	dDOS  =  ax1 * DO -  DOS * coef_DOS;
	dDIC3 =  b33 * DC3 + b11 * DIC2 - DIC3 * coef_DIC3;
	dDIC2 =  b33 * DC2 + a11 * DIC3 + b12 * DIF - DIC2 * coef_DIC2;
	dDIF  =  b33 * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO - DIF * coef_DIF;
	dDIM1 =  a44 * DIF + b55 * DIM2 - DIM1 * coef_DIM1;
	dDIM2 =  a55 * DIM1 - DIM2 * coef_DIM2;

	//Neutral Drug Bound States
	dD_O   = k_on  * O     + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS  - D_O * coef_D_O;
	dD_C1  = kc_on * C1    + a12  * D_C2 + a_33 * D_IF + b13n * D_O  - D_C1 * coef_D_C1;
	dD_C2  = kc_on * C2    + a11  * D_C3 + a_33 * D_IC2 + b12 * D_C1  - D_C2 * coef_D_C2;
	dD_C3  = kc_on * C3    + a_33 * D_IC3 + b11 * D_C2  - D_C3 * coef_D_C3;
	dD_OS  = ax2   * D_O   + ki_on * OS - D_OS * coef_D_OS;
	dD_IC3 = b_33  * D_C3  + b11  * D_IC2 + ki_on * IC3 - D_IC3 * coef_D_IC3;
	dD_IC2 = b_33  * D_C2  + a11  * D_IC3 + b12 * D_IF + ki_on * IC2 - D_IC2 * coef_D_IC2;
	dD_IF  = b_33  * D_C1  + a12  * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF - D_IF * coef_D_IF;
	dD_IM1 = a_44  * D_IF  + b_55 * D_IM2 - D_IM1 * coef_D_IM1;
	dD_IM2 = a_55  * D_IM1 - D_IM2 * coef_D_IM2;


	I_Na = G_Na* O *(V - E_Na);
}

/******************************************************************************
 WT Markov Model with Zheng's new drug binding states
 January 14, 2009
 ******************************************************************************/
void Moreno2011Cell::WT_SCN5A_Flecainide(){
	checkNaChannelStates();

	const double pKa=9.3;
	const double portion = 1/(1+ pow(10, (pH-pKa)) );
	const double diffusion=5500;

	const double drug_charged=cellenv->drug*portion;
	const double drug_neutral=cellenv->drug*(1-portion);
	const double drug_distance= -0.7;

	const double E_Na = RTF*log(Na_out/Na_in);


	//Rate Constants **********************************************************
	//WT Fits Reduced Model (no IM1, IM2)
	const double Tfactor = 1.0/(pow(Q10, (37.0-(T-273))/10.0));

	const double a11= Tfactor*8.5539/(7.4392e-2*exp(-V/17.0)+ 2.0373e-1*exp(-V/150));
	const double a12= Tfactor*8.5539/(7.4392e-2*exp(-V/15.0)+ 2.0373e-1*exp(-V/150));
	const double a13= Tfactor*8.5539/(7.4392e-2*exp(-V/12.0)+ 2.0373e-1*exp(-V/150));
	const double b11= Tfactor*7.5215e-2*exp(-V/20.3);
	const double b12= Tfactor*2.7574*exp(-(V-5)/20.3);
	const double b13= Tfactor*4.7755e-1*exp(-(V-10)/20.3);

	const double a3 = Tfactor*5.1458e-6*exp(-V/8.2471);
	const double b3=Tfactor*6.1205*exp((V)/13.542);


	const double a2= Tfactor*(13.370*exp(V/43.749));
	const double b2= ((a13*a2*a3)/(b13*b3));

	const double a4 = 0*a2;
	const double b4 = 0*a3;
	const double a5= 0*a2;
	const double b5 = 0*a3;

	const double ax = 3.4229e-2*a2;
	const double bx = 1.7898e-2*a3;

	const double ax1 =5.7839e-05 * ax;
	const double bx1 =  1.6689e-08* bx;
	const double a13c = 3.6324e-03 *a13;
	const double a22 = 1.4847e+03 *a2;
	const double b33 =  1.7352e-06* b3;
	const double a33 = 6.7505e-05 * a3;
	const double a44 =  2.4135e+00* a2;
	const double b44 =  4.9001e-02* a3;
	const double a55 = 0;
	const double b55 = 0;

	const double ax2 = 2.6126e-01 * ax;
	const double a13n = 2.6452e+00 * a13;
	const double a_22 =  4.2385e+01 * a2;
	const double b_33 = 2.1181e+00 * b3;
	const double a_44 =  1.0326e-03 * a2;
	const double b_44 = 2.1378e-02 * a3;

	const double a_55 = 0;
	const double b_55 = 0;

	const double kd0	= 11.2*(1e-6);
	const double kd_open= kd0*exp( (drug_distance*V*F) /(R*T));

	// charged drug
	const double kon	= drug_charged*diffusion;
	const double koff	= kd_open*diffusion;
	const double kcoff	= koff;
	const double kcon	= kon;

	const double b13c = (cellenv->drug ==0 || drug_charged ==0 ) ?  0 : (b13*kcon*koff*a13c)/(kon*kcoff*a13);
	const double b22  = (b13c ==0) ? 0 : (a13c*a22*a33)/(b13c*b33);

	// neutral drug
	const double k_on	= drug_neutral*diffusion;
	const double k_off	= 400* 1e-6 *diffusion;
	const double ki_on	= k_on/2;
	const double ki_off	= 5.4* 1e-6 *diffusion;
	const double kc_on	= k_on/2;
	const double kc_off	= 800* 1e-6 *diffusion;

	const double a_33 = (cellenv->drug ==0 || drug_neutral ==0 ) ? 0 : (ki_off*a3*kc_on*b_33)/(ki_on*kc_off*b3);
	const double b13n = (cellenv->drug ==0 || drug_neutral ==0 ) ? 0 : (b13*kc_on*a13n*k_off)/(kc_off*a13*k_on);
	const double b_22 = (b13n==0) ? 0 : (a_33*a13n*a_22)/(b_33*b13n);
	const double bx2  = (cellenv->drug ==0 || drug_neutral ==0)  ? 0 : (bx*k_on*ax2*ki_off)/(ax*ki_on*k_off);

	const double coef_O = (b13 + a2  + kon + k_on + ax);
	const double coef_C1 = (b12 + b3 + a13  + kcon + kc_on);
	const double coef_C2 = (b11 + b3 + a12  + kcon + kc_on);
	const double coef_C3 = (b3 + a11  + kcon + kc_on);
	const double coef_IC3 = (a11 + a3 + ki_on);
	const double coef_IC2 = (b11 + a3 + a12 + ki_on);
	const double coef_IF = (b12 + b2 + a3 + a4 + ki_on);
	const double coef_IM1 = (b4 + a5);
	const double coef_IM2 = b5;
	const double coef_OS = (bx + ki_on);

	const double coef_DO = (koff + b13c + a22 + ax1 );
	const double coef_DC1 = (kcoff + b12 + b33 + a13c );
	const double coef_DC2 = (kcoff + b11 + b33 + a12 );
	const double coef_DC3 = (kcoff+ b33 + a11 );
	const double coef_DOS = bx1;
	const double coef_DIC3 = (a11 + a33);
	const double coef_DIC2 = (a33 + b11 + a12);
	const double coef_DIF = (a33 + b12 + a44 + b22);
	const double coef_DIM1 = ( b44 + a55 );
	const double coef_DIM2 = b55 ;

	const double coef_D_O = (k_off + b13n + a_22 + ax2 );
	const double coef_D_C1 = (kc_off + b12 + b_33 + a13n );
	const double coef_D_C2 = (kc_off + b11 + b_33 + a12 );
	const double coef_D_C3 = (kc_off + b_33 + a11 );
	const double coef_D_OS = (bx2 + ki_off);
	const double coef_D_IC3 = (a_33 + a11 + ki_off);
	const double coef_D_IC2 = (a_33 + b11 + a12 + ki_off);
	const double coef_D_IF = (a_33 + a_44 + b_22 + b12 + ki_off);
	const double coef_D_IM1 = (b_44 + a_55);
	const double coef_D_IM2 = b_55;

	//Drug Free States
	dO	 = a13 * C1 + b2 * IF + koff * DO + k_off * D_O + bx * OS  - O * coef_O;
	dC1	 = a12 * C2 + a3 * IF + b13 * O  + kcoff * DC1 + kc_off * D_C1 - C1 * coef_C1;
	dC2	 = a11 * C3 + a3 * IC2 + b12 * C1  + kcoff * DC2 + kc_off * D_C2 - C2 * coef_C2;
	dC3	 =  a3 * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 - C3 * coef_C3;
	dIC3 =  b3 * C3 + b11 * IC2 + ki_off * D_IC3 - IC3 * coef_IC3;
	dIC2 = a11 * IC3 + b3 * C2 + b12 * IF + ki_off * D_IC2 - IC2 * coef_IC2;
	dIF	 = a12 * IC2 + b3 * C1 + b4 * IM1 + a2 * O + ki_off * D_IF - IF * coef_IF;
	dIM1 =  a4 * IF + b5 * IM2 - IM1 * coef_IM1;
	dIM2 =  a5 * IM1 - IM2 * coef_IM2;
	dOS	 =  ax * O + ki_off * D_OS - OS * coef_OS;

	//Charged Drug Bound States
	dDO		= kon  * O + a13c * DC1 + bx1 * DOS + b22 * DIF  - DO * coef_DO;
	dDC1	= kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO  - DC1 * coef_DC1;
	dDC2	= kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1  - DC2 * coef_DC2;
	dDC3	= kcon * C3 + b11 * DC2 + a33 * DIC3  - DC3 * coef_DC3;
	dDOS	= ax1  * DO -  DOS * coef_DOS;
	dDIC3	= b33  * DC3 + b11 * DIC2 - DIC3 * coef_DIC3;
	dDIC2	= b33  * DC2 + a11 * DIC3 + b12 * DIF - DIC2 * coef_DIC2;
	dDIF	= b33  * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO - DIF * coef_DIF;
	dDIM1	= a44  * DIF + b55 * DIM2 - DIM1 * coef_DIM1;
	dDIM2	= a55  * DIM1 - DIM2 * coef_DIM2;

	//Neutral Drug Bound States
	dD_O	=  k_on * O + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS  - D_O * coef_D_O;
	dD_C1	= kc_on * C1 + a12 * D_C2 + a_33 * D_IF + b13n * D_O  - D_C1 * coef_D_C1;
	dD_C2	= kc_on * C2 + a11 * D_C3 + a_33 * D_IC2 + b12 * D_C1  - D_C2 * coef_D_C2;
	dD_C3	= kc_on * C3 + a_33 * D_IC3 + b11 * D_C2  - D_C3 * coef_D_C3;
	dD_OS	=  ax2  * D_O + ki_on * OS - D_OS * coef_D_OS;
	dD_IC3	=  b_33 * D_C3 + b11 * D_IC2 + ki_on * IC3 - D_IC3 * coef_D_IC3;
	dD_IC2	=  b_33 * D_C2 + a11 * D_IC3 + b12 * D_IF + ki_on * IC2 - D_IC2 * coef_D_IC2;
	dD_IF	=  b_33 * D_C1 + a12 * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF - D_IF * coef_D_IF;
	dD_IM1	=  a_44 * D_IF + b_55 * D_IM2 - D_IM1 * coef_D_IM1;
	dD_IM2	=  a_55 * D_IM1 - D_IM2 * coef_D_IM2;
	I_Na = G_Na* O *(V - E_Na);
}



/*
 * Na current in heart failure cells - Rush Larsen version
 */
void Moreno2011Cell::rush_I_Na_L(const double dt) {
	double G_Na_L = (cellenv->HeartFailure) ? 5.0 * 0.0065 : 0;
	const double tau_h_L = 600;

	const double E_Na = RTF * log(Na_out / Na_in);
	const double amL = 0.32 * (V + 47.13) / (1.0 - exp(-0.1 * (V + 47.13)));
	const double bmL = 0.08 * exp(-V / 11.0);
	const double hLss = 1.0 / (1 + exp((V + 91.0) / 6.1));

	const double tau_m_L = 1. / (amL + bmL);
	const double mLss = amL / (amL + bmL);

	// Rush Larsen formulation:
	mL = mLss - (mLss - mL) * exp(-dt / tau_m_L);
	hL = hLss - (hLss - hL) * exp(-dt / tau_h_L);

	I_Na_L = G_Na_L * mL * mL * mL * hL * (V - E_Na);
}



//Fast sodium current
void Moreno2011Cell::semi_impl_I_Na(double t, const double dt){
	if (t==0)
		WT_SCN5A_Initial_Conditions();
	else {
		// Step from t to t+dt
		if (NaChanSolver == Iterative)
			semi_impl_WT_SCN5A( dt );
		else if (NaChanSolver == LUDecomp)
			semi_impl_WT_SCN5A_LU( dt );

		/* Code used to test whether the above method
		 *  - semi_impl_WT_SCN5A -
		 * is equivalent to the older - two separate - methods:
		 *  - semi_impl_WT_SCN5A_Lidocaine
		 *  - semi_impl_WT_SCN5A_Flecainide
		 * Found to be working 12.23.2014 for both Lidocaine
		 * and Flecainide 1 mM.
		 *
		// Testing: save new version
		gsl_vector* Diag    = gsl_vector_alloc (nNaChStates);
		gsl_vector* DiagInv = gsl_vector_alloc (nNaChStates);
		gsl_matrix* OffDiag = gsl_matrix_alloc (nNaChStates, nNaChStates);
		gsl_vector_memcpy(Diag,    DiagonalAn);
		gsl_vector_memcpy(DiagInv, DiagonalAnInv);
		gsl_matrix_memcpy(OffDiag, OffDiagonalAn);

		// Call old code
		if (cellenv->blocker == Lidocaine)
			semi_impl_WT_SCN5A_Lidocaine(dt);
		else // if cellenv->blocker is Flecainide || None
			semi_impl_WT_SCN5A_Flecainide(dt);

		// Compare new and old
		gsl_vector_sub (Diag, DiagonalAn);
		gsl_vector_sub (DiagInv, DiagonalAnInv);
		gsl_matrix_sub (OffDiag, OffDiagonalAn);
		double diff = gsl_blas_dasum( Diag );
		diff += gsl_blas_dasum( DiagInv );
		for (int i=0; i<nNaChStates; i++) {
		      gsl_vector_view column = gsl_matrix_column (OffDiag, i);
		      diff += gsl_blas_dnrm2 (&column.vector);
		}
		gsl_vector_free( Diag );
		gsl_vector_free( DiagInv );
		gsl_matrix_free( OffDiag );
		if ( fabs(diff) > 1e-12 ){
			std::cout << "Total diff between new and old is " << diff << std::endl;
			throw std::runtime_error("Error - difference!!!");
		}
		*/
	}
}


void Moreno2011Cell::semi_impl_WT_SCN5A(const double dt){
	NaChParams.ComputeVDepTerms( V, cellenv );
	NaChParams.set_NaCh_Matrices( OffDiagonalAn, DiagonalAn );

	const double hh2 = dt/2;

	// DiagonalAnInv = [I - dt/2 Dn ]^-1
	// The following loop should be equivalent to the 30 lines of code setting co_ constants
	for (unsigned int i=0; i<nNaChStates; i++)
		gsl_vector_set( DiagonalAnInv, i, 1.0/( 1.0 - hh2 * gsl_vector_get( DiagonalAn, i )) );

	// Compute Yn+1/2 = Yn + dt/2 An Yn
	//                = Yn + dt/2 (An - Dn) Yn + dt/2 Dn Yn
	// 1 -   Set Yn
	for (unsigned int i=0; i<nNaChStates; i++)
		gsl_vector_set( Yn, i, state[i] );
	// 2 -   Ynph <-- Dn Yn
	gsl_vector_memcpy( Ynph, Yn ); 		// Yn+1/2 <-- Yn
	gsl_vector_mul( Ynph, DiagonalAn ); // Yn+1/2 <-- Dn Yn
	// 3 -  dt/2 An Yn == dt/2 (An-Dn) Yn + dt/2 Dn Yn
	gsl_blas_dgemv( CblasNoTrans, hh2, OffDiagonalAn, Yn, hh2, Ynph );
	// 4 -  add Yn  i.e. Yn+1/2 <-- dt/2 An Yn  +   Yn
	gsl_vector_add( Ynph, Yn );


	// Set Yn+1^0 to be Yn   -- check if setting it to Yn+1/2 reduces the number of iterations required
	gsl_vector_memcpy( Ynp1i, Yn );

	int iter = 0;
	double err_sum = 0;
	do {
		// Set Yn+1 to be Ynph
		gsl_vector_memcpy( Ynp1, Ynph );

		// Yn+1^i+1 = [ I - dt/2 Dn ]^-1 ( Yn+1/2 + dt/2 A~n Yn+1^i )
		gsl_blas_dgemv( CblasNoTrans, hh2, OffDiagonalAn, Ynp1i, 1.0, Ynp1 ); // Yn+1 <-- Yn+1/2 + dt/2 A~n Yn+1^i
		gsl_vector_mul( Ynp1, DiagonalAnInv ); // Yn+1 <-- [I - dt/2 Dn ]^-1 Yn+1

		// find difference of iterates Ynp1i and Ynp1
		gsl_vector_sub( Ynp1i, Ynp1 );
		err_sum = gsl_blas_dasum( Ynp1i ); // double abs value sum

		gsl_vector_memcpy( Ynp1i, Ynp1 );
		iter++;
	} while ( err_sum > 1E-100 && iter < 100 );

	MarkovIterations = iter;
	MarkovIterationError = err_sum;

	// Once the BLAS code is working, use the following loop
	for (unsigned int i=0; i<nNaChStates; i++)
		state[i] = gsl_vector_get( Ynp1, i );

	checkNaChannelStates();

	const double E_Na = RTF*log(Na_out/Na_in);
	I_Na = G_Na*(O)*(V - E_Na);
}



void Moreno2011Cell::semi_impl_WT_SCN5A_LU(const double dt){
	NaChParams.ComputeVDepTerms( V, cellenv );
	NaChParams.set_NaCh_Matrices( OffDiagonalAn, DiagonalAn );

	const double hh2 = dt/2;

	// DiagonalAnInv = [I - dt/2 Dn ]^-1
	// The following loop should be equivalent to the 30 lines of code setting co_ constants
	for (unsigned int i=0; i<nNaChStates; i++)
		gsl_vector_set( DiagonalAnInv, i, 1.0/( 1.0 - hh2 * gsl_vector_get( DiagonalAn, i )) );

	// Compute Yn+1/2 = Yn + dt/2 An Yn
	//                = Yn + dt/2 (An - Dn) Yn + dt/2 Dn Yn
	// 1 -   Set Yn
	for (unsigned int i=0; i<nNaChStates; i++)
		gsl_vector_set( Yn, i, state[i] );
	// 2 -   Ynph <-- Dn Yn
	gsl_vector_memcpy( Ynph, Yn ); 		// Yn+1/2 <-- Yn
	gsl_vector_mul( Ynph, DiagonalAn ); // Yn+1/2 <-- Dn Yn
	// 3 -  dt/2 An Yn == dt/2 (An-Dn) Yn + dt/2 Dn Yn
	gsl_blas_dgemv( CblasNoTrans, hh2, OffDiagonalAn, Yn, hh2, Ynph );
	// 4 -  add Yn  i.e. Yn+1/2 <-- dt/2 An Yn  +   Yn
	gsl_vector_add( Ynph, Yn );
	// -------------------------------------------------------------------------

	// Compute -dt/2 A~n  -- i.e. off-diagonal terms
	gsl_matrix_scale ( OffDiagonalAn, -hh2 );
	// Add in diagonal terms I - dt/2 Dn
	for (unsigned int i=0; i<nNaChStates; i++)
		gsl_matrix_set( OffDiagonalAn, i, i, 1.0 - hh2* gsl_vector_get( DiagonalAn, i) );




	// Constrained version: add constraint: sum of states is 1:
	// This works but does not change behaviour.
	//for (unsigned int i=0; i<nNaChStates; i++)
	//	gsl_matrix_set( OffDiagonalAn,nNaChStates-1, i, 1.0);
	//gsl_vector_set( Ynph, nNaChStates-1, 1.0);




	int s;
	gsl_permutation * p = gsl_permutation_alloc (nNaChStates);
	gsl_linalg_LU_decomp(OffDiagonalAn, p, &s);
	gsl_linalg_LU_solve (OffDiagonalAn, p, Ynph, Ynp1);
	gsl_permutation_free( p );


	// Once the BLAS code is working, use the following loop
	for (unsigned int i=0; i<nNaChStates; i++)
		state[i] = gsl_vector_get( Ynp1, i );

	checkNaChannelStates();

	const double E_Na = RTF*log(Na_out/Na_in);
	I_Na = G_Na*(O)*(V - E_Na);
}

#include <iostream> // for debugging

/*
 * exponential_I_Na1 - first step in 2nd order Rush Larsen
 */
void Moreno2011Cell::exponential_I_Na1(double t, const double dt){
	if (t==0)
		WT_SCN5A_Initial_Conditions();
	else {
		checkNaChannelStates();

		NaChParams.ComputeVDepTerms( V, cellenv );
		NaChParams.set_NaCh_Matrices( OffDiagonalAn, DiagonalAn );

		// 1. Put the state into the Yn vector and add in diagonal terms  Dn to An
		for (unsigned int i=0; i<nNaChStates; i++){
			gsl_vector_set( Yn, i, state[i] );
			gsl_matrix_set( OffDiagonalAn, i, i, gsl_vector_get( DiagonalAn, i) );
		}
		// 2 -  Ynph == An Yn
		gsl_vector_set_zero( Ynp1i );
		gsl_blas_dgemv( CblasNoTrans, 1.0, OffDiagonalAn, Yn, 1.0, Ynp1i );
		// 3. exponential terms
		for (unsigned int i=0; i<nNaChStates; i++){
			const double a = gsl_vector_get( Ynp1i, i);
			const double b = gsl_vector_get( DiagonalAn, i);

			double mmst = 0;
			// The first step is a half-step; that is why we divide by 2:
			if (fabs(b) < tolerance)
				mmst = state[i] + a * dt/2;
			else
				mmst = state[i] + (a/b)*(exp(b*dt/2) - 1.0);
			gsl_vector_set( Ynph, i, mmst );
		}
		// Debugging: gsl_vector_fprintf( stdout, Ynph, "  %f");

		const double E_Na = RTF*log(Na_out/Na_in);
		I_Na = G_Na * O * (V - E_Na);
	}
}

/*
 * exponential_I_Na2 - second step in 2nd order Rush Larsen
 */
void Moreno2011Cell::exponential_I_Na2(double t, const double dt, const double Vhalf){
	if (t==0)
		WT_SCN5A_Initial_Conditions();
	else {
		// 1. The state is in the Ynph vector
		NaChParams.ComputeVDepTerms( Vhalf, cellenv );
		NaChParams.set_NaCh_Matrices( OffDiagonalAn, DiagonalAn );

		// Debugging: gsl_vector_fprintf( stdout, Ynph, "  %f");
		// Debugging: gsl_vector_fprintf( stdout, DiagonalAn, "  %f");
		// Debugging: gsl_matrix_fprintf( stdout, OffDiagonalAn, "  %f");
		// 2 add in diagonal terms Dn to An
		for (unsigned int i=0; i<nNaChStates; i++){
			gsl_matrix_set( OffDiagonalAn, i, i, gsl_vector_get( DiagonalAn, i) );
		}
		// 2 -  Ynp1 == An Ynph
		gsl_vector_set_zero( Ynp1i );
		gsl_blas_dgemv( CblasNoTrans, 1.0, OffDiagonalAn, Ynph, 1.0, Ynp1i );
		// 3. exponential terms
		for (unsigned int i=0; i<nNaChStates; i++){
			const double a = gsl_vector_get( Ynp1i, i);
			const double b = gsl_vector_get( DiagonalAn, i);

			double mmst = 0;
			if (fabs(b) < tolerance)
				mmst = state[i] + a * dt;
			else
				mmst = state[i] + (a/b)*(exp(b*dt) - 1.0);
			state[i] = mmst; // state[i] + (a/b)*(exp(b*dt) - 1.0);
		}
		checkNaChannelStates();
		const double E_Na = RTF*log(Na_out/Na_in);
		I_Na = G_Na * O * (Vhalf - E_Na);
	}
}

//L type calcium channel Ca contribution
void Moreno2011Cell::rush_I_Ca_L(const double dt){
	const double G_Ca_L = (cellenv->HeartFailure) ? 2.5 * 3.980E-5 : 3.980E-5;

	const double d_inf = 1. / (1. + exp((-8 - V) / 7.5));
	const double a_d = 1.4 / (1. + exp((-35 - V) / 13)) + 0.25;
	const double b_d = 1.4 / (1. + exp((V + 5) / 5));
	const double c_d = 1. / (1. + exp((50 - V) / 20));
	const double tau_d = a_d * b_d + c_d;

	// Rush Larsen:
	d = d_inf - (d_inf - d) * exp(-dt / tau_d);

	const double f_inf = 1. / (1. + exp((V + 20) / 7));
	const double a_f = 1102.5 * exp(-(V + 27) * (V + 27) / 225);
	const double b_f = 200. / (1 + exp((13 - V) / 10.));
	const double c_f = (180. / (1 + exp((V + 30) / 10))) + 20;
	const double tau_f = a_f + b_f + c_f;

	// Rush Larsen:
	f = f_inf - (f_inf - f) * exp(-dt / tau_f);

	const double f2_inf = 0.67 / (1. + exp((V + 35) / 7)) + 0.33;
	const double a_f2 = 600 * exp(-(V + 25) * (V + 25) / 170);
	const double b_f2 = 31 / (1. + exp((25 - V) / 10));
	const double c_f2 = 16 / (1. + exp((V + 30) / 10));
	const double tau_f2 = a_f2 + b_f2 + c_f2;

	// Rush Larsen:
	f2 = f2_inf - (f2_inf - f2) * exp(-dt / tau_f2);

	const double f_Ca_inf = 0.6 / (1 + (Ca_ss / 0.05) * (Ca_ss / 0.05)) + 0.4;
	const double tau_f_Ca = 80. / (1 + (Ca_ss / 0.05) * (Ca_ss / 0.05)) + 2.;

	// Rush Larsen:
	f_Ca = f_Ca_inf - (f_Ca_inf - f_Ca) * exp(-dt / tau_f_Ca);

	const double Vm15RT2F = z_Ca * (V - 15)/RTF;
	// Need to apply L'Hopital's rule if V is close to 15
	I_Ca_L = G_Ca_L * F * d * f * f2 * f_Ca * 2 * Vm15RT2F
			* (0.25 * exp( Vm15RT2F ) * Ca_ss - Ca_out)
			/ (exp( Vm15RT2F ) - 1.);
}

//Rapid rectifier
void Moreno2011Cell::rush_I_Kr(const double dt){
	const double G_Kr = 0.153;		// nS/pF
	const double E_Kr = RTF * log(K_out / K_in);

	const double a_xr1 = 450. / (1. + exp((-45. - V) / 10.));
	const double b_xr1 = 6. / (1. + exp((V - (-30.)) / 11.5));
	const double xr1_ss = 1. / (1. + exp((-26. - V) / 7.));
	const double tau_xr1 = a_xr1 * b_xr1;
	// Rush Larsen
	xr1 = xr1_ss - (xr1_ss - xr1) * exp(-dt / tau_xr1);

	const double a_xr2 = 3. / (1. + exp((-60. - V) / 20.));
	const double b_xr2 = 1.12 / (1. + exp((V - 60.) / 20.));
	const double xr2_ss = 1. / (1. + exp((V - (-88.)) / 24.));
	const double tau_xr2 = a_xr2 * b_xr2;
	// Rush Larsen
	xr2 = xr2_ss - (xr2_ss - xr2) * exp(-dt / tau_xr2);

	I_Kr = G_Kr * sqrt(K_out / 5.4) * xr1 * xr2 * (V - E_Kr);
}

//Slow rectifier
void Moreno2011Cell::rush_I_Ks(const double dt){
	double G_Ks = 0.392;
	switch (Cell_type) {
	case 1:	G_Ks = 0.392;	break;	// Cell type == 1 for endo
	case 2:	G_Ks = 0.098;	break;	// Cell type == 2 for M
	case 3:	G_Ks = 0.392;	break;	// Cell type == 3 for epi
	default: ;
	}

	const double PR_NaK = 0.03;
	const double E_Ks = RTF * log((K_out + PR_NaK * Na_out) / (K_in + PR_NaK * Na_in));

	const double xs_ss = 1. / (1. + exp((-5. - V) / 14.));
	const double ax_s = (1400. / (sqrt(1. + exp((5. - V) / 6))));
	const double bx_s = (1. / (1. + exp((V - 35.) / 15.)));
	const double tau_xs = (ax_s * bx_s) + 80;

	// Rush Larsen
	xs = xs_ss - (xs_ss - xs) * exp(-dt / tau_xs);
	I_Ks = G_Ks * xs * xs * (V - E_Ks);
}

//Transient outward current
void Moreno2011Cell::rush_I_to(const double dt){
	const double endoGtomax  = 0.073;
	const double McellGtomax = 0.294;
	const double epiGtomax   = 0.294;
	double G_to;
	// default values for s_inf and tau_s:
	double s_inf = 1. / (1. + exp((V + 20) / 5.));
	double tau_s = 85. * exp(-(V + 45.) * (V + 45.) / 320.)
			+ 5. / (1. + exp((V - 20.) / 5.)) + 3.;

	switch (Cell_type) {
	case 1:				// Cell type = (1) for endocardial
		G_to = (cellenv->HeartFailure) ? 0.64 * endoGtomax : endoGtomax;		// nS/pF
		s_inf = 1. / (1. + exp((V + 28) / 5.));
		tau_s = 1000. * exp(-(V + 67) * (V + 67) / 1000.) + 8.;
		break;
	case 2:				// Cell type = (2) for M cells
		G_to = (cellenv->HeartFailure) ? 0.64 * McellGtomax : McellGtomax;	// nS/pF
		// default values for s_inf and tau_s
		break;
	case 3:				// Cell type = (3) for epicardial
	default:
		G_to = (cellenv->HeartFailure) ? 0.64 * epiGtomax : epiGtomax;		// nS/pF
	}

	const double E_to = ((R * T) / (z_K * F)) * log(K_out / K_in);

	const double r_inf = 1. / (1. + exp((20 - V) / 6.));
	const double tau_r = 9.5 * exp(-(V + 40.) * (V + 40.) / 1800.) + 0.8;

	// Rush Larsen
	r = r_inf - (r_inf - r) * exp(-dt / tau_r);
	s = s_inf - (s_inf - s) * exp(-dt / tau_s);

	I_to = G_to * r * s * (V - E_to);
}

//Release from SR
void Moreno2011Cell::rush_I_rel(const double dt){
	double k_Ca_sr, k1_rel, k2_rel, dR_bar;
	const double G_rel = 0.102;	// Max. rate constant of Ca release from JSR due to overload (mM/ms)
	const double k1_rel_ = 0.15;// R to O and RI to I I_rel transition rate (mM^2/ms)
	const double k2_rel_ = 0.045;// O to I and R to RI I_rel transition rate (mM^2/ms)
	const double k3_rel = 0.060;// O to R and I to RI I_rel transition rate (ms^-1)
	const double k4_rel = 0.005;// I to O and RI to I I_rel transition rate (ms^-1)
	const double EC_sr = 1.5;		// Ca_sr half-saturation constant of k_Ca_sr
	const double max_sr = 2.5;			// Max value of k_Ca_sr (dimensionless)
	const double min_sr = 1;			// Min value of k_Ca_sr (dimensionless)

	k_Ca_sr = max_sr
			- ((max_sr - min_sr) / (1 + (EC_sr / Ca_sr) * (EC_sr / Ca_sr)));

	k1_rel = k1_rel_ / k_Ca_sr;
	k2_rel = k2_rel_ * k_Ca_sr;

	dR_bar = (-k2_rel * Ca_ss * R_bar + k4_rel * (1 - R_bar)) * dt;
	R_bar = R_bar + dR_bar;

	double OO = (k1_rel * Ca_ss * Ca_ss * R_bar) / (k3_rel + k1_rel * Ca_ss * Ca_ss);

	I_rel = G_rel * OO * (Ca_sr - Ca_ss);
}

//Dynamic [Na] in Myoplasm
void Moreno2011Cell::rush_Na_in(const double dt){
	Calculate_Na_in();
	Na_in = Na_in + dt * dNa_in;
}

//Dynamic [K] in Myoplasm
void Moreno2011Cell::rush_K_in(const double dt){
	Calculate_K_in ();
	K_in = K_in + dt * dK_in;
}

//Dynamic [Ca] in Myoplasm
// Old version of Calculate_Ca_in used in Moreno et al. 2011
//  New Ca concentration in cytoplasm
void Moreno2011Cell::rush_Ca_in(const double dt){
	const double buffer = 0.2;	//Total cytoplasmic buffer concentration (mM)
	const double K_buffer = 0.001;// Ca_in half-saturation constant for cytoplasmic buffer (mM)

	const double Ca_in_buffer = (Ca_in * buffer) / (Ca_in + K_buffer);
	dCa_in = (-(I_Ca_b + I_p_Ca - 2 * I_Na_Ca) * CAP / (V_cyto * z_Ca * F))
					+ (V_sr / V_cyto) * (I_leak - I_up) + I_tr;

	const double b_Ca = buffer - Ca_in_buffer - dCa_in * dt - Ca_in + K_buffer;
	const double c_Ca = K_buffer * (Ca_in_buffer + dCa_in * dt + Ca_in);
	Ca_in = (sqrt(b_Ca * b_Ca + 4 * c_Ca) - b_Ca) / 2;
}

//Dynamic [Ca] in JSR
//  Old version of Calculate_Ca_in used in Moreno et al. 2011
//  New Ca concentration in SR
void Moreno2011Cell::rush_Ca_sr(const double dt){
	const double buffer_sr = 10;			//Total SR buffer concentration (mM)
	const double K_buffer_sr = 0.3;	// Ca_sr half-saturation constant for SR buffer (mM)

	const double Ca_sr_buffer = (Ca_sr * buffer_sr) / (Ca_sr + K_buffer_sr);
	dCa_sr = I_up - I_rel - I_leak;

	const double b_jsr = buffer_sr - Ca_sr_buffer - dCa_sr * dt - Ca_sr + K_buffer_sr;
	const double c_jsr = K_buffer_sr * (Ca_sr_buffer + dCa_sr * dt + Ca_sr);
	Ca_sr = (sqrt(b_jsr * b_jsr + 4 * c_jsr) - b_jsr) / 2;
}

//Dynamic [Ca] in NSR
void Moreno2011Cell::rush_Ca_ss(const double dt){
	double b_Ca_ss, c_Ca_ss, dCa_ss;
	const double buffer_ss = 0.4;			//Total SS buffer concentration (mM)
	const double K_buffer_ss = 0.00025;	// Ca_ss half-saturation constant for SR buffer (mM)

	const double Ca_ss_buffer = (Ca_ss * buffer_ss) / (Ca_ss + K_buffer_ss);
	dCa_ss = dt
			* ((-I_Ca_L * CAP) / (V_ss * z_Ca * F) + (V_sr / V_ss) * I_rel
					- (V_cyto / V_ss) * I_tr);

	b_Ca_ss = buffer_ss - Ca_ss_buffer - dCa_ss - Ca_ss + K_buffer_ss;
	c_Ca_ss = K_buffer_ss * (Ca_ss_buffer + dCa_ss + Ca_ss);

	Ca_ss = (sqrt(b_Ca_ss * b_Ca_ss + 4 * c_Ca_ss) - b_Ca_ss) / 2;
}

void Moreno2011Cell::NaChannelParameters::set_NaCh_Matrices( gsl_matrix* M, gsl_vector* D){
	// The following 30 consts are the diagonal coefficients of the A matrix
	const double coef_IC3	= (a11 + a3 + ki_on);
	const double coef_IC2	= (b11 + a3 + a12 + ki_on);
	const double coef_IF	= (b12 + b2 + a3 + a4 + ki_on);
	const double coef_IM1	= (b4 + a5);
	const double coef_IM2	= b5;
	const double coef_C3	= (b3 + a11  + kcon + kc_on);
	const double coef_C2	= (b11 + b3 + a12  + kcon + kc_on);
	const double coef_O		= (b13 + a2  + kon + k_on + ax);
	const double coef_OS	= (bx + ki_on);
	const double coef_C1	= (b12 + b3 + a13  + kcon + kc_on);


	const double coef_DIC3	= (a11 + a33);
	const double coef_DIC2	= (a33 + b11 + a12);
	const double coef_DIF	= (a33 + b12 + a44 + b22);
	const double coef_DIM1	= ( b44 + a55 );
	const double coef_DIM2	= b55 ;
	const double coef_DC3	= (kcoff+ b33 + a11 );
	const double coef_DC2	= (kcoff + b11 + b33 + a12 );
	const double coef_DO	= (koff + b13c + a22 + ax1 );
	const double coef_DOS	= bx1;
	const double coef_DC1	= (kcoff + b12 + b33 + a13c );

	const double coef_D_IC3	= (a_33 + a11 + ki_off);
	const double coef_D_IC2	= (a_33 + b11 + a12 + ki_off);
	const double coef_D_IF	= (a_33 + a_44 + b_22 + b12 + ki_off);
	const double coef_D_IM1	= (b_44 + a_55);
	const double coef_D_IM2	= b_55;
	const double coef_D_C3	= (kc_off + b_33 + a11 );
	const double coef_D_C2	= (kc_off + b11 + b_33 + a12 );
	const double coef_D_O	= (k_off + b13n + a_22 + ax2 );
	const double coef_D_OS	= (bx2 + ki_off);
	const double coef_D_C1	= (kc_off + b12 + b_33 + a13n );

	// Set up diagonal coefficients
	gsl_vector_set_all( D, 0.0 );
	gsl_vector_set( D,  i3, -coef_IC3 );
	gsl_vector_set( D,  i2, -coef_IC2 );
	gsl_vector_set( D,  iF, -coef_IF  );
	gsl_vector_set( D,  im, -coef_IM1 );
	gsl_vector_set( D,  iM, -coef_IM2 );
	gsl_vector_set( D,  c3, -coef_C3  );
	gsl_vector_set( D,  c2, -coef_C2  );
	gsl_vector_set( D,  op, -coef_O	  );
	gsl_vector_set( D,  oS, -coef_OS  );
	gsl_vector_set( D,  c1, -coef_C1  );

	gsl_vector_set( D, di3, -coef_DIC3 );
	gsl_vector_set( D, di2, -coef_DIC2 );
	gsl_vector_set( D, diF, -coef_DIF  );
	gsl_vector_set( D, dim, -coef_DIM1 );
	gsl_vector_set( D, diM, -coef_DIM2 );
	gsl_vector_set( D, dc3, -coef_DC3  );
	gsl_vector_set( D, dc2, -coef_DC2  );
	gsl_vector_set( D, dop, -coef_DO   );
	gsl_vector_set( D, doS, -coef_DOS  );
	gsl_vector_set( D, dc1, -coef_DC1  );

	gsl_vector_set( D, Di3, -coef_D_IC3 );
	gsl_vector_set( D, Di2, -coef_D_IC2 );
	gsl_vector_set( D, DiF, -coef_D_IF  );
	gsl_vector_set( D, Dim, -coef_D_IM1 );
	gsl_vector_set( D, DiM, -coef_D_IM2 );
	gsl_vector_set( D, Dc3, -coef_D_C3	);
	gsl_vector_set( D, Dc2, -coef_D_C2	);
	gsl_vector_set( D, Dop, -coef_D_O	);
	gsl_vector_set( D, DoS, -coef_D_OS	);
	gsl_vector_set( D, Dc1, -coef_D_C1	);


	gsl_matrix_set_all( M, 0 );
	// Drug-free states
	gsl_matrix_set( M,  i3,   c3, b3     ); // Row i3
	gsl_matrix_set( M,  i3,   i2, b11    );
	gsl_matrix_set( M,  i3,  Di3, ki_off );
	gsl_matrix_set( M,  i2,   i3, a11    ); // Row i2
	gsl_matrix_set( M,  i2,   c2, b3     );
	gsl_matrix_set( M,  i2,   iF, b12    );
	gsl_matrix_set( M,  i2,  Di2, ki_off );
	gsl_matrix_set( M,  iF,   i2, a12    ); // Row iF
	gsl_matrix_set( M,  iF,   c1, b3     );
	gsl_matrix_set( M,  iF,   im, b4     );
	gsl_matrix_set( M,  iF,   op, a2	 );
	gsl_matrix_set( M,  iF,  DiF, ki_off );
	gsl_matrix_set( M,  im,   iF, a4	 ); // Row im
	gsl_matrix_set( M,  im,   iM, b5	 );
	gsl_matrix_set( M,  iM,   im, a5	 ); // Row iM
	gsl_matrix_set( M,  c3,   i3, a3	 ); // Row c3
	gsl_matrix_set( M,  c3,   c2, b11	 );
	gsl_matrix_set( M,  c3,  dc3, kcoff  );
	gsl_matrix_set( M,  c3,  Dc3, kc_off );
	gsl_matrix_set( M,  c2,   c3, a11	 ); // Row c2
	gsl_matrix_set( M,  c2,   i2, a3	 );
	gsl_matrix_set( M,  c2,   c1, b12	 );
	gsl_matrix_set( M,  c2,  dc2, kcoff  );
	gsl_matrix_set( M,  c2,  Dc2, kc_off );
	gsl_matrix_set( M,  op,   c1, a13	 ); // Row op
	gsl_matrix_set( M,  op,   iF, b2	 );
	gsl_matrix_set( M,  op,  dop, koff	 );
	gsl_matrix_set( M,  op,  Dop, k_off	 );
	gsl_matrix_set( M,  op,   oS, bx	 );
	gsl_matrix_set( M,  oS,   op, ax	 ); // Row oS
	gsl_matrix_set( M,  oS,  DoS, ki_off );
	gsl_matrix_set( M,  c1,   c2, a12	 ); // Row c1
	gsl_matrix_set( M,  c1,   iF, a3	 );
	gsl_matrix_set( M,  c1,   op, b13	 );
	gsl_matrix_set( M,  c1,  dc1, kcoff  );
	gsl_matrix_set( M,  c1,  Dc1, kc_off );

	// Charged-drug states
	gsl_matrix_set( M,  di3,  dc3, b33 ); // Row di3
	gsl_matrix_set( M,  di3,  di2, b11 );
	gsl_matrix_set( M,  di2,  di3, a11 ); // Row di2
	gsl_matrix_set( M,  di2,  dc2, b33 );
	gsl_matrix_set( M,  di2,  diF, b12 );
	gsl_matrix_set( M,  diF,  di2, a12 ); // Row diF
	gsl_matrix_set( M,  diF,  dc1, b33 );
	gsl_matrix_set( M,  diF,  dim, b44 );
	gsl_matrix_set( M,  diF,  dop, a22 );
	gsl_matrix_set( M,  dim,  diF, a44 ); // Row dim
	gsl_matrix_set( M,  dim,  diM, b55 );
	gsl_matrix_set( M,  diM,  dim, a55 ); // Row diM
	gsl_matrix_set( M,  dc3,  di3, a33 ); // Row dc3
	gsl_matrix_set( M,  dc3,  dc2, b11 );
	gsl_matrix_set( M,  dc3,   c3, kcon);
	gsl_matrix_set( M,  dc2,  dc3, a11 ); // Row dc2
	gsl_matrix_set( M,  dc2,  di2, a33 );
	gsl_matrix_set( M,  dc2,  dc1, b12 );
	gsl_matrix_set( M,  dc2,   c2, kcon);
	gsl_matrix_set( M,  dop,  dc1, a13c); // Row dop
	gsl_matrix_set( M,  dop,  diF, b22 );
	gsl_matrix_set( M,  dop,   op, kon );
	gsl_matrix_set( M,  dop,  doS, bx1 );
	gsl_matrix_set( M,  doS,  dop, ax1 ); // Row doS
	gsl_matrix_set( M,  dc1,  dc2, a12 ); // Row dc1
	gsl_matrix_set( M,  dc1,  diF, a33 );
	gsl_matrix_set( M,  dc1,  dop, b13c);
	gsl_matrix_set( M,  dc1,   c1, kcon);

	// Neutral-drug states
	gsl_matrix_set( M,  Di3,  Dc3, b_33 ); // Row Di3
	gsl_matrix_set( M,  Di3,  Di2, b11  );
	gsl_matrix_set( M,  Di3,   i3, ki_on);
	gsl_matrix_set( M,  Di2,  Di3, a11  ); // Row Di2
	gsl_matrix_set( M,  Di2,  Dc2, b_33 );
	gsl_matrix_set( M,  Di2,  DiF, b12  );
	gsl_matrix_set( M,  Di2,   i2, ki_on);
	gsl_matrix_set( M,  DiF,  Di2, a12  ); // Row DiF
	gsl_matrix_set( M,  DiF,  Dc1, b_33 );
	gsl_matrix_set( M,  DiF,  Dim, b_44 );
	gsl_matrix_set( M,  DiF,  Dop, a_22 );
	gsl_matrix_set( M,  DiF,   iF, ki_on);
	gsl_matrix_set( M,  Dim,  DiF, a_44 ); // Row Dim
	gsl_matrix_set( M,  Dim,  DiM, b_55 );
	gsl_matrix_set( M,  DiM,  Dim, a_55 ); // Row DiM
	gsl_matrix_set( M,  Dc3,  Di3, a_33 ); // Row Dc3
	gsl_matrix_set( M,  Dc3,  Dc2, b11  );
	gsl_matrix_set( M,  Dc3,   c3, kc_on);
	gsl_matrix_set( M,  Dc2,  Dc3, a11  ); // Row Dc2
	gsl_matrix_set( M,  Dc2,  Di2, a_33 );
	gsl_matrix_set( M,  Dc2,  Dc1, b12  );
	gsl_matrix_set( M,  Dc2,   c2, kc_on);
	gsl_matrix_set( M,  Dop,  Dc1, a13n ); // Row Dop
	gsl_matrix_set( M,  Dop,  DiF, b_22 );
	gsl_matrix_set( M,  Dop,   op, k_on );
	gsl_matrix_set( M,  Dop,  DoS, bx2  );
	gsl_matrix_set( M,  DoS,  Dop, ax2  ); // Row DoS
	gsl_matrix_set( M,  DoS,   oS, ki_on);
	gsl_matrix_set( M,  Dc1,  Dc2, a12  ); // Row Dc1
	gsl_matrix_set( M,  Dc1,  DiF, a_33 );
	gsl_matrix_set( M,  Dc1,  Dop, b13n );
	gsl_matrix_set( M,  Dc1,   c1, kc_on);
}


#include <stdexcept>

// The following code is a workaround for MinGW-32: it does not have std::to_string
// ---- start patch
#include <string>
#include <sstream>

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}
// On Win32 systems use patch::to_string instead of std::to_string
// ---- end patch

/*
 * Check sum
 */
void Moreno2011Cell::checkNaChannelStates(){
	const double tolerance = 1e-4;
	double sum = 0;
	for (unsigned int i=nNaChFirstState; i<=nNaChLastState; i++)
		if (state[i] < -tolerance){
			throw std::runtime_error(
					"Negative state value in Na Channel state " +
					std::to_string(i) );
		} else
			sum = sum + state[i];
	if ( fabs (sum - 1.0) > tolerance ) {
		throw std::runtime_error(
				"Error in Na Channel - sum of states,\t sum = " +
				patch::to_string(sum) );
	}
}



#include <iostream>

void Moreno2011Cell::semi_impl_WT_SCN5A_Lidocaine(const double dt){
	const double pKa	= 7.6;
	const double diffusion	= 500;
	const double portion	= 1/(1+ pow(10, (pH-pKa)) );
	const double drug_charged  = cellenv->drug*portion;
	const double drug_neutral  = cellenv->drug*(1-portion);
	const double drug_distance = -0.7;

	const double E_Na = RTF*log(Na_out/Na_in);

	//Rate Constants **********************************************************
	//WT Fits Reduced Model (no IM1, IM2)
	const double Tfactor = 1.0/(pow(Q10, (37.0-(T-273))/10.0));

	const double a11 = Tfactor*8.5539/(7.4392e-2*exp(-V/17.0)+ 2.0373e-1*exp(-V/150));
	const double a12 = Tfactor*8.5539/(7.4392e-2*exp(-V/15.0)+ 2.0373e-1*exp(-V/150));
	const double a13 = Tfactor*8.5539/(7.4392e-2*exp(-V/12.0)+ 2.0373e-1*exp(-V/150));
	const double b11 = Tfactor*7.5215e-2*exp(-V/20.3);
	const double b12 = Tfactor*2.7574*exp(-(V-5)/20.3);
	const double b13 = Tfactor*4.7755e-1*exp(-(V-10)/20.3);

	const double a3	= Tfactor*5.1458e-6*exp(-V/8.2471);
	const double b3	= Tfactor*6.1205*exp((V)/13.542);
	const double a2	= Tfactor*(13.370*exp(V/43.749));
	const double b2	= a13*a2*a3/(b13*b3);

	const double a4 = 0*a2;
	const double b4 = 0*a3;
	const double a5	= 0*a2;
	const double b5 = 0*a3;

	const double ax = 3.4229e-2*a2;
	const double bx = 1.7898e-2*a3;


	const double ax1 = 6.3992e-07  *ax;
	const double bx1 = 1.3511e+00  *bx;
	const double a13c= 5.6974e-03  *a13;
	const double a22 = 6.7067e-06  *a2;
	const double b33 =  1.9698e-05 *b3;
	const double a33 = 3.2976e+00  *a3;

	const double a44 = 0  *a2;
	const double b44 = 0  *a3;
	const double a55 = 0;
	const double b55 = 0;

	const double ax2  = 1.3110e-01  *ax;
	const double a13n = 8.4559e+01  *a13;
	const double a_22 = 1.7084e-05  *a2;
	const double b_33 = 4.8477e+00  *b3;

	const double a_44 = 0  *a2;
	const double b_44 = 0  *a3;
	const double a_55 = 0;
	const double b_55 = 0;


	const double kd0	= 318e-6;
	const double kd_open= kd0*exp( (drug_distance*V*F) /(R*T));

	// charged drug
	const double kon	= drug_charged*diffusion;
	const double koff	= kd_open*diffusion;
	const double kcoff	= koff;
	const double kcon	= kon;

	const double b13c	= (cellenv->drug ==0 || drug_charged ==0 ) ? 0 : b13*kcon*koff*a13c/(kon*kcoff*a13);
	const double b22	= (b13c ==0) ? 0 : (a13c*a22*a33)/(b13c*b33);

	// neutral drug
	const double k_on	= drug_neutral*diffusion;
	const double k_off	= 400e-6*diffusion;
	const double ki_on	= k_on/2;
	const double ki_off	= 3.4e-6*diffusion;
	const double kc_on	= k_on/2;
	const double kc_off	= 900e-6*diffusion;

	const double a_33 = (cellenv->drug ==0 || drug_neutral ==0 ) ? 0 : (ki_off*a3*kc_on*b_33)/(ki_on*kc_off*b3);
	const double b13n = (cellenv->drug ==0 || drug_neutral ==0 ) ? 0 : (b13*kc_on*a13n*k_off)/(kc_off*a13*k_on);
	const double b_22 = (b13n == 0) ? 0 : (a_33*a13n*a_22)/(b_33*b13n);
	const double bx2  = (cellenv->drug ==0 || drug_neutral ==0) ? 0 : (bx*k_on*ax2*ki_off)/(ax*ki_on*k_off);

	set_NaCh_Matrices( OffDiagonalAn, DiagonalAn,
			a11,  a12,	a13,
			b11,  b12,	b13,
			a3,	  b3,	a2,    b2,
			a4,	  b4,	a5,	   b5,
			ax,	  bx,	ax1,   bx1,
			a13c, a22,	b33,   a33,
			a44,  b44,	a55,   b55,
			ax2,  a13n, a_22,  b_33,
			a_44, b_44, a_55,  b_55,
			kon,  koff, kcoff, kcon,
			b13c, b22,  k_on,  k_off,
			ki_on,ki_off, kc_on, kc_off,
			a_33, b13n, b_22,  bx2 );

	const double hh	 = dt;
	const double hh2 = hh/2;

	// DiagonalAnInv = [I - dt/2 Dn ]^-1
	// The following loop should be equivalent to the 30 lines of code setting co_ constants
	for (unsigned int i=0; i<nNaChStates; i++)
		gsl_vector_set( DiagonalAnInv, i, 1.0/( 1.0 - hh2 * gsl_vector_get( DiagonalAn, i )) );

	// Compute Yn+1/2 = Yn + dt/2 An Yn
	//                = Yn + dt/2 (An - Dn) Yn + dt/2 Dn Yn
	// 1 -   Set Yn
	for (unsigned int i=0; i<nNaChStates; i++)
		gsl_vector_set( Yn, i, state[i] );
	// 2 -   Ynph <-- Dn Yn
	gsl_vector_memcpy( Ynph, Yn ); 		// Yn+1/2 <-- Yn
	gsl_vector_mul( Ynph, DiagonalAn ); // Yn+1/2 <-- Dn Yn
	// 3 -  dt/2 An Yn == dt/2 (An-Dn) Yn + dt/2 Dn Yn
	gsl_blas_dgemv( CblasNoTrans, hh2, OffDiagonalAn, Yn, hh2, Ynph );
	// 4 -  add Yn  i.e. Yn+1/2 <-- dt/2 An Yn  +   Yn
	gsl_vector_add( Ynph, Yn );


	// Set Yn+1^0 to be Yn   -- check if setting it to Yn+1/2 reduces the number of iterations required
	gsl_vector_memcpy( Ynp1i, Yn );

	int iter = 0;
	double err_sum = 0;
	do {
		// Set Yn+1 to be Ynph
		gsl_vector_memcpy( Ynp1, Ynph );

		// Yn+1^i+1 = [ I - dt/2 Dn ]^-1 ( Yn+1/2 + dt/2 A~n Yn+1^i )
		gsl_blas_dgemv( CblasNoTrans, hh2, OffDiagonalAn, Ynp1i, 1.0, Ynp1 ); // Yn+1 <-- Yn+1/2 + dt/2 A~n Yn+1^i
		gsl_vector_mul( Ynp1, DiagonalAnInv ); // Yn+1 <-- [I - dt/2 Dn ]^-1 Yn+1

		// find difference of iterates Ynp1i and Ynp1
		gsl_vector_sub( Ynp1i, Ynp1 );
		err_sum = gsl_blas_dasum( Ynp1i ); // double abs value sum

		gsl_vector_memcpy( Ynp1i, Ynp1 );
		iter++;
	} while ( err_sum > 1E-100 && iter < 100 );

	MarkovIterations = iter;
	MarkovIterationError = err_sum;

	// Once the BLAS code is working, use the following loop
	for (unsigned int i=0; i<nNaChStates; i++)
		state[i] = gsl_vector_get( Ynp1, i );

	checkNaChannelStates();

	I_Na = G_Na*(O)*(V - E_Na);
}


void Moreno2011Cell::semi_impl_WT_SCN5A_Flecainide(const double dt){
	const double pKa	= 9.3;
	const double diffusion	= 5500;
	const double portion	= 1/(1+ pow(10, (pH-pKa)) );
	const double drug_charged	= cellenv->drug*portion;
	const double drug_neutral	= cellenv->drug*(1-portion);
	const double drug_distance	= -0.7;

	const double E_Na = RTF*log(Na_out/Na_in);

	//Rate Constants **********************************************************
	//WT Fits Reduced Model (no IM1, IM2)
	const double Tfactor = 1.0/(pow(Q10, (37.0-(T-273))/10.0));

	const double a11= Tfactor*8.5539/(7.4392e-2*exp(-V/17.0)+ 2.0373e-1*exp(-V/150));
	const double a12= Tfactor*8.5539/(7.4392e-2*exp(-V/15.0)+ 2.0373e-1*exp(-V/150));
	const double a13= Tfactor*8.5539/(7.4392e-2*exp(-V/12.0)+ 2.0373e-1*exp(-V/150));
	const double b11= Tfactor*7.5215e-2*exp(-V/20.3);
	const double b12= Tfactor*2.7574*exp(-(V-5)/20.3);
	const double b13= Tfactor*4.7755e-1*exp(-(V-10)/20.3);

	const double a3 = Tfactor*5.1458e-6*exp(-V/8.2471);
	const double b3=Tfactor*6.1205*exp((V)/13.542);
	const double a2= Tfactor*(13.370*exp(V/43.749));
	const double b2= a13*a2*a3/(b13*b3);

	const double a4 = 0*a2;
	const double b4 = 0*a3;
	const double a5 = 0*a2;
	const double b5 = 0*a3;

	const double ax = 3.4229e-2*a2;
	const double bx = 1.7898e-2*a3;

	const double ax1 =5.7839e-05 * ax;
	const double bx1 =  1.6689e-08* bx;
	const double a13c = 3.6324e-03 *a13;
	const double a22 = 1.4847e+03 *a2;
	const double b33 =  1.7352e-06* b3;
	const double a33 = 6.7505e-05 * a3;

	const double a44 =  2.4135e+00* a2;
	const double b44 =  4.9001e-02* a3;
	const double a55 =  0;
	const double b55 = 0;

	const double ax2 = 2.6126e-01 * ax;
	const double a13n = 2.6452e+00 * a13;
	const double a_22 =  4.2385e+01 * a2;
	const double b_33 = 2.1181e+00 * b3;
	const double a_44 =  1.0326e-03 * a2;
	const double b_44 = 2.1378e-02 * a3;

	const double a_55 = 0;
	const double b_55 = 0;

	const double kd0=11.2*(1e-6);
	const double kd_open=kd0*exp( (drug_distance*V*F) /(R*T));

	// charged drug
	const double kon	= drug_charged*diffusion;
	const double koff	= kd_open*diffusion;
	const double kcoff	= koff;
	const double kcon	= kon;

	const double b13c	= (cellenv->drug ==0 || drug_charged ==0 ) ?  0 : b13*kcon*koff*a13c/(kon*kcoff*a13);
	const double b22	= (b13c == 0) ? 0 : a13c*a22*a33/(b13c*b33);

	// neutral drug
	const double k_on = drug_neutral*diffusion;
	const double k_off=400*(1e-6)*diffusion;
	const double ki_on=k_on/2;
	const double ki_off=5.4*(1e-6)*diffusion;
	const double kc_on=k_on/2;
	const double kc_off=800*(1e-6)*diffusion;

	const double a_33	= (cellenv->drug ==0 || drug_neutral ==0 ) ? 0 : ki_off*a3*kc_on*b_33/(ki_on*kc_off*b3);
	const double b13n	= (cellenv->drug ==0 || drug_neutral ==0)  ? 0 : b13*kc_on*a13n*k_off/(kc_off*a13*k_on);
	const double b_22	= (b13n == 0) ? 0 : a_33*a13n*a_22/(b_33*b13n);
	const double bx2	= (cellenv->drug ==0 || drug_neutral ==0) ? 0 : bx*k_on*ax2*ki_off/(ax*ki_on*k_off);

	const double hh	 = dt;
	const double hh2 = hh/2;

	set_NaCh_Matrices( OffDiagonalAn, DiagonalAn,
			a11,  a12,	a13,
			b11,  b12,	b13,
			a3,	  b3,	a2,    b2,
			a4,	  b4,	a5,	   b5,
			ax,	  bx,	ax1,   bx1,
			a13c, a22,	b33,   a33,
			a44,  b44,	a55,   b55,
			ax2,  a13n, a_22,  b_33,
			a_44, b_44, a_55,  b_55,
			kon,  koff, kcoff, kcon,
			b13c, b22,  k_on,  k_off,
			ki_on,ki_off, kc_on, kc_off,
			a_33, b13n, b_22,  bx2 );

	// DiagonalAnInv = [I - dt/2 Dn ]^-1
	// The following loop should be equivalent to the 30 lines of code setting co_ constants
	for (unsigned int i=0; i<nNaChStates; i++)
		gsl_vector_set( DiagonalAnInv, i, 1.0/( 1.0 - hh2 * gsl_vector_get( DiagonalAn, i )) );

	// Compute Yn+1/2 = Yn + dt/2 An Yn
	//                = Yn + dt/2 (An - Dn) Yn + dt/2 Dn Yn
	// 1 -   Set Yn
	for (unsigned int i=0; i<nNaChStates; i++)
		gsl_vector_set( Yn, i, state[i] );
	// 2 -   Ynph <-- Dn Yn
	gsl_vector_memcpy( Ynph, Yn ); 		// Yn+1/2 <-- Yn
	gsl_vector_mul( Ynph, DiagonalAn ); // Yn+1/2 <-- Dn Yn
	// 3 -  dt/2 An Yn == dt/2 (An-Dn) Yn + dt/2 Dn Yn
	gsl_blas_dgemv( CblasNoTrans, hh2, OffDiagonalAn, Yn, hh2, Ynph );
	// 4 -  add Yn  i.e. Yn+1/2 <-- dt/2 An Yn  +   Yn
	gsl_vector_add( Ynph, Yn );


	// Set Yn+1^0 to be Yn   -- check if setting it to Yn+1/2 reduces the number of iterations required
	gsl_vector_memcpy( Ynp1i, Yn );

	int iter = 0;
	double err_sum = 0;
	do {
		// Set Yn+1 to be Ynph
		gsl_vector_memcpy( Ynp1, Ynph );

		// Yn+1^i+1 = [ I - dt/2 Dn ]^-1 ( Yn+1/2 + dt/2 A~n Yn+1^i )
		gsl_blas_dgemv( CblasNoTrans, hh2, OffDiagonalAn, Ynp1i, 1.0, Ynp1 ); // Yn+1 <-- Yn+1/2 + dt/2 A~n Yn+1^i
		gsl_vector_mul( Ynp1, DiagonalAnInv ); // Yn+1 <-- [I - dt/2 Dn ]^-1 Yn+1

		// find difference of iterates Ynp1i and Ynp1
		gsl_vector_sub( Ynp1i, Ynp1 );
		err_sum = gsl_blas_dasum( Ynp1i ); // double abs value sum

		gsl_vector_memcpy( Ynp1i, Ynp1 );
		iter++;
	} while ( err_sum > 1E-100 && iter < 100 );

	MarkovIterations = iter;
	MarkovIterationError = err_sum;

	// Once the BLAS code is working, use the following loop
	for (unsigned int i=0; i<nNaChStates; i++)
		state[i] = gsl_vector_get( Ynp1, i );

	checkNaChannelStates();

	I_Na = G_Na*(O)*(V - E_Na);
}

void Moreno2011Cell::set_NaCh_Matrices( gsl_matrix* M, gsl_vector* D,
		const double a11,	const double a12,	const double a13,
		const double b11,	const double b12,	const double b13,
		const double a3,	const double b3,	const double a2, const double b2,
		const double a4,	const double b4,	const double a5, const double b5,
		const double ax,	const double bx,	const double ax1,const double bx1,
		const double a13c,  const double a22,	const double b33,const double a33,
		const double a44,   const double b44,	const double a55,const double b55,
		const double ax2,   const double a13n,  const double a_22, const double b_33,
		const double a_44,	const double b_44,  const double a_55, const double b_55,
		const double kon,   const double koff,  const double kcoff, const double kcon,
		const double b13c,	const double b22,	const double k_on,	const double k_off,
		const double ki_on,	const double ki_off,const double kc_on, const double kc_off,
		const double a_33,	const double b13n,	const double b_22,	const double bx2 ){

	double DiagAn[nNaChStates];
	// The following 30 consts are the diagonal coefficients of the A matrix
	const double &coef_IC3	= DiagAn[ 0] = (a11 + a3 + ki_on);
	const double &coef_IC2	= DiagAn[ 1] = (b11 + a3 + a12 + ki_on);
	const double &coef_IF	= DiagAn[ 2] = (b12 + b2 + a3 + a4 + ki_on);
	const double &coef_IM1	= DiagAn[ 3] = (b4 + a5);
	const double &coef_IM2	= DiagAn[ 4] = b5;
	const double &coef_C3	= DiagAn[ 5] = (b3 + a11  + kcon + kc_on);
	const double &coef_C2	= DiagAn[ 6] = (b11 + b3 + a12  + kcon + kc_on);
	const double &coef_O	= DiagAn[ 7] = (b13 + a2  + kon + k_on + ax);
	const double &coef_OS	= DiagAn[ 8] = (bx + ki_on);
	const double &coef_C1	= DiagAn[ 9] = (b12 + b3 + a13  + kcon + kc_on);


	const double &coef_DIC3	= DiagAn[10] = (a11 + a33);
	const double &coef_DIC2 = DiagAn[11] = (a33 + b11 + a12);
	const double &coef_DIF	= DiagAn[12] = (a33 + b12 + a44 + b22);
	const double &coef_DIM1	= DiagAn[13] = ( b44 + a55 );
	const double &coef_DIM2	= DiagAn[14] = b55 ;
	const double &coef_DC3	= DiagAn[15] = (kcoff+ b33 + a11 );
	const double &coef_DC2	= DiagAn[16] = (kcoff + b11 + b33 + a12 );
	const double &coef_DO	= DiagAn[17] = (koff + b13c + a22 + ax1 );
	const double &coef_DOS	= DiagAn[18] = bx1;
	const double &coef_DC1	= DiagAn[19] = (kcoff + b12 + b33 + a13c );

	const double &coef_D_IC3= DiagAn[20] = (a_33 + a11 + ki_off);
	const double &coef_D_IC2= DiagAn[21] = (a_33 + b11 + a12 + ki_off);
	const double &coef_D_IF	= DiagAn[22] = (a_33 + a_44 + b_22 + b12 + ki_off);
	const double &coef_D_IM1= DiagAn[23] = (b_44 + a_55);
	const double &coef_D_IM2= DiagAn[24] = b_55;
	const double &coef_D_C3	= DiagAn[25] = (kc_off + b_33 + a11 );
	const double &coef_D_C2	= DiagAn[26] = (kc_off + b11 + b_33 + a12 );
	const double &coef_D_O	= DiagAn[27] = (k_off + b13n + a_22 + ax2 );
	const double &coef_D_OS	= DiagAn[28] = (bx2 + ki_off);
	const double &coef_D_C1	= DiagAn[29] = (kc_off + b12 + b_33 + a13n );

	// Set up diagonal coefficients
	gsl_vector_set_all( D, 0.0 );
	gsl_vector_set( D,  i3, -coef_IC3 );
	gsl_vector_set( D,  i2, -coef_IC2 );
	gsl_vector_set( D,  iF, -coef_IF  );
	gsl_vector_set( D,  im, -coef_IM1 );
	gsl_vector_set( D,  iM, -coef_IM2 );
	gsl_vector_set( D,  c3, -coef_C3  );
	gsl_vector_set( D,  c2, -coef_C2  );
	gsl_vector_set( D,  op, -coef_O	  );
	gsl_vector_set( D,  oS, -coef_OS  );
	gsl_vector_set( D,  c1, -coef_C1  );

	gsl_vector_set( D, di3, -coef_DIC3 );
	gsl_vector_set( D, di2, -coef_DIC2 );
	gsl_vector_set( D, diF, -coef_DIF  );
	gsl_vector_set( D, dim, -coef_DIM1 );
	gsl_vector_set( D, diM, -coef_DIM2 );
	gsl_vector_set( D, dc3, -coef_DC3  );
	gsl_vector_set( D, dc2, -coef_DC2  );
	gsl_vector_set( D, dop, -coef_DO   );
	gsl_vector_set( D, doS, -coef_DOS  );
	gsl_vector_set( D, dc1, -coef_DC1  );

	gsl_vector_set( D, Di3, -coef_D_IC3 );
	gsl_vector_set( D, Di2, -coef_D_IC2 );
	gsl_vector_set( D, DiF, -coef_D_IF  );
	gsl_vector_set( D, Dim, -coef_D_IM1 );
	gsl_vector_set( D, DiM, -coef_D_IM2 );
	gsl_vector_set( D, Dc3, -coef_D_C3	);
	gsl_vector_set( D, Dc2, -coef_D_C2	);
	gsl_vector_set( D, Dop, -coef_D_O	);
	gsl_vector_set( D, DoS, -coef_D_OS	);
	gsl_vector_set( D, Dc1, -coef_D_C1	);


	gsl_matrix_set_all( M, 0 );
	// Drug-free states
	gsl_matrix_set( M,  i3,   c3, b3     ); // Row i3
	gsl_matrix_set( M,  i3,   i2, b11    );
	gsl_matrix_set( M,  i3,  Di3, ki_off );
	gsl_matrix_set( M,  i2,   i3, a11    ); // Row i2
	gsl_matrix_set( M,  i2,   c2, b3     );
	gsl_matrix_set( M,  i2,   iF, b12    );
	gsl_matrix_set( M,  i2,  Di2, ki_off );
	gsl_matrix_set( M,  iF,   i2, a12    ); // Row iF
	gsl_matrix_set( M,  iF,   c1, b3     );
	gsl_matrix_set( M,  iF,   im, b4     );
	gsl_matrix_set( M,  iF,   op, a2	 );
	gsl_matrix_set( M,  iF,  DiF, ki_off );
	gsl_matrix_set( M,  im,   iF, a4	 ); // Row im
	gsl_matrix_set( M,  im,   iM, b5	 );
	gsl_matrix_set( M,  iM,   im, a5	 ); // Row iM
	gsl_matrix_set( M,  c3,   i3, a3	 ); // Row c3
	gsl_matrix_set( M,  c3,   c2, b11	 );
	gsl_matrix_set( M,  c3,  dc3, kcoff  );
	gsl_matrix_set( M,  c3,  Dc3, kc_off );
	gsl_matrix_set( M,  c2,   c3, a11	 ); // Row c2
	gsl_matrix_set( M,  c2,   i2, a3	 );
	gsl_matrix_set( M,  c2,   c1, b12	 );
	gsl_matrix_set( M,  c2,  dc2, kcoff  );
	gsl_matrix_set( M,  c2,  Dc2, kc_off );
	gsl_matrix_set( M,  op,   c1, a13	 ); // Row op
	gsl_matrix_set( M,  op,   iF, b2	 );
	gsl_matrix_set( M,  op,  dop, koff	 );
	gsl_matrix_set( M,  op,  Dop, k_off	 );
	gsl_matrix_set( M,  op,   oS, bx	 );
	gsl_matrix_set( M,  oS,   op, ax	 ); // Row oS
	gsl_matrix_set( M,  oS,  DoS, ki_off );
	gsl_matrix_set( M,  c1,   c2, a12	 ); // Row c1
	gsl_matrix_set( M,  c1,   iF, a3	 );
	gsl_matrix_set( M,  c1,   op, b13	 );
	gsl_matrix_set( M,  c1,  dc1, kcoff  );
	gsl_matrix_set( M,  c1,  Dc1, kc_off );

	// Charged-drug states
	gsl_matrix_set( M,  di3,  dc3, b33 ); // Row di3
	gsl_matrix_set( M,  di3,  di2, b11 );
	gsl_matrix_set( M,  di2,  di3, a11 ); // Row di2
	gsl_matrix_set( M,  di2,  dc2, b33 );
	gsl_matrix_set( M,  di2,  diF, b12 );
	gsl_matrix_set( M,  diF,  di2, a12 ); // Row diF
	gsl_matrix_set( M,  diF,  dc1, b33 );
	gsl_matrix_set( M,  diF,  dim, b44 );
	gsl_matrix_set( M,  diF,  dop, a22 );
	gsl_matrix_set( M,  dim,  diF, a44 ); // Row dim
	gsl_matrix_set( M,  dim,  diM, b55 );
	gsl_matrix_set( M,  diM,  dim, a55 ); // Row diM
	gsl_matrix_set( M,  dc3,  di3, a33 ); // Row dc3
	gsl_matrix_set( M,  dc3,  dc2, b11 );
	gsl_matrix_set( M,  dc3,   c3, kcon);
	gsl_matrix_set( M,  dc2,  dc3, a11 ); // Row dc2
	gsl_matrix_set( M,  dc2,  di2, a33 );
	gsl_matrix_set( M,  dc2,  dc1, b12 );
	gsl_matrix_set( M,  dc2,   c2, kcon);
	gsl_matrix_set( M,  dop,  dc1, a13c); // Row dop
	gsl_matrix_set( M,  dop,  diF, b22 );
	gsl_matrix_set( M,  dop,   op, kon );
	gsl_matrix_set( M,  dop,  doS, bx1 );
	gsl_matrix_set( M,  doS,  dop, ax1 ); // Row doS
	gsl_matrix_set( M,  dc1,  dc2, a12 ); // Row dc1
	gsl_matrix_set( M,  dc1,  diF, a33 );
	gsl_matrix_set( M,  dc1,  dop, b13c);
	gsl_matrix_set( M,  dc1,   c1, kcon);

	// Neutral-drug states
	gsl_matrix_set( M,  Di3,  Dc3, b_33 ); // Row Di3
	gsl_matrix_set( M,  Di3,  Di2, b11  );
	gsl_matrix_set( M,  Di3,   i3, ki_on);
	gsl_matrix_set( M,  Di2,  Di3, a11  ); // Row Di2
	gsl_matrix_set( M,  Di2,  Dc2, b_33 );
	gsl_matrix_set( M,  Di2,  DiF, b12  );
	gsl_matrix_set( M,  Di2,   i2, ki_on);
	gsl_matrix_set( M,  DiF,  Di2, a12  ); // Row DiF
	gsl_matrix_set( M,  DiF,  Dc1, b_33 );
	gsl_matrix_set( M,  DiF,  Dim, b_44 );
	gsl_matrix_set( M,  DiF,  Dop, a_22 );
	gsl_matrix_set( M,  DiF,   iF, ki_on);
	gsl_matrix_set( M,  Dim,  DiF, a_44 ); // Row Dim
	gsl_matrix_set( M,  Dim,  DiM, b_55 );
	gsl_matrix_set( M,  DiM,  Dim, a_55 ); // Row DiM
	gsl_matrix_set( M,  Dc3,  Di3, a_33 ); // Row Dc3
	gsl_matrix_set( M,  Dc3,  Dc2, b11  );
	gsl_matrix_set( M,  Dc3,   c3, kc_on);
	gsl_matrix_set( M,  Dc2,  Dc3, a11  ); // Row Dc2
	gsl_matrix_set( M,  Dc2,  Di2, a_33 );
	gsl_matrix_set( M,  Dc2,  Dc1, b12  );
	gsl_matrix_set( M,  Dc2,   c2, kc_on);
	gsl_matrix_set( M,  Dop,  Dc1, a13n ); // Row Dop
	gsl_matrix_set( M,  Dop,  DiF, b_22 );
	gsl_matrix_set( M,  Dop,   op, k_on );
	gsl_matrix_set( M,  Dop,  DoS, bx2  );
	gsl_matrix_set( M,  DoS,  Dop, ax2  ); // Row DoS
	gsl_matrix_set( M,  DoS,   oS, ki_on);
	gsl_matrix_set( M,  Dc1,  Dc2, a12  ); // Row Dc1
	gsl_matrix_set( M,  Dc1,  DiF, a_33 );
	gsl_matrix_set( M,  Dc1,  Dop, b13n );
	gsl_matrix_set( M,  Dc1,   c1, kc_on);
}

void Moreno2011Cell::semi_impl_WT_SCN5A_Lidocaine_compare_BLAS_old(const double dt){
	const double pKa	= 7.6;
	const double diffusion	= 500;
	const double portion	= 1/(1+ pow(10, (pH-pKa)) );
	const double drug_charged  = cellenv->drug*portion;
	const double drug_neutral  = cellenv->drug*(1-portion);
	const double drug_distance = -0.7;

	const double E_Na = (R*T/F)*log(Na_out/Na_in);

	//Rate Constants **********************************************************
	const double Tfactor = 1.0/(pow(Q10, (37.0-(T-273))/10.0));

	//WT Fits Reduced Model (no IM1, IM2)
	const double a11 = Tfactor*8.5539/(7.4392e-2*exp(-V/17.0)+ 2.0373e-1*exp(-V/150));
	const double a12 = Tfactor*8.5539/(7.4392e-2*exp(-V/15.0)+ 2.0373e-1*exp(-V/150));
	const double a13 = Tfactor*8.5539/(7.4392e-2*exp(-V/12.0)+ 2.0373e-1*exp(-V/150));
	const double b11 = Tfactor*7.5215e-2*exp(-V/20.3);
	const double b12 = Tfactor*2.7574*exp(-(V-5)/20.3);
	const double b13 = Tfactor*4.7755e-1*exp(-(V-10)/20.3);

	const double a3	= Tfactor*5.1458e-6*exp(-V/8.2471);
	const double b3	= Tfactor*6.1205*exp((V)/13.542);

	const double a2	= Tfactor*(13.370*exp(V/43.749));
	const double b2	= a13*a2*a3/(b13*b3);

	const double a4 = 0*a2;
	const double b4 = 0*a3;
	const double a5	= 0*a2;
	const double b5 = 0*a3;

	const double ax = 3.4229e-2*a2;
	const double bx = 1.7898e-2*a3;


	const double ax1 = 6.3992e-07  *ax;
	const double bx1 = 1.3511e+00  *bx;
	const double a13c= 5.6974e-03  *a13;
	const double a22 = 6.7067e-06  *a2;
	const double b33 =  1.9698e-05 *b3;
	const double a33 = 3.2976e+00  *a3;

	const double a44 = 0  *a2;
	const double b44 = 0  *a3;
	const double a55 = 0;
	const double b55 = 0;

	const double ax2  = 1.3110e-01  *ax;
	const double a13n = 8.4559e+01  *a13;
	const double a_22 = 1.7084e-05  *a2;
	const double b_33 = 4.8477e+00  *b3;

	const double a_44 = 0  *a2;
	const double b_44 = 0  *a3;
	const double a_55 = 0;
	const double b_55 = 0;


	const double kd0	= 318e-6;
	const double kd_open= kd0*exp( (drug_distance*V*F) /(R*T));

	// charged drug
	const double kon	= drug_charged*diffusion;
	const double koff	= kd_open*diffusion;
	const double kcoff	= koff;
	const double kcon	= kon;

	const double b13c	= (cellenv->drug ==0 || drug_charged ==0 ) ? 0 : (b13*kcon*koff*a13c)/(kon*kcoff*a13);
	const double b22	= (b13c ==0) ? 0 : (a13c*a22*a33)/(b13c*b33);

	// neutral drug
	const double k_on	= drug_neutral*diffusion;
	const double k_off	= 400e-6*diffusion;
	const double ki_on	= k_on/2;
	const double ki_off	= 3.4e-6*diffusion;
	const double kc_on	= k_on/2;
	const double kc_off	= 900e-6*diffusion;

	const double a_33 = (cellenv->drug ==0 || drug_neutral ==0 ) ? 0 : (ki_off*a3*kc_on*b_33)/(ki_on*kc_off*b3);
	const double b13n = (cellenv->drug ==0 || drug_neutral ==0 ) ? 0 : (b13*kc_on*a13n*k_off)/(kc_off*a13*k_on);
	const double b_22 = (b13n == 0) ? 0 : (a_33*a13n*a_22)/(b_33*b13n);
	const double bx2  = (cellenv->drug ==0 || drug_neutral ==0) ? 0 : (bx*k_on*ax2*ki_off)/(ax*ki_on*k_off);



	//		// Calculation of k parameters for the Runge Kutta interations
	//
	//		//Calculation of K1 ***********************************************
	//		//Drug Free States
	//		O = hh*( a13 * C1 + b2 * IF + koff * DO + k_off * D_O + bx * OS  - O * (b13 + a2  + kon + k_on + ax));
	//		C1 = hh*(a12 * C2 + a3 * IF + b13 * O  + kcoff * DC1 + kc_off * D_C1 - C1*(b12 + b3 + a13  + kcon + kc_on));
	//		C2 = hh*(a11 * C3 + a3 * IC2 + b12 * C1  + kcoff * DC2 + kc_off * D_C2 - C2*(b11 + b3 + a12  + kcon + kc_on));
	//		C3 = hh*(a3 * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 - C3*(b3 + a11  + kcon + kc_on));
	//		IC3 = hh*(b3 * C3 + b11 * IC2 + ki_off * D_IC3 - IC3*(a11 + a3 + ki_on));
	//		IC2 = hh*(a11 * IC3 + b3 * C2 + b12 * IF + ki_off * D_IC2 - IC2*(b11 + a3 + a12 + ki_on));
	//		IF = hh*(a12 * IC2 + b3 * C1 + b4 * IM1 + a2 * O + ki_off * D_IF - IF*(b12 + b2 + a3 + a4 + ki_on));
	//		IM1 = hh*(a4 * IF + b5 * IM2 - IM1*(b4 + a5));
	//		IM2 = hh*(a5 * IM1 - IM2*(b5));
	//		OS = hh*(ax * O + ki_off * D_OS - OS*(bx + ki_on));
	//
	//		//Charged Drug Bound States
	//		DO = hh*(kon * O + a13c * DC1 + bx1 * DOS + b22 * DIF  - DO*(koff + b13c + a22 + ax1 ));
	//		DC1 = hh*(kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO  - DC1*(kcoff + b12 + b33 + a13c ));
	//		DC2 = hh*(kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1  - DC2*(kcoff + b11 + b33 + a12 ));
	//		DC3 = hh*(kcon * C3 + b11 * DC2 + a33 * DIC3  - DC3*(kcoff+ b33 + a11 ));
	//		DOS = hh*(ax1 * DO - bx1 * DOS);
	//		DIC3 = hh*(b33 * DC3 + b11 * DIC2 - DIC3*(a11 + a33));
	//		DIC2 = hh*(b33 * DC2 + a11 * DIC3 + b12 * DIF - DIC2*(a33 + b11 + a12));
	//		DIF = hh*(b33 * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO - DIF*(a33 + b12 + a44 + b22));
	//		DIM1 = hh*(a44 * DIF + b55 * DIM2 - DIM1*( b44 + a55));
	//		DIM2 = hh*(a55 * DIM1 - b55 * DIM2);
	//
	//		//Neutral Drug Bound States
	//		D_O = hh*(k_on * O + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS  - D_O*(k_off + b13n + a_22 + ax2 ));
	//		D_C1 = hh*(kc_on * C1 + a12 * D_C2 + a_33 * D_IF + b13n * D_O  - D_C1*(kc_off + b12 + b_33 + a13n ));
	//		D_C2 = hh*(kc_on * C2 + a11 * D_C3 + a_33 * D_IC2 + b12 * D_C1  - D_C2*(kc_off + b11 + b_33 + a12 ));
	//		D_C3 = hh*(kc_on * C3 + a_33 * D_IC3 + b11 * D_C2  - D_C3*(kc_off + b_33 + a11 ));
	//		D_OS = hh*(ax2 * D_O + ki_on * OS - D_OS*(bx2 + ki_off));
	//		D_IC3 = hh*(b_33 * D_C3 + b11 * D_IC2 + ki_on * IC3 - D_IC3*(a_33 + a11 + ki_off));
	//		D_IC2 = hh*(b_33 * D_C2 + a11 * D_IC3 + b12 * D_IF + ki_on * IC2 - D_IC2*(a_33 + b11 + a12 + ki_off));
	//		D_IF = hh*(b_33 * D_C1 + a12 * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF - D_IF*(a_33 + a_44 + b_22 + b12 + ki_off));
	//		D_IM1 = hh*(a_44 * D_IF + b_55 * D_IM2 - D_IM1*(b_44 + a_55));
	//		D_IM2 = hh*(a_55 * D_IM1 - b_55 * D_IM2);




	const double hh	 = dt;
	const double hh2 = hh/2;

	// The following 30 constants are the negative of the Diagonal of An
	// These can be removed when the BLAS version is working
	const double coef_IC3	= (a11 + a3 + ki_on);
	const double coef_IC2	= (b11 + a3 + a12 + ki_on);
	const double coef_IF	= (b12 + b2 + a3 + a4 + ki_on);
	const double coef_IM1	= (b4 + a5);
	const double coef_IM2	= b5;
	const double coef_C3	= (b3 + a11  + kcon + kc_on);
	const double coef_C2	= (b11 + b3 + a12  + kcon + kc_on);
	const double coef_O		= (b13 + a2  + kon + k_on + ax);
	const double coef_OS	= (bx + ki_on);
	const double coef_C1	= (b12 + b3 + a13  + kcon + kc_on);
	const double coef_DIC3	= (a11 + a33);
	const double coef_DIC2	= (a33 + b11 + a12);
	const double coef_DIF	= (a33 + b12 + a44 + b22);
	const double coef_DIM1	= ( b44 + a55 );
	const double coef_DIM2	= b55 ;
	const double coef_DC3	= (kcoff+ b33 + a11 );
	const double coef_DC2	= (kcoff + b11 + b33 + a12 );
	const double coef_DO	= (koff + b13c + a22 + ax1 );
	const double coef_DOS	= bx1;
	const double coef_DC1	= (kcoff + b12 + b33 + a13c );
	const double coef_D_IC3	= (a_33 + a11 + ki_off);
	const double coef_D_IC2	= (a_33 + b11 + a12 + ki_off);
	const double coef_D_IF	= (a_33 + a_44 + b_22 + b12 + ki_off);
	const double coef_D_IM1	= (b_44 + a_55);
	const double coef_D_IM2	= b_55;
	const double coef_D_C3	= (kc_off + b_33 + a11 );
	const double coef_D_C2	= (kc_off + b11 + b_33 + a12 );
	const double coef_D_O	= (k_off + b13n + a_22 + ax2 );
	const double coef_D_OS	= (bx2 + ki_off);
	const double coef_D_C1	= (kc_off + b12 + b_33 + a13n );

	// The following 30 consts are related to Dn the diagonal coefficients of the An matrix
	// Can be removed for the BLAS version
	const double co_IC3   = 1. / ( 1 + hh2 * ( a11 + a3 + ki_on));
	const double co_IC2   = 1. / ( 1 + hh2 * ( b11 + a3 + a12 + ki_on));
	const double co_IF	  = 1. / ( 1 + hh2 * ( b12 + b2 + a3 + a4 + ki_on));
	const double co_IM1   = 1. / ( 1 + hh2 * ( b4 + a5));
	const double co_IM2   = 1. / ( 1 + hh2 * ( b5 ));
	const double co_C3	  = 1. / ( 1 + hh2 * ( b3 + a11  + kcon + kc_on));
	const double co_C2	  = 1. / ( 1 + hh2 * ( b11 + b3 + a12  + kcon + kc_on));
	const double co_O	  = 1. / ( 1 + hh2 * ( b13 + a2  + kon + k_on + ax ) );
	const double co_OS	  = 1. / ( 1 + hh2 * ( bx + ki_on));
	const double co_C1	  = 1. / ( 1 + hh2 * ( b12 + b3 + a13  + kcon + kc_on));
	const double co_DIC3  = 1. / ( 1 + hh2 * ( a11 + a33));
	const double co_DIC2  = 1. / ( 1 + hh2 * ( a33 + b11 + a12));
	const double co_DIF   = 1. / ( 1 + hh2 * ( a33 + b12 + a44 + b22));
	const double co_DIM1  = 1. / ( 1 + hh2 * ( b44 + a55 ));
	const double co_DIM2  = 1. / ( 1 + hh2 * ( b55 ) );
	const double co_DC3	  = 1. / ( 1 + hh2 * ( kcoff+ b33 + a11 ));
	const double co_DC2   = 1. / ( 1 + hh2 * ( kcoff + b11 + b33 + a12 ));
	const double co_DO    = 1. / ( 1 + hh2 * ( koff + b13c + a22 + ax1 ));
	const double co_DOS	  = 1. / ( 1 + hh2 * ( bx1));
	const double co_DC1   = 1. / ( 1 + hh2 * ( kcoff + b12 + b33 + a13c ));
	const double co_D_IC3 = 1. / ( 1 + hh2 * ( a_33 + a11 + ki_off));
	const double co_D_IC2 = 1. / ( 1 + hh2 * ( a_33 + b11 + a12 + ki_off));
	const double co_D_IF  = 1. / ( 1 + hh2 * ( a_33 + a_44 + b_22 + b12 + ki_off));
	const double co_D_IM1 = 1. / ( 1 + hh2 * ( b_44 + a_55));
	const double co_D_IM2 = 1. / ( 1 + hh2 * ( b_55 ));
	const double co_D_C3  = 1. / ( 1 + hh2 * ( kc_off + b_33 + a11 ));
	const double co_D_C2  = 1. / ( 1 + hh2 * ( kc_off + b11 + b_33 + a12 ));
	const double co_D_O   = 1. / ( 1 + hh2 * ( k_off + b13n + a_22 + ax2 ));
	const double co_D_OS  = 1. / ( 1 + hh2 * ( bx2 + ki_off));
	const double co_D_C1  = 1. / ( 1 + hh2 * ( kc_off + b12 + b_33 + a13n ));
	double MCYn[nNaChStates];
	double &IC3_o  = MCYn[ 0], &IC2_o  = MCYn[ 1], &IF_o  = MCYn[ 2], &IM1_o  = MCYn[ 3], &IM2_o  = MCYn[ 4], &C3_o  = MCYn[ 5], &C2_o  = MCYn[ 6], &O_o  = MCYn[ 7], &OS_o  = MCYn[ 8], &C1_o  = MCYn[ 9];
	double &DIC3_o = MCYn[10], &DIC2_o = MCYn[11], &DIF_o = MCYn[12], &DIM1_o = MCYn[13], &DIM2_o = MCYn[14], &DC3_o = MCYn[15], &DC2_o = MCYn[16], &DO_o = MCYn[17], &DOS_o = MCYn[18], &DC1_o = MCYn[19];
	double &D_IC3_o= MCYn[20], &D_IC2_o= MCYn[21], &D_IF_o= MCYn[22], &D_IM1_o= MCYn[23], &D_IM2_o= MCYn[24], &D_C3_o= MCYn[25], &D_C2_o= MCYn[26], &D_O_o= MCYn[27], &D_OS_o= MCYn[28], &D_C1_o= MCYn[29];

	set_NaCh_Matrices( OffDiagonalAn, DiagonalAn,
			a11,  a12,	a13,
			b11,  b12,	b13,
			a3,	  b3,	a2,    b2,
			a4,	  b4,	a5,	   b5,
			ax,	  bx,	ax1,   bx1,
			a13c, a22,	b33,   a33,
			a44,  b44,	a55,   b55,
			ax2,  a13n, a_22,  b_33,
			a_44, b_44, a_55,  b_55,
			kon,  koff, kcoff, kcon,
			b13c, b22,  k_on,  k_off,
			ki_on,ki_off, kc_on, kc_off,
			a_33, b13n, b_22,  bx2 );
	//gsl_matrix_fprintf( stdout, OffDiagonalAn, "%f " );

	// DiagonalAnInv = [I - dt/2 Dn ]^-1
	// The following loop should be equivalent to the 30 lines of code setting co_ constants
	for (unsigned int i=0; i<nNaChStates; i++)
		gsl_vector_set( DiagonalAnInv, i, 1.0/( 1.0 - hh2 * gsl_vector_get( DiagonalAn, i )) );


	// Compute Yn+1/2 = Yn + dt/2 An Yn
	//                = Yn + dt/2 (An - Dn) Yn + dt/2 Dn Yn
	// 1 -   Set Yn
	for (unsigned int i=0; i<nNaChStates; i++)
		gsl_vector_set( Yn, i, state[i] );
	// 2 -   Ynph <-- Dn Yn
	gsl_vector_memcpy( Ynph, Yn ); 		// Yn+1/2 <-- Yn
	gsl_vector_mul( Ynph, DiagonalAn ); // Yn+1/2 <-- Dn Yn
	// 3 -  dt/2 An Yn == dt/2 (An-Dn) Yn + dt/2 Dn Yn
	gsl_blas_dgemv( CblasNoTrans, hh2, OffDiagonalAn, Yn, hh2, Ynph );
	// 4 -  add Yn  i.e. Yn+1/2 <-- dt/2 An Yn  +   Yn
	gsl_vector_add( Ynph, Yn );
	//	gsl_vector_fprintf( stdout, Ynph, "%5.4f " );

	//Drug Free States
	IC3_o	= IC3 + hh2 * (b3  * C3  + b11 * IC2 + ki_off * D_IC3 - IC3 * coef_IC3 );
	IC2_o	= IC2 + hh2 * (a11 * IC3 + b3  * C2  + b12  * IF  + ki_off * D_IC2 - IC2 * coef_IC2 );
	IF_o 	=  IF + hh2 * (a12 * IC2 + b3  * C1  + b4   * IM1 + a2 * O + ki_off * D_IF - IF * coef_IF );
	IM1_o	= IM1 + hh2 * (a4  * IF  + b5  * IM2 - IM1  * coef_IM1 );
	IM2_o	= IM2 + hh2 * (a5  * IM1 - IM2 * coef_IM2 );
	C3_o 	=  C3 + hh2 * (a3  * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 - C3 * coef_C3 );
	C2_o 	=  C2 + hh2 * (a11 * C3  + a3  * IC2 + b12  * C1  + kcoff * DC2 + kc_off * D_C2 - C2 * coef_C2 );
	O_o  	=   O + hh2 * (a13 * C1  + b2  * IF  + koff * DO  + k_off * D_O + bx * OS  - O * coef_O );
	OS_o 	=  OS + hh2 * (ax  * O   + ki_off * D_OS - OS * coef_OS );
	C1_o 	=  C1 + hh2 * (a12 * C2  + a3  * IF  + b13  * O   + kcoff * DC1 + kc_off * D_C1 - C1 * coef_C1 );

	//Charged Drug Bound States
	DIC3_o	= DIC3 + hh2 * (b33 * DC3 + b11 * DIC2 - DIC3 * coef_DIC3 );
	DIC2_o	= DIC2 + hh2 * (b33 * DC2 + a11 * DIC3 + b12 * DIF - DIC2 * coef_DIC2 );
	DIF_o	=  DIF + hh2 * (b33 * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO - DIF * coef_DIF );
	DIM1_o	= DIM1 + hh2 * (a44 * DIF + b55 * DIM2 - DIM1 * coef_DIM1 );
	DIM2_o	= DIM2 + hh2 * (a55 * DIM1 - DIM2 * coef_DIM2 );
	DC3_o	=  DC3 + hh2 * (kcon * C3 + b11 * DC2 + a33 * DIC3  - DC3 * coef_DC3 );
	DC2_o	=  DC2 + hh2 * (kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1  - DC2 * coef_DC2 );
	DO_o	=   DO + hh2 * (kon * O + a13c * DC1 + bx1 * DOS + b22 * DIF  - DO * coef_DO );
	DOS_o	=  DOS + hh2 * (ax1 * DO -  DOS * coef_DOS );
	DC1_o	=  DC1 + hh2 * (kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO  - DC1 * coef_DC1 );

	//Neutral Drug Bound States
	D_IC3_o	= D_IC3 + hh2 * (b_33  * D_C3  + b11  * D_IC2 + ki_on * IC3 - D_IC3 * coef_D_IC3 );
	D_IC2_o	= D_IC2 + hh2 * (b_33  * D_C2  + a11  * D_IC3 + b12 * D_IF + ki_on * IC2 - D_IC2 * coef_D_IC2 );
	D_IF_o	=  D_IF + hh2 * (b_33  * D_C1  + a12  * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF - D_IF * coef_D_IF );
	D_IM1_o	= D_IM1 + hh2 * (a_44  * D_IF  + b_55 * D_IM2 - D_IM1 * coef_D_IM1 );
	D_IM2_o	= D_IM2 + hh2 * (a_55  * D_IM1 - D_IM2 * coef_D_IM2 );
	D_C3_o	=  D_C3 + hh2 * (kc_on * C3    + a_33 * D_IC3 + b11 * D_C2  - D_C3 * coef_D_C3 );
	D_C2_o	=  D_C2 + hh2 * (kc_on * C2    + a11  * D_C3 + a_33 * D_IC2 + b12 * D_C1  - D_C2 * coef_D_C2 );
	D_O_o	=   D_O + hh2 * (k_on  * O     + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS  - D_O * coef_D_O );
	D_OS_o	=  D_OS + hh2 * (ax2   * D_O   + ki_on * OS - D_OS * coef_D_OS );
	D_C1_o	=  D_C1 + hh2 * (kc_on * C1    + a12  * D_C2 + a_33 * D_IF + b13n * D_O  - D_C1 * coef_D_C1 );

	// Compare the BLAS version and the hand-calculated version
	double difference = 0;
	for (unsigned int i=0; i<nNaChStates; i++)
		difference += fabs( MCYn[i] - gsl_vector_get( Ynph, i ) );
	// std::cout << difference << std::endl;

	// Set Yn+1^0 to be Yn   -- check if setting it to Yn+1/2 reduces the number of iterations required
	gsl_vector_memcpy( Ynp1i, Yn );

	int iter = 0;
	double err_sum = 1;
	double err_sum_blas = 1;
	while ( err_sum_blas > 1E-100 && iter < 100 ) {
		// Set Yn+1 to be Ynph
		gsl_vector_memcpy( Ynp1, Ynph );

		// Yn+1^i+1 = [ I - dt/2 Dn ]^-1 ( Yn+1/2 + dt/2 A~n Yn+1^i )
		gsl_blas_dgemv( CblasNoTrans, hh2, OffDiagonalAn, Ynp1i, 1.0, Ynp1 ); // Yn+1 <-- Yn+1/2 + dt/2 A~n Yn+1^i
		gsl_vector_mul( Ynp1, DiagonalAnInv ); // Yn+1 <-- [I - dt/2 Dn ]^-1 Yn+1

		// find difference of iterates Ynp1i and Ynp1
		gsl_vector_sub( Ynp1i, Ynp1 );
		err_sum_blas = gsl_blas_dasum( Ynp1i ); // double abs value sum

		gsl_vector_memcpy( Ynp1i, Ynp1 );


		//Drug Free States
		mO   = co_O  * ( O_o	+ hh2 * (a13 * C1 + b2 * IF + koff * DO + k_off * D_O + bx * OS ) );
		mC1  = co_C1 * ( C1_o	+ hh2 * (a12 * C2 + a3 * IF + b13 * O  + kcoff * DC1 + kc_off * D_C1 ) );
		mC2  = co_C2 * ( C2_o	+ hh2 * (a11 * C3 + a3 * IC2 + b12 * C1  + kcoff * DC2 + kc_off * D_C2 ) );
		mC3  = co_C3 * ( C3_o	+ hh2 * (a3 * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 ) );
		mIC3 = co_IC3* ( IC3_o	+ hh2 * (b3 * C3 + b11 * IC2 + ki_off * D_IC3  ) );
		mIC2 = co_IC2* ( IC2_o	+ hh2 * (a11 * IC3 + b3 * C2 + b12 * IF + ki_off * D_IC2  ) );
		mIF  = co_IF * ( IF_o	+ hh2 * (a12 * IC2 + b3 * C1 + b4 * IM1 + a2 * O + ki_off * D_IF  ) );
		mIM1 = co_IM1* ( IM1_o	+ hh2 * (a4 * IF + b5 * IM2  ) );
		mIM2 = co_IM2* ( IM2_o	+ hh2 * (a5 * IM1  ) );
		mOS  = co_OS * ( OS_o	+ hh2 * (ax * O + ki_off * D_OS  ) );

		//Charged Drug Bound States
		mDO	 = co_DO  * ( DO_o	+ hh2 * (kon * O + a13c * DC1 + bx1 * DOS + b22 * DIF   ) );
		mDC1 = co_DC1 * ( DC1_o	+ hh2 * (kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO  ) );
		mDC2 = co_DC2 * ( DC2_o	+ hh2 * (kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1   ) );
		mDC3 = co_DC3 * ( DC3_o	+ hh2 * (kcon * C3 + b11 * DC2 + a33 * DIC3   ) );
		mDOS = co_DOS * ( DOS_o	+ hh2 * (ax1 * DO  ) );
		mDIC3= co_DIC3* ( DIC3_o + hh2 * (b33 * DC3 + b11 * DIC2 ) );
		mDIC2= co_DIC2* ( DIC2_o + hh2 * (b33 * DC2 + a11 * DIC3 + b12 * DIF  ) );
		mDIF = co_DIF * ( DIF_o	 + hh2 * (b33 * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO ) );
		mDIM1= co_DIM1* ( DIM1_o + hh2 * (a44 * DIF + b55 * DIM2  ) );
		mDIM2= co_DIM2* ( DIM2_o + hh2 * (a55 * DIM1 ) );

		//Neutral Drug Bound States
		mD_O  = co_D_O  * ( D_O_o	+ hh2 * (k_on * O + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS   ) );
		mD_C1 = co_D_C1 * ( D_C1_o	+ hh2 * (kc_on * C1 + a12 * D_C2 + a_33 * D_IF + b13n * D_O   ) );
		mD_C2 = co_D_C2 * ( D_C2_o	+ hh2 * (kc_on * C2 + a11 * D_C3 + a_33 * D_IC2 + b12 * D_C1  ) );
		mD_C3 = co_D_C3 * ( D_C3_o	+ hh2 * (kc_on * C3 + a_33 * D_IC3 + b11 * D_C2   ) );
		mD_OS = co_D_OS * ( D_OS_o	+ hh2 * (ax2 * D_O + ki_on * OS ) );
		mD_IC3= co_D_IC3* ( D_IC3_o	+ hh2 * (b_33 * D_C3 + b11 * D_IC2 + ki_on * IC3 ) );
		mD_IC2= co_D_IC2* ( D_IC2_o	+ hh2 * (b_33 * D_C2 + a11 * D_IC3 + b12 * D_IF + ki_on * IC2  ) );
		mD_IF = co_D_IF * ( D_IF_o	+ hh2 * (b_33 * D_C1 + a12 * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF  ) );
		mD_IM1= co_D_IM1* ( D_IM1_o	+ hh2 * (a_44 * D_IF + b_55 * D_IM2  ) );
		mD_IM2= co_D_IM2* ( D_IM2_o	+ hh2 * (a_55 * D_IM1  ) );

		err_sum = fabs( IC3 - mIC3 ) + fabs( IC2 - mIC2 ) +  + fabs( IF - mIF ) +  + fabs( IM1 - mIM1 ) +  + fabs( IM2 - mIM2 ) +  + fabs( C3 - mC3 ) +  + fabs( C2 - mC2 ) +  + fabs( C1 - mC1 ) +  + fabs( O - mO ) +  + fabs( OS - mOS ) +  + fabs( DC3 - mDC3 ) +  + fabs( DC2 - mDC2 ) +  + fabs( DC1 - mDC1 ) +  + fabs( DO - mDO ) +  + fabs( DOS - mDOS ) +   + fabs( DIC3 - mDIC3 ) +  + fabs( DIC2 - mDIC2 ) +  + fabs( DIF - mDIF ) +  + fabs( DIM1 - mDIM1 ) +  + fabs( DIM2 - mDIM2 ) +  + fabs( D_C3 - mD_C3 ) +  + fabs( D_C2 - mD_C2 ) +  + fabs( D_C1 - mD_C1 ) +  + fabs( D_O - mD_O ) +  + fabs( D_OS - mD_OS ) +  + fabs( D_IC3 - mD_IC3 ) +  + fabs( D_IC2 - mD_IC2 ) +  + fabs( D_IF - mD_IF ) +  + fabs( D_IM1 - mD_IM1 ) +  + fabs( D_IM2 - mD_IM2 );

		if (fabs(err_sum - err_sum_blas) > 1e-14 ){
			std::cout << "Old err_sum: " << err_sum << " Blas err_sum: " << err_sum_blas << std::endl;
			throw std::runtime_error("Error in Na Channel BLAS computations");
		}

		IC3   = mIC3;
		IC2   = mIC2;
		IF    = mIF;
		IM1   = mIM1;
		IM2   = mIM2;
		C3    = mC3;
		C2    = mC2;
		C1    = mC1;
		O     = mO;
		OS    = mOS;
		DC3   = mDC3;
		DC2   = mDC2;
		DC1   = mDC1;
		DO    = mDO;
		DOS   = mDOS;
		DIC3  = mDIC3;
		DIC2  = mDIC2;
		DIF   = mDIF;
		DIM1  = mDIM1;
		DIM2  = mDIM2;
		D_C3  = mD_C3;
		D_C2  = mD_C2;
		D_C1  = mD_C1;
		D_O   = mD_O;
		D_OS  = mD_OS;
		D_IC3 = mD_IC3;
		D_IC2 = mD_IC2;
		D_IF  = mD_IF;
		D_IM1 = mD_IM1;
		D_IM2 = mD_IM2;

		iter++;
	}
	MarkovIterations = iter;
	MarkovIterationError = err_sum;

	// cout << iter << "\t" << err_sum << endl << endl;

	/**************************************************************************/
	//		//        MF = eye(30) - 0.5 * dt * MB;
	//		//        MB = 0.5 * dt * MB + eye(30);
	//		for (id1 = 0; id1 < 30; id1 += 1){
	//			for (id2 = 0; id2 < 30; id2 += 1) {
	//				MB[id1][id2] = 0.5 * dt * MB[id1][id2];
	//				MF[id1][id2] = -MB[id1][id2];
	//			}
	//			MB[id1][id1] = MB[id1][id1] + 1;
	//			MF[id1][id1] = MF[id1][id1] + 1;
	//		}



	//		// Compute inv(MB) : Inverse of the matrix MB
	//		// Get LU factorization first,
	//		// then use LU to compute inv(MB)
	//		id3 = 0;
	//		for (id1 = 0; id1 < 30; id1 += 1) {
	//			for (id2 = 0; id2 < 30; id2 += 1){
	//				// For Mac using Accelerate.h
	//				// Fortran is Column major, C++ is Row major
	//				//  A[id3] = MB[id2][id1];
	//
	//				// For MTJ C++ subroutine
	//				A[id3] = MB[id1][id2];
	//
	//				id3 += 1;
	//			}
	//		}
	//
	//		LWork = N;
	//		// For Mac using Accelerate.h
	//		// dgetrf_( &M, &N, A, &LDA, IPIV, &info);
	//		// dgetri_( &N, A, &LDA, IPIV, Work, &LWork, &info);
	//
	//
	//		// dgetrf( &M, &N, A, &LDA, IPIV, &info);
	//		// dgetri( &N, A, &LDA, IPIV, Work, &LWork, &info);
	//
	//		// For MTJ C++ subroutine
	//		mat_inv( M, A, invA );
	//
	//
	//		// if (info != 0) { return info;}
	//		// Let invMB = inv(MB)
	//		id3 = 0;
	//		for (id1 = 0; id1 < 30; id1 += 1) {
	//			for (id2 = 0; id2 < 30; id2 += 1){
	//				// Fortran is Column major, C++ is Row major
	//				//  invMB[id2][id1] = A[id3];
	//
	//				// For MTJ C++ subroutine
	//				 invMB[id1][id2] = invA[id3];
	//
	//				id3 += 1;
	//				// cout << invMB[id1][id2] << "\t";
	//			}
	//			// cout << endl;
	//		}
	//
	//
	//
	//		for (id1 = 0; id1 < 30; id1 += 1) {
	//			for (id2 = 0; id2 < 30; id2 += 1) {
	//				MT[id1][id2] = 0;
	//				for(id3 = 0; id3 < 30; id3 += 1) {
	//					MT[id1][id2] += invMB[id1][id3] * MF[id3][id2];
	//				}
	//			}
	//		}
	//
	//
	//
	//
	//		// Compute Vect_(n+1) = inv(MB) * Vect_n
	//		for (id1 = 0; id1 < 30; id1 += 1){
	//			nextVect[id1] = 0;
	//			for(id2 = 0; id2 < 30; id2 += 1) {
	//				// the value for the entry id2 of the vector
	//				nextVect[id1] += MT[id1][id2] * Vect[id2];
	//			} // end of id2 looping
	//		} // end of id1 looping
	//
	//		//#pragma omp critical
	//		//{
	//		//		for ( id1 = 0; id1 < 30; id1++ ) {
	//		//		cout << Vect[id1] << ", ";
	//		//		}
	//		//	cout << endl << endl;
	//		//}
	//
	//
	//		// Vect = [IC3; IC2; IF; IM1; IM2; C3; C2; C1; O; OS; DC3; DC2; DC1; DO; DOS; ...
	//		//     DIC3; DIC2; DIF; DIM1; DIM2; D_C3; D_C2; D_C1; D_O; D_OS; ...
	//		//     D_IC3; D_IC2; D_IF; D_IM1; D_IM2];
	//
	//		//        nextVect = inv(MB) * MF * Vect;


	//		IC3 = nextVect[0];
	//		IC2 = nextVect[1];
	//		IF = nextVect[2];
	//		IM1 = nextVect[3];
	//		IM2 = nextVect[4];
	//		C3 = nextVect[5];
	//		C2 = nextVect[6];
	//		C1 = nextVect[7];
	//		O = nextVect[8];
	//		OS = nextVect[9];
	//		DC3 = nextVect[10];
	//		DC2 = nextVect[11];
	//		DC1 = nextVect[12];
	//		DO = nextVect[13];
	//		DOS = nextVect[14];
	//		DIC3 = nextVect[15];
	//		DIC2 = nextVect[16];
	//		DIF = nextVect[17];
	//		DIM1 = nextVect[18];
	//		DIM2 = nextVect[19];
	//		D_C3 = nextVect[20];
	//		D_C2 = nextVect[21];
	//		D_C1 = nextVect[22];
	//		D_O = nextVect[23];
	//		D_OS = nextVect[24];
	//		D_IC3 = nextVect[25];
	//		D_IC2 = nextVect[26];
	//		D_IF = nextVect[27];
	//		D_IM1 = nextVect[28];
	//		D_IM2 = nextVect[29];


	IC3   = mIC3;
	IC2   = mIC2;
	IF    = mIF;
	IM1   = mIM1;
	IM2   = mIM2;
	C3    = mC3;
	C2    = mC2;
	C1    = mC1;
	O     = mO;
	OS    = mOS;
	DC3   = mDC3;
	DC2   = mDC2;
	DC1   = mDC1;
	DO    = mDO;
	DOS   = mDOS;
	DIC3  = mDIC3;
	DIC2  = mDIC2;
	DIF   = mDIF;
	DIM1  = mDIM1;
	DIM2  = mDIM2;
	D_C3  = mD_C3;
	D_C2  = mD_C2;
	D_C1  = mD_C1;
	D_O   = mD_O;
	D_OS  = mD_OS;
	D_IC3 = mD_IC3;
	D_IC2 = mD_IC2;
	D_IF  = mD_IF;
	D_IM1 = mD_IM1;
	D_IM2 = mD_IM2;

	// Compute the difference in states computed using BLAS vs. old
	difference = 0;
	for (unsigned int i=0; i<nNaChStates; i++) {
		double statedifference = fabs( state[i] - gsl_vector_get( Ynp1, i ) );
		if (statedifference > 1e-14 ){
			std::cout << "Diff between BLAS and old code: " << statedifference << " for state " << i << std::endl;
			throw std::runtime_error("Diff between BLAS computations vs old code");
		}
	}
	// Once the BLAS code is working, use the following loop
	for (unsigned int i=0; i<nNaChStates; i++)
		state[i] = gsl_vector_get( Ynp1, i );


	checkNaChannelStates();

	I_Na = G_Na*(O)*(V - E_Na);
}


void Moreno2011Cell::semi_impl_WT_SCN5A_Flecainide_old (const double dt){
	double O_n , OS_n , C1_n , C2_n , C3_n , IC3_n , IC2_n , IF_n , IM1_n , IM2_n ;
	double DO_n , DOS_n , DC1_n , DC2_n , DC3_n , DIC3_n , DIC2_n , DIF_n , DIM1_n , DIM2_n ;
	double D_O_n , D_OS_n , D_C1_n , D_C2_n , D_C3_n , D_IC3_n , D_IC2_n , D_IF_n , D_IM1_n , D_IM2_n ;

	double O_o , OS_o , C1_o , C2_o , C3_o , IC3_o , IC2_o , IF_o , IM1_o , IM2_o ;
	double DO_o , DOS_o , DC1_o , DC2_o , DC3_o , DIC3_o , DIC2_o , DIF_o , DIM1_o , DIM2_o ;
	double D_O_o , D_OS_o , D_C1_o , D_C2_o , D_C3_o , D_IC3_o , D_IC2_o , D_IF_o , D_IM1_o , D_IM2_o ;


//	double MB[30][30], MF[30][30], MT[30][30], invMB[30][30], nextVect[30], Vect[30];
//	int id1, id2, id3;

	// For Mac using Accelerate.h
	//	__CLPK_doublereal x,y;
	//	__CLPK_integer info, LWork;
	//	__CLPK_integer M = 30;
	//	__CLPK_integer N = 30;
	//	__CLPK_integer LDA = 30;
	//	__CLPK_integer IPIV[30];
	//	__CLPK_doublereal A[900], Work[30];

//	double x,y;
//	int info, LWork;
//	int M = 30;
//	int N = 30;
//	int LDA = 30;
//	int IPIV[30];
//	double A[900], Work[30], invA[900];


	double a11, a12, a13, a2, a3, a4, a5;
	double b11, b12, b13, b2, b3, b4, b5;
	double ax, bx, ax1, bx1, ax2, bx2;
	double a13c, b13c, a13n, b13n;
	double a22, b22, a33, b33, a44, b44, a55, b55;
	double a_22, b_22, a_33, b_33, a_44, b_44, a_55, b_55;

	double kon, k_on, koff, k_off, kcon, kc_on, kcoff, kc_off, ki_on, ki_off;
	double kd_open;
	double hh;

	hh = dt;
	const double G_Na=15.;	//27.5  //18.5

	const double Q10=3;

	const double Tfactor=1.0/(pow(Q10, (37.0-(T-273))/10.0));

	const double pH=7.4;
	const double pKa=9.3;
	const double portion = 1/(1+ pow(10, (pH-pKa)) );
	const double diffusion=5500;

	const double drug_charged=cellenv->drug*portion;
	const double drug_neutral=cellenv->drug*(1-portion);
	const double drug_distance= -0.7;

	const double E_Na = (R*T/F)*log(Na_out/Na_in);


	//Rate Constants **********************************************************

	//WT Fits Reduced Model (no IM1, IM2)
	a11= Tfactor*8.5539/(7.4392e-2*exp(-V/17.0)+ 2.0373e-1*exp(-V/150));
	a12= Tfactor*8.5539/(7.4392e-2*exp(-V/15.0)+ 2.0373e-1*exp(-V/150));
	a13= Tfactor*8.5539/(7.4392e-2*exp(-V/12.0)+ 2.0373e-1*exp(-V/150));
	b11= Tfactor*7.5215e-2*exp(-V/20.3);
	b12= Tfactor*2.7574*exp(-(V-5)/20.3);
	b13= Tfactor*4.7755e-1*exp(-(V-10)/20.3);

	a3 = Tfactor*5.1458e-6*exp(-V/8.2471);
	b3=Tfactor*6.1205*exp((V)/13.542);


	a2= Tfactor*(13.370*exp(V/43.749));
	b2= ((a13*a2*a3)/(b13*b3));

	a4 = 0*a2;
	b4 = 0*a3;
	a5= 0*a2;
	b5 = 0*a3;

	ax = 3.4229e-2*a2;
	bx = 1.7898e-2*a3;


	//************************************ Global Fitting *********************
	ax1 =5.7839e-05 * ax;
	bx1 =  1.6689e-08* bx;
	a13c = 3.6324e-03 *a13;
	a22 = 1.4847e+03 *a2;
	b33 =  1.7352e-06* b3;
	a33 = 6.7505e-05 * a3;
	a44 =  2.4135e+00* a2;
	b44 =  4.9001e-02* a3;
	a55 = b55 = 0;



	ax2 = 2.6126e-01 * ax;
	a13n = 2.6452e+00 * a13;
	a_22 =  4.2385e+01 * a2;
	b_33 = 2.1181e+00 * b3;
	a_44 =  1.0326e-03 * a2;
	b_44 = 2.1378e-02 * a3;

	a_55 = 0;
	b_55 = 0;

	const double kd0=11.2*(1e-6);
	kd_open=kd0*exp( (drug_distance*V*F) /(R*T));


	// charged drug
	kon=drug_charged*diffusion;
	koff=kd_open*diffusion;
	kcoff = koff;
	kcon = kon;

	if (cellenv->drug ==0 || drug_charged ==0 ){b13c = 0;}
	else{b13c = (b13*kcon*koff*a13c)/(kon*kcoff*a13);}

	if (b13c ==0){b22 = 0;}
	else {b22=(a13c*a22*a33)/(b13c*b33);}

	// neutral drug
	k_on = drug_neutral*diffusion;
	k_off=400*(1e-6)*diffusion;
	ki_on=k_on/2;
	ki_off=5.4*(1e-6)*diffusion;
	kc_on=k_on/2;
	kc_off=800*(1e-6)*diffusion;

	if (cellenv->drug ==0 || drug_neutral ==0 ){a_33 = 0;}
	else {a_33 = (ki_off*a3*kc_on*b_33)/(ki_on*kc_off*b3);}

	if (cellenv->drug ==0 || drug_neutral ==0){b13n = 0;}
	else {b13n = (b13*kc_on*a13n*k_off)/(kc_off*a13*k_on);}

	if (b13n==0){b_22 =0;}
	else {b_22 = (a_33*a13n*a_22)/(b_33*b13n);}

	if (cellenv->drug ==0 || drug_neutral ==0){bx2 = 0;}
	else{bx2 = (bx*k_on*ax2*ki_off)/(ax*ki_on*k_off);}


		// Backward method
		// MB * nextVect_(n+1) = Vect_n
		// where MB is a matrix, and Vect_i is a vector

		// Set up the matrix
		//  MB = zeros(30);
//		for (id1 = 0; id1 < 30; id1 += 1){
//			for(id2 = 0; id2 < 30; id2 += 1) {
//				MB[id1][id2] = 0;
//			} // end of id2 looping
//		} // end of id1 looping


		const double coef_O = (b13 + a2  + kon + k_on + ax);
		const double coef_C1 = (b12 + b3 + a13  + kcon + kc_on);
		const double coef_C2 = (b11 + b3 + a12  + kcon + kc_on);
		const double coef_C3 = (b3 + a11  + kcon + kc_on);
		const double coef_IC3 = (a11 + a3 + ki_on);
		const double coef_IC2 = (b11 + a3 + a12 + ki_on);
		const double coef_IF = (b12 + b2 + a3 + a4 + ki_on);
		const double coef_IM1 = (b4 + a5);
		const double coef_IM2 = b5;
		const double coef_OS = (bx + ki_on);

		const double coef_DO = (koff + b13c + a22 + ax1 );
		const double coef_DC1 = (kcoff + b12 + b33 + a13c );
		const double coef_DC2 = (kcoff + b11 + b33 + a12 );
		const double coef_DC3 = (kcoff+ b33 + a11 );
		const double coef_DOS = bx1;
		const double coef_DIC3 = (a11 + a33);
		const double coef_DIC2 = (a33 + b11 + a12);
		const double coef_DIF = (a33 + b12 + a44 + b22);
		const double coef_DIM1 = ( b44 + a55 );
		const double coef_DIM2 = b55 ;

		const double coef_D_O = (k_off + b13n + a_22 + ax2 );
		const double coef_D_C1 = (kc_off + b12 + b_33 + a13n );
		const double coef_D_C2 = (kc_off + b11 + b_33 + a12 );
		const double coef_D_C3 = (kc_off + b_33 + a11 );
		const double coef_D_OS = (bx2 + ki_off);
		const double coef_D_IC3 = (a_33 + a11 + ki_off);
		const double coef_D_IC2 = (a_33 + b11 + a12 + ki_off);
		const double coef_D_IF = (a_33 + a_44 + b_22 + b12 + ki_off);
		const double coef_D_IM1 = (b_44 + a_55);
		const double coef_D_IM2 = b_55;

		double hh2 = hh/2;

		const double co_O = 1. / ( 1 + hh2 * ( b13 + a2  + kon + k_on + ax ) );
		const double co_C1 = 1. / ( 1 + hh2 * ( b12 + b3 + a13  + kcon + kc_on));
		const double co_C2 = 1. / ( 1 + hh2 * ( b11 + b3 + a12  + kcon + kc_on));
		const double co_C3 = 1. / ( 1 + hh2 * ( b3 + a11  + kcon + kc_on));
		const double co_IC3 = 1. / ( 1 + hh2 * ( a11 + a3 + ki_on));
		const double co_IC2 = 1. / ( 1 + hh2 * ( b11 + a3 + a12 + ki_on));
		const double co_IF = 1. / ( 1 + hh2 * ( b12 + b2 + a3 + a4 + ki_on));
		const double co_IM1 = 1. / ( 1 + hh2 * ( b4 + a5));
		const double co_IM2 = 1. / ( 1 + hh2 * ( b5 ));
		const double co_OS = 1. / ( 1 + hh2 * ( bx + ki_on));

		const double co_DO = 1. / ( 1 + hh2 * ( koff + b13c + a22 + ax1 ));
		const double co_DC1 = 1. / ( 1 + hh2 * ( kcoff + b12 + b33 + a13c ));
		const double co_DC2 = 1. / ( 1 + hh2 * ( kcoff + b11 + b33 + a12 ));
		const double co_DC3 = 1. / ( 1 + hh2 * ( kcoff+ b33 + a11 ));
		const double co_DOS = 1. / ( 1 + hh2 * ( bx1));
		const double co_DIC3 = 1. / ( 1 + hh2 * ( a11 + a33));
		const double co_DIC2 = 1. / ( 1 + hh2 * ( a33 + b11 + a12));
		const double co_DIF = 1. / ( 1 + hh2 * ( a33 + b12 + a44 + b22));
		const double co_DIM1 = 1. / ( 1 + hh2 * (  b44 + a55 ));
		const double co_DIM2 = 1. / ( 1 + hh2 * ( b55 ) );

		const double co_D_O = 1. / ( 1 + hh2 * ( k_off + b13n + a_22 + ax2 ));
		const double co_D_C1 = 1. / ( 1 + hh2 * ( kc_off + b12 + b_33 + a13n ));
		const double co_D_C2 = 1. / ( 1 + hh2 * ( kc_off + b11 + b_33 + a12 ));
		const double co_D_C3 = 1. / ( 1 + hh2 * ( kc_off + b_33 + a11 ));
		const double co_D_OS = 1. / ( 1 + hh2 * ( bx2 + ki_off));
		const double co_D_IC3 = 1. / ( 1 + hh2 * ( a_33 + a11 + ki_off));
		const double co_D_IC2 = 1. / ( 1 + hh2 * ( a_33 + b11 + a12 + ki_off));
		const double co_D_IF = 1. / ( 1 + hh2 * ( a_33 + a_44 + b_22 + b12 + ki_off));
		const double co_D_IM1 = 1. / ( 1 + hh2 * ( b_44 + a_55));
		const double co_D_IM2 = 1. / ( 1 + hh2 * ( b_55 ));


		//Drug Free States
		O_o = O + hh2 * ( a13 * C1 + b2 * IF + koff * DO + k_off * D_O + bx * OS  - O * coef_O );
		C1_o = C1 + hh2 * (a12 * C2 + a3 * IF + b13 * O  + kcoff * DC1 + kc_off * D_C1 - C1 * coef_C1 );
		C2_o = C2 + hh2 * (a11 * C3 + a3 * IC2 + b12 * C1  + kcoff * DC2 + kc_off * D_C2 - C2 * coef_C2 );
		C3_o = C3 + hh2 * (a3 * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 - C3 * coef_C3 );
		IC3_o = IC3 + hh2 * (b3 * C3 + b11 * IC2 + ki_off * D_IC3 - IC3 * coef_IC3 );
		IC2_o = IC2 + hh2 * (a11 * IC3 + b3 * C2 + b12 * IF + ki_off * D_IC2 - IC2 * coef_IC2 );
		IF_o = IF + hh2 * (a12 * IC2 + b3 * C1 + b4 * IM1 + a2 * O + ki_off * D_IF - IF * coef_IF );
		IM1_o = IM1 + hh2 * (a4 * IF + b5 * IM2 - IM1 * coef_IM1 );
		IM2_o = IM2 + hh2 * (a5 * IM1 - IM2 * coef_IM2 );
		OS_o = OS + hh2 * (ax * O + ki_off * D_OS - OS * coef_OS );

		//Charged Drug Bound States
		DO_o = DO + hh2 * (kon * O + a13c * DC1 + bx1 * DOS + b22 * DIF  - DO * coef_DO );
		DC1_o = DC1 + hh2 * (kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO  - DC1 * coef_DC1 );
		DC2_o = DC2 + hh2 * (kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1  - DC2 * coef_DC2 );
		DC3_o = DC3 + hh2 * (kcon * C3 + b11 * DC2 + a33 * DIC3  - DC3 * coef_DC3 );
		DOS_o = DOS + hh2 * (ax1 * DO -  DOS * coef_DOS );
		DIC3_o = DIC3 + hh2 * (b33 * DC3 + b11 * DIC2 - DIC3 * coef_DIC3 );
		DIC2_o = DIC2 + hh2 * (b33 * DC2 + a11 * DIC3 + b12 * DIF - DIC2 * coef_DIC2 );
		DIF_o = DIF + hh2 * (b33 * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO - DIF * coef_DIF );
		DIM1_o = DIM1 + hh2 * (a44 * DIF + b55 * DIM2 - DIM1 * coef_DIM1 );
		DIM2_o = DIM2 + hh2 * (a55 * DIM1 - DIM2 * coef_DIM2 );

		//Neutral Drug Bound States
		D_O_o = D_O + hh2 * (k_on * O + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS  - D_O * coef_D_O );
		D_C1_o = D_C1 + hh2 * (kc_on * C1 + a12 * D_C2 + a_33 * D_IF + b13n * D_O  - D_C1 * coef_D_C1 );
		D_C2_o = D_C2 + hh2 * (kc_on * C2 + a11 * D_C3 + a_33 * D_IC2 + b12 * D_C1  - D_C2 * coef_D_C2 );
		D_C3_o = D_C3 + hh2 * (kc_on * C3 + a_33 * D_IC3 + b11 * D_C2  - D_C3 * coef_D_C3 );
		D_OS_o = D_OS + hh2 * (ax2 * D_O + ki_on * OS - D_OS * coef_D_OS );
		D_IC3_o = D_IC3 + hh2 * (b_33 * D_C3 + b11 * D_IC2 + ki_on * IC3 - D_IC3 * coef_D_IC3 );
		D_IC2_o = D_IC2 + hh2 * (b_33 * D_C2 + a11 * D_IC3 + b12 * D_IF + ki_on * IC2 - D_IC2 * coef_D_IC2 );
		D_IF_o = D_IF + hh2 * (b_33 * D_C1 + a12 * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF - D_IF * coef_D_IF );
		D_IM1_o = D_IM1 + hh2 * (a_44 * D_IF + b_55 * D_IM2 - D_IM1 * coef_D_IM1 );
		D_IM2_o = D_IM2 + hh2 * (a_55 * D_IM1 - D_IM2 * coef_D_IM2 );

		int iter = 0;
		double err_sum = 1;
		while ( err_sum > 1E-100 && iter < 100 ) {
			//Drug Free States
			O_n = co_O * ( O_o + hh2 * ( a13 * C1 + b2 * IF + koff * DO + k_off * D_O + bx * OS ) );
			C1_n = co_C1 * ( C1_o + hh2 * (a12 * C2 + a3 * IF + b13 * O  + kcoff * DC1 + kc_off * D_C1 ) );
			C2_n = co_C2 * ( C2_o + hh2 * (a11 * C3 + a3 * IC2 + b12 * C1  + kcoff * DC2 + kc_off * D_C2 ) );
			C3_n = co_C3 * ( C3_o + hh2 * (a3 * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 ) );
			IC3_n = co_IC3 * ( IC3_o + hh2 * (b3 * C3 + b11 * IC2 + ki_off * D_IC3  ) );
			IC2_n = co_IC2 * ( IC2_o + hh2 * (a11 * IC3 + b3 * C2 + b12 * IF + ki_off * D_IC2  ) );
			IF_n = co_IF * ( IF_o + hh2 * (a12 * IC2 + b3 * C1 + b4 * IM1 + a2 * O + ki_off * D_IF  ) );
			IM1_n = co_IM1 * ( IM1_o + hh2 * (a4 * IF + b5 * IM2  ) );
			IM2_n = co_IM2 * ( IM2_o + hh2 * (a5 * IM1  ) );
			OS_n = co_OS * ( OS_o + hh2 * (ax * O + ki_off * D_OS  ) );

			//Charged Drug Bound States
			DO_n = co_DO * ( DO_o + hh2 * (kon * O + a13c * DC1 + bx1 * DOS + b22 * DIF   ) );
			DC1_n = co_DC1 * ( DC1_o + hh2 * (kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO  ) );
			DC2_n = co_DC2 * ( DC2_o + hh2 * (kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1   ) );
			DC3_n = co_DC3 * ( DC3_o + hh2 * (kcon * C3 + b11 * DC2 + a33 * DIC3   ) );
			DOS_n = co_DOS * ( DOS_o + hh2 * (ax1 * DO  ) );
			DIC3_n = co_DIC3 * ( DIC3_o + hh2 * (b33 * DC3 + b11 * DIC2 ) );
			DIC2_n = co_DIC2 * ( DIC2_o + hh2 * (b33 * DC2 + a11 * DIC3 + b12 * DIF  ) );
			DIF_n = co_DIF * ( DIF_o + hh2 * (b33 * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO ) );
			DIM1_n = co_DIM1 * ( DIM1_o + hh2 * (a44 * DIF + b55 * DIM2  ) );
			DIM2_n = co_DIM2 * ( DIM2_o + hh2 * (a55 * DIM1 ) );

			//Neutral Drug Bound States
			D_O_n = co_D_O * ( D_O_o + hh2 * (k_on * O + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS   ) );
			D_C1_n = co_D_C1 * ( D_C1_o + hh2 * (kc_on * C1 + a12 * D_C2 + a_33 * D_IF + b13n * D_O   ) );
			D_C2_n = co_D_C2 * ( D_C2_o + hh2 * (kc_on * C2 + a11 * D_C3 + a_33 * D_IC2 + b12 * D_C1  ) );
			D_C3_n = co_D_C3 * ( D_C3_o + hh2 * (kc_on * C3 + a_33 * D_IC3 + b11 * D_C2   ) );
			D_OS_n = co_D_OS * ( D_OS_o + hh2 * (ax2 * D_O + ki_on * OS ) );
			D_IC3_n = co_D_IC3 * ( D_IC3_o + hh2 * (b_33 * D_C3 + b11 * D_IC2 + ki_on * IC3 ) );
			D_IC2_n = co_D_IC2 * ( D_IC2_o + hh2 * (b_33 * D_C2 + a11 * D_IC3 + b12 * D_IF + ki_on * IC2  ) );
			D_IF_n = co_D_IF * ( D_IF_o + hh2 * (b_33 * D_C1 + a12 * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF  ) );
			D_IM1_n = co_D_IM1 * ( D_IM1_o + hh2 * (a_44 * D_IF + b_55 * D_IM2  ) );
			D_IM2_n = co_D_IM2 * ( D_IM2_o + hh2 * (a_55 * D_IM1  ) );


			err_sum = fabs( IC3 - IC3_n ) + fabs( IC2 - IC2_n ) +  + fabs( IF - IF_n ) +  + fabs( IM1 - IM1_n ) +  + fabs( IM2 - IM2_n ) +  + fabs( C3 - C3_n ) +  + fabs( C2 - C2_n ) +  + fabs( C1 - C1_n ) +  + fabs( O - O_n ) +  + fabs( OS - OS_n ) +  + fabs( DC3 - DC3_n ) +  + fabs( DC2 - DC2_n ) +  + fabs( DC1 - DC1_n ) +  + fabs( DO - DO_n ) +  + fabs( DOS - DOS_n ) +   + fabs( DIC3 - DIC3_n ) +  + fabs( DIC2 - DIC2_n ) +  + fabs( DIF - DIF_n ) +  + fabs( DIM1 - DIM1_n ) +  + fabs( DIM2 - DIM2_n ) +  + fabs( D_C3 - D_C3_n ) +  + fabs( D_C2 - D_C2_n ) +  + fabs( D_C1 - D_C1_n ) +  + fabs( D_O - D_O_n ) +  + fabs( D_OS - D_OS_n ) +  + fabs( D_IC3 - D_IC3_n ) +  + fabs( D_IC2 - D_IC2_n ) +  + fabs( D_IF - D_IF_n ) +  + fabs( D_IM1 - D_IM1_n ) +  + fabs( D_IM2 - D_IM2_n );


			IC3 = IC3_n;
			IC2 = IC2_n;
			IF = IF_n;
			IM1 = IM1_n;
			IM2 = IM2_n;
			C3 = C3_n;
			C2 = C2_n;
			C1 = C1_n;
			O = O_n;
			OS = OS_n;
			DC3 = DC3_n;
			DC2 = DC2_n;
			DC1 = DC1_n;
			DO = DO_n;
			DOS = DOS_n;
			DIC3 = DIC3_n;
			DIC2 = DIC2_n;
			DIF = DIF_n;
			DIM1 = DIM1_n;
			DIM2 = DIM2_n;
			D_C3 = D_C3_n;
			D_C2 = D_C2_n;
			D_C1 = D_C1_n;
			D_O = D_O_n;
			D_OS = D_OS_n;
			D_IC3 = D_IC3_n;
			D_IC2 = D_IC2_n;
			D_IF = D_IF_n;
			D_IM1 = D_IM1_n;
			D_IM2 = D_IM2_n;


			iter++;
		}
		MarkovIterations = iter;
		MarkovIterationError = err_sum;

		// cout << iter << "\t" << err_sum << endl << endl;

		/*********************************************************************/
		//		//        MF = eye(30) - 0.5 * dt * MB;
		//		//        MB = 0.5 * dt * MB + eye(30);
		//		for (id1 = 0; id1 < 30; id1 += 1){
		//			for (id2 = 0; id2 < 30; id2 += 1) {
		//				MB[id1][id2] = 0.5 * dt * MB[id1][id2];
		//				MF[id1][id2] = -MB[id1][id2];
		//			}
		//			MB[id1][id1] = MB[id1][id1] + 1;
		//			MF[id1][id1] = MF[id1][id1] + 1;
		//		}



		//		// Compute inv(MB) : Inverse of the matrix MB
		//		// Get LU factorization first,
		//		// then use LU to compute inv(MB)
		//		id3 = 0;
		//		for (id1 = 0; id1 < 30; id1 += 1) {
		//			for (id2 = 0; id2 < 30; id2 += 1){
		//				// For Mac using Accelerate.h
		//				// Fortran is Column major, C++ is Row major
		//				//  A[id3] = MB[id2][id1];
		//
		//				// For MTJ C++ subroutine
		//				A[id3] = MB[id1][id2];
		//
		//				id3 += 1;
		//			}
		//		}
		//
		//		LWork = N;
		//		// For Mac using Accelerate.h
		//		// dgetrf_( &M, &N, A, &LDA, IPIV, &info);
		//		// dgetri_( &N, A, &LDA, IPIV, Work, &LWork, &info);
		//
		//
		//		// dgetrf( &M, &N, A, &LDA, IPIV, &info);
		//		// dgetri( &N, A, &LDA, IPIV, Work, &LWork, &info);
		//
		//		// For MTJ C++ subroutine
		//		mat_inv( M, A, invA );
		//
		//
		//		// if (info != 0) { return info;}
		//		// Let invMB = inv(MB)
		//		id3 = 0;
		//		for (id1 = 0; id1 < 30; id1 += 1) {
		//			for (id2 = 0; id2 < 30; id2 += 1){
		//				// Fortran is Column major, C++ is Row major
		//				//  invMB[id2][id1] = A[id3];
		//
		//				// For MTJ C++ subroutine
		//				 invMB[id1][id2] = invA[id3];
		//
		//				id3 += 1;
		//				// cout << invMB[id1][id2] << "\t";
		//			}
		//			// cout << endl;
		//		}
		//
		//
		//
		//		for (id1 = 0; id1 < 30; id1 += 1) {
		//			for (id2 = 0; id2 < 30; id2 += 1) {
		//				MT[id1][id2] = 0;
		//				for(id3 = 0; id3 < 30; id3 += 1) {
		//					MT[id1][id2] += invMB[id1][id3] * MF[id3][id2];
		//				}
		//			}
		//		}
		//
		//
		//
		//
		//		// Compute Vect_(n+1) = inv(MB) * Vect_n
		//		for (id1 = 0; id1 < 30; id1 += 1){
		//			nextVect[id1] = 0;
		//			for(id2 = 0; id2 < 30; id2 += 1) {
		//				// the value for the entry id2 of the vector
		//				nextVect[id1] += MT[id1][id2] * Vect[id2];
		//			} // end of id2 looping
		//		} // end of id1 looping
		//
		//		//#pragma omp critical
		//		//{
		//		//		for ( id1 = 0; id1 < 30; id1++ ) {
		//		//		cout << Vect[id1] << ", ";
		//		//		}
		//		//	cout << endl << endl;
		//		//}
		//
		//
		//		// Vect = [IC3; IC2; IF; IM1; IM2; C3; C2; C1; O; OS; DC3; DC2; DC1; DO; DOS; ...
		//		//     DIC3; DIC2; DIF; DIM1; DIM2; D_C3; D_C2; D_C1; D_O; D_OS; ...
		//		//     D_IC3; D_IC2; D_IF; D_IM1; D_IM2];
		//
		//		//        nextVect = inv(MB) * MF * Vect;


		//		IC3 = nextVect[0];
		//		IC2 = nextVect[1];
		//		IF = nextVect[2];
		//		IM1 = nextVect[3];
		//		IM2 = nextVect[4];
		//		C3 = nextVect[5];
		//		C2 = nextVect[6];
		//		C1 = nextVect[7];
		//		O = nextVect[8];
		//		OS = nextVect[9];
		//		DC3 = nextVect[10];
		//		DC2 = nextVect[11];
		//		DC1 = nextVect[12];
		//		DO = nextVect[13];
		//		DOS = nextVect[14];
		//		DIC3 = nextVect[15];
		//		DIC2 = nextVect[16];
		//		DIF = nextVect[17];
		//		DIM1 = nextVect[18];
		//		DIM2 = nextVect[19];
		//		D_C3 = nextVect[20];
		//		D_C2 = nextVect[21];
		//		D_C1 = nextVect[22];
		//		D_O = nextVect[23];
		//		D_OS = nextVect[24];
		//		D_IC3 = nextVect[25];
		//		D_IC2 = nextVect[26];
		//		D_IF = nextVect[27];
		//		D_IM1 = nextVect[28];
		//		D_IM2 = nextVect[29];

		IC3 = IC3_n;
		IC2 = IC2_n;
		IF = IF_n;
		IM1 = IM1_n;
		IM2 = IM2_n;
		C3 = C3_n;
		C2 = C2_n;
		C1 = C1_n;
		O = O_n;
		OS = OS_n;
		DC3 = DC3_n;
		DC2 = DC2_n;
		DC1 = DC1_n;
		DO = DO_n;
		DOS = DOS_n;
		DIC3 = DIC3_n;
		DIC2 = DIC2_n;
		DIF = DIF_n;
		DIM1 = DIM1_n;
		DIM2 = DIM2_n;
		D_C3 = D_C3_n;
		D_C2 = D_C2_n;
		D_C1 = D_C1_n;
		D_O = D_O_n;
		D_OS = D_OS_n;
		D_IC3 = D_IC3_n;
		D_IC2 = D_IC2_n;
		D_IF = D_IF_n;
		D_IM1 = D_IM1_n;
		D_IM2 = D_IM2_n;

		checkNaChannelStates();
		I_Na = G_Na*(O)*(V - E_Na);
}

/****************************************************************************
 * Moreno2011Cell.h - defines cell_param struct - i.e. Cell
 *
 * Based on code from Moreno et al. 2011
 *
 *  Created on: Sep 30, 2014
 *      Author: Anthony Varghese
 ***************************************************************************/

#ifndef CELL_H_
#define CELL_H_


#define M_PI 3.14159265358979323846 /* pi - needed if using -std=c++11 */
#include <cmath>
#include <gsl/gsl_blas.h>

namespace cellmodel{
	class Cells;
}

enum NaChannelMMSolver{
	NotUsed,
	Iterative,
	LUDecomp
};

class Moreno2011Cell {
public:
	Moreno2011Cell();
	virtual ~Moreno2011Cell();
	Moreno2011Cell(const Moreno2011Cell& source);
	Moreno2011Cell& operator=(const Moreno2011Cell& rhs);

	void setParams(const double* par);
	void setXpos(unsigned int p){ pos_x = p; }
	void Copy(const Moreno2011Cell& in);
	void setEnvironment( cellmodel::Cells* env );
	void setState(double* newstate) {
		for (unsigned int i = 0; i < nstates; i++)
			state[i] = newstate[i];
	}
	/*
	 * Calling method is responsible for allocating Y array
	 */
	void getState(double* Y) {
		for (unsigned int i = 0; i < nstates; i++)
			Y[i] = state[i];
	}
	/*
	 * Calling method is responsible for allocating dYdt array
	 */
	void getRates(double* dYdt) {
		for (unsigned int i = 0; i < nstates; i++)
			dYdt[i] = rates[i];
	}
	void Calculate_Points(double t);
	double CalculateCurrents(double t);
	double CalculateCurrents(double t, const double Vx);
	void Calculate_All(const double t);
	void SetInitialConditions();
	void Calculate_Reset(double t);
	void setI_stimulus(double t);
	void updatedVdt(const double Iaxnm1, const double Iaxnp1){
		dV_old = dVdt;
		I_axnm1 = Iaxnm1;
		I_axnp1 = Iaxnp1;
		I_axial = Iaxnp1 + Iaxnm1;
		dVdt = -I_axial - I_total; // Rate of change of Vm
	}
	double getV()	const {	return V;	}
	double getdVdt()const {	return dVdt;	}
	double getI_total()const {	return I_total;	}
	void   setV	  (const double Vnew)	{ V = Vnew; }
	void   setdVdt(const double dVdtnew){ dVdt = dVdtnew; }
	/*
	 * Forward Euler update of V
	 */
	void   fe_updateV(const double dt)	{ V = V + dVdt * dt; }
	double getI_K1()	const {	return I_K1; }
	double getI_Na()	const {	return I_Na; }
	double getI_axnm1() const {	return I_axnm1; }
	double getI_axnp1() const {	return I_axnp1; }
	double getI_axial() const {	return I_axial; }
	// Na Channel States
	double getNaClosed()const { return C1 + C2 + C3; }
	double getNaIC23()	const { return IC2 + IC3; }
	double getNaIF()	const { return IF; }
	double getNaOS()	const { return OS; }
	double getNaO()		const { return O; }
	double getNaDClosed() const { return DC1 + DC2 + DC3; }
	double getNaDIC23()	const { return DIC2 + DIC3; }
	double getNaDIF()	const { return DIF; }
	double getNaDOS()	const { return DOS; }
	double getNaDO()	const { return DO; }
	double getNaDIM1()	const { return DIM1; }
	double getNaD_Closed() const { return D_C1 + D_C2 + D_C3; }
	double getNaD_IC23()const { return D_IC2 + D_IC3; }
	double getNaD_IF()	const { return D_IF; }
	double getNaD_OS()	const { return D_OS; }
	double getNaD_O()	const { return D_O; }
	double getNaD_IM1()	const { return D_IM1; }

	/*
	 * computeRushLarsenStepState
	 * - take a time step using Rush Larsen for gating variables
	 *   and semi-implicit for Markov State models
	 * Calling method is responsible for allocating Y array
	 */
	void computeRushLarsenStepState( const double t, const double dt );
	/*
	 * computeRushLarsen2StepState
	 * - take a time step using Second-order Rush Larsen for all states
	 *   and semi-implicit for Markov State models
	 * Calling method is responsible for allocating Y array
	 */
	void computeRushLarsen2StepState( const double t, const double dt, const int step, const double Vhalf );

	// Indices into the state array
	static const unsigned int nNaChFirstState = 0;
	// Number of states for Na Channel Markov state model:
	static const unsigned int nNaChLastState = 29;
	static const unsigned int nNaChStates = nNaChLastState + 1;
	// start of gating variables:
	static const unsigned int ngatem = nNaChStates;
	static const unsigned int nLastState = ngatem + 20;		// last state
	// Currently, the cell transmembrane potential is the last state
	static const unsigned int nV = nLastState;
	static const unsigned int nstates = nV + 1;	// Total number of states

	// Cell States
	double state[nstates]; // state variables for integration

	double peak_slope, t_min, V_min, t_thr, V_thr, t_max, V_max, t_EAD, V_EAD;
	double t_EAD2, V_EAD2, t_90, V_90, dV_old;
	bool flag2;
	double APD_90, DI, L_EAD;

	// To measure the work involved in the iterations:
	int 	MarkovIterations;
	double	MarkovIterationError;

private:
	// Need a pointer back to the container class to
	// find out conditions like waitTime, drug conc, etc.
	cellmodel::Cells* cellenv = nullptr;

	/****************************************************************************
	 Universal Constants
	 *****************************************************************************/
	static constexpr double R = 8314.472;		// J/mol*K
	static constexpr double T = 310;			// K
	static constexpr double F = 96485.3415;
	static constexpr double RTF = R*T/F;

	// const double Cm = 2.0;					// uF/ cm^2 -- not used anywhere
	static constexpr double CAP = 0.185;		// Cellular capacitance
	static constexpr double rho = 162;			// ohm*cm

	// Channel conductances - these can be set externally
	// used to be:
	// static constexpr double G_Na = 15.0;	//27.5  //18.5
	double G_Na = 15.0;		// 15.0 in Moreno et al 2011
	double GK1	= 5.405;	// 5.405 in Moreno et al 2011

	// Ion Valences
	static constexpr double z_Na = 1;			// Valence of Na
	static constexpr double z_Ca = 2;			// Valence of Ca
	static constexpr double z_K = 1;			// Valence of K

	//Cell Geometry
	static constexpr double pi = M_PI;			// Should be 3.141592
	static constexpr double S_cg = 0.2;			// Surface to volume ratio (um^-1)


	//Intracellular volumes
	static constexpr double V_cyto=0.016404;	//16404 uL
	static constexpr double V_sr=0.001094;
	static constexpr double V_ss=0.00005468;

	// initial conditions of ion concentrations
	static constexpr double Na_out = 140;
	static constexpr double Ca_out = 2.0;
	static constexpr double K_out = 5.4;

	static constexpr double Vmin_default =	 0.0;	// mV - a "large" value so that minimums (expected to be < 0) can be detected
	static constexpr double Vthr_default = -90.0;	// mV - a "low" value so that the threshold can be detected

	bool initializedExternally = false;
	static const NaChannelMMSolver NaChanSolver = Iterative; // LUDecomp; // Iterative;

	// Arrays for each cell:
	double rates[nstates]; // rates of change of state - derivative
	// Total number of currents and auxilliary values:
	static const unsigned int naux = 26;
	double currents[naux];

	int Cell_type;
	unsigned int pos_x=0;
//	unsigned int pos_y=0, pos_z=0; // Position of cell in 3-d

	// Na Channel constants
	static constexpr double Q10  = 3;
	// This does not work in clang, ok in gcc
	//static constexpr double Tfactor = 1.0/(pow(Q10, (37.0-(T-273))/10.0));
	static constexpr double pH   = 7.4;


	// Na channel markov state model - the individual states are only needed for outputs
	double &IC3		= state[ 0], &IC2	= state[1];
	double &IF		= state[ 2], &IM1	= state[3];
	double &IM2		= state[ 4], &C3	= state[5];
	double &C2		= state[ 6], &O		= state[7];
	double &OS		= state[ 8], &C1	= state[9];

	double &DIC3	= state[10], &DIC2	= state[11];
	double &DIF 	= state[12], &DIM1	= state[13];
	double &DIM2	= state[14], &DC3 	= state[15];
	double &DC2		= state[16], &DO	= state[17];
	double &DOS		= state[18], &DC1	= state[19];

	double &D_IC3	= state[20], &D_IC2 = state[21];
	double &D_IF	= state[22], &D_IM1 = state[23];
	double &D_IM2	= state[24], &D_C3  = state[25];
	double &D_C2	= state[26], &D_O   = state[27];
	double &D_OS	= state[28], &D_C1  = state[29];
	double &V  = state[ngatem + 20];

	// Gating variables
	double &m	= state[ngatem],		&h	  = state[ngatem + 1];
	double &j   = state[ngatem + 2];
	double &d   = state[ngatem + 3],	&f	  = state[ngatem + 4];
	double &f2  = state[ngatem + 5],	&f_Ca = state[ngatem + 6];
	double &r   = state[ngatem + 7],	&s	  = state[ngatem + 8];
	double &xr1 = state[ngatem + 9],	&xr2  = state[ngatem + 10];
	double &xs  = state[ngatem + 11];
	double &R_bar = state[ngatem + 12];
	double &mL    = state[ngatem + 13], &hL	  = state[ngatem + 14];
	double &Ca_in = state[ngatem + 15], &Ca_sr = state[ngatem + 16];
	double &Ca_ss = state[ngatem + 17];
	double &Na_in = state[ngatem + 18], &K_in = state[ngatem + 19];
	double &dVdt  = rates[ngatem + 20];

	// old gating variable values before half-step in 2nd order Rush Larsen
	double oldgates[20];
	double &oldm	= oldgates[0],	&oldh		= oldgates[1];
	double &oldj	= oldgates[2];
	double &oldd	= oldgates[3],	&oldf		= oldgates[4];
	double &oldf2	= oldgates[5],	&oldf_Ca	= oldgates[6];
	double &oldr	= oldgates[7],	&olds		= oldgates[8];
	double &oldxr1	= oldgates[9],	&oldxr2 	= oldgates[10];
	double &oldxs	= oldgates[11];
	double &oldR_bar = oldgates[12];
	double &oldmL    = oldgates[13], &oldhL		= oldgates[14];
	double &oldCa_in = oldgates[15], &oldCa_sr	= oldgates[16];
	double &oldCa_ss = oldgates[17];
	double &oldNa_in = oldgates[18], &oldK_in	= oldgates[19];

	// Currents and auxilliary values - anything that may be useful to print out for plotting
	double &I_total = currents[0];
	double &I_axial = currents[1];
	double &I_axnm1 = currents[2], &I_axnp1 = currents[3];
	double &I_stim	= currents[4];
	double &I_Na	= currents[5], &I_Na_L = currents[6];
	double &mI_Na	= currents[7], &I_Na_b = currents[8];
	double &I_Ca_L	= currents[9], &I_Ca_b = currents[10];
	double &I_Kr	= currents[11], &I_Ks = currents[12];
	double &I_K1	= currents[13], &I_Kp = currents[14];
	double &I_to	= currents[15];
	double &I_Na_Ca = currents[16], &I_Na_K = currents[17];
	double &I_p_Ca	= currents[18];
	double &I_Na_ion_total	= currents[19];
	double &I_Ca_ion_total	= currents[20];
	double &I_K_ion_total	= currents[21];
	double &I_tr	= currents[22], &I_leak = currents[23];
	double &I_up	= currents[24], &I_rel = currents[25];

	// RATES
	// Na channel markov state model
	double &dIC3	= rates[ 0], &dIC2	= rates[1];
	double &dIF		= rates[ 2], &dIM1	= rates[3];
	double &dIM2	= rates[ 4], &dC3	= rates[5];
	double &dC2		= rates[ 6], &dO	= rates[7];
	double &dOS		= rates[ 8], &dC1	= rates[9];
	double &dDIC3	= rates[10], &dDIC2 = rates[11];
	double &dDIF	= rates[12], &dDIM1 = rates[13];
	double &dDIM2	= rates[14], &dDC3	= rates[15];
	double &dDC2	= rates[16], &dDO	= rates[17];
	double &dDOS	= rates[18], &dDC1	= rates[19];
	double &dD_IC3	= rates[20], &dD_IC2 = rates[21];
	double &dD_IF	= rates[22], &dD_IM1 = rates[23];
	double &dD_IM2	= rates[24], &dD_C3	= rates[25];
	double &dD_C2	= rates[26], &dD_O	= rates[27];
	double &dD_OS	= rates[28], &dD_C1	= rates[29];

	// Gating variable rates
	double &dm		= rates[ngatem     ],	&dh		= rates[ngatem +  1];
	double &dj		= rates[ngatem +  2];
	double &dd		= rates[ngatem +  3],	&df		= rates[ngatem +  4];
	double &df2		= rates[ngatem +  5],	&df_Ca	= rates[ngatem +  6];
	double &dr		= rates[ngatem +  7],	&ds		= rates[ngatem +  8];
	double &dxr1	= rates[ngatem +  9],	&dxr2	= rates[ngatem + 10];
	double &dxs		= rates[ngatem + 11];
	double &dR_bar	= rates[ngatem + 12];
	double &dmL		= rates[ngatem + 13],	&dhL	= rates[ngatem + 14];
	double &dCa_in	= rates[ngatem + 15],	&dCa_sr	= rates[ngatem + 16];
	double &dCa_ss	= rates[ngatem + 17];
	double &dNa_in	= rates[ngatem + 18],	&dK_in	= rates[ngatem + 19];


	// Cell data points - markers of cell function
	double I_total_old, dI_total, dI_total_old, t_I_total, I_total_pt;
	double V_I_total;
	bool   flagI_total, flag_EAD, flag_EAD2;
	double t_90_old;
	double CV, CV_EAD;

	// The temporary array is used for a modified-Newton iteration to solve
	// for the markov channel state variables.
	static const unsigned int	i3  = 0,  i2  = 1,	iF  = 2,  im  = 3,
								iM  = 4,  c3  = 5,	c2  = 6,  op  = 7,
								oS  = 8,  c1  = 9;
	static const unsigned int	di3 = 10, di2 = 11,	diF = 12, dim = 13,
								diM	= 14, dc3 = 15,	dc2	= 16, dop = 17,
								doS	= 18, dc1 = 19;
	static const unsigned int	Di3	= 20, Di2 = 21,	DiF	= 22, Dim = 23,
								DiM	= 24, Dc3 = 25,	Dc2	= 26, Dop = 27,
								DoS	= 28, Dc1 = 29;

	// BLAS vectors and methods to allocate, free and copy them
	gsl_vector* Yn   = nullptr;
	gsl_vector* Ynph = nullptr;
	gsl_vector* Ynp1 = nullptr;
	gsl_vector* Ynp1i= nullptr;
	gsl_vector* DiagonalAn = nullptr;
	gsl_vector* DiagonalAnInv = nullptr;
	gsl_matrix* OffDiagonalAn = nullptr;
	void allocTemps();
	void freeTemps();
	void copyTemps(const Moreno2011Cell& source);



	/*********************************************************************************************************
	 Current Function Prototypes
	 *********************************************************************************************************/
	void WT_SCN5A_Initial_Conditions();
	void WT_SCN5A_Lidocaine();
	void WT_SCN5A_Flecainide();
	void checkNaChannelStates();
	void calcI_total(double t);
	void calcI_Na(double t);	//Fast sodium current
	void calcI_Na_L();
	void calcI_K1();			//Time independent potassium current
	void calcI_Kp();			//Plateau K current, time independent, K_out insensitive
	void calcI_Na_Ca();			//Sodium calcium exchanger
	void calcI_Na_K();			//Sodium potassium pump
	void calcI_Ca_L();			//L type calcium channel Ca contribution
	void calcI_Kr();			//Rapid rectifier
	void calcI_Ks();			//Slow rectifier
	void calcI_to();			//Transient outward current
	void calcI_p_Ca();			//Sarcolemmal calcium pump
	void calcI_Ca_b();			//Calcium background current
	void calcI_Na_b();			//Sodium background current

	void calcI_tr();
	void calcI_leak();			//Leak from SR
	void calcI_up();			//Uptake to SR
	void calcI_rel();			//Release from SR

	//Dynamic Concentrations
	void Calculate_Na_in();			//Dynamic [Na] in Myoplasm
	void Calculate_K_in();			//Dynamic [K] in Myoplasm
	void Calculate_Ca_in();			//Dynamic [Ca] in Myoplasm
	void Calculate_Ca_sr();			//Dynamic [Ca] in JSR
	void Calculate_Ca_ss();			//Dynamic [Ca] in NSR

	// Semi-implicit and Rush-Larsen
	struct NaChannelParameters{
		double pKa, portion, diffusion,
				drug_charged, drug_neutral,
				drug_distance, kd0, kd_open;
		double	a11, a12, a13, b11, b12, b13,
				a3,	 b3,  a2,  b2,  a4,  b4, a5, b5,
				ax,	 bx,  ax1, bx1,
				a13c,a22, b33, a33,
				a44, b44, a55, b55,
				ax2, a13n, a_22, b_33,
				a_44,b_44, a_55, b_55,
				kon, koff, kcoff, kcon,
				b13c,b22,  k_on,  k_off,
				ki_on, ki_off, kc_on, kc_off,
				a_33, b13n,	b_22, bx2;
		void Initialize( cellmodel::Cells* cellenv );
		void ComputeVDepTerms(const double V, cellmodel::Cells* env );
		void set_NaCh_Matrices( gsl_matrix* M, gsl_vector* D);
	};
	NaChannelParameters NaChParams;
	static constexpr double tolerance = 1e-7;
	void semi_impl_I_Na(double t, const double dt);	//Fast sodium current
	void semi_impl_WT_SCN5A(const double dt);
	void semi_impl_WT_SCN5A_LU(const double dt);
	void exponential_I_Na1(double t, const double dt);
	void exponential_I_Na2(double t, const double dt, const double Vhalf);

	void rush_I_Na_L(const double dt);
	void rush_I_Ca_L(const double dt);			//L type calcium channel Ca contribution
	void rush_I_Kr(const double dt);			//Rapid rectifier
	void rush_I_Ks(const double dt);			//Slow rectifier
	void rush_I_to(const double dt);			//Transient outward current

	void rush_I_rel(const double dt);			//Release from SR
	void rush_Na_in(const double dt);			//Dynamic [Na] in Myoplasm
	void rush_K_in(const double dt);			//Dynamic [K] in Myoplasm
	void rush_Ca_in(const double dt);			//Dynamic [Ca] in Myoplasm
	void rush_Ca_sr(const double dt);			//Dynamic [Ca] in JSR
	void rush_Ca_ss(const double dt);			//Dynamic [Ca] in NSR


	// Old code - not used and not needed - new code is working
	void semi_impl_WT_SCN5A_Lidocaine(const double dt);
	void semi_impl_WT_SCN5A_Flecainide(const double dt);
	void set_NaCh_Matrices( gsl_matrix* M, gsl_vector* D,
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
			const double a_33,	const double b13n,	const double b_22,	const double bx2 );
	double MCtemp[nNaChStates];
	double &mIC3	= MCtemp[ 0], &mIC2 	= MCtemp[1];
	double &mIF		= MCtemp[ 2], &mIM1		= MCtemp[3];
	double &mIM2	= MCtemp[ 4], &mC3		= MCtemp[5];
	double &mC2		= MCtemp[ 6], &mO		= MCtemp[7];
	double &mOS		= MCtemp[ 8], &mC1		= MCtemp[9];
	double &mDIC3	= MCtemp[10], &mDIC2	= MCtemp[11];
	double &mDIF	= MCtemp[12], &mDIM1	= MCtemp[13];
	double &mDIM2	= MCtemp[14], &mDC3		= MCtemp[15];
	double &mDC2	= MCtemp[16], &mDO		= MCtemp[17];
	double &mDOS	= MCtemp[18], &mDC1		= MCtemp[19];
	double &mD_IC3	= MCtemp[20], &mD_IC2	= MCtemp[21];
	double &mD_IF	= MCtemp[22], &mD_IM1	= MCtemp[23];
	double &mD_IM2	= MCtemp[24], &mD_C3	= MCtemp[25];
	double &mD_C2	= MCtemp[26], &mD_O		= MCtemp[27];
	double &mD_OS	= MCtemp[28], &mD_C1	= MCtemp[29];
	void semi_impl_WT_SCN5A_Lidocaine_compare_BLAS_old(const double dt);
	void semi_impl_WT_SCN5A_Flecainide_old (const double dt);
};

#endif /* CELL_H_ */

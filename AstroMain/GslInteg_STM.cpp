/*% © 2014 Kathleen Howell All Rights Reserved.
% The tool ATD (Adaptive Trajectory Design) was developed at Purdue University in the 
% School of Aeronautics and Astronautics as a collaborative effort between Purdue University 
% and the NASA Goddard Space Flight Center. The project was initiated under NASA Grant Nos. NNX12AC57G and NNX13AE55G. 
% Significant contributors include Professor Kathleen Howell and 
% Ph.D. students Amanda Haapala and Tom Pavlak.

%Written by Amanda Haapala

 * call as [sol,t] = GslInteg_STM(IC,T,mu)
 * IC - 42-element rotating, barycentric initial state written as a row vector
 * T - nondimentional time interval; can be in the form [t0,tf] or [t0,t1,t2,...,tf]. 
 * mu - CR3BP mass paramter value.
*/
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <matrix.h>
#include <mex.h>   
#include <tnt.h>
#include <jama_eig.h>
#include <vector>
#include <stdio.h>

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

int nsteps; //declare number of integration steps as global variable 

int func_eoms (double t, const double f[], double fp[], void *params){
	double mu = *(double *)params;

	double x;
	double y;
    double z;
	double xdot;
	double ydot;
    double zdot;
	double d;
	double r;
    double Uxx,Uxy,Uxz,Uyy,Uyz,Uzz;
    
    double  phi[6][6];
    double  A[6][6];
    double  phidot[6][6];

    //define state vector variables
	x       =   f[0];
	y       =   f[1];
    z       =   f[2];
	xdot    =   f[3];
	ydot    =   f[4];
    zdot    =   f[5];
    
	d       =   sqrt((x+mu)*(x+mu) + y*y + z*z); 
	r       =   sqrt((x-1+mu)*(x-1+mu) + y*y + z*z);

	fp[0]   =   xdot;
	fp[1]   =   ydot;
    fp[2]   =   zdot;
	fp[3]   =   2*ydot + x - (1-mu)*(x+mu)/(d*d*d) - mu*(x-1+mu)/(r*r*r);
	fp[4]   =  -2*xdot + y - (1-mu) * y/(d*d*d) - mu*y/(r*r*r);
    fp[5]   =  -(1-mu)*z/(d*d*d) - mu*z/(r*r*r); 
        
    //define elements of STM, phi
    for (int k=0; k<6; k++){
		for (int j=0; j<6; j++){
            phi[k][j] = f[6+(6*k)+j]; 
        }
	}
    
    //compute terms for A, phidot = A*phi
    Uxx     = 1 - ((1-mu)/pow(d,3)) - (mu/(pow(r,3))) + (3*(1-mu)*((x+mu)*(x+mu))/(pow(d,5))) + (3*mu*((x-1+mu)*(x-1+mu))/(pow(r,5)));
    Uyy     = 1 - ((1-mu)/(pow(d,3))) - (mu/(pow(r,3))) + (3*(1-mu)*(y*y)/(pow(d,5))) + (3*mu*(y*y)/(pow(r,5)));
    Uzz     = -((1-mu)/(pow(d,3))) - (mu/(pow(r,3))) + (3*(1-mu)*(z*z)/(pow(d,5))) + (3*mu*(z*z)/(pow(r,5)));
    Uxy     = (3*(1-mu)*(x+mu)*y/(pow(d,5))) + (3*mu*(x-1+mu)*y/(pow(r,5)));
    Uxz     = (3*(1-mu)*(x+mu)*z/(pow(d,5))) + (3*mu*(x-1+mu)*z/(pow(r,5)));
    Uyz     = (3*(1-mu)*y*z/(pow(d,5))) + (3*mu*y*z/(pow(r,5)));

    A[0][0] = 0; A[0][1] = 0; A[0][2] = 0;
    A[1][0] = 0; A[1][1] = 0; A[1][2] = 0;
    A[2][0] = 0; A[2][1] = 0; A[2][2] = 0;
    
    A[0][3] = 1; A[0][4] = 0; A[0][5] = 0;
    A[1][3] = 0; A[1][4] = 1; A[1][5] = 0;
    A[2][3] = 0; A[2][4] = 0; A[2][5] = 1;
    
    A[3][0] = Uxx; A[3][1] = Uxy; A[3][2] = Uxz;
    A[4][0] = Uxy; A[4][1] = Uyy; A[4][2] = Uyz;
    A[5][0] = Uxz; A[5][1] = Uyz; A[5][2] = Uzz;
    
    A[3][3] = 0; A[3][4] = 2; A[3][5] = 0;
    A[4][3] = -2; A[4][4] = 0; A[4][5] = 0;
    A[5][3] = 0; A[5][4] = 0; A[5][5] = 0;

    for (int k=0; k<6; k++){
		for (int j=0; j<6; j++){
            phidot[k][j] = A[k][0]*phi[0][j]+A[k][1]*phi[1][j]+A[k][2]*phi[2][j]+A[k][3]*phi[3][j]+A[k][4]*phi[4][j]+A[k][5]*phi[5][j];
        }
	}
    
    for (int k=0; k<6; k++){
        for (int j=0; j<6; j++){
            fp[6+(6*k)+j] = phidot[k][j];
        }
	}
    
	return GSL_SUCCESS;
}
//----------------------------------------------------------------------------------------------------------------------------------

int jac (void)
{
	return GSL_SUCCESS;
}
//----------------------------------------------------------------------------------------------------------------------------------
vector<double> prop(double ic[], double time[], double mu[], int state_dim, int t_dim){   //, std::vector<std::vector<double>> state)

	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
// 	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkck;
//     const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;
//     const gsl_odeiv2_step_type * T = gsl_odeiv2_step_msadams;
    
    double mu0 = mu[0];
    int i = 0; //initialize number of integration steps
    vector<double> state; //initialize dynamic array for state data storage
//     double y[6]; //initialize state vector
    double *y = new double[state_dim];
    double *times = new double[t_dim];

    gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, state_dim);
    gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1e-12, 1e-14); //set tolerances
    gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (state_dim);
    gsl_odeiv2_system sys = {func_eoms, NULL, static_cast<size_t>(state_dim), &mu0};
//     gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, T, 1e-6, 1e-12, 1e-14);
    
    //set y = initial state
    for (int k=0; k<state_dim; k++){
        y[k] = ic[k];
    }
    
    //set times = time vector
    for (int k=0; k<t_dim; k++){
        times[k] = time[k];
    }
    
    for (int k=0;k<state_dim+1;k++){
        if (k<state_dim){
            state.push_back(y[k]);
        }
        else if (k==state_dim){
            state.push_back(times[0]);

        }
    }
    i++;
    
    if (t_dim==2){
        double t = times[0], t1 = times[1]; //define start and end times for integration
        double h = (t1>t) ? 1e-12 : -1e-12; //set step size initial guess (can be positive or negative)
        
        int sgn = (int)((t1-t)/fabs(t1-t));

//         while (fabs(t) < fabs(t1)) {
        while (sgn*t < sgn*t1) {
            int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
//             int status = gsl_odeiv2_driver_apply (d, &t, t1, y);

            if (status != GSL_SUCCESS){ //check for successful integration step
                break;}

            i++;//update number of integration steps

            //concatenate current state,time in state vector
            for (int k=0;k<state_dim+1;k++){
                if (k<state_dim){
                    state.push_back(y[k]);
                }
                else if (k==state_dim){
                    state.push_back(t);

                }
            }
        }
    }
    else {
        double t, t1;
        for (int tcount=0; tcount<t_dim-1; tcount++){
            t = time[tcount]; 
            t1 = time[tcount+1]; //define start and end times for integration
            double h = (t1>t) ? 1e-12 : -1e-12; //set step size initial guess (can be positive or negative)

            int sgn = (int)((t1-t)/fabs(t1-t));
            while (sgn*t < sgn*t1) {
                int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
//                 int status = gsl_odeiv2_driver_apply (d, &t, t1, y);

                if (status != GSL_SUCCESS){ //check for successful integration step
                    break;}
            }
            //concatenate current state,time in state vector
            for (int k=0;k<state_dim+1;k++){
                if (k<state_dim){
                    state.push_back(y[k]);
                }
                else if (k==state_dim){
                    state.push_back(t);

                }
            }
            i++;//update number of integration steps
        }

    }
    
//     else {
//         for (int tcount=0; tcount<t_dim; tcount++){
//             double t = time[tcount], t1 = time[tcount]; //define start and end times for integration
//             double h = (t1>t) ? 1e-12 : -1e-12; //set step size initial guess (can be positive or negative)
// 
//             while (fabs(t) < fabs(t1)) {
//                 int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
// 
//                 if (status != GSL_SUCCESS){ //check for successful integration step
//                     break;}
//             }
//             //concatenate current state,time in state vector
//             for (int k=0;k<state_dim+1;k++){
//                 if (k<state_dim){
//                     state.push_back(y[k]);
//                 }
//                 else if (k==state_dim){
//                     state.push_back(t);
// 
//                 }
//             }
//             i++;//update number of integration steps
//         }
//     }
    
    //free memory allocated to e, c, s
    gsl_odeiv2_evolve_free (e); 
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);
  
    nsteps = i; //set number of integrator steps globally
    
    delete(y);
    delete(times);
    
	return state; //return the array of states and times back to mex function
}
//----------------------------------------------------------------------------------------------------------------------------------

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    //declare variables
    mxArray *c_out_m, *t_out_m;
    const mwSize *dims, *time_dims;
    double *a, *c, *t, *time, *mu;
    int dim, time_dim;
    vector<double> state;

    //figure out dimensions
    dims = mxGetDimensions(prhs[0]); //get dimensions of state input vector
    time_dims = mxGetDimensions(prhs[1]); //get dimensions of time input vector
	
	if ((int)dims[0] != 1) {
		printf("input vector should be 1xn\n");
		return;
	}
    if ((int)time_dims[0] != 1) {
		printf("input vector should be 1xm\n");
		return;
	}
    
	dim = (int)dims[1]; //get dimension of state vector
    time_dim = (int)time_dims[1]; //get number of time inputs -> 2: [t0,tf], >2: [t0,t1,...,tf]
	
    //associate pointers
	a = mxGetPr(prhs[0]);
    time = mxGetPr(prhs[1]);
    mu = mxGetPr(prhs[2]);
    
    if (dim==42){
        state = prop(a,time,mu,dim,time_dim);
    }
    else {
        printf("State vector should have 42 columns");
    }

//     associate outputs
    c_out_m = plhs[0] = mxCreateDoubleMatrix(nsteps,dim,mxREAL);
    t_out_m = plhs[1] = mxCreateDoubleMatrix(nsteps,1,mxREAL);
    
    //associate pointers
    c = mxGetPr(c_out_m);
    t = mxGetPr(t_out_m);

    
	for (int k=0; k<nsteps; k++){
        for (int j=0; j<dim+1; j++){
            if (j<dim){
                c[j*nsteps+k] = state[(dim+1)*k+j];
            }
            else if (j==dim){
                t[k] = state[(dim+1)*k+j];
            }
        }
	}
    


    return ;
}
/*% © 2014 Kathleen Howell All Rights Reserved.
% The tool ATD (Adaptive Trajectory Design) was developed at Purdue University in the 
% School of Aeronautics and Astronautics as a collaborative effort between Purdue University 
% and the NASA Goddard Space Flight Center. The project was initiated under NASA Grant Nos. NNX12AC57G and NNX13AE55G. 
% Significant contributors include Professor Kathleen Howell and 
% Ph.D. students Amanda Haapala and Tom Pavlak.

 * Written by Amanda Haapala
 * Updated by Andrew Cox
 * 
 * GslInteg_CR3BP.cpp
 * call as [sol,t] = GslInteg_CR3BP(IC,T,mu)
 * IC - 6-element rotating, barycentric initial state written as a row vector
 * T - nondimentional time interval; can be in the form [t0,tf] or [t0,t1,t2,...,tf]. 
 * mu - CR3BP mass paramter value.
 * 
 * To compile/mex on OS X, run the following command in Matlab
 *  mex GslInteg_CR3BP.cpp -I../include -L../libmac -lgsl -lgslcblas
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

int nsteps = 0; //declare number of integration steps as global variable 

int func_CR3BP_eoms (double t, const double f[], double fp[], void *params){
	double mu = *(double *)params;

	double x;
	double y;
    double z;
	double xdot;
	double ydot;
    double zdot;
	double d;
	double r;
    
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
	fp[3]   =   2*ydot + x - (1-mu)*(x+mu)/pow(d,3) - mu*(x-1+mu)/pow(r,3);
	fp[4]   =  -2*xdot + y - (1-mu) * y/pow(d,3) - mu*y/pow(r,3);
    fp[5]   =  -(1-mu)*z/pow(d,3) - mu*z/pow(r,3); 
        
	return GSL_SUCCESS;
}//==============================================

int jac (){
	return GSL_SUCCESS;
}//==============================================

vector<double> prop(double ic[], double in_times[], double mu[], int state_dim, int t_dim){ 

	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;

    double mu0 = mu[0];
    int i = 0;                  //initialize number of integration steps
    vector<double> state;       //initialize dynamic array for state data storage

    // Copy IC and times into vector
    vector<double> y_vec (ic, ic + state_dim);
    vector<double> times_vec (in_times, in_times+t_dim);

    // Get pointers to the start of both vectors
    double *y = &(y_vec[0]);
    double *times = &(times_vec[0]);

    gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, state_dim);
    gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1e-12, 1e-14); //set tolerances
    gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (state_dim);
    gsl_odeiv2_system sys = {func_CR3BP_eoms, NULL, static_cast<size_t>(state_dim), &mu0};
    
    //add initial state and time to state vector
    state.insert(state.end(), y, y+state_dim);
    state.insert(state.end(), in_times[0]);
    i++;
    
    if (t_dim==2){
        double t = times[0], t1 = times[1]; //define start and end times for integration
        double h = t1 > t ? 1e-12 : -1e-12; //set step size initial guess (can be positive or negative)

        int sgn = (int)((t1-t)/fabs(t1-t));


        while (sgn*t < sgn*t1) {
            int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);

            if (status != GSL_SUCCESS){ //check for successful integration step
                printf("integration terminating");
                break;
            }

            //concatenate current state,time in state vector
            state.insert(state.end(), y, y+state_dim);
            state.insert(state.end(), t);
            i++; //update number of integration steps
        }
    } else {
        double t, t1;
        for (int tcount = 0; tcount < t_dim-1; tcount++){
            t = in_times[tcount];
            t1 = in_times[tcount+1]; //define start and end times for integration
            double h = t1 > t ? 1e-12 : -1e-12; //set step size initial guess (can be positive or negative)

            int sgn = (int)((t1-t)/fabs(t1-t));

            while (sgn*t < sgn*t1) {                
                int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
                
                //check for successful integration step
                if (status != GSL_SUCCESS){
                    break;
                }
            }

            //concatenate current state, time in state vector
            state.insert(state.end(), y, y+state_dim);
            state.insert(state.end(), t);
            i++; //update number of integration steps
        }
    }
    
    //free memory allocated to e, c, s
    gsl_odeiv2_evolve_free (e); 
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);
  
    nsteps = i; //set number of integrator steps globally
    
	return state; //return the array of states and times back to mex function
}
//----------------------------------------------------------------------------------------------------------------------------------

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    //declare variables
    mxArray *c_out_m, *t_out_m;
    const mwSize *dims, *time_dims;
    double *a, *c, *t, *times, *mu;
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
    times = mxGetPr(prhs[1]);
    mu = mxGetPr(prhs[2]);
    
    if (dim == 6){
        state = prop(a,times,mu,dim,time_dim);
    } else {
        printf("State vector should have 6 columns");
    }

    //associate outputs
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
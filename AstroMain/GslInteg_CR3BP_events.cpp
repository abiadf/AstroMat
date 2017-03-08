/*% © 2014 Kathleen Howell All Rights Reserved.
% The tool ATD (Adaptive Trajectory Design) was developed at Purdue University in the 
% School of Aeronautics and Astronautics as a collaborative effort between Purdue University 
% and the NASA Goddard Space Flight Center. The project was initiated under NASA Grant Nos. NNX12AC57G and NNX13AE55G. 
% Significant contributors include Professor Kathleen Howell and 
% Ph.D. students Amanda Haapala and Tom Pavlak.

 * Written by Amanda Haapala
 * Modified by Andrew Cox
 *
 * GslInteg_CR3BP_events.cpp
 * call as [c,t,cev,tev,iev] = GslInteg_CR3BP_events(IC,T,mu,event_id,event_state,dir,stop);
 * IC - 6-element rotating, barycentric initial state written as a row vector
 * T - nondimentional time interval; can be in the form [t0,tf] or [t0,t1,t2,...,tf]
 * mu - CR3BP mass paramter value
 * event_id - row vector that identifies the event type: 
 * 1 - x=event_state[0]; 
 * 2 - y=event_state[1];
 * 3 - z=event_state[2];
 * 4 - xdot=event_state[3];
 * 5 - ydot=event_state[4];
 * 6 - zdot=event_state[5];
 * 7 - apses relative to P1
 * 8 - apses relative to P2
 * 9 - angle relative to P1, radians
 * 10 - angle relative to P2, radians
 * 11 - radius relative to P1
 * 12 - radius relative to P2
 * event_state - row vector defining desired event states
 * dir - row vector indicating desired event crossing directions (+1 - positive crossing, -1 - negative crossing, 0 - all crossings).
 * stop - row vector indicating whether integration should terminate at the events
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
#define _USE_MATH_DEFINES

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
	fp[3]   =   2*ydot + x - (1-mu)*(x+mu)/(d*d*d) - mu*(x-1+mu)/(r*r*r);
	fp[4]   =  -2*xdot + y - (1-mu) * y/(d*d*d) - mu*y/(r*r*r);
    fp[5]   =  -(1-mu)*z/(d*d*d) - mu*z/(r*r*r); 
        
	return GSL_SUCCESS;
}//===========================================

int jac (void){ return GSL_SUCCESS; }

double eventcheck(double y[], int event_id, double event_state_j, double mu0){

    double event;
    if (event_id < 7){
        event = y[event_id-1] - event_state_j;
    } else if (event_id == 7){//P1 apse
        event = (y[0]+mu0)*y[3] + y[1]*y[4] + y[2]*y[5];
    } else if (event_id == 8){//P2 apse
        event = (y[0]-1+mu0)*y[3] + y[1]*y[4] + y[2]*y[5];
    } else if (event_id == 9){//P1 angle
        double pi = 3.141592653589793;
        double x0 = y[0]+mu0;
        double y0 = y[1];
        
        if (abs(event_state_j-pi/2)<1e-5 || abs(event_state_j+pi/2)<1e-5 || abs(event_state_j-3*pi/2)<1e-5) {
            event = x0;
        } else {
            event = y0 - x0*tan(event_state_j);
        }
                
    } else if (event_id == 10){//P2 angle
        double pi = 3.141592653589793;
        double x0 = y[0]-1+mu0;
        double y0 = y[1];
        
        if (abs(event_state_j-pi/2)<1e-5 || abs(event_state_j+pi/2)<1e-5 || abs(event_state_j-3*pi/2)<1e-5) {
            event = x0;
        }
        else {
            event = y0 - x0*tan(event_state_j);
        }
                
    } else if (event_id == 11){ // P1 Distance
        double r = sqrt(pow(y[0]+mu0,2)+pow(y[1],2)+pow(y[2],2));
        event = r - event_state_j;
    } else if (event_id == 12){ // P2 Distance
        double r = sqrt(pow(y[0]-1+mu0,2)+pow(y[1],2)+pow(y[2],2));
        event = r - event_state_j;
    }

    return event;
}//====================================================

vector<double> newtonraphson(double y[], double t, int state_dim, int event_id, double event_state_j, double mu0){

    int count = 0;
    double err = 1;
    double tol = 1e-12, h;
    double y0[6], y_in[6];
    double evn0, evn_old, evn2, evn_new, deltat, devndt;
    double t_in = t, ti, tf;
    int status, k;
    double step = 1e-6;
    vector<double> event; //initialize dynamic array for event data storage
    
    const gsl_odeiv2_step_type * T_nr = gsl_odeiv2_step_rk8pd; //

    gsl_odeiv2_step * s_nr = gsl_odeiv2_step_alloc (T_nr, state_dim);
    gsl_odeiv2_control * c_nr = gsl_odeiv2_control_y_new (1e-12, 1e-14); //set tolerances
    gsl_odeiv2_evolve * e_nr = gsl_odeiv2_evolve_alloc (state_dim);
    gsl_odeiv2_system sys_nr = {func_CR3BP_eoms, NULL, static_cast<size_t>(state_dim), &mu0};

    for (k=0; k<state_dim; k++){
        y_in[k] = y[k]; //set y0 = input state
        y0[k] = y[k];
    }
    
    while (err>tol && count<10 && err<=1){

        //define event at previous guess
        evn_old = eventcheck(y_in,event_id,event_state_j,mu0); //check event
        
        //reset y0 
        for (k=0; k<state_dim; k++){
            y0[k] = y_in[k]; //set y0 = input state
        }

        //integrate y0 for t+step and evaluate event
     	ti = t_in, tf = t_in+step; //define start and end times for integration
        h = (tf>ti) ? 1e-12 : -1e-12; //set step size initial guess (can be positive or negative)
        int sgn = (int)((tf-ti)/fabs(tf-ti));
        while (sgn*ti<sgn*tf) {
            status = gsl_odeiv2_evolve_apply (e_nr, c_nr, s_nr, &sys_nr, &ti, tf, &h, y0);

            if (status != GSL_SUCCESS){ //check for successful integration step
                break;
            }
        }
        evn2 = eventcheck(y0, event_id, event_state_j, mu0); //check event

        //reset y0 
        for (k=0; k<state_dim; k++){
            y0[k] = y_in[k]; //set y0 = input state
        }

        //integrate y0 for t-step and evaluate event
        ti = t_in, tf = t_in-step; //define start and end times for integration
        h = (tf>ti) ? 1e-12 : -1e-12; //set step size initial guess (can be positive or negative)
        sgn = (int)((tf-ti)/fabs(tf-ti));
        
        while (sgn*ti<sgn*tf) {
            status = gsl_odeiv2_evolve_apply (e_nr, c_nr, s_nr, &sys_nr, &ti, tf, &h, y0);
            if (status != GSL_SUCCESS){ //check for successful integration step
                break;
            }
        }      
        evn0 = eventcheck(y0,event_id,event_state_j,mu0); //check event

        //get new t for integration via 1st order Taylor series expansion
        devndt = (evn2-evn0)/(2*step); //compute numerical derivative
        deltat = -evn_old/devndt;
        t = t + deltat;
        
        //integrate original y0 to new t and evaluate event function again
        h = (t>t_in) ? 1e-12 : -1e-12; //set step size initial guess (can be positive or negative)
        sgn = (int)((t-t_in)/fabs(t-t_in));
        
        while (sgn*t_in<sgn*t) {
            status = gsl_odeiv2_evolve_apply (e_nr, c_nr, s_nr, &sys_nr, &t_in, t, &h, y_in);
            if (status != GSL_SUCCESS){ //check for successful integration step
                break;}
        }        
        evn_new = eventcheck(y_in,event_id,event_state_j,mu0); //check event

        err = fabs(evn_new); //exit loop when event = 0
        
        count++;
    }

    if (err > 1e-12){
        t_in = 50000;
    }

    for (k=0;k<state_dim+1;k++){
        if (k<state_dim){
            event.push_back(y_in[k]);
        }
        else {
            event.push_back(t_in);
        }
    }

    //free memory allocated to e, c, s
    gsl_odeiv2_evolve_free (e_nr); 
    gsl_odeiv2_control_free (c_nr);
    gsl_odeiv2_step_free (s_nr);
    
    return event;
}
//----------------------------------------------------------------------------------------------------------------------------------
double bisection(double y1[], double y2[], double t1, double t2, int state_dim, int event_id, double event_state_j, double mu0){
    
    int count = 0;
    double err = 1;
    double tol = 1e-12, h;
    double ti, tf;
    double evn1, evn2, evn_mid;
    double t_mid, y_mid[6];
    int status, k;
    double step = 1e-6;
    
    const gsl_odeiv2_step_type * T_nr = gsl_odeiv2_step_rk8pd; //
    
    gsl_odeiv2_step * s_nr = gsl_odeiv2_step_alloc (T_nr, state_dim);
    gsl_odeiv2_control * c_nr = gsl_odeiv2_control_y_new (1e-12, 1e-14); //set tolerances
    gsl_odeiv2_evolve * e_nr = gsl_odeiv2_evolve_alloc (state_dim);
    gsl_odeiv2_system sys_nr = {func_CR3BP_eoms, NULL, static_cast<size_t>(state_dim), &mu0};
    
    evn1 = eventcheck(y1,event_id,event_state_j,mu0); //check event
    evn2 = eventcheck(y2,event_id,event_state_j,mu0); //check event
        
    while (err>tol && count<1000){
        t_mid = (t1+t2)/2;

        //reset y_mid to start integration 
        for (k=0; k<state_dim; k++){
            y_mid[k] = y1[k]; //set y_mid = y1 state
        }
        
        //integrate y1 to t_mid to locate y_mid
        ti = t1, tf = t_mid; //define start and end times for integration
        h = (tf>ti) ? 1e-12 : -1e-12; //set step size initial guess (can be positive or negative)
        int sgn = (int)((tf-ti)/fabs(tf-ti));
        while (sgn*ti<sgn*tf) {
            status = gsl_odeiv2_evolve_apply (e_nr, c_nr, s_nr, &sys_nr, &ti, tf, &h, y_mid);

            if (status != GSL_SUCCESS){ //check for successful integration step
                break;
            }
        }

        //define event at t1, t_mid, t2
        evn1 = eventcheck(y1,event_id,event_state_j,mu0); //check event
        evn_mid = eventcheck(y_mid,event_id,event_state_j,mu0); //check event
        evn2 = eventcheck(y2,event_id,event_state_j,mu0); //check event

        if (evn1*evn2<0){
            if (evn_mid*evn1>0){//if evn_mid and evn1 same sign, set t1 to t_mid
                t1 = t_mid;
                for (k=0; k<state_dim; k++){
                    y1[k] = y_mid[k];
                }
            }
            else {//if (evn_mid*evn2>0){//if evn_mid and evn2 same sign, set t2 to t_mid
                t2 = t_mid;
                for (k=0; k<state_dim; k++){
                    y2[k] = y_mid[k];
                }            
            }
        }

        err = fabs(evn_mid); //exit loop when event = 0

        count++;
    }
       
    if (err>1e-12){
        t_mid = 50000;//if tolerance not met, set output to indicate this
    }

    //free memory allocated to e, c, s
    gsl_odeiv2_evolve_free (e_nr); 
    gsl_odeiv2_control_free (c_nr);
    gsl_odeiv2_step_free (s_nr);
    
    return t_mid;
}//===================================================================

vector<double> prop(double ic[], double time[], double mu[], int event_id[], double event_state[], int event_dir[], int event_stop[], int state_dim, int t_dim, int evn_dim){
          
	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;

    double mu0 = mu[0];
    vector<double> state; //initialize dynamic array for state data storage
    double y[6];
    double y_event[6];
    double t, t_event, ti, tf, tfin;
    int findroot, k, j;
    int quitflag = 0;
    double h, h_event;
    double t1, y1[6];//for bisection method
    vector<double> state_events, u_event;//initialize dynamic array for state data storage

    double *evn_prev = new double[evn_dim];
    double *evn = new double[evn_dim];
    int *evn_reset = new int[evn_dim];
    
    gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, state_dim);
    gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1e-12, 1e-14); //set tolerances
    gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (state_dim);
    gsl_odeiv2_system sys = {func_CR3BP_eoms, NULL, static_cast<size_t>(state_dim), &mu0};

    for (k=0; k<state_dim; k++){
        y[k] = ic[k]; //set y = initial state
    }

    //set initial state/time in state vector for output
    for (k=0;k<state_dim+2;k++){
        if (k==0){
            state.push_back(0);
        }
        else if (k>0 && k<state_dim+1){
            state.push_back(y[k-1]);
        }
        else {
            state.push_back(time[0]);
        }
    }
    
    if (t_dim==2){
        t = time[0], tf = time[1]; //define start and end times for integration
        h = (tf>t) ? 1e-12 : -1e-12; //set step size initial guess (can be positive or negative)
        int sgn = (int)((tf-t)/fabs(tf-t));
        while (sgn*t<sgn*tf && quitflag == 0) {

            for (j=0; j<evn_dim; j++){
                evn_prev[j] = eventcheck(y,event_id[j],event_state[j],mu0); //check events for initial condition
                if (fabs(evn_prev[j])<1e-12) {
                    evn_prev[j] = 0;
                }
            }
            
            t1 = t; //set t1 to time before integration
            for (k=0; k<state_dim; k++){
                y1[k] = y[k];//set y1 to state before integration
            }

            int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, tf, &h, y);

            if (status != GSL_SUCCESS){ //check for successful integration step
                quitflag = 1;
                break;
            }
            
            if (sqrt(pow(y[0]-1+mu0,2)+pow(y[1],2)+pow(y[2],2))<1e-6 || sqrt(pow(y[0]+mu0,2)+pow(y[1],2)+pow(y[2],2))<1e-6){
                quitflag = 1;
            }
            
            //check events ---------------------------------------
            for (j=0; j<evn_dim; j++){
                evn[j] = eventcheck(y,event_id[j],event_state[j],mu0); //check event
                if (fabs(evn[j])<1e-12) {
                    evn[j] = 0;
                }
            }
            
            for (j=0; j<evn_dim; j++){
                findroot = 0;
                if (evn_prev[j]*evn[j]<0) {
                     //Direction check
                     if (event_dir[j]==0) {
                         //Both Directions
                         findroot = 1;}
                     else if (event_dir[j] == -1 && evn[j]<evn_prev[j]) {
                         //Decreasing + to -
                         findroot = 1;}
                     else if (event_dir[j] == 1 && evn[j]>evn_prev[j]) {
                         //Increasing - to +
                         findroot = 1;}
                }
                if (findroot==1 && event_id[j]!=0){
                   u_event = newtonraphson(y, t, state_dim, event_id[j], event_state[j], mu0); 
                   t_event = u_event[state_dim];

                   if (t_event==50000){//if newton raphson method fails, implement bisection method
                       t_event = bisection(y1, y, t1, t, state_dim, event_id[j], event_state[j], mu0); 
                   }
                   if (t_event!=50000){
                        evn_reset[j] = 1;

                        for (k=0; k<state_dim; k++){
                            y_event[k] = u_event[k];
                        }

                        //concatenate current state,time in state vector
                        for (k=0;k<state_dim+2;k++){
                            if (k==0){
                                state.push_back(j+1);
                            }
                            else if (k>0 && k<state_dim+1){
                                state.push_back(y_event[k-1]);
                            }
                            else {
                                state.push_back(t_event);
                            }
                        }
                        if (event_stop[j] == 1){
                            t = t_event;
                            for (k=0; k<state_dim; k++){
                                y[k] = y_event[k];
                            }
                            quitflag = 1;
                        }
                   }
                }
            }
            //concatenate current state,time in state vector
            for (k=0;k<state_dim+2;k++){
                if (k==0){
                    state.push_back(0);
                }
                else if (k>0 && k<state_dim+1){
                    state.push_back(y[k-1]);
                }
                else {
                    state.push_back(t);
                }
            } 
        }
    }
    else {
        int tcount = 0;
        while (tcount<t_dim-1 && quitflag == 0){
            double t = time[tcount], tf = time[tcount+1]; //define start and end times for integration
            h = (tf>t) ? 1e-12 : -1e-12; //set step size initial guess (can be positive or negative)
            int sgn = (int)((tf-t)/fabs(tf-t));
            while (sgn*t<sgn*tf && quitflag == 0) {
                
                for (j=0; j<evn_dim; j++){
                    evn_prev[j] = eventcheck(y,event_id[j],event_state[j],mu0); //check events for initial condition
                    if (fabs(evn_prev[j])<1e-12 && evn_reset[j]==1){
                        evn_prev[j] = 0;
                    }
                }

                t1 = t; //set t1 to time before integration
                for (k=0; k<state_dim; k++){
                    y1[k] = y[k];//set y1 to state before integration
                }
                int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, tf, &h, y);

                if (status != GSL_SUCCESS){ //check for successful integration step
                    break;}
                if (sqrt(pow(y[0]-1+mu0,2)+pow(y[1],2)+pow(y[2],2))<1e-6 || sqrt(pow(y[0]+mu0,2)+pow(y[1],2)+pow(y[2],2))<1e-6){
                    quitflag = 1;
                }

                //check events ---------------------------------------
                for (j=0; j<evn_dim; j++){
                    evn[j] = eventcheck(y,event_id[j],event_state[j],mu0); //check event
                    if (fabs(evn[j])<1e-12 && evn_reset[j]==1){
                        evn[j] = 0;
                    }
                    evn_reset[j] = 0;
                }
                for (j=0; j<evn_dim; j++){
                    findroot = 0;
                    if (evn_prev[j]*evn[j]<0) {
                         //Direction check
                         if (event_dir[j]==0) {
                             //Both Directions
                             findroot = 1;}
                         else if (event_dir[j] == -1 && evn[j]<evn_prev[j]) {
                             //Decreasing + to -
                             findroot = 1;}
                         else if (event_dir[j] == 1 && evn[j]>evn_prev[j]) {
                             //Increasing - to +
                             findroot = 1;}
                    }
                    if (findroot==1 && event_id[j]!=0){

                       u_event = newtonraphson(y, t, state_dim, event_id[j], event_state[j], mu0); 
                       t_event = u_event[state_dim];
                       if (t_event==50000){//if newton raphson method fails, implement bisection method
                           t_event = bisection(y1, y, t1, t, state_dim, event_id[j], event_state[j], mu0); 
                       }
                       if (t_event!=50000){
                            evn_reset[j] = 1;
                            for (k=0; k<state_dim; k++){
                                y_event[k] = u_event[k];
                            }

                            //concatenate current state,time in state vector
                            for (k=0;k<state_dim+2;k++){
                                if (k==0){
                                    state.push_back(j+1);
                                }
                                else if (k>0 && k<state_dim+1){
                                    state.push_back(y_event[k-1]);
                                }
                                else {
                                    state.push_back(t_event);
                                }
                            }
                            if (event_stop[j] == 1){
                                t = t_event;
                                for (k=0; k<state_dim; k++){
                                    y[k] = y_event[k];
                                }
                                quitflag = 1;
                            }
                       }
                    }
                }
                
            }
            
            //concatenate current state,time in state vector
            for (int k=0;k<state_dim+2;k++){
                if (k==0){
                    state.push_back(0);
                }
                else if (k>0 && k<state_dim+1){
                    state.push_back(y[k-1]);
                }
                else {
                    state.push_back(t);
                }
            }
            tcount = tcount+1;
        }
    }

    
    //free memory allocated to e, c, s
    gsl_odeiv2_evolve_free (e); 
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);
    
    delete(evn_prev);
    delete(evn);
    delete(evn_reset);
        
    
	return state; //return the array of states and times back to mex function
}

//----------------------------------------------------------------------------------------------------------------------------------

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    //declare variables
    mxArray *c_out_m, *t_out_m, *cev_out_m, *tev_out_m, *iev_out_m;
    const mwSize *dims, *time_dims, *evn_dims;
    double *a, *time, *mu,  *evn_state; //inputs
    double *evn_id, *evn_dir, *evn_stop; //inputs
    double *c, *t, *cev, *tev, *iev; //outputs
    int dim, time_dim, evn_dim;//event_id, event_dir, event_stop;
    int k, j;
    vector<double> state; 

    //figure out dimensions
    dims = mxGetDimensions(prhs[0]); //get dimensions of state input vector
    time_dims = mxGetDimensions(prhs[1]); //get dimensions of time input vector
    evn_dims = mxGetDimensions(prhs[3]);
	
	dim = (int)dims[1]; //get dimension of state vector
    time_dim = (int)time_dims[1]; //get number of time inputs -> 2: [t0,tf], >2: [t0,t1,...,tf]
	evn_dim = (int)evn_dims[1]; //get number of events to be evaluated
    
    //associate pointers
	a = mxGetPr(prhs[0]);
    time = mxGetPr(prhs[1]);
    mu = mxGetPr(prhs[2]);
    evn_id = mxGetPr(prhs[3]);
    evn_state = mxGetPr(prhs[4]);
    evn_dir = mxGetPr(prhs[5]);
    evn_stop = mxGetPr(prhs[6]);
    
    //extract values for event_id, event_dir, event_stop
    int *event_id = new int[evn_dim];
    int *event_dir = new int[evn_dim];
    int *event_stop = new int[evn_dim];
   
    for (k=0; k<evn_dim; k++){
        event_id[k] = (int)evn_id[k];
        event_dir[k] = (int)evn_dir[k];
        event_stop[k] = (int)evn_stop[k];
    }
    
    if (dim==6){
        state = prop(a,time,mu,event_id,evn_state,event_dir,event_stop,dim,time_dim,evn_dim);
    }
    else{
        printf("State vector should have 6 columns");
    }

    int nrows = state.size()/(dim+2);
    int nelements = state.size();

    int nsteps=0, nevents=0;
    for (k=0; k<nrows; k++){
        if (state[(dim+2)*k]==0){
            nsteps++;
        }
        else {
            nevents++;
        }
	}
        
    //associate outputs
    c_out_m = plhs[0] = mxCreateDoubleMatrix(nsteps,dim,mxREAL);
    t_out_m = plhs[1] = mxCreateDoubleMatrix(nsteps,1,mxREAL);
    cev_out_m = plhs[2] = mxCreateDoubleMatrix(nevents,dim,mxREAL);
    tev_out_m = plhs[3] = mxCreateDoubleMatrix(nevents,1,mxREAL);
    iev_out_m = plhs[4] = mxCreateDoubleMatrix(nevents,1,mxREAL);

    //associate pointers
    c = mxGetPr(c_out_m);
    t = mxGetPr(t_out_m);
    cev = mxGetPr(cev_out_m);
    tev = mxGetPr(tev_out_m);
    iev = mxGetPr(iev_out_m);

    
    int k1 = 0, k2 = 0;
    for (k=0; k<nrows; k++){
        for (j=0; j<dim+2; j++){
            if (state[(dim+2)*k]==0){
                if (j>0 && j<dim+1){
                    c[(j-1)*nsteps+k1] = state[(dim+2)*k+j];
                }
                else if (j==dim+1){
                    t[k1] = state[(dim+2)*k+j];
                }
            }
            else {
                if (j==0){
                    iev[k2] = state[(dim+2)*k+j];
                }
                else if (j>0 && j<dim+1){
                    cev[(j-1)*nevents+k2] = state[(dim+2)*k+j];
                }
                else if (j==dim+1){
                    tev[k2] = state[(dim+2)*k+j];
                }
            }
        }
        if (state[(dim+2)*k]==0){
            k1++;
        }
        else {
            k2++;
        }
	}
    
    delete(event_id);
    delete(event_dir);
    delete(event_stop);
    
    return ;
}

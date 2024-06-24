#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cblas.h>
#include <openblas_config.h>
#include <immintrin.h>
#include <float.h>

#include "timer_linux.h"

void zero(double* x, int n){
    for (int i = 0; i < n; i++){
        x[i] = 0.00f;
    }
}
// Function to calculate the mean square value of a vector
double mean_square(const double* vec, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += vec[i] * vec[i];
    }
    return sum / n;
}
// Function to calculate the SNR between two vectors
double calculate_snr(const double* original, const double* estimated, int n) {
    double signal_power = mean_square(original, n);

    // Calculate noise
    double* noise = (double*)malloc(n * sizeof(double));
    if (noise == NULL) {
        fprintf(stderr, "Memory allocation failed for noise\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; i++) {
        noise[i] = original[i] - estimated[i];
    }
    
    double noise_power = mean_square(noise, n);
    free(noise);

    // Calculate SNR
    double snr = 10 * log10(signal_power / noise_power);
    return snr;
}
// Function to create a vector with n elements and k non-zero entries
void createVector(int n, int k, double *x) {
    // Initialize all elements to 0
    for (int i = 0; i < n; i++) {
        x[i] = 0.0f;
    }
    
    // Ensure we do not attempt to place more non-zero elements than the size of the vector
    if (k > n) {
        k = n;
    }

    // Randomly assign k elements to be non-zero
    int count = 0;
    while (count < k) {
        int index = rand() % n;
        if (x[index] == 0.0f) {
            x[index] = 3*(double)rand()/RAND_MAX; // Assign a non-zero value, here we use 1 for simplicity
            count++;
        }
    }
}



// GPSR_BB function definition
void GPSR_BB(double* b, double* A, int n, int m, double tau, double* x_init, double* x_out) {
    // Parameters
    int stopCriterion = 3;
    double tolA = 0.01;
    int maxiter = 10000;
    int miniter = 5;
    int init = 0;
    int enforceMonotone = 1;
    double alphamin = 1e-30;
    double alphamax = 1e30;
    int continuation = 0;
    int cont_steps = -1;

    // Initialization
    
    double* Atb = (double*)malloc(m * sizeof(double));
    double* x = (double*)malloc(m * sizeof(double));
    memcpy(x, x_init, m * sizeof(double));



    // Precompute A'*b
    cblas_dgemv(CblasRowMajor, CblasTrans, n, m, 1.0, A, m, b, 1, 0.0, Atb, 1);

    double* aux = Atb;
    double max_tau = 0;
    for ( int i = 0; i < m; i++){
        double __aux = fabs(aux[i]);
        if (max_tau<__aux) max_tau = __aux;
    }
    
    if (tau >= max_tau) {
        memset(x_out, 0, m * sizeof(double));
       
        free(Atb);
        free(x);
        return;
    }

    // Initialize u and v
    double* u = (double*)malloc(m * sizeof(double));
    double* v = (double*)malloc(m * sizeof(double));
    for (int i = 0; i < m; i++) {
        double __temp_x = x[i];
        u[i] = fmax(__temp_x, 0);
        v[i] = -fmin(__temp_x, 0);
    }

    // Define the indicator vector or matrix of nonzeros in x
    // int num_nz_x = 0;
    // for (int i = 0; i < m; i++) {
    //     if (x[i] != 0.0) num_nz_x++;
    // }

    // Store given tau for continuation procedure
    double final_tau = tau;
    int final_stopCriterion = stopCriterion;
    double final_tolA = tolA;

    // Set continuation factors
    double cont_factors = 1.0;
    if (!continuation) {
        cont_steps = 1;
    }

    int iter = 1;
    int keep_continuation = 1;
    int cont_loop = 1;

    double* resid      = (double*)malloc(n * sizeof(double));
    double* resid_base = (double*)malloc(n * sizeof(double));
    double* temp       = (double*)malloc(m * sizeof(double));
    double* term       = (double*)malloc(m * sizeof(double));
    double* gradu      = (double*)malloc(m * sizeof(double));
    double* gradv      = (double*)malloc(m * sizeof(double));
    double* du         = (double*)malloc(m * sizeof(double));
    double* dv         = (double*)malloc(m * sizeof(double));
    double* dx         = (double*)malloc(m * sizeof(double));
    double* old_u      = (double*)malloc(m * sizeof(double));
    double* old_v      = (double*)malloc(m * sizeof(double));
    double* auv        = (double*)malloc(n * sizeof(double));
    double* uvmin      = (double*)malloc(m * sizeof(double));
  
      
  
        

    // Loop for continuation
    while (keep_continuation) {
        
        cblas_dcopy(n, b, 1, resid, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, n, m, -1.0, A, m, x, 1, 1.0, resid, 1);

        tau = final_tau * cont_factors;
        if (cont_loop == cont_steps) {
            stopCriterion = final_stopCriterion;
            tolA = final_tolA;
            keep_continuation = 0;
        } else {
            stopCriterion = 1;
            tolA = 1e-5;
        }

        double alpha = 1.0;
        

        for (int i = 0; i < n; i++)
        {
            resid_base[i] = b[i] - resid[i];
        }
        
        // cblas_dcopy(n, resid, 1, resid_base, 1);

        int keep_going = 1;

        while (keep_going) {
            
            cblas_dgemv(CblasRowMajor, CblasTrans, n, m, 1.0, A, m, resid_base, 1, 0.0, temp, 1);

            double dd2 = 0;
            for (int i = 0; i < m; i++) {
                double __temp_term;
                double __temp_grad_u;
                double __temp_grad_v;
                double __temp_du;
                double __temp_dv;
                double __temp_u = u[i];
                double __temp_v = v[i];

                __temp_term = temp[i] - Atb[i];
            
                __temp_grad_u = __temp_term + tau;
                __temp_grad_v = -1*__temp_term + tau;
                

                __temp_du = fmax(__temp_u - alpha * __temp_grad_u, 0.0) - __temp_u;
                __temp_dv = fmax(__temp_v - alpha * __temp_grad_v, 0.0) - __temp_v;

                term[i]  = __temp_term;

                gradu[i] = __temp_grad_u;
                gradv[i] = __temp_grad_v;
           
                dx[i] = __temp_du - __temp_dv;

                du[i] = __temp_du;
                dv[i] = __temp_dv;

                old_u[i] = __temp_u;
                old_v[i] = __temp_v;

                dd2 += __temp_du*__temp_grad_u + __temp_dv*__temp_grad_v;
            }

            
            cblas_dgemv(CblasRowMajor, CblasNoTrans, n, m, 1.0, A, m, dx, 1, 0.0, auv, 1);
            double dGd = cblas_ddot(n, auv, 1, auv, 1);

            double lambda_val;
            if (enforceMonotone == 1) {
               
                // double lambda0 = - (cblas_ddot(m, gradu, 1, du, 1) + cblas_ddot(m, gradv, 1, dv, 1)) / (DBL_EPSILON + dGd);
                double lambda0 = -dd2 / (DBL_EPSILON + dGd);
                if (lambda0 < 0) {
                    printf("ERROR: lambda0 = %f negative. Quit\n", lambda0);
                    free(Atb);
                    free(x);
                    free(u);
                    free(v);
                    free(resid);
                    free(resid_base);
                    free(temp);
                    free(term);
                    free(gradu);
                    free(gradv);
                    free(du);
                    free(dv);
                    free(dx);
                    free(old_u);
                    free(old_v);
                    free(auv);
                    return;
                }
                lambda_val = fmin(lambda0, 1.0);
            } else {
                lambda_val = 1.0;
            }

            double dd = 0;
            for (int i = 0; i < m; i++) {
                double __temp_u;
                double __temp_v;
                double __temp_uvmin;
                double __temp_du = du[i];
                double __temp_dv = dv[i];
                // u[i] = old_u[i] + lambda_val * du[i];
                // v[i] = old_v[i] + lambda_val * dv[i];
                __temp_u = old_u[i] + lambda_val * __temp_du;
                __temp_v = old_v[i] + lambda_val * __temp_dv;
            
                __temp_uvmin = fmin(__temp_u, __temp_v);
            
                __temp_u -= __temp_uvmin;
                __temp_v -= __temp_uvmin;

                x[i] = __temp_u - __temp_v;
                u[i] = __temp_u;
                v[i] = __temp_v;

                dd += __temp_du*__temp_du + __temp_dv*__temp_dv;
            }

            for (int i = 0; i < n; i++) {
                resid[i] = b[i]-resid_base[i] - lambda_val*auv[i];
            }


            // double dd = cblas_ddot(m, du, 1, du, 1) + cblas_ddot(m, dv, 1, dv, 1);
            if (dGd <= 0) {
                printf("dGd=%f, nonpositive curvature detected\n", dGd);
                alpha = alphamax;
            } else {
                alpha = fmin(alphamax, fmax(alphamin, dd / dGd));
            }

            for (int i = 0; i < n; i++) {
                resid_base[i] +=  lambda_val*auv[i];
            }

            iter++;

            if (stopCriterion == 3) {
                double* w = (double*)malloc(2 * m * sizeof(double));
                for (int i = 0; i < m; i++) {
                    w[i] = fmin(gradu[i], old_u[i]);
                    w[i + m] = fmin(gradv[i], old_v[i]);
                }

                double criterionLCP = cblas_dasum(2 * m, w, 1);
                double den = fmax(1.0e-6, fmax(cblas_dasum(m, old_u, 1), cblas_dasum(m, old_v, 1)));
                criterionLCP /= den;
                keep_going = criterionLCP > tolA;

                free(w);
            }

            if (iter <= miniter) {
                keep_going = 1;
            } else if (iter > maxiter) {
                keep_going = 0;
            }

            
        }

        free(temp);
        free(term);
        free(gradu);
        free(gradv);
        free(du);
        free(dv);
        free(dx);
        free(old_u);
        free(old_v);
        free(auv);
        free(uvmin);
        free(resid);
        free(resid_base);
        cont_loop++;
    }

    memcpy(x_out, x, m * sizeof(double));

 
    free(Atb);
    free(x);
    free(u);
    free(v);
}



double mean(double* values, int n)
{
    int i;
    double s = 0;

    for ( i = 0; i < n; i++ )
        s += values[i];
    return s / n;
}

double stddev(double* values, int n)
{
    int i;
    double average = mean(values,n);
    double s = 0;

    for ( i = 0; i < n; i++ )
        s += (values[i] - average) * (values[i] - average);
    return sqrt(s / (n - 1));
}

void generate(int numberOfrows, int numberOfcolumns, double* A){
    int m = numberOfrows;
    for (int index_j = 0; index_j < numberOfcolumns; index_j++){
        for (int i = 0; i < m; i += 2 ){
            double x, y, rsq, f;
            do {
                x = 2.0 * rand() / (double)RAND_MAX - 1.0;
                y = 2.0 * rand() / (double)RAND_MAX - 1.0;
                rsq = x * x + y * y;
            } while ( rsq >= 1. || rsq == 0. );
            f = sqrt( -2.0 * log(rsq) / rsq );
            A[i*numberOfcolumns + index_j] = x * f;
            if(i+1 < m) {
                A[(i+1)*numberOfcolumns + index_j] = y * f;
            }
        }
    }
}

void normalize_columns(double* A, int n, int m) {
    for (int j = 0; j < m; ++j) {
        // Compute the norm of the j-th column
        double norm = cblas_dnrm2(n, &A[j], m);
        
        // Normalize the j-th column if norm is greater than 0
        if (norm > 0) {
            cblas_dscal(n, 1.0 / norm, &A[j], m);
        }
    }
}

int main(int argc, char* argv[]){
    srand((unsigned int)time(NULL));
    
    double start, med, end;

    // Variables
  
    openblas_set_num_threads(1);


    int      n      = atoi(argv[1]);
    int      m      = atoi(argv[2]);
    int      k      = atoi(argv[3]);
    double sigma    = 0.0;
    double tau     = 0.1;         
                                                                                                  
                                                                 //                                                                                 [ ]
                                                                 //                                                                                 [ ]
    double* A       = (double*)malloc(n * m * sizeof(double));   // A is a matrix with size of n*m   Sensing matrix            [                   ][ ] = [ ]
    double* b       = (double*)malloc(n * sizeof(double));       // b is a vector with size of n*1   Measurement Signal        [                   ][ ] = [ ]
    double* true_x  = (double*)malloc(m * sizeof(double));       // true_x is a vector with size of m*1  Original Sparse Signal[                   ][ ] = [ ]
    double* x       = (double*)malloc(m * sizeof(double));       // x is a vector with size of m*1  Estimated Sparse Signal                         [ ]
                                                                 //                                                                                 [ ] 
    double* noise  = (double*)malloc(n * sizeof(double));        // noise is a vector with size of n*1  Noise Signal                    A            x  = b


    

    // Create true sparse signal vector
    createVector(m, k, true_x);
    
    // Generate sensing matrix A
    generate(n, m, A);
    
    // Normalize columns of A
    normalize_columns(A, n, m);
    
    // Perform matrix-vector multiplication A * true_x to get b
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, m, 1.0, A, m, true_x, 1, 0.0, b, 1);
    
    // Add noise to vector b
    for (int i = 0; i < n; i++) {
        b[i] += sigma * ((double)rand() / RAND_MAX);
    }

    // Allocate memory for initial and output sparse signals
    double* x_init = (double*)malloc(m * sizeof(double));
    memset(x_init, 0,m);
    
    // Timing the execution of GPSR_BB algorithm
    GET_TIME(start);
    GPSR_BB(b, A, n, m, tau, x_init, x);
    GET_TIME(end);

    // Calculate SNR between true_x and x
    double snr = calculate_snr(true_x, x, m);
    printf("Serial execution time: %f seconds, SNR = %lf\n", end - start, snr);
    

    // Free the variables
    free(A);
    free(true_x);
    free(x);
    free(b);
    free(noise);
    return 0;
}

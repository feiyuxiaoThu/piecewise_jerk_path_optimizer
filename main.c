#include "qpSWIFT.h"
#include "Matrices.h"

#include <stdlib.h>
#include "matrix_gen_utils.h"

int main(int argc, char *argv[])
{

    qp_real C0 = ST_Result[0];
    qp_real C1 = ST_Result[1];
    qp_real C2 = ST_Result[2];
    qp_real C3 = ST_Result[3];
    qp_real C4 = ST_Result[4];
    qp_real C5 = ST_Result[5];

    qp_real seg_t = ST_Result[7];
    qp_real seg_s = ST_Result[6];

    qp_real T,T2,T3,T4,T5;
    T = seg_t;
    T2 = seg_t*seg_t;
    T3 = T2*seg_t;
    T4 = T3*seg_t;
    T5 = T4*seg_t;

    const qp_real s_ini = 0.0;
    const qp_real v_ini = C1; // velocity
    const qp_real s_end = seg_s;
    const qp_real v_end = C1 + 2*C2*T + 3*C3*T2 + 4*C4*T3;
    const qp_real a_end = 2*C2 + 6*C3*T + 12*C4*T2;

    printf("%lf:%lf\n", v_end, a_end);

    const qp_int nrows_A = 5;// delete a_end constraints
    const qp_int ncols_A = 6;

    enum MatFormat format = RowMajor;

    //Equality Constraint A Matrix
    qp_real A_Matrix[5][6] = {{0,0,0,0,0,0},{0,0,0,0,0,0},\
    {0,0,2,0,0,0},{1,T,T2,T3,T4,T5},\
    {0,1,2*T,3*T2,4*T3,5*T4}};

    qp_real* A = (qp_real*) malloc(nrows_A*ncols_A*sizeof(qp_real));

    for (int i = 0; i<nrows_A; i++){
        for (int j = 0; j<ncols_A; j++){
            A[j + i * ncols_A] = A_Matrix[i][j];
        }
    }
    print_dense_matrix(A, nrows_A, ncols_A, format, "sparsified dense matrix");
    // beq
    qp_real beq[5] = {0,v_ini,0,s_end,v_end}; // boundary condition at initial and seg point

    csc* csc_matrix_A  = dense_to_csc_matrix(A,3,6,RowMajor);
    print_csc_matrix(csc_matrix_A, "compressed sparse column matrix A");
    //printf("%lf:%lf\n", a_end, v_end);


    //Inequality Constraint A Matrix
    const qp_int num_t_sample = floor(seg_t/dt_sample); // 0.5 s sample t
    const qp_int nrows_G = num_t_sample;
    const qp_int ncols_G = 6;

    qp_real G_Matrix[24][6]; // max 8s 16 points

    for (int i = 0; i < nrows_G; i++){
        // upper bounder
        qp_real t_sam = i*dt_sample;
        G_Matrix[i][0] = 1;
        G_Matrix[i][1] = t_sam;
        G_Matrix[i][2] = t_sam*t_sam;
        G_Matrix[i][3] = G_Matrix[i][2] * t_sam;
        G_Matrix[i][4] = G_Matrix[i][3] * t_sam;
        G_Matrix[i][5] = G_Matrix[i][4] * t_sam;

        // lower bounder
        G_Matrix[i+num_t_sample][0] = -1;
        G_Matrix[i+num_t_sample][1] = -G_Matrix[i][1];
        G_Matrix[i+num_t_sample][2] = -G_Matrix[i][2];
        G_Matrix[i+num_t_sample][3] = -G_Matrix[i][3];
        G_Matrix[i+num_t_sample][4] = -G_Matrix[i][4];
        G_Matrix[i+num_t_sample][5] = -G_Matrix[i][5];

        // velocity always bigger than zero
        G_Matrix[i+2*num_t_sample][0] = 0;
        G_Matrix[i+2*num_t_sample][1] = -1;
        G_Matrix[i+2*num_t_sample][2] = -2*t_sam;
        G_Matrix[i+2*num_t_sample][3] = -3*t_sam*t_sam;
        G_Matrix[i+2*num_t_sample][4] = -4*t_sam*t_sam*t_sam;
        G_Matrix[i+2*num_t_sample][5] = -5*t_sam*t_sam*t_sam*t_sam;
    }


    qp_real* G = (qp_real*) malloc(3*nrows_G*ncols_G*sizeof(qp_real));

    for (int i = 0; i<3*nrows_G; i++){
        for (int j = 0; j<ncols_G; j++){
            G[j + i * ncols_G] = G_Matrix[i][j];
        }
    }


    csc* csc_matrix_G  = dense_to_csc_matrix(G,nrows_G,ncols_G,RowMajor);




    print_csc_matrix(csc_matrix_G, "compressed sparse column matrix");


    // b
    qp_real* b_ineq = (qp_real*) malloc(3*nrows_G*sizeof(qp_real));

    for (int i = 0 ; i<nrows_G ; i++){
        b_ineq[i] = vehicles_front[0]*i*dt_sample + vehicles_front[i];
        if(veh_num == 2)
        {
            b_ineq[i+num_t_sample] = vehicle_behind[0]*i*dt_sample + vehicle_behind[i];
        }
        else{
            b_ineq[i+num_t_sample] = 0.0;
            }
        b_ineq[i+2*num_t_sample] = 0.0;
    }

    qp_real Q_Matrix[6][6];
    const qp_int nrows_Q = 6;
    const qp_int ncols_Q = 6;

    for (int i = 0; i<nrows_Q; i++){
        for (int j = 0; j<ncols_Q; j++){
            Q_Matrix[i][j] = 0;
        }
    }

    for (int i = 0; i<nrows_Q; i++){
        Q_Matrix[i][i] = 1;
    }

    qp_real* Q = (qp_real*) malloc(nrows_Q*ncols_Q*sizeof(qp_real));
    for (int i = 0; i<nrows_Q; i++){
        for (int j = 0; j<ncols_Q; j++){
            Q[j + i * ncols_Q] = Q_Matrix[i][j];
        }
    }

    csc* csc_matrix_Q  = dense_to_csc_matrix(Q,nrows_Q,ncols_Q,RowMajor);



    print_csc_matrix(csc_matrix_Q, "compressed sparse column matrix");



    QP *myQP;
    //myQP = QP_SETUP(n, m, p, Pjc, Pir, Ppr, Ajc, Air, Apr, Gjc, Gir, Gpr, c, h, b, sigma_d, NULL);
    /*For only inequality constrained QP set the pointers of A matrix and b vectro to zero and p = 0 and appropraiatley sized Permut matrix*/
    /*myQP = QP_SETUP(n, m, 0 , Pjc, Pir, Ppr, NULL, NULL, NULL, Gjc, Gir, Gpr, c, h, NULL, sigma_d, NULL);  */
    myQP = QP_SETUP(6, 3*num_t_sample, 5, csc_matrix_Q->p, csc_matrix_Q->i, csc_matrix_Q->x, csc_matrix_A->p,\
    csc_matrix_A->i, csc_matrix_A->x,csc_matrix_G->p,csc_matrix_G->i,csc_matrix_G->x,c,b_ineq,beq,sigma_d,NULL);

    qp_int ExitCode = QP_SOLVE(myQP);

    if (myQP != NULL)
        printf("Setup Time     : %f ms\n", myQP->stats->tsetup * 1000.0);
    if (ExitCode == QP_OPTIMAL)
    {
        printf("Solve Time     : %f ms\n", (myQP->stats->tsolve + myQP->stats->tsetup) * 1000.0);
        printf("KKT_Solve Time : %f ms\n", myQP->stats->kkt_time * 1000.0);
        printf("LDL Time       : %f ms\n", myQP->stats->ldl_numeric * 1000.0);
        printf("Diff	       : %f ms\n", (myQP->stats->kkt_time - myQP->stats->ldl_numeric) * 1000.0);
        printf("Iterations     : %ld\n", myQP->stats->IterationCount);
        printf("Optimal Solution Found\n");
    }
    if (ExitCode == QP_MAXIT)
    {
        printf("Solve Time     : %f ms\n", myQP->stats->tsolve * 1000.0);
        printf("KKT_Solve Time : %f ms\n", myQP->stats->kkt_time * 1000.0);
        printf("LDL Time       : %f ms\n", myQP->stats->ldl_numeric * 1000.0);
        printf("Diff	       : %f ms\n", (myQP->stats->kkt_time - myQP->stats->ldl_numeric) * 1000.0);
        printf("Iterations     : %ld\n", myQP->stats->IterationCount);
        printf("Maximum Iterations reached\n");
    }

    if (ExitCode == QP_FATAL)
    {
        printf("Unknown Error Detected\n");
    }

    if (ExitCode == QP_KKTFAIL)
    {
        printf("LDL Factorization fail\n");
    }

    printf("Solution\n");

    for (int i = 0; i < n; ++i)
        printf("x[%d]:%lf\n", i, myQP->x[i]);

    QP_CLEANUP(myQP);

    free(csc_matrix_A->x);
    free(csc_matrix_A->i);
    free(csc_matrix_A->p);
    free(csc_matrix_A);

    free(csc_matrix_G->x);
    free(csc_matrix_G->i);
    free(csc_matrix_G->p);
    free(csc_matrix_G);

    free(csc_matrix_Q->x);
    free(csc_matrix_Q->i);
    free(csc_matrix_Q->p);
    free(csc_matrix_Q);

    free(A);
    free(G);
    free(b_ineq);
    free(Q);

    return 0;
}






/*! @file */
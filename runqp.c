#include "qpSWIFT.h"
//#include "Matrices.h"
#include "GlobalOptions.h"
#include <stdlib.h>
#include "matrix_gen_utils.h"

/* For reading csv files */
#include<stdio.h>
#include<stdlib.h>

int main(){
    /* Parameters settings */
    qp_real sigma_d = 0.0;
    
    /*input data */
    const qp_real ST_Results_Seg1[8] = {0, 23.8712, -0.2676, -0.4442, 0.1086, 0 , 46.8477, 2.1000 };
    const qp_real ST_Results_Seg2[8] = {1.8821, 22.2034, -0.4456, 0.0471, -0.0016, 0, 163.2100, 7.6000};
    const qp_real ST_Results_Seg3[8] = {0,0,0,0,0,0,0,0};
    const c_int veh1_e = 1; // 1 exit ; 0 not exist
    const c_float veh1_dis = -15.0; // distance
    const c_float veh1_vel = 18.0; // velocity
    const c_int veh2_e = 1;
    const c_float veh2_dis = 60.0;
    const c_float veh2_vel = 18.0;
    const c_int veh3_e = 1;
    const c_float veh3_dis = 40;
    const c_float veh3_vel = 22.0;

    /* qp_problem settings */
    qp_real t_tol = 7.0;
    qp_real v_final = 0.0;
    qp_int st_seg_num = 0;


    // get the ref from the DP results
    const qp_real* ST_seg1 = ST_Results_Seg1;
    const qp_real* ST_seg2 = ST_Results_Seg2;
    const qp_real* ST_seg3 = ST_Results_Seg3;

    qp_real v_ini = ST_seg1[1];

    if (ST_seg3[7] > 0.1){
        st_seg_num = 3;
        t_tol = (double) floor(ST_seg3[7]*10)/10.0;
        v_final = ST_seg3[1] + 2*ST_seg3[2]*t_tol + 3*ST_seg3[3]*pow(t_tol,2) + 4*ST_seg3[4]*pow(t_tol,3);
    }
    else if(ST_seg2[7] > 0.1){ 
        st_seg_num = 2;
        t_tol = (double) floor(ST_seg2[7]*10)/10.0;
        v_final = ST_seg2[1] + 2*ST_seg2[2]*t_tol + 3*ST_seg2[3]*pow(t_tol,2) + 4*ST_seg2[4]*pow(t_tol,3);
    }
    else{
        st_seg_num = 1;
        t_tol = (double) floor(ST_seg1[7]*10)/10.0;
        v_final = ST_seg1[1] + 2*ST_seg1[2]*t_tol + 3*ST_seg1[3]*pow(t_tol,2) + 4*ST_seg1[4]*pow(t_tol,3);
    }

    //Note the t_tol needs to be modified for real condition
    qp_real delta_t = t_tol/(NSTEP-1); 
    qp_real JERKLIMIT = JERKPARA * delta_t;
    // Note veh1 and veh2 must exist for Entering QP
    c_int veh1_exist = veh1_e; // Target lane behind vehicle
    c_int veh2_exist = veh2_e; // Targer lane front vehicle
    c_int veh3_exist = veh3_e; // ego lane front vehicle
    c_float veh1_s = veh1_dis;
    c_float veh2_s = veh2_dis;
    c_float veh3_s = veh3_dis;
    c_float veh1_v = veh1_vel;
    c_float veh2_v = veh2_vel;
    c_float veh3_v = veh3_vel;

    c_float veh_len = safe_zone;
    c_float ref_v = v_final; // ? 

    // S Bound
    c_float s_front[NSTEP];
    c_float s_behind[NSTEP];
    c_float s_egolanefront[NSTEP];

    //ref line
    c_float ref_s[NSTEP];

    //For Plot
    c_float ref_vel[NSTEP];
    c_float ref_a[NSTEP];

    if (st_seg_num == 1){
        for(int i = 0; i<NSTEP; i++){
            double t = i*delta_t;
            ref_s[i] = ST_seg1[0] + ST_seg1[1]*t + ST_seg1[2]*t*t + ST_seg1[3]*pow(t,3) + ST_seg1[4]*pow(t,4);
            
            //For Plot
            ref_vel[i] = ST_seg1[1] + 2*ST_seg1[2]*t + 3*ST_seg1[3]*pow(t,2) + 4*ST_seg1[4]*pow(t,3);
            ref_a[i] = 2*ST_seg1[2] + 6*ST_seg1[3]*t + 12*ST_seg1[4]*pow(t,2);
        }
    }
    else if (st_seg_num == 2){
        for(int i = 0; i<NSTEP; i++){
            double t = i*delta_t;
            if(t<=ST_seg1[7]){
                ref_s[i] = ST_seg1[0] + ST_seg1[1]*t + ST_seg1[2]*t*t + ST_seg1[3]*pow(t,3) + ST_seg1[4]*pow(t,4);
                //For Plot
                ref_vel[i] = ST_seg1[1] + 2*ST_seg1[2]*t + 3*ST_seg1[3]*pow(t,2) + 4*ST_seg1[4]*pow(t,3);
                ref_a[i] = 2*ST_seg1[2] + 6*ST_seg1[3]*t + 12*ST_seg1[4]*pow(t,2);
            }
            else{
                ref_s[i] = ST_seg2[0] + ST_seg2[1]*t + ST_seg2[2]*t*t + ST_seg2[3]*pow(t,3) + ST_seg2[4]*pow(t,4);

                //For Plot
                ref_vel[i] = ST_seg2[1] + 2*ST_seg2[2]*t + 3*ST_seg2[3]*pow(t,2) + 4*ST_seg2[4]*pow(t,3);
                ref_a[i] = 2*ST_seg2[2] + 6*ST_seg2[3]*t + 12*ST_seg2[4]*pow(t,2);
            }
        }
    }
    else{
        for(int i = 0; i<NSTEP; i++){
            double t = i*delta_t;
            if(t<=ST_seg1[7]){
                ref_s[i] = ST_seg1[0] + ST_seg1[1]*t + ST_seg1[2]*t*t + ST_seg1[3]*pow(t,3) + ST_seg1[4]*pow(t,4);

                //For Plot
                ref_vel[i] = ST_seg1[1] + 2*ST_seg1[2]*t + 3*ST_seg1[3]*pow(t,2) + 4*ST_seg1[4]*pow(t,3);
                ref_a[i] = 2*ST_seg1[2] + 6*ST_seg1[3]*t + 12*ST_seg1[4]*pow(t,2);
            }
            else if(t>=ST_seg2[7]){
                ref_s[i] = ST_seg3[0] + ST_seg3[1]*t + ST_seg3[2]*t*t + ST_seg3[3]*pow(t,3) + ST_seg3[4]*pow(t,4);

                //For Plot
                ref_vel[i] = ST_seg3[1] + 2*ST_seg3[2]*t + 3*ST_seg3[3]*pow(t,2) + 4*ST_seg3[4]*pow(t,3);
                ref_a[i] = 2*ST_seg3[2] + 6*ST_seg3[3]*t + 12*ST_seg3[4]*pow(t,2);
            }
            else{
                ref_s[i] = ST_seg2[0] + ST_seg2[1]*t + ST_seg2[2]*t*t + ST_seg2[3]*pow(t,3) + ST_seg2[4]*pow(t,4);

                //For Plot
                ref_vel[i] = ST_seg2[1] + 2*ST_seg2[2]*t + 3*ST_seg2[3]*pow(t,2) + 4*ST_seg2[4]*pow(t,3);
                ref_a[i] = 2*ST_seg2[2] + 6*ST_seg2[3]*t + 12*ST_seg2[4]*pow(t,2);
            }
        }
    }

    /*
    for(int i =0; i<NSTEP; i++){
        printf("%.2f,", ref_s[i]);
    }
    */ 

    // Now we get ref_s which is the DP result

    // The dynamic obstacles ST using constant velocity model
    for(int i = 0; i<NSTEP; i++){
        if(veh1_exist == 1){
            s_behind[i] = veh1_s + veh1_v*delta_t*i;
        }
        else{
            s_behind[i] = -10000.0; // not consider
        }

        if (veh2_exist == 1){
            s_front[i] = veh2_s + veh2_v *delta_t*i;
        }
        else{
            s_front[i] = 10000.0;
        }

        if(veh3_exist == 1){
            s_egolanefront[i] = veh3_s + veh3_v*delta_t*i;
        }
        else{
            s_egolanefront[i] = 10000.0;
        }
    }

    // Now we begin to construct the QP matrixs!

    /* Cost Function Q matrxi in dense format */
    c_float Q[NSTEP3][NSTEP3];
    for(int i = 0 ; i< 3*NSTEP; i++){
        for(int j=0; j<3*NSTEP; j++){
            Q[i][j] = 0.0;
        }
    }

    for(int i = 0; i< NSTEP; i++){
        Q[i][i] = 2*w1;
    }

    for(int i = NSTEP; i<2*NSTEP; i++){
        Q[i][i] = 2*w2;
    }

    for(int i= 2*NSTEP; i<3*NSTEP; i++){
        Q[i][i] = 2*w3;
    }

    enum MatFormat format = RowMajor;
    const c_int nrows1 = NSTEP3;
    const c_int ncols1 = NSTEP3;

    c_float* P = (c_float *)malloc(nrows1 * ncols1 * sizeof(c_float));
    for (int i_row = 0; i_row < nrows1; i_row++){
        for(int i_col = 0; i_col < ncols1; i_col++){
            P[arr_ind(i_col,i_row,nrows1,ncols1,format)] = Q[i_row][i_col];
        }
    }


    csc* csc_matrix_Q = dense_to_csc_matrix(P, NSTEP3, NSTEP3, RowMajor);
    /*
    c_int* Pjc = csc_matrix_Q->p;
    c_int* Pir = csc_matrix_Q->i;
    c_float* Ppr = csc_matrix_Q->x;

    free(Q);
    */

   free(P);

   /* Cost Function c vector */
   c_float c[NSTEP3];
   for(int i = 0; i< NSTEP; i++){
       c[i] = ref_s[i] * (-2) * w1;
   }

   for(int i = NSTEP; i<2*NSTEP; i++){
        c[i] = ref_v * (-2) * w2;
   }

   for (int i = 2*NSTEP; i< 3*NSTEP; i++){
    c[i] = 0.0;
   }

    /* Equality Constraint Aeq Matrix in dense format */
    // Equality matrix and constraints
    c_float Aeq[2*NSTEP + 1][NSTEP3];
    for(int i = 0; i < 2*NSTEP+1; i++){
        for (int j = 0; j<NSTEP3; j++){
            Aeq[i][j] = 0.0;
        }
    }

    // Part 1
    for(int i =0; i< NSTEP-1; i++){
        int j = i + NSTEP;
        // int j = NSTEP; j<2*NSTEP -1; j++
        Aeq[i][j] = -1.0;
        Aeq[i][j+1] = 1.0;
    }

    for(int i = 0; i< NSTEP -1; i++){
        int j = i+2*NSTEP;
        //int j = 2*NSTEP; j<3*NSTEP-1; j++
        Aeq[i][j] = -delta_t/2.0;
        Aeq[i][j+1] = -delta_t/2.0;
    }

    // Part 2
    for(int i = NSTEP -1; i < 2*NSTEP-2; i++){
        //int j = 0; j< NSTEP-1; j++
        int j = i-NSTEP + 1;
        Aeq[i][j] = -1.0;
        Aeq[i][j+1] = 1.0;
    }

    for(int i = NSTEP-1; i< 2*NSTEP-2; i++){
        //int j = NSTEP; j< 2*NSTEP-1; j++
        int j = i+1;
        Aeq[i][j] = -delta_t;
    }

    for(int i = NSTEP-1; i< 2*NSTEP-2; i++){
        //int j = 2*NSTEP; j<3*NSTEP-1; j++
        int j = i+NSTEP +1;
        Aeq[i][j] = -delta_t*delta_t/3.0;
        Aeq[i][j+1] = -delta_t*delta_t/6.0;
    }

    // Part 3
    Aeq[2*NSTEP - 2][0] = 1.0;
    Aeq[2*NSTEP - 1][NSTEP] = 1.0;
    Aeq[2*NSTEP][2*NSTEP] = 1.0;

    const c_int nrows2 = 2*NSTEP+1;
    const c_int ncols2 = NSTEP3;

    c_float* AA = (c_float *)malloc(nrows2*ncols2*sizeof(c_float));
    for(int i_row = 0; i_row < nrows2; i_row++){
        for(int i_col = 0; i_col < ncols2; i_col++){
            AA[arr_ind(i_col,i_row,nrows2,ncols2, format)] = Aeq[i_row][i_col];
        }
    }

    csc* csc_matrix_Aeq = dense_to_csc_matrix(AA, 2*NSTEP+1, NSTEP3, RowMajor);

    c_int* Ajc = csc_matrix_Aeq->p;
    c_int* Air = csc_matrix_Aeq->i;
    c_float* Apr = csc_matrix_Aeq->x;

    free(AA);

    c_float beq[2*NSTEP+1];

    for(int i =0; i<2*NSTEP+1; i=i+1){
        beq[i] = 0;
    }

    // Should fix the v_final and s_final ???
    beq[2*NSTEP-1] = v_ini;

    /* Inequality Constraint A matrix in dense format */
    c_float A[8*NSTEP-2][NSTEP3];

    for(int i=0; i<8*NSTEP-2; i++){
        for(int j =0; j<NSTEP3; j++){
            A[i][j] = 0.0;
        }
    }

    // A
    for(int i=0; i< NSTEP3; i++){
        A[i][i] = 1.0;
    }

    for(int i = NSTEP3; i< 4*NSTEP-1; i++){
        int j = i- NSTEP;
        A[i][j] = -1;
        A[i][j+1] = 1; 
    }

    //-A
    for(int i=4*NSTEP-1; i<7*NSTEP-1; i = i+1){
        int j = i-4*NSTEP + 1;
        A[i][j] = -1.0;
    }

    for(int i = 7*NSTEP-1; i<8*NSTEP-2; i++){
        int j = i-5*NSTEP+1;
        A[i][j] = 1.0;
        A[i][j+1] = -1.0;
    }

    const c_int nrows3 = 8*NSTEP-2;
    const c_int ncols3 = NSTEP3;
    
    c_float* G = (c_float *)malloc( ncols3 * nrows3 * sizeof(c_float)); //???
    for(c_int i_row = 0; i_row < nrows3; i_row++){
        for(c_int i_col = 0; i_col < ncols3; i_col++){
            c_int index = arr_ind(i_col,i_row,nrows3,ncols3,format);
            G[arr_ind(i_col,i_row,nrows3,ncols3,format)] = A[i_row][i_col];
        }
    }

   

    csc* csc_matrix_A = dense_to_csc_matrix(G, 8*NSTEP-2, NSTEP3, RowMajor);

    c_int* Gjc = csc_matrix_A->p;
    c_int* Gir = csc_matrix_A->i;
    c_float* Gpr = csc_matrix_A->x;

    free(G);

    c_float ul[8*NSTEP-2]; // lb < A < ub
    for(int i =0; i<8*NSTEP-2; i++){
        ul[i] = 0;
    }

    if(veh3_exist == 1){
        for(int i =0; i< NSTEP; i=i+1){
            ul[i] = MIN(s_front[i],s_egolanefront[i]);
            ul[i+4*NSTEP-1] = -s_behind[i] ; // -l_s
        }
    }
    else{ // only two vehicle in target lane
        for(int i=0; i< NSTEP; i++){
            ul[i] = s_front[i];
            ul[i+4*NSTEP-1] = -s_behind[i]; // -l_s
        }
    }

    for(int i = NSTEP;i<2*NSTEP; i++){
        ul[i] = MAXSPEED;
    }

    for(int i=2*NSTEP; i< 3*NSTEP; i++){
        ul[i] = ACCMAX;
        ul[i+4*NSTEP-1] = ACCMAX; //-l_s
    }

    for (int i = 3*NSTEP;i<4*NSTEP-1;i++){
        ul[i] =JERKLIMIT;
        ul[i+4*NSTEP-1] = JERKLIMIT;
    }




    /* Solving Process */
    QP *myQP;
   
    myQP = QP_SETUP(NSTEP3, 8*NSTEP - 2, 2*NSTEP + 1, csc_matrix_Q->p, csc_matrix_Q->i, csc_matrix_Q->x, csc_matrix_Aeq->p,\
    csc_matrix_Aeq->i, csc_matrix_Aeq->x,csc_matrix_A->p,csc_matrix_A->i,csc_matrix_A->x,c,ul,beq,sigma_d,NULL);



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

    //Write to file to visualize
    for (int i = 0; i < NSTEP; ++i){
        printf("x[%d]:%lf  ", i, myQP->x[i]);
        printf("v[%d]:%lf  ", i, myQP->x[i+NSTEP]);
        printf("a[%d]:%lf\n", i, myQP->x[i+2*NSTEP]);
    }

    FILE *QP_res = fopen("res.csv", "w+");
    if(QP_res == NULL){
        fprintf(stderr, "fopen() failed \n");
        exit(EXIT_FAILURE);
    }

    FILE *obs_info = fopen("obs.csv", "w+");
    if(obs_info == NULL){
        fprintf(stderr, "fopen() failed \n");
        exit(EXIT_FAILURE);
    }

    c_float jerk_original[NSTEP];
    c_float jerk_opt[NSTEP];

    jerk_original[0] = 0.0;
    jerk_opt[0] = 0.0;

    for(int i=1; i<NSTEP; i++){
        jerk_original[i] = (ref_a[i] - ref_a[i-1])/delta_t;
        jerk_opt[i] = (myQP->x[i+2*NSTEP-1] - myQP->x[i+2*NSTEP])/delta_t;
    }

    fprintf(QP_res, "time,original_dis,original_vel,original_acc,original_jerk,op_dis,op_vel,op_acc,op_jerk\n");
    for (int i = 0; i < NSTEP; ++i){
        fprintf(QP_res, "%f,%f,%f,%f,%f,%f,%f,%f,%f\n", i*delta_t, ref_s[i],ref_vel[i],ref_a[i],jerk_original[i] ,myQP->x[i], myQP->x[i+NSTEP], myQP->x[i+2*NSTEP],jerk_opt[i]);
    } 

    fclose(QP_res);

    /*
    c_float s_front[NSTEP];
    c_float s_behind[NSTEP];
    c_float s_egolanefront[NSTEP];
    */
    fprintf(obs_info, "time,s_front,s_behind,s_egolanefront\n");
    for (int i = 0; i < NSTEP; ++i){
        fprintf(obs_info, "%f,%f,%f,%f\n", i*delta_t, s_front[i],s_behind[i],s_egolanefront[i]);
    } 

    fclose(obs_info);


        

    QP_CLEANUP(myQP);

    free(csc_matrix_A->x);
    free(csc_matrix_A->i);
    free(csc_matrix_A->p);
    free(csc_matrix_A);

    free(csc_matrix_Aeq->x);
    free(csc_matrix_Aeq->i);
    free(csc_matrix_Aeq->p);
    free(csc_matrix_Aeq);

    free(csc_matrix_Q->x);
    free(csc_matrix_Q->i);
    free(csc_matrix_Q->p);
    free(csc_matrix_Q);


    return 0;
}
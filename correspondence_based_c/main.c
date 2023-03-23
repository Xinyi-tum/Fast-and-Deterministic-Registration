
#include "main.h"

int main(int argc,char* argv[])
{
    char* input_file=argv[1];
    char* parameter_file=argv[2];
    char* results_file=argv[3];

    clock_t start_tim,end_tim;//<---timer

    // input data

    unsigned int point_number;

    point_number=get_number(input_file);

    type_point_p* const point_p=(type_point_p*)calloc(point_number,sizeof(type_point_p));

    type_point_q* const point_q=(type_point_q*)calloc(point_number,sizeof(type_point_q));

    double v_p[3]={0};
    double v_q[3]={0};

    read_data(input_file,point_p,point_q,v_p,v_q);

    // parameter setting

    double epsilon=0.0175;

    type_tran_cube init_branch;

    read_setting(parameter_file,&epsilon,&init_branch);

    //registration

    double opt_tran[3]={0};
    double opt_R[9]={0};

    start_tim = clock();
    globally_search(point_number,point_p,point_q,v_p,v_q,&init_branch,epsilon,opt_tran );
    voting_R(point_number,point_p,point_q,v_p,v_q,opt_tran,opt_R);
    end_tim = clock();

    free(point_p);
    free(point_q);

    double tim=(double)(end_tim-start_tim)/(CLOCKS_PER_SEC/1000);

    write_result(results_file,opt_tran,opt_R,tim);

    return 0;
}

#ifndef BNBANDVOTING_H_INCLUDED
#define BNBANDVOTING_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define BIN_NUM 360

typedef struct
{
    double x;
    double y;
    double z;
} type_point_p;

typedef struct
{
    double x;
    double y;
    double z;
} type_point_q;

typedef struct
{
    double tx;
    double ty;
    double tz;
    double half_side;
} type_tran_cube;

typedef struct sub_branch//branch
{
    type_tran_cube branch_tran;
    unsigned int branch_lower;
    unsigned int branch_upper;
    struct sub_branch *previous;
    struct sub_branch *next;
} type_sub_branch;

unsigned int get_number(char *path_file);
void read_data(char *path_file,type_point_p *point_p_pointer,type_point_q *point_q_pointer,double *v_p,double *v_q);
void read_setting(char *path_file,double *epsilon,type_tran_cube *init_branch);
void branching(type_tran_cube given_branch,type_tran_cube *branch_tran_pointer);
void free_all_branches(type_sub_branch *branch_pointer);
type_sub_branch *add_sub_branch(type_sub_branch *branch_pool_start,type_sub_branch *new_branch_pointer);
type_sub_branch *delete_sub_branch(type_sub_branch *branch_pool_start,type_sub_branch *branch_pointer);
void get_bounds(unsigned int input_number,type_point_p *point_p_pointer,double *ang_gama_pointer,
                double *v_p,type_tran_cube *branch_pointer,double epsilon, unsigned int *bounds_pointer);
void globally_search(unsigned int input_number,type_point_p *point_p,type_point_q *point_q,double *v_p,double *v_q,
                     type_tran_cube *init_branch,double epsilon, double *opt_tran );
void voting_R(unsigned int input_number,type_point_p *point_p,type_point_q *point_q,double *v_p,double *v_q,double *opt_tran,double *opt_r);
void write_result(char *path_file,double *opt_t,double *opt_rot,double tim);

#endif // BNBANDVOTING_H_INCLUDED

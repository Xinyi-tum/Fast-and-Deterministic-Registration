#include "bnbandvoting.h"

unsigned int get_number(char *path_file)
{
    FILE *fp;
    unsigned int point_number=0;

    if((fp=fopen(path_file,"r"))==NULL)
    {
        printf("data file cannot be opened\n");
        exit(-1);
    }
    if(fscanf(fp,"%d",&point_number)==EOF)
    {
        printf("point number cannot be read\n");
        exit(-1);
    }
    fclose(fp);
    return point_number;
}

void read_data(char *path_file,type_point_p *point_p_pointer,type_point_q *point_q_pointer,double *v_p,double *v_q)
{
    FILE *fp;
    unsigned int point_number=0;

    if((fp=fopen(path_file,"r"))==NULL)
    {
        printf("data file cannot be opened\n");
        exit(-1);
    }

    if(fscanf(fp,"%d",&point_number)==EOF)
    {
        printf("input number cannot be read\n");
        exit(-1);
    }
    fscanf(fp,"%lf,%lf,%lf,%lf,%lf,%lf\n",&v_p[0],&v_p[1],&v_p[2],&v_q[0],&v_q[1],&v_q[2]);

    for(unsigned int i=0; i<point_number; i++)
    {
        fscanf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
               &(point_q_pointer[i].x),&(point_q_pointer[i].y),&(point_q_pointer[i].z),
               &(point_p_pointer[i].x),&(point_p_pointer[i].y),&(point_p_pointer[i].z));
    }
    fclose(fp);
    return;
}

void read_setting(char *path_file,double *epsilon,type_tran_cube *init_branch)
{
    FILE *fp;

    if((fp=fopen(path_file,"r"))==NULL)
    {
        printf("setting file cannot be opened\n");
        exit(-1);
    }
    fscanf(fp,"%lf,%lf,%lf,%lf,%lf\n",
           epsilon,&(init_branch->tx),&(init_branch->ty),&(init_branch->tz),&(init_branch->half_side));

    fclose(fp);
    return;
}

void branching(type_tran_cube given_branch,type_tran_cube *branch_tran_pointer)
{

    double new_half_side=0.5*given_branch.half_side;
    double t_x=given_branch.tx;
    double t_y=given_branch.ty;
    double t_z=given_branch.tz;

    for(unsigned int i=0; i<8; i++)
    {

        branch_tran_pointer[i].half_side=new_half_side;

        switch(i)
        {
            case 0:
                branch_tran_pointer[i].tx=t_x-new_half_side;
                branch_tran_pointer[i].ty=t_y-new_half_side;
                branch_tran_pointer[i].tz=t_z-new_half_side;
                break;
            case 1:
                branch_tran_pointer[i].tx=t_x-new_half_side;
                branch_tran_pointer[i].ty=t_y-new_half_side;
                branch_tran_pointer[i].tz=t_z+new_half_side;
                break;
            case 2:
                branch_tran_pointer[i].tx=t_x-new_half_side;
                branch_tran_pointer[i].ty=t_y+new_half_side;
                branch_tran_pointer[i].tz=t_z+new_half_side;
                break;
            case 3:
                branch_tran_pointer[i].tx=t_x-new_half_side;
                branch_tran_pointer[i].ty=t_y+new_half_side;
                branch_tran_pointer[i].tz=t_z-new_half_side;
                break;
            case 4:
                branch_tran_pointer[i].tx=t_x+new_half_side;
                branch_tran_pointer[i].ty=t_y-new_half_side;
                branch_tran_pointer[i].tz=t_z-new_half_side;
                break;
            case 5:
                branch_tran_pointer[i].tx=t_x+new_half_side;
                branch_tran_pointer[i].ty=t_y-new_half_side;
                branch_tran_pointer[i].tz=t_z+new_half_side;
                break;
            case 6:
                branch_tran_pointer[i].tx=t_x+new_half_side;
                branch_tran_pointer[i].ty=t_y+new_half_side;
                branch_tran_pointer[i].tz=t_z-new_half_side;
                break;
            case 7:
                branch_tran_pointer[i].tx=t_x+new_half_side;
                branch_tran_pointer[i].ty=t_y+new_half_side;
                branch_tran_pointer[i].tz=t_z+new_half_side;
                break;
            default:
                printf("something incorrect!");
                break;
        }
    }
    return ;
}

void free_all_branches(type_sub_branch *branch_pointer)
{
    type_sub_branch *pNode=branch_pointer;

    while(pNode!=NULL)
    {
        branch_pointer=branch_pointer->next;
        free(pNode);
        pNode=branch_pointer;
    }
    return;
}

type_sub_branch *add_sub_branch(type_sub_branch *branch_pool_start,type_sub_branch *new_branch_pointer)
{

    type_sub_branch* b_ptr=branch_pool_start;// add new branch
    if(branch_pool_start==NULL)
    {
        new_branch_pointer->previous=NULL;
        new_branch_pointer->next=NULL;
        return new_branch_pointer;
    }

    while(b_ptr!=NULL)
    {

        if((new_branch_pointer->branch_upper)<(b_ptr->branch_upper))
        {
            if(b_ptr->next==NULL)
            {
                b_ptr->next= new_branch_pointer;
                new_branch_pointer->previous=b_ptr;
                break;
            }
        }

        if((new_branch_pointer->branch_upper)>=b_ptr->branch_upper)
        {
            if(b_ptr->previous==NULL)
            {
                b_ptr->previous= new_branch_pointer;
                new_branch_pointer->next=b_ptr;
                branch_pool_start=new_branch_pointer;
                branch_pool_start->previous=NULL;
                break;
            }
            else
            {
                (b_ptr->previous)->next=new_branch_pointer;
                new_branch_pointer->previous=b_ptr->previous;
                new_branch_pointer->next=b_ptr;
                b_ptr->previous= new_branch_pointer;
                break;
            }
            break;
        }
        b_ptr=b_ptr->next;

    }
    return branch_pool_start;
}

type_sub_branch *delete_sub_branch(type_sub_branch *branch_pool_start,type_sub_branch *branch_pointer)
{
    if((branch_pointer->previous==NULL) &&(branch_pointer->next==NULL))
    {
        free(branch_pointer);
        branch_pointer=NULL;
        printf("It's empty!\n");
        getchar();
        return branch_pool_start;
    }
    if(branch_pointer->previous==NULL)
    {
        branch_pool_start=branch_pointer->next;
        (branch_pointer->next)->previous=NULL;
        free(branch_pointer);
        return branch_pool_start;
    }
    if(branch_pointer->next==NULL)//β
    {
        (branch_pointer->previous)->next=NULL;
        free(branch_pointer);
        return branch_pool_start;
    }
    (branch_pointer->previous)->next=branch_pointer->next;
    (branch_pointer->next)->previous=branch_pointer->previous;
    free(branch_pointer);
    return branch_pool_start;
}

void get_bounds(unsigned int point_number,type_point_p *point_p_pointer,double *ang_gama_pointer,
                double *v_p,type_tran_cube *branch_pointer,double epsilon, unsigned int *bounds_pointer)
{
    bounds_pointer[0]=0;
    bounds_pointer[1]=0;
    double t_c[3]= {branch_pointer->tx,branch_pointer->ty,branch_pointer->tz};
    double radius=sqrt(3)*(branch_pointer->half_side);
    double temp_dot=0;
    double tran_x=0;
    double tran_y=0;
    double tran_z=0;
    double length_3d=0;
    double ang_alpha=0;
    double diff=0;
    double delta=0;

    for(unsigned int i=0; i<point_number; i++)
    {
        tran_x=point_p_pointer[i].x+t_c[0];
        tran_y=point_p_pointer[i].y+t_c[1];
        tran_z=point_p_pointer[i].z+t_c[2];
        temp_dot=v_p[0]*tran_x+v_p[1]*tran_y+v_p[2]*tran_z;
        length_3d=sqrt(tran_x*tran_x+tran_y*tran_y+tran_z*tran_z);
        temp_dot=temp_dot/length_3d;
        ang_alpha=acos(temp_dot);
        diff=fabs(ang_alpha-ang_gama_pointer[i]);

        if(diff<=epsilon)
        {
            bounds_pointer[0]++;// lower bound +1
        }

        if(length_3d<=radius)
        {
            delta=M_PI; // case 2 in paper
        }
        else
        {
            delta=asin(radius/length_3d); // case 1 in paper
        }

        if(diff<=epsilon+delta)
        {
            bounds_pointer[1]++;// upper bound +1
        }
    }
    return;
}

void get_cross(double *a,double *b,double *c)
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
    return;
}

void make_unit(double *a)
{
    double length=0;
    length=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    a[0]=a[0]/length;
    a[1]=a[1]/length;
    a[2]=a[2]/length;
    return;
}

void make_exp_R(double *r,double angle,double *R)
{
    double sin_angle=sin(angle);
    double cos_angle=1.0-cos(angle);

    R[0]=1.0+cos_angle*(r[0]*r[0]-1.0);
    R[1]=-sin_angle*r[2]+cos_angle*r[0]*r[1];
    R[2]=sin_angle*r[1]+cos_angle*r[0]*r[2];
    R[3]=sin_angle*r[2]+cos_angle*r[0]*r[1];
    R[4]=1.0+cos_angle*(r[1]*r[1]-1.0);
    R[5]=-sin_angle*r[0]+cos_angle*r[1]*r[2];
    R[6]=-sin_angle*r[1]+cos_angle*r[0]*r[2];
    R[7]=sin_angle*r[0]+cos_angle*r[1]*r[2];
    R[8]=1.0+cos_angle*(r[2]*r[2]-1.0);
    return;
}

void multi_RR(double *a,double *b,double *c)
{
    c[0]=a[0]*b[0]+a[1]*b[3]+a[2]*b[6];
    c[1]=a[0]*b[1]+a[1]*b[4]+a[2]*b[7];
    c[2]=a[0]*b[2]+a[1]*b[5]+a[2]*b[8];
    c[3]=a[3]*b[0]+a[4]*b[3]+a[5]*b[6];
    c[4]=a[3]*b[1]+a[4]*b[4]+a[5]*b[7];
    c[5]=a[3]*b[2]+a[4]*b[5]+a[5]*b[8];
    c[6]=a[6]*b[0]+a[7]*b[3]+a[8]*b[6];
    c[7]=a[6]*b[1]+a[7]*b[4]+a[8]*b[7];
    c[8]=a[6]*b[2]+a[7]*b[5]+a[8]*b[8];
    return;
}

void globally_search(unsigned int point_number,type_point_p *point_p,type_point_q *point_q,double *v_p,double *v_q,
                     type_tran_cube *init_branch,double epsilon, double *opt_tran )
{
    double *const ang_gama =(double*)calloc(point_number,sizeof(double));
    double temp_dot=0;
    double length=0;

    for(unsigned int i=0; i<point_number; i++)
    {
        temp_dot=v_q[0]*point_q[i].x+v_q[1]*point_q[i].y+v_q[2]*point_q[i].z;
        length=sqrt(point_q[i].x*point_q[i].x+point_q[i].y*point_q[i].y+point_q[i].z*point_q[i].z);
        temp_dot=temp_dot/length;
        ang_gama[i]=acos(temp_dot);
    }

    type_sub_branch *branch_pool_start=(type_sub_branch*)malloc(sizeof(type_sub_branch));
    branch_pool_start->branch_tran=*init_branch;
    branch_pool_start->next=NULL;
    branch_pool_start->previous=NULL;
    branch_pool_start->branch_lower=0;
    branch_pool_start->branch_upper=point_number;
    type_sub_branch *best_branch=branch_pool_start;

    unsigned int global_L=0;
    unsigned int global_U=point_number;

    type_tran_cube new_translation[8]= {'\0'};
    unsigned int estimate_bounds[2]= {'\0'};

    type_sub_branch *opt_branch=best_branch;
    while(global_L<global_U)
    {
        branching(best_branch->branch_tran,new_translation);
        for(unsigned int i=0; i<8; i++)
        {
            get_bounds(point_number,point_p,ang_gama,v_p,&(new_translation[i]),epsilon, estimate_bounds);
            if(estimate_bounds[1]<global_L)
            {
                continue;
            }
            type_sub_branch *new_branch_pointer=(type_sub_branch*)malloc(sizeof(type_sub_branch));
            new_branch_pointer->previous=NULL;
            new_branch_pointer->next=NULL;
            new_branch_pointer->branch_upper=estimate_bounds[1];
            new_branch_pointer->branch_lower=estimate_bounds[0];
            new_branch_pointer->branch_tran=new_translation[i];
            if(estimate_bounds[0]>global_L)
            {
                opt_branch=new_branch_pointer;
                global_L=estimate_bounds[0];
            }
            branch_pool_start=add_sub_branch(branch_pool_start,new_branch_pointer);
        }
        branch_pool_start=delete_sub_branch(branch_pool_start,best_branch);
        best_branch=branch_pool_start;
        global_U=best_branch->branch_upper;
    }
    opt_tran[0]=opt_branch->branch_tran.tx;
    opt_tran[1]=opt_branch->branch_tran.ty;
    opt_tran[2]=opt_branch->branch_tran.tz;
    free(ang_gama);
    free_all_branches(branch_pool_start);
    return;
}

void voting_R(unsigned int point_number,type_point_p* point_p,type_point_q* point_q,
                       double *v_p,double*v_q,double*opt_tran,double *opt_R)
{
    double R_min[9]= {'\0'};
    double r_min[3]= {'\0'};
    double angle_min=0;

    make_unit(v_p);
    make_unit(v_q);
    angle_min=acos(v_p[0]*v_q[0]+v_p[1]*v_q[1]+v_p[2]*v_q[2]);
    get_cross(v_p,v_q,r_min);
    if(r_min[0]==0&&r_min[1]==0&&r_min[2]==0)
    {
        R_min[0]=1;
        R_min[1]=0;
        R_min[2]=0;
        R_min[3]=0;
        R_min[4]=1;
        R_min[5]=0;
        R_min[6]=0;
        R_min[7]=0;
        R_min[8]=1;
    }
    else
    {
        make_unit(r_min);
        make_exp_R(r_min,angle_min,R_min);
    }

    double point_3d_tran[3]= {'\0'};
    double vv_p[3]= {'\0'};
    double vv_1q[3]= {'\0'};
    double vv_2q[3]= {'\0'};
    double temp_dot=0;
    double temp_dot_q=0;
    double m_3d[3]= {'\0'};
    double n_3d[3]= {'\0'};
    unsigned int voting_pool[BIN_NUM]= {0};
    double single_bin=2*M_PI/BIN_NUM;
    double angle_extend=0;
    double r_extend[3]= {'\0'};
    double R_extend[9]= {'\0'};
    unsigned int index=0;
    double p_3d[3]= {0};
    double q_3d[3]= {0};

    for(unsigned int ii=0; ii<point_number; ii++)
    {
        p_3d[0]=point_p[ii].x;
        p_3d[1]=point_p[ii].y;
        p_3d[2]=point_p[ii].z;

        q_3d[0]=point_q[ii].x;
        q_3d[1]=point_q[ii].y;
        q_3d[2]=point_q[ii].z;

        point_3d_tran[0]=p_3d[0]+opt_tran[0];
        point_3d_tran[1]=p_3d[1]+opt_tran[1];
        point_3d_tran[2]=p_3d[2]+opt_tran[2];

        vv_p[0]=R_min[0]*point_3d_tran[0]+R_min[1]*point_3d_tran[1]+R_min[2]*point_3d_tran[2];
        vv_p[1]=R_min[3]*point_3d_tran[0]+R_min[4]*point_3d_tran[1]+R_min[5]*point_3d_tran[2];
        vv_p[2]=R_min[6]*point_3d_tran[0]+R_min[7]*point_3d_tran[1]+R_min[8]*point_3d_tran[2];

        temp_dot=v_q[0]*vv_p[0]+v_q[1]*vv_p[1]+v_q[2]*vv_p[2];
        temp_dot_q=v_q[0]*q_3d[0]+v_q[1]*q_3d[1]+v_q[2]*q_3d[2];

        vv_1q[0]=temp_dot*v_q[0];
        vv_1q[1]=temp_dot*v_q[1];
        vv_1q[2]=temp_dot*v_q[2];

        m_3d[0]=vv_p[0]-vv_1q[0];
        m_3d[1]=vv_p[1]-vv_1q[1];
        m_3d[2]=vv_p[2]-vv_1q[2];
        make_unit(m_3d);

        vv_2q[0]=temp_dot_q*v_q[0];
        vv_2q[1]=temp_dot_q*v_q[1];
        vv_2q[2]=temp_dot_q*v_q[2];

        n_3d[0]=q_3d[0]-vv_2q[0];
        n_3d[1]=q_3d[1]-vv_2q[1];
        n_3d[2]=q_3d[2]-vv_2q[2];
        make_unit(n_3d);

        angle_extend=acos(m_3d[0]*n_3d[0]+m_3d[1]*n_3d[1]+m_3d[2]*n_3d[2]);
        get_cross(vv_p,q_3d,r_extend);

        if(r_extend[0]*v_q[0]+r_extend[1]*v_q[1]+r_extend[2]*v_q[2]<0)
        {
            angle_extend=2*M_PI-angle_extend;
        }
        index=floor(angle_extend/(single_bin));//<------------counting bin
        voting_pool[index]++;
    }
    unsigned int max_count_angle=0;
    index=0;
    max_count_angle=voting_pool[0];
    for(unsigned int ii=0; ii<BIN_NUM; ii++)
    {
        if(max_count_angle<voting_pool[ii])
        {
            max_count_angle=voting_pool[ii];
            index=ii;
        }
    }
    angle_extend=single_bin*(index+0.5);
    make_exp_R(v_q,angle_extend,R_extend);
    multi_RR(R_extend,R_min,opt_R);
    return;
}

void write_result(char* path_file,double*opt_t,double* opt_rot,double tim)
{
    FILE *fp;

    if((fp=fopen(path_file,"w"))==NULL)
    {
        printf("data file cannot be opened\n");
        exit(-1);
    }
    fprintf(fp,"%.10lf\t%.10lf\t%.10lf\t",opt_rot[0],opt_rot[1],opt_rot[2]);
    fprintf(fp,"%.10lf\t%.10lf\t%.10lf\t",opt_rot[3],opt_rot[4],opt_rot[5]);
    fprintf(fp,"%.10lf\t%.10lf\t%.10lf\t",opt_rot[6],opt_rot[7],opt_rot[8]);
    fprintf(fp,"%.10lf\t%.10lf\t%.10lf\t",opt_t[0],opt_t[1],opt_t[2]);
    fprintf(fp,"%lf\t",tim);
    fclose(fp);
    return ;
}

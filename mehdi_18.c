#include "udf.h"
#include "math.h"
#include "sg_udms.h"
#include "sg.h"
#include "stdio.h"
#include "mem.h"
#include "dpm.h"
#include "surf.h"

double Total_Lenght;
double xx, xx_shadow;
double heat_flux, htc, temp_infinity;
double HFG, zigma;
double M_m, PI_num;
double R_univ, R_R;
double T_ref, P_ref;
double P_op1,P_op2;
double M_v1, M_v2;
double M_l1, M_l2;
double Porousity, K_eff, wick_volume, wick_liquid_density_0;
double CP_L,CP_S,RO_L,RO_S;
double Y_wick1,Y_wick2,Y_wick3,Y_wick4;
double K_wick1,K_wick2,K_wick3;
double D_wick1,D_wick2,D_wick3;
double Porousity1,Porousity2,Porousity3;
double Visc_Re_1,Visc_Re_2,Visc_Re_3;
double Inter_Re_1,Inter_Re_2,Inter_Re_3;
double temp_init;

double aa1,aa2,aa3,aa4=0.0;
double bb1,bb2,bb3,bb4=0.0;
double x[ND_ND];
double A[ND_ND],es[ND_ND],A_by_es,dr0[ND_ND], ds;
double As[ND_ND],ess[ND_ND],A_by_ess,dr0s[ND_ND], ds_shadow;
double temp_cell, temp_cell_shadow;
double temp_face, temp_face_shadow;
double k_cell, k_cell_shadow;
double ro_cell, ro_cell_shadow;
double cp_cell, cp_cell_shadow;
double P_cell, P_cell_shadow;
double v_face, v_face_shadow;
double mdot,mdot_balance, P_int,P_vapor_0;
double u_r_mass,u_r_m_p,u_r_m_p1,u_r_m_p2,u_r_temp;
double vel_diff,vel_max,vel_max_p;
double vel_max_V_vapor,vel_max_U_vapor,vel_max_U_wick,x_seperation,x_seperation1,x_seperation2,x_vel_max_U_vapor,x_vel_max_U_wick,T_wall_max,T_wall_min;
double Q_vapor,Q_vapor_e,Q_vapor_c,Q_wick,W_wall,Q_out,Q_in;
double Ru_vapor,Ru_vapor_e,Ru_vapor_c,Ru_vapor_max,Ru_vapor_min, Ru_wick;
double P_vapor_max,P_vapor_min,P_wick_max,P_wick_min;
double mdot_error,velocity_error, tmp_error;
double k_ds; k_ds_shadow; 
double alfa,beta, alfa_max;
double m_dot[1000],V_inter_v[1000],V_inter_w[1000],T_int_w_v[1000];

int iii_wick_vapor[1000],iii_vapor_wick[1000];
int zone_ID, N_iteration,N_time,N_print,N_print_time,N_print_iter,iii,jjj,i,j;

int	marz_wick_vapor_ID=189;
int	marz_vapor_wick_ID=228;
int	marz_wall_wick_ID=215;
int	marz_wick_wall_ID=213;

int	wall_heating_ID=220;
int	wall_insulated_ID=219;
int	wall_cooling_ID=218;

int	vapor_core_ID=8;
int	wick_core_ID=14;
int	wall_a_ID=56;

Domain *dd;

face_t ff;
face_t ff_shadow;
Thread *tt;
Thread *t0;
Thread *t0_shadow;
Thread *tt_shadow;
cell_t c0;
cell_t c0_shadow;

face_t ff_wall_wick [1000];
face_t ff_wick_wall [1000];
face_t ff_wick_vapor [1000];
face_t ff_vapor_wick [1000];

DEFINE_INIT(init_parameter, domain)
{

	FILE *fp0;
	FILE *fp1;
	FILE *fp2;
	FILE *fp3;
	fp0 = fopen ("Data_0Transient.txt", "w");
	fp1 = fopen ("Data_1Interface.txt", "w"); 
	fp2 = fopen ("Data_2Wall.txt", "w"); 
	fp3= fopen ("Data_3UnderRelax.txt", "w");
	fprintf (fp0,"TIME            P_op2            T_wall_max      T_wall_min      Vap_max_vel     Q_out           Q_in            Q_vapor          Q_vapor_e       Q_vapor_c        mdot_balance     M_v2           M_l2            vel_max_inte    Wick_max_vel    x_seperation    x_max_U_vapor   x_max_U_wick    P_drop_Vapor    P_drop_Wick       Ru_vapor        Ru_vapor_e      Ru_vapor_c      Ru_vapor_max    Ru_vapor_min\n"); 
	fclose (fp0);
	fclose (fp1);
	fclose (fp2);
	fclose (fp3);

	//reading solid thermal properties from a wall-face
	zone_ID = marz_wall_wick_ID;
	tt = Lookup_Thread(domain,zone_ID);
	begin_f_loop (ff,tt)   
	{
		c0 = F_C0(ff,tt); 
		t0= F_C0_THREAD(ff,tt); 
		RO_S=C_R_M1(c0,t0);    ///???????
		CP_S=C_CP(c0,t0);
	}
	end_f_loop (ff,tt)
		bb4=bb4;
	//reading liquid thermal properties from a wick-face
	zone_ID = marz_wick_wall_ID;
	tt = Lookup_Thread(domain,zone_ID);
	begin_f_loop (ff,tt)   
	{
		c0 = F_C0(ff,tt); 
		t0= F_C0_THREAD(ff,tt); 
		RO_L=C_R_M1(c0,t0);   ////????
		CP_L=C_CP(c0,t0);
	}
	end_f_loop (ff,tt)
		bb4=bb4;

	Total_Lenght=0.370;

//heat_flux=17.42084554;   //Q=30
//heat_flux=34.39316262;   //Q=60
heat_flux=51.48508728;   //Q=90
//heat_flux=68.35573791;   //Q=120
//heat_flux=84.03031263;   //Q=150

//htc=985.21;
//htc=1153.97;
htc=1163.47;
//htc=1154.64;
//htc=1217.54;


	temp_infinity=21+273.15;
	temp_init=temp_infinity;
	HFG=2406*1E3;    //@ T=40 C
	zigma=0.03;
	M_m=18.015;
	PI_num=3.141592653589;
	R_univ=8314.40;
	R_R=R_univ/M_m; 

	wick_liquid_density_0=992.45;

	Y_wick1=(5.55)*0.001;
	Y_wick2=(5.55-0.28)*0.001;
	Y_wick3=(5.55-0.280-0.2)*0.001;
	Y_wick4=(5.55-0.280-0.2-0.1)*0.001;  //for non-ideal case
	K_wick1=1.72;
	K_wick2=62.507;
	K_wick3=62.507;     //for non-ideal case
	D_wick1=K_wick1/CP_L;
	D_wick2=K_wick2/CP_L;
	D_wick3=K_wick3/CP_L;
	Porousity1=0.713;
	Porousity2=0.707;
	Porousity3=0.707; //for non-ideal case
	Visc_Re_1=3.331E+08;
	Visc_Re_2=1.161E+10;
	Visc_Re_3=1.161E+10; //for non-ideal case
	Inter_Re_1=8.663E+03;
	Inter_Re_2=5.178E+04;
	Inter_Re_3=5.178E+04; //for non-ideal case

	P_op1=2490; 
	
	P_ref=P_op1;
	T_ref=temp_init;
	P_ref=P_op1;
	u_r_m_p2=-1.0; ///?????
	u_r_m_p1=-1.0;
	u_r_m_p=u_r_m_p2;
	u_r_mass=pow(10.0,u_r_m_p);
	u_r_temp=8.0E-1;
	vel_diff=-10.0;
	vel_max=0.0;
	vel_max_p=30.0;
	N_iteration=0;
	N_time=1;
	N_print=1;
	N_print_time=1;
	N_print_iter=100;
	P_vapor_0=0.0;
	/// tempreature initial condition for all the domains
	thread_loop_c (t0,domain)
	{
		begin_c_loop_all (c0,t0)   
		{
			C_UDSI(c0,t0,0)=temp_init;			
		}
		end_c_loop_all (c0,t0)
	}
	/// mass of vapor calcaulation
	zone_ID = vapor_core_ID;
	t0 = Lookup_Thread(domain,zone_ID);
	aa1=0.0;
	iii=0;
	begin_c_loop (c0,t0)   
	{
		iii=iii+1;
		C_CENTROID(x,c0,t0);
		aa1=aa1+C_VOLUME(c0,t0)/C_UDSI(c0,t0,0); 
	}
	end_c_loop (c0,t0)
		bb4=bb4;
	M_v1=aa1*P_op1/R_R;
	/// mass of liquid calcaulation

	zone_ID = wick_core_ID;
	t0 = Lookup_Thread(domain,zone_ID);
	wick_volume=0.0;
	begin_c_loop (c0,t0)   
	{
		C_CENTROID(x,c0,t0);

		if ((x[1] <= Y_wick1) && (x[1] >= Y_wick2))
		{
			wick_volume=wick_volume+C_VOLUME(c0,t0)*Porousity1; 
		}
		else if ((x[1] < Y_wick2) && (x[1] >= Y_wick3))
		{
			wick_volume=wick_volume+C_VOLUME(c0,t0)*Porousity2; 
		}
		else 
		{
			wick_volume=wick_volume+C_VOLUME(c0,t0)*Porousity3; 
		}

	}
	end_c_loop (c0,t0)
		bb4=bb4;
	M_l1=wick_liquid_density_0*wick_volume;

	P_op2=P_op1;
	M_v2=M_v1;
	M_l2=M_l1;

	iii=0;
	while (iii<1000)
	{
		m_dot[iii]=0;
		T_int_w_v[iii]=temp_init;
		V_inter_v[iii]=0.0;
		V_inter_w[iii]=0.0;
		iii=iii+1;
	}

	Message("P_op=%e   M_v=%e   M_l=%e \n",P_op1,M_v1,M_l1);
}
DEFINE_INIT(shadow_reading, domain)
{

	zone_ID = marz_wall_wick_ID;
	tt= Lookup_Thread(domain,zone_ID);
	zone_ID = marz_wick_wall_ID;
	tt_shadow = Lookup_Thread(domain,zone_ID);
	iii=-1;
	begin_f_loop(ff, tt)
	{
		iii=iii+1;
		F_CENTROID(x,ff,tt);
		xx = x[0];
		jjj=-1;
		begin_f_loop(ff_shadow, tt_shadow)
		{
			jjj=jjj+1;
			F_CENTROID(x,ff_shadow, tt_shadow);
			xx_shadow = x[0];
			if (fabs((xx_shadow-xx)/Total_Lenght) < 1.0E-6)
			{
				ff_wall_wick [iii]=ff_shadow;
			}
		}
		end_f_loop(ff_shadow, tt_shadow)
	}
	end_f_loop(ff, tt)
		bb4=bb4;
	//===================================
	zone_ID = marz_wick_wall_ID;
	tt= Lookup_Thread(domain,zone_ID);
	zone_ID = marz_wall_wick_ID;
	tt_shadow = Lookup_Thread(domain,zone_ID);
	iii=-1;
	begin_f_loop(ff, tt)
	{
		iii=iii+1;
		F_CENTROID(x,ff,tt);
		xx = x[0];
		begin_f_loop(ff_shadow, tt_shadow)
		{
			F_CENTROID(x,ff_shadow, tt_shadow);
			xx_shadow = x[0];

			if (fabs((xx_shadow-xx)/Total_Lenght) < 1.0E-6)
			{
				ff_wick_wall [iii]=ff_shadow;
			}
		}
		end_f_loop(ff_shadow, tt_shadow)
	}
	end_f_loop(ff, tt)
		bb4=bb4;
	//===================================		
	//=================================== Wick
	zone_ID = marz_wick_vapor_ID;
	tt= Lookup_Thread(domain,zone_ID);
	zone_ID = marz_vapor_wick_ID;
	tt_shadow = Lookup_Thread(domain,zone_ID);
	iii=-1;
	begin_f_loop(ff, tt)
	{
		iii=iii+1;
		F_CENTROID(x,ff,tt);
		xx = x[0];
		jjj=-1;
		begin_f_loop(ff_shadow, tt_shadow)
		{
			jjj=jjj+1;
			F_CENTROID(x,ff_shadow, tt_shadow);
			xx_shadow = x[0];
			if (fabs((xx_shadow-xx)/Total_Lenght) < 1.0E-6)
			{
				ff_wick_vapor [iii]=ff_shadow;
				iii_wick_vapor[jjj]=iii;
			}
		}
		end_f_loop(ff_shadow, tt_shadow)
	}
	end_f_loop(ff, tt)
		bb4=bb4;
	//===================================  Vapor
	zone_ID = marz_vapor_wick_ID;
	tt= Lookup_Thread(domain,zone_ID);
	zone_ID = marz_wick_vapor_ID;
	tt_shadow = Lookup_Thread(domain,zone_ID);
	iii=-1;
	begin_f_loop(ff, tt)
	{
		iii=iii+1;
		F_CENTROID(x,ff,tt);
		xx = x[0];
		jjj=-1;
		begin_f_loop(ff_shadow, tt_shadow)
		{
			jjj=jjj+1;
			F_CENTROID(x,ff_shadow, tt_shadow);
			xx_shadow = x[0];

			if (fabs((xx_shadow-xx)/Total_Lenght) < 1.0E-6)
			{
				ff_vapor_wick [iii]=ff_shadow;
				iii_vapor_wick[jjj]=iii;
			}
		}
		end_f_loop(ff_shadow, tt_shadow)
	}
	end_f_loop(ff, tt)
		bb4=bb4;
}
DEFINE_ADJUST(parameter_update, domain)
{

	FILE *fp;
	double P_op2_old,P_op2_new;
	double AAA,BBB,CCC,DDD;
	double NV_VEC(f_area);
	double d_area;
	double d_t;  ///time step

	d_t=CURRENT_TIMESTEP; /// ???
	//d_t=RP_Get_Real("physical-time-step");
	N_iteration=N_iteration+1;

	//====================================================
	/// mass balance at the interface vapor side
	zone_ID = marz_vapor_wick_ID;
	tt = Lookup_Thread(domain,zone_ID);
	aa1=0.0;
	iii=-1;
	begin_f_loop (ff,tt)   
	{
		iii=iii+1;
		c0 = F_C0(ff,tt); 
		t0= F_C0_THREAD(ff,tt); 
		v_face =F_V(ff,tt); 
		ro_cell=C_R(c0,t0);   

		F_AREA(f_area,ff,tt);
		d_area = NV_MAG(f_area);

		aa1=aa1+v_face*ro_cell*d_area; 
		//aa1=aa1+m_dot[iii]*d_area; 

	}
	end_f_loop (ff,tt)
		bb4=bb4;
	mdot_balance=aa1;
	M_v2=M_v1-d_t*(aa1);  //  ++?--?????????????
	M_l2=M_l1+d_t*(aa1);  // ++?--?????????????

	////???????????or or or or It could be calcuated using the vapor domain, just like initial calculation, could it be
	//====================================================
	//====================================================
	zone_ID = vapor_core_ID;
	t0 = Lookup_Thread(domain,zone_ID);
	aa1=0.0;
	begin_c_loop (c0,t0)   
	{
		aa1=aa1+C_VOLUME(c0,t0)/C_UDSI(c0,t0,0); 
	}
	end_c_loop (c0,t0)
		bb4=bb4;
	P_op2_old=M_v2/(aa1/R_R);
	BBB=R_R/aa1;

	zone_ID = marz_vapor_wick_ID;
	tt = Lookup_Thread(domain,zone_ID);
	CCC=0.0;
	DDD=0.0;
	iii=-1;
	begin_f_loop (ff,tt)   
	{
		iii=iii+1;
		c0 = F_C0(ff,tt);
		t0= F_C0_THREAD(ff,tt);
		v_face = F_V(ff,tt);
		ro_cell=C_R(c0,t0);
		temp_cell=C_UDSI(c0,t0,0);
		temp_face=T_int_w_v [iii];
		//temp_face=F_UDSI(c0,t0,0);
		P_int=P_ref*exp(HFG/R_R*(1.0/T_ref-1.0/temp_face));
		P_cell=C_P(c0,t0)-P_vapor_0;
		F_AREA(f_area,ff,tt);
		d_area = NV_MAG(f_area);

		CCC=CCC+d_area*(P_cell/sqrt(temp_cell)-P_int/sqrt(temp_face));
		DDD=DDD+d_area*(1/sqrt(temp_cell)); 
	}
	end_f_loop (ff,tt)
		bb4=bb4;
	AAA=(2.0*zigma/(2.0-zigma))*(1/sqrt(2*PI_num*R_R));
	P_op2_new=(BBB*(M_v1-d_t*AAA*CCC))/(1+BBB*d_t*AAA*DDD);
	//====================================================
	if ( fabs(P_op2_old-P_op2_new) > 1.01)
	{
		u_r_m_p=u_r_m_p2;
		u_r_mass=pow(10.0,u_r_m_p);
	}
	//P_op2=P_op2+u_r_mass*(P_op2_new-P_op2);
	P_op2=P_op2_new;
	//====================================================
	//====================================================
	//if ((floor(N_TIME/N_print_time)*N_print_time == N_TIME) && (floor(N_iteration/N_print_iter)*N_print_iter == N_iteration))
	if (floor(N_iteration/N_print_iter)*N_print_iter == N_iteration)
	{
		fp = fopen ("Data_3UnderRelax.txt", "a"); 
		fprintf (fp, "%5d    %E   %E   %E   %E   %E   %E   %E\n", N_iteration, vel_max, vel_diff, u_r_mass, fabs(P_op2_old-P_op2_new),mdot_error,velocity_error,tmp_error); 
		fclose (fp);
	}

	//====================================================
	//====================================================
	if (floor(N_iteration/5)*5 == N_iteration)
	{
		u_r_m_p=u_r_m_p+0.05*(u_r_m_p1-u_r_m_p);
		u_r_mass=pow(10.0,u_r_m_p);
		vel_max_p=vel_max;
	}
	vel_diff=0.0;
	vel_max=0.0;
	//====================================================
	//====================================================
}
DEFINE_ADJUST(boundray_condition, domain)
{

	double mdot_error2,velocity_error2,tmp_error2;

	mdot_error2=mdot_error;
	velocity_error2=velocity_error;
	tmp_error2=tmp_error;

	mdot_error=0.0;
	velocity_error=0.0;
	tmp_error=0.0;
	zone_ID = marz_vapor_wick_ID;  //Vapor
	tt = Lookup_Thread(domain,zone_ID);	
	zone_ID = marz_wick_vapor_ID;  //Wick
	tt_shadow = Lookup_Thread(domain,zone_ID);
	alfa_max=0.0;
	iii=-1;
	begin_f_loop(ff, tt)
	{
		iii=iii+1;
		F_CENTROID(x,ff,tt);
		xx = x[0];
		ff_shadow=ff_vapor_wick [iii];

		c0 = F_C0(ff,tt); 
		t0 = F_C0_THREAD(ff,tt);
		temp_cell=C_UDSI(c0,t0,0);
		temp_face=F_UDSI(ff,tt,0);
		k_cell=C_K_L(c0,t0); 
		ro_cell=C_R(c0,t0);   
		cp_cell=C_CP(c0,t0);   
		v_face=F_V(ff,tt);
		P_cell=C_P(c0,t0)-P_vapor_0;
		BOUNDARY_FACE_GEOMETRY(ff,tt,A,ds,es,A_by_es,dr0);

		c0_shadow = F_C0(ff_shadow,tt_shadow); 
		t0_shadow = F_C0_THREAD(ff_shadow,tt_shadow);
		temp_cell_shadow=C_UDSI(c0_shadow,t0_shadow,0);
		temp_face_shadow=F_UDSI(ff_shadow,tt_shadow,0);
		k_cell_shadow=C_K_L(c0_shadow,t0_shadow);
		ro_cell_shadow=C_R(c0_shadow,t0_shadow);   
		cp_cell_shadow=C_CP(c0_shadow,t0_shadow);   
		v_face_shadow=F_V(ff_shadow,tt_shadow);
		BOUNDARY_FACE_GEOMETRY(ff_shadow, tt_shadow,A,ds_shadow,ess,A_by_ess,dr0s);

		////======== mass transfer
		temp_face=T_int_w_v [iii];
		aa1=HFG/R_R*(1.0/T_ref-1.0/temp_face);
		P_int=P_ref*exp(aa1);
		aa1=2.0*zigma/(2.0-zigma);
		aa2=1/sqrt(2*PI_num*R_R);
		aa3=((P_op2+P_cell)/sqrt(temp_cell)-P_int/sqrt(temp_face)); 
		mdot=aa1*aa2*aa3;  
		bb1=mdot/ro_cell;
		velocity_error=velocity_error+fabs(bb1-V_inter_v[iii]);//?????
		V_inter_v[iii] = V_inter_v[iii]+u_r_mass*(bb1-V_inter_v[iii]);   
		V_inter_w[iii]= V_inter_v[iii]*ro_cell/ro_cell_shadow;
		mdot_error=mdot_error+fabs(mdot-m_dot[iii]);  ///????
		m_dot[iii]=m_dot[iii]+(u_r_mass/10)*(mdot-m_dot[iii]);

		if (vel_diff < fabs(fabs(bb1)-fabs(v_face)) )
		{
			vel_diff=fabs(fabs(bb1)-fabs(v_face));
		}
		if (fabs(vel_max) < fabs(bb1) )
		{
			vel_max=fabs(bb1);
		}
				////======== heat transfer
		
		k_ds=k_cell/ds;
		k_ds_shadow=k_cell_shadow/ds_shadow;

		aa1=temp_cell_shadow*k_ds_shadow;
		aa2=temp_cell*k_ds;
		aa3=2.0*zigma/(2.0-zigma)*1/sqrt(2*PI_num*R_R)*((P_op2+P_cell)/sqrt(temp_cell))*HFG;
		bb1=k_ds_shadow;
		bb2=k_ds;
		bb3=2.0*zigma/(2.0-zigma)*1/sqrt(2*PI_num*R_R)*(P_int/sqrt(temp_face))*HFG/temp_face;

		tmp_error=tmp_error+fabs(((aa1+aa2+aa3)/(bb1+bb2+bb3)-T_int_w_v [iii]));///??????
		T_int_w_v [iii] =T_int_w_v [iii]+u_r_mass*((aa1+aa2+aa3)/(bb1+bb2+bb3)-T_int_w_v [iii]);
	}
	end_f_loop(ff, tt)
		bb4=bb4;
	/*
	if (mdot_error2*1.15 < mdot_error)
	if ( mdot_error > 5.0E-5)
	{
	u_r_m_p=u_r_m_p2;
	u_r_mass=pow(10.0,u_r_m_p);
	}
	if (velocity_error > 1.0)
	{
	u_r_m_p=u_r_m_p2;
	u_r_mass=pow(10.0,u_r_m_p);
	}*/

	if (fabs(tmp_error2-tmp_error)/tmp_error > 0.2  )

	{
		u_r_m_p=u_r_m_p2;
		u_r_mass=pow(10.0,u_r_m_p);
	}

}
DEFINE_EXECUTE_AT_END(execute_at_end_4)
{
	FILE *fp1; 
	FILE *fp2; 
	//====================================================	
	if (floor(N_TIME/N_print)*N_print == N_TIME)
	{
		fp1 = fopen ("Data_0Transient.txt", "a"); 
		fprintf (fp1, "%E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E   %E\n", CURRENT_TIME, P_op2,T_wall_max,T_wall_min,vel_max_U_vapor,Q_out,Q_in,Q_vapor,Q_vapor_e,Q_vapor_c,mdot_balance,M_v2,M_l2,vel_max_V_vapor,vel_max_U_wick,x_seperation,x_vel_max_U_vapor,x_vel_max_U_wick,P_vapor_max-P_vapor_min,P_wick_max-P_wick_min,Ru_vapor,Ru_vapor_e,Ru_vapor_c,Ru_vapor_max,Ru_vapor_min);
		fclose (fp1);
	}
	//====================================================
	fp2 = fopen ("Data_3UnderRelax.txt", "a"); 
	fprintf (fp2, "TIME=%E \n",CURRENT_TIME); 
	fclose (fp2);
	//====================================================
	Ru_wick=M_l2/(wick_volume);
	P_op1=P_op2;
	M_v1=M_v2;
	M_l1=M_l2;
	N_iteration=0;
	N_time=N_time+1;
	u_r_m_p=u_r_m_p2;
	u_r_mass=pow(10.0,u_r_m_p);
}
DEFINE_EXECUTE_AT_END(interface_printout_1)
{
	FILE *fp2; 
	double data_interface [1000][10];
	double NV_VEC(f_area);
	double d_area;
	double ru;
	//====================================================
	//x_seperation on the wall of vapor-wick
	//Pressure Zero Point
	// Velocities at the Interface
	for (iii=0;iii<1000;iii=iii+1)
	{
		for (jjj=0;jjj<10;jjj=jjj+1)
		{
			data_interface [iii][jjj]=12345.0;
		}
	}

	dd=Get_Domain(1);
	zone_ID = marz_vapor_wick_ID;   /// marz_vapor_wick_zone_ID
	tt = Lookup_Thread(dd,zone_ID);
	zone_ID = marz_wick_vapor_ID;
	tt_shadow = Lookup_Thread(dd,zone_ID);
	P_vapor_max=-1.0E15;
	P_vapor_min=1.0E15;
	P_wick_max=-1.0E15;
	P_wick_min=1.0E15;
	iii=-1;
	begin_f_loop (ff,tt)   
	{
		iii=iii+1;
		c0 = F_C0(ff,tt); 
		t0 = F_C0_THREAD(ff,tt); 
		ff_shadow=ff_vapor_wick [iii];
		c0_shadow= F_C0(ff_shadow,tt_shadow); 
		t0_shadow= F_C0_THREAD(ff_shadow,tt_shadow); 

		F_CENTROID(x,ff,tt);
		aa1=C_P(c0,t0);

		temp_cell= C_UDSI(c0, t0, 0);
		ru=P_op2/R_R/temp_cell;

		F_AREA(f_area,ff,tt);
		d_area = NV_MAG(f_area);

		data_interface [iii][0]=x[0];
		data_interface [iii][1]=F_V(ff,tt);
		data_interface [iii][2]=F_V(ff_shadow,tt_shadow);
		data_interface [iii][3]=aa1;
		data_interface [iii][4]=m_dot[iii]*d_area;
		data_interface [iii][5]=ru;
		aa1=C_P(c0,t0);
		bb1=C_P(c0_shadow,t0_shadow);
		if (aa1>P_vapor_max)
		{
			P_vapor_max=aa1;
		}
		if (aa1<P_vapor_min)
		{
			P_vapor_min=aa1;
		}
		if (bb1>P_wick_max)
		{
			P_wick_max=bb1;
		}
		if (bb1<P_wick_min)
		{
			P_wick_min=bb1;
		}
	}
	end_f_loop (ff,tt)
		bb4=bb4;
	for (i=0;i<1000;i=i+1)
	{
		for (iii=0; iii<1000-1;iii=iii+1)
		{
			if (data_interface [iii][0]>data_interface [iii+1][0])
			{
				for (jjj=0;jjj<10;jjj=jjj+1)
				{
					bb4=data_interface [iii+1][jjj];
					data_interface [iii+1][jjj]=data_interface [iii][jjj];
					data_interface [iii][jjj]=bb4;
				}
			}
		}
	}
	vel_max_V_vapor=-1.0;	///?????
	x_seperation=0.0;
	P_vapor_0=0.0;
	Q_vapor_e=0.0;
	Q_vapor_c=0.0;
	Ru_vapor_e=0.0;
	Ru_vapor_c=0.0;
	i=0;
	j=0;
	for (iii=0;iii<1000;iii=iii+1)
	{
		if (data_interface [iii][0]<99.0)
		{
			if (vel_max_V_vapor < fabs(data_interface [iii][1]))
			{
				vel_max_V_vapor=fabs(data_interface [iii][1]);
			}
			if (data_interface [iii][4]<0)
			{
				i=i+1;
				Q_vapor_e=Q_vapor_e+data_interface [iii][4];
				Ru_vapor_e=Ru_vapor_e+data_interface [iii][5];
			}
			if (data_interface [iii][4]>0)
			{
				j=j+1;
				Q_vapor_c= Q_vapor_c+data_interface [iii][4];
				Ru_vapor_c=Ru_vapor_c+data_interface [iii][5];
			}
		}

		if (data_interface [iii][1]*data_interface [iii+1][1]< 0.0)
		{
			x_seperation1=data_interface [iii][0];   //x1
			aa2=data_interface [iii][1];   //v1
			aa3=data_interface [iii][3];   //P1
			x_seperation2=data_interface [iii+1][0];   //x2
			bb2=data_interface [iii+1][1];   //v2
			bb3=data_interface [iii+1][3];   //P2

			bb4=(bb2-aa2)/(x_seperation2-x_seperation1);  //(v2-v1)/(x2-x1)
			x_seperation=1/bb4*(0.0-aa2)+x_seperation1; // y-y1=a(x-x1)  ==>> (y-y1)*(1/a)+x1=x
			bb4=(bb3-aa3)/(x_seperation2-x_seperation1); //(p2-p1)/(x2-x1)
			P_vapor_0=bb4*(x_seperation-x_seperation1)+aa3;  // y-y1=a(x-x1)  ==>> y=a(x-x1)+y1

		}

	}
	Q_vapor_e=Q_vapor_e*HFG;
	Q_vapor_c=-Q_vapor_c*HFG;
	Ru_vapor_e=Ru_vapor_e/i;
	Ru_vapor_c=-Ru_vapor_c/j;
	//====================================================
	//====================================================
	if (floor(N_TIME/N_print_time)*N_print_time == N_TIME) 
	{
		fp2 = fopen ("Data_1Interface.txt", "a"); 
		fprintf (fp2, "Time= %E \n", CURRENT_TIME); 
		for (iii=0;iii<1000;iii=iii+1)
		{
			if (data_interface [iii][0]<99.0)
			{
				fprintf (fp2, "%E   %E   %E   %E   %E   %E\n",data_interface [iii][0],data_interface [iii][1],data_interface [iii][2],data_interface [iii][3]-P_vapor_0,data_interface [iii][4],data_interface [iii][5]); 
			}
		}
		fclose (fp2);
	}
}
DEFINE_EXECUTE_AT_END(vapor_wick_domain_2)
{

	double data_interface [1000][10];
	double NV_VEC(f_area);
	double d_area;
	double ru;
	//====================================================
	//Vapor max velocity in x direction and its position at the vapor domain
	vel_max_U_vapor=-1.0;  ///?????? sign
	Ru_vapor=0.0;
	Ru_vapor_max=-1.0E15;
	Ru_vapor_min=1.0E15;

	zone_ID = vapor_core_ID;
	t0 = Lookup_Thread(dd,zone_ID);
	iii=-1;
	begin_c_loop (c0,t0)   
	{
		iii=iii+1;
		if (vel_max_U_vapor < fabs(C_U(c0,t0)))
		{
			vel_max_U_vapor=fabs(C_U(c0,t0));
			C_CENTROID(x,c0,t0);
			x_vel_max_U_vapor=x[0];
		}

		temp_cell= C_UDSI(c0, t0, 0);
		ru=P_op2/R_R/temp_cell;
		if (Ru_vapor_max < ru )
		{
			Ru_vapor_max = ru;
		}
		if (Ru_vapor_min > ru )
		{
			Ru_vapor_min = ru;
		}
		Ru_vapor=Ru_vapor+ru;
	}
	end_c_loop (c0,t0)
		bb4=bb4;
	Ru_vapor=Ru_vapor/(iii+1);

	//====================================================
	//Wick max velocity in x direction and its position at the wick domain

	vel_max_U_wick=-1.0;
	zone_ID = wick_core_ID;
	t0 = Lookup_Thread(dd,zone_ID);
	begin_c_loop (c0,t0)   
	{
		if (vel_max_U_wick < fabs(C_U(c0,t0)))
		{
			vel_max_U_wick=fabs(C_U(c0,t0));
			C_CENTROID(x,c0,t0);
			x_vel_max_U_wick=x[0];
		}
	}
	end_c_loop (c0,t0)
		bb4=bb4;
	//====================================================
	//====================================================  Q going through the vapor core through mass flux at seperation point
	Q_vapor=0.0;
	zone_ID = vapor_core_ID;
	t0 = Lookup_Thread(dd,zone_ID);

	for (iii=0;iii<1000;iii=iii+1)
	{
		for (jjj=0;jjj<10;jjj=jjj+1)
		{
			data_interface [iii][jjj]=12345.0;
		}
	}

	iii=-1;
	begin_c_loop (c0,t0)   
	{
		c_face_loop(c0, t0, i)   
		{
			ff = C_FACE(c0,t0,i);
			tt = C_FACE_THREAD(c0,t0,i);
			F_CENTROID(x,ff,tt);
			aa1=x[0];
			aa2=x[1];
			C_CENTROID(x,c0,t0);
			bb1=x[0];
			bb2=x[1];
			if ((aa1 < x_seperation2) && (aa1 > x_seperation1) && (fabs(aa2-bb2)<1.0E-7))
			{
				iii=iii+1;
				F_AREA(f_area,ff,tt);
				d_area = NV_MAG(f_area);
				data_interface [iii][0]=bb2;
				data_interface [iii][1]=bb1;
				data_interface [iii][2]=d_area;
				data_interface [iii][3]=F_U(c0,t0);
				data_interface [iii][4]=C_R(c0,t0);
			}
		}
	}
	end_c_loop (c0,t0)
		bb4=bb4;
	for (i=0;i<1000;i=i+1)
	{
		for (iii=0; iii<1000-1;iii=iii+1)
		{
			if (data_interface [iii][0]>data_interface [iii+1][0])
			{
				for (jjj=0;jjj<10;jjj=jjj+1)
				{
					bb4=data_interface [iii+1][jjj];
					data_interface [iii+1][jjj]=data_interface [iii][jjj];
					data_interface [iii][jjj]=bb4;
				}
			}
		}
	}


	for (iii=0;iii<1000;iii=iii+2)
	{
		if ((data_interface [iii][0]< 99.0) && (data_interface [iii][1]< 99.0))
		{
			aa1=(data_interface [iii+1][3]-data_interface [iii][3])/(data_interface [iii+1][1]-data_interface [iii][1]);  //(U2-U1)/(x2-x1)
			bb1=aa1*(x_seperation-data_interface [iii][1])+data_interface [iii][3];// y-y1=a(x-x1)  ==>> y=a(x-x1)+y1
			data_interface [iii][5]=bb1;
			data_interface [iii][6]=bb1*data_interface [iii][2]*data_interface [iii][4];			
			Q_vapor=Q_vapor+data_interface [iii][6];
		}
	}
	Q_vapor=Q_vapor*HFG;
}
DEFINE_EXECUTE_AT_END(wall_printout_3)
{

	FILE *fp3;
	double data_wall [1000][10];
	double NV_VEC(f_area);
	double d_area;
	double ru;
	//====================================================
	T_wall_max=-10.0E10;
	T_wall_min=10.0E10;
	Q_in=0.0;
	Q_out=0.0;
	zone_ID = wall_heating_ID; //heating_wall_zone_ID
	tt = Lookup_Thread(dd,zone_ID);
	begin_f_loop (ff,tt)   
	{
		c0 = F_C0(ff,tt); 
		t0= F_C0_THREAD(ff,tt); 
		temp_cell=C_UDSI(c0,t0,0);
		temp_face=F_UDSI(ff,tt,0);
		k_cell=C_K_L(c0,t0);   
		BOUNDARY_FACE_GEOMETRY(ff,tt,A,ds,es,A_by_es,dr0);
		F_AREA(f_area,ff,tt);
		d_area = NV_MAG(f_area);
		Q_in=Q_in+(k_cell*(temp_cell-temp_face)/ds)*d_area;
		if (T_wall_max < F_UDSI(ff,tt,0))
		{
			T_wall_max=F_UDSI(ff,tt,0);
		}
	}
	end_f_loop (ff,tt)
		bb4=bb4;
	zone_ID = wall_cooling_ID; //cooling_wall_zone_ID
	tt = Lookup_Thread(dd,zone_ID);
	begin_f_loop (ff,tt)   
	{
		c0 = F_C0(ff,tt); 
		t0= F_C0_THREAD(ff,tt); 
		temp_cell=C_UDSI(c0,t0,0);
		temp_face=F_UDSI(ff,tt,0);
		k_cell=C_K_L(c0,t0);   
		BOUNDARY_FACE_GEOMETRY(ff,tt,A,ds,es,A_by_es,dr0);
		F_AREA(f_area,ff,tt);
		d_area = NV_MAG(f_area);
		Q_out=Q_out+(k_cell*(temp_cell-temp_face)/ds)*d_area;
		if (T_wall_min > F_UDSI(ff,tt,0))
		{
			T_wall_min=F_UDSI(ff,tt,0);
		}
	}
	end_f_loop (ff,tt)
		bb4=bb4;
	//====================================================
	for (iii=0;iii<1000;iii=iii+1)
	{
		for (jjj=0;jjj<10;jjj=jjj+1)
		{
			data_wall [iii][jjj]=12345.0;
		}
	}
	//=====
	zone_ID = wall_heating_ID;
	tt = Lookup_Thread(dd,zone_ID);
	iii=-1;
	begin_f_loop (ff,tt)   
	{
		iii=iii+1;
		F_CENTROID(x,ff,tt);
		xx=x[0];
		data_wall [iii][0]=xx;data_wall [iii][1]=F_UDSI(ff,tt,0);
	}
	end_f_loop (ff,tt)
		bb4=bb4;
	zone_ID = wall_insulated_ID;
	tt = Lookup_Thread(dd,zone_ID);
	begin_f_loop (ff,tt)   
	{
		iii=iii+1;
		F_CENTROID(x,ff,tt);
		xx=x[0];
		data_wall [iii][0]=xx;data_wall [iii][1]=F_UDSI(ff,tt,0);
	}
	end_f_loop (ff,tt)
		bb4=bb4;
	zone_ID = wall_cooling_ID;
	tt = Lookup_Thread(dd,zone_ID);
	begin_f_loop (ff,tt)   
	{
		iii=iii+1;
		F_CENTROID(x,ff,tt);
		xx=x[0];
		data_wall [iii][0]=xx;data_wall [iii][1]=F_UDSI(ff,tt,0);
	}
	end_f_loop (ff,tt)
		bb4=bb4;
	for (i=0;i<1000;i=i+1)
	{
		for (iii=0; iii<1000-1;iii=iii+1)
		{
			if (data_wall [iii][0]>data_wall [iii+1][0])
			{
				for (jjj=0;jjj<2;jjj=jjj+1)
				{
					bb4=data_wall [iii+1][jjj];
					data_wall [iii+1][jjj]=data_wall [iii][jjj];
					data_wall [iii][jjj]=bb4;
				}
			}
		}
	}

	//======================= print
	if (floor(N_TIME/N_print_time)*N_print_time == N_TIME)
	{
		fp3 = fopen ("Data_2Wall.txt", "a"); 
		fprintf (fp3, "Time= %E \n", CURRENT_TIME); 
		for (iii=0;iii<1000;iii=iii+1)
		{
			if (data_wall [iii][0]<99.0)
			{
				fprintf (fp3, "%E   %E\n",data_wall [iii][0],data_wall [iii][1]); 
			}
		}
		fclose (fp3);
	}
}
DEFINE_PROFILE(velo_vapor_wick, tt, ii) 
{
	iii=-1;
	begin_f_loop(ff, tt)
	{
		iii=iii+1;
		F_PROFILE(ff, tt, ii) = -V_inter_v[iii];   
	}
	end_f_loop(ff, tt)
}  //=============  END OF THE UDF
DEFINE_PROFILE(velo_wick_vapor, tt, ii) 
{
	iii=-1;
	begin_f_loop(ff, tt)
	{
		iii=iii+1;
		F_PROFILE(ff, tt, ii) = V_inter_w[iii_wick_vapor[iii]];
	}
	end_f_loop(ff, tt)
}
DEFINE_PROFILE(tmp_wick_vapor, tt, ii) 
{
	iii=-1;
	begin_f_loop(ff, tt)
	{
		iii=iii+1;
		F_PROFILE(ff, tt, ii) = T_int_w_v [iii_vapor_wick[iii]];
	}
	end_f_loop(ff, tt)
} //=============  END OF THE UDF
DEFINE_PROFILE(tmp_vapor_wick, tt, ii) 
{
	iii=-1;
	begin_f_loop(ff, tt)
	{
		iii=iii+1;
		F_PROFILE(ff, tt, ii) =T_int_w_v [iii];
	}
	end_f_loop(ff, tt)

} //=============  END OF THE UDF
DEFINE_PROFILE(tmp_wall_wick, tt, ii) 
{
	dd=Get_Domain(1);
	zone_ID = marz_wick_wall_ID;
	tt_shadow = Lookup_Thread(dd,zone_ID);
	iii=-1;

	begin_f_loop(ff, tt)
	{
		iii=iii+1;

		c0 = F_C0(ff,tt); 
		t0 = F_C0_THREAD(ff,tt);
		temp_cell=F_UDSI(c0,t0,0);
		temp_face=F_UDSI(ff,tt,0);
		k_cell=C_K_L(c0,t0);   

		BOUNDARY_FACE_GEOMETRY(ff,tt,A,ds,es,A_by_es,dr0);

		ff_shadow=ff_wall_wick [iii];

		c0_shadow = F_C0(ff_shadow,tt_shadow); 
		t0_shadow = F_C0_THREAD(ff_shadow,tt_shadow);
		temp_cell_shadow=C_UDSI(c0_shadow,t0_shadow,0);
		temp_face_shadow=F_UDSI(ff_shadow,tt_shadow,0);
		k_cell_shadow=C_K_L(c0_shadow,t0_shadow);

		BOUNDARY_FACE_GEOMETRY(ff_shadow, tt_shadow,A,ds_shadow,ess,A_by_ess,dr0s);

		k_ds=k_cell/ds;
		k_ds_shadow=k_cell_shadow/ds_shadow;
		aa1 = (k_ds*temp_cell+k_ds_shadow*temp_cell_shadow)/(k_ds+k_ds_shadow);
		F_PROFILE(ff, tt, ii) =aa1;
	}
	end_f_loop(ff, tt)
		bb4=bb4;
}
DEFINE_PROFILE(tmp_wick_wall, tt, ii) 
{

	dd=Get_Domain(1);
	zone_ID = marz_wall_wick_ID;
	tt_shadow = Lookup_Thread(dd,zone_ID);
	iii=-1;

	begin_f_loop(ff, tt)
	{
		iii=iii+1;

		c0 = F_C0(ff,tt); 
		t0 = F_C0_THREAD(ff,tt);
		temp_cell=C_UDSI(c0,t0,0);
		temp_face=F_UDSI(ff,tt,0);
		k_cell=C_K_L(c0,t0);   //??????

		BOUNDARY_FACE_GEOMETRY(ff,tt,A,ds,es,A_by_es,dr0);

		ff_shadow=ff_wick_wall [iii];

		c0_shadow = F_C0(ff_shadow,tt_shadow); 
		t0_shadow = F_C0_THREAD(ff_shadow,tt_shadow);
		temp_cell_shadow=C_UDSI(c0_shadow,t0_shadow,0);
		temp_face_shadow=F_UDSI(ff_shadow,tt_shadow,0);
		k_cell_shadow=C_K_L(c0_shadow,t0_shadow);   //???

		BOUNDARY_FACE_GEOMETRY(ff_shadow, tt_shadow,A,ds_shadow,ess,A_by_ess,dr0s);

		k_ds=k_cell/ds;
		k_ds_shadow=k_cell_shadow/ds_shadow;
		aa1 = (k_ds*temp_cell+k_ds_shadow*temp_cell_shadow)/(k_ds+k_ds_shadow);
		F_PROFILE(ff, tt, ii) =aa1;
	}
	end_f_loop(ff, tt)
		bb4=bb4;
} 
DEFINE_PROFILE(tmp_wall_cooling, tt, ii) 
{
	iii=-1;
	begin_f_loop(ff, tt)
	{
		iii=iii+1;
		c0 = F_C0(ff,tt); 
		t0 = F_C0_THREAD(ff,tt);
		temp_cell=C_UDSI(c0,t0,0);
		temp_face=F_UDSI(ff,tt,0);
		k_cell=C_K_L(c0,t0);   //?????????
		BOUNDARY_FACE_GEOMETRY(ff,tt,A,ds,es,A_by_es,dr0);
		aa1=(temp_cell+htc*ds/k_cell*temp_infinity)/(1.0+htc*ds/k_cell);
		//aa2=-htc/k_cell*(temp_face-temp_infinity);
		F_PROFILE(ff, tt, ii) =aa1;

	}
	end_f_loop(ff, tt)
		bb4=bb4;
}
DEFINE_PROPERTY(vapor_density, c0, t0)
{
	double ru;
	//temp_cell= C_UDSI_M1(c0, t0, 0);
	//ru=P_op1/R_R/temp_cell;
	temp_cell= C_UDSI(c0, t0, 0);
	ru=P_op2/R_R/temp_cell;
	//Message("P_op=%f R_R=%f temp_cell=%f RU=%f \n", P_op1,R_R,temp_cell,ru);
	return ru;
}
DEFINE_PROPERTY(liquid_density, c0, t0)
{
	double ru;
	ru=M_l2/(wick_volume);
	//Message("M_l1=%f  RU=%f \n", M_l1,ru);
	return ru;
}
DEFINE_PROPERTY(conductivity_water, c0, t0)
{
	double kkk;
	C_CENTROID(x,c0,t0);

	if ((x[1] <= Y_wick1) && (x[1] >= Y_wick2))
	{
		kkk=K_wick1;

	}
	else if ((x[1] < Y_wick2) && (x[1] >= Y_wick3))
	{
		kkk=K_wick2;

	}
	else 
	{
		kkk=K_wick3;

	}

	return kkk;
}
DEFINE_DIFFUSIVITY(diffusivity_water, c0, t0, i)
{
	double DDD;
	C_CENTROID(x,c0,t0);

	if ((x[1] <= Y_wick1) && (x[1] >= Y_wick2))
	{
		DDD=D_wick1;
	}
	else if ((x[1] < Y_wick2) && (x[1] >= Y_wick3))
	{
		DDD=D_wick2;
	}
	else 
	{
		DDD=D_wick3;
	}

	return DDD;
}
DEFINE_PROFILE(Porosity_water, t0, ii)
{
	double ppp;
	begin_c_loop(c0, t0)
	{
		C_CENTROID(x,c0,t0);

		if ((x[1] <= Y_wick1) && (x[1] >= Y_wick2))
		{
			ppp=Porousity1;
		}
		else if ((x[1] < Y_wick2) && (x[1] >= Y_wick3))
		{
			ppp=Porousity2;
		}
		else 
		{
			ppp=Porousity3;
		}

		C_PROFILE(c0, t0, ii) =ppp;

	}
	end_c_loop(ff, tt)
}
DEFINE_PROFILE(viscous_resistance, t0, ii)
{
	double v_re;
	begin_c_loop(c0, t0)
	{
		C_CENTROID(x,c0,t0);

		if ((x[1] <= Y_wick1) && (x[1] >= Y_wick2))
		{
			v_re=Visc_Re_1;
		}
		else if ((x[1] < Y_wick2) && (x[1] >= Y_wick3))
		{
			v_re=Visc_Re_2;
		}
		else 
		{
			v_re=Visc_Re_3;
		}

		C_PROFILE(c0, t0, ii) =v_re;

	}
	end_c_loop(ff, tt)
}
DEFINE_PROFILE(intertial_resistance, t0, ii)
{
	double int_re;
	begin_c_loop(c0, t0)
	{
		C_CENTROID(x,c0,t0);

		if ((x[1] <= Y_wick1) && (x[1] >= Y_wick2))
		{
			int_re=Inter_Re_1;
		}
		else if ((x[1] < Y_wick2) && (x[1] >= Y_wick3))
		{
			int_re=Inter_Re_2;
		}
		else 
		{
			int_re=Inter_Re_3;
		}

		C_PROFILE(c0, t0, ii) =int_re;

	}
	end_c_loop(ff, tt)
}
DEFINE_UDS_UNSTEADY(wick_uds_unsteady,c,t,i,apu,su)
{
	real physical_dt, vol, rho, rho_old, phi_old;
	physical_dt = RP_Get_Real("physical-time-step");
	zone_ID = THREAD_ID (t);

	vol = C_VOLUME(c,t);

	C_CENTROID(x,c,t);

	if ((x[1] <= Y_wick1) && (x[1] >= Y_wick2))
	{
		Porousity=Porousity1;
	}
	else if ((x[1] < Y_wick2) && (x[1] >= Y_wick3))
	{
		Porousity=Porousity2;
	}
	else 
	{
		Porousity=Porousity3;
	}

	RO_L=Ru_wick;
	if (zone_ID == wick_core_ID)
	{
		aa1=(1.0-Porousity)*RO_S*CP_S+Porousity*RO_L*CP_L;
		aa2=aa1/CP_L;
		rho = aa2;
		rho_old=aa2;
	}
	else
	{
		rho = C_R(c,t);
		rho_old=C_R_M1(c,t);
	}

	*apu = -rho*vol/physical_dt;   
	phi_old = C_STORAGE_R(c,t,SV_UDSI_M1(i));
	*su  = rho_old*vol*phi_old/physical_dt;  
}
DEFINE_PROFILE(flux_heat_dif_cond, tt, ii)
{
	begin_f_loop(ff, tt)
	{
		F_PROFILE(ff, tt, ii) =heat_flux;
	}
	end_f_loop(ff, tt)
}
DEFINE_DELTAT(mydeltat,d)
{
	real time_step;
	real flow_time = CURRENT_TIME;
	time_step=0.01*pow(1.20,N_time-1);
	time_step=floor(time_step*1000)/1000;
	if (time_step >= 1)
	{
		time_step=1.0;
	}
	
	return time_step;
} 
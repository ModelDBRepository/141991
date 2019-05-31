#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

//#####################################################################
/* Runge Kutta integrator from numerical recipies plus improvements */
/* void *deri(int n,double h[],double D[],double t);  */
/* function argument not tested yet */

void rk4(void deri(int , double [], double [], double ), \
		 double h[], int n, double t, double dt)
{
#define naux 60 
	
	int i;
	double k1[naux],k2[naux],k3[naux],k4[naux],h0[naux];
	double dt2, dt6;
	
	dt2=dt/2.;
	dt6=dt/6.;
	
	for (i = 0 ; i<n; i++)
		h0[i] = h[i];
	
	deri(n,h0,k1,t);
	for (i =0 ; i<n ; i++)
		h0[i]=h[i]+dt2*k1[i];
	
	deri(n,h0,k2,t+dt2);
	for (i =0 ; i<n ; i++)
		h0[i]=h[i]+dt2*k2[i];
	
	deri(n,h0,k3,t+dt2);
	for (i =0 ; i<n ; i++)
		h0[i]=h[i]+dt*k3[i];
	
	deri(n,h0,k4,t+dt);
	for (i = 0 ; i<n ; i++)
	{h0[i]=h[i]+dt*k4[i];};
	
	for (i =0; i<n ; i++)
		h[i]=h[i]+dt6*(2.*(k2[i]+k3[i])+k1[i]+k4[i]);
	
	return;
}
/*End of Runge Kutta Code*/
//#####################################################################



//#####################################################################
/*Global variables definitions*/

#define Nneurons 32
#define PI 3.1415926535897932384626433832795028841971 
#define Ntrials 5
#define Npoints 20
#define Nspikes 25
#define TRUE 1
#define FALSE 0
#define Ndiscret 1000
#define N_steps_sta 1000

double t_q,R_r,R_st,t_r,R_st,t_st,R_ss,phi,phi1,phi2,phi01,phi02,freq,nu,nu1,nu2,tau_d,t_k; 
double Rss,RA,C0,C1,S0,S1;
double V,a,b,c,w,z,nn,p,m,h,r,x,y,xx,yy,xxx,yyy,xxxx,yyyy;
double cur_term,IA,ILT,IHT,INA,IH,ILK,syn_term,syn1ex,syn2ex,syn1in,syn2in,CM;
double a_inf,b_inf,c_inf,w_inf,z_inf,n_inf,p_inf,m_inf,h_inf,r_inf,tau_a,tau_b,tau_c,tau_w,tau_z,tau_n,tau_p,tau_m,tau_h,tau_r;
double tauex1,tauex2,tauinh1,tauinh2,GEXC1,GEXC2,GINH1,GINH2;

/* Neuron model with synaptic inputs equations */

void takens(int n, double v[], double dv[], double t)
{
 V=v[0];
 a=v[1];
 b=v[2];
 c=v[3];
 w=v[4];
 z=v[5];
 nn=v[6];
 p=v[7];
 m=v[8];
 h=v[9];
 r=v[10];
 x=v[11];
 y=v[12];
 xx=v[13];
 yy=v[14];
 xxx=v[15];
 yyy=v[16];
 xxxx=v[17];
 yyyy=v[18];
                                                                                            
dv[0]=0.001*(cur_term-IA-ILT-IHT-INA-IH-ILK-syn_term)/CM;
dv[1]=(a_inf-v[1])/tau_a;
dv[2]=(b_inf-v[2])/tau_b;
dv[3]=(c_inf-v[3])/tau_c;
dv[4]=(w_inf-v[4])/tau_w;
dv[5]=(z_inf-v[5])/tau_z;
dv[6]=(n_inf-v[6])/tau_n;
dv[7]=(p_inf-v[7])/tau_p;
dv[8]=(m_inf-v[8])/tau_m;
dv[9]=(h_inf-v[9])/tau_h;
dv[10]=(r_inf-v[10])/tau_r;
dv[11]=v[12];
dv[12]= - 2.*v[12]/tauex1 - 1.*v[11]/(tauex1*tauex1);
dv[13]=v[14];
dv[14]= - 2.*v[14]/tauex2 - 1.*v[13]/(tauex2*tauex2);
dv[15]=v[16];
dv[16]= - 2.*v[16]/tauinh1 - 1.*v[15]/(tauinh1*tauinh1);
dv[17]=v[18];
dv[18]= - 2.*v[18]/tauinh2 - 1.*v[17]/(tauinh2*tauinh2);

return;
}


/* This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input as a pointer to the function to be 
 integrated between limits a and b, also input. When called with	n=1, the routine returns the crudest estimate of integral from a 
 to b of f(x)dx. Subsequent calls with n=2,3,...(in that sequential order) will improve the accuracy by adding 2n-2 additional interior 
 points. */

#define FUNC(x) ((*func)(x))

double trapzd(double (*func)(double), double a, double b, int n){		

	double x,tnm,sum,del;
	static double s;
	int it,j;
	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
	for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm; 						//This is the spacing of the points to be added.
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm); 				//This replaces s by its refined value.
		return s;
	}
}


/* This routine returns the integral of the function func from a to b. The parameters EPS can be set to the desired fractional accuracy and 
 JMAX so that 2 to the power JMAX-1 is the maximum allowed number of steps. Integration is performed by the trapezoidal rule. */ 

#define EPS 1.0e-3
#define JMAX 20
double qtrap(double (*func)(double), double a, double b){ 
double trapzd(double (*func)(double), double a, double b, int n);
int j;
double s,olds=0.0; 					//Initial value of olds is arbitrary.

	for (j=1;j<=JMAX;j++) {
		s=trapzd(func,a,b,j);
		if (j > 5) 				//Avoid spurious early convergence.
		if (fabs(s-olds) < EPS*fabs(olds)||(s == 0.0 && olds == 0.0)) return s;
		olds=s;
	}
	return 0.0; 					//Never get here.
}

/*#############################################################################################################
Definition of the spontaneuos figing rate functions to be evaluated by a time dependent poission process scheme. 
##############################################################################################################*/ 

double lambdaIpsi(double t){
	return(( R_ss*(exp(phi1*sin(2.*PI*freq*(t)))/exp(phi1))+phi01)*nu1*(1.-0.55*exp(-(t-t_k)/(0.8))-0.45*exp(-(t-t_k)/(25.)))); //0.0008 & 0.09 //Vs=0.88: 0.0025 & 0.13
}

double lambdaContra(double t){
	return(( R_ss*(exp(phi2*sin(2.*PI*freq*(t)))/exp(phi2))+phi02)*nu2*(1.-0.55*exp(-(t-t_k)/(0.8))-0.45*exp(-(t-t_k)/(25.))));//*(1.-exp(nu*(t-t-1))));//*C0*exp(-(t)/S0)); //0.0008 & 0.09 //Vs=0.88: 0.0025 & 0.13
}

/*Other functions:
double lambda(double t){
	return( (1.-exp(-(t+t_k)/t_q))*( R_r*exp(-(t+t_k)/t_r) + R_st*exp(-(t+t_k)/t_st) + R_ss*(exp(phi*sin(2.*PI*freq*(t+t_k)))/exp(phi)))*(1.-exp(-nu*(t))) ); //0.0008 & 0.09 //Vs=0.88: 0.0025 & 0.13
}
double lambda1(double t){
	return( R_ss*(exp(sin(2.*PI*freq*(t+t_k)))/(exp(phi)))*(1.-exp(-nu*(t))));
}
*/


/*########################################
 Local variables definition and main MSO model code
######################################## */

int main(){
FILE *prnt0, *prnt1, *prnt2, *prnt3, *prnt4,*prnt5;
double v[19], e1[Nneurons], e2[Nneurons],e3[Nneurons],old1[Nneurons],old2[Nneurons],old3[Nneurons];
double I0, IEXT,q0p17,q3p03;
double startex,EEXC,startinh,EINH,ZETA,PHI;
double EK,ENA,EH,ELK,GNA,GHT,GLT,GA,GH,GLK,VREST,RREST,TAUM,VTH,SM,GE22,GE38;
double t,t_int,dt,tTOT,T,ITDmax;
double sigma1,sigma2,sigma3,d1,d2,d3,v1,v2,v3,integral;
double u,u1,u2,R,Phi,nu0;
int i,j,k,l,o,h,aa,pp;
int N,decision,Npeak,N_over_trials;
int Nsmallepsp1,Nsmallepsp2,Nsmallepsp3,k1,k2,k3;
double tis,R1x,R1y,R2x,R2y,R3x,R3y,R1,R2,R3,y2;
double x[Nneurons],y[Nneurons],yy[Nneurons];
double Gtotal, gA,gHT,gLT,gNA,gH,gLK;
double V0,w0,z0,n0,p0,m0,h0,r0;	
double max,twenty,eighty,t1,t2,t3,t4,peak,v0;	
double aver_epsp[Ndiscret],aver_epsg_ips[Ndiscret],aver_epsg_con[Ndiscret],aver_ipsg[Ndiscret],sta[N_steps_sta],memory[N_steps_sta];

//#######################	
//Output Files definition
//#######################
	
prnt1=fopen("Model2_syn-events-ipsi.dat","wt");
prnt2=fopen("Model2_syn-events-contra.dat","wt");	
prnt3=fopen("Model2_train_EPSGs-EPSPs.dat","wt");
prnt4=fopen("Model2_ITD-curve.dat","w");
prnt5=fopen("Model2_STA_epsg_epsp_iklt.dat","wt");

	/*Paramter values for lambda functions*/
	
	freq=.5;//in KHz
	phi=7.70 - 3.69*freq + 0.40*freq*freq;//3.79 - 1.08*freq + 0.07*freq*freq;
	t_q=0.2;
	t_r=3.;
	t_st=10.;
	R_r=600.;
	R_st=200.;
	R_ss=2296 - 977*freq + 105*freq*freq;
	
	nu=(10./8.)*(R_r+R_st+2.*R_ss);
	
	nu1=0.006;//0.0065;//0.0018;
	nu2	=0.0009;//0.0042;//0.0005;//0.001;
	phi1=32.7525;//7.8606;//phi*1.45;//phi1=phi*10.,nu1=0.01875,VS=0.879479,%=1.00875,Slope=4.82404,HW=0.586,GEXC=1.15//phi1=phi*10.,nu1=0.003,VS=0.889077,%=0.34125,Slope=4.92063,HW=0.602
	phi2=1.5483;//1.7865;//1.0719;//phi*.25;//phi2=phi*.1,nu2=0.00105,VS=0.236316,%=1.03,Slope=1.72839,HW=0.32//phi2=phi*.25,nu2=0.00185,VS=0.482553,%=0.97125,Slope=2.09915,HW=0.422//phi2=phi*.2,nu2=0.00035;VS=0.483799,%=0.35125,Slope=1.64128,HW=0.436
	phi01=0.0;//-0.05;
	phi02=0.;//-0.1;//0.0;	
	
//	32.7525 1.7865  0       0       0.005   0.0042  0.196   -0.24
	
	
//#######################
//Experimental setting for MSO neuron cell    
//#######################
I0=0.0;
IEXT=0.0;
q0p17=0.17;//1.0;
q3p03=3.03;//1.0;
	

//#######################	
//parameters for syn_term:
//#######################	
	GEXC1=3.;//3.;//2.5;//.9;//.5;//.50;//1.35;
	GEXC2=4.5;//4.5;//4.8;//.5;//.50;//1.35;
GINH1=0.0;//.70;
GINH2=0.0;
syn1ex=0.;
syn2ex=0.;
syn1in=0.;
syn2in=0.;
tauex1=0.10;
tauex2=0.10;
tauinh1=.10;
tauinh2=.10;
EEXC=0.0;
EINH=-90.;//-75.;
ZETA=0.5;
PHI=0.85;
                                                                                            
// Type-II and IA=0
CM=0.012;
EK=-70.0;
ENA=55.0;
EH=-43.0;
ELK=-65.0;
GNA=1000.0;
GHT=150.0;
GLT=200.0;
GA=0.0;
GH=20.0;
GLK=2.0;

//#######################	
//Simulations time values.
//#######################	
T=(1./(freq*1.));
tau_d=.75;
tTOT=T*(Nspikes*1.+2.);
dt=T/Ndiscret;
ITDmax=.6;	
srand48(time(0));	

    o=0;
    for(o=0;o<=Npoints;o++){
		l=0;
		N_over_trials=0;
		for(l=0;l<Ntrials;l++){
			
			//#####################################
			//Inizialization of dynamical variables
			//#####################################
			v[0]=-63.628399;
			v[1]=0.;
			v[2]=0.;
			v[3]=0.;
			v[4]=0.51221395;
			v[5]=0.66181272;
			v[6]=0.0077282726;
			v[7]=0.0011447769;
			v[8]=0.025057629;
			v[9]=0.44309759;
			v[10]=0.14586952;
			v[11]=0.;
			v[12]=0.;
			v[13]=0.;
			v[14]=0.;
			v[15]=0.;
			v[16]=0.;
			v[17]=0.;
			v[18]=0.;
		
			V0=v[0];
			w0=v[4];	
			z0=v[5];
			n0=v[6];
			p0=v[7];
			m0=v[8];
			h0=v[9];
			r0=v[10];
	
			N=0;
			v1=0.;
			v2=0.;
			v3=0.;
		
			d1=o*ITDmax*T/(1.*Npoints);
			d2=(ITDmax/2.)*T;
			d3=0.0;		

			Nsmallepsp1=0;
			Nsmallepsp2=0;
			Nsmallepsp3=0;
			k=0;
			k1=0;
			k2=0;
			k3=0;
			R1x=0.;
			R1y=0.;
			R2x=0.;
			R2y=0.;
			R3x=0.;
			R3y=0.;
		
			i=0;
			for(i=0;i<Nneurons;i++){
        		e1[i]=drand48()*(T/5.);
        		e2[i]=drand48()*(T/5.);
        		e3[i]=drand48()*(T/5.);
        		old1[i]=0.;
        		old2[i]=0.;
        		old3[i]=0.;
			}
			
			for(aa=0;aa<Ndiscret;aa++){
				aver_epsp[aa]=v[0];
				aver_epsg_ips[aa]=v[11];
				aver_epsg_con[aa]=v[13];
				aver_ipsg[aa]=v[17];
			}
			
			t=0.;
			aa=0;

			//%%%%%%//Main simulation time loop//%%%%%%%%%
			
			while(t<tTOT){
				Npeak=floor(t/T);
				
				if(fabs(fmod(t,T)-2.)<dt)aa=0;
				
				if(o==Npoints){
					aver_epsp[aa]=aver_epsp[aa]+v[0];
					aver_epsg_ips[aa]=aver_epsg_ips[aa]+v[11];
					aver_epsg_con[aa]=aver_epsg_con[aa]+v[13];
					aver_ipsg[aa]=aver_ipsg[aa]+v[17];
				}
			
				//Decision to make an event on the ipsilateral input
				i=0;
				for(i=0;i<Nneurons;i++){
					decision=TRUE;
					if((fabs(e1[i]-t+d1)<0.5*dt)){
						u=drand48();
						t_k=e1[i];
						t_int=tau_d;
						decision=FALSE;
						while(decision == FALSE){
							if(t_int>tTOT)decision=TRUE;
							
							integral=qtrap(&lambdaIpsi,t_k,t_k+t_int);
							
							if(((log(u)+integral)>=0.)){
								old1[i]=e1[i];
								v[12]=v[12]+GEXC1*exp(1.)/tauex1;
								e1[i]=t_int+old1[i];
								fprintf(prnt1,"%lg\t%d\t%lg\t%lg\n",t,i,e1[i],lambdaIpsi);
								decision=TRUE;
							}
							
							t_int=t_int+10.*dt;
							
						}
						
						k1=k1+1;
						Nsmallepsp1=k1;
						
						if(k1>1){
							R1x=R1x+cos(e1[i]*2.*PI*freq);
							R1y=R1y+sin(e1[i]*2.*PI*freq);
						}
						
					}
				}
				
				//Decision to make an event on the contralateral input
				j=0;
				for(j=0;j<Nneurons;j++){
					decision=TRUE;
					if((fabs(e2[j]-t+d2)<0.5*dt)){
						u=drand48();
						t_k=e2[j];
						t_int=tau_d;
						decision=FALSE;
						while(decision == FALSE){
							if(t_int>tTOT)decision=TRUE;
							
							integral=qtrap(&lambdaContra,t_k,t_k+t_int);
							
							if((log(u)+integral)>=0. ){
								old2[j]=e2[j];
								v[14]=v[14]+GEXC2*exp(1.)/tauex2;
								e2[j]=t_int+old2[j];
								fprintf(prnt2,"%lg\t%d\t%lg\t%lg\n",t,j,e2[j],lambdaContra);
								decision=TRUE;
							}
							
							t_int=t_int+10.*dt;
							
						}
						k2=k2+1;
						Nsmallepsp2=k2;
						
						if(k2>1){
							R2x=R2x+cos(e2[j]*freq*2.*PI);
							R2y=R2y+sin(e2[j]*freq*2.*PI);
						}
					}
				}


					//#######################
					//MSO Point Neuron Model	
					//#######################
					// Fast transient K+ current
					a_inf=1./(pow((1.0+exp(-(v[0]+31.0)/6.0)),0.25));
					b_inf=1./(pow((1.0+exp((v[0]+66.0)/7.0)),0.5));
					c_inf=b_inf;
					tau_a=q0p17*(100.0/(7.0*exp((v[0]+60.0)/14.0)+29.0*exp(-(v[0]+60.0)/24.0))+0.1);
					tau_b=q0p17*(1000.0/(14.0*exp((v[0]+60.0)/27.0)+29.0*exp(-(v[0]+60.0)/24.0))+1.0);
					tau_c=q0p17*(90.0/(1.0+exp(-(v[0]+66.0)/17.0))+10.0);
																						
					// Low threshold K+ current
					w_inf=1./(pow((1.0+exp(-(v[0]+48.0)/6.0)),0.25));
					z_inf=(1.0-ZETA)/(1.0+exp((v[0]+71.0)/10.0))+ZETA;
					tau_w=q0p17*(90.0/(6.0*exp((v[0]+60.0)/6.0)+16.0*exp(-(v[0]+60.0)/45.0))+0.5);
					tau_z=q0p17*(1000.0/(exp((v[0]+60.0)/20.0)+exp(-(v[0]+60.0)/8.0))+50.0);
                                                                                            
					// High threshold K+ current
					n_inf=1./(pow((1.0+exp(-(v[0]+15.0)/5.0)),0.5));
					p_inf=1.0/(1.0+exp(-(v[0]+23.0)/6.0));
					tau_n=q0p17*(100.0/(11.0*exp((v[0]+60.0)/24.0)+21.0*exp(-(v[0]+60.0)/23.0))+0.7);
					tau_p=q0p17*(100.0/(4.0*exp((v[0]+60.0)/32.0)+5.0*exp(-(v[0]+60.0)/22.0))+5.0);
                                                                                            
					// Fast Na+ current
					m_inf=1.0/(1.0+exp(-(v[0]+38.0)/7.0));
					h_inf=1.0/(1.0+exp((v[0]+65.0)/6.0));
					tau_m=q0p17*(10.0/(5.0*exp((v[0]+60.0)/18.0)+36.0*exp(-(v[0]+60.0)/25.0))+0.04);
					tau_h=q0p17*(100.0/(7.0*exp((v[0]+60.0)/11.0)+10.0*exp(-(v[0]+60.0)/25.0))+0.6);
                                                                                            
					// Hyperpolarization-activated cation current
					r_inf=1.0/(1.0+exp((v[0]+76.0)/7.0));
					tau_r=q0p17*(100000./(237.0*exp((v[0]+60.0)/12.0)+17.0*exp(-(v[0]+60.0)/14.0))+25.0);
                                                                                            
					// Current equations
					IA=0.0;//GA*(v[1]*v[1]*v[1]*v[1])*v[2]*v[3]*(v[0]-EK)
					ILT=1.3*q3p03*GLT*(v[4]*v[4]*v[4]*v[4])*v[5]*(v[0]-EK);//q3p03*GLT*(w0*w0*w0*w0)*z0*(v[0]-EK);//q3p03*GLT*(v[4]*v[4]*v[4]*v[4])*v[5]*(v[0]-EK);
					IHT=q3p03*GHT*(PHI*v[6]*v[6]+(1.0-PHI)*v[7])*(v[0]-EK);
					INA=q3p03*GNA*(v[8]*v[8]*v[8])*v[9]*(v[0]-ENA);
					IH=1.29*q3p03*GH*v[10]*(v[0]-EH);
					ILK=q3p03*GLK*(v[0]-ELK);
                                                                                            
					cur_term=I0;

					syn1ex = v[11]*(v[0]-EEXC);
					syn2ex = v[13]*(v[0]-EEXC);
					syn1in = v[15]*(v[0]-EINH);
					syn2in = v[17]*(v[0]-EINH);
					syn_term= syn1ex + syn2ex + syn1in +syn2in;

					rk4(takens,v,19,t,dt);
			
				//####################//Creating the spike triggered average//################
				pp=0;
				for(pp=0;pp<N_steps_sta;pp++){
					memory[N_steps_sta-pp]=memory[N_steps_sta - pp - 1];
				}	
				memory[1]=v[0];
				
				v3=v2;
				v2=v1;
				v1=v[0];
				if(((v2-v3)>0.)&&((v1-v2)<0.)&&(v2>-20.)){
					N=N+1;
					printf("%d\n",N);
					pp=0;
					for(pp=0;pp<N_steps_sta;pp++){
						sta[pp]=sta[pp]+memory[N_steps_sta-pp];
					}
				}
				
				
				if(o==Npoints){
					fprintf(prnt3,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",t,v[0],sin(2.*PI*freq*(t)),syn_term,v[11],v[13]);
				}
				
				t=t+dt;
				aa=aa+1;

		}
		fprintf(prnt2,"\n");
		
		//####################//Creating the Vector strength for each input//################
			
		N_over_trials=N_over_trials+N;	
			
		R1 = R1+(sqrt(R1x*R1x+R1y*R1y)/(1.*(Nsmallepsp1-1)));
		R2 = R2+(sqrt(R2x*R2x+R2y*R2y)/(1.*(Nsmallepsp2-1)));
		R3 = R3+(sqrt(R3x*R3x+R3y*R3y)/(1.*(Nsmallepsp3-1)));

	
	}
	
		prnt3=fopen("Model2_ITD-curve.dat","a");
		fprintf(prnt4,"%lg\t%lg\n",d1-(ITDmax/2.)*T,((1.*N_over_trials)/(1.*Ntrials))/Nspikes);
		fclose(prnt4);	
		}
		
		//####################//Creating the average slope and halfwidth//################
	
		//REMEMBER TO PUT A "-" SIGN IN FRONT OF THE VECTOR VALUE IN CASE OF IPSP
		
		v0=(aver_epsp[1]/((1.*Nspikes+3.)));//-(aver_epsp[1]/((1.*Nspikes+3.)));
		v3=v2=v1=v0;
		max=0.;
		for(aa=1;aa<(Ndiscret-1);aa++){
			fprintf(prnt5,"%lg\t%lg\t%lg\t%lg\t%lg\n",aa*dt,aver_epsp[aa]/(1.*Nspikes+3.),aver_epsg_ips[aa]/(1.*Nspikes+3.),aver_epsg_con[aa]/(1.*Nspikes+3.),aver_ipsg[aa]/(1.*Nspikes+3.));
			v3=v2;
			v2=v1;
			v1=(aver_epsp[aa]/((1.*Nspikes+3.)));//-(aver_epsp[aa]/((1.*Nspikes+3.)));
			if((v3<v2)&&(v1<v2)&&(max<(v2-v0))){
				peak=aa*dt;
				max=v2-v0;
			}
		}
		
		for(aa=1;aa<(Ndiscret-1);aa++){
			v3=v2;
			v2=v1;
			v1=(aver_epsp[aa]/((1.*Nspikes+3.)));//-(aver_epsp[aa]/((1.*Nspikes+3.)));
			if(((v3-v0)<(max/2.))&&((max/2.)<(v2-v0))&&(aa*dt<(peak))){t1=aa*dt;}
			if(((v2-v0)<(max/2.))&&((max/2.)<(v3-v0))&&(aa*dt<(peak))){t2=aa*dt;}
			if(((v3-v0)<(max*0.2))&&((max*0.2)<(v2-v0))&&(aa*dt<(peak))){t3=aa*dt;twenty=v2-v0;}
			if(((v3-v0)<(max*0.8))&&((max*0.8)<(v2-v0))&&(aa*dt<(peak))){t4=aa*dt;eighty=v2-v0;}	
		}		
    
    R1 = sqrt( (R1x/(1.*Nsmallepsp1-1.))*(R1x/(1.*Nsmallepsp1-1.)) + (R1y/(1.*Nsmallepsp1-1.))*(R1y/(1.*Nsmallepsp1-1.)) );
    R2 = sqrt( (R2x/(1.*Nsmallepsp2-1.))*(R2x/(1.*Nsmallepsp2-1.)) + (R2y/(1.*Nsmallepsp2-1.))*(R2y/(1.*Nsmallepsp2-1.)) );
    R3 = sqrt( (R3x/(1.*Nsmallepsp3-1.))*(R3x/(1.*Nsmallepsp3-1.)) + (R3y/(1.*Nsmallepsp3-1.))*(R3y/(1.*Nsmallepsp3-1.)) );	
	
	//####################//Saving all the results on terminal//################
    
	printf("		Phi	R_ss	VS	Prcnt.	\n");
	printf("ipsi_exc	");
    printf("%lg\t%lg\t%lg\t%lg\n",phi1,R_ss*nu1,R1,(1.*Nsmallepsp1)/(1.*Nneurons)/Nspikes);
	printf("cont_exc	");
    printf("%lg\t%lg\t%lg\t%lg\n",phi2,R_ss*nu2,R2,(1.*Nsmallepsp2)/(1.*Nneurons)/Nspikes);
	printf("cont_inh	");
    printf("%lg\t%lg\t%lg\t%lg\n",phi,R_ss,R3,(1.*Nsmallepsp3)/(1.*Nneurons)/Nspikes);		
	
	printf("\n");
	
	fclose(prnt1);
	fclose(prnt2);
	fclose(prnt3);
	fclose(prnt4);
	fclose(prnt5);
	
	Gtotal=(q3p03*GLT*(w0*w0*w0*w0)*z0+q3p03*GHT*(PHI*n0*n0+(1.0-PHI)*p0)+q3p03*GNA*(m0*m0*m0)*h0+q3p03*GH*r0+q3p03*GLK); 
		
	printf("membrane time constant= ");
	printf("%lg",(CM/Gtotal)*1000. );
	printf("msec \n");
		
	printf("\n");
	printf("Slope	Norm-Slope	Halfwidth	\n");
	printf("%lg\t%lg\t%lg\n",(eighty-twenty)/(t4-t3),((eighty-twenty)/(t4-t3))/max,t1-t2);


exit(0);

}




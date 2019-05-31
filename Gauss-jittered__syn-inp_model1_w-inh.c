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


#define Nneurons_exc 128 //256
#define Nneurons_inh 128
#define PI 3.1415926535897932384626433832795028841971 
#define Ntrials 20 //20
#define Npoints 60 //40
#define Nspikes 40 //20 //10
#define N_steps_sta 10000
#define Ndiscret 1000

double V,a,b,c,w,z,nn,p,m,h,r,x,y,xx,yy,xxx,yyy,xxxx,yyyy;
double cur_term,IA,ILT,IHT,INA,IH,ILK,syn_term,syn1ex,syn2ex,syn1in,syn2in,CM;
double a_inf,b_inf,c_inf,w_inf,z_inf,n_inf,p_inf,m_inf,h_inf,r_inf,tau_a,tau_b,tau_c,tau_w,tau_z,tau_n,tau_p,tau_m,tau_h,tau_r;
double tauex1,tauex2,tauinh1,tauinh2,GEXC1,GEXC2,GINH1,GINH2;

void takens(int n, double v[], double dv[], double t)
/* Takens equations */
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
dv[12]= - 2.*v[12]/tauex1 - v[11]/(tauex1*tauex1);
dv[13]=v[14];
dv[14]= - 2.*v[14]/tauex2 - v[13]/(tauex2*tauex2);
dv[15]=v[16];
dv[16]= - 2.*v[16]/tauinh1 - v[15]/(tauinh1*tauinh1);
dv[17]=v[18];
dv[18]= - 2.*v[18]/tauinh2 - v[17]/(tauinh2*tauinh2);

return;
}

int main(){
FILE *prnt1, *prnt2, *prnt3, *prnt4,*prnt5;
double v[19], e1[Nneurons_exc], e2[Nneurons_exc],e3[Nneurons_inh],e4[Nneurons_inh];
double sta[N_steps_sta],sta_stdev[N_steps_sta],memory[N_steps_sta],iklt_memory[N_steps_sta],gklt_memory[N_steps_sta],epsp_gklt_ave[N_steps_sta],uper[N_steps_sta],lower[N_steps_sta],epsp_ave[N_steps_sta],epsp_ave_lower[N_steps_sta],epsp_ave_upper[N_steps_sta],epsp_stdev[N_steps_sta],iklt_ave[N_steps_sta],gklt_ave[N_steps_sta],ipsp_ave[N_steps_sta];
double syn1[N_steps_sta],syn2[N_steps_sta],syn_ex1[N_steps_sta], syn_ex2[N_steps_sta],syn_spike_ex1[N_steps_sta],syn_spike_ex2[N_steps_sta];
int old1[Nneurons_exc],old2[Nneurons_exc],old3[Nneurons_inh],old4[Nneurons_inh];
double T_exc_ipsi[Nneurons_exc][Nspikes],T_exc_contra[Nneurons_exc][Nspikes],T_inh_ipsi[Nneurons_inh][Nspikes],T_inh_contra[Nneurons_inh][Nspikes];
double A_exc_ipsi[Nneurons_exc][Nspikes],A_exc_contra[Nneurons_exc][Nspikes],A_inh_ipsi[Nneurons_inh][Nspikes],A_inh_contra[Nneurons_inh][Nspikes];
double I0, IEXT,q0p17,q3p03;
double startex,EEXC,startinh,EINH,ZETA,PHI,tauex,GEXC;
double EK,ENA,EH,ELK,GNA,GHT,GLT,GA,GH,GLK,VREST,RREST,TAUM,VTH,SM,GE22,GE38;
double t,dt,tTOT,T,freq,v0,ITDmax;
double u1,u2,R,Phi;
double sigmaT1,sigmaT2,sigmaT3,sigmaT4,sigmaA1,sigmaA2,sigmaA3,sigmaA4,v1,v2,v3,integral;
double release_now,inh_ipsi_delay,exc_contra_delay,inh_contra_delay,exc_ipsi_delay;
double prob_max_ipsi_exc,prob_max_contra_exc,prob_max_ipsi_inh,prob_max_contra_inh,threshold1,threshold2,threshold3,threshold4,lambda1,lambda2,lambda3,lambda4,nu1,nu2,nu3,nu4;
int i,j,k,l,o,h,p,n_cycles;
int N,decision,Npeak,N_over_trials,P;
int Nsmallepsp1,Nsmallepsp2,Nsmallepsp3,Nsmallepsp4,k1,k2,k3,k4;
double R1x,R1y,R2x,R2y,R3x,R3y,R4x,R4y,R1,R2,R3,R4,y2,R_ss;
double Gtotal, gA,gHT,gLT,gNA,gH,gLK;
double V0,w0,z0,n0,p0,m0,h0,r0;


prnt1=fopen("Model1_STA_epsg_epsp_iklt_w-inh.dat","w");//ITD_w-IKLT_09_-099_09_500.dat
fclose(prnt1);
prnt2=fopen("Model1_train_EPSGs-EPSPs_w-inh.dat","wt");
prnt3=fopen("Model1_ITD-curve_w-inh.dat","w");//ITD_w-IKLT_w-INH_09_-099_09_500.dat
fclose(prnt3);

srand48(time(0));

//*********************************//
//Definition os the local constants//
//*********************************//																						
freq=.5;//in KHz

I0=0.0;
IEXT=0.0;
q0p17=0.17;//1.0;
q3p03=3.03;//1.0;
                                                                                            
//parameters for syn_term:
syn1ex=0.;
syn2ex=0.;
syn1in=0.;
syn2in=0.;
tauex1=0.1;
tauex2=0.1;
tauinh1=0.4;
tauinh2=0.4;
GEXC1=2.25;//2.25//0.7;//15.;//0.72;//.55;//.50;//1.35;
GEXC2=2.25;//2.25//(GEXC1*tauex1)/tauex2;
GINH1=0.0;//0.05;//0.05;//.70;
GINH2=1.4;//1.2;//0.40;//0.17;//0.05;//.7;

EEXC=0.0;
EINH=-90.;//-75.0;
ZETA=0.5;
PHI=0.85;

// Type-II and IA=0
CM=0.012;// in nF
EK=-70.0;
ENA=55.0;
EH=-43.0;
ELK=-65.0;
GNA=1000.0;// in nS
GHT=150.0;
GLT=200.;//200.0;
GA=0.0;
GH=20.;//20.0;
GLK=2.;//2.0;

T=(1./(freq*1.));
tTOT=T*(Nspikes*1.+2.);
dt=T/Ndiscret;
ITDmax=1.5;
sigmaT1=.1*T;
sigmaT2=.1*T;
sigmaT3=.1*T;
sigmaT4=.1*T;
sigmaA1=.1*GEXC;
sigmaA2=.1*GEXC;
sigmaA3=.1*GINH1;
sigmaA4=.1*GINH2;

R1=0.;
R2=0.;
R3=0.;
R4=0.;

//************//
//Loop in ITDs//
//************//
o=0;
for(o=0;o<=Npoints;o++){
	//**************//
	//Loop in trials//
	//**************//
	for(p=0;p<N_steps_sta;p++){
		sta[p]=0.;
		sta_stdev[p]=0.;
		epsp_ave[p]=0.;
		epsp_stdev[p]=0.;
		uper[p]=-120.;
		epsp_ave_upper[p]=-120.;
		lower[p]=60.;
		epsp_ave_lower[p]=60.;
		iklt_ave[p]=0.;
		gklt_ave[p]=0.;
		epsp_gklt_ave[p]=0.;
		syn_ex1[p]=0.; 
		syn_ex2[p]=0.;
		syn_spike_ex1[p]=0.;
		syn_spike_ex2[p]=0.;
	}
	N_over_trials=0;
	Npeak=0;
	
	h=0;
	for(h=0;h<Ntrials;h++){
		//********************************//
		//Initianlization of the variables//
		//********************************//
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
		P=0;
		v1=0.;
		v2=0.;
		v3=0.;
		for(p=0;p<N_steps_sta;p++){
			memory[p]=-63.628399;
			syn1[p]=0.;
			syn2[p]=0.;
			iklt_memory[p]=q3p03*GLT*(v[4]*v[4]*v[4]*v[4])*v[5]*(-63.628399-EK);
			gklt_memory[p]=q3p03*GLT*(v[4]*v[4]*v[4]*v[4])*v[5];
		}
		
		Nsmallepsp1=0;
		Nsmallepsp2=0;
		Nsmallepsp3=0;
		Nsmallepsp4=0;
		k=0;
		k1=0;
		k2=0;
		k3=0;
		k4=0;
		R1x=0.;
		R1y=0.;
		R2x=0.;
		R2y=0.;
		R3x=0.;
		R3y=0.;
		R4x=0.;
		R4y=0.;
		
		i=0;
		for(i=0;i<Nneurons_exc;i++){
			e1[i]=0.;
			e2[i]=0.;
			old1[i]=0;
			old2[i]=0;
		}
		
		i=0;
		for(i=0;i<Nneurons_inh;i++){
			e3[i]=0.;
			e4[i]=0.;
			old3[i]=0;
			old4[i]=0;
		}
		
		//*********************************************//
		//EPSPs and IPSPs time delay generation for ITD//		
		//*********************************************//				
		
		exc_ipsi_delay=o*ITDmax*T/(1.*Npoints);
		inh_ipsi_delay=o*ITDmax*T/(1.*Npoints)+0.2;
		inh_contra_delay=(ITDmax/2.)*T-0.6;
		exc_contra_delay=(ITDmax/2.)*T;
				
		//************//
		//Loop in time//
		//************//
		release_now=0.;
		t=0.;
		while(t<tTOT){
			
			//******************************************************************************//
			// lambda  //  nu  //  R=Vect_strength //  % of synapses active pre input-cycle // 
			// 0.75		  95       0.95					0.2									//	
			// 0.9		  24	   0.98					0.2									//
			// 0.55		  230      0.9					0.2									//
 			// 0.05		  800      0.8					0.2									//
			// -0.35	  1300     0.7					0.2									//	
			// -0.75	  1750     0.6					0.2									//	
			// -0.99	  2500     0.5					0.2									//	
			//
			//******************************************************************************//
		
			lambda1=0.975;
			lambda2=0.975;
			lambda3=0.975;
			lambda4=0.975;
			threshold1=0.9;//0.55;//0.75;//0.90;
			threshold2=-.99;//0.55;//0.75;//-0.99;
			threshold3=-.99;//-0.75;//0.9;
			threshold4=-.99;//-0.75;//0.9;//.75
			nu1=24;
			nu2=2500;
			nu3=2500;
			nu4=2500;
			
			prob_max_ipsi_exc=sin(2.*PI*(t - exc_ipsi_delay)*freq);
			prob_max_contra_exc=sin(2.*PI*(t - exc_contra_delay)*freq);
			prob_max_ipsi_inh=sin(2.*PI*(t - inh_ipsi_delay)*freq);
			prob_max_contra_inh=sin(2.*PI*(t - inh_contra_delay)*freq);
			if(t<exc_contra_delay){prob_max_ipsi_exc=-1.;}
			if(t<exc_ipsi_delay){prob_max_contra_exc=-1.;}
			if(t<inh_ipsi_delay){prob_max_ipsi_inh=-1.;}
			if(t<inh_contra_delay){prob_max_contra_inh=-1.;}
			
		
			//IPSILATERAL EPSPs
				i=0;
				for(i=0;i<Nneurons_exc;i++){
					release_now=(nu1*(2.*drand48())-1.);//24 for 0.90//95 for 0.75///230 for 0.55//800 for 0.05//1250 -0.35//1700 for -0.75// 2300 for -0.99 //2500 for -0.999 //  
					if((threshold1 < release_now )&&(release_now < prob_max_ipsi_exc)){
						if((e1[i])<t){
							v[12]=v[12]+GEXC1*exp(1.)/tauex1;
							old1[i]=e1[i];
							e1[i]=t;
							//fprintf(prnt4,"%lg\t%d\n",e1[i],i);
							k1=k1+1;
							Nsmallepsp1=k1;
							if(k1>1){
								R1x=R1x+cos(e1[i]*2.*PI*freq);
								R1y=R1y+sin(e1[i]*2.*PI*freq);
							}
						}
					}
				}
			
			//CONTRALATERAL EPSPs	
				i=0.;
				for(i=0;i<Nneurons_exc;i++){
					release_now=(nu2*(2.*drand48())-1.);
					if((threshold2 < release_now )&&(release_now < prob_max_contra_exc)){
						if((e2[i])<t){
							v[14]=v[14]+GEXC2*exp(1.)/tauex2;
							old2[i]=e2[i];
							e2[i]=t;
							//fprintf(prnt4,"%lg\t%d\n",e2[i],i);
							k2=k2+1;
							Nsmallepsp2=k2;
							if(k2>1){
								R2x=R2x+cos(e2[i]*2.*PI*freq);
								R2y=R2y+sin(e2[i]*2.*PI*freq);
							}
						}
					}
				}
			
			//IPSILATERAL IPSPs	
				i=0.;
				for(i=0;i<Nneurons_inh;i++){
					release_now=(nu3*(2.*drand48())-1.);
					if((threshold3 < release_now )&&(release_now < prob_max_ipsi_inh)){
						if((e3[i])<t){
							v[16]=v[16]+GINH1*exp(1.)/tauinh2;
							old3[i]=e3[i];
							e3[i]=t;
							//fprintf(prnt4,"%lg\t%d\n",e3[i],i);
							k3=k3+1;
							Nsmallepsp3=k3;
							if(k3>1){
								R3x=R3x+cos(e3[i]*2.*PI*freq);
								R3y=R3y+sin(e3[i]*2.*PI*freq);
							}
						}
					}
				}
						
			//CONTRALATERAL IPSPs	
			i=0.;
				for(i=0;i<Nneurons_inh;i++){
					release_now=(nu4*(2.*drand48())-1.);//80//1700//2200
					if((threshold4 < release_now )&&(release_now < prob_max_contra_inh)){
						if((e4[i])<t){
							v[18]=v[18]+GINH2*exp(1.)/tauinh2;
							old4[i]=e4[i];
							e4[i]=t;
							//fprintf(prnt4,"%lg\t%d\n",e4[i],i);
							k4=k4+1;
							Nsmallepsp4=k4;
							if(k4>1){
								R4x=R4x+cos(e4[i]*2.*PI*freq);
								R4y=R4y+sin(e4[i]*2.*PI*freq);
							}
						}
					}
				}	
				
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
			tau_w=0.8*q0p17*(100.0/(6.0*exp((v[0]+60.0)/6.0)+16.0*exp(-(v[0]+60.0)/45.0))+1.5);
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
			gA=0.0;//GA*(v[1]*v[1]*v[1]*v[1])*v[2]*v[3]*(v[0]-EK)
			gLT=1.3*q3p03*GLT*(v[4]*v[4]*v[4]*v[4])*v[5];//1.6*q3p03*GLT*(v[4]*v[4]*v[4]*v[4])*v[5];//1.6*q3p03*GLT*(w0*w0*w0*w0)*z0;//1.6*q3p03*GLT*(v[4]*v[4]*v[4]*v[4])*v[5];
			gHT=q3p03*GHT*(PHI*v[6]*v[6]+(1.0-PHI)*v[7]);
			gNA=q3p03*GNA*(v[8]*v[8]*v[8])*v[9];
			gH=1.3*q3p03*GH*v[10];//1.56*q3p03*GH*r0;//1.56*q3p03*GH*v[10];
			gLK=q3p03*GLK;
            
			IA=gA*(v[0]-EK);
			ILT=gLT*(v[0]-EK);
			IHT=gHT*(v[0]-EK);
			INA=gNA*(v[0]-ENA);
			IH=gH*(v[0]-EH);
			ILK=gLK*(v[0]-ELK);
                                                                                            
			cur_term=I0;
			
			syn1ex = v[11]*(v[0]-EEXC);
			syn2ex = v[13]*(v[0]-EEXC);
			syn1in = v[15]*(v[0]-EINH);
			syn2in = v[17]*(v[0]-EINH);
			syn_term= syn1ex + syn2ex + syn1in +syn2in;

			rk4(takens,v,19,t,dt);
			
			//******************************************************************************************************//
			//Loop to compute Spike triggered voltage, synaptic conductance, IKLT current, IKLT conductance averages//
			//******************************************************************************************************//
			
			
			for(p=0;p<(N_steps_sta-1);p++){
				memory[N_steps_sta - p -1] = memory[N_steps_sta - p - 2];
				syn1[N_steps_sta - p -1]=syn1[N_steps_sta - p - 2];
				syn2[N_steps_sta - p -1]=syn2[N_steps_sta - p - 2];
				iklt_memory[N_steps_sta - p -1] = iklt_memory[N_steps_sta - p - 2];
				gklt_memory[N_steps_sta - p -1] = gklt_memory[N_steps_sta - p - 2];
			}	
			memory[0]=v[0];
			iklt_memory[0]=ILT;
			gklt_memory[0]=q3p03*GLT*(v[4]*v[4]*v[4]*v[4])*v[5];
			syn1[0]=syn1ex;
			syn2[0]=syn2ex;
				
			v3=v2;
			v2=v1;
			v1=v[0];
			
			if(((v2-v3)>0.)&&((v1-v2)<0.)){
				P=P+1;
				p=0;
				for(p=0;p<N_steps_sta;p++){
					epsp_ave[p]=epsp_ave[p]+memory[N_steps_sta-p-1];
					epsp_stdev[p]=epsp_stdev[p]+memory[N_steps_sta-p-1]*memory[N_steps_sta-p-1];
					epsp_gklt_ave[p]=epsp_gklt_ave[p]+gklt_memory[N_steps_sta - p -1];
					syn_ex1[p]=syn_ex1[p] + syn1[N_steps_sta-p-1]; 
					syn_ex2[p]=syn_ex2[p] + syn2[N_steps_sta-p-1];
					if(memory[N_steps_sta-p-1]<epsp_ave_lower[p]){epsp_ave_lower[p]=memory[N_steps_sta-p-1];}
					if(epsp_ave_upper[p]<memory[N_steps_sta-p-1]){epsp_ave_upper[p]=memory[N_steps_sta-p-1];}
				}	
			}
			if(((v2-v3)>0.)&&((v1-v2)<0.)&&(v2>-25.)){
				N=N+1;
				p=0;
				for(p=0;p<N_steps_sta;p++){
					sta[p]=sta[p] + memory[N_steps_sta-p-1];
					sta_stdev[p]=sta_stdev[p] + memory[N_steps_sta-p-1]*memory[N_steps_sta-p-1];
					iklt_ave[p]=iklt_ave[p] + iklt_memory[N_steps_sta-p-1];
					gklt_ave[p]=gklt_ave[p] + gklt_memory[N_steps_sta-p-1];
					syn_spike_ex1[p]=syn_spike_ex1[p] + syn1[N_steps_sta-p-1]; 
					syn_spike_ex2[p]=syn_spike_ex2[p] + syn2[N_steps_sta-p-1];
					if(memory[N_steps_sta-p-1]<lower[p]){lower[p]=memory[N_steps_sta-p-1];}
					if(uper[p]<memory[N_steps_sta-p-1]){uper[p]=memory[N_steps_sta-p-1];}
				}
				
			}
			
			if(o==Npoints){
				fprintf(prnt2,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",t,v[0],sin(2.*PI*freq*(t)),syn_term,v[11],v[13]);
			}
			
			t=t+dt;
			
		}
		
		//***********************************//
		//END OF THE LOOP FOR THE RK INTEGRATOR
		//***********************************//
		
		N_over_trials=N_over_trials+N;
		Npeak=Npeak+P;
		
		R1 = R1+(sqrt(R1x*R1x+R1y*R1y)/(1.*(Nsmallepsp1-1)));
		R2 = R2+(sqrt(R2x*R2x+R2y*R2y)/(1.*(Nsmallepsp2-1)));
		R3 = R3+(sqrt(R3x*R3x+R3y*R3y)/(1.*(Nsmallepsp3-1)));
		R4 = R4+(sqrt(R4x*R4x+R4y*R4y)/(1.*(Nsmallepsp4-1)));

	}
	
	prnt1=fopen("Model1_STA_epsg_epsp_iklt_w-inh.dat","a");
	p=0;
	for(p=0;p<N_steps_sta;p++){
		fprintf(prnt1,"%lg\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",exc_ipsi_delay-(ITDmax/2.)*T,p,1.*sta[p]/(1.*N_over_trials),lower[p],uper[p],epsp_ave[p]/(1.*Npeak),epsp_ave_lower[p],epsp_ave_upper[p],sqrt((epsp_stdev[p]/(1.*Npeak))-(epsp_ave[p]/(1.*Npeak))*(epsp_ave[p]/(1.*Npeak))),sqrt((sta_stdev[p]/(1.*N_over_trials))-(sta[p]/(1.*N_over_trials))*(sta[p]/(1.*N_over_trials))),iklt_ave[p]/(1.*N_over_trials),gklt_ave[p]/(1.*N_over_trials),epsp_gklt_ave[p]/(1.*Npeak),syn_ex1[p]/(1.*Npeak),syn_ex2[p]/(1.*Npeak),syn_spike_ex1[p]/(1.*N_over_trials),syn_spike_ex2[p]/(1.*N_over_trials));
	}
	
	fprintf(prnt1,"\n");	
	fclose(prnt1);
	
	prnt3=fopen("Model1_ITD-curve_w-inh.dat","a");
	fprintf(prnt3,"%lg\t%lg\n",exc_ipsi_delay-(ITDmax/2.)*T,(1.*N_over_trials)/(1.*Ntrials));
	fclose(prnt3);
	
}	
	    
    R1 = sqrt( (R1x/(1.*Nsmallepsp1-1.))*(R1x/(1.*Nsmallepsp1-1.)) + (R1y/(1.*Nsmallepsp1-1.))*(R1y/(1.*Nsmallepsp1-1.)) );
    R2 = sqrt( (R2x/(1.*Nsmallepsp2-1.))*(R2x/(1.*Nsmallepsp2-1.)) + (R2y/(1.*Nsmallepsp2-1.))*(R2y/(1.*Nsmallepsp2-1.)) );
    R3 = sqrt( (R3x/(1.*Nsmallepsp3-1.))*(R3x/(1.*Nsmallepsp3-1.)) + (R3y/(1.*Nsmallepsp3-1.))*(R3y/(1.*Nsmallepsp3-1.)) );		
    R4 = sqrt( (R4x/(1.*Nsmallepsp4-1.))*(R4x/(1.*Nsmallepsp4-1.)) + (R4y/(1.*Nsmallepsp4-1.))*(R4y/(1.*Nsmallepsp4-1.)) );
	 
	printf("Exc Ipsi Vect-Stren= ");
	printf("%lg\t%lg\n",R1,(1.*Nsmallepsp1)/(1.*Nneurons_exc)/Nspikes);
	printf("Exc Contra Vect-Stren= ");
    printf("%lg\t%lg\n",R2,(1.*Nsmallepsp2)/(1.*Nneurons_exc)/Nspikes);
	printf("Inh Ipsi Vect-Stren= ");
	;printf("%lg\t%lg\n",R3,(1.*Nsmallepsp3)/(1.*Nneurons_inh)/Nspikes);
	printf("Inh Contra Vect-Stren= ");
	printf("%lg\t%lg\n",R4,(1.*Nsmallepsp4)/(1.*Nneurons_inh)/Nspikes);		

fclose(prnt1);
fclose(prnt2);
fclose(prnt3);

Gtotal=(q3p03*GLT*(w0*w0*w0*w0)*z0+q3p03*GHT*(PHI*n0*n0+(1.0-PHI)*p0)+q3p03*GNA*(m0*m0*m0)*h0+q3p03*GH*r0+q3p03*GLK); 

printf("membrane time constant= ");
printf("%lg",(CM/Gtotal)*1000. );
printf("msec \n");
printf("membrane resistance= ");
printf("%lg",(1./Gtotal));
printf("Ohm/cm^2 \n");

exit(0);

}

#include "pair_airebo.h"
#include "declarations.h"
#include "pair_tersoff.h"

//**************************************
//***** Global variables*********
//**************************************
double LNB;
double LCC;
double Rcut_List; 
int MaxNBond;
int MaxNeighb;
int MaxList;
double Normal_cuttoff;
int Frame_limit;
double shift;
int N_periodic;
int print_terms;
int ILP_flag, N_ILP_Layers;
int N_layers;
//*********list
double Rcut_List_long,Rcut_List_short;
int *List_long,*NList_long;int MaxList_long;
//**************************************
//***** Tapper and cutt-off variables***
//**************************************
double  swa,swb,non_bond_cut, swa2,swa3,swb2,swb3,Tap_7,Tap_6,Tap_5,Tap_4,Tap_3,Tap_2,Tap_1,Tap_0;
double  swa_Rep,swb_Rep,non_bond_cut_Rep, swa2_Rep,swa3_Rep,swb2_Rep,swb3_Rep,Tap_7_Rep,Tap_6_Rep,Tap_5_Rep,Tap_4_Rep,Tap_3_Rep,Tap_2_Rep,Tap_1_Rep,Tap_0_Rep;

double  alpha0N,alpha0C,alpha0B,alpha0H,C60N,C60B,C60H,C60C,VeffN,VeffB,VeffH,VeffC,alphaN;
double  alphaB,alphaH,alphaC,NN_C6,BB_C6,HH_C6,CC_C6,BN_C6,BH_C6,NH_C6,CN_C6,CB_C6,CH_C6,r0N,r0B,r0H,r0C,r_eff_NN,r_eff_BB,r_eff_HH,r_eff_CC,r_eff_BN,r_eff_BH,r_eff_NH;
double  r_eff_CN,r_eff_CB,r_eff_CH,d_TS,Sr_TS;
double  CC_C_vdW,BN_C_vdW,NN_C_vdW,BB_C_vdW,HH_C_vdW,CB_C_vdW,CN_C_vdW,CH_C_vdW,BH_C_vdW,NH_C_vdW;
double  NN_rvdW_long,BN_rvdW_long,BB_rvdW_long,BH_rvdW_long,NH_rvdW_long,HH_rvdW_long,CH_rvdW_long,CB_rvdW_long,CN_rvdW_long,CC_rvdW_long;
double  H_Epsilon_long,B_Epsilon_long,N_Epsilon_long,BH_Epsilon_long,NH_Epsilon_long,BN_Epsilon_long,CH_Epsilon_long,CN_Epsilon_long,C_Epsilon_long,CB_Epsilon_long;
double  H_alpha_long,B_alpha_long,N_alpha_long,BH_alpha_long,NH_alpha_long,BN_alpha_long,CH_alpha_long,CN_alpha_long,C_alpha_long,CB_alpha_long;
double  H_Trans,B_Trans,N_Trans,BH_Trans,NH_Trans,BN_Trans,CH_Trans,CN_Trans,C_Trans,CB_Trans;
double  original_H_gamma,original_B_gamma,original_N_gamma,original_C_gamma,H_Etha,B_Etha,N_Etha,C_Etha;
double  H_Kai,N_Kai,B_Kai,C_Kai,H_gamma,B_gamma,N_gamma,C_gamma,Kappa,Kappa_eV_kcal;
// ***********************Tersoff Params*********************
double BN_Tersoff_c,BN_Tersoff_d,BN_Tersoff_sqr_c,BN_Tersoff_sqr_d,BN_Tersoff_h,BN_Tersoff_n,BN_Tersoff_A,BN_Tersoff_B,BN_Tersoff_beta,BN_Tersoff_lambda1,BN_Tersoff_lambda2,BN_Tersoff_lambda3,BN_Tersoff_D,BN_Tersoff_R;
double NN_Tersoff_c,NN_Tersoff_d,NN_Tersoff_sqr_c,NN_Tersoff_sqr_d,NN_Tersoff_h,NN_Tersoff_n,NN_Tersoff_A,NN_Tersoff_B,NN_Tersoff_beta,NN_Tersoff_lambda1,NN_Tersoff_lambda2,NN_Tersoff_lambda3,NN_Tersoff_D,NN_Tersoff_R;
double BB_Tersoff_c,BB_Tersoff_d,BB_Tersoff_sqr_c,BB_Tersoff_sqr_d,BB_Tersoff_h,BB_Tersoff_n,BB_Tersoff_A,BB_Tersoff_B,BB_Tersoff_beta,BB_Tersoff_lambda1,BB_Tersoff_lambda2,BB_Tersoff_lambda3,BB_Tersoff_D,BB_Tersoff_R;
double CC_Tersoff_c,CC_Tersoff_d,CC_Tersoff_sqr_c,CC_Tersoff_sqr_d,CC_Tersoff_h,CC_Tersoff_n,CC_Tersoff_A,CC_Tersoff_B,CC_Tersoff_beta,CC_Tersoff_lambda1,CC_Tersoff_lambda2,CC_Tersoff_lambda3,CC_Tersoff_D,CC_Tersoff_R;
double CB_Tersoff_c,CB_Tersoff_d,CB_Tersoff_sqr_c,CB_Tersoff_sqr_d,CB_Tersoff_h,CB_Tersoff_n,CB_Tersoff_A,CB_Tersoff_B,CB_Tersoff_beta,CB_Tersoff_lambda1,CB_Tersoff_lambda2,CB_Tersoff_lambda3,CB_Tersoff_D,CB_Tersoff_R;
double CN_Tersoff_c,CN_Tersoff_d,CN_Tersoff_sqr_c,CN_Tersoff_sqr_d,CN_Tersoff_h,CN_Tersoff_n,CN_Tersoff_A,CN_Tersoff_B,CN_Tersoff_beta,CN_Tersoff_lambda1,CN_Tersoff_lambda2,CN_Tersoff_lambda3,CN_Tersoff_D,CN_Tersoff_R;
//**************************************
//*****Tersoff Global variables*********
//**************************************
int *map;                        // 0 (C), 1 (H), or -1 (NULL) for each type
int me;
int ljflag,torflag;              // 0/1 if LJ,torsion terms included
double cutlj;                    // user-specified LJ cutoff
double cutljrebosq;              // cut for when to compute
double **cutljsq;                // LJ cutoffs for C,H types
double **lj1,**lj2,**lj3,**lj4;  // pre-computed LJ coeffs for C,H types
double cut3rebo;                 // maximum distance for 3rd REBO neigh
int maxlocal;                    // size of numneigh, firstneigh arrays
int pgsize;                      // size of neighbor page
int oneatom;                     // max # of neighbors for one atom
int *REBO_numneigh;              // # of pair neighbors for each atomint **REBO_firstneigh;           // ptr to 1st neighbor of each atom
double *closestdistsq;           // closest owned atom dist to each ghost
double *nC,*nH;                  // sum of weighting fns with REBO neighs
double smin,Nmin,Nmax,NCmin,NCmax,thmin,thmax;
double rcmin[2][2],rcmax[2][2],rcmaxsq[2][2],rcmaxp[2][2];
double Q[2][2],alpha[2][2],A[2][2],rho[2][2],BIJc[2][2][3],Beta[2][2][3],A_ij[2][2];
double rcLJmin[2][2],rcLJmax[2][2],rcLJmaxsq[2][2],bLJmin[2][2],bLJmax[2][2];
double epsilon[2][2],sigma[2][2],epsilonT[2][2];
  // spline coefficients

double gCdom[5],gC1[4][6],gC2[4][6],gHdom[4],gH[3][6];
double pCCdom[2][2],pCHdom[2][2],pCC[4][4][16],pCH[4][4][16];
double piCCdom[3][2],piCHdom[3][2],piHHdom[3][2];
double piCC[4][4][9][64],piCH[4][4][9][64],piHH[4][4][9][64];
double Tijdom[3][2],Tijc[4][4][9][64];

  // spline knot values

double PCCf[5][5],PCCdfdx[5][5],PCCdfdy[5][5],PCHf[5][5];
double PCHdfdx[5][5],PCHdfdy[5][5];
double piCCf[5][5][11],piCCdfdx[5][5][11];
double piCCdfdy[5][5][11],piCCdfdz[5][5][11];
double piCHf[5][5][11],piCHdfdx[5][5][11];
double piCHdfdy[5][5][11],piCHdfdz[5][5][11];
double piHHf[5][5][11],piHHdfdx[5][5][11];
double piHHdfdy[5][5][11],piHHdfdz[5][5][11];
double Tf[5][5][10],Tdfdx[5][5][10],Tdfdy[5][5][10],Tdfdz[5][5][10];


/********************************************************************/
/* Routines required for the calculation of the potential energy.   */
/* Written by Itai Leven June 2013.                                 */
/********************************************************************/

/**************************************************/
/*                    E_{bond}                    */
/**************************************************/
struct timeb t_initial,t_final,caliper[3];

int main()
{
  /********** Variable declarations **********/
  
  time_t  t0; /* time_t is defined on <time.h> and <sys/types.h> as long */
  clock_t c0; /* clock_t is defined on <time.h> and <sys/types.h> as int */
  int i,j,atomi,Interlayer,Periodic_flag, Coulomb_flag,Lock_COM_flag,LRI_flag;
  int CalcType, Natoms, NConserve, NTraj,nthreads,layer1,layer2;
  int *NBond, *Bonded,*NBond_once, *Bonded_once, *BondNum,*List,*NList,Max_NBond,Max_NNeighb,Max_NList;
  int *Normal_atom, *NNeighb, *Neighb,p,Input_flag,potential_flag;
  int BN_flag,Ncar_hyd,N_init,N_init_ILP;  
  double dt,tini, tend, T, OptTol;
  double *rij;
  double  *angle,*rij_long,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,Term;           //Coulomb term array.
  double *Bmat, *qvec, *dqvec, *Avec, *Fc; //Coulomb term arrays.
  double Energy,Energy_;
  char atomtype[3];
  char Simulation[50];
  point L;
  point *R0;
  point *Fbond;
  long idum;
  Randomize();
  idum = -random();
  FILE *logfile = fopen("log.dat","w");
  AtomParamsStruct AtomParams[NAtomParams];
  BondParamsStruct BondParams[NAtomParams*NAtomParams];
  
  /********** Initialization **********/
  
  BN_flag=0;
  print_terms=0;
  
  read_ILP_parameters();
  read_Tersoff_parameters();  
  read_RUN_parameters(tini,tend,dt,T,NConserve,NTraj,OptTol,CalcType,Interlayer,Periodic_flag,nthreads,potential_flag,Input_flag,Lock_COM_flag,LRI_flag,layer1,layer2);
  
  if(CalcType == 0)strcpy(Simulation,"MD simulation");
  else if(CalcType == 1)strcpy(Simulation,"Single point calculation");
  else if(CalcType == 4)strcpy(Simulation,"FIRE Optimization");
  else if(CalcType == 5)strcpy(Simulation,"Global Registry Index + LRI calculations");
  else {
    cerr<<"Wrong value for CalcType in main() exiting"<<endl;
    exit(0);
  }
  
  fprintf(logfile,"Starting %s\n",Simulation);
  
    //**Temp**
  MaxList_long   = MaxList;
  Rcut_List_long = Rcut_List;
  Rcut_List_short = 5;
 if(N_ILP_Layers==0)N_ILP_Layers=100000000;
  //********
  
  Init_Atom_Params(AtomParams);
  Init_Bond_Params(BondParams);
  
  Natoms              = GetNatoms("Coords.xyz",Input_flag);
  atom *particle      = new atom[Natoms];
  atom *temp_particle = new atom[Natoms];
  NBond               = new int[Natoms];
  NBond_once          = new int[Natoms];
  BondNum             = new int[1];
  NNeighb             = new int[Natoms];
  Neighb              = new int[Natoms*MaxNeighb];
  List                = new int[Natoms*MaxList];
  NList               = new int[Natoms];
  List_long           = new int[Natoms*MaxList_long];
  NList_long          = new int[Natoms];
  Fc                  = new double[Natoms*MaxNeighb];
  nC                  = new double[Natoms];
  nH                  = new double[Natoms];
  Bonded              = new int[Natoms*MaxNBond];
  Bonded_once         = new int[Natoms*MaxNBond];
  rij                 = new double[1];
  rij_long            = new double[1];
  //rij                 = new double[Natoms*MaxNBond];
  //rij_long            = new double[Natoms*MaxNeighb];
  Fbond               = new point[Natoms];
  R0                  = new point[Natoms];
  Normal_struct *Normal = new Normal_struct[Natoms];
  Normal_atom = new int[Natoms*3];
  Bmat             = new double[1];
  qvec             = new double[Natoms+1];
  dqvec            = new double[Natoms+1];
  Avec             = new double[Natoms+1];
  angle            = new double[1];
  spline_init();
  
  ReadCoor("Coords.xyz",Natoms,particle,Interlayer,L,Input_flag,potential_flag,temp_particle,Ncar_hyd,BN_flag);
  read_file("CH.airebo");
  
  PBC_Init(particle,Natoms,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L,R0);
  
  N_layers=0;
  
  for (atomi=0 ; atomi < Natoms ; atomi++){
    R0[atomi].x = particle[atomi].r[0];
    R0[atomi].y = particle[atomi].r[1];
    R0[atomi].z = particle[atomi].r[2];
    particle[atomi].r[0] -= L.x * floor((particle[atomi].r[0] - Xmin) / (L.x + TINY));
    particle[atomi].r[1] -= L.y * floor((particle[atomi].r[1] - Ymin) / (L.y + TINY));
    particle[atomi].r[2] -= L.z * floor((particle[atomi].r[2] - Zmin) / (L.z + TINY));
    if(particle[atomi].layer > N_layers)N_layers=particle[atomi].layer;
  }
  
  N_layers++;
  
  if (nthreads != 0)
    omp_set_num_threads(nthreads);
  
  int milestone(0), rc;
  rc = ftime(&caliper[milestone]);
  //cerr<<"hey"<<endl;
    
  //*****************************
  fprintf(logfile,"\nNumber of atoms in simulation= %i\n",Natoms);
  fprintf(logfile,"\nNumber of layers in simulation= %i\n",N_layers);
  
  UpdateNeighbList(Natoms,particle,List,NList,R0,L);
  
  calc_rij(NBond,Bonded,BondNum,Natoms,particle,L,rij,Fc,BondParams,Neighb,NNeighb,rij_long,NList,List,NBond_once,Bonded_once,potential_flag);
  
  //calc_rij_long(Natoms,particle,L,Neighb,NNeighb);
  //cerr<<"yoyoyo1"<<endl;
  calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rij_long,NList,List,R0,NBond_once,Bonded_once,potential_flag);
  
  if (Interlayer == 1)
    {
      set_Normal_atom(rij,Natoms,Normal,particle,Normal_atom,L,NList,List);//The normal is calculating according to an atoms three nearest neighbors, Sp2 bonding is required appart from hydrogen atoms.
      calc_Normal(Natoms,Normal_atom,particle,Normal,L);
    }
  else if (Interlayer == 0){
#pragma omp parallel for
    for(atomi=0; atomi < Natoms; atomi++){
      particle[atomi].layer = 0;
      Normal[atomi].x = 0.0;
      Normal[atomi].y = 0.0;
      Normal[atomi].z = 0.0;
    }
  }
  else if(Interlayer == 2){
#pragma omp parallel for  
    for(atomi=0; atomi < Natoms; atomi++){
      Normal[atomi].x = 0.0;
      Normal[atomi].y = 0.0;
      Normal[atomi].z = 1.0;
    }
  }

  //Calc_Charge(particle,Natoms,Bmat,qvec,Avec,AtomParams,BondParams,rij,rij_long,logfile,L);
  /********** Calculate **********/
  //cerr<<"yoyoyo"<<endl;
  if(potential_flag==2 || potential_flag==3)N_init=0;
  else N_init=Ncar_hyd;
  
  if(potential_flag==1 || potential_flag==3)N_init_ILP=0;
  else N_init_ILP=Ncar_hyd;

  if(CalcType == 0){
    vel(particle,T,Natoms,&idum); 
    Propagate(Natoms,NBond,Bonded,particle,AtomParams,BondParams,BondNum,rij,L,dt,tini,tend,T,NConserve,NTraj,Bmat,qvec,Avec,dqvec,rij_long,logfile,t0,c0,angle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,potential_flag,R0,NList,List,NBond_once,Bonded_once,Lock_COM_flag,BN_flag,Ncar_hyd,N_init,N_init_ILP);
  }
  else if(CalcType == 1){
    
    FILE *TrajFile = fopen("Single_point_traj.xyz","w");
    
    
    point *Fterm = new point[Natoms];
        
    //calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rij_long,NList,List,R0,NBond_once,Bonded_once,potential_flag);
    
    //if(Interlayer == 1)
    //calc_Normal(Natoms,Normal_atom,particle,Normal,L);

    print_terms=1;
    
    Energy = Calc_Potential(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,rij_long,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,potential_flag,Energy_,NBond_once,Bonded_once,BN_flag,Ncar_hyd,N_init,N_init_ILP);
    
    fprintf(logfile,"\nPotential ToT_Energy= %.16f\nEnergy_per_atom= %.16f\n",Energy,Energy/Natoms);

    if(LRI_flag==0)Print_Traj(TrajFile,Natoms,0,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
    
    fclose(TrajFile);
    delete [] Fterm;
  }
  else if(CalcType == 2){
    //steepest_descent(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,NTraj,NConserve,OptTol,Bmat,qvec,Avec,dqvec,rij_long,logfile,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,Normal,Normal_atom,Interlayer,Fc,Periodic_flag,Neighb,NNeighb,potential_flag,NList,List,R0); 
  }
  else if(CalcType == 3){
    //conj_grad(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,OptTol,Bmat,qvec,Avec,dqvec,rij_long,logfile,NConserve, NTraj,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,Normal,Normal_atom,Interlayer,Fc,Periodic_flag,Neighb,NNeighb,potential_flag,NList,List,R0);
  }
  else if(CalcType ==4){

    fprintf(logfile,"\nInitial Energy:\n");
    print_terms=1;

    //calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rij_long,NList,List,R0,NBond_once,Bonded_once,potential_flag);
    
    //if(Interlayer == 1)
    //calc_Normal(Natoms,Normal_atom,particle,Normal,L);
    
    Energy = Calc_Potential(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,rij_long,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,potential_flag,Energy_,NBond_once,Bonded_once,BN_flag,Ncar_hyd,N_init,N_init_ILP);

    print_terms=0;
    
    vel(particle,T,Natoms,&idum); 
    Quench( Natoms,NBond,Bonded,particle,AtomParams,BondParams,BondNum,rij,L,dt,tini,tend,T,NConserve,NTraj,Bmat,qvec,dqvec,Avec,rij_long,logfile,t0,c0,angle,OptTol,Normal,Normal_atom,Interlayer,Fc,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,Neighb,NNeighb,Periodic_flag,potential_flag,NList,List,R0,NBond_once,Bonded_once,Lock_COM_flag,BN_flag,Ncar_hyd,N_init,N_init_ILP);
    
    fprintf(logfile,"\nFinal Energy:\n");

    calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rij_long,NList,List,R0,NBond_once,Bonded_once,potential_flag);
    
    if(Interlayer == 1)
      calc_Normal(Natoms,Normal_atom,particle,Normal,L);
    
    print_terms=1;
    
    Energy = Calc_Potential(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,rij_long,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,potential_flag,Energy_,NBond_once,Bonded_once,BN_flag,Ncar_hyd,N_init,N_init_ILP);
    
    fprintf(logfile,"\nToT_Energy= %.16f\nEnergy_per_atom= %.16f\n",Energy,Energy/Natoms);
    
  }
  else if(CalcType==5);
  else{
    cerr<<"Unrecognized requested calculation type! Ending session.\n";
    exit(0);
  }
  
  if (LRI_flag==1 || CalcType==5){
    double RI;
    FILE *TrajFile = fopen("Single_point_traj.xyz","w");
    FILE *TrajFile_RI    = fopen("Traj-RI.xyz","w");
    FILE *TrajFile_Inter = fopen("Traj-Inter.xyz","w");
    FILE *TrajFile_Intra = fopen("Traj-Intra.xyz","w");
    //LRI calculations
    calc_rij_LRI(Natoms,particle,L,rij,BondParams,Neighb,NNeighb,rij_long,NList_long,List_long,layer1,layer2);
    calc_RI(Natoms,rij_long,particle,BondParams,Normal,Neighb,NNeighb,L,RI,Normal_atom,layer1);
    Print_Traj_LRI(TrajFile,TrajFile_RI,TrajFile_Inter,TrajFile_Intra,Natoms,0,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L,layer1);
    fprintf(logfile,"Global Registry Index= %.16f\n",RI);
    fclose(TrajFile_RI);
    fclose(TrajFile_Intra);
    fclose(TrajFile_Inter);
    fclose(TrajFile); 
  }
   
  
  milestone++;
  rc = ftime(&caliper[milestone]);
  
  double delta_t(0);
  int hours,minutes,seconds,time;
  for (int i=0;i<milestone;++i)
    {
      delta_t = (caliper[i+1].time-caliper[i].time);
      delta_t += (caliper[i+1].millitm-caliper[i].millitm)/1000.0;
      //printf("time=%-.3f seconds using %i threads\n",delta_t,omp_get_max_threads());
      time    = int(delta_t);
      hours   = (time/3600);
      time    = time%3600;
      minutes = time/60;
      seconds  = time%60;
      fprintf(logfile,"time=%i:%i:%i h:m:s using %i threads\n",hours,minutes,seconds,omp_get_max_threads());
    }      

  fprintf(logfile,"\nNormal Termination\n");
  //**********************************************
  
  fclose(logfile);
  
  delete [] particle;
 
  delete [] temp_particle;
  
  delete [] NBond;
  delete [] Bonded;
  delete [] BondNum;
  delete [] rij;
  delete [] Fbond;
  delete [] angle;
  
  delete [] rij_long;
  
  delete [] NBond_once;
  delete [] Bonded_once;
  
  // Deallocate Coulomb term related arrays
  //delete [] Fterm;
  delete [] Bmat;
  delete [] qvec;
  delete [] dqvec;
  delete [] Avec;
  
}



double dFc_drn(int atomi,int atomj,int atomn,double r_ij,double Tersoff_R,double Tersoff_D,double dr_ij_drn){
  
  double dFc_ij_drn;
  
  if(r_ij < (Tersoff_R-Tersoff_D) || r_ij > (Tersoff_R+Tersoff_D)) dFc_ij_drn=0.0;
  else dFc_ij_drn = -0.5*cos(PIE/2*(r_ij-Tersoff_R)/Tersoff_D)*PIE*0.5/Tersoff_D*dr_ij_drn;
  
  return (dFc_ij_drn);
}

void calc_dE_REBO(point *F_REBO,double *angle,int *BondNum,int Natoms,point L,BondParamsStruct *BondParams,int *Bonded,atom *particle,double *Fc,int *NBond,double &E_REBO_,int *NBond_once,int *Bonded_once,int Ncar_hyd)
{
  double E_REBO=0.0;

    int index,atomk,typek,typej,typei,index2;
    double g,Theta_ijk,b_ij,r_ij,r_jk,Tersoff_lambda1,Tersoff_lambda2,Tersoff_lambda3,Tersoff_c,Tersoff_d,Tersoff_h,Tersoff_n,A,B,beta;
    double Xsi,ETersoff,dTheta_ijk_drn,dg_drn,dXsi_drn,Fc_ij,Fc_jk,Tersoff_R,Tersoff_D;
    double dr_ij_drn,dr_jk_drn,dFc_ij_drn,dFc_jk_drn,Fa,Fr,dFa_drn,dFr_drn,db_ij_drn,n,dE_Trsoff_ij_drn;
    double Tersoff_sqr_c,Tersoff_sqr_d,sqr_Tersoff_h_Cos_Theta;
    double inside_acos,cos_ijk;
    int i,j,k,m,ii,itype,jtype,atomi,atomj;
    int itag ,jtag;
    double delx,dely,delz,evdwl,fpair,xtmp,ytmp,ztmp;
    double rsq,rij,wij;
    double Qij,Aij,alphaij,VR,pre,dVRdi,VA,term,bij,dVAdi,dVA;
    double dwij,del[3];
    int vflag_atom=0;
    evdwl = 0.0;
    
#pragma omp for //reduction(+:E_REBO)
    for (atomi = 0; atomi < Ncar_hyd; atomi++) {
      itag=atomi;
      itype = particle[atomi].type;
      xtmp = particle[atomi].r[0];
      ytmp = particle[atomi].r[1];
      ztmp = particle[atomi].r[2];
      for (j = 0; j < NBond_once[atomi]; j++) {
	
	atomj = Bonded_once[atomi*MaxNBond + j];
	
	jtype = particle[atomj].type;
	
	delx = R_PBC(particle,atomi,atomj,0,L.x);
	dely = R_PBC(particle,atomi,atomj,1,L.y);	
	delz = R_PBC(particle,atomi,atomj,2,L.z);
	
	rsq = delx*delx + dely*dely + delz*delz;
	rij = sqrt(rsq);
	
	wij = Sp(rij,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
	
	if (wij <= TOL) continue;
	
	Qij = Q[itype][jtype];
	Aij = A_ij[itype][jtype];
	alphaij = alpha[itype][jtype];
	VR = wij*(1.0+(Qij/rij)) * Aij*exp(-alphaij*rij);
	pre = wij*Aij * exp(-alphaij*rij);
	dVRdi = pre * ((-alphaij)-(Qij/rsq)-(Qij*alphaij/rij));
	dVRdi += VR/wij * dwij;
	
	VA = dVA = 0.0;
	
	for (m = 0; m < 3; m++) {
	  term = -wij * BIJc[itype][jtype][m] * exp(-Beta[itype][jtype][m]*rij);
	  VA += term;
	  dVA += -Beta[itype][jtype][m] * term;
	}
	
	dVA += VA/wij * dwij;
	
	del[0] = delx;
	del[1] = dely;
	del[2] = delz;
	
	bij = bondorder(atomi,atomj,del,rij,VA,F_REBO,vflag_atom,particle,Bonded,BondNum,NBond,L);
	
	dVAdi = bij*dVA;
	
	fpair = -(dVRdi+dVAdi) / rij;
	
	
	F_REBO[atomi].x += delx*fpair;
	F_REBO[atomi].y += dely*fpair;
	F_REBO[atomi].z += delz*fpair;
	F_REBO[atomj].x -= delx*fpair;
	F_REBO[atomj].y -= dely*fpair;
	F_REBO[atomj].z -= delz*fpair;
	E_REBO_+=VR + bij*VA;
	//printf("Force VA=%.16f VR=%.16f bij=%.16f E_REBO=%.16f atomi=%i atomj=%i \n",VA,VR,bij,VR + bij*VA,atomi,atomj);
	//E_REBO_+=VR;
	// E_REBO_+=VR;
	
      }
    }
  
    //E_REBO_=E_REBO;
}


double bondorder(int i, int j, double rij[3], double rijmag, double VA, point *f, int vflag_atom,atom *particle,int *Bonded,int *BondNum,int *NBond,point L){

  int atomi,atomj,k,n,l,atomk,atoml,atomn,atom1,atom2,atom3,atom4;
  int itype,jtype,ktype,ltype,ntype;
  double rik[3],rjl[3],rkn[3],rji[3],rki[3],rlj[3],rknmag,dNki,dwjl,bij;
  double NijC,NijH,NjiC,NjiH,wik,dwik,dwkn,wjl;
  double rikmag,rjlmag,cosjik,cosijl,g,tmp2,tmp3;
  double Etmp,pij,tmp,wij,dwij,NconjtmpI,NconjtmpJ,Nki,Nlj,dS;
  double lamdajik,lamdaijl,dgdc,dgdN,pji,Nijconj,piRC;
  double dcosjikdri[3],dcosijldri[3],dcosjikdrk[3];
  double dN2[2],dN3[3];
  double dcosjikdrj[3],dcosijldrj[3],dcosijldrl[3];
  double Tij;
  double r32[3],r32mag,cos321,r43[3],r13[3];
  double dNlj;
  double om1234,rln[3];
  double rlnmag,dwln,r23[3],r23mag,r21[3],r21mag;
  double w21,dw21,r34[3],r34mag,cos234,w34,dw34;
  double cross321[3],cross234[3],prefactor,SpN;
  double fcijpc,fcikpc,fcjlpc,fcjkpc,fcilpc;
  double dt2dik[3],dt2djl[3],dt2dij[3],aa,aaa1,aaa2,at2,cw,cwnum,cwnom;
  double sin321,sin234,rr,rijrik,rijrjl,rjk2,rik2,ril2,rjl2;
  double dctik,dctjk,dctjl,dctij,dctji,dctil,rik2i,rjl2i,sink2i,sinl2i;
  double rjk[3],ril[3],dt1dik,dt1djk,dt1djl,dt1dil,dt1dij;
  double F23[3],F12[3],F34[3],F31[3],F24[3],fi[3],fj[3],fk[3],fl[3];
  double f1[3],f2[3],f3[3],f4[4];
  double dcut321,PijS,PjiS;
  double rij2,tspjik,dtsjik,tspijl,dtsijl,costmp;
  int *REBO_neighs,*REBO_neighs_i,*REBO_neighs_j,*REBO_neighs_k,*REBO_neighs_l;

  atomi = i;
  atomj = j;
  itype = particle[atomi].type;
  jtype = particle[atomj].type;

  wij = Sp(rijmag,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
  
  NijC = nC[i]-(wij*kronecker(jtype,0));
  NijH = nH[i]-(wij*kronecker(jtype,1));
  NjiC = nC[j]-(wij*kronecker(itype,0));
  NjiH = nH[j]-(wij*kronecker(itype,1));
  
  bij = 0.0;
  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  dgdc = 0.0;
  dgdN = 0.0;
  NconjtmpI = 0.0;
  NconjtmpJ = 0.0;
  Etmp = 0.0;
  
  //REBO_neighs = REBO_firstneigh[i];
  
    for(k=0 ; k < NBond[atomi] ; k++){ // Go over all j neighbors k != i.
      //atomk = REBO_neighs[k];
      atomk = Bonded[atomi*MaxNBond + k];
    if (atomk != atomj) {
      ktype = particle[atomk].type;
      //ktype = particle[atomk].type;
      rik[0] = R_PBC(particle,atomi,atomk,0,L.x);
      rik[1] = R_PBC(particle,atomi,atomk,1,L.y);
      rik[2] = R_PBC(particle,atomi,atomk,2,L.z);
      
      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      lamdajik = 4.0*kronecker(itype,1) *
        ((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag));
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dS);
      Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] -
        (wik*kronecker(itype,1));
      cosjik = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2])) /
        (rijmag*rikmag);
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosjik,(NijC+NijH),itype,&dgdc,&dgdN);
      Etmp = Etmp+(wik*g*exp(lamdajik));
      
      tmp3 = tmp3+(wik*dgdN*exp(lamdajik));
      NconjtmpI = NconjtmpI+(kronecker(ktype,0)*wik*Sp(Nki,Nmin,Nmax,dS)); //parameter Nmin Nmax
    }
  }

  PijS = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  PijS = PijSpline(NijC,NijH,itype,jtype,dN2);
  pij = pow(1.0+Etmp+PijS,-0.5);
  tmp = -0.5*cube(pij);
  //printf("pij=%.16f Etmp=%.16f PijS=%.16f \n",pij,Etmp,PijS);
  // pij forces
  
  //REBO_neighs = REBO_firstneigh[i];
  for (k = 0; k < NBond[atomi]; k++) {
    //atomk = REBO_neighs[k];
    atomk = Bonded[atomi*MaxNBond + k];
    if (atomk != atomj) {
      ktype = particle[atomk].type;
      rik[0] = R_PBC(particle,atomi,atomk,0,L.x);
      rik[1] = R_PBC(particle,atomi,atomk,1,L.y);
      rik[2] = R_PBC(particle,atomi,atomk,2,L.z);
      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      lamdajik = 4.0*kronecker(itype,1) *
        ((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag));
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
      cosjik = (rij[0]*rik[0] + rij[1]*rik[1] + rij[2]*rik[2]) /
        (rijmag*rikmag);
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      dcosjikdri[0] = ((rij[0]+rik[0])/(rijmag*rikmag)) -
        (cosjik*((rij[0]/(rijmag*rijmag))+(rik[0]/(rikmag*rikmag))));
      dcosjikdri[1] = ((rij[1]+rik[1])/(rijmag*rikmag)) -
        (cosjik*((rij[1]/(rijmag*rijmag))+(rik[1]/(rikmag*rikmag))));
      dcosjikdri[2] = ((rij[2]+rik[2])/(rijmag*rikmag)) -
        (cosjik*((rij[2]/(rijmag*rijmag))+(rik[2]/(rikmag*rikmag))));
      dcosjikdrk[0] = (-rij[0]/(rijmag*rikmag)) +
        (cosjik*(rik[0]/(rikmag*rikmag)));
      dcosjikdrk[1] = (-rij[1]/(rijmag*rikmag)) +
        (cosjik*(rik[1]/(rikmag*rikmag)));
      dcosjikdrk[2] = (-rij[2]/(rijmag*rikmag)) +
        (cosjik*(rik[2]/(rikmag*rikmag)));
      dcosjikdrj[0] = (-rik[0]/(rijmag*rikmag)) +
        (cosjik*(rij[0]/(rijmag*rijmag)));
      dcosjikdrj[1] = (-rik[1]/(rijmag*rikmag)) +
        (cosjik*(rij[1]/(rijmag*rijmag)));
      dcosjikdrj[2] = (-rik[2]/(rijmag*rikmag)) +
        (cosjik*(rij[2]/(rijmag*rijmag)));

      g = gSpline(cosjik,(NijC+NijH),itype,&dgdc,&dgdN);
      tmp2 = VA*.5*(tmp*wik*dgdc*exp(lamdajik));
      fj[0] = -tmp2*dcosjikdrj[0];
      fj[1] = -tmp2*dcosjikdrj[1];
      fj[2] = -tmp2*dcosjikdrj[2];
      fi[0] = -tmp2*dcosjikdri[0];
      fi[1] = -tmp2*dcosjikdri[1];
      fi[2] = -tmp2*dcosjikdri[2];
      fk[0] = -tmp2*dcosjikdrk[0];
      fk[1] = -tmp2*dcosjikdrk[1];
      fk[2] = -tmp2*dcosjikdrk[2];

      tmp2 = VA*.5*(tmp*wik*g*exp(lamdajik)*4.0*kronecker(itype,1));
      fj[0] -= tmp2*(-rij[0]/rijmag);
      fj[1] -= tmp2*(-rij[1]/rijmag);
      fj[2] -= tmp2*(-rij[2]/rijmag);
      fi[0] -= tmp2*((-rik[0]/rikmag)+(rij[0]/rijmag));
      fi[1] -= tmp2*((-rik[1]/rikmag)+(rij[1]/rijmag));
      fi[2] -= tmp2*((-rik[2]/rikmag)+(rij[2]/rijmag));
      fk[0] -= tmp2*(rik[0]/rikmag);
      fk[1] -= tmp2*(rik[1]/rikmag);
      fk[2] -= tmp2*(rik[2]/rikmag);

      // coordination forces

      // dwik forces

      tmp2 = VA*.5*(tmp*dwik*g*exp(lamdajik))/rikmag;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];

      // PIJ forces

      tmp2 = VA*.5*(tmp*dN2[ktype]*dwik)/rikmag;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];

      // dgdN forces

      tmp2 = VA*.5*(tmp*tmp3*dwik)/rikmag;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];

      f[atomi].x += fi[0]; f[atomi].y += fi[1]; f[atomi].z += fi[2];
      f[atomj].x += fj[0]; f[atomj].y += fj[1]; f[atomj].z += fj[2];
      f[atomk].x += fk[0]; f[atomk].y += fk[1]; f[atomk].z += fk[2];

      /*
      if (vflag_atom) {
        rji[0] = -rij[0]; rji[1] = -rij[1]; rji[2] = -rij[2];
        rki[0] = -rik[0]; rki[1] = -rik[1]; rki[2] = -rik[2];
        v_tally3(atomi,atomj,atomk,fj,fk,rji,rki);
	}
      */
    }
  }

  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  Etmp = 0.0;

  //REBO_neighs = REBO_firstneigh[j];
  
  for (l = 0; l < NBond[atomj]; l++) {
    atoml = Bonded[atomj*MaxNBond + l];
    if (atoml != atomi) {
      ltype = particle[atoml].type;
      rjl[0] = R_PBC(particle,atomj,atoml,0,L.x);
      rjl[1] = R_PBC(particle,atomj,atoml,1,L.y);
      rjl[2] = R_PBC(particle,atomj,atoml,2,L.z);
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      lamdaijl = 4.0*kronecker(jtype,1) *
        ((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dS);
      Nlj = nC[atoml]-(wjl*kronecker(jtype,0)) +
        nH[atoml]-(wjl*kronecker(jtype,1));
      cosijl = -1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2])) /
        (rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosijl,NjiC+NjiH,jtype,&dgdc,&dgdN);
      Etmp = Etmp+(wjl*g*exp(lamdaijl));
      //printf("Etmp=%.16f wjl=%.16f g=%.16f lamdaijl=%.16f\n",Etmp,wjl,g,lamdaijl);
      tmp3 = tmp3+(wjl*dgdN*exp(lamdaijl));
      NconjtmpJ = NconjtmpJ+(kronecker(ltype,0)*wjl*Sp(Nlj,Nmin,Nmax,dS));
    }
  }

  PjiS = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  //NjiC =1.0;
  //NjiH =1.0;
  PjiS = PijSpline(NjiC,NjiH,jtype,itype,dN2);
  //printf("PjiS=%.16f NjiC=%.16f NjiH=%.16f dN2[0]=%.16f dN2[1]=%.16f   jtype=%i itype=%i\n",PjiS,NjiC,NjiH,dN2[0],dN2[1],jtype,itype);
  pji = pow(1.0+Etmp+PjiS,-0.5);
  tmp = -0.5*cube(pji);
  //printf("pji=%.16f Etmp=%.16f PjiS=%.16f \n",pji,Etmp,PjiS);
  
  //REBO_neighs = REBO_firstneigh[j];
  for (l = 0; l < NBond[atomj]; l++) {
    atoml = Bonded[atomj*MaxNBond + l];
    if (atoml != atomi) {
      ltype = particle[atoml].type;
      
      rjl[0] = R_PBC(particle,atomj,atoml,0,L.x);
      rjl[1] = R_PBC(particle,atomj,atoml,1,L.y);
      rjl[2] = R_PBC(particle,atomj,atoml,2,L.z);
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      lamdaijl = 4.0*kronecker(jtype,1) *
        ((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
      cosijl = (-1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2]))) /
        (rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      dcosijldri[0] = (-rjl[0]/(rijmag*rjlmag)) -
        (cosijl*rij[0]/(rijmag*rijmag));
      dcosijldri[1] = (-rjl[1]/(rijmag*rjlmag)) -
        (cosijl*rij[1]/(rijmag*rijmag));
      dcosijldri[2] = (-rjl[2]/(rijmag*rjlmag)) -
        (cosijl*rij[2]/(rijmag*rijmag));
      dcosijldrj[0] = ((-rij[0]+rjl[0])/(rijmag*rjlmag)) +
        (cosijl*((rij[0]/square(rijmag))-(rjl[0]/(rjlmag*rjlmag))));
      dcosijldrj[1] = ((-rij[1]+rjl[1])/(rijmag*rjlmag)) +
        (cosijl*((rij[1]/square(rijmag))-(rjl[1]/(rjlmag*rjlmag))));
      dcosijldrj[2] = ((-rij[2]+rjl[2])/(rijmag*rjlmag)) +
        (cosijl*((rij[2]/square(rijmag))-(rjl[2]/(rjlmag*rjlmag))));
      dcosijldrl[0] = (rij[0]/(rijmag*rjlmag))+(cosijl*rjl[0]/(rjlmag*rjlmag));
      dcosijldrl[1] = (rij[1]/(rijmag*rjlmag))+(cosijl*rjl[1]/(rjlmag*rjlmag));
      dcosijldrl[2] = (rij[2]/(rijmag*rjlmag))+(cosijl*rjl[2]/(rjlmag*rjlmag));

      // evaluate splines g and derivatives dg

      g = gSpline(cosijl,NjiC+NjiH,jtype,&dgdc,&dgdN);
      tmp2 = VA*.5*(tmp*wjl*dgdc*exp(lamdaijl));
      fi[0] = -tmp2*dcosijldri[0];
      fi[1] = -tmp2*dcosijldri[1];
      fi[2] = -tmp2*dcosijldri[2];
      fj[0] = -tmp2*dcosijldrj[0];
      fj[1] = -tmp2*dcosijldrj[1];
      fj[2] = -tmp2*dcosijldrj[2];
      fl[0] = -tmp2*dcosijldrl[0];
      fl[1] = -tmp2*dcosijldrl[1];
      fl[2] = -tmp2*dcosijldrl[2];

      tmp2 = VA*.5*(tmp*wjl*g*exp(lamdaijl)*4.0*kronecker(jtype,1));
      fi[0] -= tmp2*(rij[0]/rijmag);
      fi[1] -= tmp2*(rij[1]/rijmag);
      fi[2] -= tmp2*(rij[2]/rijmag);
      fj[0] -= tmp2*((-rjl[0]/rjlmag)-(rij[0]/rijmag));
      fj[1] -= tmp2*((-rjl[1]/rjlmag)-(rij[1]/rijmag));
      fj[2] -= tmp2*((-rjl[2]/rjlmag)-(rij[2]/rijmag));
      fl[0] -= tmp2*(rjl[0]/rjlmag);
      fl[1] -= tmp2*(rjl[1]/rjlmag);
      fl[2] -= tmp2*(rjl[2]/rjlmag);

      // coordination forces

      // dwik forces

      tmp2 = VA*.5*(tmp*dwjl*g*exp(lamdaijl))/rjlmag;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];

      // PIJ forces

      tmp2 = VA*.5*(tmp*dN2[ltype]*dwjl)/rjlmag;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];

      // dgdN forces

      tmp2 = VA*.5*(tmp*tmp3*dwjl)/rjlmag;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];

      f[atomi].x += fi[0]; f[atomi].y += fi[1]; f[atomi].z += fi[2];
      f[atomj].x += fj[0]; f[atomj].y += fj[1]; f[atomj].z += fj[2];
      f[atoml].x += fl[0]; f[atoml].y += fl[1]; f[atoml].z += fl[2];

      /*
      if (vflag_atom) {
        rlj[0] = -rjl[0]; rlj[1] = -rjl[1]; rlj[2] = -rjl[2];
        v_tally3(atomi,atomj,atoml,fi,fl,rij,rlj);
	}
      */
    }
  }

  // evaluate Nij conj

  Nijconj = 1.0+(NconjtmpI*NconjtmpI)+(NconjtmpJ*NconjtmpJ);
  piRC = piRCSpline(NijC+NijH,NjiC+NjiH,Nijconj,itype,jtype,dN3);

  // piRC forces
  
  //REBO_neighs_i = REBO_firstneigh[i];
  for (k = 0; k < NBond[atomi]; k++) {
    atomk = Bonded[atomi*MaxNBond + k];
    if (atomk !=atomj) {
      ktype = particle[atomk].type;
      
      rik[0] = R_PBC(particle,atomi,atomk,0,L.x);
      rik[1] = R_PBC(particle,atomi,atomk,1,L.y);
      rik[2] = R_PBC(particle,atomi,atomk,2,L.z);
      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
      Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] -
        (wik*kronecker(itype,1));
      SpN = Sp(Nki,Nmin,Nmax,dNki);

      tmp2 = VA*dN3[0]*dwik/rikmag;

      f[atomi].x -= tmp2*rik[0];
      f[atomi].y -= tmp2*rik[1];
      f[atomi].z -= tmp2*rik[2];
      f[atomk].x += tmp2*rik[0];
      f[atomk].y += tmp2*rik[1];
      f[atomk].z += tmp2*rik[2];

      //if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

	tmp2 = VA*dN3[2]*(2.0*NconjtmpI*dwik*SpN)/rikmag;

      f[atomi].x -= tmp2*rik[0];
      f[atomi].y -= tmp2*rik[1];
      f[atomi].z -= tmp2*rik[2];
      f[atomk].x += tmp2*rik[0];
      f[atomk].y += tmp2*rik[1];
      f[atomk].z += tmp2*rik[2];

      //if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);
      atomk = Bonded[atomi*MaxNBond + k];
      if (fabs(dNki) > TOL) {
        //REBO_neighs_k = REBO_firstneigh[atomk];
        for (n = 0; n < NBond[atomk]; n++) {
          atomn = Bonded[atomk*MaxNBond + n];
          if (atomn != atomi) {
            ntype = particle[atomn].type;
	    rkn[0] = R_PBC(particle,atomk,atomn,0,L.x);
	    rkn[1] = R_PBC(particle,atomk,atomn,1,L.y);
	    rkn[2] = R_PBC(particle,atomk,atomn,2,L.z);
            rknmag = sqrt((rkn[0]*rkn[0])+(rkn[1]*rkn[1])+(rkn[2]*rkn[2]));
            Sp(rknmag,rcmin[ktype][ntype],rcmax[ktype][ntype],dwkn);

            tmp2 = VA*dN3[2]*(2.0*NconjtmpI*wik*dNki*dwkn)/rknmag;

            f[atomk].x -= tmp2*rkn[0];
            f[atomk].y -= tmp2*rkn[1];
            f[atomk].z -= tmp2*rkn[2];
            f[atomn].x += tmp2*rkn[0];
            f[atomn].y += tmp2*rkn[1];
            f[atomn].z += tmp2*rkn[2];

            //if (vflag_atom) v_tally2(atomk,atomn,-tmp2,rkn);
          }
        }
      }
    }
  }

  // piRC forces

  //REBO_neighs = REBO_firstneigh[atomj];
  for (l = 0; l < NBond[atomj]; l++) {
    atoml = Bonded[atomj*MaxNBond + l];
    if (atoml !=atomi) {
      ltype = particle[atoml].type;
     
      rjl[0] = R_PBC(particle,atomj,atoml,0,L.x);
      rjl[1] = R_PBC(particle,atomj,atoml,1,L.y);
      rjl[2] = R_PBC(particle,atomj,atoml,2,L.z);
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
      Nlj = nC[atoml]-(wjl*kronecker(jtype,0))+nH[atoml] -
        (wjl*kronecker(jtype,1));
      SpN = Sp(Nlj,Nmin,Nmax,dNlj);

      tmp2 = VA*dN3[1]*dwjl/rjlmag;

      f[atomj].x -= tmp2*rjl[0];
      f[atomj].y -= tmp2*rjl[1];
      f[atomj].z -= tmp2*rjl[2];
      f[atoml].x += tmp2*rjl[0];
      f[atoml].y += tmp2*rjl[1];
      f[atoml].z += tmp2*rjl[2];

      //if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

      tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*dwjl*SpN)/rjlmag;

      f[atomj].x -= tmp2*rjl[0];
      f[atomj].y -= tmp2*rjl[1];
      f[atomj].z -= tmp2*rjl[2];
      f[atoml].x += tmp2*rjl[0];
      f[atoml].y += tmp2*rjl[1];
      f[atoml].z += tmp2*rjl[2];

      //if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

      if (fabs(dNlj) > TOL) {
        //REBO_neighs_l = REBO_firstneigh[atoml];
        for (n = 0; n <NBond[atoml]; n++) {
	  atomn = Bonded[atoml*MaxNBond + n];
          //atomn = REBO_neighs_l[n];
          if (atomn != atomj) {
            ntype = particle[atomn].type;
	    rln[0] = R_PBC(particle,atoml,atomn,0,L.x);
	    rln[1] = R_PBC(particle,atoml,atomn,1,L.y);
	    rln[2] = R_PBC(particle,atoml,atomn,2,L.z);
            rlnmag = sqrt((rln[0]*rln[0])+(rln[1]*rln[1])+(rln[2]*rln[2]));
            Sp(rlnmag,rcmin[ltype][ntype],rcmax[ltype][ntype],dwln);

            tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*wjl*dNlj*dwln)/rlnmag;

            f[atoml].x -= tmp2*rln[0];
            f[atoml].y -= tmp2*rln[1];
            f[atoml].z -= tmp2*rln[2];
            f[atomn].x += tmp2*rln[0];
            f[atomn].y += tmp2*rln[1];
            f[atomn].z += tmp2*rln[2];
	  
            //if (vflag_atom) v_tally2(atoml,atomn,-tmp2,rln);
          }
        }
      }
    }
  }

  Tij = 0.0;
  dN3[0] = 0.0;
  dN3[1] = 0.0;
  dN3[2] = 0.0;
  if (itype == 0 && jtype == 0)
    Tij=TijSpline((NijC+NijH),(NjiC+NjiH),Nijconj,dN3);
  Etmp = 0.0;

  if (fabs(Tij) > TOL) {
    atom2 = atomi;
    atom3 = atomj;
    r32[0] = R_PBC(particle,atom3,atom2,0,L.x);
    r32[1] = R_PBC(particle,atom3,atom2,1,L.y);
    r32[2] = R_PBC(particle,atom3,atom2,2,L.z);
    r32mag = sqrt((r32[0]*r32[0])+(r32[1]*r32[1])+(r32[2]*r32[2]));
    r23[0] = -r32[0];
    r23[1] = -r32[1];
    r23[2] = -r32[2];
    r23mag = r32mag;
    //REBO_neighs_i = REBO_firstneigh[i];
    for (k = 0; k < NBond[atomi]; k++) {
      atomk = Bonded[atomi*MaxNBond + k];
      atom1 = atomk;
      ktype = particle[atomk].type;
      if (atomk != atomj) {
	r21[0] = R_PBC(particle,atom2,atom1,0,L.x);
	r21[1] = R_PBC(particle,atom2,atom1,1,L.y);
	r21[2] = R_PBC(particle,atom2,atom1,2,L.z);
        r21mag = sqrt(r21[0]*r21[0] + r21[1]*r21[1] + r21[2]*r21[2]);
        cos321 = -1.0*((r21[0]*r32[0])+(r21[1]*r32[1])+(r21[2]*r32[2])) /
          (r21mag*r32mag);
        cos321 = MIN(cos321,1.0);
        cos321 = MAX(cos321,-1.0);
        Sp2(cos321,thmin,thmax,dcut321);
        sin321 = sqrt(1.0 - cos321*cos321);
        sink2i = 1.0/(sin321*sin321);
        rik2i = 1.0/(r21mag*r21mag);
        if (sin321 != 0.0) {
          rr = (r23mag*r23mag)-(r21mag*r21mag);
          //rjk[0] = r21[0]-r23[0];
          //rjk[1] = r21[1]-r23[1];
          //rjk[2] = r21[2]-r23[2];

	  rjk[0] = R_PBC(particle,atomj,atomk,0,L.x);
          rjk[1] = R_PBC(particle,atomj,atomk,1,L.y);
          rjk[2] = R_PBC(particle,atomj,atomk,2,L.z);
	  
          rjk2 = (rjk[0]*rjk[0])+(rjk[1]*rjk[1])+(rjk[2]*rjk[2]);
          rijrik = 2.0*r23mag*r21mag;
          rik2 = r21mag*r21mag;
          dctik = (-rr+rjk2)/(rijrik*rik2);
          dctij = (rr+rjk2)/(rijrik*r23mag*r23mag);
          dctjk = -2.0/rijrik;
          w21 = Sp(r21mag,rcmin[itype][ktype],rcmaxp[itype][ktype],dw21);
          rijmag = r32mag;
          rikmag = r21mag;
          rij2 = r32mag*r32mag;
          rik2 = r21mag*r21mag;
          costmp = 0.5*(rij2+rik2-rjk2)/rijmag/rikmag;
          tspjik = Sp2(costmp,thmin,thmax,dtsjik);
          dtsjik = -dtsjik;

          //REBO_neighs_j = REBO_firstneigh[j];
          for (l = 0; l < NBond[atomj]; l++) {
            //atoml = REBO_neighs_j[l];
	    atoml = Bonded[atomj*MaxNBond + l];
            atom4 = atoml;
            ltype = particle[atoml].type;
            if (!(atoml == atomi || atoml == atomk)) {
	      r34[0] = R_PBC(particle,atom3,atom4,0,L.x);
	      r34[1] = R_PBC(particle,atom3,atom4,1,L.y);
	      r34[2] = R_PBC(particle,atom3,atom4,2,L.z);
              r34mag = sqrt((r34[0]*r34[0])+(r34[1]*r34[1])+(r34[2]*r34[2]));
              cos234 = (r32[0]*r34[0] + r32[1]*r34[1] + r32[2]*r34[2]) /
                (r32mag*r34mag);
              cos234 = MIN(cos234,1.0);
              cos234 = MAX(cos234,-1.0);
              sin234 = sqrt(1.0 - cos234*cos234);
              sinl2i = 1.0/(sin234*sin234);
              rjl2i = 1.0/(r34mag*r34mag);

              if (sin234 != 0.0) {
                w34 = Sp(r34mag,rcmin[jtype][ltype],rcmaxp[jtype][ltype],dw34);
                rr = (r23mag*r23mag)-(r34mag*r34mag);
                //ril[0] = r23[0]+r34[0];
                //ril[1] = r23[1]+r34[1];
                //ril[2] = r23[2]+r34[2];

		ril[0] = R_PBC(particle,atomi,atoml,0,L.x);
                ril[1] = R_PBC(particle,atomi,atoml,1,L.y);
                ril[2] = R_PBC(particle,atomi,atoml,2,L.z);
		
                ril2 = (ril[0]*ril[0])+(ril[1]*ril[1])+(ril[2]*ril[2]);
                rijrjl = 2.0*r23mag*r34mag;
                rjl2 = r34mag*r34mag;
                dctjl = (-rr+ril2)/(rijrjl*rjl2);
                dctji = (rr+ril2)/(rijrjl*r23mag*r23mag);
                dctil = -2.0/rijrjl;
                rjlmag = r34mag;
                rjl2 = r34mag*r34mag;
                costmp = 0.5*(rij2+rjl2-ril2)/rijmag/rjlmag;
                tspijl = Sp2(costmp,thmin,thmax,dtsijl);
                dtsijl = -dtsijl;
                prefactor = VA*Tij;

                cross321[0] = (r32[1]*r21[2])-(r32[2]*r21[1]);
                cross321[1] = (r32[2]*r21[0])-(r32[0]*r21[2]);
                cross321[2] = (r32[0]*r21[1])-(r32[1]*r21[0]);
                cross234[0] = (r23[1]*r34[2])-(r23[2]*r34[1]);
                cross234[1] = (r23[2]*r34[0])-(r23[0]*r34[2]);
                cross234[2] = (r23[0]*r34[1])-(r23[1]*r34[0]);

                cwnum = (cross321[0]*cross234[0]) +
                  (cross321[1]*cross234[1]) + (cross321[2]*cross234[2]);
                cwnom = r21mag*r34mag*r23mag*r23mag*sin321*sin234;
                om1234 = cwnum/cwnom;
                cw = om1234;

		Etmp += ( (1.0 - sqr(om1234) )*w21*w34)*(1.0-tspjik)*(1.0-tspijl);
		
                dt1dik = (rik2i)-(dctik*sink2i*cos321);
                dt1djk = (-dctjk*sink2i*cos321);
                dt1djl = (rjl2i)-(dctjl*sinl2i*cos234);
                dt1dil = (-dctil*sinl2i*cos234);
                dt1dij = (2.0/(r23mag*r23mag))-(dctij*sink2i*cos321) -
                  (dctji*sinl2i*cos234);

                dt2dik[0] = (-r23[2]*cross234[1])+(r23[1]*cross234[2]);
                dt2dik[1] = (-r23[0]*cross234[2])+(r23[2]*cross234[0]);
                dt2dik[2] = (-r23[1]*cross234[0])+(r23[0]*cross234[1]);

                dt2djl[0] = (-r23[1]*cross321[2])+(r23[2]*cross321[1]);
                dt2djl[1] = (-r23[2]*cross321[0])+(r23[0]*cross321[2]);
                dt2djl[2] = (-r23[0]*cross321[1])+(r23[1]*cross321[0]);

                dt2dij[0] = (r21[2]*cross234[1])-(r34[2]*cross321[1]) -
                  (r21[1]*cross234[2])+(r34[1]*cross321[2]);
                dt2dij[1] = (r21[0]*cross234[2])-(r34[0]*cross321[2]) -
                  (r21[2]*cross234[0])+(r34[2]*cross321[0]);
                dt2dij[2] = (r21[1]*cross234[0])-(r34[1]*cross321[0]) -
                  (r21[0]*cross234[1])+(r34[0]*cross321[1]);

                aa = (prefactor*2.0*cw/cwnom)*w21*w34 *
                  (1.0-tspjik)*(1.0-tspijl);
                aaa1 = -prefactor*(1.0-square(om1234)) *
                  (1.0-tspjik)*(1.0-tspijl);
                aaa2 = aaa1*w21*w34;
                at2 = aa*cwnum;

                fcijpc = (-dt1dij*at2)+(aaa2*dtsjik*dctij*(1.0-tspijl)) +
                  (aaa2*dtsijl*dctji*(1.0-tspjik));
                fcikpc = (-dt1dik*at2)+(aaa2*dtsjik*dctik*(1.0-tspijl));
                fcjlpc = (-dt1djl*at2)+(aaa2*dtsijl*dctjl*(1.0-tspjik));
                fcjkpc = (-dt1djk*at2)+(aaa2*dtsjik*dctjk*(1.0-tspijl));
                fcilpc = (-dt1dil*at2)+(aaa2*dtsijl*dctil*(1.0-tspjik));

                F23[0] = (fcijpc*r23[0])+(aa*dt2dij[0]);
                F23[1] = (fcijpc*r23[1])+(aa*dt2dij[1]);
                F23[2] = (fcijpc*r23[2])+(aa*dt2dij[2]);

                F12[0] = (fcikpc*r21[0])+(aa*dt2dik[0]);
                F12[1] = (fcikpc*r21[1])+(aa*dt2dik[1]);
                F12[2] = (fcikpc*r21[2])+(aa*dt2dik[2]);

                F34[0] = (fcjlpc*r34[0])+(aa*dt2djl[0]);
                F34[1] = (fcjlpc*r34[1])+(aa*dt2djl[1]);
                F34[2] = (fcjlpc*r34[2])+(aa*dt2djl[2]);

                F31[0] = (fcjkpc*rjk[0]);
                F31[1] = (fcjkpc*rjk[1]);
                F31[2] = (fcjkpc*rjk[2]);

                F24[0] = (fcilpc*ril[0]);
                F24[1] = (fcilpc*ril[1]);
                F24[2] = (fcilpc*ril[2]);

                f1[0] = -F12[0]-F31[0];
                f1[1] = -F12[1]-F31[1];
                f1[2] = -F12[2]-F31[2];
                f2[0] = F23[0]+F12[0]+F24[0];
                f2[1] = F23[1]+F12[1]+F24[1];
                f2[2] = F23[2]+F12[2]+F24[2];
                f3[0] = -F23[0]+F34[0]+F31[0];
                f3[1] = -F23[1]+F34[1]+F31[1];
                f3[2] = -F23[2]+F34[2]+F31[2];
                f4[0] = -F34[0]-F24[0];
                f4[1] = -F34[1]-F24[1];
                f4[2] = -F34[2]-F24[2];

                // coordination forces

                tmp2 = VA*Tij*((1.0-(om1234*om1234))) *
                  (1.0-tspjik)*(1.0-tspijl)*dw21*w34/r21mag;
                f2[0] -= tmp2*r21[0];
                f2[1] -= tmp2*r21[1];
                f2[2] -= tmp2*r21[2];
                f1[0] += tmp2*r21[0];
                f1[1] += tmp2*r21[1];
                f1[2] += tmp2*r21[2];

                tmp2 = VA*Tij*((1.0-(om1234*om1234))) *
                  (1.0-tspjik)*(1.0-tspijl)*w21*dw34/r34mag;
                f3[0] -= tmp2*r34[0];
                f3[1] -= tmp2*r34[1];
                f3[2] -= tmp2*r34[2];
                f4[0] += tmp2*r34[0];
                f4[1] += tmp2*r34[1];
                f4[2] += tmp2*r34[2];

                f[atom1].x += f1[0]; f[atom1].y += f1[1];
                f[atom1].z += f1[2];
                f[atom2].x += f2[0]; f[atom2].y += f2[1];
                f[atom2].z += f2[2];
                f[atom3].x += f3[0]; f[atom3].y += f3[1];
                f[atom3].z += f3[2];
                f[atom4].x += f4[0]; f[atom4].y += f4[1];
                f[atom4].z += f4[2];

		/*
                if (vflag_atom) {
                  r13[0] = -rjk[0]; r13[1] = -rjk[1]; r13[2] = -rjk[2];
                  r43[0] = -r34[0]; r43[1] = -r34[1]; r43[2] = -r34[2];
                  v_tally4(atom1,atom2,atom3,atom4,f1,f2,f4,r13,r23,r43);
		  }
		*/
              }
            }
          }
        }
      }
    }

    // Tij forces now that we have Etmp

    //REBO_neighs = REBO_firstneigh[i];
    for (k = 0; k < NBond[atomi]; k++) {
      atomk = Bonded[atomi*MaxNBond + k];
      //atomk = REBO_neighs[k];
      if (atomk != atomj) {
        ktype = particle[atomk].type;
	rik[0] = R_PBC(particle,atomi,atomk,0,L.x);
	rik[1] = R_PBC(particle,atomi,atomk,1,L.y);
	rik[2] = R_PBC(particle,atomi,atomk,2,L.z);
        rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
        wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
        Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] -
          (wik*kronecker(itype,1));
        SpN = Sp(Nki,Nmin,Nmax,dNki);

        tmp2 = VA*dN3[0]*dwik*Etmp/rikmag;

        f[atomi].x -= tmp2*rik[0];
        f[atomi].y -= tmp2*rik[1];
        f[atomi].z -= tmp2*rik[2];
        f[atomk].x += tmp2*rik[0];
        f[atomk].y += tmp2*rik[1];
        f[atomk].z += tmp2*rik[2];

        //if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

        tmp2 = VA*dN3[2]*(2.0*NconjtmpI*dwik*SpN)*Etmp/rikmag;

        f[atomi].x -= tmp2*rik[0];
        f[atomi].y -= tmp2*rik[1];
        f[atomi].z -= tmp2*rik[2];
        f[atomk].x += tmp2*rik[0];
        f[atomk].y += tmp2*rik[1];
        f[atomk].z += tmp2*rik[2];

        //if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

        if (fabs(dNki) > TOL) {
          //REBO_neighs_k = REBO_firstneigh[atomk];
          for (n = 0; n < NBond[atomk]; n++) {
            atomn = Bonded[atomk*MaxNBond + n];
            ntype = particle[atomn].type;
            if (atomn != atomi) {
	      rkn[0] = R_PBC(particle,atomk,atomn,0,L.x);
	      rkn[1] = R_PBC(particle,atomk,atomn,1,L.y);
	      rkn[2] = R_PBC(particle,atomk,atomn,2,L.z);
              rknmag = sqrt((rkn[0]*rkn[0])+(rkn[1]*rkn[1])+(rkn[2]*rkn[2]));
              Sp(rknmag,rcmin[ktype][ntype],rcmax[ktype][ntype],dwkn);

              tmp2 = VA*dN3[2]*(2.0*NconjtmpI*wik*dNki*dwkn)*Etmp/rknmag;

              f[atomk].x -= tmp2*rkn[0];
              f[atomk].y -= tmp2*rkn[1];
              f[atomk].z -= tmp2*rkn[2];
              f[atomn].x += tmp2*rkn[0];
              f[atomn].y += tmp2*rkn[1];
              f[atomn].z += tmp2*rkn[2];

              //if (vflag_atom) v_tally2(atomk,atomn,-tmp2,rkn);
            }
          }
        }
      }
    }

    // Tij forces

    //REBO_neighs = REBO_firstneigh[j];
    for (l = 0; l < NBond[atomj]; l++) {
      atoml = Bonded[atomj*MaxNBond + l];
      if (atoml != atomi) {
        ltype = particle[atoml].type;
      	rjl[0] = R_PBC(particle,atomj,atoml,0,L.x);
	rjl[1] = R_PBC(particle,atomj,atoml,1,L.y);
	rjl[2] = R_PBC(particle,atomj,atoml,2,L.z);
        rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
        wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
        Nlj = nC[atoml]-(wjl*kronecker(jtype,0))+nH[atoml] -
          (wjl*kronecker(jtype,1));
        SpN = Sp(Nlj,Nmin,Nmax,dNlj);

        tmp2 = VA*dN3[1]*dwjl*Etmp/rjlmag;

        f[atomj].x -= tmp2*rjl[0];
        f[atomj].y -= tmp2*rjl[1];
        f[atomj].z -= tmp2*rjl[2];
        f[atoml].x += tmp2*rjl[0];
        f[atoml].y += tmp2*rjl[1];
        f[atoml].z += tmp2*rjl[2];

        //if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

        tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*dwjl*SpN)*Etmp/rjlmag;

        f[atomj].x -= tmp2*rjl[0];
        f[atomj].y -= tmp2*rjl[1];
        f[atomj].z -= tmp2*rjl[2];
        f[atoml].x += tmp2*rjl[0];
        f[atoml].y += tmp2*rjl[1];
        f[atoml].z += tmp2*rjl[2];
	
        //if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

        if (fabs(dNlj) > TOL) {
          //REBO_neighs_l = REBO_firstneigh[atoml];
          for (n = 0; n < NBond[atoml]; n++) {
	    atomn = Bonded[atoml*MaxNBond + n];
            //atomn = REBO_neighs_l[n];
            ntype = particle[atomn].type;
            if (atomn !=atomj) {
	      rln[0] = R_PBC(particle,atoml,atomn,0,L.x);
	      rln[1] = R_PBC(particle,atoml,atomn,1,L.y);
	      rln[2] = R_PBC(particle,atoml,atomn,2,L.z);
              rlnmag = sqrt((rln[0]*rln[0])+(rln[1]*rln[1])+(rln[2]*rln[2]));
              Sp(rlnmag,rcmin[ltype][ntype],rcmax[ltype][ntype],dwln);
	      
              tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*wjl*dNlj*dwln)*Etmp/rlnmag;
	
	      f[atoml].x -= tmp2*rln[0];
              f[atoml].y -= tmp2*rln[1];
              f[atoml].z -= tmp2*rln[2];
              f[atomn].x += tmp2*rln[0];
              f[atomn].y += tmp2*rln[1];
              f[atomn].z += tmp2*rln[2];
	
              //if (vflag_atom) v_tally2(atoml,atomn,-tmp2,rln);
            }
          }
        }
      }
    }
  }

  bij = (0.5*(pij+pji))+piRC+(Tij*Etmp);
  //  if(atomi==3434 && atomj ==3433)printf("bij=%.16f pij=%.16f pji=%.16f piRC=%.16f Tij=%.16f Etmp=%.16f\n",bij,pij,pji,piRC,Tij,Etmp);
  return bij;
 
 }


double gSpline(double costh, double Nij, int typei, double *dgdc, double *dgdN)
{
  double coeffs[6],dS,g1,g2,dg1,dg2,cut,g;
  int i,j;

  i = 0;
  j = 0;
  g = 0.0;
  cut = 0.0;
  dS = 0.0;
  dg1 = 0.0;
    dg2 = 0.0;
   *dgdc = 0.0;
   *dgdN = 0.0;

  // central atom is Carbon

   if (typei == 0) {
    if (costh < gCdom[0]) costh = gCdom[0]; //gCdom[i] parameter
    if (costh > gCdom[4]) costh = gCdom[4];
    if (Nij >= NCmax) {
      for (i = 0; i < 4; i++) {
        if (costh >= gCdom[i] && costh <= gCdom[i+1])  {
          for (j = 0; j < 6; j++) coeffs[j] = gC2[i][j]; //gC2[i] parameter
         } 
      }
      g2 = Sp5th(costh,coeffs,&dg2);
      g = g2;
      *dgdc = dg2;
      *dgdN = 0.0;
    }
    if (Nij <= NCmin) {                        //NCmin parameter
      for (i = 0; i < 4; i++) {
        if (costh >= gCdom[i] && costh <= gCdom[i+1]) {
          for (j = 0; j < 6; j++) coeffs[j] = gC1[i][j];
        }
      }
      g1 = Sp5th(costh,coeffs,&dg1);
      g = g1;
      *dgdc = dg1;
      *dgdN = 0.0;
    }
    if (Nij > NCmin && Nij < NCmax) { //NCmax parameter
      for (i = 0; i < 4; i++) {
        if (costh >= gCdom[i] && costh <= gCdom[i+1]) {
          for (j = 0; j < 6; j++) coeffs[j] = gC1[i][j];
        }
      }
      g1 = Sp5th(costh,coeffs,&dg1);
      for (i = 0; i < 4; i++) {
        if (costh >= gCdom[i] && costh <= gCdom[i+1]) {
          for (j = 0; j < 6; j++) coeffs[j] = gC2[i][j];
        }
      }
      g2 = Sp5th(costh,coeffs,&dg2);
      cut = Sp(Nij,NCmin,NCmax,dS);
      g = g2+cut*(g1-g2);
      *dgdc = dg2+(cut*(dg1-dg2));
      *dgdN = dS*(g1-g2);
    }
  }

  // central atom is Hydrogen

  if (typei == 1) {
    if (costh < gHdom[0]) costh = gHdom[0]; //gHdom parameter
    if (costh > gHdom[3]) costh = gHdom[3];
    for (i = 0; i < 3; i++) {
      if (costh >= gHdom[i] && costh <= gHdom[i+1]) {
        for (j = 0; j < 6; j++) coeffs[j] = gH[i][j];
      }
    }
    g = Sp5th(costh,coeffs,&dg1);
    *dgdN = 0.0;
    *dgdc = dg1;
  }

  return g;
}

double PijSpline(double NijC, double NijH, int typei, int typej,double dN2[2])
{
  int x,y,i,done;
  double Pij,coeffs[16];
  
  for (i = 0; i < 16; i++) coeffs[i]=0.0;
  
  x = 0;
  y = 0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  done = 0;

  // if inputs are out of bounds set them back to a point in bounds

  if (typei == 0 && typej == 0) {
    if (NijC < pCCdom[0][0]) NijC=pCCdom[0][0];
    if (NijC > pCCdom[0][1]) NijC=pCCdom[0][1];
    if (NijH < pCCdom[1][0]) NijH=pCCdom[1][0];
    if (NijH > pCCdom[1][1]) NijH=pCCdom[1][1];
    //printf("NijC=%.16f NijH=%.16f\n",NijC,NijH);
    if (fabs(NijC-floor(NijC)) < TOL && fabs(NijH-floor(NijH)) < TOL) {
      Pij = PCCf[(int) NijC][(int) NijH];
      dN2[0] = PCCdfdx[(int) NijC][(int) NijH];
      dN2[1] = PCCdfdy[(int) NijC][(int) NijH];
      done = 1;
      //printf("Pij=%.16f NijC=%.16f NijH=%.16f dN[0]=%.16f dN[1]=%.16f \n",Pij,NijC,NijH,dN2[0],dN2[1]);
    }
    if (done == 0) {
      x = (int) (floor(NijC));
      y = (int) (floor(NijH));
      for (i = 0; i<16; i++) coeffs[i] = pCC[x][y][i];
      Pij = Spbicubic(NijC,NijH,coeffs,dN2);
      //printf("Pij=%.16f NijC=%.16f NijH=%.16f coeffs=%.16f dN2[0]==%.16f dN2[1]=%.16f \n",Pij,NijC,NijH,coeffs,dN2[0],dN2[1]);
    }
  }

  // if inputs are out of bounds set them back to a point in bounds

   if (typei == 0 && typej == 1){
     if (NijC < pCHdom[0][0]) NijC=pCHdom[0][0];
     if (NijC > pCHdom[0][1]) NijC=pCHdom[0][1];
      if (NijH < pCHdom[1][0]) NijH=pCHdom[1][0];
      if (NijH > pCHdom[1][1]) NijH=pCHdom[1][1];

    if (fabs(NijC-floor(NijC)) < TOL && fabs(NijH-floor(NijH)) < TOL) {
      Pij = PCHf[(int) NijC][(int) NijH];
      dN2[0] = PCHdfdx[(int) NijC][(int) NijH];
      dN2[1] = PCHdfdy[(int) NijC][(int) NijH];
      done = 1;
    }
    if (done == 0) {
      x = (int) (floor(NijC));
      y = (int) (floor(NijH));
      for (i = 0; i<16; i++) coeffs[i] = pCH[x][y][i];
      Pij = Spbicubic(NijC,NijH,coeffs,dN2);
    }
  }

  if (typei == 1 && typej == 0) {
    Pij = 0.0;
    dN2[0] = 0.0;
    dN2[1] = 0.0;
  }


  if (typei == 1 && typej == 1) {
    Pij = 0.0;
    dN2[0] = 0.0;
    dN2[1] = 0.0;
  }
  return Pij;
}
double Sp5th(double x, double coeffs[6], double *df)
{
  double f, d;
  const double x2 = x*x;
  const double x3 = x2*x;
  
  f  = coeffs[0];
  f += coeffs[1]*x;
  d  = coeffs[1];
  f += coeffs[2]*x2;
  d += 2.0*coeffs[2]*x;
  f += coeffs[3]*x3;
  d += 3.0*coeffs[3]*x2;
  f += coeffs[4]*x2*x2;
  d += 4.0*coeffs[4]*x3;
  f += coeffs[5]*x2*x3;
  d += 5.0*coeffs[5]*x2*x2;
  
  *df = d;
  return f;
}
double Spbicubic(double x, double y,double coeffs[16], double df[2])
{
  double f,xn,yn,xn1,yn1,c;
  int i,j;

  f = 0.0;
  df[0] = 0.0;
  df[1] = 0.0;

  xn = 1.0;
  for (i = 0; i < 4; i++) {
    yn = 1.0;
    for (j = 0; j < 4; j++) {
      c = coeffs[i*4+j];

      f += c*xn*yn;
      if (i > 0) df[0] += c * ((double) i) * xn1 * yn;
      if (j > 0) df[1] += c * ((double) j) * xn * yn1;

      yn1 = yn;
      yn *= y;
    }
    xn1 = xn;
    xn *= x;
  }

  return f;
}

double piRCSpline(double Nij, double Nji, double Nijconj,int typei, int typej, double dN3[3])
{
  int x,y,z,i,done;
  double piRC,coeffs[64];
  x=0;
  y=0;
  z=0;
  i=0;

  done=0;

  for (i=0; i<64; i++) coeffs[i]=0.0;

  if (typei==0 && typej==0) {
    //if the inputs are out of bounds set them back to a point in bounds
    if (Nij<piCCdom[0][0]) Nij=piCCdom[0][0];
    if (Nij>piCCdom[0][1]) Nij=piCCdom[0][1];
    if (Nji<piCCdom[1][0]) Nji=piCCdom[1][0];
    if (Nji>piCCdom[1][1]) Nji=piCCdom[1][1];
    if (Nijconj<piCCdom[2][0]) Nijconj=piCCdom[2][0];
    if (Nijconj>piCCdom[2][1]) Nijconj=piCCdom[2][1];

    if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL &&
        fabs(Nijconj-floor(Nijconj))<TOL) {
      piRC=piCCf[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[0]=piCCdfdx[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[1]=piCCdfdy[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[2]=piCCdfdz[(int) Nij][(int) Nji][(int) Nijconj];
      done=1;
    }

    if (done==0) {
      for (i=0; i<piCCdom[0][1]; i++)
        if (Nij>=(double) i && Nij<=(double) i+1 || Nij==(double) i) x=i;
      for (i=0; i<piCCdom[1][1]; i++)
        if (Nji>=(double) i && Nji<=(double) i+1 || Nji==(double) i) y=i;
      for (i=0; i<piCCdom[2][1]; i++)
        if (Nijconj>=(double) i && Nijconj<=(double) i+1 ||
            Nijconj==(double) i) z=i;

      for (i=0; i<64; i++) coeffs[i]=piCC[x][y][z][i];
      piRC=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
    }
  }


  // CH interaction

  if (typei==0 && typej==1 || typei==1 && typej==0) {
    // if the inputs are out of bounds set them back to a point in bounds

    if (Nij<piCHdom[0][0] || Nij>piCHdom[0][1] ||
        Nji<piCHdom[1][0] || Nji>piCHdom[1][1] ||
        Nijconj<piCHdom[2][0] || Nijconj>piCHdom[2][1]) {
      if (Nij<piCHdom[0][0]) Nij=piCHdom[0][0];
      if (Nij>piCHdom[0][1]) Nij=piCHdom[0][1];
      if (Nji<piCHdom[1][0]) Nji=piCHdom[1][0];
      if (Nji>piCHdom[1][1]) Nji=piCHdom[1][1];
      if (Nijconj<piCHdom[2][0]) Nijconj=piCHdom[2][0];
      if (Nijconj>piCHdom[2][1]) Nijconj=piCHdom[2][1];
    }

    if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL &&
        fabs(Nijconj-floor(Nijconj))<TOL) {
      piRC=piCHf[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[0]=piCHdfdx[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[1]=piCHdfdy[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[2]=piCHdfdz[(int) Nij][(int) Nji][(int) Nijconj];
      done=1;
    }

    if (done==0) {
      for (i=0; i<piCHdom[0][1]; i++)
        if (Nij>=i && Nij<=i+1) x=i;
      for (i=0; i<piCHdom[1][1]; i++)
        if (Nji>=i && Nji<=i+1) y=i;
      for (i=0; i<piCHdom[2][1]; i++)
        if (Nijconj>=i && Nijconj<=i+1) z=i;

      for (i=0; i<64; i++) coeffs[i]=piCH[x][y][z][i];
      piRC=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
    }
  }

  if (typei==1 && typej==1) {
    if (Nij<piHHdom[0][0] || Nij>piHHdom[0][1] ||
        Nji<piHHdom[1][0] || Nji>piHHdom[1][1] ||
        Nijconj<piHHdom[2][0] || Nijconj>piHHdom[2][1]) {
      Nij=0.0;
      Nji=0.0;
      Nijconj=0.0;
    }
    if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL &&
        fabs(Nijconj-floor(Nijconj))<TOL) {
      piRC=piHHf[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[0]=piHHdfdx[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[1]=piHHdfdy[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[2]=piHHdfdz[(int) Nij][(int) Nji][(int) Nijconj];
      done=1;
    }
    if (done==0) {
      for (i=0; i<piHHdom[0][1]; i++)
        if (Nij>=i && Nij<=i+1) x=i;
      for (i=0; i<piHHdom[1][1]; i++)
        if (Nji>=i && Nji<=i+1) y=i;
      for (i=0; i<piHHdom[2][1]; i++)
        if (Nijconj>=i && Nijconj<=i+1) z=i;

      for (i=0; i<64; i++) coeffs[i]=piHH[x][y][z][i];
      piRC=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
    }
  }

  return piRC;
}

double TijSpline(double Nij, double Nji,double Nijconj, double dN3[3])
{
  int x,y,z,i,done;
  double Tijf,coeffs[64];

  x=0;
  y=0;
  z=0;
  i=0;
  Tijf=0.0;
  done=0;
  for (i=0; i<64; i++) coeffs[i]=0.0;

  //if the inputs are out of bounds set them back to a point in bounds

  if (Nij<Tijdom[0][0]) Nij=Tijdom[0][0];
  if (Nij>Tijdom[0][1]) Nij=Tijdom[0][1];
  if (Nji<Tijdom[1][0]) Nji=Tijdom[1][0];
  if (Nji>Tijdom[1][1]) Nji=Tijdom[1][1];
  if (Nijconj<Tijdom[2][0]) Nijconj=Tijdom[2][0];
  if (Nijconj>Tijdom[2][1]) Nijconj=Tijdom[2][1];

  if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL &&
      fabs(Nijconj-floor(Nijconj))<TOL) {
    Tijf=Tf[(int) Nij][(int) Nji][(int) Nijconj];
    dN3[0]=Tdfdx[(int) Nij][(int) Nji][(int) Nijconj];
    dN3[1]=Tdfdy[(int) Nij][(int) Nji][(int) Nijconj];
    dN3[2]=Tdfdz[(int) Nij][(int) Nji][(int) Nijconj];
    done=1;
  }

  if (done==0) {
    for (i=0; i<Tijdom[0][1]; i++)
      if (Nij>=i && Nij<=i+1) x=i;
    for (i=0; i<Tijdom[1][1]; i++)
      if (Nji>=i && Nji<=i+1) y=i;
    for (i=0; i<Tijdom[2][1]; i++)
      if (Nijconj>=i && Nijconj<=i+1) z=i;

    for (i=0; i<64; i++) coeffs[i]=Tijc[x][y][z][i];
    Tijf=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
  }

  return Tijf;
}

double Sptricubic(double x, double y, double z,
                              double coeffs[64], double df[3])
{
  double f,ir,jr,kr,xn,yn,zn,xn1,yn1,zn1,c;
  int i,j,k;

  f = 0.0;
  df[0] = 0.0;
  df[1] = 0.0;
  df[2] = 0.0;

  xn = 1.0;
  for (i = 0; i < 4; i++) {
    ir = (double) i;
    yn = 1.0;
    for (j = 0; j < 4; j++) {
      jr = (double) j;
      zn = 1.0;
      for (k = 0; k < 4; k++) {
        kr = (double) k;
        c = coeffs[16*i+4*j+k];
        f += c*xn*yn*zn;
        if (i > 0) df[0] += c * ir * xn1 * yn * zn;
        if (j > 0) df[1] += c * jr * xn * yn1 * zn;
        if (k > 0) df[2] += c * kr * xn * yn * zn1;
        zn1 = zn;
        zn *= z;
      }
      yn1 = yn;
      yn *= y;
    }
    xn1 = xn;
    xn *= x;
  }

  return f;
}

void calc_force_arrays(int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *angle, double *Fc,int *Neighb,int *NNeighb,double *rij_long,int *NList,int *List,point *R0,int *NBond_once,int *Bonded_once,int potential_flag)
{
  Check_List(particle,R0,Natoms,List,NList,L);
 
  //calc_rij(NBond,Bonded,BondNum,Natoms,particle,L,rij,Fc,BondParams,Neighb,NNeighb,rij_long,NList,List,NBond_once,Bonded_once,potential_flag);
  calc_rij_long(Natoms,particle,L,Neighb,NNeighb);
   
  calc_rij_short(NBond,Bonded,BondNum,Natoms,particle,L,rij,Fc,BondParams,Neighb,NNeighb,rij_long,NList,List,NBond_once,Bonded_once,potential_flag);
 
  //calc_angles(angle,Natoms,NBond,Bonded,rij,particle,L);
}

/****************************************************************/
/* Calculate the total force array                              */

void Calc_Force(point *Fterm, int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rij_long, FILE* logfile,double *angle,Normal_struct *Normal,int *Normal_atom,int Interlayer,double *Fc,int *Neighb,int *NNeighb,int potential_flag,double &E_tot,int *NBond_once,int *Bonded_once,int BN_flag,int Ncar_hyd,int N_init,int N_init_ILP)
{
  E_tot     = 0.0;
  
#pragma omp parallel
  {
    int atomi;
    double E_TORSION,E_KC,E_REBO,E_Tersoff,Ecoulomb;
    double EvdW,ERep;
    clock_t start , diff;
    int msec;
    //struct timeval t1, t2;
    int elapsedTime;
    point *Fterm_;
    
    Fterm_ = new point[Natoms];
    
    
    dr_ij *dN = new dr_ij[1];
    
    E_TORSION = 0.0; 
    E_REBO    = 0.0;
    E_KC      = 0.0;
    
    E_Tersoff = 0.0;
    EvdW      = 0.0;
    ERep      = 0.0;
    Ecoulomb  = 0.0;
    
    // Calculate the bond force term
    
    //#pragma omp for
    for(atomi=0 ; atomi < Natoms ; atomi++){
      Fterm_[atomi].x=0.0;
      Fterm_[atomi].y=0.0;
      Fterm_[atomi].z=0.0;
    }
    
#pragma omp for
    for(atomi=0 ; atomi < Natoms ; atomi++){
      particle[atomi].Fx=0.0;
      particle[atomi].Fy=0.0;
      particle[atomi].Fz=0.0;
    }
    //start = clock();
    
    if (potential_flag==0 || potential_flag==1)
      calc_dE_REBO(Fterm_,angle,BondNum,Natoms,L,BondParams,Bonded,particle,Fc,NBond,E_REBO,NBond_once,Bonded_once,Ncar_hyd);
    
    
    if ((potential_flag==0 || potential_flag==1))
      TORSION_(particle,Fterm_,BondNum,Natoms,L,Bonded,NBond,E_TORSION,NBond_once,Bonded_once,Ncar_hyd);
    
    //diff = clock() - start;
    //msec = diff * 1000 / CLOCKS_PER_SEC;
    
    //printf("Torsion Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
    
    
    
    //if (potential_flag==3 || potential_flag==4 || potential_flag==5){
    // calc_F_Tersoff(angle,BondNum,Natoms,rij,L,BondParams,Bonded,particle,Fc,Fterm,NBond,E_Tersoff);
    
    if (potential_flag==2 || potential_flag==3 || BN_flag==1){
      //gettimeofday(&t1, NULL);
      calc_dE_Tersoff(Fterm_,angle,BondNum,Natoms,L,BondParams,Bonded,particle,Fc,NBond,E_Tersoff,NBond_once,Bonded_once,N_init);
      //gettimeofday(&t2, NULL);
      //elapsedTime = (t2.tv_sec - t1.tv_sec)*1000;      // sec to ms
      //elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000;   // us to ms
      //fprintf(logfile,"Tersoff forces Time taken %d seconds %d milliseconds \n",elapsedTime/1000 ,elapsedTime%1000);
    }
    
    
    if ((potential_flag==0 || potential_flag==2) && Interlayer > 0){
      //gettimeofday(&t1, NULL);
      if(ILP_flag==0)
	calc_F_KC(Fterm_,Natoms,particle,rij_long,BondParams,L,Normal,Normal_atom,Interlayer,Neighb,NNeighb,dN,E_KC,Ncar_hyd);
      else if(ILP_flag==1)
	calc_F_KC(Fterm_,Natoms,particle,rij_long,BondParams,L,Normal,Normal_atom,Interlayer,Neighb,NNeighb,dN,E_KC,Natoms);
      
      //gettimeofday(&t2, NULL);
      //elapsedTime = (t2.tv_sec - t1.tv_sec)*1000;      // sec to ms
      //elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000;   // us to ms
      //fprintf(logfile,"KC forces Time taken %d seconds %d milliseconds \n",elapsedTime/1000 ,elapsedTime%1000);
    }
  
  if (Interlayer > 0 && ILP_flag==0 &&(potential_flag==1 || potential_flag==3 || BN_flag==1)){
    //gettimeofday(&t1, NULL);
    calc_F_vdW(Fterm_,Natoms,particle,rij_long,BondParams,L,Normal,Normal_atom,Interlayer,Neighb,NNeighb,EvdW,N_init_ILP);
    //gettimeofday(&t2, NULL);
    //elapsedTime = (t2.tv_sec - t1.tv_sec)*1000;      // sec to ms
    //elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000;   // us to ms
    //fprintf(logfile,"vdW forces Time taken %d seconds %d milliseconds \n",elapsedTime/1000 ,elapsedTime%1000);
    //gettimeofday(&t1, NULL);
    calc_F_Rep(Fterm_,Natoms,particle,rij_long,BondParams,L,Normal,Normal_atom,Interlayer,Neighb,NNeighb,dN,ERep,N_init_ILP);
    //cerr<<"yoyo5"<<endl;
    //gettimeofday(&t2, NULL);
    //elapsedTime = (t2.tv_sec - t1.tv_sec)*1000;      // sec to ms
    //elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000;   // us to ms
    //fprintf(logfile,"Rep forces Time taken %d seconds %d milliseconds \n",elapsedTime/1000 ,elapsedTime%1000);
    //gettimeofday(&t1, NULL);
    calc_F_coul(Fterm_,particle,Natoms,Bmat,qvec,Avec,dqvec,AtomParams,BondParams,rij,rij_long,logfile,L,Neighb,NNeighb,potential_flag,Ecoulomb,N_init_ILP);
    //gettimeofday(&t2, NULL);
    //elapsedTime = (t2.tv_sec - t1.tv_sec)*1000;      // sec to ms
    //elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000;   // us to ms
    //fprintf(logfile,"Coulomd forces Time taken %d seconds %d milliseconds \n",elapsedTime/1000 ,elapsedTime%1000);
  }
  
  //exit(0);
  //diff = clock() - start;
  //msec = diff * 1000 / CLOCKS_PER_SEC;
  //printf("Inter forces Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
#pragma omp critical
  {
    for(atomi=0 ; atomi < Natoms ; atomi++){
      particle[atomi].Fx += Fterm_[atomi].x;
      particle[atomi].Fy += Fterm_[atomi].y;
      particle[atomi].Fz += Fterm_[atomi].z;
      //if(particle[atomi].type==0)printf("atomi=%i Fx=%.16f Fy=%.16f Fz=%.16f\n",atomi,particle[atomi].Fx,particle[atomi].Fy,particle[atomi].Fz);
    }
    
    E_tot += E_REBO + E_TORSION + E_KC + E_Tersoff + EvdW + ERep + Ecoulomb;
  }
  //fprintf(logfile,"Force E_REBO=%.16f E_Torsion=%.16f E_Tersoff=%.16f E_KC=%.16f  EvdW=%.16f ERep=%.16f Ecoulomb=%.16f\n",E_REBO,E_TORSION,E_Tersoff,E_KC,EvdW,ERep,Ecoulomb);

    delete [] dN;
    delete [] Fterm_;
  }
}

void calc_F_coul(point *Fcoul, atom *particle, int Natoms, double *Bmat, double *qvec, double *Avec, double *dqvec, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, double *rij, double *rij_long, FILE* logfile, point L,int *Neighb,int *NNeighb,int potential_flag,double &Ecoulomb,int N_init_ILP)
{
  Ecoulomb = 0.0;
  double Ecoulomb_ = 0.0;
  //#pragma omp parallel
  // {
  //int atomk;%
  //#pragma omp for
  // for(atomk=0 ; atomk < Natoms ; atomk++)
  //Fcoul[atomk].x = Fcoul[atomk].y = Fcoul[atomk].z = 0.0;
  //**
  //if(potential_flag == 1){
  
  //}
  // }
  //**
 
    int atomi, atomj, atomk, index,j;
    double term1,term2;
    double Tap, dTap, r_ij;
    double rijShiled,shield,gamma,pow_gamma;
    double drij_dr_k_x,drij_dr_k_y,drij_dr_k_z;
    double dTerm;
#pragma omp for schedule (dynamic,10) //reduction(+:Ecoulomb_)  //private(dqvec,Avec,atomi, atomj, atomk, index,j,term1,term2,Tap, dTap, r_ij,rijShiled,shield,gamma,drij_dr_k_x,drij_dr_k_y,drij_dr_k_z)     
    for(atomi=N_init_ILP ; atomi < Natoms ; atomi++){
      
      for (j=0 ; j < NNeighb[atomi]; j++){
	
	atomj = Neighb[atomi*MaxNeighb + j] ;

	if(particle[atomi].type ==0 || particle[atomj].type ==0)continue;
	
	//if(atomk != atomj && atomk !=atomi && potential_flag == 0)continue;
	
	//index = atomi*Natoms + atomj;
	
	r_ij  = R_ij(particle,atomi,atomj, L);
	//r_ij = R_ij_(particle[atomi].r[0],particle[atomj].r[0],particle[atomi].r[1],particle[atomj].r[1],particle[atomi].r[2],particle[atomj].r[2],L);
	//r_ij = rij_long[atomi*MaxNeighb + j];
	Tap  = calc_Tap(r_ij);
	dTap = calc_dTap(r_ij);
	gamma = BondParams[particle[atomi].type*NAtomParams + particle[atomj].type].original_gamma;
	pow_gamma = BondParams[particle[atomi].type*NAtomParams + particle[atomj].type].original_gamma;
	shield    = pow(r_ij,3.0) + pow(1.0/gamma,3.0);
	//shield    = pow(r_ij,3.0) + pow_gamma;
	rijShiled = 1.0/pow(shield,0.333333333333);
	term1=- particle[atomi].charge * particle[atomj].charge * fourth(rijShiled) * sqr(r_ij);
	term2= particle[atomi].charge * particle[atomj].charge * rijShiled;
	//#pragma omp atomic
	Ecoulomb_ += Tap*(term2);
	
	drij_dr_k_x = calc_dr_ij_drk(atomi,atomj,atomi,0,r_ij,particle,Natoms,L);
	drij_dr_k_y = calc_dr_ij_drk(atomi,atomj,atomi,1,r_ij,particle,Natoms,L);
	drij_dr_k_z = calc_dr_ij_drk(atomi,atomj,atomi,2,r_ij,particle,Natoms,L);
	
	//dTerm = (Tap * term1 + dTap * term2)*Kappa_eV_kcal;
	dTerm = (Tap * term1 + dTap * term2)*Kappa;

	Fcoul[atomi].x -= drij_dr_k_x*dTerm;

	Fcoul[atomi].y -= drij_dr_k_y*dTerm;

	Fcoul[atomi].z -= drij_dr_k_z*dTerm;

	Fcoul[atomj].x += drij_dr_k_x*dTerm;

	Fcoul[atomj].y += drij_dr_k_y*dTerm;

	Fcoul[atomj].z += drij_dr_k_z*dTerm;
	/*
	  if(potential_flag == 0){
	  Fcoul[atomk].x -= drij_dr_k_x*(Tap * term1 + dTap * term2);
	  Fcoul[atomk].y -= drij_dr_k_y*(Tap * term1 + dTap * term2);
	  Fcoul[atomk].z -= drij_dr_k_z*(Tap * term1 + dTap * term2);
	  }
	  else{
	  
	  }
	*/
      }
      //Fcoul[atomi].x *= (Kappa * eV2kcalmole);
      //Fcoul[atomi].y *= (Kappa * eV2kcalmole);
      //Fcoul[atomi].z *= (Kappa * eV2kcalmole);
    }
    //}
  //Ecoulomb *= (Kappa * eV2kcalmole);
  Ecoulomb_ *= (Kappa);
  Ecoulomb = Ecoulomb_;
}

void calc_F_vdW(point *FvdW, int Natoms, atom *particle, double *rij_long, BondParamsStruct *BondParams, point L,Normal_struct *Normal,int *Normal_atom,int Interlayer,int *Neighb,int *NNeighb,double &EvdW,int N_init_ILP)
{
  EvdW=0.0;
  double EvdW_=0.0;
 
    int atomk, atomi, typek, typei,atomj,typej, index,j;
    double alpha, rvdW, term1, term2,_term2, dEvdW_drij;
    double dTap, r_ij, Tap,C6,term6,term7;
    double R,dr_ij_drk,C6_r_ij;
    double reff,r_ik,r_jk;
    int n,dir;
    
    //#pragma omp for 
    //for(atomk=0 ; atomk < Natoms ; atomk++)
    //FvdW[atomk].x = FvdW[atomk].y = FvdW[atomk].z = 0.0;
    
#pragma omp for schedule (dynamic,10) //reduction(+:EvdW_)
    for(atomi=N_init_ILP ; atomi < Natoms ; atomi++){
      typei = particle[atomi].type;
      for (j=0 ; j < NNeighb[atomi]; j++){
	atomk = atomj = Neighb[atomi*MaxNeighb + j] ;	  
	typek = typej = particle[atomj].type;
	
	index = typej*NAtomParams + typei;
	r_ij  = R_ij(particle,atomi,atomj, L);
	//r_ij=R_ij_(particle[atomi].r[0],particle[atomj].r[0],particle[atomi].r[1],particle[atomj].r[1],particle[atomi].r[2],particle[atomj].r[2],L);
	//r_ij = rij_long[atomi*MaxNeighb + j];
	Tap   = calc_Tap(r_ij);
	dTap  = calc_dTap(r_ij);  
	
	//interlayer vdW
	reff    = BondParams[index].r_eff;	    
	R       = BondParams[index].rvdW_long;
	C6      = BondParams[index].C6;	
	term2   = exp(-d_TS*(r_ij/(Sr_TS*reff) - 1.0));
	_term2  = 1.0/(1.0 + term2);
	
	C6_r_ij = C6/(pow(r_ij,6.0)+0.1);
	
	term6   = (-(pow(_term2,2.0)*term2*(-d_TS/(Sr_TS*reff))*C6_r_ij ) + _term2*(-1)*C6*pow(pow(r_ij,6)+0.1,-2)*6*pow(r_ij,5));
	term7   = (-(_term2)*C6_r_ij)*dTap;
	
	//#pragma omp critical
	EvdW_ += Tap*(-(1.0/(1.0 + term2 ))*C6_r_ij);
	//********dirivative in the x direction******
	
	dir = 0;
	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	
	dEvdW_drij = Tap*(- (term6*dr_ij_drk)) + dr_ij_drk*term7;
	//#pragma omp atomic	
	FvdW[atomk].x -= dEvdW_drij;
	//#pragma omp atomic	
	FvdW[atomi].x += dEvdW_drij;
	
	//********dirivative in the y direction******
	
	dir = 1;
	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	
	dEvdW_drij = Tap*(- (term6*dr_ij_drk)) + dr_ij_drk*term7;
	
	//#pragma omp atomic	
	FvdW[atomk].y -= dEvdW_drij;
	//#pragma omp atomic	
	FvdW[atomi].y += dEvdW_drij;
	
	//********dirivative in the z direction******
	
	dir = 2;
	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	
	dEvdW_drij = Tap*(- (term6*dr_ij_drk)) + dr_ij_drk*term7;
	
	//#pragma omp atomic	
	FvdW[atomk].z -= dEvdW_drij;
	//#pragma omp atomic	
	FvdW[atomi].z += dEvdW_drij;
      }
    }
  
  EvdW=EvdW_;
}

void calc_F_Rep(point *FvdW, int Natoms, atom *particle, double *rij_long, BondParamsStruct *BondParams, point L,Normal_struct *Normal,int *Normal_atom,int Interlayer,int *Neighb,int *NNeighb,dr_ij *dN1,double &ERep, int N_init_ILP)
{
  ERep=0.0;
  double ERep_=0.0;

    int atomk, atomi, typek, typei,atomj,typej, index,j;
    double alpha, rvdW, Epsilon, term1, dEvdW_drij;
    double dTap, r_ij, Tap,Epsilon_C_vdW_exp,term3,term4,term5,term7;
    double dP_ij_drk,dP_ji_drk,R,Pij,Pji,gamma,dr_ij_drk,exp_term1;
    double exp_Pij_gamma, exp_Pji_gamma,Pij_gamma,Pji_gamma,reff,r_ik,r_jk;
     double rx_ij,ry_ij,rz_ij,rx_ji,ry_ji,rz_ji,inside_acos_ij,inside_acos_ji,Pij_term1,Pji_term1,Pij_term,Pji_term,dr_ji_drk;
    double angle_ij, sqrt_inside_ij, sin_angle_ij,cos_angle_ij,angle_ji,sqrt_inside_ji,sin_angle_ji,cos_angle_ji,dinside_acos,dangle;
    double dP_ij_term,dP_ji_term,d_Norm_length,Pow_normal_i,Pow_normal_j;
    double C_vdW,Nori_L,Norj_L;
    double temp_b,temp_c,dpij_flag,dpji_flag;
    int n,dir;
    dP_ij_drk = 0.0;
    dP_ji_drk = 0.0;
    dr_ij *dN = new dr_ij[1];
        
#pragma omp for schedule(dynamic,10) //reduction(+:ERep_)
    for(atomi=N_init_ILP ; atomi < Natoms ; atomi++){
      typei = particle[atomi].type;
      for (j=0 ; j < NNeighb[atomi]; j++){
	
	atomk = atomj = Neighb[atomi*MaxNeighb + j] ;	  
	typek = typej = particle[atomj].type;
	//cerr<<"typej= "<<typej<<" typei= "<<typei<<endl; 
	index = typej*NAtomParams + typei;
	//r_ij  = R_ij(particle,atomi,atomj, L);

	rx_ij = R_PBC(particle,atomi,atomj,0,L.x);
	ry_ij = R_PBC(particle,atomi,atomj,1,L.y);
	rz_ij = R_PBC(particle,atomi,atomj,2,L.z);
	
	rx_ji = -rx_ij;
	ry_ji = -ry_ij;
	rz_ji = -rz_ij;
	
	r_ij =sqrt(rx_ij*rx_ij + ry_ij*ry_ij + rz_ij*rz_ij);
	//r_ij=R_ij_(particle[atomi].r[0],particle[atomj].r[0],particle[atomi].r[1],particle[atomj].r[1],particle[atomi].r[2],particle[atomj].r[2],L);
	if(r_ij > non_bond_cut_Rep)continue;
		
	Tap   = calc_Tap_Rep(r_ij);
	dTap  = calc_dTap_Rep(r_ij);  
	//r_ij2 = r_ij*r_ij;
	reff    = BondParams[index].r_eff;	    
	alpha   = BondParams[index].alpha_long;
	Epsilon = BondParams[index].Epsilon_long;
	R       = BondParams[index].rvdW_long;
	gamma   = BondParams[index].Trans;
    	C_vdW   = BondParams[index].C_vdW;
	//****
	dP_ij_drk=0.0;
	dP_ji_drk=0.0;
	Pij_term1=0.0;
	Pji_term1=0.0;
	dN[0].r[0]=0;
	dN[0].r[1]=0;
	dN[0].r[2]=0;
	dpij_flag = 1.0;
	dpji_flag = 1.0;

	
	
	Nori_L=rx_ij*Normal[atomi].x + ry_ij*Normal[atomi].y + rz_ij*Normal[atomi].z;
	Norj_L=rx_ji*Normal[atomj].x + ry_ji*Normal[atomj].y + rz_ji*Normal[atomj].z;

	inside_acos_ij = (Nori_L) / (r_ij);
	inside_acos_ji = (Norj_L) / (r_ij);

	if(((r_ij*r_ij-sqr(Nori_L)) <= 0) || Interlayer==2){
	  Pij=0.0;
	  dpij_flag = 0.0;
	}
	else Pij=sqrt(r_ij*r_ij - sqr(Nori_L));

	if(((r_ij*r_ij-sqr(Norj_L)) <= 0) || Interlayer==2){
	  Pji=0.0;
	  dpji_flag = 0.0;
	}
	else Pji=sqrt(r_ij*r_ij - sqr(Norj_L));
	
	Pij_gamma = Pij/gamma; //new
	Pji_gamma = Pji/gamma; //new
	exp_Pij_gamma = exp(-pow((Pij_gamma),2.0)); //new
	exp_Pji_gamma = exp(-pow((Pji_gamma),2.0)); //new
	exp_term1 = exp(alpha*(1.0-r_ij/R));

	Epsilon_C_vdW_exp = Epsilon + C_vdW*(exp_Pij_gamma + exp_Pji_gamma);
	//printf("exp_term1=%.16f \n",exp_term1);
	//printf("Epsilon_C_vdW_exp=%.16f \n",Epsilon_C_vdW_exp);
	//#pragma omp critical
	
	//ERep_ += (exp_term1*(Epsilon_C_vdW_exp) );
	ERep_ += Tap*(exp_term1*(Epsilon_C_vdW_exp) );
	
	term3   = exp_term1*(-alpha/R)*(Epsilon_C_vdW_exp);
	term4   = exp_term1*C_vdW*(exp_Pij_gamma*(-2.0)/sqr(gamma))*dpij_flag;
	term5   = exp_term1*C_vdW*(exp_Pji_gamma*(-2.0)/sqr(gamma))*dpji_flag;
	term7   = ((exp_term1*(Epsilon_C_vdW_exp)))*dTap;

	//********dirivative in the x direction******
	
	dir = 0;
	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	
	dP_ij_drk=r_ij*dr_ij_drk+Nori_L*Normal[atomi].x;
	dP_ji_drk=r_ij*dr_ij_drk-Norj_L*Normal[atomj].x;
	
	dEvdW_drij = Tap*(term3*dr_ij_drk + term4*dP_ij_drk + term5*dP_ji_drk) + dr_ij_drk*term7;
	

	FvdW[atomk].x -= dEvdW_drij;

	FvdW[atomi].x += dEvdW_drij;
	
	//********dirivative in the y direction******
	
	dir = 1;
	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);

	dP_ij_drk=r_ij*dr_ij_drk+Nori_L*Normal[atomi].y;
	dP_ji_drk=r_ij*dr_ij_drk-Norj_L*Normal[atomj].y;
	
	dEvdW_drij = Tap*(term3*dr_ij_drk + term4*dP_ij_drk + term5*dP_ji_drk ) + dr_ij_drk*term7;


	FvdW[atomk].y -= dEvdW_drij;

	FvdW[atomi].y += dEvdW_drij;
	
	//********dirivative in the z direction******
	
	dir = 2;

	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);

	dP_ij_drk=r_ij*dr_ij_drk+Nori_L*Normal[atomi].z;
	dP_ji_drk=r_ij*dr_ij_drk-Norj_L*Normal[atomj].z;
	
	dEvdW_drij = Tap*(term3*dr_ij_drk + term4*dP_ij_drk + term5*dP_ji_drk ) + dr_ij_drk*term7;

	FvdW[atomk].z -= dEvdW_drij;

	FvdW[atomi].z += dEvdW_drij;

	//****Normal[atomi].vector_2y
	
	if(((r_ij*r_ij-sqr(Nori_L)) <= 0) || Interlayer==2);
	else{
	  
	  temp_b = Tap*term4*(-Nori_L)*( (ry_ij*Normal[atomi].vector_2z - rz_ij*Normal[atomi].vector_2y) - Nori_L*( Normal[atomi].y*Normal[atomi].vector_2z - Normal[atomi].z*Normal[atomi].vector_2y))/(Normal[atomi].Normal_length);
	  temp_c = Tap*term4*(-Nori_L)*((-ry_ij*Normal[atomi].vector_1z + rz_ij*Normal[atomi].vector_1y) - Nori_L*(-Normal[atomi].y*Normal[atomi].vector_1z + Normal[atomi].z*Normal[atomi].vector_1y))/(Normal[atomi].Normal_length);
	  

	  FvdW[Normal_atom[atomi*3 + 0]].x+=( temp_b + temp_c );

	  FvdW[Normal_atom[atomi*3 + 1]].x-=temp_b;

	  FvdW[Normal_atom[atomi*3 + 2]].x-=temp_c;
	  
	  temp_b = Tap*term4*(-Nori_L)*((-rx_ij*Normal[atomi].vector_2z + rz_ij*Normal[atomi].vector_2x) - Nori_L*(-Normal[atomi].x*Normal[atomi].vector_2z + Normal[atomi].z*Normal[atomi].vector_2x))/(Normal[atomi].Normal_length);
	  temp_c = Tap*term4*(-Nori_L)*(( rx_ij*Normal[atomi].vector_1z - rz_ij*Normal[atomi].vector_1x) - Nori_L*( Normal[atomi].x*Normal[atomi].vector_1z - Normal[atomi].z*Normal[atomi].vector_1x))/(Normal[atomi].Normal_length);
	  

	  FvdW[Normal_atom[atomi*3 + 0]].y+=( temp_b + temp_c );

	  FvdW[Normal_atom[atomi*3 + 1]].y-=temp_b;

	  FvdW[Normal_atom[atomi*3 + 2]].y-=temp_c;
	  
	  temp_b = Tap*term4*(-Nori_L)*(( rx_ij*Normal[atomi].vector_2y - ry_ij*Normal[atomi].vector_2x) - Nori_L*( Normal[atomi].x*Normal[atomi].vector_2y - Normal[atomi].y*Normal[atomi].vector_2x))/(Normal[atomi].Normal_length);
	  temp_c = Tap*term4*(-Nori_L)*((-rx_ij*Normal[atomi].vector_1y + ry_ij*Normal[atomi].vector_1x) - Nori_L*(-Normal[atomi].x*Normal[atomi].vector_1y + Normal[atomi].y*Normal[atomi].vector_1x))/(Normal[atomi].Normal_length);
	  

	  FvdW[Normal_atom[atomi*3 + 0]].z+=( temp_b + temp_c );

	  FvdW[Normal_atom[atomi*3 + 1]].z-=temp_b;

	  FvdW[Normal_atom[atomi*3 + 2]].z-=temp_c;
	}
	//****************************************
	//****************************************
	
	if(((r_ij*r_ij-sqr(Norj_L)) <= 0) || Interlayer==2);
	else{
	  
	  temp_b = Tap*term5*(-Norj_L)*(( ry_ji*Normal[atomj].vector_2z - rz_ji*Normal[atomj].vector_2y) - Norj_L*( Normal[atomj].y*Normal[atomj].vector_2z - Normal[atomj].z*Normal[atomj].vector_2y))/(Normal[atomj].Normal_length);
	  temp_c = Tap*term5*(-Norj_L)*((-ry_ji*Normal[atomj].vector_1z + rz_ji*Normal[atomj].vector_1y) - Norj_L*(-Normal[atomj].y*Normal[atomj].vector_1z + Normal[atomj].z*Normal[atomj].vector_1y))/(Normal[atomj].Normal_length);

	  FvdW[Normal_atom[atomj*3 + 0]].x+=( temp_b + temp_c );

	  FvdW[Normal_atom[atomj*3 + 1]].x-=temp_b;

	  FvdW[Normal_atom[atomj*3 + 2]].x-=temp_c;
	  
	  temp_b = Tap*term5*(-Norj_L)*((-rx_ji*Normal[atomj].vector_2z + rz_ji*Normal[atomj].vector_2x) - Norj_L*(-Normal[atomj].x*Normal[atomj].vector_2z + Normal[atomj].z*Normal[atomj].vector_2x))/(Normal[atomj].Normal_length);
	  temp_c = Tap*term5*(-Norj_L)*(( rx_ji*Normal[atomj].vector_1z - rz_ji*Normal[atomj].vector_1x) - Norj_L*( Normal[atomj].x*Normal[atomj].vector_1z - Normal[atomj].z*Normal[atomj].vector_1x))/(Normal[atomj].Normal_length);

	  FvdW[Normal_atom[atomj*3 + 0]].y+=( temp_b + temp_c );

	  FvdW[Normal_atom[atomj*3 + 1]].y-=temp_b;

	  FvdW[Normal_atom[atomj*3 + 2]].y-=temp_c;
	  
	  temp_b = Tap*term5*(-Norj_L)*(( rx_ji*Normal[atomj].vector_2y - ry_ji*Normal[atomj].vector_2x) - Norj_L*( Normal[atomj].x*Normal[atomj].vector_2y - Normal[atomj].y*Normal[atomj].vector_2x))/(Normal[atomj].Normal_length);
	  temp_c = Tap*term5*(-Norj_L)*((-rx_ji*Normal[atomj].vector_1y + ry_ji*Normal[atomj].vector_1x) - Norj_L*(-Normal[atomj].x*Normal[atomj].vector_1y + Normal[atomj].y*Normal[atomj].vector_1x))/(Normal[atomj].Normal_length);

	  FvdW[Normal_atom[atomj*3 + 0]].z+=( temp_b + temp_c );

	  FvdW[Normal_atom[atomj*3 + 1]].z-=temp_b;

	  FvdW[Normal_atom[atomj*3 + 2]].z-=temp_c;
	}
      }
    }
    delete [] dN;
    //}
  ERep=ERep_;
}

/* ---------------------------------------------------------------------- */

void calc_dE_Tersoff(point *F_REBO,double *angle,int *BondNum,int Natoms,point L,BondParamsStruct *BondParams,int *Bonded,atom *particle,double *Fc,int *NBond,double &E_Tersoff_,int *NBond_once,int *Bonded_once,int N_init)
{
  
  double evdwl = 0.0;


    int i,j,k,ii,jj,kk,inum,jnum;
    int itype,jtype,ktype,iparam_ij,iparam_ijk;
    //tagint itag,jtag;
    double xtmp,ytmp,ztmp,delx,dely,delz,fpair;
    double rsq,rsq1,rsq2;
    double delr1[3],delr2[3],fi[3],fj[3],fk[3];
    double zeta_ij,prefactor;
    double c1,c2,c3,c4,bigr,biga,cut,cutsq,gamma,d,c,h,beta,bigb,bigd,powern,lam1,lam2,lam3;
    int powermint,atomi,itag,atomj,jtag,atomk,ktag;
    double shift_,rsq_tmp; 
    
#pragma omp for //reduction(+:evdwl)
    for (atomi = N_init; atomi < Natoms; atomi++) {
      
      itag=atomi;
      itype = particle[atomi].type;
      xtmp = particle[atomi].r[0];
      ytmp = particle[atomi].r[1];
      ztmp = particle[atomi].r[2];
      
      for (j = 0; j < NBond_once[atomi]; j++) {
	atomj = Bonded_once[atomi*MaxNBond + j];
	jtag = atomj;
	jtype = particle[atomj].type;
      
      powern  = BondParams[itype*NAtomParams + jtype].Tersoff_n;
      bigr    = BondParams[itype*NAtomParams + jtype].Tersoff_R;
      bigd    = BondParams[itype*NAtomParams + jtype].Tersoff_D;
      biga    = BondParams[itype*NAtomParams + jtype].Tersoff_A;
      cut   = bigr + bigd;
      cutsq = cut*cut;
      c1 = pow(2.0*powern*1.0e-16,-1.0/powern);
      c2 = pow(2.0*powern*1.0e-8,-1.0/powern);
      c3 = 1.0/c2;
      c4 = 1.0/c1;
      lam1  = BondParams[itype*NAtomParams + jtype].Tersoff_lambda1;
      delx = R_PBC(particle,atomi,atomj,0,L.x);
      dely = R_PBC(particle,atomi,atomj,1,L.y);	
      delz = R_PBC(particle,atomi,atomj,2,L.z);
      
      rsq = delx*delx + dely*dely + delz*delz;

      rsq_tmp = sqrt(rsq);
      
      if(itype == 0 && jtype == 0){
	rsq += shift*shift + 2*sqrt(rsq)*shift;
	shift_=shift;
      }
      else shift_=0.0;
      
      if (rsq > cutsq) continue;

	repulsive(rsq,fpair,evdwl,lam1,biga,bigr,bigd);
      
	F_REBO[atomi].x += delx*fpair/rsq_tmp;

	F_REBO[atomi].y += dely*fpair/rsq_tmp;

	F_REBO[atomi].z += delz*fpair/rsq_tmp;

	F_REBO[atomj].x -= delx*fpair/rsq_tmp;

	F_REBO[atomj].y -= dely*fpair/rsq_tmp;

	F_REBO[atomj].z -= delz*fpair/rsq_tmp;
	//}
    }
    for (j = 0; j < NBond[atomi]; j++) {
      
      atomj = Bonded[atomi*MaxNBond + j];
      jtag = atomj;
      jtype = particle[atomj].type;
      delr1[0] = R_PBC(particle,atomj,atomi,0,L.x);
      delr1[1] = R_PBC(particle,atomj,atomi,1,L.y);	
      delr1[2] = R_PBC(particle,atomj,atomi,2,L.z);
      
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      rsq_tmp = sqrt(rsq1);
      
      if(itype == 0 && jtype == 0){
	rsq1 += shift*shift + 2.0*rsq_tmp*shift;
	shift_=shift;
      }
      else shift_=0.0;
      
      bigr    = BondParams[itype*NAtomParams + jtype].Tersoff_R;
      bigd    = BondParams[itype*NAtomParams + jtype].Tersoff_D;
      cut   = bigr + bigd;
      cutsq = cut*cut;
      
      if (rsq1 > cutsq) continue;
      zeta_ij = 0.0;
      for (k = 0; k < NBond[atomi]; k++) {
	if (j == k) continue;
	atomk = Bonded[atomi*MaxNBond + k];
	ktag = atomk;
	ktype = particle[atomk].type;
	delr2[0] = R_PBC(particle,atomk,atomi,0,L.x);
	delr2[1] = R_PBC(particle,atomk,atomi,1,L.y);	
	delr2[2] = R_PBC(particle,atomk,atomi,2,L.z);
	
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
	
	if(itype == 0 && ktype == 0){
	  rsq2 += shift*shift + 2.0*sqrt(rsq2)*shift;
	  shift_=shift;
	}
	else shift_=0.0;
	
	bigr    = BondParams[itype*NAtomParams + ktype].Tersoff_R;
	bigd    = BondParams[itype*NAtomParams + ktype].Tersoff_D;
	cut   = bigr + bigd;
	cutsq = cut*cut;
	
        if (rsq2 > cutsq) continue;

	gamma     = 1.0;
	powermint = 3;
	lam1  = BondParams[itype*NAtomParams + ktype].Tersoff_lambda1;
	lam3  = BondParams[itype*NAtomParams + ktype].Tersoff_lambda3;
	c  = BondParams[itype*NAtomParams + ktype].Tersoff_c;
	h  = BondParams[itype*NAtomParams + ktype].Tersoff_h;
	d  = BondParams[itype*NAtomParams + ktype].Tersoff_d;
	
	zeta_ij += zeta(rsq1,rsq2,delr1,delr2,powermint,lam1,lam3,c,d,h,gamma, bigr,bigd,shift_);
      }
      
      // pairwise force due to zeta
      lam2  = BondParams[itype*NAtomParams + jtype].Tersoff_lambda2;
      beta    = BondParams[itype*NAtomParams + jtype].Tersoff_beta;
      powern  = BondParams[itype*NAtomParams + jtype].Tersoff_n;
      bigb    = BondParams[itype*NAtomParams + jtype].Tersoff_B;
      
      c1 = pow(2.0*powern*1.0e-16,-1.0/powern);
      c2 = pow(2.0*powern*1.0e-8,-1.0/powern);
      c3 = 1.0/c2;
      c4 = 1.0/c1;
      //#pragma omp critical
      //{
      //#pragma omp atomic
      force_zeta(rsq1,zeta_ij,fpair,prefactor,evdwl,c1,c2,c3,c4,powern,beta,bigr,bigb,bigd,lam2);

      F_REBO[atomi].x += delr1[0]*fpair/rsq_tmp;

      F_REBO[atomi].y += delr1[1]*fpair/rsq_tmp;

      F_REBO[atomi].z += delr1[2]*fpair/rsq_tmp;

      F_REBO[atomj].x -= delr1[0]*fpair/rsq_tmp;

      F_REBO[atomj].y -= delr1[1]*fpair/rsq_tmp;

      F_REBO[atomj].z -= delr1[2]*fpair/rsq_tmp;
    //}
      for (k = 0; k < NBond[atomi]; k++) {
	if (j == k) continue;
	
	atomk = Bonded[atomi*MaxNBond + k];
	ktag = atomk;
	ktype = particle[atomk].type;
	delr2[0] = R_PBC(particle,atomk,atomi,0,L.x);
	delr2[1] = R_PBC(particle,atomk,atomi,1,L.y);	
	delr2[2] = R_PBC(particle,atomk,atomi,2,L.z);
	
	rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
	
	if(itype == 0 && jtype == 0){
	  rsq2 += shift*shift + 2*sqrt(rsq2)*shift;
	  shift_=shift;
	}
	else shift_=0.0;
	
	bigr    = BondParams[itype*NAtomParams + ktype].Tersoff_R;
	bigd    = BondParams[itype*NAtomParams + ktype].Tersoff_D;
	cut   = bigr + bigd;
	cutsq = cut*cut;
	
        if (rsq2 > cutsq) continue;
	
	gamma     = 1.0;
	powermint = 3;
	c  = BondParams[itype*NAtomParams + ktype].Tersoff_c;
	h  = BondParams[itype*NAtomParams + ktype].Tersoff_h;
	d  = BondParams[itype*NAtomParams + ktype].Tersoff_d;
	lam3  = BondParams[itype*NAtomParams + ktype].Tersoff_lambda3;
	//#pragma omp critical
	//{
	//#pragma omp atomic
	attractive(prefactor,rsq1,rsq2,delr1,delr2,fi,fj,fk,powermint,lam3,c,d,h,gamma,bigr,bigd,shift_);


	F_REBO[atomi].x += fi[0];

	F_REBO[atomi].y += fi[1];

	F_REBO[atomi].z += fi[2];

	F_REBO[atomj].x += fj[0];

	F_REBO[atomj].y += fj[1];

	F_REBO[atomj].z += fj[2];

	F_REBO[atomk].x += fk[0];

	F_REBO[atomk].y += fk[1];

	F_REBO[atomk].z += fk[2];
	//}
      }
    }
  }
  
  E_Tersoff_= evdwl;
  
}

void repulsive(double rsq, double &fforce, double &eng,double lam1,double biga,double bigr,double bigd)
{
  double r,tmp_fc,tmp_fc_d,tmp_exp;
  
  r        = sqrt(rsq);
  tmp_fc   = ters_fc(r,bigr,bigd);//ters_fc(r,param);
  tmp_fc_d = ters_fc_d(r,bigr,bigd); //ters_fc_d(r,param);
  tmp_exp  = exp(-lam1 * r);
  //fforce   = -biga * tmp_exp * (tmp_fc_d - tmp_fc*lam1) / r;
  fforce   = -biga * tmp_exp * (tmp_fc_d - tmp_fc*lam1) ;
  eng      += tmp_fc * biga * tmp_exp;
  
  //printf("r= %.16f tmp_fc= %.16f tmp_fc_d= %.16f tmp_exp= %.16f fforce= %.16f eng= %.16f\n",r,tmp_fc,tmp_fc_d,tmp_exp,fforce,eng);
}

/* ---------------------------------------------------------------------- */

double zeta(double rsqij, double rsqik,double *delrij, double *delrik,int powermint,double lam1,double lam3,double c,double d,double h,double gamma,double bigr,double bigd,double shift_)
{
  double rij,rik,costheta,arg,ex_delr;
  
  rij = sqrt(rsqij);
  rik = sqrt(rsqik);
  costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] + delrij[2]*delrik[2]) / ((rij-shift_)*(rik-shift_));
  //costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] + delrij[2]*delrik[2]) / ((rij)*(rik));

if (powermint == 3) arg = pow(lam3 * (rij-rik),3.0);
  else arg = lam3 * (rij-rik);
  
  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);
  
  return ters_fc(rik,bigr,bigd) * ters_gijk(costheta,c,d,h,gamma) * ex_delr;
}

/* ---------------------------------------------------------------------- */

void force_zeta(double rsq, double zeta_ij,double &fforce, double &prefactor, double &eng,double c1,double c2,double c3,double c4,double powern,double beta,double bigr,double bigb,double bigd,double lam2)
{
  double r,fa,fa_d,bij;
  
  r = sqrt(rsq);
  //fa = ters_fa(r,param);
  fa= ters_fa(r,bigr,bigb,bigd,lam2);
  //fa_d = ters_fa_d(r,param);
  fa_d = ters_fa_d(r,bigr,bigd,lam2,bigb);
  //bij = ters_bij(zeta_ij,param);
  bij =  ters_bij(zeta_ij,c1,c2,c3,c4,powern, beta);
  //fforce = 0.5*bij*fa_d / r;
  fforce = 0.5*bij*fa_d;
  //prefactor = -0.5*fa * ters_bij_d(zeta_ij,param);
  prefactor = -0.5*fa *  ters_bij_d(zeta_ij,c1,c2,c3,c4,beta,powern);
  //if (eflag) eng = 0.5*bij*fa;
  eng += 0.5*bij*fa;
  //printf("r= %.16f fa= %.16f fa_d= %.16f bij= %.16f fforce= %.16f prefactor= %.16f eng= %.16f\n",r,fa,fa_d,bij,fforce,prefactor,eng);
  
}

/* ----------------------------------------------------------------------
   attractive term
   use param_ij cutoff for rij test
   use param_ijk cutoff for rik test
   ------------------------------------------------------------------------- */

void attractive(double prefactor,double rsqij, double rsqik,double *delrij, double *delrik,double *fi, double *fj, double *fk,int powermint,double lam3,double c,double d,double h,double gamma,double bigr,double bigd,double shift_)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;
  double rij_hat_[3],rik_hat_[3];
  
  rij = sqrt(rsqij);
  rijinv = 1.0/(rij-shift_);
  vec3_scale(rijinv,delrij,rij_hat);
  vec3_scale(1.0/rij,delrij,rij_hat_);
 
  rik = sqrt(rsqik);
  rikinv = 1.0/(rik-shift_);
  vec3_scale(rikinv,delrik,rik_hat);
  vec3_scale(1.0/rik,delrik,rik_hat_);
  //printf("rikinv= %.16f rik= %.16f delrik= %.16f  rik_hat0= %.16f rik_hat1= %.16f rik_hat2= %.16f\n",rikinv,rik,delrik,rik_hat[0],rik_hat[1],rij_hat[2]);
  //printf("cor_theta= %.16f rij_hat1= %.16f rik_hat1= %.16f \n",cos_theta,rij_hat[1],rik_hat[1]);
  //printf("cor_theta= %.16f rij_hat2= %.16f rik_h/rsq_tmpat2= %.16f \n",cos_theta,rij_hat[2],rik_hat[2]);
  ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,powermint,lam3,c,d,h,gamma,bigr,bigd,shift_,rij_hat_,rik_hat_);
  //ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,param);
}

/* ---------------------------------------------------------------------- */

double ters_fc(double r,double bigr,double bigd)
{
  double ters_R = bigr;
  double ters_D = bigd;

  if (r < ters_R-ters_D) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  return 0.5*(1.0 - sin(MY_PI2*(r - ters_R)/ters_D));
}

/* ---------------------------------------------------------------------- */

double ters_fc_d(double r,double bigr,double bigd)
{
  double ters_R = bigr;
  double ters_D = bigd;

  if (r < ters_R-ters_D) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  return -(MY_PI4/ters_D) * cos(MY_PI2*(r - ters_R)/ters_D);
}

/* ---------------------------------------------------------------------- */

double ters_fa(double r,double bigr,double bigb,double bigd,double lam2)
{
  if (r > bigr + bigd) return 0.0;
  return -bigb * exp(-lam2 * r) * ters_fc(r,bigr,bigd);
}

/* ---------------------------------------------------------------------- */

double ters_fa_d(double r,double bigr,double bigd,double lam2,double bigb)
{
  if (r > bigr + bigd) return 0.0;
  return bigb * exp(-lam2 * r) *
    (lam2 * ters_fc(r,bigr,bigd) - ters_fc_d(r,bigr,bigd));
}

/* ---------------------------------------------------------------------- */

double ters_bij(double zeta, double c1,double c2,double c3,double c4,double powern,double beta)
{
  double tmp = beta * zeta;
  if (tmp > c1) return 1.0/sqrt(tmp);
  if (tmp > c2)
    return (1.0 - pow(tmp,-powern) / (2.0*powern))/sqrt(tmp);
  if (tmp < c4) return 1.0;
  if (tmp < c3)
    return 1.0 - pow(tmp,powern)/(2.0*powern);
  return pow(1.0 + pow(tmp,powern), -1.0/(2.0*powern));
}

/* ---------------------------------------------------------------------- */

double ters_bij_d(double zeta,double c1,double c2,double c3,double c4,double beta,double powern)
{
  double tmp = beta * zeta;
  if (tmp > c1) return beta * -0.5*pow(tmp,-1.5);
  if (tmp > c2)
    return beta * (-0.5*pow(tmp,-1.5) *
                          (1.0 - 0.5*(1.0 +  1.0/(2.0*powern)) *
                           pow(tmp,-powern)));
  if (tmp < c4) return 0.0;
  if (tmp < c3)
    return -0.5*beta * pow(tmp,powern-1.0);

  double tmp_n = pow(tmp,powern);
  return -0.5 * pow(1.0+tmp_n, -1.0-(1.0/(2.0*powern)))*tmp_n / zeta;
}

/* ---------------------------------------------------------------------- */

void ters_zetaterm_d(double prefactor,double *rij_hat, double rij,double *rik_hat, double rik,double *dri, double *drj, double *drk,int powermint,double lam3,double c,double d,double h,double gamma,double bigr,double bigd,double shift_,double *rij_hat_,double *rik_hat_)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];
  
  //fc = ters_fc(rik,param);
  fc = ters_fc(rik,bigr,bigd) ;
  //dfc = ters_fc_d(rik,param);
  dfc = ters_fc_d(rik,bigr,bigd);
  if (powermint == 3) tmp = pow(lam3 * (rij-rik),3.0);
  else tmp = lam3 * (rij-rik);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (powermint == 3)
    ex_delr_d = 3.0*pow(lam3,3.0) * pow(rij-rik,2.0)*ex_delr;
  else ex_delr_d = lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  
  //gijk = ters_gijk(cos_theta,param);
  gijk   = ters_gijk(cos_theta,c,d,h,gamma);
  //printf("gijk= %.16f cor_theta= %.16f c= %.16f d=%.16f h= %.16f",gijk,cos_theta,c,d,h);
  gijk_d = ters_gijk_d(cos_theta,c,d,h,gamma);
  //gijk_d = ters_gijk_d(cos_theta,param);
  //costheta_d(rij_hat,rij-shift_,rik_hat,rik-shift_,dcosdri,dcosdrj,dcosdrk,shift_,rij_hat_,rik_hat_);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk,shift_,rij_hat_,rik_hat_);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);
  //printf("dri[0]= %.16f dfc= %.16f gijk= %.16f ex_delr= %.16f rik_hat= %.16f",dri[0],dfc,gijk,ex_delr,rik_hat);
  vec3_scale(-dfc*gijk*ex_delr,rik_hat,dri);
  
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdri,dri,dri);
  vec3_scaleadd(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  vec3_scale(prefactor,dri,dri);
  //printf("dri[0]= %.16f ",dri[0]);
  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  vec3_scale(fc*gijk_d*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);
  
  
}
/*
void read_file_Tersoff(char *file)
{
  int params_per_line = 17;
  char **words = new char*[params_per_line+1];

  //memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  //if (comm->me == 0) {
    //fp = force->open_potential(file);
    //if (fp == NULL) {
      //char str[128];
      //sprintf(str,"Cannot open Tersoff potential file %s",file);
      //error->one(FLERR,str);
    //}
  //}

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  int n,nwords,ielement,jelement,kelement;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    //if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
      //}
    //MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    //MPI_Bcast(&n,1,MPI_INT,0,world);
    //MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      //MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      //MPI_Bcast(&n,1,MPI_INT,0,world);
      //MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in Tersoff potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    // ielement,jelement,kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next line

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    for (kelement = 0; kelement < nelements; kelement++)
      if (strcmp(words[2],elements[kelement]) == 0) break;
    if (kelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                          "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].kelement = kelement;
    params[nparams].powerm = atof(words[3]);
    params[nparams].gamma = atof(words[4]);
    params[nparams].lam3 = atof(words[5]);
    params[nparams].c = atof(words[6]);
    params[nparams].d = atof(words[7]);
    params[nparams].h = atof(words[8]);
    params[nparams].powern = atof(words[9]);
    params[nparams].beta = atof(words[10]);
    params[nparams].lam2 = atof(words[11]);
    params[nparams].bigb = atof(words[12]);
    params[nparams].bigr = atof(words[13]);
    params[nparams].bigd = atof(words[14]);
    params[nparams].lam1 = atof(words[15]);
    params[nparams].biga = atof(words[16]);

    // currently only allow m exponent of 1 or 3

    params[nparams].powermint = int(params[nparams].powerm);

    if (params[nparams].c < 0.0 || params[nparams].d < 0.0 ||
        params[nparams].powern < 0.0 || params[nparams].beta < 0.0 ||
        params[nparams].lam2 < 0.0 || params[nparams].bigb < 0.0 ||
        params[nparams].bigr < 0.0 ||params[nparams].bigd < 0.0 ||
        params[nparams].bigd > params[nparams].bigr ||
        params[nparams].lam1 < 0.0 || params[nparams].biga < 0.0 ||
        params[nparams].powerm - params[nparams].powermint != 0.0 ||
        (params[nparams].powermint != 3 && params[nparams].powermint != 1) ||
        params[nparams].gamma < 0.0)
      error->all(FLERR,"Illegal Tersoff parameter");

    nparams++;
  }

  delete [] words;
}
*/
 /* ---------------------------------------------------------------------- */

void costheta_d(double *rij_hat, double rij,double *rik_hat, double rik,double *dri, double *drj, double *drk,double shift_,double *rij_hat_,double *rik_hat_)
{
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk
  
  double cos_theta = vec3_dot(rij_hat,rik_hat);
  
  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(1.0/(rij-shift_),drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(1.0/(rik-shift_),drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale(-1.0,dri,dri);
}





double calc_dE_Tersoff_ij_drn(int atomi,int atomj,int atomn,int dir,int i,double *angle,int *BondNum,int Natoms,double *rij,point L,BondParamsStruct *BondParams,int *Bonded,atom *particle,double *Fc,int *NBond,double &Etmp){
  Etmp=0.0;
  int index,atomk,typek,typej,typei,k,index2;
  double g,Theta_ijk,b_ij,r_ij,r_jk,Tersoff_lambda1,Tersoff_lambda2,Tersoff_lambda3,Tersoff_c,Tersoff_d,Tersoff_h,Tersoff_n,A,B,beta;
  double Xsi,ETersoff,dTheta_ijk_drn,dg_drn,dXsi_drn,Fc_ij,Fc_jk,Tersoff_R,Tersoff_D;
  double dr_ij_drn,dr_jk_drn,dFc_ij_drn,dFc_jk_drn,Fa,Fr,dFa_drn,dFr_drn,db_ij_drn,n,dE_Tersoff_ij_drn;
  double Tersoff_sqr_c,Tersoff_sqr_d,sqr_Tersoff_h_Cos_Theta;
  //i=BondNum[atomj*Natoms + atomi];
  double inside_acos,cos_ijk;
  dr_ij rji,rjk;
 
  //r_ij       = rij[atomj*MaxNBond + i];
  //r_ij = R_ij_(particle[atomi].r[0],particle[atomj].r[0],particle[atomi].r[1],particle[atomj].r[1],particle[atomi].r[2],particle[atomj].r[2],L);
  r_ij       = R_ij(particle,atomi,atomj, L);
  typei      = particle[atomi].type;
  typej      = particle[atomj].type;
  index      = typei*NAtomParams + typej;
  
  Tersoff_n  = BondParams[index].Tersoff_n;
  A          = BondParams[index].Tersoff_A;
  B          = BondParams[index].Tersoff_B;
  beta       = BondParams[index].Tersoff_beta;
  Tersoff_R  = BondParams[index].Tersoff_R;
  Tersoff_D  = BondParams[index].Tersoff_D;
  Tersoff_lambda1 = BondParams[index].Tersoff_lambda1;
  Tersoff_lambda2 = BondParams[index].Tersoff_lambda2;
  Tersoff_lambda3 = BondParams[index].Tersoff_lambda3;
  Xsi=0;
  dXsi_drn=0;
  Fc_ij = Fc_(r_ij,Tersoff_R,Tersoff_D);
  //Fc_ij = Fc[atomj*MaxNBond + i];
  dr_ij_drn = calc_dr_ij_drk(atomi,atomj,atomn,dir,r_ij,particle,Natoms,L);
  dFc_ij_drn = dFc_drn(atomi,atomj,atomn,r_ij,Tersoff_R,Tersoff_D,dr_ij_drn);
  
  for(k=0 ; k < NBond[atomj] ; k++){ // Go over all j neighbors k != i.
    atomk = Bonded[atomj*MaxNBond + k];
    if(atomk!=atomi){
      
      typek      = particle[atomk].type;
      index2     = typek*NAtomParams + typej;
      Tersoff_c  = BondParams[index2].Tersoff_c;
      Tersoff_d  = BondParams[index2].Tersoff_d;
      Tersoff_sqr_c  = BondParams[index2].Tersoff_sqr_c;
      Tersoff_sqr_d  = BondParams[index2].Tersoff_sqr_d;
      
      Tersoff_h  = BondParams[index2].Tersoff_h;
      Tersoff_R  = BondParams[index2].Tersoff_R;
      Tersoff_D  = BondParams[index2].Tersoff_D;
      //r_jk  = R_ij(particle,atomk,atomj, L);
      //r_jk      = rij[atomj*MaxNBond + k];
      r_jk = R_ij(particle,atomj,atomk, L);
      //r_jk=R_ij_(particle[atomj].r[0],particle[atomk].r[0],particle[atomj].r[1],particle[atomk].r[1],particle[atomj].r[2],particle[atomk].r[2],L);
      Fc_jk =  Fc_(r_jk,Tersoff_R,Tersoff_D);
      //Fc_jk = Fc[atomj*MaxNBond + k];
      dr_jk_drn = calc_dr_ij_drk(atomj,atomk,atomn,dir,r_jk,particle,Natoms,L);
      
      Theta_ijk = calc_val_angle(atomi,atomj,atomk,Natoms,NBond,Bonded,r_jk,r_ij,particle,L,rji,rjk,inside_acos);
      cos_ijk=cos(Theta_ijk);
      sqr_Tersoff_h_Cos_Theta = sqr(Tersoff_h-cos_ijk);
      
      dTheta_ijk_drn = calc_dthetaijk_drn(atomi,atomj,atomk,atomn,dir,particle,angle,BondNum,Natoms,rij,L,Theta_ijk,r_jk,r_ij,rji,rjk,cos_ijk);
      
      g        = 1+Tersoff_sqr_c/Tersoff_sqr_d-Tersoff_sqr_c/(Tersoff_sqr_d+sqr_Tersoff_h_Cos_Theta);
      
      
      dg_drn = Tersoff_sqr_c*pow(Tersoff_sqr_d + sqr_Tersoff_h_Cos_Theta,-2)*2*(Tersoff_h-cos_ijk)*sin(Theta_ijk)*dTheta_ijk_drn;
      
      dFc_jk_drn = dFc_drn(atomj,atomk,atomn,r_jk,Tersoff_R,Tersoff_D,dr_jk_drn);
      Xsi  += Fc_jk*g;
      dXsi_drn +=dFc_jk_drn*g + Fc_jk*dg_drn;
    }
  }
  
  //**
  //Xsi=1;
  //dXsi_drn=0;
  //**
  
  Fa      = -B*exp(-Tersoff_lambda2*r_ij);
  Fr      = A*exp(-Tersoff_lambda1*r_ij);
  b_ij    = pow(1+pow(beta*Xsi,Tersoff_n),-1/(2*Tersoff_n));
  dFa_drn = -Fa*Tersoff_lambda2*dr_ij_drn;
  dFr_drn = -Fr*Tersoff_lambda1*dr_ij_drn;
  
  if(Xsi==0)db_ij_drn=0;
  else db_ij_drn = -pow(beta,Tersoff_n)*pow(Xsi,Tersoff_n-1)/(2*pow(1+pow(beta*Xsi,Tersoff_n),1+1/(2*Tersoff_n)))*dXsi_drn;
    
  dE_Tersoff_ij_drn = dFc_ij_drn*(Fr + b_ij*Fa) + Fc_ij*(dFr_drn + db_ij_drn*Fa + b_ij*dFa_drn);
  /*
  cerr<<"atomi= "<<atomi<<"atomj= "<<atomj<<"atomn= "<<atomn<<endl;
  cerr<<"typei= "<<typei<<"typej= "<<typej<<endl;
  cerr<<"dE_Tersoff_ij_drn= "<<dE_Tersoff_ij_drn<<endl;
  cerr<<"dFc_ij_drn= "<<dFc_ij_drn<<endl;
  cerr<<"Fr= "<<Fr<<endl;
  cerr<<"Fa= "<<Fa<<endl;
  cerr<<"b_ij= "<<b_ij<<endl;
  cerr<<"Fc_ij= "<<Fc_ij<<endl;
  cerr<<"dFr_drn= "<<dFr_drn<<endl;
  cerr<<"db_ij_drn= "<<db_ij_drn<<endl;
  cerr<<"dFa_drn= "<<dFa_drn<<endl;
  cerr<<"Tersoff_lambda2= "<<Tersoff_lambda2<<endl;
  cerr<<"Tersoff_lambda1= "<<Tersoff_lambda1<<endl;
  cerr<<"dr_ij_drn= "<<dr_ij_drn<<endl;
  */
  if(dir==0 && (atomn=atomj))
    Etmp  = 0.5*Fc_ij*(Fr + b_ij*Fa);
  
  return(0.5*dE_Tersoff_ij_drn);
}

void calc_F_Tersoff(double *angle,int *BondNum,int Natoms,double *rij,point L,BondParamsStruct *BondParams,int *Bonded,atom *particle,double *Fc,point *FTersoff,int *NBond,double &ETersoff){

  ETersoff = 0.0;
  double ETersoff_=0.0;
  //#pragma omp parallel
  //{
    
    int i, j,k, atomk, atomi, atomj;
    //int hey=0;
    double Etmp;
    
    //#pragma omp for //reduction(+:ETersoff_)             
    for(atomk=0 ; atomk < Natoms ; atomk++){
      
      FTersoff[atomk].x = FTersoff[atomk].y = FTersoff[atomk].z = 0.0;
      for(i=0 ; i < NBond[atomk] ; i++){
	atomi = Bonded[atomk*MaxNBond + i];
	atomj = atomk;
	FTersoff[atomk].x -= calc_dE_Tersoff_ij_drn(atomi,atomj,atomk,0,i,angle,BondNum,Natoms,rij,L,BondParams,Bonded,particle,Fc,NBond,Etmp);
	//#pragma omp critical
	ETersoff_+=Etmp;
	FTersoff[atomk].y -= calc_dE_Tersoff_ij_drn(atomi,atomj,atomk,1,i,angle,BondNum,Natoms,rij,L,BondParams,Bonded,particle,Fc,NBond,Etmp);
	FTersoff[atomk].z -= calc_dE_Tersoff_ij_drn(atomi,atomj,atomk,2,i,angle,BondNum,Natoms,rij,L,BondParams,Bonded,particle,Fc,NBond,Etmp);
	//#pragma omp critical
	//ETersoff_+=Etmp;
	for(j=0 ; j < NBond[atomi] ; j++){
	  atomj=Bonded[atomi*MaxNBond + j];
	  FTersoff[atomk].x -= calc_dE_Tersoff_ij_drn(atomj,atomi,atomk,0,j,angle,BondNum,Natoms,rij,L,BondParams,Bonded,particle,Fc,NBond,Etmp);
	  FTersoff[atomk].y -= calc_dE_Tersoff_ij_drn(atomj,atomi,atomk,1,j,angle,BondNum,Natoms,rij,L,BondParams,Bonded,particle,Fc,NBond,Etmp);
	  FTersoff[atomk].z -= calc_dE_Tersoff_ij_drn(atomj,atomi,atomk,2,j,angle,BondNum,Natoms,rij,L,BondParams,Bonded,particle,Fc,NBond,Etmp);
	  
	}
      }
    }
    //}
  ETersoff=ETersoff_;
}



void read_file(char *filename)
{
  int i,j,k,l,limit;
  char s[MAXLINE];

  // REBO Parameters (AIREBO)

  double rcmin_CC,rcmin_CH,rcmin_HH,rcmax_CC,rcmax_CH,
    rcmax_HH,rcmaxp_CC,rcmaxp_CH,rcmaxp_HH;
  double Q_CC,Q_CH,Q_HH,alpha_CC,alpha_CH,alpha_HH,A_CC,A_CH,A_HH;
  double BIJc_CC1,BIJc_CC2,BIJc_CC3,BIJc_CH1,BIJc_CH2,BIJc_CH3,
    BIJc_HH1,BIJc_HH2,BIJc_HH3;
  double Beta_CC1,Beta_CC2,Beta_CC3,Beta_CH1,Beta_CH2,Beta_CH3,Beta_HH1,Beta_HH2,Beta_HH3;
  double rho_CC,rho_CH,rho_HH;

  // LJ Parameters (AIREBO)
  
  double rcLJmin_CC,rcLJmin_CH,rcLJmin_HH,rcLJmax_CC,rcLJmax_CH,
    rcLJmax_HH,bLJmin_CC;
  double bLJmin_CH,bLJmin_HH,bLJmax_CC,bLJmax_CH,bLJmax_HH,
    epsilon_CC,epsilon_CH,epsilon_HH;
  double sigma_CC,sigma_CH,sigma_HH,epsilonT_CCCC,epsilonT_CCCH,epsilonT_HCCH;
  
  // read file on proc 0
  
  
  //FILE *fp = open_potential(filename);
  FILE *fp;

  /*
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open AIREBO potential file %s",filename);
    error->one(FLERR,str);
  }
  */
  if ((fp = fopen(filename,"r"))== NULL) {
    char str[128];
    sprintf(str,"Cannot open AIREBO potental file %s",filename);
    exit(0);
  }
  
  // skip initial comment lines
  
  while (1) {
    fgets(s,MAXLINE,fp);
    if (s[0] != '#') break;
  }
    
    // read parameters
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmin_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmin_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmin_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmax_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmax_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmax_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmaxp_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmaxp_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmaxp_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&smin);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Nmin);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Nmax);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&NCmin);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&NCmax);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Q_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Q_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Q_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&alpha_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&alpha_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&alpha_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&A_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&A_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&A_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CC1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CC2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CC3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_HH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_HH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_HH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CC1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CC2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CC3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_HH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_HH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_HH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rho_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rho_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rho_HH);
    //cerr<<"BIJc_CH1 ="<<BIJc_CH1<<endl;
    //cerr<<"rcmin_CC ="<<rcmin_CC<<endl;
    // LJ parameters

    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmin_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmin_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmin_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmax_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmax_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmax_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmin_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmin_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmin_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmax_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmax_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmax_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilon_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilon_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilon_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&sigma_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&sigma_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&sigma_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilonT_CCCC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilonT_CCCH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilonT_HCCH);

    // gC spline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);

    // number-1 = # of domains for the spline

    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);
  
    for (i = 0; i < limit; i++) {
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lg",&gCdom[i]);
    }
    fgets(s,MAXLINE,fp);
    for (i = 0; i < limit-1; i++) {
      for (j = 0; j < 6; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&gC1[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    for (i = 0; i < limit-1; i++) {
      for (j = 0; j < 6; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&gC2[i][j]);
      }
    }
  
    // gH spline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);

    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);

    for (i = 0; i < limit; i++) {
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lg",&gHdom[i]);
    }

    fgets(s,MAXLINE,fp);

    for (i = 0; i < limit-1; i++) {
      for (j = 0; j < 6; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&gH[i][j]);
      }
    }

    // pCC spline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);

    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);

    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/2; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&pCCdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);

    for (i = 0; i < (int) pCCdom[0][1]; i++) {
      for (j = 0; j < (int) pCCdom[1][1]; j++) {
        for (k = 0; k < 16; k++) {
          fgets(s,MAXLINE,fp);
          sscanf(s,"%lg",&pCC[i][j][k]);
        }
      }
    }

    // pCH spline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);

    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/2; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&pCHdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);

    for (i = 0; i < (int) pCHdom[0][1]; i++) {
      for (j = 0; j < (int) pCHdom[1][1]; j++) {
        for (k = 0; k < 16; k++) {
          fgets(s,MAXLINE,fp);
          sscanf(s,"%lg",&pCH[i][j][k]);
        }
      }
    }

    // piCC cpline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);

    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);

    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&piCCdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);

    for (i = 0; i < (int) piCCdom[0][1]; i++) {
      for (j = 0; j < (int) piCCdom[1][1]; j++) {
        for (k = 0; k < (int) piCCdom[2][1]; k++) {
          for (l = 0; l < 64; l = l+1) {
            fgets(s,MAXLINE,fp);
            sscanf(s,"%lg",&piCC[i][j][k][l]);
          }
        }
      }
    }
  
    // piCH spline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);

    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);

    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&piCHdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);

    for (i = 0; i < (int) piCHdom[0][1]; i++) {
      for (j = 0; j < (int) piCHdom[1][1]; j++) {
        for (k = 0; k < (int) piCHdom[2][1]; k++) {
          for (l = 0; l < 64; l = l+1) {
            fgets(s,MAXLINE,fp);
            sscanf(s,"%lg",&piCH[i][j][k][l]);
          }
        }
      }
    }

    // piHH spline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);

    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);

    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&piHHdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);

    for (i = 0; i < (int) piHHdom[0][1]; i++) {
      for (j = 0; j < (int) piHHdom[1][1]; j++) {
        for (k = 0; k < (int) piHHdom[2][1]; k++) {
          for (l = 0; l < 64; l = l+1) {
            fgets(s,MAXLINE,fp);
            sscanf(s,"%lg",&piHH[i][j][k][l]);
          }
        }
      }
    }

    // Tij spline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);

    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);

    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&Tijdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);

    for (i = 0; i < (int) Tijdom[0][1]; i++) {
      for (j = 0; j < (int) Tijdom[1][1]; j++) {
        for (k = 0; k < (int) Tijdom[2][1]; k++) {
          for (l = 0; l < 64; l = l+1) {
            fgets(s,MAXLINE,fp);
            sscanf(s,"%lg",&Tijc[i][j][k][l]);
          }
        }
      }
    }

    fclose(fp);
  

  // store read-in values in arrays

  // REBO

    rcmin[0][0] = rcmin_CC;
    
    rcmin[0][1] = rcmin_CH;
    rcmin[1][0] = rcmin[0][1];
    rcmin[1][1] = rcmin_HH;

    rcmax[0][0] = rcmax_CC;
    rcmax[0][1] = rcmax_CH;
    rcmax[1][0] = rcmax[0][1];
    rcmax[1][1] = rcmax_HH;

    rcmaxsq[0][0] = rcmax[0][0]*rcmax[0][0];
    rcmaxsq[1][0] = rcmax[1][0]*rcmax[1][0];
    rcmaxsq[0][1] = rcmax[0][1]*rcmax[0][1];
    rcmaxsq[1][1] = rcmax[1][1]*rcmax[1][1];

    rcmaxp[0][0] = rcmaxp_CC;
    rcmaxp[0][1] = rcmaxp_CH;
    rcmaxp[1][0] = rcmaxp[0][1];
    rcmaxp[1][1] = rcmaxp_HH;

    Q[0][0] = Q_CC;
    Q[0][1] = Q_CH;
    Q[1][0] = Q[0][1];
    Q[1][1] = Q_HH;

    alpha[0][0] = alpha_CC;
    alpha[0][1] = alpha_CH;
    alpha[1][0] = alpha[0][1];
    alpha[1][1] = alpha_HH;

    A_ij[0][0] = A_CC;
    A_ij[0][1] = A_CH;
    A_ij[1][0] = A_ij[0][1];
    A_ij[1][1] = A_HH;

    rho[0][0] = rho_CC;
    rho[0][1] = rho_CH;
    rho[1][0] = rho[0][1];
    rho[1][1] = rho_HH;

    BIJc[0][0][0] = BIJc_CC1;
    BIJc[0][0][1] = BIJc_CC2;
    BIJc[0][0][2] = BIJc_CC3;
    BIJc[0][1][0] = BIJc_CH1;
    BIJc[0][1][1] = BIJc_CH2;
    BIJc[0][1][2] = BIJc_CH3;
    BIJc[1][0][0] = BIJc_CH1;
    BIJc[1][0][1] = BIJc_CH2;
    BIJc[1][0][2] = BIJc_CH3;
    BIJc[1][1][0] = BIJc_HH1;
    BIJc[1][1][1] = BIJc_HH2;
    BIJc[1][1][2] = BIJc_HH3;

    Beta[0][0][0] = Beta_CC1;
    Beta[0][0][1] = Beta_CC2;
    Beta[0][0][2] = Beta_CC3;
    Beta[0][1][0] = Beta_CH1;
    Beta[0][1][1] = Beta_CH2;
    Beta[0][1][2] = Beta_CH3;
    Beta[1][0][0] = Beta_CH1;
    Beta[1][0][1] = Beta_CH2;
    Beta[1][0][2] = Beta_CH3;
    Beta[1][1][0] = Beta_HH1;
    Beta[1][1][1] = Beta_HH2;
    Beta[1][1][2] = Beta_HH3;

    // LJ

    rcLJmin[0][0] = rcLJmin_CC;
    rcLJmin[0][1] = rcLJmin_CH;
    rcLJmin[1][0] = rcLJmin[0][1];
    rcLJmin[1][1] = rcLJmin_HH;

    rcLJmax[0][0] = rcLJmax_CC;
    rcLJmax[0][1] = rcLJmax_CH;
    rcLJmax[1][0] = rcLJmax[0][1];
    rcLJmax[1][1] = rcLJmax_HH;

    rcLJmaxsq[0][0] = rcLJmax[0][0]*rcLJmax[0][0];
    rcLJmaxsq[1][0] = rcLJmax[1][0]*rcLJmax[1][0];
    rcLJmaxsq[0][1] = rcLJmax[0][1]*rcLJmax[0][1];
    rcLJmaxsq[1][1] = rcLJmax[1][1]*rcLJmax[1][1];

    bLJmin[0][0] = bLJmin_CC;
    bLJmin[0][1] = bLJmin_CH;
    bLJmin[1][0] = bLJmin[0][1];
    bLJmin[1][1] = bLJmin_HH;

    bLJmax[0][0] = bLJmax_CC;
    bLJmax[0][1] = bLJmax_CH;
    bLJmax[1][0] = bLJmax[0][1];
    bLJmax[1][1] = bLJmax_HH;

    epsilon[0][0] = epsilon_CC;
    epsilon[0][1] = epsilon_CH;
    epsilon[1][0] = epsilon[0][1];
    epsilon[1][1] = epsilon_HH;

    sigma[0][0] = sigma_CC;
    sigma[0][1] = sigma_CH;
    sigma[1][0] = sigma[0][1];
    sigma[1][1] = sigma_HH;

    // torsional

    thmin = -1.0;
    thmax = -0.995;
    epsilonT[0][0] = epsilonT_CCCC;
    epsilonT[0][1] = epsilonT_CCCH;
    epsilonT[1][0] = epsilonT[0][1];
    epsilonT[1][1] = epsilonT_HCCH;
  


}
void scale_coords(int dir,int Natoms,point &L,atom *particle, double factor){
  int i;
  
  for(i=0;i<Natoms;i++) particle[i].r[dir] *= factor;
  
  if (dir==0)
    L.x *= factor;
  else if (dir==1)
    L.y *= factor;
  else
    L.z *= factor;
  //  *pbc *= factor;
  return;
}


void relax_cell(int Natoms,int *NBond,AtomParamsStruct *AtomParams,BondParamsStruct *BondParams,int *Bonded,double *Bmat,double *qvec,double *Avec,double *rij,double *rij_long,point &L, atom *particle,Normal_struct *Normal,double *Fc,FILE* logfile,double *angle,int *BondNum,int *Neighb,int *NNeighb,int potential_flag,int Interlayer,int *Normal_atom,int *NList,int *List,point *R0,int *NBond_once,int *Bonded_once,int BN_flag,int Ncar_hyd,int N_init,int N_init_ILP){
  double origEnergy,origEnergy_;
  double newEnergy,newEnergy_;
  double scaleup = 1.001;
  double scaledown = 1./scaleup;
  ///static const char message[]="\rrelax box size %c: %.3f    ";
  int nchange=0, maxchange=1000,dir;
  //try to increase the volume
  //printf(message,dir,pow(scaleup,++nchange)); fflush(stdout);

  calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rij_long,NList,List,R0,NBond_once,Bonded_once, potential_flag);

  if(Interlayer == 1)
    calc_Normal(Natoms,Normal_atom,particle,Normal,L);
  
  origEnergy = Calc_Potential(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,rij_long,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,potential_flag,origEnergy_,NBond_once,Bonded_once,BN_flag,Ncar_hyd,N_init,N_init_ILP);
  
  for(dir=0;dir<3;dir++){
    scale_coords(dir,Natoms,L,particle,scaleup);
    //newEnergy=E_interlayer()+E_intralayer();
    calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rij_long,NList,List,R0,NBond_once,Bonded_once, potential_flag);
    if(Interlayer == 1)
      calc_Normal(Natoms,Normal_atom,particle,Normal,L);
    
    newEnergy = Calc_Potential(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,rij_long,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,potential_flag,newEnergy_,NBond_once,Bonded_once,BN_flag,Ncar_hyd,N_init,N_init_ILP);
    
    if(newEnergy < origEnergy){ //keep increasing the volume
      origEnergy = newEnergy;
      while(nchange<maxchange){
	//printf(message,dir,pow(scaleup,++nchange)); fflush(stdout);
	nchange++;
	scale_coords(dir,Natoms,L,particle,scaleup);
	calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rij_long,NList,List,R0,NBond_once,Bonded_once, potential_flag);
	if(Interlayer == 1)
	  calc_Normal(Natoms,Normal_atom,particle,Normal,L);
	
	newEnergy = Calc_Potential(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,rij_long,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,potential_flag,newEnergy_,NBond_once,Bonded_once,BN_flag,Ncar_hyd,N_init,N_init_ILP);
	
	if(newEnergy < origEnergy){
	  origEnergy = newEnergy;
	}
	else{ //step back and exit
	  //printf(message,dir,pow(scaleup,--nchange)); fflush(stdout);
	  nchange--;
	  scale_coords(dir,Natoms,L,particle,scaledown);
	  break;
	}
      }
    }
    else { //try to reduce the volume
      origEnergy = newEnergy;
      while(nchange>-maxchange){
	//printf(message,dir,pow(scaleup,--nchange)); fflush(stdout);
	nchange--;
	scale_coords(dir,Natoms,L,particle,scaledown);
	calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rij_long,NList,List,R0,NBond_once,Bonded_once, potential_flag);
	if(Interlayer == 1)
	  calc_Normal(Natoms,Normal_atom,particle,Normal,L);
	
	newEnergy = Calc_Potential(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,rij_long,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,potential_flag,newEnergy_,NBond_once,Bonded_once,BN_flag,Ncar_hyd,N_init,N_init_ILP);
	
	if(newEnergy < origEnergy){
	  origEnergy = newEnergy;
	}
	else{ //step back and exit
	  //printf(message,dir,pow(scaleup,++nchange)); fflush(stdout);
	  nchange++;
	  scale_coords(dir,Natoms,L,particle,scaleup);
	  break;
	}
      }
    }
  }
  //printf("\n"); fflush(stdout);
  return;
}

void TORSION_(atom *particle,point *f,int *BondNum,int Natoms,point L,int *Bonded,int *NBond,double &E_TORSION_,int *NBond_once,int *Bonded_once,int  Ncar_hyd)
{
  double E_TORSION=0.0;

    int i,j,k,l,ii,inum;
    //tagint itag,jtag;
    int itag,jtag;
    double evdwl,fpair,xtmp,ytmp,ztmp;
    double cos321;
    double w21,dw21,cos234,w34,dw34;
    double cross321[3],cross321mag,cross234[3],cross234mag;
    double w23,dw23,cw2,ekijl,Ec;
    double cw,cwnum,cwnom;
    double rij,rij2,rik,rjl,tspjik,dtsjik,tspijl,dtsijl,costmp,fcpc;
    double sin321,sin234,rjk2,rik2,ril2,rjl2;
    double rjk,ril;
    double Vtors;
    double dndij[3],tmpvec[3],dndik[3],dndjl[3];
    double dcidij,dcidik,dcidjk,dcjdji,dcjdjl,dcjdil;
    double dsidij,dsidik,dsidjk,dsjdji,dsjdjl,dsjdil;
    double dxidij,dxidik,dxidjk,dxjdji,dxjdjl,dxjdil;
    double ddndij,ddndik,ddndjk,ddndjl,ddndil,dcwddn,dcwdn,dvpdcw,Ftmp[3];
    double del32[3],rsq,r32,del23[3],del21[3],r21;
    double deljk[3],del34[3],delil[3],delkl[3],r23,r34;
    double fi[3],fj[3],fk[3],fl[3];
    int itype,jtype,ktype,ltype,kk,ll,jj;
    int *ilist,*REBO_neighs_i,*REBO_neighs_j;
    
    //double **x = atom->x;
    //double **f = atom->f;
    //int *type = atom->type;
    //tagint *tag = atom->tag;
    
    //inum = list->inum;
    //ilist = list->ilist;
    //map[type[k]]
    //REBO_numneigh[i]
    //for (j = 0; j < NBond[atomi]; j++) {
    // atomj = Bonded[atomi*MaxNBond + j];
    // itag = atomj;
    //}
    
#pragma omp for //reduction(+:E_TORSION)
    for (ii = 0; ii < Ncar_hyd; ii++) {
      //cerr<<"hey1"<<endl;
      //i = ilist[ii];
      i=ii;
      //itag = tag[i];
      itag = ii;
      itype = particle[i].type;
      if (itype != 0) continue;
      xtmp = particle[i].r[0];
      ytmp = particle[i].r[1];
      ztmp = particle[i].r[2];
      //REBO_neighs_i = REBO_firstneigh[i];
      //cerr<<"hey2"<<endl;
      for (jj = 0; jj < NBond_once[ii]; jj++) {
	//j = REBO_neighs_i[jj];
	j = Bonded_once[i*MaxNBond + jj];
	//jtag = tag[j];
	jtag = j;
	/*
	  if (itag > jtag) {
	  if ((itag+jtag) % 2 == 0) continue;
	  } else if (itag < jtag) {
	  if ((itag+jtag) % 2 == 1) continue;
	  } else {
	  if (particle[j].r[2] < ztmp) continue;
	  if (particle[j].r[2] == ztmp && particle[j].r[1] < ytmp) continue;
	  if (particle[j].r[2] == ztmp && particle[j].r[1] == ytmp && particle[j].r[0] < xtmp) continue;
	  }
	*/
      jtype = particle[j].type;
      if (jtype != 0) continue;
      
      del32[0]  = R_PBC(particle,j,i,0,L.x);
      del32[1]  = R_PBC(particle,j,i,1,L.y);
      del32[2]  = R_PBC(particle,j,i,2,L.z);
      rsq = del32[0]*del32[0] + del32[1]*del32[1] + del32[2]*del32[2];
      r32 = sqrt(rsq);
      del23[0] = -del32[0];
      del23[1] = -del32[1];
      del23[2] = -del32[2];
      r23 = r32;
      w23 = Sp(r23,rcmin[itype][jtype],rcmax[itype][jtype],dw23);
      //cerr<<"hey3"<<endl;
      for (kk = 0; kk < NBond[i]; kk++) {
        //k = REBO_neighs_i[kk];
	k = Bonded[i*MaxNBond + kk];
        ktype = particle[k].type;
        if (k == j) continue;
	del21[0] = R_PBC(particle,i,k,0,L.x);
	del21[1] = R_PBC(particle,i,k,1,L.y);
	del21[2] = R_PBC(particle,i,k,2,L.z);
        rsq = del21[0]*del21[0] + del21[1]*del21[1] + del21[2]*del21[2];
        r21 = sqrt(rsq);
        cos321 = - ((del21[0]*del32[0]) + (del21[1]*del32[1]) +
                    (del21[2]*del32[2])) / (r21*r32);
        cos321 = MIN(cos321,1.0);
        cos321 = MAX(cos321,-1.0);
        sin321 = sqrt(1.0 - cos321*cos321);
        if (sin321 < TOL) continue;

        //deljk[0] = del21[0]-del23[0];
        //deljk[1] = del21[1]-del23[1];
        //deljk[2] = del21[2]-del23[2];
	
	deljk[0] = R_PBC(particle,j,k,0,L.x);
        deljk[1] = R_PBC(particle,j,k,1,L.y);
        deljk[2] = R_PBC(particle,j,k,2,L.z);
        rjk2 = deljk[0]*deljk[0] + deljk[1]*deljk[1] + deljk[2]*deljk[2];
        rjk=sqrt(rjk2);
        rik2 = r21*r21;
        w21 = Sp(r21,rcmin[itype][ktype],rcmax[itype][ktype],dw21);

        rij = r32;
        rik = r21;
        rij2 = r32*r32;
        rik2 = r21*r21;
        costmp = 0.5*(rij2+rik2-rjk2)/rij/rik;
        tspjik = Sp2(costmp,thmin,thmax,dtsjik);
        dtsjik = -dtsjik;
	//cerr<<"hey4"<<endl;
        //REBO_neighs_j = REBO_firstneigh[j];
        for (ll = 0; ll <NBond[j]; ll++) {
          //l = REBO_neighs_j[ll];
	  l = Bonded[j*MaxNBond + ll];
          ltype = particle[l].type;
          if (l == i || l == k) continue;
	  del34[0] = R_PBC(particle,j,l,0,L.x);
	  del34[1] = R_PBC(particle,j,l,1,L.y);
	  del34[2] = R_PBC(particle,j,l,2,L.z);
          rsq = del34[0]*del34[0] + del34[1]*del34[1] + del34[2]*del34[2];
          r34 = sqrt(rsq);
          cos234 = (del32[0]*del34[0] + del32[1]*del34[1] +
                    del32[2]*del34[2]) / (r32*r34);
          cos234 = MIN(cos234,1.0);
          cos234 = MAX(cos234,-1.0);
          sin234 = sqrt(1.0 - cos234*cos234);
          if (sin234 < TOL) continue;
          w34 = Sp(r34,rcmin[jtype][ltype],rcmax[jtype][ltype],dw34);
          //delil[0] = del23[0] + del34[0];
          //delil[1] = del23[1] + del34[1];
          //delil[2] = del23[2] + del34[2];
	  delil[0] = R_PBC(particle,i,l,0,L.x);
	  delil[1] = R_PBC(particle,i,l,1,L.y);
	  delil[2] = R_PBC(particle,i,l,2,L.z);
          ril2 = delil[0]*delil[0] + delil[1]*delil[1] + delil[2]*delil[2];
          ril=sqrt(ril2);
          rjl2 = r34*r34;

          rjl = r34;
          rjl2 = r34*r34;
          costmp = 0.5*(rij2+rjl2-ril2)/rij/rjl;
          tspijl = Sp2(costmp,thmin,thmax,dtsijl);
          dtsijl = -dtsijl; //need minus sign
          cross321[0] = (del32[1]*del21[2])-(del32[2]*del21[1]);
          cross321[1] = (del32[2]*del21[0])-(del32[0]*del21[2]);
          cross321[2] = (del32[0]*del21[1])-(del32[1]*del21[0]);
          cross321mag = sqrt(cross321[0]*cross321[0]+
                             cross321[1]*cross321[1]+
                             cross321[2]*cross321[2]);
          cross234[0] = (del23[1]*del34[2])-(del23[2]*del34[1]);
          cross234[1] = (del23[2]*del34[0])-(del23[0]*del34[2]);
          cross234[2] = (del23[0]*del34[1])-(del23[1]*del34[0]);
          cross234mag = sqrt(cross234[0]*cross234[0]+
                             cross234[1]*cross234[1]+
                             cross234[2]*cross234[2]);
          cwnum = (cross321[0]*cross234[0]) +
            (cross321[1]*cross234[1])+(cross321[2]*cross234[2]);
          cwnom = r21*r34*r32*r32*sin321*sin234;
          cw = cwnum/cwnom;

          cw2 = (.5*(1.0-cw));
          ekijl = epsilonT[ktype][ltype];
          Ec = 256.0*ekijl/405.0;
          Vtors = (Ec*(powint(cw2,5)))-(ekijl/10.0);

          //if (eflag) evdwl = Vtors*w21*w23*w34*(1.0-tspjik)*(1.0-tspijl);

	  //evdwl =

	  E_TORSION += Vtors*w21*w23*w34*(1.0-tspjik)*(1.0-tspijl);
	  //printf("E_TORSION= %.16f\n", E_TORSION);
          dndij[0] = (cross234[1]*del21[2])-(cross234[2]*del21[1]);
          dndij[1] = (cross234[2]*del21[0])-(cross234[0]*del21[2]);
          dndij[2] = (cross234[0]*del21[1])-(cross234[1]*del21[0]);

          tmpvec[0] = (del34[1]*cross321[2])-(del34[2]*cross321[1]);
          tmpvec[1] = (del34[2]*cross321[0])-(del34[0]*cross321[2]);
          tmpvec[2] = (del34[0]*cross321[1])-(del34[1]*cross321[0]);

          dndij[0] = dndij[0]+tmpvec[0];
          dndij[1] = dndij[1]+tmpvec[1];
          dndij[2] = dndij[2]+tmpvec[2];

          dndik[0] = (del23[1]*cross234[2])-(del23[2]*cross234[1]);
          dndik[1] = (del23[2]*cross234[0])-(del23[0]*cross234[2]);
          dndik[2] = (del23[0]*cross234[1])-(del23[1]*cross234[0]);

          dndjl[0] = (cross321[1]*del23[2])-(cross321[2]*del23[1]);
          dndjl[1] = (cross321[2]*del23[0])-(cross321[0]*del23[2]);
          dndjl[2] = (cross321[0]*del23[1])-(cross321[1]*del23[0]);

          dcidij = ((r23*r23)-(r21*r21)+(rjk*rjk))/(2.0*r23*r23*r21);
          dcidik = ((r21*r21)-(r23*r23)+(rjk*rjk))/(2.0*r23*r21*r21);
          dcidjk = (-rjk)/(r23*r21);
          dcjdji = ((r23*r23)-(r34*r34)+(ril*ril))/(2.0*r23*r23*r34);
          dcjdjl = ((r34*r34)-(r23*r23)+(ril*ril))/(2.0*r23*r34*r34);
          dcjdil = (-ril)/(r23*r34);

          dsidij = (-cos321/sin321)*dcidij;
          dsidik = (-cos321/sin321)*dcidik;
          dsidjk = (-cos321/sin321)*dcidjk;

          dsjdji = (-cos234/sin234)*dcjdji;
          dsjdjl = (-cos234/sin234)*dcjdjl;
          dsjdil = (-cos234/sin234)*dcjdil;

          dxidij = (r21*sin321)+(r23*r21*dsidij);
          dxidik = (r23*sin321)+(r23*r21*dsidik);
          dxidjk = (r23*r21*dsidjk);

          dxjdji = (r34*sin234)+(r23*r34*dsjdji);
          dxjdjl = (r23*sin234)+(r23*r34*dsjdjl);
          dxjdil = (r23*r34*dsjdil);

          ddndij = (dxidij*cross234mag)+(cross321mag*dxjdji);
          ddndik = dxidik*cross234mag;
          ddndjk = dxidjk*cross234mag;
          ddndjl = cross321mag*dxjdjl;
          ddndil = cross321mag*dxjdil;
          dcwddn = -cwnum/(cwnom*cwnom);
          dcwdn = 1.0/cwnom;
          dvpdcw = (-1.0)*Ec*(-.5)*5.0*powint(cw2,4) *
            w23*w21*w34*(1.0-tspjik)*(1.0-tspijl);

          Ftmp[0] = dvpdcw*((dcwdn*dndij[0])+(dcwddn*ddndij*del23[0]/r23));
          Ftmp[1] = dvpdcw*((dcwdn*dndij[1])+(dcwddn*ddndij*del23[1]/r23));
          Ftmp[2] = dvpdcw*((dcwdn*dndij[2])+(dcwddn*ddndij*del23[2]/r23));
          fi[0] = Ftmp[0];
          fi[1] = Ftmp[1];
          fi[2] = Ftmp[2];
          fj[0] = -Ftmp[0];
          fj[1] = -Ftmp[1];
          fj[2] = -Ftmp[2];

          Ftmp[0] = dvpdcw*((dcwdn*dndik[0])+(dcwddn*ddndik*del21[0]/r21));
          Ftmp[1] = dvpdcw*((dcwdn*dndik[1])+(dcwddn*ddndik*del21[1]/r21));
          Ftmp[2] = dvpdcw*((dcwdn*dndik[2])+(dcwddn*ddndik*del21[2]/r21));
          fi[0] += Ftmp[0];
          fi[1] += Ftmp[1];
          fi[2] += Ftmp[2];
          fk[0] = -Ftmp[0];
          fk[1] = -Ftmp[1];
          fk[2] = -Ftmp[2];

          Ftmp[0] = (dvpdcw*dcwddn*ddndjk*deljk[0])/rjk;
          Ftmp[1] = (dvpdcw*dcwddn*ddndjk*deljk[1])/rjk;
          Ftmp[2] = (dvpdcw*dcwddn*ddndjk*deljk[2])/rjk;
          fj[0] += Ftmp[0];
          fj[1] += Ftmp[1];
          fj[2] += Ftmp[2];
          fk[0] -= Ftmp[0];
          fk[1] -= Ftmp[1];
          fk[2] -= Ftmp[2];

          Ftmp[0] = dvpdcw*((dcwdn*dndjl[0])+(dcwddn*ddndjl*del34[0]/r34));
          Ftmp[1] = dvpdcw*((dcwdn*dndjl[1])+(dcwddn*ddndjl*del34[1]/r34));
          Ftmp[2] = dvpdcw*((dcwdn*dndjl[2])+(dcwddn*ddndjl*del34[2]/r34));
          fj[0] += Ftmp[0];
          fj[1] += Ftmp[1];
          fj[2] += Ftmp[2];
          fl[0] = -Ftmp[0];
          fl[1] = -Ftmp[1];
          fl[2] = -Ftmp[2];

          Ftmp[0] = (dvpdcw*dcwddn*ddndil*delil[0])/ril;
          Ftmp[1] = (dvpdcw*dcwddn*ddndil*delil[1])/ril;
          Ftmp[2] = (dvpdcw*dcwddn*ddndil*delil[2])/ril;
          fi[0] += Ftmp[0];
          fi[1] += Ftmp[1];
          fi[2] += Ftmp[2];
          fl[0] -= Ftmp[0];
          fl[1] -= Ftmp[1];
          fl[2] -= Ftmp[2];

          // coordination forces

          fpair = Vtors*dw21*w23*w34*(1.0-tspjik)*(1.0-tspijl) / r21;
          fi[0] -= del21[0]*fpair;
          fi[1] -= del21[1]*fpair;
          fi[2] -= del21[2]*fpair;
          fk[0] += del21[0]*fpair;
          fk[1] += del21[1]*fpair;
          fk[2] += del21[2]*fpair;

          fpair = Vtors*w21*dw23*w34*(1.0-tspjik)*(1.0-tspijl) / r23;
          fi[0] -= del23[0]*fpair;
          fi[1] -= del23[1]*fpair;
          fi[2] -= del23[2]*fpair;
          fj[0] += del23[0]*fpair;
          fj[1] += del23[1]*fpair;
          fj[2] += del23[2]*fpair;

          fpair = Vtors*w21*w23*dw34*(1.0-tspjik)*(1.0-tspijl) / r34;
          fj[0] -= del34[0]*fpair;
          fj[1] -= del34[1]*fpair;
          fj[2] -= del34[2]*fpair;
          fl[0] += del34[0]*fpair;
          fl[1] += del34[1]*fpair;
          fl[2] += del34[2]*fpair;

          // additional cut off function forces

          fcpc = -Vtors*w21*w23*w34*dtsjik*(1.0-tspijl);
          fpair = fcpc*dcidij/rij;
          fi[0] += fpair*del23[0];
          fi[1] += fpair*del23[1];
          fi[2] += fpair*del23[2];
          fj[0] -= fpair*del23[0];
          fj[1] -= fpair*del23[1];
          fj[2] -= fpair*del23[2];

          fpair = fcpc*dcidik/rik;
          fi[0] += fpair*del21[0];
          fi[1] += fpair*del21[1];
          fi[2] += fpair*del21[2];
          fk[0] -= fpair*del21[0];
          fk[1] -= fpair*del21[1];
          fk[2] -= fpair*del21[2];

          fpair = fcpc*dcidjk/rjk;
          fj[0] += fpair*deljk[0];
          fj[1] += fpair*deljk[1];
          fj[2] += fpair*deljk[2];
          fk[0] -= fpair*deljk[0];
          fk[1] -= fpair*deljk[1];
          fk[2] -= fpair*deljk[2];

          fcpc = -Vtors*w21*w23*w34*(1.0-tspjik)*dtsijl;
          fpair = fcpc*dcjdji/rij;
          fi[0] += fpair*del23[0];
          fi[1] += fpair*del23[1];
          fi[2] += fpair*del23[2];
          fj[0] -= fpair*del23[0];
          fj[1] -= fpair*del23[1];
          fj[2] -= fpair*del23[2];

          fpair = fcpc*dcjdjl/rjl;
          fj[0] += fpair*del34[0];
          fj[1] += fpair*del34[1];
          fj[2] += fpair*del34[2];
          fl[0] -= fpair*del34[0];
          fl[1] -= fpair*del34[1];
          fl[2] -= fpair*del34[2];

          fpair = fcpc*dcjdil/ril;
          fi[0] += fpair*delil[0];
          fi[1] += fpair*delil[1];
          fi[2] += fpair*delil[2];
          fl[0] -= fpair*delil[0];
          fl[1] -= fpair*delil[1];
          fl[2] -= fpair*delil[2];

          // sum per-atom forces into atom force array

          f[i].x += fi[0]; f[i].y += fi[1]; f[i].z += fi[2];
          f[j].x += fj[0]; f[j].y += fj[1]; f[j].z += fj[2];
          f[k].x += fk[0]; f[k].y += fk[1]; f[k].z += fk[2];
          f[l].x += fl[0]; f[l].y += fl[1]; f[l].z += fl[2];

          //if (evflag) {
	  // delkl[0] = delil[0] - del21[0];
	  // delkl[1] = delil[1] - del21[1];
	  // delkl[2] = delil[2] - del21[2];
	  // ev_tally4(i,j,k,l,evdwl,fi,fj,fk,delil,del34,delkl);
	  // }
	  
        }
      }
    }
      //}
  }
  //#pragma omp critical
  E_TORSION_=E_TORSION;
}

/****************************************************************/
double calc_dP_ij_dr(int atomi,int atomj, int atomk,int dir,atom *particle,int Natoms,double *rij,double r_ij,Normal_struct *Normal,point L,int *Normal_atom,int Interlayer,double rx,double ry,double rz,double inside_acos,double acos_inside_acos,double sqrt_inside,double sin_angle,double cos_angle,double Pij_term,double dR_ij,dr_ij *dN)
{
  
  double angle,dPij,dinside_acos,dangle;
  //double dR_ij;
    
  //dR_ij = calc_dr_ij_drk(atomi,atomj,atomk,dir,r_ij,particle,Natoms,L);
  
  //rx = particle[atomi].r[0] - particle[atomj].r[0];
  //ry = particle[atomi].r[1] - particle[atomj].r[1];
  //rz = particle[atomi].r[2] - particle[atomj].r[2];
  
  //rx -= L.x * rint(rx / (L.x + TINY));
  //ry -= L.y * rint(ry / (L.y + TINY));
  //rz -= L.z * rint(rz / (L.z + TINY));
  
  //inside_acos = (rx*Normal[atomi].x + ry*Normal[atomi].y + rz*Normal[atomi].z) / (r_ij);
  if (fabs(inside_acos) < 1.0){
    if (atomk == atomi){
      if (dir == 0)dinside_acos = (Normal[atomi].x)/(r_ij) - Pij_term*dR_ij;
      if (dir == 1)dinside_acos = (Normal[atomi].y)/(r_ij) - Pij_term*dR_ij;
      if (dir == 2)dinside_acos = (Normal[atomi].z)/(r_ij) - Pij_term*dR_ij;
    }
    else if (atomk == atomj){
      if (dir == 0)dinside_acos = (-Normal[atomi].x)/(r_ij) - Pij_term*dR_ij;
      if (dir == 1)dinside_acos = (-Normal[atomi].y)/(r_ij) - Pij_term*dR_ij;
      if (dir == 2)dinside_acos = (-Normal[atomi].z)/(r_ij) - Pij_term*dR_ij;
    }
    else{
      if (Interlayer == 1)
	calc_dNormal_k(atomi,atomk,dir,dN,particle,Normal,Normal_atom,L);
      dinside_acos = (rx*dN[0].r[0] + ry*dN[0].r[1] + rz*dN[0].r[2])/r_ij ;
    }
    dangle = -1.0/sqrt_inside*dinside_acos;
    dPij = dR_ij*sin_angle + r_ij*cos_angle*dangle;  
  }
  else {
    dPij = 0.0;
  }
    
  return(dPij);
  //return(dN[0].r[0]+dN[0].r[1]+dN[0].r[2]);
}

//*********************************

void calc_F_KC(point *FvdW, int Natoms, atom *particle, double *rij_long, BondParamsStruct *BondParams, point L,Normal_struct *Normal,int *Normal_atom,int Interlayer,int *Neighb,int *NNeighb,dr_ij *dN,double &ERep_,int Ncar_hyd)
{
   ERep_=0;
  
  double ERep=0.0;

    
    int atomk, atomi, typek, typei,atomj,typej, index,j;
    double alpha, rvdW, Epsilon, term1, dEvdW_drij;
    double dTap, r_ij, Tap,Epsilon_C_vdW_exp,term3,term4,term5,term7;
    double dP_ij_drk,dP_ji_drk,R,Pij,Pji,gamma,dr_ij_drk,exp_term1;
    double exp_Pij_gamma, exp_Pji_gamma,Pij_gamma,Pji_gamma,reff,r_ik,r_jk;
    double rx_ij,ry_ij,rz_ij,rx_ji,ry_ji,rz_ji,inside_acos_ij,inside_acos_ji,Pij_term1,Pji_term1,Pij_term,Pji_term,dr_ji_drk;
    double angle_ij, sqrt_inside_ij, sin_angle_ij,cos_angle_ij,angle_ji,sqrt_inside_ji,sin_angle_ji,cos_angle_ji,dinside_acos,dangle;
    double C0_kc,C2_kc,C4_kc,A_kc,fji,fij,term6,dfji,dfij,A_term;
    double ERep_Tmp,dP_ij_term,dP_ji_term,d_Norm_length,Pow_normal_i,Pow_normal_j;
    double temp_b,temp_c,Nori_L,Norj_L,Fasolina;
    int n,dir,dpij_flag,dpji_flag;
 
    dP_ij_drk = 0.0;
    dP_ji_drk = 0.0;

    C0_kc      = 15.71/1000;
    C2_kc      = 12.29/1000;
    C4_kc      = 4.933/1000;
    alpha   = 3.629;
    Epsilon = 3.03/1000;
    R       = 3.34;
    gamma   = 0.578;
    A_kc    = 10.238/1000;
    
#pragma omp for schedule(dynamic,10) //reduction(+:ERep)
    for(atomi=0 ; atomi < Ncar_hyd ; atomi++){
      typei = particle[atomi].type;
      if(typei==1)continue;
      for (j=0 ; j < NNeighb[atomi]; j++){
	atomk = atomj = Neighb[atomi*MaxNeighb + j] ;	  
	typek = typej = particle[atomj].type;	
	if(typek==1)continue;

	if(typei==0&&typej==0)Fasolina=1.0;
	else if((typei==0 && typej==7) || typei==7 && typej==0)Fasolina=2.0;
	else if((typei==0 && typej==5) || typei==5 && typej==0)Fasolina=0.6;
	else{
	  cerr<<"wrong atom type in F_KC atomi= "<<atomi<<"atomj= "<<atomj<< " exiting"<<endl;
	  exit(0);
	}
	
	index = typej*NAtomParams + typei;
	//r_ij  = R_ij(particle,atomi,atomj, L);
	//r_ij = rij_long[atomi*MaxNeighb + j];
	rx_ij = R_PBC(particle,atomi,atomj,0,L.x);
	ry_ij = R_PBC(particle,atomi,atomj,1,L.y);	
	rz_ij = R_PBC(particle,atomi,atomj,2,L.z);
	
	rx_ji = -rx_ij;
	ry_ji = -ry_ij;
	rz_ji = -rz_ij;

	r_ij = sqrt(rx_ij*rx_ij + ry_ij*ry_ij + rz_ij*rz_ij);
	//r_ij = R_ij_(particle[atomi].r[0],particle[atomj].r[0],particle[atomi].r[1],particle[atomj].r[1],particle[atomi].r[2],particle[atomj].r[2],L);
	Tap   = calc_Tap(r_ij);
	dTap  = calc_dTap(r_ij);  
	
	reff    = BondParams[index].r_eff;	    
	
	dP_ij_drk=0.0;
	dP_ji_drk=0.0;
	Pij_term1=0.0;
	Pji_term1=0.0;
	dN[0].r[0]=0;
	dN[0].r[1]=0;
	dN[0].r[2]=0;
	dpij_flag = 1;
	dpji_flag = 1;
	
	
	
	Nori_L=rx_ij*Normal[atomi].x + ry_ij*Normal[atomi].y + rz_ij*Normal[atomi].z;
	Norj_L=rx_ji*Normal[atomj].x + ry_ji*Normal[atomj].y + rz_ji*Normal[atomj].z;

	inside_acos_ij = (Nori_L) / (r_ij);
	inside_acos_ji = (Norj_L) / (r_ij);
	
	if(((r_ij*r_ij-sqr(Nori_L)) <= 0) || Interlayer==2){
	  Pij=0.0+TINY;
	  dpij_flag = 0;
	}
	else Pij=sqrt(r_ij*r_ij - sqr(Nori_L));
	
	if(((r_ij*r_ij-sqr(Norj_L)) <= 0) || Interlayer==2){
	  Pji=0.0+TINY;
	  dpji_flag = 0;
	}
	else Pji=sqrt(r_ij*r_ij - sqr(Norj_L));
	//****
	
	Pij_gamma = Pij/gamma*dpij_flag; //new
	Pji_gamma = Pji/gamma*dpji_flag; //new
	exp_Pij_gamma = exp(-pow((Pij_gamma),2.0)); //new
	exp_Pji_gamma = exp(-pow((Pji_gamma),2.0)); //new
	exp_term1 = exp(-alpha*(r_ij-R));
	
	fij = C0_kc+ C2_kc*pow((Pij_gamma),2.0)+ C4_kc*pow((Pij_gamma),4.0);
	fji = C0_kc+ C2_kc*pow((Pji_gamma),2.0)+ C4_kc*pow((Pji_gamma),4.0);
	
	fij*=exp_Pij_gamma;
	fji*=exp_Pji_gamma;
	A_term=A_kc*pow((r_ij/R),-6);
	Epsilon_C_vdW_exp = Epsilon + fij + fji;
	ERep_Tmp = (exp_term1*(Epsilon_C_vdW_exp) - A_term);
	ERep += Fasolina*Tap*(ERep_Tmp);
	
	dfij    = 2.0*C2_kc*Pij_gamma/gamma+4.0*C4_kc*pow(Pij_gamma,3)/gamma;
        dfji    = 2.0*C2_kc*Pji_gamma/gamma+4.0*C4_kc*pow(Pji_gamma,3)/gamma;
	//***
	term3   = exp_term1*(-alpha)*(Epsilon_C_vdW_exp);
	//term4   = exp_term1*((-2.0)*Pij_gamma/gamma*fij+exp_Pij_gamma*dfij)*dpij_flag;
	//term5   = exp_term1*((-2.0)*Pji_gamma/gamma*fji+exp_Pji_gamma*dfji)*dpji_flag;
	term4   = exp_term1*((-2.0)*Pij_gamma/(gamma)*fij+exp_Pij_gamma*dfij)/Pij*dpij_flag;
	term5   = exp_term1*((-2.0)*Pji_gamma/(gamma)*fji+exp_Pji_gamma*dfji)/Pji*dpji_flag;
	term6   = 6*A_kc*pow(r_ij/R,-7)/R;

	//********dirivative in the x direction******
	
	dir = 0;
	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	
	dP_ij_drk=r_ij*dr_ij_drk+Nori_L*Normal[atomi].x;
	dP_ji_drk=r_ij*dr_ij_drk-Norj_L*Normal[atomj].x;
	
	dEvdW_drij = Fasolina*(Tap*(term3*dr_ij_drk + term4*dP_ij_drk + term5*dP_ji_drk + term6*dr_ij_drk) + dTap*(ERep_Tmp)*dr_ij_drk);
	
	FvdW[atomk].x -= dEvdW_drij;
	FvdW[atomi].x += dEvdW_drij;
	
	//********dirivative in the y direction******
	
	dir = 1;
	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	
	dP_ij_drk=r_ij*dr_ij_drk+Nori_L*Normal[atomi].y;
	dP_ji_drk=r_ij*dr_ij_drk-Norj_L*Normal[atomj].y;
	
	dEvdW_drij = Fasolina*(Tap*(term3*dr_ij_drk + term4*dP_ij_drk + term5*dP_ji_drk + term6*dr_ij_drk) + dTap*(ERep_Tmp)*dr_ij_drk);
	FvdW[atomk].y -= dEvdW_drij;
	FvdW[atomi].y += dEvdW_drij;

	//********dirivative in the z direction******
	
	dir = 2;
	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	
	dP_ij_drk=r_ij*dr_ij_drk+Nori_L*Normal[atomi].z;
	dP_ji_drk=r_ij*dr_ij_drk-Norj_L*Normal[atomj].z;
	
	dEvdW_drij = Fasolina*(Tap*(term3*dr_ij_drk + term4*dP_ij_drk + term5*dP_ji_drk + term6*dr_ij_drk) + dTap*(ERep_Tmp)*dr_ij_drk);
	
	FvdW[atomk].z -= dEvdW_drij;
	FvdW[atomi].z += dEvdW_drij;
	
	if(((r_ij*r_ij-sqr(Nori_L)) <= 0) || Interlayer==2);
	else{
	  
	  temp_b = Tap*term4*(-Nori_L)*( (ry_ij*Normal[atomi].vector_2z - rz_ij*Normal[atomi].vector_2y) - Nori_L*( Normal[atomi].y*Normal[atomi].vector_2z - Normal[atomi].z*Normal[atomi].vector_2y))/(Normal[atomi].Normal_length);
	  temp_c = Tap*term4*(-Nori_L)*((-ry_ij*Normal[atomi].vector_1z + rz_ij*Normal[atomi].vector_1y) - Nori_L*(-Normal[atomi].y*Normal[atomi].vector_1z + Normal[atomi].z*Normal[atomi].vector_1y))/(Normal[atomi].Normal_length);
	  

	  FvdW[Normal_atom[atomi*3 + 0]].x+=Fasolina*( temp_b + temp_c );

	  FvdW[Normal_atom[atomi*3 + 1]].x-=Fasolina*temp_b;

	  FvdW[Normal_atom[atomi*3 + 2]].x-=Fasolina*temp_c;
	  
	  temp_b = Tap*term4*(-Nori_L)*((-rx_ij*Normal[atomi].vector_2z + rz_ij*Normal[atomi].vector_2x) - Nori_L*(-Normal[atomi].x*Normal[atomi].vector_2z + Normal[atomi].z*Normal[atomi].vector_2x))/(Normal[atomi].Normal_length);
	  temp_c = Tap*term4*(-Nori_L)*(( rx_ij*Normal[atomi].vector_1z - rz_ij*Normal[atomi].vector_1x) - Nori_L*( Normal[atomi].x*Normal[atomi].vector_1z - Normal[atomi].z*Normal[atomi].vector_1x))/(Normal[atomi].Normal_length);
	  

	  FvdW[Normal_atom[atomi*3 + 0]].y+=Fasolina*( temp_b + temp_c );

	  FvdW[Normal_atom[atomi*3 + 1]].y-=Fasolina*temp_b;

	  FvdW[Normal_atom[atomi*3 + 2]].y-=Fasolina*temp_c;
	  
	  temp_b = Tap*term4*(-Nori_L)*(( rx_ij*Normal[atomi].vector_2y - ry_ij*Normal[atomi].vector_2x) - Nori_L*( Normal[atomi].x*Normal[atomi].vector_2y - Normal[atomi].y*Normal[atomi].vector_2x))/(Normal[atomi].Normal_length);
	  temp_c = Tap*term4*(-Nori_L)*((-rx_ij*Normal[atomi].vector_1y + ry_ij*Normal[atomi].vector_1x) - Nori_L*(-Normal[atomi].x*Normal[atomi].vector_1y + Normal[atomi].y*Normal[atomi].vector_1x))/(Normal[atomi].Normal_length);
	  

	  FvdW[Normal_atom[atomi*3 + 0]].z+=Fasolina*( temp_b + temp_c );

	  FvdW[Normal_atom[atomi*3 + 1]].z-=Fasolina*temp_b;

	  FvdW[Normal_atom[atomi*3 + 2]].z-=Fasolina*temp_c;
	}
	//****************************************
	//****************************************
	
	if(((r_ij*r_ij-sqr(Norj_L)) <= 0) || Interlayer==2);
	else{
	  
	  temp_b = Tap*term5*(-Norj_L)*(( ry_ji*Normal[atomj].vector_2z - rz_ji*Normal[atomj].vector_2y) - Norj_L*( Normal[atomj].y*Normal[atomj].vector_2z - Normal[atomj].z*Normal[atomj].vector_2y))/(Normal[atomj].Normal_length);
	  temp_c = Tap*term5*(-Norj_L)*((-ry_ji*Normal[atomj].vector_1z + rz_ji*Normal[atomj].vector_1y) - Norj_L*(-Normal[atomj].y*Normal[atomj].vector_1z + Normal[atomj].z*Normal[atomj].vector_1y))/(Normal[atomj].Normal_length);

	  FvdW[Normal_atom[atomj*3 + 0]].x+=Fasolina*( temp_b + temp_c );

	  FvdW[Normal_atom[atomj*3 + 1]].x-=Fasolina*temp_b;

	  FvdW[Normal_atom[atomj*3 + 2]].x-=Fasolina*temp_c;
	  
	  temp_b = Tap*term5*(-Norj_L)*((-rx_ji*Normal[atomj].vector_2z + rz_ji*Normal[atomj].vector_2x) - Norj_L*(-Normal[atomj].x*Normal[atomj].vector_2z + Normal[atomj].z*Normal[atomj].vector_2x))/(Normal[atomj].Normal_length);
	  temp_c = Tap*term5*(-Norj_L)*(( rx_ji*Normal[atomj].vector_1z - rz_ji*Normal[atomj].vector_1x) - Norj_L*( Normal[atomj].x*Normal[atomj].vector_1z - Normal[atomj].z*Normal[atomj].vector_1x))/(Normal[atomj].Normal_length);

	  FvdW[Normal_atom[atomj*3 + 0]].y+=Fasolina*( temp_b + temp_c );

	  FvdW[Normal_atom[atomj*3 + 1]].y-=Fasolina*temp_b;

	  FvdW[Normal_atom[atomj*3 + 2]].y-=Fasolina*temp_c;
	  
	  temp_b = Tap*term5*(-Norj_L)*(( rx_ji*Normal[atomj].vector_2y - ry_ji*Normal[atomj].vector_2x) - Norj_L*( Normal[atomj].x*Normal[atomj].vector_2y - Normal[atomj].y*Normal[atomj].vector_2x))/(Normal[atomj].Normal_length);
	  temp_c = Tap*term5*(-Norj_L)*((-rx_ji*Normal[atomj].vector_1y + ry_ji*Normal[atomj].vector_1x) - Norj_L*(-Normal[atomj].x*Normal[atomj].vector_1y + Normal[atomj].y*Normal[atomj].vector_1x))/(Normal[atomj].Normal_length);

	  FvdW[Normal_atom[atomj*3 + 0]].z+=Fasolina*( temp_b + temp_c );

	  FvdW[Normal_atom[atomj*3 + 1]].z-=Fasolina*temp_b;

	  FvdW[Normal_atom[atomj*3 + 2]].z-=Fasolina*temp_c;
	}
      }
    }
  
  ERep_=ERep; 
}
double calc_dTap(double r_ij)
{
  double dTap;
  
  dTap = 7*Tap_7 * r_ij + 6*Tap_6;
  dTap = dTap * r_ij + 5*Tap_5;
  dTap = dTap * r_ij + 4*Tap_4;
  dTap = dTap * r_ij + 3*Tap_3;
  dTap = dTap * r_ij + 2*Tap_2;
  dTap = dTap * r_ij + Tap_1;
    
  return(dTap);
}
double calc_Tap(double r_ij)
{
  double Tap;
  
  Tap = 0.0;
  
  Tap = Tap_7 * r_ij + Tap_6;
  Tap = Tap * r_ij   + Tap_5;
  Tap = Tap * r_ij   + Tap_4;
  Tap = Tap * r_ij   + Tap_3;
  Tap = Tap * r_ij   + Tap_2;
  Tap = Tap * r_ij   + Tap_1;
  Tap = Tap * r_ij   + Tap_0;
  
  return(Tap);
}
double calc_dTap_Rep(double r_ij)
{
  double dTap;
  
  dTap = 7*Tap_7_Rep * r_ij + 6*Tap_6_Rep;
  dTap = dTap * r_ij + 5*Tap_5_Rep;
  dTap = dTap * r_ij + 4*Tap_4_Rep;
  dTap = dTap * r_ij + 3*Tap_3_Rep;
  dTap = dTap * r_ij + 2*Tap_2_Rep;
  dTap = dTap * r_ij + Tap_1_Rep;
    
  return(dTap);
}
double calc_Tap_Rep(double r_ij)
{
  double Tap;
  
  Tap = 0.0;
  
  Tap = Tap_7_Rep * r_ij + Tap_6_Rep;
  Tap = Tap * r_ij   + Tap_5_Rep;
  Tap = Tap * r_ij   + Tap_4_Rep;
  Tap = Tap * r_ij   + Tap_3_Rep;
  Tap = Tap * r_ij   + Tap_2_Rep;
  Tap = Tap * r_ij   + Tap_1_Rep;
  Tap = Tap * r_ij   + Tap_0_Rep;
  
  return(Tap);
}
void spline_init()
{
  int i,j,k;

  for (i = 0; i < 5; i++) {
    for (j = 0; j < 5; j++) {
      PCCf[i][j] = 0.0;
      PCCdfdx[i][j] = 0.0;
      PCCdfdy[i][j] = 0.0;
      PCHf[i][j] = 0.0;
      PCHdfdx[i][j] = 0.0;
      PCHdfdy[i][j] = 0.0;
    }
  }

  PCCf[0][2] = -0.00050;
  PCCf[0][3] = 0.0161253646;
  PCCf[1][1] = -0.010960;
  PCCf[1][2] = 0.00632624824;
  PCCf[2][0] = -0.0276030;
  PCCf[2][1] = 0.00317953083;

  PCHf[0][1] = 0.209336733;
  PCHf[0][2] = -0.0644496154;
  PCHf[0][3] = -0.303927546;
  PCHf[1][0] = 0.010;
  PCHf[1][1] = -0.125123401;
  PCHf[1][2] = -0.298905246;
  PCHf[2][0] = -0.122042146;
  PCHf[2][1] = -0.300529172;
  PCHf[3][0] = -0.307584705;

  for (i = 0; i < 5; i++) {
    for (j = 0; j < 5; j++) {
      for (k = 0; k < 10; k++) {
        piCCf[i][j][k] = 0.0;
        piCCdfdx[i][j][k] = 0.0;
        piCCdfdy[i][j][k] = 0.0;
        piCCdfdz[i][j][k] = 0.0;
        piCHf[i][j][k] = 0.0;
        piCHdfdx[i][j][k] = 0.0;
        piCHdfdy[i][j][k] = 0.0;
        piCHdfdz[i][j][k] = 0.0;
        piHHf[i][j][k] = 0.0;
        piHHdfdx[i][j][k] = 0.0;
        piHHdfdy[i][j][k] = 0.0;
        piHHdfdz[i][j][k] = 0.0;
        Tf[i][j][k] = 0.0;
        Tdfdx[i][j][k] = 0.0;
        Tdfdy[i][j][k] = 0.0;
        Tdfdz[i][j][k] = 0.0;
      }
    }
  }

  for (i = 3; i < 10; i++) piCCf[0][0][i] = 0.0049586079;
  piCCf[1][0][1] = 0.021693495;
  piCCf[0][1][1] = 0.021693495;
  for (i = 2; i < 10; i++) piCCf[1][0][i] = 0.0049586079;
  for (i = 2; i < 10; i++) piCCf[0][1][i] = 0.0049586079;
  piCCf[1][1][1] = 0.05250;
  piCCf[1][1][2] = -0.002088750;
  for (i = 3; i < 10; i++) piCCf[1][1][i] = -0.00804280;
  piCCf[2][0][1] = 0.024698831850;
  piCCf[0][2][1] = 0.024698831850;
  piCCf[2][0][2] = -0.00597133450;
  piCCf[0][2][2] = -0.00597133450;
  for (i = 3; i < 10; i++) piCCf[2][0][i] = 0.0049586079;
  for (i = 3; i < 10; i++) piCCf[0][2][i] = 0.0049586079;
  piCCf[2][1][1] = 0.00482478490;
  piCCf[1][2][1] = 0.00482478490;
  piCCf[2][1][2] = 0.0150;
  piCCf[1][2][2] = 0.0150;
  piCCf[2][1][3] = -0.010;
  piCCf[1][2][3] = -0.010;
  piCCf[2][1][4] = -0.01168893870;
  piCCf[1][2][4] = -0.01168893870;
  piCCf[2][1][5] = -0.013377877400;
  piCCf[1][2][5] = -0.013377877400;
  piCCf[2][1][6] = -0.015066816000;
  piCCf[1][2][6] = -0.015066816000;
  for (i = 7; i < 10; i++) piCCf[2][1][i] = -0.015066816000;
  for (i = 7; i < 10; i++) piCCf[1][2][i] = -0.015066816000;
  piCCf[2][2][1] = 0.0472247850;
  piCCf[2][2][2] = 0.0110;
  piCCf[2][2][3] = 0.0198529350;
  piCCf[2][2][4] = 0.01654411250;
  piCCf[2][2][5] = 0.013235290;
  piCCf[2][2][6] = 0.00992646749999 ;
  piCCf[2][2][7] = 0.006617644999;
  piCCf[2][2][8] = 0.00330882250;
  piCCf[3][0][1] = -0.05989946750;
  piCCf[0][3][1] = -0.05989946750;
  piCCf[3][0][2] = -0.05989946750;
  piCCf[0][3][2] = -0.05989946750;
  for (i = 3; i < 10; i++) piCCf[3][0][i] = 0.0049586079;
  for (i = 3; i < 10; i++) piCCf[0][3][i] = 0.0049586079;
  piCCf[3][1][2] = -0.0624183760;
  piCCf[1][3][2] = -0.0624183760;
  for (i = 3; i < 10; i++) piCCf[3][1][i] = -0.0624183760;
  for (i = 3; i < 10; i++) piCCf[1][3][i] = -0.0624183760;
  piCCf[3][2][1] = -0.02235469150;
  piCCf[2][3][1] = -0.02235469150;
  for (i = 2; i < 10; i++) piCCf[3][2][i] = -0.02235469150;
  for (i = 2; i < 10; i++) piCCf[2][3][i] = -0.02235469150;

  piCCdfdx[2][1][1] = -0.026250;
  piCCdfdx[2][1][5] = -0.0271880;
  piCCdfdx[2][1][6] = -0.0271880;
  for (i = 7; i < 10; i++) piCCdfdx[2][1][i] = -0.0271880;
  piCCdfdx[1][3][2] = 0.0187723882;
  for (i = 2; i < 10; i++) piCCdfdx[2][3][i] = 0.031209;

  piCCdfdy[1][2][1] = -0.026250;
  piCCdfdy[1][2][5] = -0.0271880;
  piCCdfdy[1][2][6] = -0.0271880;
  for (i = 7; i < 10; i++) piCCdfdy[1][2][i] = -0.0271880;
  piCCdfdy[3][1][2] = 0.0187723882;
  for (i = 2; i < 10; i++) piCCdfdy[3][2][i] = 0.031209;

  piCCdfdz[1][1][2] = -0.0302715;
  piCCdfdz[2][1][4] = -0.0100220;
  piCCdfdz[1][2][4] = -0.0100220;
  piCCdfdz[2][1][5] = -0.0100220;
  piCCdfdz[1][2][5] = -0.0100220;
  for (i = 4; i < 9; i++) piCCdfdz[2][2][i] = -0.0033090;

  //  make top end of piCC flat instead of zero
  i = 4;
  for (j = 0; j < 4; j++){
      for (k = 1; k < 11; k++){
          piCCf[i][j][k] = piCCf[i-1][j][k];
      }
  }
  for (i = 0; i < 4; i++){ // also enforces some symmetry
      for (j = i+1; j < 5; j++){
          for (k = 1; k < 11; k++){
              piCCf[i][j][k] = piCCf[j][i][k];
          }
      }
  }
  for (k = 1; k < 11; k++) piCCf[4][4][k] = piCCf[3][4][k];
  k = 10;
  for (i = 0; i < 5; i++){
      for (j = 0; j < 5; j++){
      piCCf[i][j][k] = piCCf[i][j][k-1];
      }
  }

  piCHf[1][1][1] = -0.050;
  piCHf[1][1][2] = -0.050;
  piCHf[1][1][3] = -0.30;
  for (i = 4; i < 10; i++) piCHf[1][1][i] = -0.050;
  for (i = 5; i < 10; i++) piCHf[2][0][i] = -0.004523893758064;
  for (i = 5; i < 10; i++) piCHf[0][2][i] = -0.004523893758064;
  piCHf[2][1][2] = -0.250;
  piCHf[1][2][2] = -0.250;
  piCHf[2][1][3] = -0.250;
  piCHf[1][2][3] = -0.250;
  piCHf[3][1][1] = -0.10;
  piCHf[1][3][1] = -0.10;
  piCHf[3][1][2] = -0.125;
  piCHf[1][3][2] = -0.125;
  piCHf[3][1][3] = -0.125;
  piCHf[1][3][3] = -0.125;
  for (i = 4; i < 10; i++) piCHf[3][1][i] = -0.10;
  for (i = 4; i < 10; i++) piCHf[1][3][i] = -0.10;

  // make top end of piCH flat instead of zero
 // also enforces some symmetry

  i = 4;
  for (j = 0; j < 4; j++){
      for (k = 1; k < 11; k++){
          piCHf[i][j][k] = piCHf[i-1][j][k];
      }
  }
  for (i = 0; i < 4; i++){
      for (j = i+1; j < 5; j++){
          for (k = 1; k < 11; k++){
              piCHf[i][j][k] = piCHf[j][i][k];
          }
      }
  }
  for (k = 1; k < 11; k++) piCHf[4][4][k] = piCHf[3][4][k];
  k = 10;
  for (i = 0; i < 5; i++){
      for (j = 0; j < 5; j++){
      piCHf[i][j][k] = piCHf[i][j][k-1];
      }
  }

  piHHf[1][1][1] = 0.124915958;

  Tf[2][2][1] = -0.035140;
  for (i = 2; i < 10; i++) Tf[2][2][i] = -0.0040480;
}
double Calc_Potential(int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *Bmat, double *qvec, double *Avec,double *rij_long, FILE* logfile,double *angle,Normal_struct *Normal,int *Normal_atom,int Interlayer,double *Fc,int *Neighb,int *NNeighb,int potential_flag,double &E_tot_,int *NBond_once,int *Bonded_once,int BN_flag,int Ncar_hyd,int N_init,int N_init_ILP)
{
  double E_tot;
  int atomi;
  double E_TORSION,E_KC,E_REBO,ETersoff;
  double E_TORSION_,E_KC_,E_REBO_;
  double RI, Interlayer_dis;
  double EDisp,ERep,Ecoulomb;
  
  dr_ij *dN = new dr_ij[1];
  
  E_TORSION = 0.0; 
  E_REBO    = 0.0;
  E_KC      = 0.0;
  E_tot     = 0.0;
  ETersoff  = 0.0;
  ERep      = 0.0;
  EDisp     = 0.0;
  Ecoulomb  = 0.0;
  
  if (potential_flag==0 || potential_flag==1)
    E_REBO = calc_dE_REBO_pot(angle,BondNum,Natoms,L,BondParams,Bonded,particle,Fc,NBond,E_REBO_,NBond_once,Bonded_once,Ncar_hyd);
  
  if (potential_flag==0 || potential_flag==1)
    E_TORSION = TORSION_pot(particle,BondNum,Natoms,L,Bonded,NBond,E_TORSION_,NBond_once,Bonded_once,Ncar_hyd);
  
  if (potential_flag==2 || potential_flag==3 || BN_flag==1)
    calc_E_Tersoff_Pot(angle,BondNum,Natoms,L,BondParams,Bonded,particle,Fc,NBond,ETersoff,NBond_once,Bonded_once,N_init);
  
  if ((potential_flag==0 || potential_flag==2) && Interlayer > 0){
    
    if(ILP_flag==0)
	E_KC = calc_F_Rep_pot(Natoms,particle,rij_long,BondParams,L,Normal,Normal_atom,Interlayer,Neighb,NNeighb,dN,E_KC_,Ncar_hyd);
      else if(ILP_flag==1)
	E_KC = calc_F_Rep_pot(Natoms,particle,rij_long,BondParams,L,Normal,Normal_atom,Interlayer,Neighb,NNeighb,dN,E_KC_,Natoms);
  }
    
  if (Interlayer > 0 && ILP_flag == 0 && (potential_flag==1 || potential_flag==3 || BN_flag==1)){
    calc_EvdW(Natoms,rij_long,particle,BondParams,Normal,Neighb,NNeighb,L,ERep,EDisp,RI,Interlayer_dis,N_init_ILP);
    calc_Ecoulomb(Natoms,particle,rij,Bmat,qvec,Avec,AtomParams,BondParams,rij_long,logfile,Neighb,NNeighb,potential_flag,L,Ecoulomb,N_init_ILP);
  }
  
  delete [] dN;
  
  //fprintf(logfile,"pot E_REBO=%.16f E_Torsion=%.16f E_Tersoff=%.16f E_KC=%.16f  EvdW=%.16f ERep=%.16f Ecoulomb=%.16f\n",E_REBO,E_TORSION,ETersoff,E_KC,EDisp,ERep,Ecoulomb);
  
  E_tot = E_REBO + E_TORSION + E_KC + ETersoff + ERep + EDisp + Ecoulomb;
  
  if(print_terms==1 )
    fprintf(logfile,"E_tot=%.16f E_AIREBO=%.16f E_Tersoff=%.16f E_KC=%.16f  E_ILP=%.16f E_Coulomb=%.16f\n",E_tot,E_REBO+E_TORSION,ETersoff,E_KC,EDisp+ERep,Ecoulomb);
  //cerr<<"E_REBO= "<<E_REBO<<" E_Torsion= "<<E_TORSION<<endl;
  //cerr<<"ERep= "<<ERep<<" EDisp= "<<EDisp<<" ETersoff= "<<ETersoff<<endl;
  //cerr<<"EKC= "<<E_KC<<endl;
  
  return(E_tot);
  
  
}

double E_Tersoff(int Natoms,int *Bonded,atom *particle,double *rij,double *Fc, BondParamsStruct *BondParams,double *angle,int *NBond,point L){
  
  double ETersoff=0.0;
  
  int atomi,atomj,atomk,i,k,typej,typei,typek,index,index2,NBond_j;
  double g,Theta_ijk,b_ij,r_ij,r_jk,Tersoff_lambda1,Tersoff_lambda2,Tersoff_lambda3,Tersoff_c,Tersoff_d,Tersoff_h,Tersoff_n,A,B,beta;
  double Xsi,Tersoff_R,Tersoff_D;
  double Fc_jk,Fc_ij;
  dr_ij rji;
  dr_ij rjk;
  double inside_acos;
  
  //#pragma omp for reduction(+:ETersoff)
  for(atomj=0 ; atomj < Natoms ; atomj++){ // Go over all central atoms j.
    //if(particle[atomi].type==0 && potential_flag==6)continue;
    typej  = particle[atomj].type;
    NBond_j=NBond[atomj];
    for(i=0 ; i < NBond_j ; i++){ // Go over all j neighbors i.
      atomi = Bonded[atomj*MaxNBond + i];
      typei = particle[atomi].type;

      //if(particle[atomj].layer != particle[atomi].layer || particle[atomi].layer)continue;
      
      index = typei*NAtomParams + typej;
      r_ij = R_ij(particle,atomi,atomj, L);
      //r_ij=R_ij_(particle[atomi].r[0],particle[atomj].r[0],particle[atomi].r[1],particle[atomj].r[1],particle[atomi].r[2],particle[atomj].r[2],L);
      //r_ij = rij[atomj*MaxNBond + i];
      Tersoff_n   = BondParams[index].Tersoff_n;
      A     = BondParams[index].Tersoff_A;
      B     = BondParams[index].Tersoff_B;
      beta  = BondParams[index].Tersoff_beta;
      Tersoff_lambda1 = BondParams[index].Tersoff_lambda1;
      Tersoff_lambda2 = BondParams[index].Tersoff_lambda2;
      Tersoff_R = BondParams[index].Tersoff_R;
      Tersoff_D = BondParams[index].Tersoff_D;
      //Fc_ij = Fc[atomj*MaxNBond + i];
      Fc_ij = Fc_(r_ij,Tersoff_R,Tersoff_D);
      Xsi=0.0;
      
      for(k=0 ; k < NBond_j ; k++){ // Go over all j neighbors k != i.
	if (k!=i){
	  atomk     = Bonded[atomj*MaxNBond + k];
	  typek     = particle[atomk].type;
	  index2    = typek*NAtomParams + typej;
	  Tersoff_c = BondParams[index2].Tersoff_c;
	  Tersoff_d = BondParams[index2].Tersoff_d;
	  Tersoff_h = BondParams[index2].Tersoff_h;
	  Tersoff_R = BondParams[index2].Tersoff_R;
	  Tersoff_D = BondParams[index2].Tersoff_D;
	  //r_jk      = rij[atomj*MaxNBond + k];
	  r_jk = R_ij(particle,atomj,atomk, L);
	  //r_jk=R_ij_(particle[atomj].r[0],particle[atomk].r[0],particle[atomj].r[1],particle[atomk].r[1],particle[atomj].r[2],particle[atomk].r[2],L);
	  
	  Theta_ijk = calc_val_angle(atomi,atomj,atomk,Natoms,NBond,Bonded,r_jk,r_ij,particle,L,rji,rjk,inside_acos);
	  g     = 1+sqr(Tersoff_c)/sqr(Tersoff_d)-sqr(Tersoff_c)/(sqr(Tersoff_d)+sqr(Tersoff_h-cos(Theta_ijk)));
	  //Fc_jk = Fc[atomj*MaxNBond + k];
	  Fc_jk = Fc_(r_jk,Tersoff_R,Tersoff_D);
	  Xsi  += Fc_jk*g;
	}
      }
      b_ij      = pow(1+pow(beta*Xsi,Tersoff_n),-1/(2*Tersoff_n));
      ETersoff  += Fc_ij*(A*exp(-Tersoff_lambda1*r_ij) - b_ij*B*exp(-Tersoff_lambda2*r_ij));
    }
  }
  
  ETersoff*=0.5;
  return(ETersoff);
  
}

double calc_dE_REBO_pot(double *angle,int *BondNum,int Natoms,point L,BondParamsStruct *BondParams,int *Bonded,atom *particle,double *Fc,int *NBond,double &E_REBO_,int *NBond_once,int *Bonded_once,int Ncar_hyd)
{
   double E_REBO;
  E_REBO=0.0;
#pragma omp parallel
  {
    int index,atomk,typek,typej,typei,index2;
    double g,Theta_ijk,b_ij,r_ij,r_jk,Tersoff_lambda1,Tersoff_lambda2,Tersoff_lambda3,Tersoff_c,Tersoff_d,Tersoff_h,Tersoff_n,A,B,beta;
    double Xsi,ETersoff,dTheta_ijk_drn,dg_drn,dXsi_drn,Fc_ij,Fc_jk,Tersoff_R,Tersoff_D;
    double dr_ij_drn,dr_jk_drn,dFc_ij_drn,dFc_jk_drn,Fa,Fr,dFa_drn,dFr_drn,db_ij_drn,n,dE_Tersoff_ij_drn;
    double Tersoff_sqr_c,Tersoff_sqr_d,sqr_Tersoff_h_Cos_Theta;
    double inside_acos,cos_ijk;
    int i,j,k,m,ii,itype,jtype,atomi,atomj;
    int itag ,jtag;
    double delx,dely,delz,evdwl,fpair,xtmp,ytmp,ztmp;
    double rsq,rij,wij;
    double Qij,Aij,alphaij,VR,pre,dVRdi,VA,term,bij,dVAdi,dVA;
    double dwij,del[3];
    
    int vflag_atom=0;
    evdwl = 0.0;
    
    // two-body interactions from REBO neighbor list, skip half of them
#pragma omp for reduction(+:E_REBO)
    for (atomi = 0; atomi < Ncar_hyd; atomi++) {
      itag=atomi;
      itype = particle[atomi].type;
      xtmp = particle[atomi].r[0];
      ytmp = particle[atomi].r[1];
      ztmp = particle[atomi].r[2];
      //REBO_neighs = REBO_firstneigh[i];
      
      for (j = 0; j < NBond_once[atomi]; j++) {
	atomj = Bonded_once[atomi*MaxNBond + j];

	jtype = particle[atomj].type;

	delx = R_PBC(particle,atomi,atomj,0,L.x);
	dely = R_PBC(particle,atomi,atomj,1,L.y);	
	delz = R_PBC(particle,atomi,atomj,2,L.z);
	
	rsq = delx*delx + dely*dely + delz*delz;
	rij = sqrt(rsq);
	
	wij = Sp(rij,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
	
	if (wij <= TOL) continue;
	Qij = Q[itype][jtype];
	
	Aij = A_ij[itype][jtype];
	
	alphaij = alpha[itype][jtype];
	
	VR = wij*(1.0+(Qij/rij)) * Aij*exp(-alphaij*rij);

	pre = wij*Aij * exp(-alphaij*rij);
		
	VA = dVA = 0.0;
	
	for (m = 0; m < 3; m++) {
	  term = -wij * BIJc[itype][jtype][m] * exp(-Beta[itype][jtype][m]*rij);
	  VA += term;
	}
	
	del[0] = delx;
	del[1] = dely;
	del[2] = delz;
	
	bij = bondorder_pot(atomi,atomj,del,rij,VA,vflag_atom,particle,Bonded,BondNum,NBond,L);
	//#pragma omp atomic
	E_REBO+=VR + bij*VA;
	//E_REBO+=VR;
     	
      }
    }
  }
  //cerr<<"E_REBO_03"<<E_REBO<<endl;
  return(E_REBO);
}


double bondorder_pot(int i, int j, double rij[3], double rijmag, double VA, int vflag_atom,atom *particle,int *Bonded,int *BondNum,int *NBond,point L){

   int atomi,atomj,k,n,l,atomk,atoml,atomn,atom1,atom2,atom3,atom4;
  int itype,jtype,ktype,ltype,ntype;
  double rik[3],rjl[3],rkn[3],rji[3],rki[3],rlj[3],rknmag,dNki,dwjl,bij;
  double NijC,NijH,NjiC,NjiH,wik,dwik,dwkn,wjl;
  double rikmag,rjlmag,cosjik,cosijl,g,tmp2,tmp3;
  double Etmp,pij,tmp,wij,dwij,NconjtmpI,NconjtmpJ,Nki,Nlj,dS;
  double lamdajik,lamdaijl,dgdc,dgdN,pji,Nijconj,piRC;
  double dcosjikdri[3],dcosijldri[3],dcosjikdrk[3];
  double dN2[2],dN3[3];
  double dcosjikdrj[3],dcosijldrj[3],dcosijldrl[3];
  double Tij;
  double r32[3],r32mag,cos321,r43[3],r13[3];
  double dNlj;
  double om1234,rln[3];
  double rlnmag,dwln,r23[3],r23mag,r21[3],r21mag;
  double w21,dw21,r34[3],r34mag,cos234,w34,dw34;
  double cross321[3],cross234[3],prefactor,SpN;
  double fcijpc,fcikpc,fcjlpc,fcjkpc,fcilpc;
  double dt2dik[3],dt2djl[3],dt2dij[3],aa,aaa1,aaa2,at2,cw,cwnum,cwnom;
  double sin321,sin234,rr,rijrik,rijrjl,rjk2,rik2,ril2,rjl2;
  double dctik,dctjk,dctjl,dctij,dctji,dctil,rik2i,rjl2i,sink2i,sinl2i;
  double rjk[3],ril[3],dt1dik,dt1djk,dt1djl,dt1dil,dt1dij;
  double F23[3],F12[3],F34[3],F31[3],F24[3],fi[3],fj[3],fk[3],fl[3];
  double f1[3],f2[3],f3[3],f4[4];
  double dcut321,PijS,PjiS;
  double rij2,tspjik,dtsjik,tspijl,dtsijl,costmp;
  int *REBO_neighs,*REBO_neighs_i,*REBO_neighs_j,*REBO_neighs_k,*REBO_neighs_l;

  atomi = i;
  atomj = j;
  itype = particle[atomi].type;
  jtype = particle[atomj].type;

  wij = Sp(rijmag,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
  NijC = nC[i]-(wij*kronecker(jtype,0));
  NijH = nH[i]-(wij*kronecker(jtype,1));
  NjiC = nC[j]-(wij*kronecker(itype,0));
  NjiH = nH[j]-(wij*kronecker(itype,1));

  bij = 0.0;
  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  dgdc = 0.0;
  dgdN = 0.0;
  NconjtmpI = 0.0;
  NconjtmpJ = 0.0;
  Etmp = 0.0;
  
  for(k=0 ; k < NBond[atomi] ; k++){ // Go over all j neighbors k != i.
      //atomk = REBO_neighs[k];
      atomk = Bonded[atomi*MaxNBond + k];
    if (atomk != atomj) {
      ktype = particle[atomk].type;
      //ktype = particle[atomk].type;
      rik[0] = R_PBC(particle,atomi,atomk,0,L.x);
      rik[1] = R_PBC(particle,atomi,atomk,1,L.y);
      rik[2] = R_PBC(particle,atomi,atomk,2,L.z);
      
      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      lamdajik = 4.0*kronecker(itype,1) *
        ((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag));
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dS);
      Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] -
        (wik*kronecker(itype,1));
      cosjik = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2])) /
        (rijmag*rikmag);
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosjik,(NijC+NijH),itype,&dgdc,&dgdN);
      Etmp = Etmp+(wik*g*exp(lamdajik));
      tmp3 = tmp3+(wik*dgdN*exp(lamdajik));
      NconjtmpI = NconjtmpI+(kronecker(ktype,0)*wik*Sp(Nki,Nmin,Nmax,dS)); //parameter Nmin Nmax
    }
  }

  PijS = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  PijS = PijSpline(NijC,NijH,itype,jtype,dN2);
  pij = pow(1.0+Etmp+PijS,-0.5);
  tmp = -0.5*cube(pij);
  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  Etmp = 0.0;
  
  for (l = 0; l < NBond[atomj]; l++) {
    atoml = Bonded[atomj*MaxNBond + l];
    if (atoml != atomi) {
      ltype = particle[atoml].type;
      rjl[0] = R_PBC(particle,atomj,atoml,0,L.x);
      rjl[1] = R_PBC(particle,atomj,atoml,1,L.y);
      rjl[2] = R_PBC(particle,atomj,atoml,2,L.z);
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      lamdaijl = 4.0*kronecker(jtype,1) *
        ((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dS);
      Nlj = nC[atoml]-(wjl*kronecker(jtype,0)) +
        nH[atoml]-(wjl*kronecker(jtype,1));
      cosijl = -1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2])) /
        (rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosijl,NjiC+NjiH,jtype,&dgdc,&dgdN);
      Etmp = Etmp+(wjl*g*exp(lamdaijl));
      tmp3 = tmp3+(wjl*dgdN*exp(lamdaijl));
      NconjtmpJ = NconjtmpJ+(kronecker(ltype,0)*wjl*Sp(Nlj,Nmin,Nmax,dS));
    }
  }

  PjiS = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  
  PjiS = PijSpline(NjiC,NjiH,jtype,itype,dN2);
  
  pji = pow(1.0+Etmp+PjiS,-0.5);
  tmp = -0.5*cube(pji);
  // evaluate Nij conj
  
  Nijconj = 1.0+(NconjtmpI*NconjtmpI)+(NconjtmpJ*NconjtmpJ);
  piRC = piRCSpline(NijC+NijH,NjiC+NjiH,Nijconj,itype,jtype,dN3);
  
  // piRC forces
  Tij = 0.0;
  dN3[0] = 0.0;
  dN3[1] = 0.0;
  dN3[2] = 0.0;
  if (itype == 0 && jtype == 0)
    Tij=TijSpline((NijC+NijH),(NjiC+NjiH),Nijconj,dN3);
  Etmp = 0.0;

  if (fabs(Tij) > TOL) {
    atom2 = atomi;
    atom3 = atomj;
    r32[0] = R_PBC(particle,atom3,atom2,0,L.x);
    r32[1] = R_PBC(particle,atom3,atom2,1,L.y);
    r32[2] = R_PBC(particle,atom3,atom2,2,L.z);
    r32mag = sqrt((r32[0]*r32[0])+(r32[1]*r32[1])+(r32[2]*r32[2]));
    r23[0] = -r32[0];
    r23[1] = -r32[1];
    r23[2] = -r32[2];
    r23mag = r32mag;
    //REBO_neighs_i = REBO_firstneigh[i];
    for (k = 0; k < NBond[atomi]; k++) {
      atomk = Bonded[atomi*MaxNBond + k];
      atom1 = atomk;
      ktype = particle[atomk].type;
      if (atomk != atomj) {
	r21[0] = R_PBC(particle,atom2,atom1,0,L.x);
	r21[1] = R_PBC(particle,atom2,atom1,1,L.y);
	r21[2] = R_PBC(particle,atom2,atom1,2,L.z);
        r21mag = sqrt(r21[0]*r21[0] + r21[1]*r21[1] + r21[2]*r21[2]);
        cos321 = -1.0*((r21[0]*r32[0])+(r21[1]*r32[1])+(r21[2]*r32[2])) /
          (r21mag*r32mag);
        cos321 = MIN(cos321,1.0);
        cos321 = MAX(cos321,-1.0);
        Sp2(cos321,thmin,thmax,dcut321);
        sin321 = sqrt(1.0 - cos321*cos321);
        sink2i = 1.0/(sin321*sin321);
        rik2i = 1.0/(r21mag*r21mag);
        if (sin321 != 0.0) {
          rr = (r23mag*r23mag)-(r21mag*r21mag);
          //rjk[0] = r21[0]-r23[0];
          //rjk[1] = r21[1]-r23[1];
          //rjk[2] = r21[2]-r23[2];

	  rjk[0] = R_PBC(particle,atomj,atomk,0,L.x);
	  rjk[1] = R_PBC(particle,atomj,atomk,1,L.y);
	  rjk[2] = R_PBC(particle,atomj,atomk,2,L.z);
	      
          rjk2 = (rjk[0]*rjk[0])+(rjk[1]*rjk[1])+(rjk[2]*rjk[2]);
          rijrik = 2.0*r23mag*r21mag;
          rik2 = r21mag*r21mag;
          dctik = (-rr+rjk2)/(rijrik*rik2);
          dctij = (rr+rjk2)/(rijrik*r23mag*r23mag);
          dctjk = -2.0/rijrik;
          w21 = Sp(r21mag,rcmin[itype][ktype],rcmaxp[itype][ktype],dw21);
          rijmag = r32mag;
          rikmag = r21mag;
          rij2 = r32mag*r32mag;
          rik2 = r21mag*r21mag;
          costmp = 0.5*(rij2+rik2-rjk2)/rijmag/rikmag;
          tspjik = Sp2(costmp,thmin,thmax,dtsjik);
          dtsjik = -dtsjik;

          //REBO_neighs_j = REBO_firstneigh[j];
          for (l = 0; l < NBond[atomj]; l++) {
            //atoml = REBO_neighs_j[l];
	    atoml = Bonded[atomj*MaxNBond + l];
            atom4 = atoml;
            ltype = particle[atoml].type;
            if (!(atoml == atomi || atoml == atomk)) {
	      r34[0] = R_PBC(particle,atom3,atom4,0,L.x);
	      r34[1] = R_PBC(particle,atom3,atom4,1,L.y);
	      r34[2] = R_PBC(particle,atom3,atom4,2,L.z);
              r34mag = sqrt((r34[0]*r34[0])+(r34[1]*r34[1])+(r34[2]*r34[2]));
              cos234 = (r32[0]*r34[0] + r32[1]*r34[1] + r32[2]*r34[2]) /
                (r32mag*r34mag);
              cos234 = MIN(cos234,1.0);
              cos234 = MAX(cos234,-1.0);
              sin234 = sqrt(1.0 - cos234*cos234);
              sinl2i = 1.0/(sin234*sin234);
              rjl2i = 1.0/(r34mag*r34mag);

              if (sin234 != 0.0) {
                w34 = Sp(r34mag,rcmin[jtype][ltype],rcmaxp[jtype][ltype],dw34);
                rr = (r23mag*r23mag)-(r34mag*r34mag);
                //ril[0] = r23[0]+r34[0];
                //ril[1] = r23[1]+r34[1];
                //ril[2] = r23[2]+r34[2];
		ril[0] = R_PBC(particle,atomi,atoml,0,L.x);
		ril[1] = R_PBC(particle,atomi,atoml,1,L.y);
		ril[2] = R_PBC(particle,atomi,atoml,2,L.z);
                ril2 = (ril[0]*ril[0])+(ril[1]*ril[1])+(ril[2]*ril[2]);
                rijrjl = 2.0*r23mag*r34mag;
                rjl2 = r34mag*r34mag;
                dctjl = (-rr+ril2)/(rijrjl*rjl2);
                dctji = (rr+ril2)/(rijrjl*r23mag*r23mag);
                dctil = -2.0/rijrjl;
                rjlmag = r34mag;
                rjl2 = r34mag*r34mag;
                costmp = 0.5*(rij2+rjl2-ril2)/rijmag/rjlmag;
                tspijl = Sp2(costmp,thmin,thmax,dtsijl);
                dtsijl = -dtsijl;
                prefactor = VA*Tij;

                cross321[0] = (r32[1]*r21[2])-(r32[2]*r21[1]);
                cross321[1] = (r32[2]*r21[0])-(r32[0]*r21[2]);
                cross321[2] = (r32[0]*r21[1])-(r32[1]*r21[0]);
                cross234[0] = (r23[1]*r34[2])-(r23[2]*r34[1]);
                cross234[1] = (r23[2]*r34[0])-(r23[0]*r34[2]);
                cross234[2] = (r23[0]*r34[1])-(r23[1]*r34[0]);

                cwnum = (cross321[0]*cross234[0]) +
                  (cross321[1]*cross234[1]) + (cross321[2]*cross234[2]);
                cwnom = r21mag*r34mag*r23mag*r23mag*sin321*sin234;
                om1234 = cwnum/cwnom;
                cw = om1234;
                Etmp += ((1.0-square(om1234))*w21*w34) *
                  (1.0-tspjik)*(1.0-tspijl);
		
	      }
            }
          }
	  
        }
      }
    }
  }
  bij = (0.5*(pij+pji))+piRC+(Tij*Etmp);
  
  return bij;
}
double TORSION_pot(atom *particle,int *BondNum,int Natoms,point L,int *Bonded,int *NBond,double &E_TORSION_,int *NBond_once,int *Bonded_once,int Ncar_hyd)
{
  double E_TORSION;
  E_TORSION=0;
#pragma omp parallel
  {
    int i,j,k,l,ii,inum;
    //tagint itag,jtag;
    int itag,jtag;
    double evdwl,fpair,xtmp,ytmp,ztmp;
    double cos321;
    double w21,dw21,cos234,w34,dw34;
    double cross321[3],cross321mag,cross234[3],cross234mag;
    double w23,dw23,cw2,ekijl,Ec;
    double cw,cwnum,cwnom;
    double rij,rij2,rik,rjl,tspjik,dtsjik,tspijl,dtsijl,costmp,fcpc;
    double sin321,sin234,rjk2,rik2,ril2,rjl2;
    double rjk,ril;
    double Vtors;
    double dndij[3],tmpvec[3],dndik[3],dndjl[3];
    double dcidij,dcidik,dcidjk,dcjdji,dcjdjl,dcjdil;
    double dsidij,dsidik,dsidjk,dsjdji,dsjdjl,dsjdil;
    double dxidij,dxidik,dxidjk,dxjdji,dxjdjl,dxjdil;
    double ddndij,ddndik,ddndjk,ddndjl,ddndil,dcwddn,dcwdn,dvpdcw,Ftmp[3];
    double del32[3],rsq,r32,del23[3],del21[3],r21;
    double deljk[3],del34[3],delil[3],delkl[3],r23,r34;
    double fi[3],fj[3],fk[3],fl[3];
    int itype,jtype,ktype,ltype,kk,ll,jj;
    int *ilist,*REBO_neighs_i,*REBO_neighs_j;
    
#pragma omp for reduction(+:E_TORSION)
    for (ii = 0; ii < Ncar_hyd; ii++) {
      i=ii;
      
      itag = ii;
      itype = particle[i].type;
      if (itype != 0) continue;
      xtmp = particle[i].r[0];
      ytmp = particle[i].r[1];
      ztmp = particle[i].r[2];
      
      for (jj = 0; jj < NBond_once[ii]; jj++) {
	
	j = Bonded_once[i*MaxNBond + jj];
	jtag = j;
	
	jtype = particle[j].type;
	if (jtype != 0) continue;
	
	del32[0] = R_PBC(particle,j,i,0,L.x);
	del32[1]  = R_PBC(particle,j,i,1,L.y);
	del32[2]  = R_PBC(particle,j,i,2,L.z);
	rsq = del32[0]*del32[0] + del32[1]*del32[1] + del32[2]*del32[2];
	r32 = sqrt(rsq);
	del23[0] = -del32[0];
	del23[1] = -del32[1];
	del23[2] = -del32[2];
	r23 = r32;
	w23 = Sp(r23,rcmin[itype][jtype],rcmax[itype][jtype],dw23);
	//cerr<<"hey3"<<endl;
	for (kk = 0; kk < NBond[i]; kk++) {
	  //k = REBO_neighs_i[kk];
	  k = Bonded[i*MaxNBond + kk];
	  ktype = particle[k].type;
	  if (k == j) continue;
	  del21[0] = R_PBC(particle,i,k,0,L.x);
	  del21[1] = R_PBC(particle,i,k,1,L.y);
	  del21[2] = R_PBC(particle,i,k,2,L.z);
	  rsq = del21[0]*del21[0] + del21[1]*del21[1] + del21[2]*del21[2];
	  r21 = sqrt(rsq);
	  cos321 = - ((del21[0]*del32[0]) + (del21[1]*del32[1]) +
		      (del21[2]*del32[2])) / (r21*r32);
	  cos321 = MIN(cos321,1.0);
	  cos321 = MAX(cos321,-1.0);
	  sin321 = sqrt(1.0 - cos321*cos321);
	  if (sin321 < TOL) continue;
	  
	  //deljk[0] = del21[0]-del23[0];
	  //deljk[1] = del21[1]-del23[1];
	  //deljk[2] = del21[2]-del23[2];
	  
	  deljk[0] = R_PBC(particle,j,k,0,L.x);
	  deljk[1] = R_PBC(particle,j,k,1,L.y);
	  deljk[2] = R_PBC(particle,j,k,2,L.z);
	  
	  rjk2 = deljk[0]*deljk[0] + deljk[1]*deljk[1] + deljk[2]*deljk[2];
	  rjk=sqrt(rjk2);
	  rik2 = r21*r21;
	  w21 = Sp(r21,rcmin[itype][ktype],rcmax[itype][ktype],dw21);
	  
	  rij = r32;
	  rik = r21;
	  rij2 = r32*r32;
	  rik2 = r21*r21;
	  costmp = 0.5*(rij2+rik2-rjk2)/rij/rik;
	  tspjik = Sp2(costmp,thmin,thmax,dtsjik);
	  dtsjik = -dtsjik;
	  //REBO_neighs_j = REBO_firstneigh[j];
	  for (ll = 0; ll <NBond[j]; ll++) {
	    //l = REBO_neighs_j[ll];
	    l = Bonded[j*MaxNBond + ll];
	    ltype = particle[l].type;
	    if (l == i || l == k) continue;
	    del34[0] = R_PBC(particle,j,l,0,L.x);
	    del34[1] = R_PBC(particle,j,l,1,L.y);
	    del34[2] = R_PBC(particle,j,l,2,L.z);
	    rsq = del34[0]*del34[0] + del34[1]*del34[1] + del34[2]*del34[2];
	    r34 = sqrt(rsq);
	    cos234 = (del32[0]*del34[0] + del32[1]*del34[1] +
		      del32[2]*del34[2]) / (r32*r34);
	    cos234 = MIN(cos234,1.0);
	    cos234 = MAX(cos234,-1.0);
	    sin234 = sqrt(1.0 - cos234*cos234);
	    if (sin234 < TOL) continue;
	    w34 = Sp(r34,rcmin[jtype][ltype],rcmax[jtype][ltype],dw34);
	    
	    delil[0] = R_PBC(particle,i,l,0,L.x);
	    delil[1] = R_PBC(particle,i,l,1,L.y);
	    delil[2] = R_PBC(particle,i,l,2,L.z);
	    
	    ril2 = delil[0]*delil[0] + delil[1]*delil[1] + delil[2]*delil[2];
	    ril=sqrt(ril2);
	    rjl2 = r34*r34;
	    
	    rjl = r34;
	    rjl2 = r34*r34;
	    costmp = 0.5*(rij2+rjl2-ril2)/rij/rjl;
	    tspijl = Sp2(costmp,thmin,thmax,dtsijl);
	    dtsijl = -dtsijl; //need minus sign
	    cross321[0] = (del32[1]*del21[2])-(del32[2]*del21[1]);
	    cross321[1] = (del32[2]*del21[0])-(del32[0]*del21[2]);
	    cross321[2] = (del32[0]*del21[1])-(del32[1]*del21[0]);
	    cross321mag = sqrt(cross321[0]*cross321[0]+
			       cross321[1]*cross321[1]+
			       cross321[2]*cross321[2]);
	    cross234[0] = (del23[1]*del34[2])-(del23[2]*del34[1]);
	    cross234[1] = (del23[2]*del34[0])-(del23[0]*del34[2]);
	    cross234[2] = (del23[0]*del34[1])-(del23[1]*del34[0]);
	    cross234mag = sqrt(cross234[0]*cross234[0]+
			       cross234[1]*cross234[1]+
			       cross234[2]*cross234[2]);
	    cwnum = (cross321[0]*cross234[0]) +
	      (cross321[1]*cross234[1])+(cross321[2]*cross234[2]);
	    cwnom = r21*r34*r32*r32*sin321*sin234;
	    cw = cwnum/cwnom;
	    
	    cw2 = (.5*(1.0-cw));
	    ekijl = epsilonT[ktype][ltype];
	    Ec = 256.0*ekijl/405.0;
	    Vtors = (Ec*(powint(cw2,5)))-(ekijl/10.0);
	    
	    //#pragma omp atomic
	    E_TORSION += Vtors*w21*w23*w34*(1.0-tspjik)*(1.0-tspijl);
	    
	  }
	}
      }
    }
  }
  return(E_TORSION);
}
double calc_F_Rep_pot(int Natoms, atom *particle, double *rij_long, BondParamsStruct *BondParams, point L,Normal_struct *Normal,int *Normal_atom,int Interlayer,int *Neighb,int *NNeighb,dr_ij *dN,double &ERep_,int Ncar_hyd)
{
  
  double ERep=0.0;  
#pragma omp parallel
  {
    int atomk, atomi, typek, typei,atomj,typej, index,j;
    double alpha, rvdW, Epsilon, term1, dEvdW_drij;
    double dTap, r_ij, Tap,Epsilon_C_vdW_exp,term3,term4,term5,term7;
    double dP_ij_drk,dP_ji_drk,R,Pij,Pji,gamma,dr_ij_drk,exp_term1;
    double exp_Pij_gamma, exp_Pji_gamma,Pij_gamma,Pji_gamma,reff,r_ik,r_jk;
    double rx_ij,ry_ij,rz_ij,rx_ji,ry_ji,rz_ji,inside_acos_ij,inside_acos_ji,Pij_term1,Pji_term1,Pij_term,Pji_term,dr_ji_drk;
    double angle_ij, sqrt_inside_ij, sin_angle_ij,cos_angle_ij,angle_ji,sqrt_inside_ji,sin_angle_ji,cos_angle_ji,dinside_acos,dangle;
    double C0_kc,C2_kc,C4_kc,A_kc,fji,fij,term6,dfji,dfij,A_term;
    double Nori_L,Norj_L,Fasolina;
    int n,dir;
    dP_ij_drk = 0.0;
    dP_ji_drk = 0.0;
    
    C0_kc      = 15.71/1000;
    C2_kc      = 12.29/1000;
    C4_kc      = 4.933/1000;
    alpha   = 3.629;
    Epsilon = 3.03/1000;
    R       = 3.34;
    gamma   = 0.578;
    A_kc    = 10.238/1000;
    
#pragma omp for schedule(dynamic,15) reduction(+:ERep)
    for(atomi=0 ; atomi < Ncar_hyd ; atomi++){
      typei = particle[atomi].type;
      if(typei==1)continue;
      for (j=0 ; j < NNeighb[atomi]; j++){
	atomk = atomj = Neighb[atomi*MaxNeighb + j] ;	  
	typek = typej = particle[atomj].type;
	if(typej==1)continue;
	index = typej*NAtomParams + typei;

	if(typei==0&&typej==0)Fasolina=1.0;
	else if((typei==0 && typej==7) || typei==7 && typej==0)Fasolina=2.0;
	else if((typei==0 && typej==5) || typei==5 && typej==0)Fasolina=0.6;
	else{
	  cerr<<"wrong atom type in F_KC atomi= "<<atomi<<"atomj= "<<atomj<< " exiting"<<endl;
	  exit(0);
	}
	

	rx_ij = R_PBC(particle,atomi,atomj,0,L.x);
	ry_ij = R_PBC(particle,atomi,atomj,1,L.y);	
	rz_ij = R_PBC(particle,atomi,atomj,2,L.z);
	
	rx_ji = -rx_ij;
	ry_ji = -ry_ij;
	rz_ji = -rz_ij;

	r_ij = sqrt(rx_ij*rx_ij + ry_ij*ry_ij + rz_ij*rz_ij);
	
	//r_ij=R_ij_(particle[atomi].r[0],particle[atomj].r[0],particle[atomi].r[1],particle[atomj].r[1],particle[atomi].r[2],particle[atomj].r[2],L);
	Tap   = calc_Tap(r_ij);
	
	reff    = BondParams[index].r_eff;	    
	
	Pij_term1=0.0;
	Pji_term1=0.0;
	dN[0].r[0]=0;
	dN[0].r[1]=0;
	dN[0].r[2]=0;
		
	Nori_L=rx_ij*Normal[atomi].x + ry_ij*Normal[atomi].y + rz_ij*Normal[atomi].z;
	Norj_L=rx_ji*Normal[atomj].x + ry_ji*Normal[atomj].y + rz_ji*Normal[atomj].z;
	
	inside_acos_ij = (Nori_L) / (r_ij);
	inside_acos_ji = (Norj_L) / (r_ij);
	
	if((r_ij*r_ij-sqr(Nori_L)) <= 0)Pij=0.0;
	else Pij=sqrt(r_ij*r_ij - sqr(Nori_L));
	
	if((r_ij*r_ij-sqr(Norj_L)) <= 0)Pji=0.0;
	else Pji=sqrt(r_ij*r_ij - sqr(Norj_L));
	
	Pij_gamma = Pij/gamma; //new
	Pji_gamma = Pji/gamma; //new
	exp_Pij_gamma = exp(-pow((Pij_gamma),2.0)); //new
	exp_Pji_gamma = exp(-pow((Pji_gamma),2.0)); //new
	exp_term1 = exp(-alpha*(r_ij-R));
	
	fij=0.0;
	fji=0.0;
	
	fij = C0_kc+ C2_kc*pow((Pij_gamma),2.0)+ C4_kc*pow((Pij_gamma),4.0);
	fji = C0_kc+ C2_kc*pow((Pji_gamma),2.0)+ C4_kc*pow((Pji_gamma),4.0);
	
	fij*=exp_Pij_gamma;
	fji*=exp_Pji_gamma;
	
	A_term=A_kc*pow((r_ij/R),-6);
	
	Epsilon_C_vdW_exp = Epsilon + fij + fji;
	//#pragma omp atomic
	
	ERep += Fasolina*Tap*(exp_term1*(Epsilon_C_vdW_exp) - A_term);
			
      }
    }
  }
  return(ERep);
}

void calc_EvdW(int Natoms, double *rij_long, atom *particle, BondParamsStruct *BondParams,Normal_struct *Normal,int *Neighb,int *NNeighb,point L,double &ERep,double &EDisp,double &RI,double &Interlayer_dis,int N_init_ILP)
{
  
  double ERep_,EDisp_;
  ERep_ = EDisp_ = 0.0;
  
#pragma omp parallel
  {
    double EvdW;
    double Pavg;
    int atomi, atomj, typei, typej, index,j;
    double lambdaW, alpha, rvdW, Epsilon, f13, term,Inner_shield;
    double ecore, acore, rcore, Tap, r_ij,C6,C12;
    double gamma,Pij,Pji,angle,rx,ry,rz;
    double inside_acos,r_eff;
    double R;
    double rx_ij,ry_ij,rz_ij,rx_ji,ry_ji,rz_ji,inside_acos_ij,inside_acos_ji,Pij_term1,Pji_term1,Pij_term,Pji_term,dr_ji_drk;
    double angle_ij, sqrt_inside_ij, sin_angle_ij,angle_ji,sqrt_inside_ji,sin_angle_ji,cos_angle_ji,dinside_acos,dangle;
    double C_vdW,term2,C6_r_ij,Pij_gamma,Pji_gamma,exp_Pij_gamma,exp_Pji_gamma;
    double Nori_L,Norj_L,exp_term1,Epsilon_C_vdW_exp;
    //EvdW = ERep_ = EDisp_ = 0.0;

       
#pragma omp for reduction(+:ERep_,EDisp_)
    for(atomi=N_init_ILP ; atomi < Natoms ; atomi++){
      typei = particle[atomi].type;
      
      for(j=0 ; j < NNeighb[atomi] ; j++){
	atomj = Neighb[atomi*MaxNeighb + j] ;
	
	typej   = particle[atomj].type;
	index   = typei*NAtomParams + typej;
	//r_ij    = R_ij(particle,atomi,atomj, L);
	//r_ij=R_ij_(particle[atomi].r[0],particle[atomj].r[0],particle[atomi].r[1],particle[atomj].r[1],particle[atomi].r[2],particle[atomj].r[2],L);
	
	//r_ij    = rij_long[atomi*MaxNeighb + j];
	rx_ij = R_PBC(particle,atomi,atomj,0,L.x);
	ry_ij = R_PBC(particle,atomi,atomj,1,L.y);
	rz_ij = R_PBC(particle,atomi,atomj,2,L.z);
	
	rx_ji = -rx_ij;
	ry_ji = -ry_ij;
	rz_ji = -rz_ij;
	
	r_ij =sqrt(rx_ij*rx_ij + ry_ij*ry_ij + rz_ij*rz_ij);
	
	Tap     = calc_Tap(r_ij);
	r_eff   = BondParams[index].r_eff;
	alpha   = BondParams[index].alpha_long;
	Epsilon = BondParams[index].Epsilon_long;
	R       = BondParams[index].rvdW_long;
	C6      = BondParams[index].C6;
	gamma   = BondParams[index].Trans;
	C_vdW   = BondParams[index].C_vdW;
	//Pij   = Calc_Pij(atomi,atomj,particle,Normal,r_ij,L);
	//Pji   = Calc_Pij(atomj,atomi,particle,Normal,r_ij,L);
	
	
	
	Nori_L=rx_ij*Normal[atomi].x + ry_ij*Normal[atomi].y + rz_ij*Normal[atomi].z;
	Norj_L=rx_ji*Normal[atomj].x + ry_ji*Normal[atomj].y + rz_ji*Normal[atomj].z;
	
	inside_acos_ij = (Nori_L) / (r_ij);
	inside_acos_ji = (Norj_L) / (r_ij);
	
	if((r_ij*r_ij-sqr(Nori_L)) <= 0)Pij=0.0;
	else Pij=sqrt(r_ij*r_ij - sqr(Nori_L));
	
	if((r_ij*r_ij-sqr(Norj_L)) <= 0)Pji=0.0;
	else Pji=sqrt(r_ij*r_ij - sqr(Norj_L));
	
	Pij_gamma = Pij/gamma; //new
	Pji_gamma = Pji/gamma; //new
	exp_Pij_gamma = exp(-pow((Pij_gamma),2.0)); //new
	exp_Pji_gamma = exp(-pow((Pji_gamma),2.0)); //new
	exp_term1 = exp(alpha*(1.0-r_ij/R));
	Epsilon_C_vdW_exp = Epsilon + C_vdW*(exp_Pij_gamma + exp_Pji_gamma);
	
	C6_r_ij = C6/(pow(r_ij,6.0)+0.1);
	term2   = exp(-d_TS*(r_ij/(Sr_TS*r_eff) - 1.0));
	
	//printf("exp_term1=%.16f \n",exp_term1);
	//printf("Epsilon_C_vdW_exp=%.16f \n",Epsilon_C_vdW_exp);
	
	//ERep_ += Tap;
	//ERep_ += (exp_term1*(Epsilon_C_vdW_exp) );
	//if(particle[atomi].layer==(particle[atomj].layer+1) || particle[atomi].layer==(particle[atomj].layer-1) || particle[atomi].layer==(particle[atomj].layer+2) || particle[atomi].layer==(particle[atomj].layer-2) || particle[atomi].layer==(particle[atomj].layer+3) || particle[atomi].layer==(particle[atomj].layer-3) ){

	  if(r_ij < non_bond_cut_Rep)ERep_ += Tap*(exp_term1*(Epsilon_C_vdW_exp) );
	  EDisp_ += Tap*(-(1.0/(1.0 + term2 ))*C6_r_ij);

	//ERep_ +=Tap*(exp(alpha*(1.0-r_ij/R))*(Epsilon + C_vdW*(exp(-pow((Pij/gamma),2.0)) + exp(-pow((Pji/gamma),2.0))) ));
	//EDisp_+=Tap*(-(1.0/(1.0 + exp(-d_TS*(r_ij/(Sr_TS*r_eff)- 1.0))))*C6/(pow(r_ij,6.0)+0.1));
      }
    }
  }
  ERep = ERep_;
  EDisp= EDisp_;
}

double Calc_Pij(int atomi,int atomj,atom *particle,Normal_struct *Normal,double r_ij,point L)
{
  double rx,ry,rz,inside_acos;
  double angle;

  rx = R_PBC(particle,atomi,atomj,0,L.x);
  ry = R_PBC(particle,atomi,atomj,1,L.y);
  rz = R_PBC(particle,atomi,atomj,2,L.z);
  
  inside_acos = (rx*Normal[atomi].x + ry*Normal[atomi].y + rz*Normal[atomi].z) / (r_ij);

  //if(atomj==640)cerr<<"rx= "<<rx<<"ry= "<<ry<<"rz= "<<rz<<"inside_acos= "<<inside_acos<<"Normal[atomi].x= "<<Normal[atomi].x<<"Normal[atomi].y= "<<Normal[atomi].y<<"Normal[atomi].z= "<<Normal[atomi].z<<endl;;
  
  if (fabs(inside_acos) < 1.0){
    angle = acos(inside_acos);
  }
  else if (inside_acos >= 1.0) angle = 0.0;
  else angle = PIE;

  return(r_ij*sin(angle));
  //return(Normal[atomi].x+Normal[atomi].y+Normal[atomi].z);
}
double calc_Ecoulomb(int Natoms, atom *particle, double *rij, double *Bmat, double *qvec, double *Avec, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, double *rij_long, FILE* logfile,int *Neighb,int *NNeighb,int potential_flag,point L,double &Ecoulomb,int N_init_ILP)
{
  //int atomi, atomj;
  //double Ecoulomb, Tap, r_ij;
  //double rijShield,shield,gamma;
  int atomi;
  // double Ecoulomb;
  
  Ecoulomb=0.0;
  
  //if (potential_flag==1)
  //Calc_Charge(particle,Natoms,Bmat,qvec,Avec,AtomParams,BondParams,rij,rij_long,logfile,L);
  
  //****
  /*
    for (atomi=0;atomi < Natoms; atomi++){
    //if (particle[i].type==5)qvec[i]=0.467;
    //else if (particle[i].type==7)qvec[i]=-0.467;
    //else qvec[i]=0.0;
    //cerr<<"atomi= "<<atomi<<" charge ="<<particle[atomi].charge<<endl;
    //printf("atomi=%i charge=%.16f\n",i,qvec[i]);
    particle[atomi].charge = particle[atomi].charge;
    }
  */ 
  //printf("\n\n");
  
  //exit(0);  
  //***
  
  //#pragma omp parallel 
  //{
  int atomj,j,i;
  double Tap, r_ij;
  double rijShield,shield,gamma;
  //#pragma omp for
  for (atomi = 0;atomi < Natoms;atomi++){
    if (fabs(particle[atomi].charge) > 1.5){
      cerr<<"Error in Ecoulomb!!! charge on atom["<<atomi<<"]"<<" is "<<particle[atomi].charge<<", bigger then 1.5 - ending session"<<endl;
      exit(0);
    }
  }
  //#pragma omp for reduction(+:Ecoulomb)
  for(atomi=N_init_ILP ; atomi < Natoms ; atomi++){
    for(j=0 ; j < NNeighb[atomi] ; j++){
      atomj = Neighb[atomi*MaxNeighb + j] ;
      r_ij = R_ij(particle,atomi,atomj, L);
      //r_ij   = R_ij_(particle[atomi].r[0],particle[atomj].r[0],particle[atomi].r[1],particle[atomj].r[1],particle[atomi].r[2],particle[atomj].r[2],L);
      //r_ij = rij_long[atomi*MaxNeighb + j];
      
      gamma = BondParams[particle[atomi].type*NAtomParams + particle[atomj].type].original_gamma;
      shield    = pow(r_ij,3.0) + pow(1.0/gamma,3.0);
      rijShield = 1.0/pow(shield,0.333333333333);
      Tap = calc_Tap(r_ij);
      
      //cerr<<"Ecoulomb= "<<Ecoulomb<<"rijShield= "<<rijShield<<"particle[atomj].charge= "<<particle[atomj].charge<<"particle[atomi].charge= "<<particle[atomi].charge<<endl;
      Ecoulomb += Tap*(particle[atomi].charge * particle[atomj].charge * rijShield);
      //Ecoulomb += (particle[atomi].charge * particle[atomj].charge/rij[atomi*Natoms + atomj]);
      
    }
  }
  //}
  //Ecoulomb *= (Kappa * eV2kcalmole);
  Ecoulomb *= (Kappa);
  
  //return(Ecoulomb);
}

void calc_E_Tersoff_Pot(double *angle,int *BondNum,int Natoms,point L,BondParamsStruct *BondParams,int *Bonded,atom *particle,double *Fc,int *NBond,double &E_Tersoff_,int *NBond_once,int *Bonded_once,int N_init)
{
  int i,j,k,ii,jj,kk,inum,jnum;
  int itype,jtype,ktype,iparam_ij,iparam_ijk;
  //tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij,prefactor;
  double c1,c2,c3,c4,bigr,biga,cut,cutsq,gamma,d,c,h,beta,bigb,bigd,powern,lam1,lam2,lam3;
  int powermint,atomi,itag,atomj,jtag,atomk,ktag;
  double shift_;
  //int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  
  for (atomi = N_init; atomi < Natoms; atomi++) {
   
    itag=atomi;
    itype = particle[atomi].type;
   
    xtmp = particle[atomi].r[0];
    ytmp = particle[atomi].r[1];
    ztmp = particle[atomi].r[2];
   
    // two-body interactions, skip half of them
    
    //jlist = firstneigh[i];
    //jnum = numneigh[i];
    for (j = 0; j < NBond_once[atomi]; j++) {
      //for (j = 0; j < NBond[atomi]; j++) {
      atomj = Bonded_once[atomi*MaxNBond + j];
      
      jtag = atomj;
      
      jtype = particle[atomj].type;
            
      powern  = BondParams[itype*NAtomParams + jtype].Tersoff_n;
      bigr    = BondParams[itype*NAtomParams + jtype].Tersoff_R;
      bigd    = BondParams[itype*NAtomParams + jtype].Tersoff_D;
      biga    = BondParams[itype*NAtomParams + jtype].Tersoff_A;
      cut   = bigr + bigd;
      cutsq = cut*cut;
      c1 = pow(2.0*powern*1.0e-16,-1.0/powern);
      c2 = pow(2.0*powern*1.0e-8,-1.0/powern);
      c3 = 1.0/c2;
      c4 = 1.0/c1;
      lam1  = BondParams[itype*NAtomParams + jtype].Tersoff_lambda1;
      delx = R_PBC(particle,atomi,atomj,0,L.x);
      dely = R_PBC(particle,atomi,atomj,1,L.y);	
      delz = R_PBC(particle,atomi,atomj,2,L.z);
      
      rsq = delx*delx + dely*dely + delz*delz;
      if(itype == 0 && jtype == 0)rsq += shift*shift + 2*sqrt(rsq)*shift;
      //iparam_ij = elem2param[itype][jtype][jtype];
      
      if (rsq > cutsq) continue;
      
      repulsive_pot(rsq,evdwl,lam1,biga,bigr,bigd);
           
      //printf("atomi= %i atomj=%i evdw1= %.16f fpair= %.16f F_REBO[atomi].x= %.16f F_REBO[atomj].x= %.16f rsq=%.16f \n",atomi,atomj,evdw1,fpair,F_REBO[atomi].x,F_REBO[atomj].x,rsq);
      //printf("itype= %i jtype=%i F_REBO[atomi].y= %.16f F_REBO[atomj].y= %.16f F_REBO[atomi].z= %.16f F_REBO[atomj].z= %.16f\n",itype,jtype,F_REBO[atomi].y,F_REBO[atomj].y,F_REBO[atomi].z,F_REBO[atomj].z);
      //if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
    }
    for (j = 0; j < NBond[atomi]; j++) {
      //for (jj = 0; jj < jnum; jj++) {
      //atomj = Bonded_once[atomi*MaxNBond + j];
      atomj = Bonded[atomi*MaxNBond + j];
      jtag = atomj;
      //jtag = tag[j];
      jtype = particle[atomj].type;
      
      delr1[0] = R_PBC(particle,atomj,atomi,0,L.x);
      delr1[1] = R_PBC(particle,atomj,atomi,1,L.y);	
      delr1[2] = R_PBC(particle,atomj,atomi,2,L.z);
            
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if(itype == 0 && jtype == 0)rsq1 += shift*shift + 2*sqrt(rsq1)*shift;
      bigr    = BondParams[itype*NAtomParams + jtype].Tersoff_R;
      bigd    = BondParams[itype*NAtomParams + jtype].Tersoff_D;
      cut   = bigr + bigd;
      cutsq = cut*cut;
      
      if (rsq1 > cutsq) continue;

      // accumulate bondorder zeta for each i-j interaction via loop over k

      zeta_ij = 0.0;
      
      //for (k = 0; k < NBond_once[atomi]; k++) {
      for (k = 0; k < NBond[atomi]; k++) {
	//for (kk = 0; kk < jnum; kk++) {
        if (j == k) continue;
	//atomk = Bonded_once[atomi*MaxNBond + k];
	atomk = Bonded[atomi*MaxNBond + k];
	ktag = atomk;
	//jtag = tag[j];
	ktype = particle[atomk].type;
	
     	delr2[0] = R_PBC(particle,atomk,atomi,0,L.x);
	delr2[1] = R_PBC(particle,atomk,atomi,1,L.y);	
	delr2[2] = R_PBC(particle,atomk,atomi,2,L.z);
	
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
	if(itype == 0 && jtype == 0){
	  rsq2 += shift*shift + 2*sqrt(rsq2)*shift;
	  shift_=shift;
	}
	else shift_=0.0;
	
	bigr    = BondParams[itype*NAtomParams + ktype].Tersoff_R;
	bigd    = BondParams[itype*NAtomParams + ktype].Tersoff_D;
	cut   = bigr + bigd;
	cutsq = cut*cut;
	
        if (rsq2 > cutsq) continue;
	
	gamma     = 1.0;
	powermint = 3;
	lam1  = BondParams[itype*NAtomParams + ktype].Tersoff_lambda1;
	lam3  = BondParams[itype*NAtomParams + ktype].Tersoff_lambda3;
	c  = BondParams[itype*NAtomParams + ktype].Tersoff_c;
	h  = BondParams[itype*NAtomParams + ktype].Tersoff_h;
	d  = BondParams[itype*NAtomParams + ktype].Tersoff_d;
	
	zeta_ij += zeta(rsq1,rsq2,delr1,delr2,powermint,lam1,lam3,c,d,h,gamma, bigr,bigd,shift_);
	//zeta_ij += zeta(&params[iparam_ijk],rsq1,rsq2,delr1,delr2);
      }
      
      // pairwise force due to zeta
      
      //force_zeta(&params[iparam_ij],rsq1,zeta_ij,fpair,prefactor,eflag,evdwl);
      lam2  = BondParams[itype*NAtomParams + jtype].Tersoff_lambda2;
      beta    = BondParams[itype*NAtomParams + jtype].Tersoff_beta;
      powern  = BondParams[itype*NAtomParams + jtype].Tersoff_n;
      bigb    = BondParams[itype*NAtomParams + jtype].Tersoff_B;
      
      c1 = pow(2.0*powern*1.0e-16,-1.0/powern);
      c2 = pow(2.0*powern*1.0e-8,-1.0/powern);
      c3 = 1.0/c2;
      c4 = 1.0/c1;
      
      force_zeta_pot(rsq1,zeta_ij,prefactor,evdwl,c1,c2,c3,c4,powern,beta,bigr,bigb,bigd,lam2);
      //printf("atomi= %i atomj=%i evdw1= %.16f fpair= %.16f rsq=%.16f \n",atomi,atomj,evdw1,fpair,rsq);
      //printf("atomi= %i atomj= %i fpair= %.16f \n",atomi,atomj,fpair);
      //if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);

      // attractive term via loop over k
      //for (k = 0; k < NBond_once[atomi]; k++) {
     
    }
  }
  //exit(0);
  E_Tersoff_= evdwl;
}
void repulsive_pot(double rsq, double &eng,double lam1,double biga,double bigr,double bigd)
{
  double r,tmp_fc,tmp_fc_d,tmp_exp;
  
  r        = sqrt(rsq);
  tmp_fc   = ters_fc(r,bigr,bigd);//ters_fc(r,param);
  //tmp_fc_d = ters_fc_d(r,bigr,bigd); //ters_fc_d(r,param);
  tmp_exp  = exp(-lam1 * r);
  //fforce   = -biga * tmp_exp * (tmp_fc_d - tmp_fc*lam1) / r;
  //fforce   = -biga * tmp_exp * (tmp_fc_d - tmp_fc*lam1);
  eng      += tmp_fc * biga * tmp_exp;
  
  //printf("r= %.16f tmp_fc= %.16f tmp_fc_d= %.16f tmp_exp= %.16f fforce= %.16f eng= %.16f\n",r,tmp_fc,tmp_fc_d,tmp_exp,fforce,eng);
}


/* ---------------------------------------------------------------------- */

void force_zeta_pot(double rsq, double zeta_ij,double &prefactor, double &eng,double c1,double c2,double c3,double c4,double powern,double beta,double bigr,double bigb,double bigd,double lam2)
{
  double r,fa,fa_d,bij;
  
  r = sqrt(rsq);
  //fa = ters_fa(r,param);
  fa= ters_fa(r,bigr,bigb,bigd,lam2);
  //fa_d = ters_fa_d(r,param);
  //fa_d = ters_fa_d(r,bigr,bigd,lam2,bigb);
  //bij = ters_bij(zeta_ij,param);
  bij =  ters_bij(zeta_ij,c1,c2,c3,c4,powern, beta);
  //fforce = 0.5*bij*fa_d / r;
  //prefactor = -0.5*fa * ters_bij_d(zeta_ij,param);
  //prefactor = -0.5*fa *  ters_bij_d(zeta_ij,c1,c2,c3,c4,beta,powern);
  //if (eflag) eng = 0.5*bij*fa;
  eng += 0.5*bij*fa;
  //printf("r= %.16f fa= %.16f fa_d= %.16f bij= %.16f fforce= %.16f prefactor= %.16f eng= %.16f\n",r,fa,fa_d,bij,fforce,prefactor,eng);
  
}
void calc_RI(int Natoms, double *rij_long, atom *particle, BondParamsStruct *BondParams,Normal_struct *Normal,int *Neighb,int *NNeighb,point L,double &RI,int *Normal_atom,int layer1)
{
  
  double EvdW;
  double SNN, SBB, SNB,SBC,SNC,SCC;
  double SNN_Tmp, SBB_Tmp, SNB_Tmp;
  double SNN2,SNB2,SBB2,SNN1,SNB1,SBB1,RegInd;
  int NClayer1,NBlayer1,NNlayer1,flag=0;
  int atomi, atomj, typei, typej, index,j,k,atomk,typek;
  double lambdaW, alpha, rvdW, Epsilon, f13, term,Inner_shield;
  double ecore, acore, rcore, Tap, r_ij,C6,C12;
  double gamma,Pij,Pji,angle,rx,ry,rz;
  double inside_acos,r_eff;
  double R,r_kj;
  double rx_ij,ry_ij,rz_ij,rx_ji,ry_ji,rz_ji,rx_kj,ry_kj,rz_kj,rx_jk,ry_jk,rz_jk,inside_acos_ij,inside_acos_ji,Pij_term1,Pji_term1,Pij_term,Pji_term,dr_ji_drk;
  double Interlayer_dis_Tmp,Interlayer_dis,nterlayer_dis;
  double Intralayer_dis,r_ik,Norj_L,Nori_L,Nork_L,Pkj,Pjk,inside_acos_kj,inside_acos_jk;
  double SNC_Tmp, SBC_Tmp,SNC1, SBC1,SNC2,SBC2,SCC_Tmp,SCC1,SCC2;
  double rC,rB,rN;
  int LRI_type,tmp_typei,tmp_typej;
  EvdW = 0;
  NClayer1=NBlayer1=NNlayer1=0;
  SNN = SBB = SNB = SBC = SNC = SCC = 0.0;
  
  for(atomi=0;atomi<Natoms;atomi++){
    typei = particle[atomi].type;
    if(particle[atomi].layer == layer1){
      if(typei == 0)NClayer1++;
      else if(typei == 5)NBlayer1++;
      else if(typei == 7)NNlayer1++;
    }
    //cerr<<"atomi= "<<atomi<<endl;
    //cerr<<"typei= "<<typei<<endl;
    //cerr<<"layeri= "<<particle[atomi].layer<<endl;
    for(j=0 ; j < NNeighb[atomi] ; j++){
      atomj = Neighb[atomi*MaxNeighb + j] ;
      if(particle[atomj].layer!=particle[atomi].layer){
	tmp_typej = particle[atomj].type;
	tmp_typei = particle[atomi].type;
      }
      //cerr<<"atomj= "<<atomj<<endl;
      //cerr<<"typej= "<<typej<<endl;
      //cerr<<"layerj= "<<particle[atomj].layer<<endl;
    }
  }
  //cerr<<"NBlayer1= "<<NBlayer1<<endl;
  //cerr<<"NNlayer1= "<<NNlayer1<<endl;
  //cerr<<"NClayer1= "<<NClayer1<<endl;
  if(tmp_typei==0 && tmp_typej==0){ //homogeneous graphene
    LRI_type=0;
    rC=0.5*LCC;
  }
  else if(tmp_typei==0 || tmp_typej==0){ //heterogeneous graphene/hBN
    LRI_type=1;
    rC=0.5*LCC;
    rB=0.2*LNB;
    rN=0.4*LNB;
  }
  else if((tmp_typei==5 && tmp_typej==7) || (tmp_typei==7 && tmp_typej==5) || (tmp_typei==5 && tmp_typej==5) || (tmp_typei==7 && tmp_typej==7)){//homogeneous hBN
    LRI_type=2;
    rB=0.15*LNB;
    rN=0.5*LNB;
  }
  else {
    cerr<<"error in calc_RI typei= "<<tmp_typei<<" typej= "<<tmp_typej<<" atom types"<<endl;
    exit(0);
  }
  //cerr<<"LRI_type= "<<LRI_type<<"typei= "<<tmp_typei<<"typej= "<<tmp_typej<<endl;
  
  SNC1=0.0;
  SBC1=CircleOverlap(0,rB,rC);
  SNC2=CircleOverlap(0,rN,rC);
  SBC2=CircleOverlap(0,rB,rC);
  SNB1=2.0*CircleOverlap(0,rN,rB);
  SBB1=0.0;
  SNN1=0.0;
  SNB2=0.0;
  SBB2=CircleOverlap(0,rB,rB);
  SNN2=CircleOverlap(0,rN,rN);
  SCC2 = 2.0*CircleOverlap(0,rC,rC);
  //SCC2 = 2.0*CircleOverlap(0,rC,rC);
  SCC1 = 2.0*CircleOverlap(2.0*sin(15*PIE/180)*LCC,rC,rC);
  
  //#pragma omp for reduction(+:EvdW)
  for(atomi=0 ; atomi < Natoms ; atomi++){
    typei = particle[atomi].type;
    Interlayer_dis = 1000;
    SNC_Tmp= SBC_Tmp= SCC_Tmp = SNB_Tmp = SNN_Tmp = SBB_Tmp =0.0;
    for(j=0 ; j < NNeighb[atomi] ; j++){
      atomj = Neighb[atomi*MaxNeighb + j] ;
      Intralayer_dis=0.0;
      typej   = particle[atomj].type;
      index   = typei*NAtomParams + typej;
      //r_ij    = R_ij(particle,atomi,atomj, L);
      
      rx_ij = R_PBC(particle,atomi,atomj,0,L.x);
      ry_ij = R_PBC(particle,atomi,atomj,1,L.y);	
      rz_ij = R_PBC(particle,atomi,atomj,2,L.z);
      
      rx_ji = -rx_ij;
      ry_ji = -ry_ij;
      rz_ji = -rz_ij;
      
      r_ij = sqrt(rx_ij*rx_ij + ry_ij*ry_ij + rz_ij*rz_ij);
      
      Nori_L=rx_ij*Normal[atomi].x + ry_ij*Normal[atomi].y + rz_ij*Normal[atomi].z;
      Norj_L=rx_ji*Normal[atomj].x + ry_ji*Normal[atomj].y + rz_ji*Normal[atomj].z;
      
      inside_acos_ij = (Nori_L) / (r_ij);
      inside_acos_ji = (Norj_L) / (r_ij);
      
      if((r_ij*r_ij-sqr(Nori_L)) <= 0)Pij=0.0;
      else Pij=sqrt(r_ij*r_ij - sqr(Nori_L));
      
      if((r_ij*r_ij-sqr(Norj_L)) <= 0)Pji=0.0;
      else Pji=sqrt(r_ij*r_ij - sqr(Norj_L));
      
      if(r_ij < 5.0 ){
	Interlayer_dis_Tmp = (sqrt(r_ij*r_ij -  Pij*Pij ) + sqrt(r_ij*r_ij -  Pji*Pji ))/2;
	if(Interlayer_dis > Interlayer_dis_Tmp)Interlayer_dis=Interlayer_dis_Tmp;
      }
      
      if(typei==0  && typej==0){
	SCC_Tmp += 0.5*(CircleOverlap(Pij,rC,rC)+CircleOverlap(Pji,rC,rC));
	SCC+=0.5*(CircleOverlap(Pij,rC,rC)+CircleOverlap(Pji,rC,rC));
      }
      else if((typei==7  && typej==0) || (typei==0  && typej==7)){
	SNC_Tmp += 0.5*(CircleOverlap(Pij,rC,rN)+CircleOverlap(Pji,rC,rN));
	SNC+=0.5*(CircleOverlap(Pij,rC,rN)+CircleOverlap(Pji,rC,rN));
      }
      else if((typei==5  && typej==0) || (typei==0  && typej==5)){
	SBC_Tmp += 0.5*(CircleOverlap(Pij,rB,rC)+CircleOverlap(Pji,rB,rC));
	SBC+=0.5*(CircleOverlap(Pij,rB,rC)+CircleOverlap(Pji,rB,rC));
      }
      else if(typei==7  && typej==7){
	SNN_Tmp += 0.5*(CircleOverlap(Pij,rN,rN)+CircleOverlap(Pji,rN,rN));
	SNN+=0.5*(CircleOverlap(Pij,rN,rN)+CircleOverlap(Pji,rN,rN));
      }
      else if(typei==5  && typej==5){
	SBB_Tmp += 0.5*(CircleOverlap(Pij,rB,rB)+CircleOverlap(Pji,rB,rB));
	SBB+= 0.5*(CircleOverlap(Pij,rB,rB)+CircleOverlap(Pji,rB,rB));
      }
      else if((typei==7 && typej==5) || (typei==5 && typej==7) ) {
	SNB_Tmp += 0.5*(CircleOverlap(Pij,rN,rB)+CircleOverlap(Pji,rN,rB));
	SNB+= 0.5*(CircleOverlap(Pij,rN,rB)+CircleOverlap(Pji,rN,rB));
      }
      else{
	cerr<<"Wrong atom type! Ending session.\n";
	exit(0);
      }
            
      //cerr<<"SCC_Tmp= "<<SCC_Tmp<<endl;
      for(k=0; k < 3; k++){
	atomk =Normal_atom[atomi*3 + k];
	typek = particle[atomk].type;
	r_ik  = R_ij(particle,atomi,atomk, L);
	rx_kj = R_PBC(particle,atomk,atomj,0,L.x);
	ry_kj = R_PBC(particle,atomk,atomj,1,L.y);	
	rz_kj = R_PBC(particle,atomk,atomj,2,L.z);
	
	rx_jk = -rx_kj;
	ry_jk = -ry_kj;
	rz_jk = -rz_kj;
	
	r_kj = sqrt(rx_kj*rx_kj + ry_kj*ry_kj + rz_kj*rz_kj);
	
	Nork_L=rx_kj*Normal[atomk].x + ry_kj*Normal[atomk].y + rz_kj*Normal[atomk].z;
	Norj_L=rx_jk*Normal[atomj].x + ry_jk*Normal[atomj].y + rz_jk*Normal[atomj].z;
	
	inside_acos_kj = (Nork_L) / (r_kj);
	inside_acos_jk = (Norj_L) / (r_kj);
	
	if((r_kj*r_kj-sqr(Nork_L)) <= 0)Pkj=0.0;
	else Pkj=sqrt(r_kj*r_kj - sqr(Nork_L));
	
	if((r_kj*r_kj-sqr(Norj_L)) <= 0)Pjk=0.0;
	else Pjk=sqrt(r_kj*r_kj - sqr(Norj_L));
	
	Intralayer_dis+=r_ik/3.0;
	
	if(typek==0  && typej==0)      SCC_Tmp += 1.0/3.0*0.5*(CircleOverlap(Pkj,rC,rC)+CircleOverlap(Pjk,rC,rC));
	else if(typek==7  && typej==0) SNC_Tmp += 1.0/3.0*0.5*(CircleOverlap(Pkj,rC,rN)+CircleOverlap(Pjk,rC,rN));
	else if(typek==0 && typej==7)  SNC_Tmp += 1.0/3.0*0.5*(CircleOverlap(Pkj,rC,rN)+CircleOverlap(Pjk,rC,rN));
	else if(typek==5 && typej==0)  SBC_Tmp += 1.0/3.0*0.5*(CircleOverlap(Pkj,rB,rC)+CircleOverlap(Pjk,rB,rC));
	else if(typek==0 && typej==5)  SBC_Tmp += 1.0/3.0*0.5*(CircleOverlap(Pkj,rB,rC)+CircleOverlap(Pjk,rB,rC));
	else if(typek==7 && typej==7)  SNN_Tmp += 1.0/3.0*0.5*(CircleOverlap(Pkj,rN,rN)+CircleOverlap(Pjk,rN,rN));
	else if(typek==5 && typej==5)  SBB_Tmp += 1.0/3.0*0.5*(CircleOverlap(Pkj,rB,rB)+CircleOverlap(Pjk,rB,rB));
	else if(typek==7 && typej==5)  SNB_Tmp += 1.0/3.0*0.5*(CircleOverlap(Pkj,rN,rB)+CircleOverlap(Pjk,rN,rB));
	else if(typek==5 && typej==7)  SNB_Tmp += 1.0/3.0*0.5*(CircleOverlap(Pkj,rN,rB)+CircleOverlap(Pjk,rN,rB));
	else{
	  cerr<<"Wrong atom type! Ending session.\n";
	  exit(0);
	}
      }
      
    }
    
    //particle[atomi].Inter = -(2.0*(Interlayer_dis-3.3)/(3.55-3.3)-1);
    particle[atomi].Inter = Interlayer_dis;
    particle[atomi].Intra = Intralayer_dis;//graphene on h-BN
    
    if(LRI_type==0)particle[atomi].RI         = -(2.0*(SCC_Tmp-SCC1)/(SCC2-SCC1)-1.0) ; //graphene/graphene
    else if(LRI_type==1)particle[atomi].RI    = -(2.0*( (SNC_Tmp -SNC1) + (SBC_Tmp-SBC1) ) / ((SNC2-SNC1) + (SBC2-SBC1) ) -1.0);//graphene/h-BN
    else if(LRI_type==2)particle[atomi].RI    = -(2.0*( ((SNN_Tmp-SNN1) + (SBB_Tmp-SBB1) - (SNB_Tmp-SNB1)) / ( ((SNN2-SNN1) + (SBB2-SBB1) - (SNB2-SNB1)) )) - 1.0); //h-BN/h-BN
    else{
      cerr<<"Wrong LRI_type in cal_RI! Ending session.\n";
      exit(0);
    }
  }
  
  SNC1=0.0;
  SBC1=CircleOverlap(0,rB,rC)*(NBlayer1  + NClayer1*0.5);
  SNC2=CircleOverlap(0,rN,rC)*(NNlayer1  + NClayer1*0.5);
  SBC2=CircleOverlap(0,rB,rC)*(NBlayer1  + NClayer1*0.5);
  SNB1=CircleOverlap(0,rN,rB)*(NBlayer1  + NNlayer1);
  SBB1=0.0;
  SNN1=0.0;
  SNB2=0.0;
  SBB2=CircleOverlap(0,rB,rB)*(NBlayer1);
  SNN2=CircleOverlap(0,rN,rN)*(NNlayer1);
  SCC2 = CircleOverlap(0,rC,rC)*NClayer1;
  SCC1 = CircleOverlap(0,rC,rC)*NClayer1*0.5;
  
  if(LRI_type==0)RI    = (SCC-SCC1)/(SCC2-SCC1) ; //graphene/graphene
  else if(LRI_type==1)RI    = ((SNC -SNC1) + (SBC-SBC1) ) / ((SNC2-SNC1) + (SBC2-SBC1) );//graphene/h-BN
  else if(LRI_type==2)RI    =   ((SNN-SNN1) + (SBB-SBB1) - (SNB-SNB1)) /( (SNN2-SNN1) + (SBB2-SBB1) - (SNB2-SNB1) ) ; //h-BN/h-BN
  
  //return(EvdW);
  
}
void calc_rij_LRI(int Natoms, atom *particle, point L, double *rij,BondParamsStruct *BondParams,int *Neighb,int *NNeighb,double *rij_long,int *NList_long,int *List_long,int layer1, int layer2)
{
  int atomi, atomj,typei,typej,j;
  double rx, ry, rz,r_ij;
  double Tersoff_R,Tersoff_D;
  int NBondj, NBondi,NBondi_once;
  double cutoff;
  double dS,rsq;
  
  
  //#pragma omp parallel for private(atomj)
  for(atomi=0 ; atomi < Natoms ; atomi++){
    NNeighb[atomi] = 0;
    for(atomj=0 ; atomj < MaxNeighb ; atomj++){
      Neighb[atomi*MaxNeighb + atomj] = -1;
    }
  }
  
  //#pragma omp parallel for private(rx,ry,rz,typei,typej,Tersoff_R,Tersoff_D,atomj,r_ij,j,atomi,NBondj, NBondi, NBondi_once,dS,cutoff) 
  for(atomi=0 ; atomi < Natoms ; atomi++){
    
    //rij[atomi*Natoms + atomi] = 0.0;
    typei = particle[atomi].type;
    for(j=0 ; j < NList_long[atomi] ; j++){
      atomj=List_long[atomi*MaxList_long + j];
      r_ij=R_ij(particle,atomi,atomj,L);
      if((r_ij < non_bond_cut)){
	if(particle[atomj].layer==layer1 && particle[atomi].layer==layer2){
	  Neighb[atomj*MaxNeighb + NNeighb[atomj]]   = atomi;
	  NNeighb[atomj]++;
	}
	else if(particle[atomi].layer==layer1 && particle[atomj].layer==layer2){
	  Neighb[atomi*MaxNeighb + NNeighb[atomi]]   = atomj;
	  NNeighb[atomi]++;
	}
	if(NNeighb[atomi] == MaxNeighb){
	  cerr<<"Exceeded maximum number of neighbors for atom "<<atomi<<"! Ending session.\n";
	  exit(0);
	}
      }
    }      
  }
}

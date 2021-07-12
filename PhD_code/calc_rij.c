#include "pair_airebo.h"
#include "declarations.h"
#include "pair_tersoff.h"
/***********************************************************/
/* Calculate and fill the following arrays:                */
/*                                                         */
/* rij[Natoms*Natoms] - interatom distance array.          */
/*                                                         */
/* NBond[Natoms]      - number of atoms bonded to each     */
/*                      atom - two atoms are consedered to */
/*			be bonded if BO_ij != 0 namely if  */
/*                      their distance is less than 3A     */
/*			(BOCutoff).                        */
/*                                                         */
/* Continue explanation ...                                */
/***********************************************************/

void calc_rij(int *NBond, int *Bonded, int *BondNum, int Natoms, atom *particle, point L, double *rij,double *Fc,BondParamsStruct *BondParams,int *Neighb,int *NNeighb,double *rij_long,int *NList,int *List,int *NBond_once,int *Bonded_once,int potential_flag)
{
  
#pragma omp parallel
  {
    int atomi, atomj,typei,typej,j;
    double r_ij;
    double Tersoff_R,Tersoff_D;
    double cutoff;
        
    int *NBond_,*NBond_once_;
    int *Bonded_,*Bonded_once_;
    
    NBond_= new int[Natoms];
    NBond_once_=  new int[Natoms];
    Bonded_= new int[Natoms*MaxNBond];
    Bonded_once_= new int[Natoms*MaxNBond];
    
#pragma omp for  
    for(atomi=0 ; atomi < Natoms ; atomi++){
      NBond[atomi]   = 0;
      NBond_once[atomi]   = 0;
      for(atomj=0 ; atomj < MaxNBond ; atomj++){
	Bonded[atomi*MaxNBond + atomj] = -1;
	Bonded_once[atomi*MaxNBond + atomj] = -1;
      }
    }
    
    for(atomi=0 ; atomi < Natoms ; atomi++){
      NBond_[atomi]   = 0;
      NBond_once_[atomi]   = 0;
      for(atomj=0 ; atomj < MaxNBond ; atomj++){
	Bonded_[atomi*MaxNBond + atomj] = -1;
	Bonded_once_[atomi*MaxNBond + atomj] = -1;
      }
    }
    
#pragma omp for private(typei,typej,Tersoff_R,Tersoff_D,atomj,r_ij,j,atomi,cutoff) 
    for(atomi=0 ; atomi < Natoms ; atomi++){
      
      typei = particle[atomi].type;
      for(j=0 ; j < NList[atomi] ; j++){
	atomj=List[atomi*MaxList + j];
	r_ij  = R_ij(particle,atomi,atomj, L);
	//r_ij=R_ij_(particle[atomi].r[0],particle[atomj].r[0],particle[atomi].r[1],particle[atomj].r[1],particle[atomi].r[2],particle[atomj].r[2],L);
	
	typej=particle[atomj].type;
	
	Tersoff_R   = BondParams[typei*NAtomParams + typej].Tersoff_R;
	Tersoff_D   = BondParams[typei*NAtomParams + typej].Tersoff_D;
	
	if ((potential_flag == 0 || potential_flag == 1) && (particle[atomi].type == 0 || particle[atomi].type == 1)) cutoff = sqrt(rcmaxsq[typei][typej]);
	else{
	  cutoff = Tersoff_R+Tersoff_D;
	  if(particle[atomj].type == 0 && typei ==0)r_ij+=shift;
	}
	
	if(r_ij < cutoff){
		  
	  Bonded_once_[atomi*MaxNBond + NBond_once_[atomi]] = atomj;
	  Bonded_[atomi*MaxNBond + NBond_[atomi]] = atomj;
	  Bonded_[atomj*MaxNBond + NBond_[atomj]] = atomi;
	    
	  NBond_[atomi]++;
	  NBond_once_[atomi]++;
	  NBond_[atomj]++;
	}
      }
    }
    
#pragma omp critical
    {
      for(atomi=0 ; atomi < Natoms ; atomi++){
	NBond[atomi]+=NBond_[atomi];
	NBond_once[atomi]+=NBond_once_[atomi];
	
	if(NBond[atomi] == MaxNBond ){
	  cerr<<"Exceeded maximum number of bonded atoms for atom "<<atomi<<"! Ending session.\n";
	  exit(0);
	}
	
	for(atomj=(NBond[atomi]-NBond_[atomi]) ; atomj < NBond[atomi] ; atomj++) Bonded[atomi*MaxNBond + atomj] = Bonded_[atomi*MaxNBond + (atomj-(NBond[atomi]-NBond_[atomi]) )];
	for(atomj=(NBond_once[atomi]-NBond_once_[atomi]) ; atomj < NBond_once[atomi] ; atomj++) Bonded_once[atomi*MaxNBond + atomj] = Bonded_once_[atomi*MaxNBond + (atomj-(NBond_once[atomi]-NBond_once_[atomi]) )];
      }
    }
    
    
    
    delete [] NBond_;
    delete [] NBond_once_;
    delete [] Bonded_once_;
    delete [] Bonded_;
  }
}

void calc_rij_long(int Natoms, atom *particle, point L,int *Neighb,int *NNeighb){
#pragma omp parallel
  {
    int typei,typej,atomj,atomi,j;
    double r_ij;
    
    int *NNeighb_,*Neighb_;
    
    NNeighb_= new int[Natoms];
    Neighb_=  new int[Natoms*MaxNeighb];
    
    for(atomi=0 ; atomi < Natoms ; atomi++){
      NNeighb_[atomi]=0;
      for(atomj=0 ; atomj < MaxNeighb ; atomj++){
	Neighb_[atomi*MaxNeighb + atomj] = -1;
      }
    }
    
#pragma omp for
    for(atomi=0 ; atomi < Natoms ; atomi++){
      NNeighb[atomi]=0;
      for(atomj=0 ; atomj < MaxNeighb ; atomj++){
	Neighb[atomi*MaxNeighb + atomj] = -1;
      }
    }
    
    
#pragma omp for private(typei,typej,atomj,r_ij,j,atomi) 
    for(atomi=0 ; atomi < Natoms ; atomi++){
      
      typei = particle[atomi].type;
      for(j=0 ; j < NList_long[atomi] ; j++){
	atomj=List_long[atomi*MaxList_long + j];
	typej = particle[atomj].type;
	//r_ij=R_ij_(particle[atomi].r[0],particle[atomj].r[0],particle[atomi].r[1],particle[atomj].r[1],particle[atomi].r[2],particle[atomj].r[2],L);
	r_ij  = R_ij(particle,atomi,atomj, L);
	//#pragma omp critical
	//{				
	if((r_ij < non_bond_cut) && (fabs(particle[atomi].layer-particle[atomj].layer) <= N_ILP_Layers || fabs(particle[atomi].layer-particle[atomj].layer+N_layers) <= N_ILP_Layers || fabs(particle[atomi].layer-particle[atomj].layer-N_layers) <= N_ILP_Layers)){
	  //if(typej==5 || typej==7){
	  //Neighb_[atomj*MaxNeighb + NNeighb_[atomj]]   = atomi;
	  //NNeighb_[atomj]++;
	  //}
	  //else{
	  Neighb_[atomi*MaxNeighb + NNeighb_[atomi]]   = atomj;
	  NNeighb_[atomi]++;
	  //}
	  //}
	}
      }      
    }
    
#pragma omp critical
    {
      for(atomi=0 ; atomi < Natoms ; atomi++){
	NNeighb[atomi]+=NNeighb_[atomi];
	if(NNeighb[atomi] == MaxNeighb){
	  cerr<<"Exceeded maximum number of neighbors for atom "<<atomi<<"! Ending session.\n";
	  exit(0);
	}
	for(atomj=(NNeighb[atomi]-NNeighb_[atomi]) ; atomj < NNeighb[atomi] ; atomj++)Neighb[atomi*MaxNeighb + atomj] = Neighb_[atomi*MaxNeighb + (atomj-(NNeighb[atomi]-NNeighb_[atomi]) )];
      }
    }
    delete [] NNeighb_;
    delete [] Neighb_; 
  }
}

void calc_rij_short(int *NBond, int *Bonded, int *BondNum, int Natoms, atom *particle, point L, double *rij,double *Fc,BondParamsStruct *BondParams,int *Neighb,int *NNeighb,double *rij_long,int *NList,int *List,int *NBond_once,int *Bonded_once,int potential_flag)
{
  
#pragma omp parallel
  {
    int atomi, atomj,typei,typej,j;
    double r_ij;
    double Tersoff_R,Tersoff_D;
    int NBondj, NBondi,NBondi_once;
    double cutoff;
    double dS,rsq;
    
#pragma omp  for 
    for(atomi=0 ; atomi < Natoms ; atomi++){
      nC[atomi]=0;
      nH[atomi]=0;
    }
    
    
    double *nC_,*nH_;
    nC_= new double[Natoms];
    nH_= new double[Natoms];
    
    for(atomi=0 ; atomi < Natoms ; atomi++){
      nC_[atomi]=0.0;
      nH_[atomi]=0.0;
    }
    
#pragma omp for //private(typei,typej,Tersoff_R,Tersoff_D,atomj,r_ij,j,atomi,NBondj, NBondi, NBondi_once,dS,cutoff) 
    for(atomi=0 ; atomi < Natoms ; atomi++){
      
      typei = particle[atomi].type;
      for(j=0 ; j < NList[atomi] ; j++){
	atomj=List[atomi*MaxList + j];
	r_ij  = R_ij(particle,atomi,atomj, L);
	//r_ij=R_ij_(particle[atomi].r[0],particle[atomj].r[0],particle[atomi].r[1],particle[atomj].r[1],particle[atomi].r[2],particle[atomj].r[2],L);
	
	//if(r_ij < 3.0){
	
	typej=particle[atomj].type;
	
	Tersoff_R   = BondParams[typei*NAtomParams + typej].Tersoff_R;
	Tersoff_D   = BondParams[typei*NAtomParams + typej].Tersoff_D;
	
	if ((potential_flag == 0 || potential_flag == 1) && (particle[atomi].type == 0 || particle[atomi].type == 1)) cutoff = sqrt(rcmaxsq[typei][typej]);
	else cutoff = Tersoff_R+Tersoff_D;
	
	if(r_ij < cutoff){
	  
	  if ((potential_flag == 0 || potential_flag == 1) && (particle[atomi].type == 0 || particle[atomi].type == 1)){
	    if (typej == 0){
	      nC_[atomi] += Sp(r_ij,rcmin[typei][typej],rcmax[typei][typej],dS);
	      nC_[atomj] += Sp(r_ij,rcmin[typej][typei],rcmax[typej][typei],dS);
	    }
	    else{
	      nH_[atomi] += Sp(r_ij,rcmin[typei][typej],rcmax[typei][typej],dS);
	      nH_[atomj] += Sp(r_ij,rcmin[typej][typei],rcmax[typej][typei],dS);
	    }
	  }
	  else {
	    if((particle[atomi].type == 0 && particle[atomj].type == 0))
	      Fc[atomi*MaxNBond + NBond[atomi]]  = Fc[atomj*MaxNBond + NBond[atomj]]  = Fc_(r_ij+shift,Tersoff_R,Tersoff_D);
	    else Fc[atomi*MaxNBond + NBond[atomi]]  = Fc[atomj*MaxNBond + NBond[atomj]]  = Fc_(r_ij,Tersoff_R,Tersoff_D);
	  }
	  
	}
	//}
      }      
    }
    
#pragma omp critical
    {
      for(atomi=0 ; atomi < Natoms ; atomi++){
	nC[atomi]+=nC_[atomi];
	nH[atomi]+=nH_[atomi];
      }
    }
    
    delete [] nC_;
    delete [] nH_; 
  }
}




void PBC_Init(atom* particle,int Natoms,double &Xmin,double&Ymin,double&Zmin,double&Xmax,double&Ymax,double&Zmax,point L,point *R0){

  double Term;
  int atomi;
  point COM;
  
  COM.x=COM.y=COM.z=0;
  
  for (atomi=0 ; atomi < Natoms ; atomi++){

    COM.x += particle[atomi].r[0]/Natoms;
    COM.y += particle[atomi].r[1]/Natoms;
    COM.z += particle[atomi].r[2]/Natoms;

  }

  Xmin = COM.x-0.5*L.x;
  Xmax = COM.x+0.5*L.x;
  Ymin = COM.y-0.5*L.y;
  Ymax = COM.y+0.5*L.y;
  Zmin = COM.z-0.5*L.z;
  Zmax = COM.z+0.5*L.z;
  /*
  Term = 0.5*(L.x - (Xmax - Xmin));
  Xmax += Term;
  Xmin -= Term;
  Term = 0.5*(L.y - (Ymax - Ymin));
  Ymax += Term;
  Ymin -= Term;
  Term = 0.5*(L.z - (Zmax - Zmin));
  Zmax += Term;
  Zmin -= Term;
  */
}
double calc_val_angle(int atomi,int atomj, int atomk, int Natoms, int *NBond, int *Bonded,double r_jk, double r_ij, atom *particle,point L,dr_ij &rji,dr_ij &rjk,double inside_acos)
{
  // This routine calculates and fills the angles array.  For each
  // atom j all the neighboring atom pairs (i and k) are considered
  // and the angle ijk is calculated, where j is the central atom. The
  // calculation is done using the scalar product formula and the acos
  // function.
  //#pragma omp parallel
  //{
  double rx,ry,rz,angle;
  
  angle=0;
  
  //point rji, rj;k
  
  rji.r[0] = R_PBC(particle,atomi,atomj,0,L.x);
  rji.r[1] = R_PBC(particle,atomi,atomj,1,L.y);
  rji.r[2] = R_PBC(particle,atomi,atomj,2,L.z);
  
  rjk.r[0] = R_PBC(particle,atomk,atomj,0,L.x);
  rjk.r[1] = R_PBC(particle,atomk,atomj,1,L.y);
  rjk.r[2] = R_PBC(particle,atomk,atomj,2,L.z);
  
  inside_acos = ((rji.r[0]*rjk.r[0] + rji.r[1]*rjk.r[1] + rji.r[2]*rjk.r[2]) / (r_ij* r_jk));
  
  if (fabs(inside_acos) < 1.0){
    angle = acos((rji.r[0]*rjk.r[0] + rji.r[1]*rjk.r[1] + rji.r[2]*rjk.r[2]) / (r_ij* r_jk));
  }
  else if (inside_acos >= 1.0) angle = 0.0;
  else angle = PIE;

  return(angle);
}
/*********************************************************************************************************/

/*********************************************************************************************************/

void calc_Normal(int Natoms,int *Normal_atom,atom *particle,Normal_struct *Normal,point L)
{
  //#pragma omp parallel
  //{
  int atomi,atomj,atomk,atoml;
  double vector_1x,vector_1y,vector_1z;
  double vector_2x,vector_2y,vector_2z;
  double Normal_length;
  
  //#pragma omp for
#pragma omp parallel for private(atomi,atomk,atomj,vector_1x,vector_1y,vector_1z,vector_2x,vector_2y,vector_2z,Normal_length)
  for(atomi=0; atomi < Natoms; atomi++){
    Normal[atomi].x = 0.0;
    Normal[atomi].y = 0.0;
    Normal[atomi].z = 0.0;
    
    atomj=Normal_atom[atomi*3 + 1];
    atomk=Normal_atom[atomi*3 + 0];
    atoml=Normal_atom[atomi*3 + 2];
    
    Normal[atomi].vector_1x = R_PBC(particle,atomk,atomj,0,L.x);
    Normal[atomi].vector_1y = R_PBC(particle,atomk,atomj,1,L.y);
    Normal[atomi].vector_1z = R_PBC(particle,atomk,atomj,2,L.z);
    
    Normal[atomi].vector_2x = R_PBC(particle,atomk,atoml,0,L.x);
    Normal[atomi].vector_2y = R_PBC(particle,atomk,atoml,1,L.y);
    Normal[atomi].vector_2z = R_PBC(particle,atomk,atoml,2,L.z);
    
    Normal[atomi].not_norm_x = Normal[atomi].vector_1y*Normal[atomi].vector_2z - Normal[atomi].vector_2y*Normal[atomi].vector_1z;
    Normal[atomi].not_norm_y = Normal[atomi].vector_1z*Normal[atomi].vector_2x - Normal[atomi].vector_1x*Normal[atomi].vector_2z;
    Normal[atomi].not_norm_z = Normal[atomi].vector_1x*Normal[atomi].vector_2y - Normal[atomi].vector_1y*Normal[atomi].vector_2x;
    
    Normal[atomi].Normal_length   = sqrt(sqr(Normal[atomi].not_norm_x)+sqr(Normal[atomi].not_norm_y)+sqr(Normal[atomi].not_norm_z));
    
    Normal[atomi].x = Normal[atomi].not_norm_x/(Normal[atomi].Normal_length);
    Normal[atomi].y = Normal[atomi].not_norm_y/(Normal[atomi].Normal_length);
    Normal[atomi].z = Normal[atomi].not_norm_z/(Normal[atomi].Normal_length);
    
  }
  
  //}
}
void calc_dNormal_k(int &atomi,int &atomn,int &dir,dr_ij *dN,atom *particle,Normal_struct *Normal,int *Normal_atom,point &L)
{
  double vector_1x,vector_1y,vector_1z;
  double vector_2x,vector_2y,vector_2z;
  int atomj,atomk,atoml;
  double Normal_length,d_Norm_length;
  double Normalx,Normaly,Normalz;
  
  atomj=Normal_atom[atomi*3 + 0];
  atomk=Normal_atom[atomi*3 + 1];
  atoml=Normal_atom[atomi*3 + 2];
  
  if(particle[atomi].type==1){
    dN[0].r[0] = 0.0;
    dN[0].r[1] = 0.0;
    dN[0].r[2] = 0.0;
  }
  else {
    if (atomn == atomj){
      if (dir == 0){
	dN[0].r[0] = 0.0;
	dN[0].r[1] = -Normal[atomi].vector_2z;
	dN[0].r[2] = Normal[atomi].vector_2y;
      }
      else if(dir == 1){
	dN[0].r[0] = Normal[atomi].vector_2z;
	dN[0].r[1] = 0.0;
	dN[0].r[2] =-Normal[atomi].vector_2x;
      }
      else if(dir == 2){
	dN[0].r[0] =-Normal[atomi].vector_2y;
	dN[0].r[1] = Normal[atomi].vector_2x;
	dN[0].r[2] = 0.0;
      }
    }
    else if (atomn == atomk){
      if (dir == 0){
	dN[0].r[0] = 0.0;
	dN[0].r[1] =-Normal[atomi].vector_1z + Normal[atomi].vector_2z;
	dN[0].r[2] =-Normal[atomi].vector_2y + Normal[atomi].vector_1y;
      }
      else if(dir == 1){
	dN[0].r[0] = -Normal[atomi].vector_2z + Normal[atomi].vector_1z;
	dN[0].r[1] = 0.0;
	dN[0].r[2] = -Normal[atomi].vector_1x + Normal[atomi].vector_2x;
      }
      else if(dir == 2){
	dN[0].r[0] =-Normal[atomi].vector_1y + Normal[atomi].vector_2y;
	dN[0].r[1] =-Normal[atomi].vector_2x + Normal[atomi].vector_1x;
	dN[0].r[2] = 0.0;
      }
    }
    else if (atomn == atoml){
      if (dir == 0){
	dN[0].r[0] = 0.0;
	dN[0].r[1] =Normal[atomi].vector_1z;
	dN[0].r[2] =-Normal[atomi].vector_1y;
      }
      else if(dir == 1){
	dN[0].r[0] =-Normal[atomi].vector_1z;
	dN[0].r[1] = 0.0;
	dN[0].r[2] =Normal[atomi].vector_1x;
      }
      else if(dir == 2){
	dN[0].r[0] =Normal[atomi].vector_1y;
	dN[0].r[1] =-Normal[atomi].vector_1x;
	dN[0].r[2] = 0.0;
      }
    }
    else{
      dN[0].r[0] =0.0;
      dN[0].r[1] =0.0;
      dN[0].r[2] =0.0;
    }
    d_Norm_length = pow(Normal[atomi].Normal_length,-3.0)*(Normal[atomi].not_norm_x*dN[0].r[0] + Normal[atomi].not_norm_y*dN[0].r[1] + Normal[atomi].not_norm_z*dN[0].r[2]);
    dN[0].r[0]= dN[0].r[0]/Normal[atomi].Normal_length-Normal[atomi].not_norm_x*d_Norm_length;
    dN[0].r[1]= dN[0].r[1]/Normal[atomi].Normal_length-Normal[atomi].not_norm_y*d_Norm_length;
    dN[0].r[2]= dN[0].r[2]/Normal[atomi].Normal_length-Normal[atomi].not_norm_z*d_Norm_length;
  }
  
}
void set_Normal_atom(double *rij,int Natoms,Normal_struct *Normal,atom *particle,int *Normal_atom,point L,int *NList,int *List)
{
  //This function fills the Normal_atom array of the 3 nearest neighbors surrounding an atom which we use to calculate its normal. 
  double rx,ry,rz;
  int n,atomj,atomi,j;
  double r_ij,R;

  n=0;
  R=Normal_cuttoff;
  
  for(atomi=0; atomi < Natoms; atomi++){
    //R=Normal_cuttoff;
    //n=0;
    for(atomj = 0;atomj < Natoms;atomj++){
    //cerr<<"NList[atomi]= "<<NList[atomi]<<endl;
    //for(j=0 ; j < NList[atomi] ; j++){
      //atomj=List[atomi*MaxList + j];

      //r_ij=R_ij_(particle[atomi].r[0],particle[atomj].r[0],particle[atomi].r[1],particle[atomj].r[1],particle[atomi].r[2],particle[atomj].r[2],L);
      
      if (atomi == atomj || (particle[atomi].layer!=particle[atomj].layer))continue;
      
      //r_ij=R_ij(particle,atomi,atomj, L); //warning!!! gives nans ?? found the problem
      rx = R_PBC(particle,atomi,atomj,0,L.x);
      ry = R_PBC(particle,atomi,atomj,1,L.y);
      rz = R_PBC(particle,atomi,atomj,2,L.z);
      
      r_ij = sqr(rx) + sqr(ry) + sqr(rz);
      
      if(r_ij > R*R)continue;
      else
	{
	  Normal_atom[atomi*3+n] = atomj;
	  n++;
	}
      
      if(n==3)break;
    
    }
    
    if(n!=3){
      
      atomi-=1;
      R = R + 0.2;
    }
    else R=Normal_cuttoff;
    
    n=0;
       
  }
  
}
/*
double R_ij(atom *particle,int atomi,int atomj, point L){
  
  double rx,ry,rz;
  
  //rx = R_PBC_(particle[atomi].r[0],particle[atomj].r[0],L.x);
  //ry = R_PBC_(particle[atomi].r[1],particle[atomj].r[1],L.y);
  //rz = R_PBC_(particle[atomi].r[2],particle[atomj].r[2],L.z);

  rx = R_PBC(particle,atomi,atomj,0,L.x);
  ry = R_PBC(particle,atomi,atomj,1,L.y);
  rz = R_PBC(particle,atomi,atomj,2,L.z);
  
  return(sqrt(sqr(rx) + sqr(ry) + sqr(rz)));
  }
*/
/*
double R_PBC(atom *particle,int atomi,int atomj,int dir,double L){
  double R;
  R = particle[atomi].r[dir] - particle[atomj].r[dir];
  
  R -= L * rint(R / (L + TINY));
  //R -= L * rint(R / (L));
  
  return(R);
  }

double Fc_(double r_ij,double Tersoff_R,double Tersoff_D)
{
  double Fc;
  
  if(r_ij < Tersoff_R+Tersoff_D){
    if(r_ij < Tersoff_R-Tersoff_D){
      Fc = 1;
    }
    else{
      Fc = 1/2-1/2*sin(PIE*0.5*(r_ij - Tersoff_R)/Tersoff_D);
    }
  }
  else Fc=0;
  
  return(Fc);
}
*/
void UpdateNeighbList(int Dim, atom *particle, int *List, int *NList, point *R0, point L)
{
  
#pragma omp parallel 
  {
    int i, j;
    double R, rx, ry, rz;
    int *List_long_;
    int *NList_long_;
    int *List_;
    int *NList_;
    
      List_long_  = new int[Dim*MaxList_long];
      NList_long_ = new int[Dim];
      List_      = new int[Dim*MaxList];
      NList_      = new int[Dim];
    
#pragma omp for
    for(i=0 ; i < Dim ;i++) NList[i] = NList_long[i] =0;
#pragma omp for 
    for(i=0 ; i < Dim*MaxList ;i++) List[i] = -1;
#pragma omp for 
    for(i=0 ; i < Dim*MaxList_long ;i++) List_long[i] = -1;
    
    for(i=0 ; i < Dim ;i++) NList_[i] = NList_long_[i] =0;
    for(i=0 ; i < Dim*MaxList ;i++) List_[i] = -1;
    for(i=0 ; i < Dim*MaxList_long ;i++) List_long_[i] = -1;
    
#pragma omp for private(j,i,rx,ry,rz,R)
    for(i=0 ; i < Dim ;i++){
      for(j=0 ; j < i ; j++){
	rx = R_PBC(particle,i,j,0,L.x);
	ry = R_PBC(particle,i,j,1,L.y);
	rz = R_PBC(particle,i,j,2,L.z);
	
	R  = sqrt(sqr(rx) + sqr(ry) + sqr(rz));
	
	if(R < Rcut_List_long && (particle[i].layer != particle[j].layer)){
	  List_long_[i*MaxList_long + NList_long_[i]] = j;
	  NList_long_[i]++;
	  //List_long[i*MaxList_long + NList_long[i]] = j;
	  //NList_long[i]++;
	}
	else if (R < Rcut_List_short && (particle[i].layer == particle[j].layer)){
	  List_[i*MaxList + NList_[i]] = j;
	  NList_[i]++;
	  //List[i*MaxList + NList[i]] = j;
	  //NList[i]++;
	}
      }
      R0[i].x = particle[i].r[0];
      R0[i].y = particle[i].r[1];
      R0[i].z = particle[i].r[2];
    }
    
#pragma omp critical
    {
      for(i=0 ; i < Dim ;i++) {
	NList[i] += NList_[i];
	NList_long[i] += NList_long_[i];
	for(j= (NList[i]- NList_[i]); j < NList[i];j++) List[i*MaxList + j] = List_[i*MaxList + j-(NList[i]- NList_[i])];
	for(j= (NList_long[i]- NList_long_[i]); j < NList_long[i];j++) List_long[i*MaxList_long + j] = List_long_[i*MaxList_long +j-(NList_long[i]- NList_long_[i])];
      }
    }
    
    for(i=0 ; i < Dim ;i++){
      if(NList[i] > MaxList){
	cerr<<"Number of neighbors of atom in List "<<i<<" exceeded maximum number of allowed neighbors in list! ("<<MaxList<<") ending session.\n";
	exit(0);
      }
      if(NList_long[i] > MaxList_long){
	cerr<<"Number of neighbors of atom in List "<<i<<" exceeded maximum number of allowed neighbors in list! ("<<MaxList_long<<") ending session.\n";
	exit(0);
      }
    }
    
    delete [] List_long_;

    delete [] NList_long_;
    delete [] List_;
    delete [] NList_;
  }
}

void Check_List(atom *particle, point *R0,int Natoms, int *List, int *NList, point L){
  int UpdatedList, i;
  double dx1, dy1, dz1, r1, R,R1, R2;
  R1=0;
  R2=0;
  
#pragma omp parallel for private(R,dx1,dy1,dz1)
  for(i=0; i<Natoms; i++){
    //for(UpdatedList=0, i=0; ((i<Natoms) && !UpdatedList) ; i++){
    dx1 = R_PBC_(particle[i].r[0],R0[i].x,L.x);
    dy1 = R_PBC_(particle[i].r[1],R0[i].y,L.y);
    dz1 = R_PBC_(particle[i].r[2],R0[i].z,L.z);

    R = sqrt(sqr(dx1) + sqr(dy1) + sqr(dz1));
    
    if(R > R1){
      R2=R1;
      R1=R;
    }
    
    if ((R1 + R2) > (Rcut_List_long-non_bond_cut)){
      UpdatedList = 1;
    }
  }
  
  if ((R1 + R2) > (Rcut_List_long-non_bond_cut)){
    UpdateNeighbList(Natoms,particle,List,NList,R0,L);
  }
  
}

double calc_dr_ij_drk(int atomi, int atomj, int atomk, int dir, double r_ij, atom *particle, int Natoms, point L)
{
  double Length;

  //r_ij  = R_ij(particle,atomi,atomj, L);
  
  if (dir == 0)    Length = L.x;
  else if (dir ==1)Length = L.y;
  else Length = L.z;
  
  if(atomk == atomi)      return( (particle[atomi].r[dir] - particle[atomj].r[dir] - (Length * rint( (particle[atomi].r[dir] - particle[atomj].r[dir]) / (Length + TINY)))) / r_ij);
  else if(atomk == atomj) return( (particle[atomj].r[dir] - particle[atomi].r[dir] - (Length * rint( (particle[atomj].r[dir] - particle[atomi].r[dir]) / (Length + TINY)))) / r_ij);
  else                    return(0.0);
}

double calc_dthetaijk_drn(int atomi, int atomj, int atomk, int atomn, int dir, atom *particle, double *angle, int *BondNum, int Natoms, double *rij, point L,double Thetaijk,double r_jk, double r_ij,dr_ij &rji,dr_ij &rjk,double cosijk)
{
  int indexji, indexjk;
  double ri, rj, rk, rjk_1, rji_1;
  double Length;
  
  indexji = atomj*Natoms + atomi;
  indexjk = atomj*Natoms + atomk;
  
  if((fabs(Thetaijk) > EPS) && (fabs(Thetaijk - PIE) > EPS)){
    if(atomn == atomi){
      
      rjk_1 = 1.0 /r_jk;
      rji_1 = 1.0 /r_ij;
      
      return( (1.0 / sqrt(1.0 - sqr(cosijk))) * rji_1 * (rji.r[dir] * rji_1 * cosijk - rjk.r[dir] * rjk_1) );
    }
    else if(atomn == atomk){
      
      rjk_1 = 1.0 / r_jk;
      rji_1 = 1.0 / r_ij;
      
      return( (1.0 / sqrt(1.0 - sqr(cosijk))) * rjk_1 * (rjk.r[dir] * rjk_1 * cosijk -  rji.r[dir] * rji_1) );
    }
    else if(atomn == atomj){
      
      rjk_1 = 1.0 / r_jk;
      rji_1 = 1.0 / r_ij;
      
      return( - (1.0 / sqrt(1.0 - sqr(cosijk))) * (rji_1 * (rji.r[dir] * rji_1 * cosijk - rjk.r[dir] * rjk_1) + rjk_1 * (rjk.r[dir] * rjk_1 * cosijk - rji.r[dir] * rji_1)) );
    }
    else return(0.0);
  }
  else return(0.0);
}



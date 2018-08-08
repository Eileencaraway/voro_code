#include <iostream>
#include <fstream>
#include "math.h"
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <vector>
#include <utility>
#include "voro++.hh"
//#include "src\voro++.cc"

using namespace std;
using namespace voro;
 
int const N = 1000; //number of particles
//double L = powl(N,1/3.0), HL=L/2.0;//box size and half of the box size
double L=9.56, HL=L/2.0;
double Rcut = 1.5, Rcut2 = Rcut*Rcut; //Rcut is determined from the first minimum of g(r)
// number of boxes the simulation volume is divided into, used by voro++; they claim is not relevant for the performance
int n_x, n_y, n_z;
double tol = 1E-8;
const double PI = 4*atan(1.0);

// class "particle" containg details on the particle and on the properties of its voronoi cell
class particle{
  public:
  double R[3]; // position [x,y,z]
  double q6;
  double Surface_areas;
  double Volume;
  double Nfaces;

  int num_neighbors;
  int neighbor[100];

  void add_neighbor(int m){ neighbor[num_neighbors] = m; num_neighbors++;}
  void reset_neighbors(){ num_neighbors = 0; }
  void set_VSN(double V, double S, int NF){ Volume = V; Surface_areas = S; Nfaces = NF; }
}Particles[N];

void determine_neighbors_using_cutoff(); 
void determine_neighbors_using_Voronoi();
double round(double r);
void ReadConf();//read the configuration file
void Eval_Q6_local();
void Eval_Q6_global(double &q6global);
void qlm(int &ql, int&qm, int &n, double &meanr, double &meani);
void outputdata(double &q6global);

int main() {
	 char buffer[100];
	 int i,j,step,s_witch, count, m;
     double q6global;
	 //===============
 //// //======read Trajectory
   ReadConf();
   //determine_neighbors_using_cutoff();
   // Create a pre-container class to import the input file and guess the best computational grid size to use.
   pre_container pcon(0, L, 0, L, 0, L,true,true,true);
   pcon.import("s_dumpfile.0100000000.txt");
   pcon.guess_optimal(n_x,n_y,n_z);
   //perform the Voronoi tessellation
   determine_neighbors_using_Voronoi();
   Eval_Q6_local();
   Eval_Q6_global(q6global);
   //out put information
   outputdata(q6global);
	return 0;
}

// it is better to use two vectors, for neigh and nnn, as the energy only involves neigh
void ReadConf(){
  int i;
    ifstream ConfFile;
	ConfFile.open("s_dumpfile.0100000000.txt");
	for(int n = 0; n < N; n++){
		ConfFile>>i>>Particles[n].R[0] >>Particles[n].R[1] >>Particles[n].R[2];
	}
	ConfFile.close();
}

void outputdata(double &q6global){
  int i;
   ofstream LOP,GOP, local_density;
   LOP.open("Q6_local.dat");
   GOP.open("Q6_global.dat");
   local_density.open("local_density.dat");
   for(int n = 0; n < N; n++){
	   LOP<<n<<"     "<<Particles[n].q6<<endl;
	   local_density<<n<<"     "<<1.0/Particles[n].Volume<<endl;
	}
   GOP<<"Q6_global:"<<endl;
   GOP<<q6global<<endl;
   LOP.close();
   GOP.close();
   local_density.close();
}

//determine the first shell neighbors
void determine_neighbors_using_cutoff(){
	double x,y,z,x1,y1,z1;
	double dr2,dx,dy,dz;
  for(int n = 0; n < N; n++){
	  Particles[n].reset_neighbors();
      x = Particles[n].R[0]; 
      y = Particles[n].R[1]; 
	  z = Particles[n].R[2]; 
	  for(int i = 0; i < N; i++){
		if(i!=n){
		  x1 = Particles[i].R[0]; 
          y1 = Particles[i].R[1]; 
	      z1 = Particles[i].R[2]; 
		  dx = x-x1;
		  dy = y-y1;
		  dz = z-z1;
		  //consider the periodic boundary
		  dx -=  L*round(dx/L);
		  dy -= L*round(dy/L);
		  dz -= L*round(dz/L);
		  dr2 = dx*dx + dy*dy + dz*dz;
		  if(dr2<Rcut2) Particles[n].add_neighbor(i);
		}
	  }
  } 
}

void determine_neighbors_using_Voronoi(){
  double x,y, z, V, S, NF;
  unsigned int i,j;
  int id;
  voronoicell_neighbor c;
  vector<int> neigh,f_vert;
  vector<double> f_area;
  vector<double> v;
  //// open files
  FILE *fp;
  //if (draw) fp=safe_fopen("tessellation.dat", "w");
  container conVoro(0, L,0, L,0, L,n_x,n_y,n_z,true,true,true,8);
  for(int n = 0; n < N; n++){
    x = Particles[n].R[0]; 
    y = Particles[n].R[1]; 
	z = Particles[n].R[2]; 
	Particles[n].reset_neighbors();
	conVoro.put(n,x,y,z);
  } 
  // Loop over all particles in the container and compute each Voronoi cell
  c_loop_all clVoro(conVoro);
  	if(clVoro.start()) do if(conVoro.compute_cell(c,clVoro)) {
		id=clVoro.pid();
		// Gather information about the computed Voronoi cell
		c.neighbors(neigh);
		c.face_areas(f_area);
		S = 0;
		for(i=0;i<neigh.size();i++) {
            S+=f_area[i];
		}
		//=========================
		V = c.volume();
	    NF = neigh.size();
		Particles[id].set_VSN(V,S,NF);
		//update neighbor information
		 Particles[id].num_neighbors = neigh.size();
		 for(i=0;i<neigh.size();i++) {
			Particles[id].neighbor[i] = neigh[i];
		}
	} while (clVoro.inc());
  conVoro.clear();
}

//====
double round(double r){
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

//function
int fact(int n){
    int  value;
	value = 1;
	for (int i =1; i<n+1;i++) value *=i;
	return value;
}

//function plgndr
double plgndr(int l, int m, double x){
	double pmm, somx2,fact,pmmp1,pll;
	//=================
	pmm = 1.0;	
	if(m>0){
		somx2 = sqrt((1.0 - x)*(1.0 + x));
	    fact = 1.0;
		for (int i=1;i<m+1;i++){
			pmm = -pmm*fact*somx2;
		    fact+=2.0;
		}
	}
	if(l==m) return pmm;
	else {
		pmmp1 = x*(2*m + 1)*pmm;
		if (l==(m+1)) return pmmp1;
		else{
			for(int ll=m+2; ll<l+1; ll++){
				pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm=pmmp1;
                pmmp1=pll;
			}
			return pll;
		}
	}
}

//===function qlm
void qlm(int &ql, int&qm, int &n, double &meanr, double &meani){
	int m;
   float const eps = 1e-6;
   double sumr = 0.0, sumi=0.0;
   double xmn,ymn,zmn,rijsq,rij,rxysq,rxy,costheta,phi;
   double constant, rpart, ipart, rylm, iylm;
   ofstream testfile;
   //========
   for (int i=0; i<Particles[n].num_neighbors;i++){
	   m = Particles[n].neighbor[i];
       xmn = Particles[m].R[0] - Particles[n].R[0];
	   ymn = Particles[m].R[1] - Particles[n].R[1];
	   zmn = Particles[m].R[2] - Particles[n].R[2];
	   //====periodic boundary condition
	   xmn -=  L*round(xmn/L);
	   ymn -=  L*round(ymn/L);
	   zmn -=  L*round(zmn/L);
	   rijsq = xmn*xmn +  ymn*ymn + zmn*zmn;
	   rij = sqrt(rijsq);
	   rxysq = xmn*xmn +  ymn*ymn;
	   rxy = sqrt(rxysq);
	   if (rxy<eps) rxy += eps;
	   costheta = zmn/rij;
	   if (ymn>=0.0) phi = acos(xmn/rxy);
	   else phi = 2.0*PI - acos(xmn/rxy);
	   if(qm>=0) {
           constant = sqrt((2*ql + 1)*fact(ql - qm)/4.0/PI/fact(ql+qm));
		   rpart = cos(qm*phi);
		   ipart = sin(qm*phi);
		   rylm = constant*plgndr(ql, qm, costheta)*rpart;
		   iylm = constant*plgndr(ql, qm, costheta)*ipart;
	   }
	   else{
		   constant = sqrt((2*ql+1)*fact(ql+qm)/4.0/PI/fact(ql-qm));
		   rpart = cos(-qm*phi)*powl(-1.0, qm);
		   ipart = sin(-qm*phi)*powl(-1.0, qm);
		   rylm = constant *plgndr(ql,-qm,costheta)*rpart;
		   iylm = constant *plgndr(ql,-qm,costheta)*ipart;
	   }
	   sumr += rylm;
	   sumi += iylm;
   }
   meanr = sumr/Particles[n].num_neighbors;
   meani = sumi/Particles[n].num_neighbors;
}

////===The local order parameters
void Eval_Q6_local(){
	int l;
	int mtr, m1,m2,m3,neij;
	double sum, meansq,sumrw,sumiw;
	double a, b, rp,ip;
	double a1,a2,a3,b1,b2,b3;
	double rqlocal[N][20], iqlocal[N][20], sumq[N];
	double rqvector[N][20], iqvector[N][20];
	l = 6;
	//ofstream q6file;
	//q6file.open("q6_local.dat");
	for(int i=0;i<N;i++){
		sum =0.0;
		for (int m=-l;m<l+1;m++){
			mtr = m + l;
		    qlm(l,m,i,a,b);
			rqlocal[i][mtr] = a;
			iqlocal[i][mtr] = b;
			meansq = a*a + b*b;
			sum += meansq;
		}
		sumq[i] = sqrt(sum);
		Particles[i].q6 = sqrt(4.0*PI*sum/(2*l + 1));
		//cout<<Particles[i].q6<<endl;
		//q6file<<Particles[i].RL[0]<<"    "<<Particles[i].RL[1]<<"     "<<Particles[i].RL[2]<<"   "<<Particles[i].q6<<endl;
	}
}

////===The global order parameters
void Eval_Q6_global(double &q6global){
	int l =6;
	int mtr,sumbond;
	double sumr,sumi, meansq,sumsq;
	double rq, iq;
	//================
	sumbond = 0;
	sumsq = 0.0;
	for(int i=0;i<N;i++){
		sumbond += Particles[i].num_neighbors;
	}
	//========
	for (int m=-6;m<7;m++){
		sumr = 0.0;
		sumi= 0.0;
		for(int i=0;i<N;i++){
		    qlm(l,m,i,rq,iq);
			//cout<<rq<<endl;
			sumr += rq*Particles[i].num_neighbors;
			sumi += iq*Particles[i].num_neighbors;
		}
		sumr/=sumbond;
		sumi/=sumbond;
		sumsq += (sumr*sumr + sumi*sumi);
	}
	q6global = sqrt(4.0*PI*sumsq/13.0);
}

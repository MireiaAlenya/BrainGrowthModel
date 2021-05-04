#include <omp.h>
#include "vema.h"
#include "eig3.h"
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <chrono>
#include <ctime>
using namespace std;

#define PATH_DIR "../Res.growth/"
#define MESH_FILE "../../Brains_LargerStudy/Meshes/F-s020-w26-820k-gt.mesh"
//#define LABEL_FILE "../../Meshes/fet-020-820k-labelfile.txt"
#define BETA 0.01
#define GROWTH_RELATIVE 1.829 
#define NAME "n.020"

#define ID 20
#define GA 26.4


class Tetra{
public:
    int n1, n2, n3, n4;
    Matrix G;
    Tetra(): n1(0), n2(0), n3(0), n4(0) {}
    Tetra(int n1_, int n2_, int n3_, int n4_): n1(n1_), n2(n2_), n3(n3_), n4(n4_) {}
};
class Face{
public:
    int n1, n2, n3;
    Face(): n1(0), n2(0), n3(0) {}
    Face(int n1_, int n2_, int n3_): n1(n1_), n2(n2_), n3(n3_) {}
};

void createNNLtriangle(vector<int>*, Vector*, vector<Face>&, int*, int, int, double, double, double);
Vector closestPointTriangle(Vector&, Vector&, Vector&, Vector&, double&, double&, double&);
void Eigensystem(Matrix A, Matrix& V, double d[3]);
void dist2surf(Vector*, int*, int, int, int*, double*);
void writeTXT(Vector*, vector<Face>&, Tetra*, int*, int, int, int, double);
void writeSTL(Vector*, vector<Face>&, int*, int*, int, int, int, int, double, int);
void writeMESH(Vector*, vector<Face>&, Tetra*, int*, int, int, int, int, int, double, int);
//void writeVTK(Vector*, vector<Face>&, Tetra*, Vector*, Vector*, Vector*, Vector*, Vector*, Vector*, Vector*, Vector*, double*, double*, double*, double*, double*, double*, Vector*, Vector*, int*, int, int, int, int, int, double, int);
void writeVTK(Vector*, vector<Face>&, Tetra*, Vector*, Vector*, Vector*, Vector*, Vector*, Vector*, Vector*, Vector*, double*, double*, double*, double*, double*, double*, Vector*, Vector*, int*, int, int, int, int, int, double, int);

int main(int argc, char* argv[]) {
	std::cout << "Try implement new relation between GA-time modifying growth tensor computation. Version 1."<< std::endl;
	
	std::cout<<"Initializing"<<std::endl;
	
	auto start = std::chrono::system_clock::now();
	auto start0 = std::chrono::system_clock::now();

	// ----------------------------------------------------------------------
    
    for (int i=0; i<argc; i++){
		cout<<"argv"<< i <<": " << argv[i] << endl;
	}
	// DEFINE CORTICAL THICKNESS
	int hname;
	double Hi; // Thickness of growing layer defined in the sh file.
	if (!strcmp(argv[1], "1")) { Hi = 0.0125; hname = 1;}  		//{ H = 0.74; }
	else if (!strcmp(argv[1], "2")) { Hi = 0.0251; hname = 2;}  	//{ H = 1.48; }
 	else if (!strcmp(argv[1], "3")) { Hi = 0.0336; hname = 3;} 	//{ H = 1.98; }
	else if (!strcmp(argv[1], "4")) { Hi = 0.0421; hname = 4;}  	//{ H = 2.48; }
	else if (!strcmp(argv[1], "5")) { Hi = 0.0505; hname = 5;} 	//{ H = 2.98; }
	else if (!strcmp(argv[1], "6")) { Hi = 0.1010; hname = 6;}  	//{ H = 5.96; }
	
	std::cout<< "Hi: "<<Hi <<std::endl;
    // ----------------------------------------------------------------------
    
    // CREATE FOLDER 
    int f = mkdir(PATH_DIR, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    char pd[50];
	sprintf(pd,"%s/%s", PATH_DIR, NAME);
    f = mkdir(pd, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    char inpa[10], inpb[10], inpc[10], inpd[10], inpe[10];
    
    // ----------------------------------------------------------------------

    // READ INPUT MESH
    ifstream filu;
    filu.open(MESH_FILE);
    std::cout<<MESH_FILE<<std::endl;
	// Read nodes
    filu >> inpa;
    int nn = atoi(inpa);
    Vector* Ut = new Vector[nn]();
    Vector* Ut0 = new Vector[nn]();
	Vector* Ut1 = new Vector[nn]();
    for (int i = 0; i < nn; i++) {
        filu >> inpa >> inpb >> inpc;
        Ut0[i] = Vector(atof(inpa), atof(inpb), atof(inpc));
        Ut[i] = Ut0[i];
		Ut1[i] = Ut0[i];
    }
    // Read elements
    filu >> inpa;
    int ne = atoi(inpa);
    Tetra* tets = new Tetra[ne]();
	Tetra* tets_list = new Tetra[ne]();
    for (int i = 0; i < ne; i++) {
        filu >> inpa >> inpb >> inpc >> inpd >> inpe;
        tets[i] = Tetra(atoi(inpb)-1, atoi(inpc)-1, atoi(inpe)-1, atoi(inpd)-1);
    }
    // Read faces
    filu >> inpa;
    int nf = atoi(inpa);
    vector<Face> faces;
    for (int i = 0; i < nf; i++) {
        filu >> inpa >> inpb >> inpc >> inpd;
        faces.push_back(Face(atoi(inpb)-1, atoi(inpc)-1, atoi(inpd)-1));
    }
    filu.close();
	
	// ----------------------------------------------------------------------
	/*
	// READ INPUT LABELS (PARCELLATION)
	ifstream filabel;
	filabel.open(LABEL_FILE);
	std::cout<<LABEL_FILE<<std::endl;
	
	//read values of parcellation lobes
	filabel >> inpa;
	int nlab = atoi(inpa);
	double*	drawem = new double[nlab]();
	for (int i = 0; i < nlab; i++) {
        filabel >> inpa;
        drawem[i] = atoi(inpa);
    }	
	filabel.close();
	*/
    // ----------------------------------------------------------------------

    // Find center of mass and dimension of the mesh
	double maxx = -1e9, minx = 1e9, maxy = -1e9, miny = 1e9, maxz = -1e9, minz = 1e9;
	Vector cog;
	for (int i = 0; i < nn; i++) {
		maxx = max(maxx, Ut[i].x); minx = min(minx, Ut[i].x);
		maxy = max(maxy, Ut[i].y); miny = min(miny, Ut[i].y);
		maxz = max(maxz, Ut[i].z); minz = min(minz, Ut[i].z);
		cog += Ut[i];
	}
	double lx = abs(maxx-minx);
	double ly = abs(maxy-miny);
	double lz = abs(maxz-minz);
	cog /= nn; // ***** center 
		
	double maxd = max(max(max(abs(maxx-cog.x), abs(minx-cog.x)), max(abs(maxy-cog.y), abs(miny-cog.y))), max(abs(maxz-cog.z), abs(minz-cog.z)));
	
	cout << "maxx: " << maxx << setw(11) << "minx: "<< minx << setw(11) << "maxy: " << maxy << setw(11) << "miny: " << miny << setw(11) << "maxz: " << maxz << setw(11) << "minz: " << minz << endl;
	cout << "lx: " << lx << setw(11) << "ly: "<< ly << setw(11) << "lz: " << lz << endl;
	cout << "cog.x: " << cog.x << setw(11) << "cog.y: " << cog.y << setw(11) << "cog.z: " << cog.z << setw(11) << "maxd: " << maxd << endl;
    
    // ----------------------------------------------------------------------

	// NORMALISE MESH: Change mesh information by values normalised 
    for (int i = 0; i < nn; i++){
        Ut0[i].x = -(Ut[i].y - cog.y)/maxd;  
        Ut0[i].y = (Ut[i].x - cog.x)/maxd;
        Ut0[i].z = -(Ut[i].z - cog.z)/maxd;
        Ut[i] = Ut0[i];
        Ut1[i] = Ut0[i];
    } 

	std::cout<< "Mesh normalisation done." <<std::endl;

	// Find new extremal values and center of mesh
	double nmaxx = -1e9, nminx = 1e9, nmaxy = -1e9, nminy = 1e9, nmaxz = -1e9, nminz = 1e9;
	Vector ncog;
	for (int i = 0; i < nn; i++) {
		nmaxx = max(nmaxx, Ut[i].x); nminx = min(nminx, Ut[i].x);
		nmaxy = max(nmaxy, Ut[i].y); nminy = min(nminy, Ut[i].y);
		nmaxz = max(nmaxz, Ut[i].z); nminz = min(nminz, Ut[i].z);
		ncog += Ut[i];
	}
	double nlx = abs(nmaxx-nminx);
	double nly = abs(nmaxy-nminy);
	double nlz = abs(nmaxz-nminz);
	ncog /= nn;                             // ***** center 
	// nmaxd = normalised max distance between the furthest point of the mesh and the center of coordinates	
	double nmaxd = max(max(max(abs(nmaxx-ncog.x), abs(nminx-ncog.x)), max(abs(nmaxy-ncog.y), abs(nminy-ncog.y))), max(abs(nmaxz-ncog.z), abs(nminz-ncog.z)));
	double mpy;    
	double at;
	
	// Midplane position
    mpy = (nminy+nmaxy);
	std::cout << "-----" << std::endl;
	cout << "nmaxx: " << nmaxx << setw(11) << "nminx: " << nminx << setw(11) << "nmaxy: " << nmaxy << setw(11) << "nminy: " << nminy << setw(11) << "nmaxz: " << nmaxz << setw(11) << "nminz: " << nminz << endl;
	cout << "nlx: " << nlx << setw(11) << "nly: " << nly << setw(11) << "nlz: " << nlz << endl;
	cout << "ncog.x: "<< ncog.x << setw(11) << "ncog.y: " << ncog.y << setw(11) << "ncog.z: " << ncog.z << setw(11) << "nmaxd: " << nmaxd << endl;
	std::cout << "-----" << std::endl;
	cout << "mpy: " << mpy << endl; 

    // ----------------------------------------------------------------------

	// MESH SPACING DEFINITION. 
	// Check minimum and maximum edge lengths at the surface
	double nmine = 1e9, nmaxe = 0.0, nave = 0.0;
	for (int i = 0; i < nf; i++) {
		nmine = min((Ut[faces[i].n2] - Ut[faces[i].n1]).length(), nmine);
		nmine = min((Ut[faces[i].n3] - Ut[faces[i].n1]).length(), nmine);
		nmine = min((Ut[faces[i].n3] - Ut[faces[i].n2]).length(), nmine);
		nmaxe = max((Ut[faces[i].n2] - Ut[faces[i].n1]).length(), nmaxe);
		nmaxe = max((Ut[faces[i].n3] - Ut[faces[i].n1]).length(), nmaxe);
		nmaxe = max((Ut[faces[i].n3] - Ut[faces[i].n2]).length(), nmaxe);
		nave += (Ut[faces[i].n3] - Ut[faces[i].n2]).length() + (Ut[faces[i].n3] - Ut[faces[i].n1]).length() + (Ut[faces[i].n2] - Ut[faces[i].n1]).length();
	}
	nave /= 3.0*nf; //***** average value of edge length
    std::cout<<"Normalised edge lengths at the surface: "<<std::endl;
	double a = (nave + nmine)/2;
	cout << "mine: " << nmine  << setw(11) << "ave: " << nave << setw(11) << "maxe: " << nmaxe << setw(11) << "a: " << a << endl;
	
    // ----------------------------------------------------------------------

    // Determine surface nodes and index maps (SN = global indices of the surface nodes, SNb = surface indices of global nodes)
    double gr_average = 0.0;
    int nsn = 0;                                                    // Number of nodes at the surface
    int SNb[nn];                                                    // Nodal index map from full mesh to surface 
    for (int i = 0; i < nn; i++) SNb[i] = 0;                        // Initialization SNb with all 0s
    for (int i = 0; i < nf; i++) { SNb[faces[i].n1] = 1; SNb[faces[i].n2] = 1; SNb[faces[i].n3] = 1; }      //*** nodes in surface triangles 
    for (int i = 0; i < nn; i++) if (SNb[i] == 1) nsn++;            //***** nsn = 50943
    int SN[nsn];                                                    // Nodal index map from surface to full mesh
    int p = 0;                                                      // Iterator
    for (int i = 0; i < nn; i++) if (SNb[i] == 1) { SN[p] = i; SNb[i] = p; p++; }
    // **** p = 50943, if the node is surfacal, SNb[nodeIndex] is the index of its surface triangle
	// **** SN saves index of surface triangles 

    // ----------------------------------------------------------------------

	std::cout<< "ne: "<< ne << std::endl;
	std::cout<< "nf: " << nf << std::endl;
	std::cout<< "nn: " << nn <<std::endl;
	std::cout<< "nsn: " << nsn <<std::endl;
	
	const int id = ID;
    int spw = 0;                            // output counter for stl paraview
	double H = Hi;                    
    
    const double mu = 1.0;                  // Shear modulus
    const double mug = 1.0;                       // Shear modulus of gray matter
	const double muw = 1.167;                     // Shear modulus of white matter
    const double K = 5.0;                 	// Bulk modulus
    const double rho = 0.005;              	// Mass density
    const double gamma = 0.5;              	// Damping constant
	const double kc = 10.0*K; 	            // Contact stiffness
	const double fn_ct = 2.0; 	            // Contact stiffness
	
	const double hs = 0.6*a;            	// Thickness of contact proximity skin
    const double hc = 0.2*a;            	// Thickness of repulsive skin
    const double bw = 3.2; 		            // Width of a bounding box, centered at origin, that encloses the whole geometry even after growth 
	const double mw = 8*a; 		            // Width of a cell in the linked cell algorithm for proximity detection
	
    double ga = GA; 						// relative ga to t
	double lastga = 0;                      // used to print output each value of integer GA
    int step = 0;							//initial step
    double t = 0.0;                    		// Current time 
    //double tt = 0.00006926*(ga*ga*ga) - 0.00665*(ga*ga) + 0.250*ga - 3.0189;  // 27.3 GA --> t0 = 0.25911434142
	double t0 = 0.987*exp(-exp(-0.134*(ga-29.433)));			// X. W. Approach time-GA, based on GI (instead of volume).
	t = t0;                          		//current time (adapted to each mesh)
    std::cout << "Time starting from t0=t= " << t << std::endl;
	double GADIF = GA - 22.0;
    double gaprint; 
    const int di = 500;               		// Output data once every di steps
	const int dicout = 500;                // Print values once every dicout steps
	const int dimesh = 1500; 		   		// mesh output once every dimesh steps  
    const double dt = 0.05*sqrt(rho*a*a/K); // Time step
	
	std::cout << "a: " << a << std::endl;
    
	std::cout << "mug: " << mug << std::endl;
	std::cout << "muw: " << muw << std::endl;
	std::cout << "K: " << K << std::endl;
	std::cout << "rho: " << rho << std::endl;
	std::cout << "gamma: " << gamma << std::endl;	
	std::cout << "kc: " << kc << std::endl;	
	std::cout << "fn_ct: " << fn_ct << std::endl;
	
	std::cout << "dt: " << dt << std::endl;
	std::cout << "di: " << di << std::endl;
	std::cout << "dicout: " << dicout << std::endl;
	std::cout << "dimesh: " << dimesh<< std::endl;
   
	Vector* Vt = new Vector[nn]();      	// Velocities
    Vector* Ft = new Vector[nn]();      	// Forces
	
	int countgm = 0;
    int countwm = 0;
	
    double* Vn0 = new double[nn];       	// Nodal volumes in reference state
    double* Vn = new double[nn];        	// Nodal volumes in deformed state
    double* Vn_init = new double[nn];       // Deformed nodal volumes
    Matrix I;                           	// Unit matrix
    double Ue, Ue_xw, Area;                 // Elastic energy, surface area      
    std::vector<int> NNLt[nsn];         	// Triangle-proximity lists for surface nodes  
    double maxDist;                     	// Maximum displacement since the last update of the proximity list
    Vector* Utold = new Vector[nsn];    	// Stores positions when proximity list is updated
    double ub, vb, wb;                  	// Barycentric coordinates
    int tp;                             	// Loop iterator

    std::cout << "bound box: " << bw << setw(11) << "  Width of a cell: " << mw << std::endl;
	std::cout << "hc: " << hc << setw(11) << "  hs: " << hs << setw(11) << "  kc: " << kc <<std::endl;
    double totalsteps = (1-t) / dt;
	std::cout << "totalsteps: " << totalsteps << std::endl;

    // ----------------------------------------------------------------------
	
	double* mu_vec = new double[ne];     	// vector mu
	double* gm_vec = new double[ne];     	// vector gm
	double* gr =new double[nn];
	double* gr_vtk = new double[nn];
	double* contacts = new double[nn];
	double* contacts_nsn = new double[nsn];
	//double* draw = new double[nn]; 
	double* Hilab = new double[nn]; 
	//double* Hilabne = new double[ne]; 
	double* Hlab = new double[nn]; 
	//double* Hlabne = new double[ne]; 
	
	// Mark non-growing areas
	for (int i = 0; i < nn; i++) {
		Vector qp = Ut0[i];
		// x: 0.05 = offset center of the ellipsoid. 0.714 defines its long-axis dimension. Reduction -> more skewed ellipsoid. 
		// z: -0.05 offset of vertical axis. If +0.25, ellipsoid in the top part of the cortex. 
		double rqp = Vector((qp.x+0.01)*0.714, qp.y, qp.z-0.05).length();
		if ( rqp < 0.6 ) {
            gr[i] = max(1.0 - 10.0*(0.6-rqp), 0.0);	            // Ellipsoid for white matter
            countwm++;
			if (gr[i] == 0.0) {gr_vtk[i] = 0.0;} 
			else {gr_vtk[i] = 0.5;}
		}
		else {                                                  // Grey matter indicator
            gr[i] = 1.0;
            countgm++;
			gr_vtk[i] = 1.0; 
		}
	} 
	std::cout << "-------------------" << std::endl;
	cout << "Number: " << 1.0 << setw(11) << "Grey matter nodes: " << countgm << "\n";
    cout << "Number: " << 0.0 << setw(11) << "White matter nodes: " << countwm << "\n";

	/* ----------------------------------------------------------------------
	int count_parietal = 0;
	int count_occipital = 0; 
	int count_other = 0; 
	double val_parietal = 1.0;
	double val_occipital = 1.0;
	double val_other = 1.0;
    // LABEL NODES according to input metric file
	for (int i = 0; i < nn; i++) {
		Vector qp = Ut0[i];
		if ( drawem[i] == 38 || drawem[i] == 39 ) {
            draw[i] = val_parietal;	            // larger if node is of parietal zone
			Hilab[i] = Hi*val_parietal;
            count_parietal++;
		}
		else if ( drawem[i] == 22 || drawem[i] == 23 ) {
			draw[i] = val_occipital;
			Hilab[i] = Hi*val_occipital;
			count_occipital++;
		}
		else {                                                 		// Grey matter indicator
            draw[i] = val_other;
            Hilab[i] = Hi*val_other;
            count_other++;

		}
	} 

	std::cout << "-------------------" << std::endl;
	cout << "Number: " << val_parietal << setw(11) << " Parietal nodes: " << count_parietal << "\n";
	cout << "Number: " << val_occipital << setw(11) << " Occipital nodes: " << count_occipital << "\n";
    cout << "Number: " << val_other << setw(11) << " Other region nodes: " << count_other << "\n";
	std::cout << "-------------------" << std::endl;
	cout << "Hi par: " << Hi*val_parietal << setw(11) << " Parietal nodes: " << count_parietal << "\n";
	cout << "Hi occ: " << Hi*val_occipital << setw(11) << " Occipital nodes: " << count_occipital << "\n";
    cout << "Hi oth: " << Hi*val_other << setw(11) << " Other region nodes: " << count_other << "\n";
	std::cout << "-------------------" << std::endl;
	
	 ----------------------------------------------------------------------*/
	
    // COMPUTE VOLUME
    for (int i = 0; i < nn; i++) {  Vn_init[i] = 0.0; }
	
	for (int i = 0; i < ne; i++) {
		int n1 = tets[i].n1;
		int n2 = tets[i].n2;
		int n3 = tets[i].n3;
		int n4 = tets[i].n4;	

		// Deformed
		Vector x1_init = Ut[n2] - Ut[n1];
		Vector x2_init = Ut[n3] - Ut[n1];
		Vector x3_init = Ut[n4] - Ut[n1];
		Matrix A_init = Matrix(x1_init, x2_init, x3_init);
		double vol_init = A_init.det()/6.0;
		Vn_init[n1] += vol_init/4.0;
		Vn_init[n2] += vol_init/4.0;
		Vn_init[n3] += vol_init/4.0;
		Vn_init[n4] += vol_init/4.0;	
	}
	
	double Volume_init = 0.0; // Volume
	for (int i = 0; i < nn; i++) Volume_init += Vn_init[i];
	
	std::cout << "Initial volume is: " << Volume_init << std::endl;
	std::cout << "Initial volume is: " << -Volume_init*maxd*maxd*maxd << std::endl;
    
	// ----------------------------------------------------------------------

    // Find closest surface nodes (csn) and distances to them (d2s)
    // Determine nearest surface nodes to nodes, distances to surface, and surface normals - these are needed to set up the growth of the gray matter
    int csn[nn];                            	// Closest surface nodes
    double* d2s = new double[nn];           	// Distances between nodes and their respective closest surface nodes
    dist2surf(Ut, SN, nn, nsn, csn, d2s);   	// Finds the nearest surface nodes (csn) to nodes and distances to them (d2s)
    Vector* N0 = new Vector[nsn];           	// Normals in reference state
    Vector Ntmp;        
	Vector* mat_G1 = new Vector[ne];			// Growth tensor values for print in VTK file
	Vector* mat_G2 = new Vector[ne]; 
	Vector* mat_G3 = new Vector[ne];
	Vector* mat_G_diag = new Vector[ne];
    
	Vector* mat_T1 = new Vector[ne];			// Stress tensor values for print in VTK file
	Vector* mat_T2 = new Vector[ne]; 
	Vector* mat_T3 = new Vector[ne];
	Vector* mat_T_diag = new Vector[ne];
	
    // RECALCULATE NORMALS 
	for (int i = 0; i < nsn; i++) {
		   N0[i].x = 0.0;
		   N0[i].y = 0.0;
		   N0[i].z = 0.0;     
	}
	
    for (int i = 0; i < nf; i++) {          	// Find normals
        Ntmp = (Ut0[faces[i].n2] - Ut0[faces[i].n1]).cross(Ut0[faces[i].n3] - Ut0[faces[i].n1]);
        N0[SNb[faces[i].n1]] += Ntmp;
        N0[SNb[faces[i].n2]] += Ntmp;
        N0[SNb[faces[i].n3]] += Ntmp;
    }

    for (int i = 0; i < nsn; i++) N0[i].normalize();  // *** normalize divided by sqrt(x*x + y*y + z*z)

    // ----------------------------------------------------------------------

	double outputS = 0, tmdistS = 0, tpenS = 0, tnodS = 0, tdefS = 0, tmpyS = 0, tprintS = 0, tnewtonS = 0;
    const double eps = 0.1;
	const double k = 0.0;
	
    cout.precision(3);
    ofstream dat;
    ofstream datt;
    ofstream dats1;
    ofstream dats2;
    ofstream dats3;
    ofstream datp;
    ofstream datgatime;
    dat.precision(5);
    datt.precision(5);
    dats1.precision(5);
    dats2.precision(5);
    dats3.precision(5);
    datp.precision(5);
    datgatime.precision(5);
	char datfn[100];
	char datfnt[100];
	char datfns1[100];
	char datfns2[100];
	char datfns3[100];
	char datfnp[100];
	char datfngatime[100];
	sprintf(datfn,"%s/%s/%s_h%d_data.dat", PATH_DIR, NAME, NAME, hname);
	sprintf(datfnt,"%s/%s/%s_h%d_timing.dat", PATH_DIR, NAME, NAME, hname);
	sprintf(datfns1,"%s/%s/%s_h%d_step1.dat", PATH_DIR, NAME, NAME, hname);
	sprintf(datfns2,"%s/%s/%s_h%d_step2.dat", PATH_DIR, NAME, NAME, hname);
	sprintf(datfns3,"%s/%s/%s_h%d_step3.dat", PATH_DIR, NAME, NAME, hname);
	sprintf(datfnp,"%s/%s/%s_h%d_print.dat", PATH_DIR, NAME, NAME, hname);
	sprintf(datfngatime,"%s/%s/%s_h%d_gatime.dat", PATH_DIR, NAME, NAME, hname);
	datt.open(datfnt);
	datgatime.open(datfngatime);
	dats1.open(datfns1);
	dats2.open(datfns2);
	dats3.open(datfns3);
	datp.open(datfnp);
	dat.open(datfn);
	
	dat << ID << endl;
	dat << GA << endl; 
	dat << GADIF << endl;
	dat << MESH_FILE << endl;
	dat << Hi << endl;
	dat << GROWTH_RELATIVE << endl;
	dat << "LINEAL" << endl;
	
    std::cout << "----- end preparation, start simulation loop -----" << std::endl;
	datgatime << "LINEAL" << endl;
	datgatime << "ga: " << setw(16) << "time: " << setw(16) << "step: " << setw(16) << "at: " << endl;
    // ----------------------------------------------------------------------
   
	int cpen = 0; 				// counter for penalty force
	int cmdist = 0; 			// counter for cmdist calculation (updating lists)

	
    //while (t < 1.1) { // Main loop
    while (ga < 46.1) { // Main loop
       	
		//LINEAL
        at = GROWTH_RELATIVE*((t-t0)/(1-t0)); 		// relative growth 
		
		// GOMPERTZ 1st:
		//at = GROWTH_RELATIVE*exp(-exp(-6.6*(t-0.43)));
		// GOMPERTZ 2nd:
		//at = GROWTH_RELATIVE*exp(-exp(-7.5*(t-0.19)));
		// LOGISTIC MODEL:
		//at = GROWTH_RELATIVE/(1+exp(-50*(t-0.9)))
		
		H = Hi + BETA*t;
		//for (int i = 0; i < nn; i ++) {
			//Hlab[i] = Hilab[i] + BETA*t; 
		//}
		if (step%di == 0) {datgatime << ga << setw(16) << t << setw(16) << step << setw(16) << at << endl; }
        // ----------------------------------------------------------------------
        // just to print values with GA <-> t consistency.
        //ga = 33.82802*(pow(t,5)) - 75.50047*(pow(t,4)) + 40.05540*(pow(t,3)) + 6.26943*(pow(t,2)) + 17.29106*t + 22.00818;      //t=1 --> 43.94, t = 1.05 --> 44.8478
		// time is what is in every equation. ga is just for printing issues!!! This ga(t) here, does not have any effect on constitutive or model equations.
		ga = -(log(log(0.987/t))/0.134)+29.433; 
        // ----------------------------------------------------------------------
        
		// CONTACT PROCESSING   
		auto count_0 = std::chrono::system_clock::now();
		
		maxDist = 0.0;
        #pragma omp parallel for reduction(max:maxDist)
		for (int i = 0; i < nsn; i++ ) { 
			// Find the maximum displacement since the last update of the proximity list
			maxDist = max(maxDist, (Ut[SN[i]] - Utold[i]).length());
		}		
		if (maxDist > 0.5*(hs-hc)) { 
			// Generates point-triangle proximity lists (NNLt[nsn]) using the linked cell algorithm
			createNNLtriangle(NNLt, Ut, faces, SN, nsn, nf, hs, bw, mw); 
			
			// Update proximity list
			for (int i = 0; i < nsn; i++) Utold[i] = Ut[SN[i]];
			cmdist++;
		}   
		
		for (int i = 0; i < nsn; i++) { contacts_nsn[i] = 0; }
        
		auto count_1 = std::chrono::system_clock::now();
        
		#pragma omp parallel for private(tp)
        for (int i = 0; i < nsn; i++) { 		                        // Determine the closest triangles to points
            for (tp = 0; tp < NNLt[i].size(); tp++) {
                int pt = SN[i];
                int tri = NNLt[i][tp];
                Vector cc = closestPointTriangle(Ut[pt], Ut[faces[tri].n1], Ut[faces[tri].n2], Ut[faces[tri].n3], ub, vb, wb) - Ut[pt];
                double rc = cc.length();
                if (rc < hc && gr[pt] + gr[faces[tri].n1] > 0.0) {      // Calculate contact force if within the contact range
                    cc.normalize();
					Vector Ntri = (Ut[faces[tri].n2] - Ut[faces[tri].n1]).cross(Ut[faces[tri].n3] - Ut[faces[tri].n1]); // **** triangle normal 
					Ntri.normalize();
					Vector fn = cc*(rc-hc)/hc*kc*a*a;                   // kc = 10.0*K Contact stiffness
					if (fn.dot(Ntri) < 0.0) fn -= Ntri*fn.dot(Ntri)*fn_ct;
					Ft[faces[tri].n1] -= fn*ub;
					Ft[faces[tri].n2] -= fn*vb;
					Ft[faces[tri].n3] -= fn*wb;
					Ft[pt] += fn;
                    contacts_nsn[i] = contacts_nsn[i] + 1;	
					cpen++;
                }
            }
        }
		
		dats1 << step << endl;
		auto count_2 = std::chrono::system_clock::now();
		
		for (int i = 0; i < nn; i++) { 
			contacts[i] = contacts_nsn[SNb[i]]; 
		}
		
		// ----------------------------------------------------------------------

        // Nodal volumes for nodal pressure calculation
        #pragma omp parallel for
        for (int i = 0; i < nn; i++) { Vn0[i] = 0.0; Vn[i] = 0.0; }
        
		// Calculation of Vn0[i] and Vn[i] for each tetrahedron
        #pragma omp parallel for
        for (int i = 0; i < ne; i++) {
            int n1 = tets[i].n1;
            int n2 = tets[i].n2;
            int n3 = tets[i].n3;
            int n4 = tets[i].n4;
            
            // Undeformed
            Vector xr1 = Ut0[n2] - Ut0[n1];
            Vector xr2 = Ut0[n3] - Ut0[n1];
            Vector xr3 = Ut0[n4] - Ut0[n1];
            Matrix Ar = Matrix(xr1, xr2, xr3);  // Reference state
            Ar = tets[i].G.prod(Ar);
            
            double vol0 = Ar.det()/6.0;         // Reference volume
            Vn0[n1] += vol0/4.0;
            Vn0[n2] += vol0/4.0;
            Vn0[n3] += vol0/4.0;
            Vn0[n4] += vol0/4.0;
            
            // Deformed
            Vector x1 = Ut[n2] - Ut[n1];
            Vector x2 = Ut[n3] - Ut[n1];
            Vector x3 = Ut[n4] - Ut[n1];
            Matrix A = Matrix(x1, x2, x3);      // Deformed state
            double vol = A.det()/6.0;           // Deformed volume
            Vn[n1] += vol/4.0;
            Vn[n2] += vol/4.0;
            Vn[n3] += vol/4.0;
            Vn[n4] += vol/4.0;
        }
		auto count_3 = std::chrono::system_clock::now();

        // ----------------------------------------------------------------------

        // DEFORMATIONS
        Ue = 0.0;   
		int ninv = 0;
        #pragma omp parallel for reduction(+:Ue, ninv)
		for (int i = 0; i < ne; i++) {// ***** loop tets
			// Nodal indices
			int n1 = tets[i].n1;
			int n2 = tets[i].n2;
			int n3 = tets[i].n3;
			int n4 = tets[i].n4;
			
			// Gray and white matter indicators
			//double gr_average = 0.25*(gr[n1]+gr[n2]+gr[n3]+gr[n4]);
			//H = 0.25*(Hlab[n1]+Hlab[n2]+Hlab[n3]+Hlab[n4]);
			double gm = 1.0/(1.0 + exp(10.0*(0.25*(d2s[n1]+d2s[n2]+d2s[n3]+d2s[n4])/H - 1.0))) * 0.25*(gr[n1]+gr[n2]+gr[n3]+gr[n4]); //***0.25 is divided by 4, get mean
			//gm = gm * 0.25*(draw[n1] + draw[n2] + draw[n3] + draw[n4]); // * regional growth approximation
			double wm = 1.0 - gm;
			double mu = muw*wm + mug*gm; // modulus of white matter and gray matter
			
			// Basis vector of reference state
			Vector xr1 = Ut0[n2] - Ut0[n1];
			Vector xr2 = Ut0[n3] - Ut0[n1];
			Vector xr3 = Ut0[n4] - Ut0[n1];
			Matrix Ar = Matrix(xr1, xr2, xr3); // Reference state
			Ar = tets[i].G.prod(Ar); // Apply growth to reference state
			
			// Deformed basis vectors
			Vector x1 = Ut[n2] - Ut[n1];
			Vector x2 = Ut[n3] - Ut[n1];
			Vector x3 = Ut[n4] - Ut[n1];

			// Undeformed normals
			xr1 = Vector(Ar.a, Ar.d, Ar.g);
			xr2 = Vector(Ar.b, Ar.e, Ar.h);
			xr3 = Vector(Ar.c, Ar.f, Ar.i);
			Vector N1 = xr3.cross(xr1);
			Vector N2 = xr2.cross(xr3);
			Vector N3 = xr1.cross(xr2);
			Vector N4 = (xr2 - xr3).cross(xr1 - xr3);

			Matrix A = Matrix(x1, x2, x3); // Deformed state
			double vol = A.det()/6.0; // Deformed volume
			Matrix F = A.prod(Ar.inv()); // Deformation gradient
			Matrix B = F.prod(F.trans()); // Left Cauchy-Green strain tensor //***** trans() --- transpose 
			double J = F.det(); // Relative volume change
			
			int flagJ = 0; 
			if (J < 0 || flagJ == 1) { 
				if (flagJ == 0) { datp << "-.-.-.-.- NEGATIVE JACOBIAN HERE, PRINTING INFO: -.-.-.-.-.-" << endl; }
				datp << "step = " << step << endl; 
				datp << "ne = " << i << endl; 
				datp << "J = " << J << endl; 
				flagJ = 1; 
			}
			
			double J1 = Vn[n1]/Vn0[n1];
			double J2 = Vn[n2]/Vn0[n2];
			double J3 = Vn[n3]/Vn0[n3];
			double J4 = Vn[n4]/Vn0[n4];
			double Ja = (J1 + J2 + J3 + J4)/4.0; // Averaged nodal volume change
						
			double powJ23, W;
			Matrix P;
			Matrix S;
			if (B.EV().z >= eps*eps && J > 0.0) { // No need for SVD
			
				//powJ23 = 1.0 + 2.0/3.0*(J - 1.0) - 1.0/9.0*(J-1.0)*(J-1.0); // Approximate pow(J, 2/3)
				powJ23 = pow(J, 2.0/3.0); // **** J^2/3
				Matrix T = (B - I*B.trace()/3.0)*mu/(J*powJ23) + I*K*(Ja-1.0); // **** 
				S = T; 
 				P = T.prod(F.trans().inv())*J; // 
				W = 0.5*mu*(B.trace()/powJ23 - 3.0) + 0.5*K*( (J1-1.0)*(J1-1.0) + (J2-1.0)*(J2-1.0) + (J3-1.0)*(J3-1.0) + (J4-1.0)*(J4-1.0) )*0.25;
			
			} else { // Needs SVD
				Matrix C = F.trans().prod(F);
				Matrix V;
				double eva[3];
				Eigensystem(C, V, eva);
			
				double l1 = sqrt(eva[0]);
				double l2 = sqrt(eva[1]);
				double l3 = sqrt(eva[2]);
				
				if (V.det() < 0.0) { V.a = -V.a; V.d = -V.d; V.g = -V.g; } // loop1++;}
			
				Matrix Fdi;
				if (l1 >= 1e-25) {Fdi.a = 1.0/l1; } // vlowl1++; }
				Fdi.e = 1.0/l2;
				Fdi.i = 1.0/l3;
			
				Matrix U = F.prod(V.prod(Fdi));
			
				if (l1 < 1e-25) {
					//loop2++;
					U.a = U.e*U.i - U.h*U.f;
					U.d = U.h*U.c - U.b*U.i;
					U.g = U.b*U.f - U.e*U.c;
				}
				if (F.det() < 0.0) {
					ninv++;
					l1 = -l1;
					U.a = -U.a; U.d = -U.d; U.g = -U.g;
				}
				
				Matrix Pd;
				double pow23 = pow(eps*l2*l3, 2.0/3.0);
				Pd.a = mu/3.0*(2.0*eps - l2*l2/eps - l3*l3/eps)/pow23 + k*(l1-eps) + K*(Ja-1.0)*l2*l3;
				Pd.e = mu/3.0*(-eps*eps/l2 + 2.0*l2 - l3*l3/l2)/pow23 + mu/9.0*(-4.0*eps/l2 - 4.0/eps*l2 + 2.0/eps/l2*l3*l3)/pow23*(l1-eps) + K*(Ja-1.0)*l1*l3;
				Pd.i = mu/3.0*(-eps*eps/l3 - l2*l2/l3 + 2.0*l3)/pow23 + mu/9.0*(-4.0*eps/l3 + 2.0/eps*l2*l2/l3 - 4.0/eps*l3)/pow23*(l1-eps) + K*(Ja-1.0)*l1*l2;
				W = 0.5*mu*((eps*eps + l2*l2 + l3*l3)/pow23 - 3.0) + mu/3.0*(2.0*eps - l2*l2/eps - l3*l3/eps)/pow23*(l1-eps) + 0.5*k*(l1-eps)*(l1-eps) + 0.5*K*((J1-1.0)*(J1-1.0) + (J2-1.0)*(J2-1.0) + (J3-1.0)*(J3-1.0) + (J4-1.0)*(J4-1.0))/4.0;
				
				P = U.prod(Pd.prod(V.trans()));
				
				//S = P.prod(F.trans())*J.inv();
			}
				
			if (J*J > 1e-50) Ue += W*vol/J; // Increment total elastic energy
			
			// Apply forces to nodes
			Ft[n1] += P.prod(N1 + N2 + N3)/6.0;
			Ft[n2] += P.prod(N1 + N3 + N4)/6.0;
			Ft[n3] += P.prod(N2 + N3 + N4)/6.0;
			Ft[n4] += P.prod(N1 + N2 + N4)/6.0;	
				
			// Growth
			Vector Ns = (N0[csn[n1]] + N0[csn[n2]] + N0[csn[n3]] + N0[csn[n4]]); // Surface normal
			Ns.normalize();
			
			tets[i].G = I + (I - Matrix(Ns.x*Ns.x, Ns.x*Ns.y, Ns.x*Ns.z, Ns.x*Ns.y, Ns.y*Ns.y, Ns.y*Ns.z, Ns.x*Ns.z, Ns.y*Ns.z, Ns.z*Ns.z))*gm*at;
		
			mat_G1[i] = Vector(tets[i].G.a, tets[i].G.b, tets[i].G.c);
			mat_G2[i] = Vector(tets[i].G.d, tets[i].G.e, tets[i].G.f);
			mat_G3[i] = Vector(tets[i].G.g, tets[i].G.h, tets[i].G.i);
			mat_G_diag[i] = Vector(tets[i].G.a, tets[i].G.e, tets[i].G.i);
			
			mat_T1[i] = Vector(P.a, P.b, P.c);
			mat_T2[i] = Vector(P.d, P.e, P.f);
			mat_T3[i] = Vector(P.g, P.h, P.i);
			mat_T_diag[i] = Vector(P.a, P.e, P.i);
			
			/*
			if ( i == 59893 || i == 13873 ) {
			//if (i > 10 && i < 16) {
				if (step > 1280 && step < 1310) {
			
					dat11 << "---------------------------------" << std::endl;
					dat11 << "step: " << step << std::endl;
					dat11 << "element: " << i << std::endl;
					dat11 << "gm: " << gm << std::endl;
					//dat11 << "elementSN: " << in << std::endl;
					
					dat11 << "---------" << std::endl;
					dat11 << "n1: " << n1 << std::endl;
					dat11 << "n2: " << n2 << std::endl;
					dat11 << "n3: " << n3 << std::endl;
					dat11 << "n4: " << n4 << std::endl;
					dat11 << "---------" << std::endl;
					dat11 << "x1: " << x1.x << " " << x1.y << " " << x1.z << std::endl;
					dat11 << "x2: " << x2.x << " " << x2.y << " " << x2.z << std::endl;
					dat11 << "x3: " << x3.x << " " << x3.y << " " << x3.z << std::endl;
					dat11 << "---------" << std::endl;
					dat11 << "xr1: " << xr1.x << " " << xr1.y << " " << xr1.z << std::endl;
					dat11 << "xr2: " << xr2.x << " " << xr2.y << " " << xr2.z << std::endl;
					dat11 << "xr3: " << xr3.x << " " << xr3.y << " " << xr3.z << std::endl;
					dat11 << "---------" << std::endl;					
					

					dat11 << "A: " << A.a << " " << A.b << " " << A.c << std::endl;
					dat11 << "A: " << A.d << " " << A.e << " " << A.f << std::endl;
					dat11 << "A: " << A.g << " " << A.h << " " << A.i << std::endl;
					
					dat11 << "---------" << std::endl;
					dat11 << "Ar: " << Ar.a << " " << Ar.b << " " << Ar.c << std::endl;
					dat11 << "Ar: " << Ar.d << " " << Ar.e << " " << Ar.f << std::endl;
					dat11 << "Ar: " << Ar.g << " " << Ar.h << " " << Ar.i << std::endl;
					dat11 << "---------" << std::endl;					
					dat11 << "Ftn1: " << Ft[n1].x << " " << Ft[n1].y << " " << Ft[n1].z << std::endl;
					dat11 << "Ftn2: " << Ft[n2].x << " " << Ft[n2].y << " " << Ft[n2].z << std::endl;
					dat11 << "Ftn3: " << Ft[n3].x << " " << Ft[n3].y << " " << Ft[n3].z << std::endl;
					dat11 << "Ftn4: " << Ft[n4].x << " " << Ft[n4].y << " " << Ft[n4].z << std::endl;
					
					dat11 << "---------" << std::endl;
					
					dat11 << "N1: " << N1.x << " " << N1.y << " " << N1.z << std::endl;
					dat11 << "N2: " << N2.x << " " << N2.y << " " << N2.z << std::endl;
					dat11 << "N3: " << N3.x << " " << N3.y << " " << N3.z << std::endl;
					dat11 << "N4: " << N4.x << " " << N4.y << " " << N4.z << std::endl;
					dat11 << "---------" << std::endl;
					dat11 << "Ns: " << Ns.x << " " << Ns.y << " " << Ns.z << std::endl;
					
					
					dat11 << "---------" << std::endl;
					
					dat11 << "DefGradient: " << F.a << " " << F.b << " " << F.c << std::endl;
					dat11 << "DefGradient: " << F.d << " " << F.e << " " << F.f << std::endl;
					dat11 << "DefGradient: " << F.g << " " << F.h << " " << F.i << std::endl;
					dat11 << "---------" << std::endl;
					dat11 << "LeftCG: " << B.a << " " << B.b << " " << B.c << std::endl;
					dat11 << "LeftCG: " << B.d << " " << B.e << " " << B.f << std::endl;
					dat11 << "LeftCG: " << B.g << " " << B.h << " " << B.i << std::endl;
					dat11 << "---------" << std::endl;
					dat11 << "Stress: " << T.a << " " << T.b << " " << T.c << std::endl;
					dat11 << "Stress: " << T.d << " " << T.e << " " << T.f << std::endl;
					dat11 << "Stress: " << T.g << " " << T.h << " " << T.i << std::endl;
					dat11 << "---------" << std::endl;

					dat11 << "Stress S: " << Ts.a << " " << Ts.b << " " << Ts.c << std::endl;
					dat11 << "Stress S: " << Ts.d << " " << Ts.e << " " << Ts.f << std::endl;
					dat11 << "Stress S: " << Ts.g << " " << Ts.h << " " << Ts.i << std::endl;
					
					dat11 << "---------" << std::endl;
					
					dat11 << "Stress V: " << Tv.a << " " << Tv.b << " " << Tv.c << std::endl;
					dat11 << "Stress V: " << Tv.d << " " << Tv.e << " " << Tv.f << std::endl;
					dat11 << "Stress V: " << Tv.g << " " << Tv.h << " " << Tv.i << std::endl;
					
					
					dat11 << "---------" << std::endl;
					
					dat11 << "GrowthTensor: " << tets[i].G.a << " " << tets[i].G.b << " " << tets[i].G.c << std::endl;
					dat11 << "GrowthTensor: " << tets[i].G.d << " " << tets[i].G.e << " " << tets[i].G.f << std::endl;
					dat11 << "GrowthTensor: " << tets[i].G.g << " " << tets[i].G.h << " " << tets[i].G.i << std::endl;
					
					dat11 << "---------" << std::endl;
					
					dat11 << "Ut0: " << Ut0[i].x << " " << Ut0[i].y << " " << Ut0[i].z << std::endl;
					
					dat11 << "---------" << std::endl;
					
					dat11 << "Utn1: " << Ut[n1].x << " " << Ut[n1].y << " " << Ut[n1].z << std::endl;
					dat11 << "Utn2: " << Ut[n2].x << " " << Ut[n2].y << " " << Ut[n2].z << std::endl;
					dat11 << "Utn3: " << Ut[n3].x << " " << Ut[n3].y << " " << Ut[n3].z << std::endl;
					dat11 << "Utn4: " << Ut[n4].x << " " << Ut[n4].y << " " << Ut[n4].z << std::endl;
					
					dat11 << "---------" << std::endl;
					
					dat11 << "Vtn1: " << Vt[n1].x << " " << Vt[n1].y << " " << Vt[n1].z << std::endl;
					dat11 << "Vtn2: " << Vt[n2].x << " " << Vt[n2].y << " " << Vt[n2].z << std::endl;
					dat11 << "Vtn3: " << Vt[n3].x << " " << Vt[n3].y << " " << Vt[n3].z << std::endl;
					dat11 << "Vtn4: " << Vt[n4].x << " " << Vt[n4].y << " " << Vt[n4].z << std::endl;
									
					dat11 << "---------" << std::endl;
					dat11 << "Ft: " << Ft[i].x << " " << Ft[i].y << " " << Ft[i].z << std::endl;
					dat11 << "Vt: " << Vt[i].x << " " << Vt[i].y << " " << Vt[i].z << std::endl;
					dat11 << "Ut: " << Ut[i].x << " " << Ut[i].y << " " << Ut[i].z << std::endl;
					dat11 << "---------" << std::endl;

					dat11 << "Vn0: " << Vn0[n1] << " " << Vn0[n2] << " " << Vn0[n3] << " " << Vn0[n4] << std::endl;
					dat11 << "Vn: " << Vn[n1] << " " << Vn[n2] << " " << Vn[n3] << " " << Vn[n4] << std::endl;
					dat11 << "---------" << std::endl;
					dat11 << "Ja: " << Ja << std::endl;
					dat11 << "J: " << J << std::endl;
					dat11 << "J1: " << J1 << std::endl;
					dat11 << "J2: " << J2 << std::endl;
					dat11 << "J3: " << J3 << std::endl;
					dat11 << "J4: " << J4 << std::endl;
					dat11 << "---------" << std::endl;
					dat11 << "vol: " << vol << std::endl;
				
					dat11 << "---------" << std::endl;
					dat11 << "us: " << us << std::endl;
					dat11 << "uv: " << uv << std::endl;
					dat11 << "W: " << W << std::endl;
					dat11 << "Ue: " << Ue << std::endl;
					dat11 << "Ue_xw: " << Ue_xw << std::endl;
					dat11 << "---------" << std::endl;
					dat11 << "powJ23: " << powJ23 << std::endl;
					dat11 << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-" << std::endl;
				}						
			}*/
			

        }
		
		auto count_4 = std::chrono::system_clock::now();

		dats2 << step << endl;
        // ----------------------------------------------------------------------
        
        // MIDPLANE COMPUTATION: Accounts for nodes crossing midplane (hemispheres) during deformation 
		int mpyL = 0; 
		int mpyR = 0;
		#pragma omp parallel for 
		for (int i = 0; i < nsn; i++) {
			int pt = SN[i];
			if ( Ut0[pt].y < mpy - 0.5*a && Ut[pt].y > mpy ) {          // if node crosses from left to right the mpy (Ut0 --> Ut), penalising force applied
				Ft[pt].y -= (mpy - Ut[pt].y)/hc*a*a*K;
                mpyL++; 
			}
			if ( Ut0[pt].y > mpy + 0.5*a && Ut[pt].y < mpy ) {          // if node crosses from right to left the mpy (Ut --> Ut0), penalising force applied
				Ft[pt].y -= (mpy - Ut[pt].y)/hc*a*a*K;
                mpyR++;
			}
		}
        auto count_5 = std::chrono::system_clock::now();

        // ----------------------------------------------------------------------
		double ymin = 1.0; 
		double ymax = -1.0; 
		double xmin = 1.0;
		double xmax = -1.0;
		double Volume, Uk;
		
        // Output
        if (step == 0 || (step%di == 0) || (floor(ga) != floor(lastga))) { 
            int ga_print = (int)floor(ga); 
			Uk = 0.0; // Kinetic energy
			for (int i = 0; i < nn; i++) Uk += 0.5*(Vn[i]*rho)*(Vt[i].dot(Vt[i]));
            Area = 0.0;
            for (int i = 0; i < nf; i++) {
                Vector N = (Ut[faces[i].n2]-Ut[faces[i].n1]).cross(Ut[faces[i].n3]-Ut[faces[i].n1]);
                Area += 0.5*N.length()*(gr[faces[i].n1] + gr[faces[i].n2] + gr[faces[i].n3])/3.0;
            }      

			Volume = 0.0; // Volume
			for (int i = 0; i < nn; i++) Volume += Vn[i];
			
            // Find new length
			for (int i = 0; i < nsn; i++) {
				xmin = min(xmin, Ut[SN[i]].x);
				xmax = max(xmax, Ut[SN[i]].x);
				ymin = min(ymin, Ut[SN[i]].y);
				ymax = max(ymax, Ut[SN[i]].y);
			}
			
			// Desnormalize mesh values using the same maxd factor used at the beginning of the codes. Original dimensions are be recovered.
			for (int i = 0; i < nn; i++){
				Ut1[i].x = (Ut[i].y)*maxd;
				Ut1[i].y = -(Ut[i].x)*maxd; 
				Ut1[i].z = -(Ut[i].z)*maxd;
				Ut1[i] = Ut1[i];
			}
			
			auto output0 = std::chrono::system_clock::now();
			std::chrono::duration<double> outputI = output0-start0;
            if (step == 0) {
			dat << setw(13) << "step" << setw(13) << "t" << setw(13) << "output" << setw(13) << "ga" << setw(13) << "Uk" << setw(13) << "Ue" << setw(13) << "Area" << setw(13) << "Volume" << setw(13) << "xmax-xmin" << setw(13) << "ymax-ymin" << endl;
			}
			dat << setw(13) << step << setw(13) << t << setw(13) << outputI.count() << setw(13) << ga << setw(13) << Uk << setw(13) << Ue << setw(13) << Area << setw(13) << Volume << setw(13) << xmax-xmin << setw(13) << ymax-ymin << endl;	
			
			//writeSTL(Ut1, faces, SN, SNb, nsn, step, hname, id, ga, spw);
            if (step == 0 || (floor(ga) != floor(lastga))) {
				spw++; 
				//if (step%dimesh == 0) {
				//writeMESH(Ut1, faces, tets, SNb, nn, ne, step, hname, id, ga, spw);
				writeSTL(Ut1, faces, SN, SNb, nsn, step, id, hname, ga, spw);
				// writeTXT(Ut1, faces, tets, SNb, nn, hname, id, ga);
				writeTXT(Ut1, faces, tets, SN, nsn, hname, id, ga);
				writeVTK(Ut1, faces, tets, mat_T1, mat_T2, mat_T3, mat_T_diag, mat_G1, mat_G2, mat_G3, mat_G_diag, contacts, mu_vec, gm_vec, gr, gr_vtk, d2s, Vt, Ft, SNb, nn, ne, step, hname, id, ga, spw);
			}
		}

        
		start0 = std::chrono::system_clock::now();
		auto count_6 = std::chrono::system_clock::now();
        
		// Newton dynamics
        #pragma omp parallel for 
        for (int i = 0; i < nn; i++) {
            Ft[i] -= Vt[i]*gamma*Vn0[i];
            Vt[i] += Ft[i]/(Vn0[i]*rho)*dt;
            Ut[i] += Vt[i]*dt;		
			Ft[i].clear();
		}
        
		// print values and define time stamp to control behavior
		auto count_7 = std::chrono::system_clock::now();		
		auto outputT = std::chrono::system_clock::now();
		
		double output = std::chrono::duration<double>(outputT - start).count();
		double tmdist = std::chrono::duration<double>(count_1 - count_0).count();
		double tpenforce = std::chrono::duration<double>(count_2-count_1).count();
		double tnodvol = std::chrono::duration<double>(count_3-count_2).count();
		double tdeform = std::chrono::duration<double>(count_4-count_3).count();
		double tmpy = std::chrono::duration<double>(count_5-count_4).count();
		double tprint = std::chrono::duration<double>(count_6-count_5).count();
        double tnewton = std::chrono::duration<double>(count_7-count_6).count();
		
		outputS += output;
		tmdistS += tmdist;
		tpenS += tpenforce;
		tnodS += tnodvol;
		tdefS += tdeform;
        tmpyS += tmpy; 
		tprintS += tprint;
		tnewtonS += tnewton;
		
		if (step == 0 || (step%dicout == 0) || (floor(ga) != floor(lastga))) { 
			if (step == 0) {
			cout  << setw(11) << "step" << setw(11) << "t" << setw(11) << "ga" << setw(11) <<"count(loop)" << setw(11) << "sumcount()" << setw(11) << "Uk" << setw(11) << "Ue" << setw(11) << 
				"Area" << setw(11) << "Volume" << setw(11) << "xmax-xmin" << setw(11) << "ymax-ymin" << endl; 
				
			datt << setw(13) << "step" << setw(13) << "t" << setw(13) << "ga" << setw(13) << "Time.Total" << setw(13) << "Time.Loop" << setw(13) << "Count.PenF" << setw(13) << "Count.MDist" << setw(13) << "Count.mpyR" << setw(13) << "Count.mpyL" << setw(13) <<
				"Time.MDist" << setw(13) << "Time.FnCalc." << setw(13) << "Time.NodalF" << setw(13) << "Time.Deform" << setw(13) << "Time.mpy" << setw(13) << "Time.Print" << setw(13) << "Time.Newton" << endl;
			}

		cout  << setw(11) << step << setw(11) << t << setw(11) << ga << setw(11) << output << setw(11) << outputS << setw(11) << Uk << setw(11) << Ue << setw(11) << Area << setw(11) << Volume << 
			setw(11) << xmax-xmin << setw(11) << ymax-ymin << endl;
		
		datt << setw(13) << step << setw(13) << t << setw(13) << ga << setw(13) << output << setw(13) << outputS << setw(13) << cpen << setw(13) << cmdist << setw(13) << mpyR << setw(13) << mpyL << 
			setw(13) << tmdistS << setw(13) << tpenS << setw(13) << tnodS << setw(13) << tdefS << setw(13) << tmpyS << setw(13) << tprintS << setw(13) << tnewtonS << endl;      
        }

		start = std::chrono::system_clock::now();
		

		dats3 << step << endl;
		
		t += dt;
        step++;
		lastga = ga; 
	}       
    dat.close();
    datgatime.close();
    datt.close();
    dats1.close();
    dats2.close();
    dats3.close();
    datp.close();
    return 0;
}

// ----------------------------------------------------------------------

// Generates point-triangle proximity lists using the linked cell algorithm
void createNNLtriangle(vector<int>* NNLt, Vector* Ut, vector<Face>& faces, int* SN, int nsn, int nf, double hs, double bw, double mw) {
	int mx = max(1, (int)(bw/mw)); // ** = 40 cells bw=3.2, mw=0.08
    //vector<int> head(mx*mx*mx, -1);
	//vector<int> list(nf);
	std::vector<int> head(mx*mx*mx, -1); //****** mx*mx*mx cells nomber, size mx*mx*mx vector with all values are -1, 40*40*40 = 64000
	std::vector<int> list(nf); // **** nf = 101882
	int xa, ya, za, xi, yi, zi;
	double ub, vb, wb;
	int pt, tri;
	Vector cog;
	Vector nn1;
	for (int i = 0; i < nf; i++) { // Divide triangle faces into cells, i index of face
		//std::cout<<"i is "<<i<<std::endl;
		//Vector cog = (Ut[faces[i].n1] + Ut[faces[i].n2] + Ut[faces[i].n3])/3.0;
        cog = (Ut[faces[i].n1] + Ut[faces[i].n2] + Ut[faces[i].n3])/3.0;
		//nn1 = Ut[faces[i].n1];
		//std::cout<<nn1.x<<" "<<nn1.y<<" "<<nn1.z<<std::endl;
		//std::cout<<"faces[i].n1 "<<faces[i].n1<<"faces[i].n2 "<<faces[i].n2<<"faces[i].n3 "<<faces[i].n3<<std::endl;
		int xa = (int)((cog.x + 0.5*bw)/bw*mx);
		int ya = (int)((cog.y + 0.5*bw)/bw*mx);
		int za = (int)((cog.z + 0.5*bw)/bw*mx);
	    //std::cout<<"cog.x "<<cog.x<<"cog.y "<<cog.y<<"cog.z "<<cog.z<<std::endl;
        int tmp =  mx*mx*za + mx*ya + xa; // *** 1641838 > 64000
	    //std::cout<<"tmp is "<<tmp<<std::endl;
		list[i]=head[mx*mx*za + mx*ya + xa];
		head[mx*mx*za + mx*ya + xa] = i;
	}

    // --------------------

    #pragma omp parallel for
    for (int i = 0; i < nsn; i++) { // Search cells around each point and build proximity list
        int pt = SN[i];
        NNLt[i].clear();
		//std::cout<<"CP"<<i<<std::endl;
        int xa = (int)((Ut[pt].x + 0.5*bw)/bw*mx);
        int ya = (int)((Ut[pt].y + 0.5*bw)/bw*mx);
        int za = (int)((Ut[pt].z + 0.5*bw)/bw*mx);
		//std::cout<<"cp1"<<std::endl;

		for (int xi = max(0, xa-1); xi <= min(mx-1, xa+1); xi++)		// *** Browse head list
		for (int yi = max(0, ya-1); yi <= min(mx-1, ya+1); yi++)
		for (int zi = max(0, za-1); zi <= min(mx-1, za+1); zi++) {
			int tri = head[mx*mx*zi + mx*yi + xi];
			while (tri != -1) {
				if ( pt != faces[tri].n1 && pt != faces[tri].n2 && pt != faces[tri].n3 ) {				
					if ( (closestPointTriangle(Ut[pt], Ut[faces[tri].n1], Ut[faces[tri].n2], Ut[faces[tri].n3], ub, vb, wb) - Ut[pt]).length() < hs) {
						NNLt[i].push_back(tri);
					}
				}
				tri = list[tri];
			}
		}
    }
}

// ----------------------------------------------------------------------


void Eigensystem(Matrix A, Matrix& V, double d[3]) {

	double A_[3][3];
	double V_[3][3];

	A_[0][0] = A.a; A_[0][1] = A.b; A_[0][2] = A.c;
	A_[1][0] = A.d; A_[1][1] = A.e; A_[1][2] = A.f;
	A_[2][0] = A.g; A_[2][1] = A.h; A_[2][2] = A.i;
	
	eigen_decomposition(A_, V_, d);

	V.a = V_[0][0]; V.b = V_[0][1]; V.c = V_[0][2];
	V.d = V_[1][0]; V.e = V_[1][1]; V.f = V_[1][2];
	V.g = V_[2][0]; V.h = V_[2][1]; V.i = V_[2][2];
}

// ----------------------------------------------------------------------

// Returns the closest point of triangle abc to point p ***** a or b or c, if not, pt projection through the barycenter inside the triangle 
Vector closestPointTriangle(Vector& p, Vector& a, Vector& b, Vector& c, double& u, double& v, double& w) {
    
    Vector ab = b - a;
    Vector ac = c - a;
    Vector ap = p - a;
    double d1 = ab.dot(ap);
    double d2 = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0) {
        u = 1.0;
        v = 0.0;
        w = 0.0;
        return a;
    }
    Vector bp = p - b;
    double d3 = ab.dot(bp);
    double d4 = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3) {
        u = 0.0;
        v = 1.0;
        w = 0.0;
        return b;
    }
    double vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        v = d1 / (d1 - d3);
        u = 1.0 - v;
        w = 0.0;
        return a + ab * v;
    }
    Vector cp = p - c;
    double d5 = ab.dot(cp);
    double d6 = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6) {
        u = 0.0;
        v = 0.0;
        w = 1.0;
        return c;
    }
    double vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        w = d2 / (d2 - d6);
        u = 1.0 - w;
        v = 0.0;    
        return a + ac * w;
    }
    double va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        u = 0.0;
        v = 1.0 - w;
        return b + (c - b) * w;
    }
    double denom = 1.0 / (va + vb + vc);
    v = vb * denom;
    w = vc * denom;
    u = 1.0 - v - w;
    return a + ab * v + ac * w;
}

// ----------------------------------------------------------------------

// Finds the nearest surface nodes (csn) to nodes and distances to them (d2s)
void dist2surf(Vector* Ut, int* SN, int nn, int nsn, int* csn, double* d2s) {
	int p, j;

	#pragma omp parallel for private(p, j)
	for (int i = 0; i < nn; i++) {
	    
    	double d2min = 1e9;
		for (j = 0; j < nsn; j++) {
			double d2 = (Ut[SN[j]] - Ut[i]).dot(Ut[SN[j]] - Ut[i]);
			if (d2 < d2min) {
				d2min = d2;
				p = j;
			}
		}
		csn[i] = p; //****** = SN[j]] ? to get index of the nearest surface node 
		d2s[i] = sqrt(d2min);
	}
}

// ----------------------------------------------------------------------

void writeSTL(Vector* Ut, vector<Face>& faces, int* SN, int* SNb, int nsn, int step, int id, int hname, double ga, int spw){
    char stlname[110];

	sprintf(stlname, "%s/%s/SD%d_h%d-GA_%.0f.stl", PATH_DIR, NAME, id, hname, ga);
	
	ofstream stl(stlname);
	stl.setf(ios::fixed);
	stl.precision(10);    
	// Normals
	Vector* N = new Vector[nsn];
	Vector Ntmp;
	stl << "solid " << "\r\n"; //'"'  << MESH_FILE << '"' << "\n";
	for (int i = 0; i < faces.size(); i++) {
		Ntmp = (Ut[faces[i].n2] - Ut[faces[i].n1]).cross(Ut[faces[i].n3] - Ut[faces[i].n1]);
		N[SNb[faces[i].n1]] += Ntmp;
		N[SNb[faces[i].n2]] += Ntmp;
		N[SNb[faces[i].n3]] += Ntmp;
		// FACES
		stl << "facet normal " << N[SNb[faces[i].n1]].x << " " << N[SNb[faces[i].n2]].y << " " << N[SNb[faces[i].n3]].z << "\r\n";
		stl << "outer loop" << "\r\n";
		// VERTICES
		stl << "vertex " << Ut[SN[faces[i].n1]].x << " " << Ut[SN[faces[i].n1]].y << " " << Ut[SN[faces[i].n1]].z << "\r\n";
		stl << "vertex " << Ut[SN[faces[i].n2]].x << " " << Ut[SN[faces[i].n2]].y << " " << Ut[SN[faces[i].n2]].z << "\r\n";
		stl << "vertex " << Ut[SN[faces[i].n3]].x << " " << Ut[SN[faces[i].n3]].y << " " << Ut[SN[faces[i].n3]].z << "\r\n";

		stl << "endloop" << "\r\n";
		stl << "endfacet" << "\r\n";
	}
	stl << "endsolid " << "\r\n";  //'"' << MESH_FILE << '"' << "\n";
	stl.close();
}

// ----------------------------------------------------------------------

void writeMESH(Vector* Ut, vector<Face>& faces, Tetra* tets, int* SNb, int nn, int ne, int step, int hname, int id, double ga, int spw){
   char meshname[100];
   sprintf(meshname, "%s/%s/D0%d_h%d_s%d-GA_%.0f.mesh", PATH_DIR, NAME, id, hname, step, ga);
   ofstream mesh(meshname);
   mesh.setf(ios::fixed);
   mesh.precision(6);    // Nodes
   mesh << nn << "\n";
   for (int i = 0; i < nn; i++) {
	mesh << Ut[i].x << " " << Ut[i].y << " " << Ut[i].z << endl;
   }   
   mesh << ne << "\n";
   for (int i = 0; i < ne; i++) {
	mesh << 1 << " " << tets[i].n1 + 1 << " " << tets[i].n2 + 1 << " " << tets[i].n3 + 1 << " " << tets[i].n4 + 1 << endl;
    }
   mesh << faces.size() << "\n";
   for (int i = 0; i < faces.size(); i++) {
       mesh << 1 << " " << SNb[faces[i].n1]+1 << " " << SNb[faces[i].n2]+1 << " " << SNb[faces[i].n3]+1 << "\n";
   }   mesh.close();
}

// ----------------------------------------------------------------------


void writeTXT(Vector* Ut, vector<Face>& faces, Tetra* tets, int* SN, int nsn, int hname, int id, double ga){
	
	char txtname[100];
	sprintf(txtname, "%s/%s/SD%d_h%d-GA_%.0f.txt", PATH_DIR, NAME, id, hname, ga);

	ofstream txt(txtname);
	txt.setf(ios::fixed);
	txt.precision(6);  
	
	txt << nsn << "\n";      
	for (int i = 0; i < nsn; i++) {
		txt << Ut[i].x << " " << Ut[i].y << " " << Ut[i].z << "\n";
	} 
	
	txt << faces.size() << "\n";
	for (int i = 0; i < faces.size(); i++) {
       txt << SN[faces[i].n1]+1 << " " << SN[faces[i].n2]+1 << " " << SN[faces[i].n3]+1 << "\n";
	}   
	txt.close();	
}


// ----------------------------------------------------------------------

void writeVTK(Vector* Ut1, vector<Face>& faces, Tetra* tets, Vector* mat_T1, Vector* mat_T2, Vector* mat_T3, Vector* mat_T_diag, Vector* mat_G1, Vector* mat_G2, Vector* mat_G3, Vector* mat_G_diag, double* contacts, double* mu_vec, double* gm_vec, double* gr, double* gr_vtk, double* d2s, Vector* Vt, Vector* Ft, int* SNb, int nn, int ne, int step, int hname, int id, double ga, int spw){
	char vtkname[100];
	sprintf(vtkname, "%s/%s/SD%d_h%d-GA_%.0f.vtk", PATH_DIR, NAME, id, hname, ga);

	ofstream vtk(vtkname);
	vtk.setf(ios::fixed);
	vtk.precision(6);  

	vtk << "# vtk DataFile Version 4.2" << "\n";
	vtk << "vtk simulation output" << "\n";
	vtk << "ASCII" << "\n";
	vtk << "DATASET UNSTRUCTURED_GRID" << "\n";
	vtk << "POINTS " << nn << " float" << "\n";                       				// Nodes
	for (int i = 0; i < nn; i++) {
		vtk << Ut1[i].x << " " << Ut1[i].y << " " << Ut1[i].z << "\n";
	} 
	
	vtk << "\n";
	vtk << "CELLS " << ne << " " << ne*5 << "\n";                       		// Tets
	for (int i = 0; i < ne; i++) {
		vtk << 4 << " " << tets[i].n1 << " " << tets[i].n2 << " " << tets[i].n3 << " " << tets[i].n4 << "\n";
    }
	
	vtk << "\n";
	vtk << "CELL_TYPES " << ne << "\n";                       					// Tets
	for (int i = 0; i < ne; i++) {
	    vtk << 10 << "\n";
    }
	/*
	vtk << "\n";
	vtk << "POINT_DATA " << nn << "\n";        
	vtk << "SCALARS drawem float 1" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (int i = 0; i < nn; i++) {
		vtk << draw[i] << "\n";
	}
	*/
	vtk << "\n";  	
	vtk << "POINT_DATA " << nn << "\n"; 
	vtk << "SCALARS contacts float 1" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (int i = 0; i < nn; i++) {
		if (isnan(contacts[i]) == 1 ) { contacts[i] = 0.000; }
		//if (contacts[i] == 'nan') { contacts[i] = 0.000; }
		vtk << contacts[i] << "\n";
	} 	
	
	vtk << "\n";      
	vtk << "SCALARS gr_indicator float 1" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (int i = 0; i < nn; i++) {
		vtk << gr[i] << "\n";
	} 
	/*
	vtk << "\n";
	vtk << "SCALARS gr_binary float 1" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (int i = 0; i < nn; i++) {
		vtk << gr_vtk[i] << "\n";
	} 
	
	vtk << "\n";
	vtk << "SCALARS d2s float 1" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (int i = 0; i < nn; i++) {
		vtk << d2s[i] << "\n";
	} 
	*/ 
	
	vtk << "\n";
	vtk << "VECTORS velocities float" << "\n";
	for (int i = 0; i < nn; i++) {
		vtk << Vt[i].x << " " << Vt[i].y << " " << Vt[i].z << "\n";
	} 
	
	vtk << "\n";
	vtk << "VECTORS forces float" << "\n";
	for (int i = 0; i < nn; i++) {
		vtk << Ft[i].x << " " << Ft[i].y << " " << Ft[i].z << "\n";
	} 
	
	vtk << "\n";
	vtk << "CELL_DATA " << ne << "\n";   
	
	/*
	vtk << "SCALARS gm float 1" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (int i = 0; i < ne; i++) {
		vtk << gm_vec[i] << "\n";
	}
	
	vtk << "\n";
	vtk << "SCALARS mu float 1" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (int i = 0; i < ne; i++) {
		vtk << mu_vec[i] << "\n";
	}

	vtk << "\n"; 
	vtk << "VECTORS growth_diag float" << "\n";
	for (int i = 0; i < ne; i++) {
		vtk << mat_G_diag[i].x << " " << mat_G_diag[i].y << " " << mat_G_diag[i].z << "\n";
	}
	*/	
	vtk << "\n"; 
	vtk << "VECTORS stress_diag float" << "\n";
	for (int i = 0; i < ne; i++) {
		vtk << mat_T_diag[i].x << " " << mat_T_diag[i].y << " " << mat_T_diag[i].z << "\n";
	}
	/*
	vtk << "\n"; 
	vtk << "VECTORS growth_G1 float" << "\n";
	for (int i = 0; i < ne; i++) {
		vtk << mat_G1[i].x << " " << mat_G1[i].y << " " << mat_G1[i].z << "\n";
	}
		
	vtk << "\n"; 
	vtk << "VECTORS growth_version float" << "\n";
	for (int i = 0; i < ne; i++) {
		vtk << mat_G1[i].x << " " << mat_G2[i].x << " " << mat_G3[i].x << "\n";
	}
	
	vtk << "\n";   
	vtk << "TENSORS growth float" << "\n";
	for (int i = 0; i < ne; i++) {
		vtk << mat_G1[i].x << " " << mat_G1[i].y << " " << mat_G1[i].z << "\n";
		vtk << mat_G2[i].x << " " << mat_G2[i].y << " " << mat_G2[i].z << "\n";
		vtk << mat_G3[i].x << " " << mat_G3[i].y << " " << mat_G3[i].z << "\n\n";
	}
	
	vtk << "\n";   
	vtk << "TENSORS stress float" << "\n";
	for (int i = 0; i < ne; i++) {
		vtk << mat_T1[i].x << " " << mat_T1[i].y << " " << mat_T1[i].z << "\n";
		vtk << mat_T2[i].x << " " << mat_T2[i].y << " " << mat_T2[i].z << "\n";
		vtk << mat_T3[i].x << " " << mat_T3[i].y << " " << mat_T3[i].z << "\n\n";
	}
	
	
	
	vtk << "CELL_DATA" << ne << "\n";        
	vtk << "VECTORS forces float 1" << "\n";
	//vtk << "LOOKUP_TABLE default" << "\n";
	for (int i = 0; i < nn; i++) {
		vtk << Ft[i].x << " " << Ft[i].y << " " << Ft[i].z << "\n";
	} 
	*/
	vtk.close();
}

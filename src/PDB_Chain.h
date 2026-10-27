#pragma once
#include "GeometryUtils.h"
#include "PDB_Residue.h"
#include "PdbUtils.h"
#include <vector>
#include <string>
#include <sstream>
using namespace std;


//--------- type_def -----------//
typedef struct _Phi_Psi_Omega
{
	void operator=(const struct _Phi_Psi_Omega & phi_psi_omega)
	{
		phi = phi_psi_omega.phi;
		psi = phi_psi_omega.psi;
		omega = phi_psi_omega.omega;
	}
	double phi;
	double psi;
	double omega;
}Phi_Psi_Omega;
typedef struct _Theta_Tau
{
	void operator=(const struct _Theta_Tau & theta_tau)
	{
		theta = theta_tau.theta;
		tau = theta_tau.tau;
	}
	double theta;
	double tau;
}Theta_Tau;



//====class_PDB_Chain====// 
class PDB_Chain
{
public:
	PDB_Chain();
	virtual ~PDB_Chain();
	PDB_Chain(const PDB_Chain & chain);
	PDB_Chain & operator =(const PDB_Chain & chain);

//======== variables ========//
public:
	char chain_id;
	int length;			
	vector<PDB_Residue> residues;
	// [0].phi and [n-1].psi, [n-1].omega are meaningless.
	vector<Phi_Psi_Omega> phi_psi_omegas;
	// [0], [1], [2].tau are meaningless
	vector<Theta_Tau> theta_taus;
	//other sequence
	string SSE;       // SSE structure
	string CLE;       // CLE structure
	string SEQ;       // of amino acid

//======== functions ========//
public:
	//init & calc
	int initialize_simple(int len, char id);
	void allocate(void);
	void allocate(int length);

	//------- main process part --------// (virtual)
	//phi_psi <-> XYZ
	virtual int refold_from_phi_psi_omega(void); //from phi_psi_omega to 3D-coordinates
	virtual int refold_from_theta_tau(void);     //from theta_tau to 3D-coordinates
	virtual int calculate_phi_psi_omega(void);   //3D-coordinates to phi_psi_omega
	virtual int calculate_theta_tau(void);       //3D-coordinates to theta_tau
	//calculate CLE,SSE and SEQ
	virtual int calculate_CLE(void);  //calculate Conformational Letter (17 digit)
	virtual int calculate_SSE(void);  //calculate DSSP Secondary Structure (8 digit)

	//--- get and set ---//
	//[main_part]
	int get_length() const;
	int get_residue(int i, PDB_Residue & residue) const;
	int set_residue(int i, const PDB_Residue & residue);
	//[vice_part]
	//normal
	char get_chain_id() const;
	//special

	//--- printf related ---//
};
//====class_PDB_Chain====//over 

#include "PDB_Residue.h"
#include "StringUtils.h"
#include "PDB_Chain.h"
#include "PdbUtils.h"
#include <iomanip>
using namespace std;


//-------- PDB_Chain_Define ---------// 
PDB_Chain::PDB_Chain()
{
	chain_id='A';
	length=0;
	SSE="";
	CLE="";
	SEQ="";
}
PDB_Chain::~PDB_Chain()
{
}
PDB_Chain::PDB_Chain(const PDB_Chain & chain)
{
	this->chain_id = chain.chain_id;
	this->length = chain.length;
	this->residues = chain.residues;
	this->phi_psi_omegas = chain.phi_psi_omegas;
	this->theta_taus = chain.theta_taus;
	this->SSE.assign(chain.SSE);
	this->CLE.assign(chain.CLE);
	this->SEQ.assign(chain.SEQ);
}
PDB_Chain &  PDB_Chain::operator =(const PDB_Chain & chain)
{
	if(this==&chain)return *this;
	this->chain_id = chain.chain_id;
	this->length = chain.length;
	this->residues = chain.residues;
	this->phi_psi_omegas = chain.phi_psi_omegas;
	this->theta_taus = chain.theta_taus;
	this->SSE.assign(chain.SSE);
	this->CLE.assign(chain.CLE);
	this->SEQ.assign(chain.SEQ);
	return *this;
}

//=================== main process part ====================//
//[init & calc]
//init
int PDB_Chain::initialize_simple(int len, char id)
{
	length = len;
	allocate(len);
	chain_id = id;
	return 0;
}
void PDB_Chain::allocate(void)
{
	residues.resize( length );
	phi_psi_omegas.resize( length );
	theta_taus.resize( length );
}
void PDB_Chain::allocate( int n )
{
	length = n;
	allocate();
}

//============ virtual functions =============//
//phi_psi <-> XYZ
int PDB_Chain::refold_from_phi_psi_omega(void)
{
	return UNSUPPORT_ERROR;
}
int PDB_Chain::refold_from_theta_tau(void)
{
	return UNSUPPORT_ERROR;
}
int PDB_Chain::calculate_phi_psi_omega(void)
{
	return UNSUPPORT_ERROR;
}
int PDB_Chain::calculate_theta_tau(void)
{
	return UNSUPPORT_ERROR;
}
//calculate CLE,SSE and SEQ
int PDB_Chain::calculate_CLE(void)
{
	return UNSUPPORT_ERROR;
}
int PDB_Chain::calculate_SSE(void)
{
	return UNSUPPORT_ERROR;
}

/*------- Learning -------
//The difference between "at" and operator[] is that, 
//"at" does range checks. Otherwise, they are the same.
-----------------------------*/
//================================= get and set ============================//
//[main_part]
int PDB_Chain::get_length() const
{
	return length;
}
//residue
int PDB_Chain::get_residue(int i, PDB_Residue & residue) const
{
	if(i < 0 || i >= (int)residues.size())
	{
		debug_info(DEBUG_WARNING, "@PDB_Chain::get_residue, index = %d, excced %d\n",
				i, residues.size());
		return -1;
	}
	residue = residues.at(i);
	return 0;
}
int PDB_Chain::set_residue(int i, const PDB_Residue & residue)
{
	if(i < 0 || i >= (int)residues.size())
	{
		debug_info(DEBUG_WARNING, "@PDB_Chain::set_residue, index = %d, excced %d\n",
				i, residues.size());
		return -1;
	}
	residues.at(i) = residue;
	return 0;
}
//phi_psi_omega
//theta_tau


//[vice_part]
char PDB_Chain::get_chain_id() const
{
	return chain_id;
}
//sequence set & get (they are different with get_SEQ!!!)

//================================= printf related ============================//



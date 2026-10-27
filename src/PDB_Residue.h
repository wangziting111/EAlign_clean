#pragma once
#include "ResidueGeometry.h"
#include <string>
#include <vector>
using namespace std;

 //====class_PDB_Residue====// 
class PDB_Residue 
{
public:
	PDB_Residue(void);
	~PDB_Residue(void);
	PDB_Residue(const PDB_Residue & pdb_residue);
	PDB_Residue & operator =(const PDB_Residue & pdb_residue);

//======== variables ========//
private:
	char AA_residue;               //aminoacid type 'A->Z'
	string PDB_residue_number;     //number in PDB_Residue, exactly 6 chars
	ResidueGeometry backbone;   //backbone (N,CA,C,O)
	ResidueGeometry sidechain;  //sidechain (maximal sidechain is 10 for W)

public:
	vector <XYZ> hydro_rec;        //hydrogen record

//======== functions ========//
public:
	//init 
	int PDB_residue_backbone_initialize(char c);
	int PDB_residue_sidechain_initialize(char c);
	//check 
	int PDB_residue_backbone_check(int check = 4); //[N,CA,C,O](backbone)
	int PDB_residue_sidechain_check(void);
	int PDB_residue_CB_check(void);

	//[backbone] (set & get)
	//(atom_name mode)
	int get_backbone_atom( const string &atom_name, XYZ &atom, int & numb, double & R, double & tmpr ) const;
	int get_backbone_atom( const string &atom_name, XYZ &atom ) const;
	int get_backbone_atom( const string &atom_name, int &numb ) const;
	//(index mode)
	int set_backbone_atom( int index, const XYZ &atom, int numb = 0, double R = 1.0, double tmpr = 20.0 );
	int get_backbone_atom( int index, XYZ &atom, int & numb, double & R, double & tmpr ) const;
	int get_backbone_atom( int index, XYZ &atom ) const;
	int get_backbone_atom( int index, int &numb ) const;
	//[sidechain] (set & get)
	//(atom_name mode)
	int get_sidechain_atom( const string &atom_name, XYZ &atom, int & numb, double & R, double & tmpr ) const;
	int get_sidechain_atom( const string &atom_name, XYZ &atom ) const;
	int get_sidechain_atom( const string &atom_name, int &numb ) const;
	//(index mode)
	int set_sidechain_atom( int index, const XYZ &atom, int numb = 0, double R = 1.0, double tmpr = 20.0 );
	int get_sidechain_atom( int index, XYZ &atom, int & numb, double & R, double & tmpr ) const;
	int get_sidechain_atom( int index, XYZ &atom ) const;
	int get_sidechain_atom( int index, int &numb ) const;

	//[set and get] (others)
	int get_backbone_totnum() const;
	int get_sidechain_totnum() const;
	int get_backbone_part_index(int i) const;
	int get_sidechain_part_index(int i) const;
	char get_AA() const;
	int set_PDB_residue_number(const string & res_num);
	int get_PDB_residue_number(string &res_num) const;
	
	//-> PDB_Residue to XYZ_array + sidechain_count //__130830__//
};
//====class_PDB_Residue====//over 

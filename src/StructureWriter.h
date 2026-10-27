#pragma once
#include "PDB_Residue.h" 
#include "PdbUtils.h"
#include "GeometryUtils.h" 

//====class: StructureWriter====//
//=> output structure
class StructureWriter
{
public:
	StructureWriter(void);
	~StructureWriter(void);

//---- macros ----//
public:
	int OutPrt;       // whether output printf    // (default:no)

//--------- functions ---------//
public:
	//---PART_IV:  output_related
	// Referenced by an inherited routine that is discarded from the final binary.
	int Output_PDB(FILE *fp,int len,XYZ *mol,char *AMI,char *ind,char Chain_ID);
	int Output_PDB_III(FILE *fp,int len,PDB_Residue *mol,char Chain_ID,
		int OutType=1,int OutMode=1,int OutGlys=1,int OutLast=1);
};

#pragma once
#include "ConformationLetters.h"
#include "TMAligner.h"
#include "SequenceScorer.h"
#include "DynamicProgramming.h"

class EnvironmentScorer : virtual public TMAligner, virtual public SequenceScorer
{
public:
	EnvironmentScorer(int num=PROT_MAX_NUM,int CLESUM=1);
	~EnvironmentScorer(void);
	int Environment_maximal;

//---- variables ----//
public:
	//macro
	int Environment_Do_TMalign;       // whether use TMalign to refine (default: no)
	//normal parameter
	double Environment_Distance;      // the distance cutoff (default: 8.5)
	int Environment_mol_extend;       // determine the initial rotmat (default: 5)
	int Environment_neib_extend;      // determine the addtional neibor (default: 0)
	//DynaProg parameter
	double Environment_AMI_GapOpen;   // the gap open penalty for AMI
	double Environment_AMI_GapExtend; // the gap extend penalty for AMI
	double Environment_CLE_GapOpen;   // the gap open penalty for CLE
	double Environment_CLE_GapExtend; // the gap extend penalty for CLE
	//input
	int Environment_moln1;
	int Environment_moln2;
	XYZ *Environment_mol1;
	XYZ *Environment_mol2;
	int *Environment_ami1;
	int *Environment_ami2;
	int *Environment_cle1;
	int *Environment_cle2;
	//temp
	vector <vector <int> > Environment_Neib1;
	vector <vector <int> > Environment_Neib2;
	double *Environment_DP_sco;  //DynaProg matrix
	XYZ *Environment_tmp1;
	XYZ *Environment_tmp2;
	int *Environment_alit;
	int *Environment_ali2;
	char *Environment_cle;

//---- functions ----//
public:
	//init function
	void InitializeEnvironmentScorer(int maxlen);
	void ReleaseEnvironmentScorer(void);
	//vice function
	//misc function
	//main function
	//Out_Type: [0] mol, [1] ami, [2] cle
};

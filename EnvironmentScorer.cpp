#include "EnvironmentScorer.h"

//------ constructor -------//
EnvironmentScorer::EnvironmentScorer(int num,int CLESUM)
:TMAligner(num),SequenceScorer(CLESUM)
{
	//init
	Environment_maximal=num;
	InitializeEnvironmentScorer(Environment_maximal);
	//macro
	Environment_Do_TMalign=0;      // whether use TMalign to refine (default: no)
	//normal parameter
	Environment_Distance=8.5;      // the distance cutoff (default: 8.5)
	Environment_mol_extend=4;      // determine the initial rotmat (default: -l,+l)
	Environment_neib_extend=0;     // determine the addtional neibor (default: -k,+k)
	//DynaProg parameter
	Environment_AMI_GapOpen=-10;   // the gap open penalty for AMI
	Environment_AMI_GapExtend=-1;  // the gap extend penalty for AMI
	Environment_CLE_GapOpen=-100;  // the gap open penalty for CLE
	Environment_CLE_GapExtend=-10; // the gap extend penalty for CLE
}
EnvironmentScorer::~EnvironmentScorer(void)
{
	ReleaseEnvironmentScorer();
}

//--------- init ---------//
void EnvironmentScorer::InitializeEnvironmentScorer(int maxlen)
{
	Environment_DP_sco=new double[maxlen*maxlen];
	//input
	Environment_mol1=new XYZ[maxlen];
	Environment_mol2=new XYZ[maxlen];
	Environment_ami1=new int[maxlen];
	Environment_ami2=new int[maxlen];
	Environment_cle1=new int[maxlen];
	Environment_cle2=new int[maxlen];
	//temp
	Environment_tmp1=new XYZ[maxlen];
	Environment_tmp2=new XYZ[maxlen];
	Environment_alit=new int[maxlen];
	Environment_ali2=new int[maxlen];
	Environment_cle=new char[maxlen+1];
}
void EnvironmentScorer::ReleaseEnvironmentScorer(void)
{
	delete [] Environment_DP_sco;
	//input
	delete [] Environment_mol1;
	delete [] Environment_mol2;
	delete [] Environment_ami1;
	delete [] Environment_ami2;
	delete [] Environment_cle1;
	delete [] Environment_cle2;
	//temp
	delete [] Environment_tmp1;
	delete [] Environment_tmp2;
	delete [] Environment_alit;
	delete [] Environment_ali2;
	delete [] Environment_cle;
}

//--------- vice ---------//

//--------- misc ---------//

//------------------ main_function --------------//
//input original structures
//user may assign his own Neib_Definition for this input
//Out_Type: [0] mol, [1] ami, [2] cle

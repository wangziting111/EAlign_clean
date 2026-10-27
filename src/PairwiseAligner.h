#pragma once
#include "SeedAligner.h"
#include "EnvironmentScorer.h"
#include "TMRefiner.h"
using namespace std;

//====class: PairwiseAligner====//
//=> end-to-end pairwise structure alignment
// Using Zheng Wei-Mou's method
class PairwiseAligner :  public EnvironmentScorer, public SeedAligner, public TMRefiner
{
public:
	PairwiseAligner(int num=PROT_MAX_NUM,int CLESUM=1);
	~PairwiseAligner(void);

	//---------- data structure ------------//
	int Kill_Frag;    //whether kill frag or not
	int Kill_Frag_Para_0;  //-> min_small_len, 25
	int Kill_Frag_Para_1;  //-> gap_len, 25
	int Kill_Frag_Para_2;  //-> nongap_len, 15
	int Pairwise_Kill_Gap;  //whether kill gap or not
	vector <Align_Record> Pairwise_align_tot;

	//----------- functions -----------//
	//single match kill
	int IN_Ali_To_Cor(int *Fragment_Cor,int thres,int moln1,int moln2,int *ali1,int *ali2);
	void IN_Elongation(double *rotmat_,double thres,int cutoff,int *ali1,int *ali2,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2);
	void Pairwise_Align_Refine(XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int *ali2_in,int *ali2_out,double *rotmat_);
	double Fragment_Make_Score(Align_Record & align_record,XYZ *mol1,XYZ *mol2,int moln1,int moln2,double *rotmat,int *ali2);
	void Fragment_Single_Match_Kill(XYZ *mol1,XYZ *mol2,int moln1,int moln2,vector <Align_Record> &tot);
	//main function
	void Pairwise_Get_Initial_CLE(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2); //CLE dynamic programming
	void Alimeng_Add(vector <Align_Record> &align_tot,int moln2,double tms);
	double Pairwise_Align_Total(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_,
		int norm_len,double norm_d0,double *MAXSCO=0); //MAXSCO should at least double[8]
	double Pairwise_Align_Normal(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_,
		int norm_len,double norm_d0,double *MAXSCO=0); //MAXSCO should at least double[8]
};

#pragma once
#include "PdbUtils.h"
#include "FragmentAligner.h"
#include "SequenceScorer.h"
#include <algorithm>
#include <vector>
using namespace std;

//====class: SeedAligner====//
//=> conformation-letter and 3Di seed generation
// Using Zheng Wei-Mou's method
class SeedAligner : public FragmentAligner, virtual public SequenceScorer
{
public:
	SeedAligner(int num=PROT_MAX_NUM,int CLESUM=1);
	~SeedAligner(void);

//---- variables ----//
public:
	//macro
	int SFP_Strategy;    // 1:CLE, 2:GEN, 3:AA, 4:AA+CLE, 5:AA+3Di, 6:3Di
	int REF_Strategy;    //-> using refinement strategy
	//parameter
	int SFP_H_Len;
	int SFP_L_Len;
	int SFP_H_Thres;
	int SFP_L_Thres;
	//additional
	int ZM_Upper_Max;    //-> default is 500
	int ZM_Lower_Max;    //-> default is 50
	int ZM_TopL;         //-> retain TopL candidate for TMscore refinement
	double Refine_Cut;   //-> TM_refine cutoff (default: 0.3)
	//structure
	char *IN_CLE1;
	char *IN_CLE2;
	char *IN_AMI1;
	char *IN_AMI2;
	const int *IN_3DI1;
	const int *IN_3DI2;
	int *IN_DUM1;
	int *IN_DUM2;

//---- functions ----//
public:
	//[input]
	void Fragment_Input_Func(char *c1,char *c2,char *a1,char *a2,
		XYZ *m1,XYZ *m2,int n1,int n2);
	void Fragment_3Di_Input(const int *three_di1,const int *three_di2);
	//[process]
	void Fragment_Second_Init(void);
	double Fragment_Second_Refine(vector <Align_Record> &tot);
	//[main]
	void Fragment_Main_Init(int SFP_Strategy);
	double Fragment_Main_Func(vector <Align_Record> &tot);
	//[mask]
	int Check_Mask(vector <SFP_Record> &in, vector <SFP_Record> &out, int *mas1, int*mas2);
	void Mask_CLE_BreakPoint(char *CLE, int *mask, int len);
};

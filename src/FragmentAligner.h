#pragma once
#include <vector>
#include "AlignmentRefiner.h"
#include "FragmentChainer.h"
#include "AlignmentPipeline.h"
#include "TMAligner.h"
using namespace std;

//====class: FragmentAligner====//
//=> CLeFAPS main algorithm
class FragmentAligner : public AlignmentRefiner, virtual public FragmentChainer, virtual public AlignmentPipeline, virtual public TMAligner
{
public:
	FragmentAligner(int num=PROT_MAX_NUM);
	~FragmentAligner(void);
	int Fragment_Maximal;

//---- variables ----//
public:
	//parameter
	int ANG_CUT;            // vector angle threshold
	double Tail_Check;      // tail alignment threshold
	double Comp_Check;      // comparison threshold
	int Comp_Check_Vari;    // comparison partial (default: 0)
	int Geo_Sparse_DP;      // sparse DP score: 0=CombinedScore, 1=geometry only
	int SparseDP_filter;
	//additional
	double CUR_MaxJ_thres;
	double CUR_MaxJ;        // TopJ maximal
	double CUR_MaxK;        // TopK maximal
	int CUR_INDEX;          // current index
	int INDEX_CUR;          // index current
	double TM_d0;           // d0 to calculate TMscore
	int TM_smaller;         // smaller to calculate TMscore
	//temp
	int *tali1,*tali2;      // tali1 & tali2
	int *Fragment_Cor;           // Fragment_Cor (rigid frag correspondence)  [once]
	int *Fragment_Cor_temp;      // Fragment_Cor_temp (rigid frag correspondence) [once]
	double *center_mat;     // SFP_High rotmat
	int *head_score;        // SFP_High score
	int *head_index;        // SFP_High index 
	short *valid_idx;       // SFP_Low index
	vector <int> align_rec; // record alignment //__110121__//

//---- functions ----//
public:
	//[init function]
	void Fragment_Init(int maxlen);
	void Fragment_Dele(void);
	//[vice process]
	void Elongation(double *rotmat_,double thres,int cutoff,int *ali1,int *ali2,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2);                             //elongation
	void Shrinking(double *rotmat_,double thres,int cutoff,int *ali1,int *ali2,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2);                             //shrinking
	double Alg_Comp_Simp(vector <int> &ali_m1,vector <int> &ali_m2,int part=0); //alignment comparison

//---- virtual functions ----//
public:
	//---- input/output functions ------//io
	void Pipeline_Input(XYZ *m1,XYZ *m2,int n1,int n2) override;
	//---- make_socre virtual functions ----//score
	//given alignment and rotmat, return scores (maybe TMscore,RMSD and LALI, etc)
	double Pipeline_Make_Score(Align_Record & align_record) override;
	//given the align_record, make the secondary score (maybe LALI*TMscore, etc)
	double Pipeline_Get_Score(Align_Record &align) override;
	
	//---- virtual functions -----//partial
	//0. Initialize fragment alignment
	void Other_Fragment_Init(void) override;
	//1. ZoomIn Add (ret_rot might be updated)
	int ZoomIn_Add(int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP,double *ret_rot) override;
	//2 .Refinement (ret_rot might be updated)
	int Refinement(int recur,double FIN_CUT,double *ret_rot) override;
	//3. Kill NonLinear (ret_rot won't be updated)
	int Kill_NonLinear(double *ret_rot) override;
	//4. Final Stage (combine Shrink & Elong)
	int Final_Stage(double *ret_rot) override;
	//5. TMscore refinement (ret_rot might be updated)
	int Final_Refine(double *ret_rot) override;

	//---- virtual functions -----//global
	//0. Init functions
	void Global_Pipeline_Init(void) override;
	//1. Select best pivot
	int Select_Best_Pivot(int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP,int TopJ) override;
	double Pipeline_Full(XYZ *m1,XYZ *m2,int n1,int n2,
		vector <SFP_Record> &SFP_High,vector <SFP_Record> &SFP_Low,
		vector <Align_Record> &align_records) override;
	//2. Break Check
	int Pipeline_Break_Check(int cur,int count,int totnum) override;
	//3. Get Index
	int Pipeline_Get_Index(int cur) override;
	//4. Head Check
	int Pipeline_Head_Check(SFP_Record &SFP) override;
	//5. Get initial rotmat
	void Pipeline_Initial_Rotmat(int index,double *rotmat) override;
	//6. Tail Check
	int Pipeline_Tail_Check(Align_Record &cur) override;
	//7. Comparison Check
	int Pipeline_Comparison_Check(Align_Record &cur,vector <Align_Record> &tot) override;
	//8. Terminal Process
	void Pipeline_Terminal_Process(Align_Record &cur) override;
};

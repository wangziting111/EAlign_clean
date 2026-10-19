#pragma once
#include "AlignmentGeometry.h"



//====class: AlignmentRefiner====//
//=> one-to-multi correspondence (ali3)
class AlignmentRefiner : virtual public AlignmentGeometry
{
public:
	AlignmentRefiner(int num=PROT_MAX_NUM);
	~AlignmentRefiner(void);
	int Refine_Maximal;

//Ali3:DataStructure//
//	ali3[i][*][*] the i-th index of mol2(smaller)
//	ali3[i][0][0] the total correspondece of i-th index(default:0)
//	ali3[i][0][1] the minimal correspondence of i-th index(default:0,then k)
//  ali3[i][k][*] the k-th correspondence of the i-th index
//	ali3[i][k][0] the k-th correspondence of position in mol1(bigger)
//	ali3[i][k][1] the k-th correspondence of distance in mol1(bigger)(integer)
public:
	int *ali3;             //ali3 universal correspondence
	int ali3_TOT;          //ali3's volumn (default:6)
//pre_sco:DataStructure//
//	pre_sco[i][*][*] the i-th index of mol2(smaller)
//	pre_sco[i][0][*] the null state
//	pre_sco[i][k][*] the k-th state
//	pre_sco[i][x][0] the dynamic programming score
//	pre_sco[i][x][1] the correspondence in mol1
public:
	int *pre_sco;          // ali3 dynamic programming data_structure
	int *pre_temp;         // ali3 dynamic programming data_structure(temp)
public:
	int BAK_GAP;           // default 50   (Refine_DynaProg Bak_Gap penalty)
	int FOR_GAP;           // default 200  (Refine_DynaProg For_Gap penalty)
	int SCALE;             // default 100  (Refine_DynaProg Scaling Factor)


//====== ali3 function_related ======//
public:
	//--create & init--//
	void Refine_Create(int totlen);                           // ali3_Create
	void Refine_Delete(void);                                 // ali3_Delete
	void Refine_Init(int *ali3,int moln2);				    // ali3_Init (complex)
	//---- ali3_kabsch ----//
	double ali3_kabsch(double *rotmat_,double &rmsd,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int *ali3);                                         // ali3_rotation (complex)
	//--main_usage--//
	void Refine_Add(int ii,int jj,int len,double *rotmat_,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int *ali3,int *zali1=0,int *zali2=0);               // ali3_Add (complex)
	void Refine_Cor_Add(int *Fragment_Cor,int Range,double *rotmat_,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int *ali3,int *zali1=0,int *zali2=0);               // ali3_Cor_Add (complex)
	//--output_part--//


//---- ali3 DynamicProgramming--//
public:
	//--main_part--//
	//vector_score
	int Refine_Vector_Score_Forward(int ii,int jj,double *rotmat_,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int SCALE);                                                 // ali3_Vector_Forward_score (complex)
	int Refine_Vector_Score_Backward(int ii,int jj,double *rotmat_,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int SCALE);                                                 // ali3_Vector_Backward_score (complex)
	//ali3_DynaProg
	void Refine_DynaProg(double *rotmat_,double beta,double d_0,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int FOR_GAP,int BAK_GAP,int SCALE,int *ali3);               // ali3_Dynamic_Programming (complex)
	//ali3_TraceBack
	double Refine_TraceBack(double *rotmat_,double d_0,
		XYZ *mol1,XYZ *mol2,int moln1,int moln2,
		int *ali3,int *ali1,int *ali2);                             // ali3_Trace_Back (complex)
	//--output_part--//
};
//====class: AlignmentRefiner====//over

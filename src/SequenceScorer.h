///////////////////////////////////////////////////////////////////////////
//  CLePAPS:                                                             //
//  Conformational Letter based Pairwise Alignment of Protein Structure  //
//                                                                       //
//  Authors: S.Wang & WM.Zheng                                           //
///////////////////////////////////////////////////////////////////////////
/*
 
Copyright (c)  2007-2008  Institute of Theoretical Physics,Academia Sinica
All Rights Reserved
 
Permission to use, copy, modify and distribute any part of this CLePAPS
software for educational, research and non-profit purposes, without fee,
and without a written agreement is hereby granted, provided that the above
copyright notice, this paragraph and the following three paragraphs appear
in all copies.
 
Those desiring to incorporate this CLePAPS Software into commercial products
or use for commercial purposes should contact Institute of Theoretical Physics,
Academia Sinica,Beijing 100190, China. 
EMAIL: zheng@itp.ac.cn, wangsheng@itp.ac.cn
 
IN NO EVENT SHALL ITP(INSTITUTE OF THEORETICAL PHYSICS) BE LIABLE TO ANY PARTY 
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
LOST PROFITS, ARISING OUT OF THE USE OF THIS CLEPAPS SOFTWARE, EVEN IF ITP HAS 
BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
THE CLEPAPS SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND ITP HAS NO 
OBLIGATION TO PROVIDE MAINTENANCE,SUPPORT, UPDATES, ENHANCEMENTS, OR 
MODIFICATIONS.  ITP MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF 
ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT
THE USE OF THE CLEPAPS SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR
OTHER RIGHTS.
*/

#pragma once
#include <vector>
#include "AlignmentTypes.h"
#include "GeometryUtils.h"


//==Definition==//
//--alphabet_number--//
#ifndef CLE_NUM
#define CLE_NUM 18
#endif
#ifndef AMI_NUM
#define AMI_NUM 21
#endif

//--average_CLE_boundaryr--//
#ifndef MAX_AVECLE
#define MAX_AVECLE 10
#endif
#ifndef MIN_AVECLE
#define MIN_AVECLE 0
#endif
//==Definition==//over


//-----extern data_structure ------//
extern int Ori_CLESUM[CLE_NUM][CLE_NUM];
extern int Ori_BLOSUM[AMI_NUM][AMI_NUM];
extern int Blo_AA_Map[AMI_NUM];
extern int Ori_AA_Map[26];
extern int Ori_Anti_Map[AMI_NUM];
extern const char Ami_Clus[AMI_NUM+1];


//====class: SequenceScorer====//
//=> Conformational_Letter & Amino_Acid string process
class SequenceScorer
{
public:
	SequenceScorer(int version=1);
	~SequenceScorer(void);

//___data_structure__//
public:
	double AMI_Weight;
	double CLE_Weight;
	int **Gen_CLESUM;        // Gen_CLESUM
//[1]Glo_CLESUM_realted//__init_part__//
public:
	void Gen_CLESUM_Create(int Version=1);           // Create_GLO_Matrix (default:ori_version)
	void Gen_CLESUM_Create_Ori(void);                // Ori_Version //
	void Gen_CLESUM_Create_Glo(void);                // Glo_Version //


//==============Function============//(related)
public:
	//transform
	void AMI_transform(const char *AMI,int *out);
	void CLE_transform(const char *CLE,int *out);
	void AMI_CLE_transform_Ori(const char *AMI,const char *CLE,int *out);  //AMI+CKE->OUT
	void AMI_CLE_transform(const char *AMI,const char *CLE,int *GEN);      //AMI+CLE->GEN

	//---------------- Normal_SFP (fixed length) --------------//
	//calc universal
	double Universal_Calc(int ii,int jj,int len,
		int *IN1,int *IN2,int mollen1,int mollen2,int InType); // 1:CLE, 2:GEN, 3:AA, 4:AA+CLE, 5:AA+3Di, 6:3Di
	//check universal
	int Universal_Check(int *IN,int start,int len,int moln,int InType);      //check broken chain universal
	//---------------- Seed_SFP (flexible length) --------------//__090420__//
	//seed explosion
	int Seed_Explosion(int temp[4],int ic,int limi,int &bk,int &fw,
		int *GEN1,int *IN2,int mollen1,int mollen2,int InType);
	void Positive_Kmer_SFP_Seed(vector <SFP_Record> &SFP_List,int core_length,
		int *IN1,int *IN2,int mollen1,int mollen2,int isChange,int InType,int &num);
	void Universal_SFP_Seed_II(vector <SFP_Record> &SFP_List1,vector <SFP_Record> &SFP_List2,
		int winlen1,int winlen2,int ic1,int ic2,int extend1,int extend2,int step1,int step2,
		int *IN1,int *IN2,int mollen1,int mollen2,int isChange,int InType,int &num1,int &num2);
};
//====class_StruCode====//over

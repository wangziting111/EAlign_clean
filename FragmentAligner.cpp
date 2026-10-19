#include "FragmentAligner.h"
#include "SpatialAlignmentProposal.h"
#include <cmath>
#include <queue>
#include <stdexcept>

// 0: rank R0 by spatial support; 1: Top-6 weighted Ali3-Kabsch;
// 2: 1-NN ordinary Kabsch; 3: rank R0 by evolutionary seed score only.
#ifndef ALI3_POSE_REFINEMENT_MODE
#define ALI3_POSE_REFINEMENT_MODE 1
#endif

#if ALI3_POSE_REFINEMENT_MODE < 0 || ALI3_POSE_REFINEMENT_MODE > 3
#error "ALI3_POSE_REFINEMENT_MODE must be 0, 1, 2, or 3"
#endif

namespace
{
double ScoreTransformedQueryAgainstTarget(const XYZ *query, int query_length,
	const XYZ *target, int target_length, const double *query_to_target,
	double threshold, double d0)
{
	if (query == 0 || target == 0 || query_to_target == 0 ||
		query_length <= 0 || target_length <= 0 || threshold <= 0.0 || d0 <= 0.0) {
		return -1.0;
	}
	const double threshold_square = threshold * threshold;
	const double d0_square = d0 * d0;
	double score = 0.0;
	for (int i = 0; i < query_length; ++i) {
		const XYZ &source = query[i];
		XYZ transformed;
		transformed.X = source.X * query_to_target[0] +
			source.Y * query_to_target[1] + source.Z * query_to_target[2] +
			query_to_target[9];
		transformed.Y = source.X * query_to_target[3] +
			source.Y * query_to_target[4] + source.Z * query_to_target[5] +
			query_to_target[10];
		transformed.Z = source.X * query_to_target[6] +
			source.Y * query_to_target[7] + source.Z * query_to_target[8] +
			query_to_target[11];

		double best_distance_square = threshold_square;
		bool found = false;
		for (int j = 0; j < target_length; ++j) {
			const double distance_square = target[j].distance_square(transformed);
			if (distance_square <= threshold_square &&
				(!found || distance_square < best_distance_square)) {
				best_distance_square = distance_square;
				found = true;
			}
		}
		if (found) {
			score += 1.0 / (1.0 + best_distance_square / d0_square);
		}
	}
	return score;
}
}
//------ constructor -------//
FragmentAligner::FragmentAligner(int num)
:AlignmentGeometry(num),FragmentChainer(num),TMAligner(num),
 AlignmentRefiner(num)
{
	//init
	Fragment_Maximal=num;
	Fragment_Init(num);
	CUR_MaxJ_thres=0;
	//parameters
	ANG_CUT=60;         //
	Tail_Check=0.9;     //whether do tail check
	Comp_Check=0.9;     //whether do comparison check
	Comp_Check_Vari=0;  //comparison check (partial count)
	Geo_Sparse_DP=1;    //default: use Geo-distance in sparse DP
	SparseDP_filter=0; //default: use full Deepscore DP after Ali3-kabsch
}
FragmentAligner::~FragmentAligner(void)
{
	//dele
	Fragment_Dele();
}

//=================== process functions ============//
//--------- init ---------//
void FragmentAligner::Fragment_Init(int maxlen)
{
	//-- ali3_init --//
	ali1=new int[maxlen];
	ali2=new int[maxlen];
	tali1=new int[maxlen];
	tali2=new int[maxlen];
	//-- head_index --//
	valid_idx=new short[maxlen*maxlen];
	center_mat=new double[12*maxlen*10];
	head_score=new int[maxlen*10];
	head_index=new int[maxlen*10];
	//-- AFP & ali ---//
	Fragment_Cor=new int[4*maxlen];
	Fragment_Cor_temp=new int[4*maxlen];
	//-- pipeline_rotmat --//
	pipeline_rotmat=new double[12];
}
void FragmentAligner::Fragment_Dele(void)
{
	//-- ali3_dele --//
	delete [] ali1;
	delete [] ali2;
	delete [] tali1;
	delete [] tali2;
	//-- head_index --//
	delete [] center_mat;
	delete [] head_score;
	delete [] head_index;
	//-- AFP & ali ---//
	delete [] Fragment_Cor;
	delete [] Fragment_Cor_temp;
	//-- pipeline_rotmat --//
	delete [] pipeline_rotmat;
}

//--------- vice process ---------//
//calculate simple TMscore
//----------- Elongation and Shrinking -----------//
void FragmentAligner::Elongation(double *rotmat_,double thres,int cutoff,int *ali1,int *ali2,
						 XYZ *mol1,XYZ *mol2,int moln1,int moln2)
{
	int i;
	int ii,jj;
	int index1,index2;
	double Dou_Thres=thres*thres;
	int success;
	XYZ temp1,temp2,temp;
	int cut1,cut2,cut;


    //forward_elongation
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]<0)continue;
		//forward
		index1 = ali2[i];
		index2 = i;

#ifdef DEBUG
	overrange_debug_test(index1,moln1);
	overrange_debug_test(index2,moln2);
	lessrange_debug_test(index1,0);
	lessrange_debug_test(index2,0);
#endif

		for(;;)
		{
			success=0; //default:fail
			index1++;
			index2++;
			if(index1>moln1-1 || index2>moln2-1) break;
			jj=ali1[index1];
			ii=ali2[index2];
			if(jj!=-1 && ii!=-1) break;
			else if(jj!=-1 && ii==-1)
			{
				Ali_Cache_Point(mol1,index1,temp,rotmat_);
				cut1=Refine_Vector_Score_Forward(index1,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				cut2=Refine_Vector_Score_Forward(index1,jj,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				if((mol2[index2].distance_square(temp) < mol2[jj].distance_square(temp)) || cut1>cut2)
				{
					ali2[jj]=-1;
					ali1[index1]=index2;
					ali2[index2]=index1;
					success=1; //now success
				}
			}
			else if(jj==-1 && ii!=-1)
			{
				Ali_Cache_Point(mol1,index1,temp1,rotmat_);
				Ali_Cache_Point(mol1,ii,temp2,rotmat_);
				cut1=Refine_Vector_Score_Forward(index1,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				cut2=Refine_Vector_Score_Forward(ii,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				if((mol2[index2].distance_square(temp1) < mol2[index2].distance_square(temp2)) || cut1>cut2)
				{
					ali1[ii]=-1;
					ali1[index1]=index2;
					ali2[index2]=index1;
					success=1; //now success
				}
			}
			else
			{
				Ali_Cache_Point(mol1,index1,temp,rotmat_);
				cut=Refine_Vector_Score_Forward(index1,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				if(mol2[index2].distance_square(temp) < Dou_Thres || cut>cutoff)// distance² < FIN_CUT²
				{
					ali1[index1]=index2;
					ali2[index2]=index1;
					success=1; //now success
				}
			}

			if(success==0)break; //failed
			i++;
		}//end of FOR(;;)
	}//end of FOR(i)


	//backward_elongation
	for(i=moln2-1;i>=0;i--)
	{
		if(ali2[i]<0)continue;
		//forward
		index1 = ali2[i];
		index2 = i;

#ifdef DEBUG
	overrange_debug_test(index1,moln1);
	overrange_debug_test(index2,moln2);
	lessrange_debug_test(index1,0);
	lessrange_debug_test(index2,0);
#endif

		for(;;)
		{
			success=0; //default:fail
			index1--;
			index2--;
			if(index1<0 || index2<0) break;
			jj=ali1[index1];
			ii=ali2[index2];
			if(jj!=-1 && ii!=-1) break;
			else if(jj!=-1 && ii==-1)
			{
				Ali_Cache_Point(mol1,index1,temp,rotmat_);
				cut1=Refine_Vector_Score_Backward(index1,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				cut2=Refine_Vector_Score_Backward(index1,jj,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				if((mol2[index2].distance_square(temp) < mol2[jj].distance_square(temp)) || cut1>cut2)
				{
					ali2[jj]=-1;
					ali1[index1]=index2;
					ali2[index2]=index1;
					success=1; //now success
				}
			}
			else if(jj==-1 && ii!=-1)
			{
				Ali_Cache_Point(mol1,index1,temp1,rotmat_);
				Ali_Cache_Point(mol1,ii,temp2,rotmat_);
				cut1=Refine_Vector_Score_Backward(index1,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				cut2=Refine_Vector_Score_Backward(ii,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				if((mol2[index2].distance_square(temp1) < mol2[index2].distance_square(temp2)) || cut1>cut2)
				{
					ali1[ii]=-1;
					ali1[index1]=index2;
					ali2[index2]=index1;
					success=1; //now success
				}
			}
			else
			{
				Ali_Cache_Point(mol1,index1,temp,rotmat_);
				cut=Refine_Vector_Score_Backward(index1,index2,rotmat_,mol1,mol2,moln1,moln2,SCALE);
				if(mol2[index2].distance_square(temp) < Dou_Thres || cut>cutoff)
				{
					ali1[index1]=index2;
					ali2[index2]=index1;
					success=1; //now success
				}
			}

			if(success==0)break; //failed
			i--;
		}//end of FOR(;;)
	}//end of FOR(i)
}
void FragmentAligner::Shrinking(double *rotmat_,double thres,int cutoff,int *ali1,int *ali2,
						XYZ *mol1,XYZ *mol2,int moln1,int moln2)
{
	int i;
	int ii,jj;
	int cut;
	double dist;
	XYZ temp;
	double Dou_Thres=thres*thres;
	double SQURE_DIST=MAX_DIST*MAX_DIST;

	for(i=0;i<moln2;i++)
	{
		ii=ali2[i];
		jj=i;
		if(ii==-1)continue;

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

		Ali_Cache_Point(mol1,ii,temp,rotmat_);
		dist=mol2[jj].distance_square(temp);
		if(dist>SQURE_DIST)
		{
			ali2[jj]=-1;
			ali1[ii]=-1;
		}
		else if(dist>Dou_Thres)
		{
			cut=Refine_Vector_Score_Forward(ii,jj,rotmat_,mol1,mol2,moln1,moln2,SCALE);
			if(cut<cutoff)
			{
				ali2[jj]=-1;
				ali1[ii]=-1;
			}
		}
	}
}
//------------- Select Best Pivot --------------//
//given a certain SFP, return its' TMScorer against other SFPs
//----------- AlignmentMetricAccuracy ----------------//
double FragmentAligner::Alg_Comp_Simp(vector <int> &ali_m1,vector <int> &ali_m2,int part) //using ali2
{
	int i;
	int ii1,ii2;
	int count=0;
	int parti=0;
	int total=0;
	for(i=0;i<moln2;i++)
	{
		ii1=ali_m1[i];
		ii2=ali_m2[i];
		if(ii1==-1 || ii2==-1)continue;
		if(ii1==ii2)count++;
		else
		{
			if(abs(ii1-ii2)<=part)parti++;
		}
		total++;
	}
	if(total==0)return 0.0;
	else return 1.0*(count+parti)/total;
}

//==================== virtual function ====================//-> part[0]: io
void FragmentAligner::Pipeline_Input(XYZ *m1,XYZ *m2,int n1,int n2)
{
	mol1=m1;
	mol2=m2;
	moln1=n1;
	moln2=n2;
	Other_Fragment_Init();
}
//==================== virtual function ====================//-> part[1]: score
//given alignment and rotmat, return scores (maybe TMscore,RMSD and LALI, etc)
double FragmentAligner::Pipeline_Make_Score(Align_Record & align_record)
{
	int i;
	int ii,jj;
	double dist2;
	double ori_d;
	double sco;
	double rms;
	double tms;
	int count;
	int lali;
	double rmsd;
	//init equal
	align_record.alignment.resize(moln2);
	for(i=0;i<moln2;i++)align_record.alignment[i]=ali2[i];
	align_record.rotmat.resize(12);
	for(i=0;i<12;i++)align_record.rotmat[i]=pipeline_rotmat[i];
	//calculate
	double d_0=TM_d0;
	XYZ temp=0.0;
	ori_d=d_0*d_0;
	sco=0.0;
	rms=0.0;
	count=0;
	for(i=0;i<moln2;i++)
	{
		ii=ali2[i];
		jj=i;
		if(ii==-1)continue;

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

		rot_point(mol1[ii],temp,pipeline_rotmat);
		dist2=mol2[jj].distance_square(temp);
		sco+=1.0/(1.0+dist2/ori_d);
		rms+=dist2;
		count++;
	}
	//evaluate
	rmsd=sqrt(1.0*rms/count);
	lali=count;
	tms=1.0*sco/TM_smaller;
	//final
	align_record.RMSD=rmsd;
	align_record.lali=lali;
	align_record.TMsco=tms;
	align_record.main_sco=tms;
	return tms;
}
//given the align_record, make the secondary score (maybe LALI*TMscore, etc)
double FragmentAligner::Pipeline_Get_Score(Align_Record &align)
{
	double tms=align.main_sco;
	return tms;
}

//==================== virtual function ====================//-> part[2]: partial
//0. Init Other FragmentAligner
void FragmentAligner::Other_Fragment_Init(void)
{
	//STEP[0] initialization
	int mol_len=moln1<=moln2?moln1:moln2;   //smaller_length
	double d_0=Calc_TM_d0_Simp(mol_len);
	TM_smaller=mol_len;
	TM_d0=d_0;
	//self_adaptive[1]->RMSD_CUTOFF
	INI_CUT=3*d_0;
	INI_CUT=INI_CUT<MAX_DIST?INI_CUT:MAX_DIST;
	INI_CUT=INI_CUT>MIN_DIST?INI_CUT:MIN_DIST;
	FIN_CUT=d_0;
	FIN_CUT=FIN_CUT<MAX_DIST?FIN_CUT:MAX_DIST;
	FIN_CUT=FIN_CUT>MIN_DIST?FIN_CUT:MIN_DIST;
	//init Maximal
	CUR_MaxJ=0.0;
	CUR_MaxK=0.0;
}
//1. ZoomIn Add (ret_rot might be updated)
int FragmentAligner::ZoomIn_Add(int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP,double *ret_rot)
{
	int i,j;
	int ii,jj;
	int winlen;
	int QQ_temp_num;
	int QQ_cur;
	int Mx_K;
	double ret_val;
	double thres,minus;
	int neo_ii,neo_jj,neo_len;
	double score;
	double RMSD;
	double d_0=d0;
	//tot_init
	QQ_cur=0;
	QQ_temp_num=(int)SFP.size();
	memset(valid_idx,1,QQ_temp_num*sizeof(short));
	//rigid_frag_adding
	minus=(INI_CUT-FIN_CUT)/recur;
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	for(i=0;i<recur;i++)
	{
		//init
		Refine_Init(ali3,moln2);
		if(i<recur-1)Mx_K=(int)(QQ_temp_num/pow(2.0,i+1));
		else Mx_K=(int)(QQ_temp_num/pow(2.0,i));
		thres=INI_CUT-(i+1)*minus;
		//check_frag
		for(j=0;j<QQ_cur+Mx_K;j++)
		{
			if(valid_idx[j]==-1)continue;  // -1 means killed
			ii=SFP[j].ii;
			jj=SFP[j].jj;
			winlen=SFP[j].winlen;
			calc_frag_dist_thres(ii,jj,winlen,ret_rot,mol1,mol2,moln1,moln2,thres,neo_ii,neo_jj,neo_len,score);//剪切出满足距离阈值的最长连续子片段:neo_ii,最长有效子片段在 mol1 中的起点;neo_jj,最长有效子片段在 mol2 中的起点;neo_len,长度
			if(neo_len==0)valid_idx[j]=-1;
			else Refine_Add(neo_ii,neo_jj,neo_len,ret_rot,mol1,mol2,moln1,moln2,ali3);
		}
		QQ_cur+=Mx_K;
		//rot_mol
		ret_val=ali3_kabsch(ret_rot,RMSD,mol1,mol2,moln1,moln2,ali3); //apply ALI3(weight) rotation,基于ali3-pairs更新ret_rot
		if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
		if(ret_val<0.0 || RMSD<0.0)return -1; //incorrect
	}
	//Refine_DynaProg
	Refine_DynaProg(ret_rot,0.0,d_0,mol1,mol2,moln1,moln2,FOR_GAP,BAK_GAP,SCALE,ali3);
	return (int)Refine_TraceBack(ret_rot,d_0,mol1,mol2,moln1,moln2,ali3,ali1,ali2);
}
//2 .Refinement (ret_rot might be updated)
int FragmentAligner::Refinement(int recur,double FIN_CUT,double *ret_rot)
{
	int k;
	double ret_val;
	int score,bak_score;
	int success;
	int count;
	double rot_bak[12];
	double RMSD;
	int Range=4;
	double DP_BETA=0.5;//考虑了forward vector dot similarity
	double d_0=d0;

	//process
	count=0;
	success=1; //default:success
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	bak_score=(int)Refine_TraceBack(ret_rot,d_0,mol1,mol2,moln1,moln2,ali3,ali1,ali2); // TM-kernel 累积值
	score=bak_score;
	for(k=0;k<recur;k++)
	{
		//Back
		if(success==1)
		{
			memcpy(rot_bak,ret_rot,12*sizeof(double));
			memcpy(tali1,ali1,moln1*sizeof(int));
			memcpy(tali2,ali2,moln2*sizeof(int));
		}
		if(count==2)return score; //[refine_better]//__081130__//
		//Elongation & Collect
		// 从当前已对齐 residue pair 出发，同时向前和向后沿*对角线*尝试扩展,条件为distance² < FIN_CUT²
		Elongation(ret_rot,FIN_CUT,INT_MAX_NUM,ali1,ali2,mol1,mol2,moln1,moln2);
		Ali_To_Cor(Fragment_Cor,1,moln1,moln2,ali1,ali2);//阈值为 1，所以任何长度至少为 1 的连续匹配都会保留
		//Rotation	
		Refine_Cor_Add(Fragment_Cor,0,ret_rot,mol1,mol2,moln1,moln2,ali3);
		ret_val=ali3_kabsch(ret_rot,RMSD,mol1,mol2,moln1,moln2,ali3); //apply ALI3(weight) rotation: 更新 pose
		if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
		if(ret_val<0.0 || RMSD<0.0)return -1; //incorrect
		//Replace
		Refine_Cor_Add(Fragment_Cor,Range,ret_rot,mol1,mol2,moln1,moln2,ali3);//Range=4,允许当前 alignment 在局部范围内发生约 ±4 residue 的 register correction
		Refine_DynaProg(ret_rot,DP_BETA,d_0,mol1,mol2,moln1,moln2,FOR_GAP,BAK_GAP,SCALE,ali3);//DP_BETA=0.5,s(i,j)=0.5sdistance​+0.5sbackbone-direction​
		score=(int)Refine_TraceBack(ret_rot,d_0,mol1,mol2,moln1,moln2,ali3,ali1,ali2);
		//Check
		if(score>bak_score)
		{
			bak_score=score;
			success=1;
			count=0;
		}
		else
		{
			if(score<bak_score*0.95)
			{
				memcpy(ret_rot,rot_bak,12*sizeof(double));
				memcpy(ali1,tali1,moln1*sizeof(int));
				memcpy(ali2,tali2,moln2*sizeof(int));
				if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
				return (bak_score);
			}
			success=0;
			count++;
		}
	}
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	return score; //[refine_better]//__081130__//
}
//3. Kill NonLinear (ret_rot won't be updated)
int FragmentAligner::Kill_NonLinear(double *ret_rot)
{
	Fragment_Kill_NonLinear(Fragment_Cor,Fragment_Cor_temp,ret_rot);
	return 1;
}
//4. Final Stage (combine Shrink & Elong)
int FragmentAligner::Final_Stage(double *ret_rot)
{
	//--Final_Statis--//
	double RMSD;
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	Elongation(ret_rot,FIN_CUT,INT_MAX_NUM,ali1,ali2,mol1,mol2,moln1,moln2);
	Shrinking(ret_rot,FIN_CUT,ANG_CUT,ali1,ali2,mol1,mol2,moln1,moln2);
	Final_Rot_Ali(ret_rot,RMSD,mol1,mol2,moln1,moln2,ali1,ali2);
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	Elongation(ret_rot,FIN_CUT,ANG_CUT,ali1,ali2,mol1,mol2,moln1,moln2);
	return 1;
}
//5. TMscore refinement (ret_rot might be updated)
int FragmentAligner::Final_Refine(double *ret_rot)
{
	Calc_TM_Align(mol1,mol2,moln1,moln2,ali2,ali2);    //apply TMscore dynamic programming
	TM_Align_TM_Score_Simp(mol1,mol2,moln1,moln2,ali2,RMSD,LALI); //use simple version
	EqualArray(ret_rot,finmat,12);
	return 1;
}

//==================== virtual function ====================//-> part[3]: global
//0. Init functions
void FragmentAligner::Global_Pipeline_Init(void)
{
	int i;
	align_rec.resize(moln1*moln2);
	for(i=0;i<moln1*moln2;i++)align_rec[i]=0;
}
//1. Select best pivot
int FragmentAligner::Select_Best_Pivot(int recur,double INI_CUT,double FIN_CUT,vector <SFP_Record> &SFP,int TopJ)
{
	int i;
	int relnum=(int)SFP.size();
	if(relnum==0)return 0;
	int totnum;
	totnum=TopJ;
	totnum=totnum<Fragment_Maximal?totnum:Fragment_Maximal;
	totnum=totnum<relnum?totnum:relnum;
	if (INI_CUT != FIN_CUT) {
    	throw std::runtime_error("INI_CUT must equal FIN_CUT");
	}
	const double threshold=FIN_CUT;
	const double score_scale=1000.0;
	(void)recur;
	(void)INI_CUT;
	//real_start
	for(i=0;i<totnum;i++)
	{
		const SFP_Record &pivot=SFP[i];
		double *query_to_target=center_mat+12*i;
		if(make_center(pivot.ii,pivot.jj,pivot.winlen,query_to_target,
			mol1,mol2,moln1,moln2)<0.0)
		{
			head_score[i]=-1;
			continue;
		}

		const double score=ScoreTransformedQueryAgainstTarget(
			mol1,moln1,mol2,moln2,query_to_target,threshold,d0);
		head_score[i]=(int)(score_scale*score);
	}
	//fast_sort
	fast_sort.fast_sort_1(head_score,head_index,totnum);
	int smaller=moln1<moln2?moln1:moln2;	
	CUR_MaxJ=1.0*head_score[head_index[0]]/(score_scale*smaller);
	//fast_check
	if(FAST_Chk==1)
	{
//		printf("%d\n",FAST_Chk);
//		printf("%lf\n",CUR_MaxJ);
		if(CUR_MaxJ<CUR_MaxJ_thres)return 0;
	}
	return totnum;
}

namespace
{
	const int CE_FRAGMENT_LENGTH=8;
	const int CE_INTERNAL_DISTANCE_COUNT=
		(CE_FRAGMENT_LENGTH-1)*(CE_FRAGMENT_LENGTH-2)/2;
	const double CE_Fragment_DISTANCE_CUTOFF=3.0;

	struct CEFragmentDescriptor
	{
		int start;
		double distances[CE_INTERNAL_DISTANCE_COUNT];
	};

	struct CEAFP
	{
		int mol1_start;
		int mol2_start;
		double distance;
	};

	struct CEAFPBetter
	{
		bool operator ()(const CEAFP &left,const CEAFP &right) const
		{
			if(left.distance!=right.distance)
				return left.distance<right.distance;
			if(left.mol1_start!=right.mol1_start)
				return left.mol1_start<right.mol1_start;
			return left.mol2_start<right.mol2_start;
		}
	};

	struct GeometryPoseSeed
	{
		int mol1_start;
		int mol2_start;
		double support_score;
		double rotmat[12];
	};

	bool GeometryPoseScoreGreater(const GeometryPoseSeed &left,
		const GeometryPoseSeed &right)
	{
		if(left.support_score!=right.support_score)
			return left.support_score>right.support_score;
		if(left.mol1_start!=right.mol1_start)
			return left.mol1_start<right.mol1_start;
		return left.mol2_start<right.mol2_start;
	}

	struct SpatialPoseCandidate
	{
		explicit SpatialPoseCandidate(int candidate_count)
		:pivot_index(-1),mol1_start(-1),mol2_start(-1),proposal_score(-1.0),
		 proposal(candidate_count) {}
		int pivot_index;
		int mol1_start;
		int mol2_start;
		double proposal_score;
		double rotmat[12];
		SpatialAlignmentProposal proposal;//!
	};

	bool SpatialPoseScoreGreater(const SpatialPoseCandidate &left,
		const SpatialPoseCandidate &right)
	{
		if(left.proposal_score!=right.proposal_score)
			return left.proposal_score>right.proposal_score;
		if(left.pivot_index!=right.pivot_index)
			return left.pivot_index<right.pivot_index;
		if(left.mol1_start!=right.mol1_start)
			return left.mol1_start<right.mol1_start;
		return left.mol2_start<right.mol2_start;
	}

	double PointDistance(const XYZ &left,const XYZ &right)
	{
		return std::sqrt(left.distance_square(right));
	}

	vector <CEFragmentDescriptor> BuildCEFragmentDescriptors(
		const XYZ *mol,int moln)
	{
		vector <CEFragmentDescriptor> descriptors;
		if(mol==NULL || moln<CE_FRAGMENT_LENGTH)return descriptors;
		descriptors.reserve(moln-CE_FRAGMENT_LENGTH+1);
		for(int start=0;start+CE_FRAGMENT_LENGTH<=moln;start++)
		{
			CEFragmentDescriptor descriptor;
			descriptor.start=start;
			int distance_index=0;
			for(int left=0;left<CE_FRAGMENT_LENGTH;left++)
				for(int right=left+2;right<CE_FRAGMENT_LENGTH;right++)
					descriptor.distances[distance_index++]=PointDistance(
						mol[start+left],mol[start+right]);
			descriptors.push_back(descriptor);
		}
		return descriptors;
	}

	vector <CEAFP> SelectTopCEAFPs(const XYZ *mol1,int moln1,
		const XYZ *mol2,int moln2,int max_pairs)
	{
		vector <CEAFP> selected;
		if(max_pairs<=0)return selected;
		const vector <CEFragmentDescriptor> first=
			BuildCEFragmentDescriptors(mol1,moln1);
		const vector <CEFragmentDescriptor> second=
			BuildCEFragmentDescriptors(mol2,moln2);
		CEAFPBetter better;
		priority_queue <CEAFP,vector <CEAFP>,CEAFPBetter> best_afps;
		const double normalization=2.0/
			(CE_FRAGMENT_LENGTH*CE_FRAGMENT_LENGTH);
		for(size_t i=0;i<first.size();i++)
		{
			for(size_t j=0;j<second.size();j++)
			{
				double distance_sum=0.0;
				for(int k=0;k<CE_INTERNAL_DISTANCE_COUNT;k++)
					distance_sum+=std::fabs(
						first[i].distances[k]-second[j].distances[k]);
				const double ce_distance=normalization*distance_sum;
				if(ce_distance>=CE_Fragment_DISTANCE_CUTOFF)continue;
				const CEAFP afp={first[i].start,second[j].start,ce_distance};
				if((int)best_afps.size()<max_pairs)best_afps.push(afp);
				else if(better(afp,best_afps.top()))
				{
					best_afps.pop();
					best_afps.push(afp);
				}
			}
		}
		selected.reserve(best_afps.size());
		while(!best_afps.empty())
		{
			selected.push_back(best_afps.top());
			best_afps.pop();
		}
		stable_sort(selected.begin(),selected.end(),better);
		return selected;
	}
	// 给定pose,对mol2每个residue找到mol1中nearest residue,计算距离累加
	double ScoreGeometryPose(const XYZ *mol1,int moln1,const XYZ *mol2,int moln2,
		const double *rotmat,double threshold,double d0)
	{
		if(mol1==NULL || mol2==NULL || rotmat==NULL || moln1<=0 || moln2<=0 ||
			threshold<=0.0 || d0<=0.0)return -1.0;
		const double threshold_square=threshold*threshold;
		const double d0_square=d0*d0;
		vector <XYZ> transformed_mol1(moln1);
		for(int i=0;i<moln1;i++)
		{
			const XYZ &source=mol1[i];
			XYZ &target=transformed_mol1[i];
			target.X=source.X*rotmat[0]+source.Y*rotmat[1]+
				source.Z*rotmat[2]+rotmat[9];
			target.Y=source.X*rotmat[3]+source.Y*rotmat[4]+
				source.Z*rotmat[5]+rotmat[10];
			target.Z=source.X*rotmat[6]+source.Y*rotmat[7]+
				source.Z*rotmat[8]+rotmat[11];
		}
		double score=0.0;
		for(int j=0;j<moln2;j++)
		{
			double best_distance_square=threshold_square;
			bool found=false;
			for(int i=0;i<moln1;i++)
			{
				const double distance_square=
					mol2[j].distance_square(transformed_mol1[i]);
				if(distance_square<=threshold_square &&
					(!found || distance_square<best_distance_square))
				{
					best_distance_square=distance_square;
					found=true;
				}
			}
			if(found)score+=1.0/(1.0+best_distance_square/d0_square);
		}
		return score;
	}

	void BuildSparseCombinedScores(TMAligner &scorer,
		const SpatialAlignmentProposal &proposal,double *rotmat,
		int moln1,int moln2,double d0,vector <double> &scores)
	{
		const int capacity=proposal.CandidateCapacity();
		scores.assign(moln2*capacity,0.0);
		for(int j=0;j<moln2;j++)
		{
			for(int k=0;k<proposal.CandidateCount(j);k++)
			{
				const SpatialNeighbor &candidate=proposal.Candidate(j,k);
				scores[j*capacity+k]=scorer.TM_Align_Get_Pair_Score(
					candidate.index,j,candidate.distance_square,d0,rotmat,
					moln1,moln2);
			}
		}
	}
}

double FragmentAligner::Pipeline_Full(XYZ *m1,XYZ *m2,int n1,int n2,
	vector <SFP_Record> &SFP_High,vector <SFP_Record> &SFP_Low,
	vector <Align_Record> &align_records)
{
	(void)SFP_Low;
	const int smaller=n1<n2?n1:n2;
	if(smaller==0)return -1.0;
	const int top_l_limit=32;
	const int top_k_limit=10;
	const double radius=INI_CUT>0.0?INI_CUT:5.0;
#if ALI3_POSE_REFINEMENT_MODE == 1
	const int ali3_size=n2*(ali3_TOT+1)*2;
#endif
	align_records.clear();
	pipeline_MAX.main_sco=0.0;

	// Evolution branch: High-SFP rotations followed by staged Ali3 refinement.
	Pipeline_Input(m1,m2,n1,n2);
	Other_Fragment_Init();
	TM_Precompute_Vector_Cache(m1,m2,n1,n2);
#if ALI3_POSE_REFINEMENT_MODE == 3
	// SFP_High is already stably sorted by descending 3Di seed score.
	const int available=(int)SFP_High.size();
#else
	const int available=Select_Best_Pivot(1,INI_CUT,INI_CUT,SFP_High,ZM_TopJ);
#endif
	Global_Pipeline_Init();
	const int pose_shortlist_limit=
		(ALI3_POSE_REFINEMENT_MODE==0 || ALI3_POSE_REFINEMENT_MODE==3)?
		top_k_limit:top_l_limit;
	const int evolution_top_l=available<pose_shortlist_limit?
		available:pose_shortlist_limit;
	vector <SpatialPoseCandidate> evolution_poses;
	evolution_poses.reserve(evolution_top_l);
	for(int rank=0;rank<evolution_top_l;rank++)
	{
#if ALI3_POSE_REFINEMENT_MODE == 3
		const int pivot_index=rank;
#else
		const int pivot_index=head_index[rank];
		if(head_score[pivot_index]<0)continue;
#endif
		const int proposal_capacity=ALI3_POSE_REFINEMENT_MODE==2?1:ali3_TOT;
		SpatialPoseCandidate candidate(proposal_capacity);
		candidate.pivot_index=pivot_index;

#if ALI3_POSE_REFINEMENT_MODE == 3
		const SFP_Record &seed=SFP_High[pivot_index];
		if(make_center(seed.ii,seed.jj,seed.winlen,candidate.rotmat,
			m1,m2,n1,n2)<0.0)continue;
		candidate.proposal_score=seed.score;
#else
		EqualArray(candidate.rotmat,center_mat+12*pivot_index,12);

#if ALI3_POSE_REFINEMENT_MODE == 0
		// Select_Best_Pivot already computed the R0 nearest-neighbor support.
		candidate.proposal_score=head_score[pivot_index]/1000.0;
#else
		if(!candidate.proposal.BuildFromCoordinates(m1,n1,m2,n2,
			candidate.rotmat,radius))continue;
		double proposal_rmsd=0.0;

#if ALI3_POSE_REFINEMENT_MODE == 1
		vector <int> proposal_ali3(ali3_size);
		if(!candidate.proposal.WriteLegacyAli3(&proposal_ali3[0],ali3_TOT))continue;
		if(ali3_kabsch(candidate.rotmat,proposal_rmsd,m1,m2,n1,n2,
			&proposal_ali3[0])<0.0 || proposal_rmsd<0.0)continue;
#else
		// Hard assignment: one nearest mol1 residue for each mol2 residue.
		vector <XYZ> nearest_mol1;
		vector <XYZ> nearest_mol2;
		nearest_mol1.reserve(n2);
		nearest_mol2.reserve(n2);
		for(int j=0;j<n2;j++)
		{
			if(candidate.proposal.CandidateCount(j)==0)continue;
			const int i=candidate.proposal.Candidate(j,0).index;
			nearest_mol1.push_back(m1[i]);
			nearest_mol2.push_back(m2[j]);
		}
		const int correspondence_count=(int)nearest_mol1.size();
		if(correspondence_count<pipeline_MIN_LEN)continue;
		proposal_rmsd=kabsch(&nearest_mol2[0],&nearest_mol1[0],
			correspondence_count,candidate.rotmat);
		if(proposal_rmsd<0.0)continue;
#endif

		if(!candidate.proposal.BuildFromCoordinates(m1,n1,m2,n2,
			candidate.rotmat,radius))continue;
		candidate.proposal_score=candidate.proposal.BestNeighborScore(d0);
#endif
#endif
		evolution_poses.push_back(candidate);
	}
	stable_sort(evolution_poses.begin(),evolution_poses.end(),SpatialPoseScoreGreater);
	const int evolution_selected=(int)evolution_poses.size()<top_k_limit?
		(int)evolution_poses.size():top_k_limit;
	for(int rank=0;rank<evolution_selected;rank++)
	{
		double sparse_score=0.0;
		vector <int> coarse_alignment;
		SpatialPoseCandidate &candidate=evolution_poses[rank];
		if(SparseDP_filter==1)
		{
			EqualArray(pipeline_rotmat,candidate.rotmat,12);
#if ALI3_POSE_REFINEMENT_MODE == 0 || ALI3_POSE_REFINEMENT_MODE == 3
			// Direct-R0 full-DP modes do not otherwise need a proposal.
			if(!candidate.proposal.BuildFromCoordinates(m1,n1,m2,n2,
				candidate.rotmat,radius))continue;
#endif
			if(Geo_Sparse_DP==1)
			{
				sparse_score=candidate.proposal.TraceOrderPreservingAlignment(
					d0,n1,coarse_alignment);
			}
			else
			{
				vector <double> sparse_combined_scores;
				BuildSparseCombinedScores(*this,candidate.proposal,candidate.rotmat,
					n1,n2,d0,sparse_combined_scores);
				sparse_score=candidate.proposal.TraceOrderPreservingAlignment(
					sparse_combined_scores,n1,coarse_alignment);
			}
		}
		else //default
		{
			//full dp + traceback
			EqualArray(pipeline_rotmat,candidate.rotmat,12);
			coarse_alignment.assign(n2,-1);
			rot_mol(m1,TM_tmp1,n1,pipeline_rotmat);
			TM_Align_Get_Ali_coarse_fullDP(TM_tmp1,m2,n1,n2,&coarse_alignment[0]);
			for(int j=0;j<n2;j++)
			{
				const int i=coarse_alignment[j];
				if(i>=0 && i<n1)sparse_score+=TM_DP_sco[i*n2+j];
			}
		}
		if(sparse_score<=0.0)continue;
		Ali_Init(n1,n2,ali1,ali2);
		for(int j=0;j<n2;j++)ali2[j]=coarse_alignment[j];
		Ali2_To_Ali1(n1,n2,ali1,ali2);
		if(Final_Rot_Ali(pipeline_rotmat,RMSD,m1,m2,n1,n2,ali1,ali2)<pipeline_MIN_LEN)
			continue;
		if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*n1);

		Align_Record record;
		Pipeline_Make_Score(record);
		record.index.clear();
		record.index.push_back(candidate.pivot_index);
		align_records.push_back(record);
	}


	if(align_records.empty())return -1.0;
	stable_sort(align_records.rbegin(),align_records.rend());
	pipeline_MAX=align_records[0];
	return pipeline_MAX.main_sco;
}
//2. Break Check
int FragmentAligner::Pipeline_Break_Check(int cur,int count,int totnum)
{
	if(count>=ZM_TopK)return 0;
	else return 1;
}
//3. Get Index
int FragmentAligner::Pipeline_Get_Index(int cur)
{
	CUR_INDEX=cur;
	INDEX_CUR=head_index[cur];
	return head_index[cur];
}
//4. Head Check
int FragmentAligner::Pipeline_Head_Check(SFP_Record &SFP)
{
	int ii=SFP.ii;
	int jj=SFP.jj;
	int winlen=SFP.winlen;

#ifdef DEBUG
	overrange_debug_test(ii+winlen-1,moln1);
	overrange_debug_test(jj+winlen-1,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

	int k;
	for(k=0;k<winlen;k++)
	{
		if(align_rec.at((ii+k)*moln2+(jj+k))==0)return 1; //un-possessed!!
	}
	return 0;
}
//5. Get initial rotmat
void FragmentAligner::Pipeline_Initial_Rotmat(int index,double *rotmat)
{
	EqualArray(rotmat,center_mat+12*index,12);
}
//6. Tail Check
int FragmentAligner::Pipeline_Tail_Check(Align_Record &cur)
{
	//init_check
	if(cur.score>pipeline_MAX.score)return 1;
	//process
	int i;
	int ii,jj;
	int count=0;
	int match=0;
	for(i=0;i<moln2;i++)
	{
		ii=cur.alignment[i];
		jj=i;
		if(ii==-1)continue;

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

		count++;
		if(align_rec.at((ii)*moln2+(jj))==1)match++; //possessed!!
	}
	if(count==0)return 0;
	if(1.0*match/count>Tail_Check)return 0;
	return 1;
}
//7. Comparison Check
int FragmentAligner::Pipeline_Comparison_Check(Align_Record &cur,vector <Align_Record> &tot)
{
	//init_check
	if(cur.score>pipeline_MAX.score)return 1;
	//process
	int i;
	double score;
	int totnum=(int)tot.size();
	for(i=0;i<totnum;i++)
	{
		score=Alg_Comp_Simp(cur.alignment,tot[i].alignment,Comp_Check_Vari);
		if(score>Comp_Check)return 0;
	}
	return 1;
}
//8. Terminal Process
void FragmentAligner::Pipeline_Terminal_Process(Align_Record &cur)
{
	int i;
	int ii,jj;
	for(i=0;i<moln2;i++)
	{
		ii=cur.alignment[i];
		jj=i;
		if(ii==-1)continue;

#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif

		align_rec.at((ii)*moln2+(jj))=1; //possess!!
	}
	//terminal_add
	cur.index.clear();
	cur.index.push_back(INDEX_CUR);
	cur.index.push_back(CUR_INDEX);
}

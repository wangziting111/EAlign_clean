#include "TMAligner.h"

//------ constructor -------//
TMAligner::TMAligner(int num)
:TMScorer(num)
{
	TM_maximal=num;
	InitializeTMAlignment(TM_maximal);
	//parameter
	TM_NOTM=0;       // 0 for NOT use NOTM (i.e., use TMscore to superimpose)
	TM_GAP_TYPE=0;   // 0 for ordinary; 1 for four-state; -1 for faster (dynaprog bound)//__110730__//
	TM_DIST_CUT=1;
	TM_DistCut=10.0; // equal to d8 in TMscore
	//DynaProg_bound//__110730__//
	TM_bound_neib=4;
	//TMali_Weight//__110830__//
	TM_Vect_Score=-1;  //vect score COMBINE !!
	TM_Wei_Score=0;    //wei score NO !!//?
	TMali_Weight=0;
	//TM_CB
	TM_cb1=0;
	TM_cb2=0;
	//TM_SCORE_TYPE
	TM_Score_Type=1;  //-> 0 for TMscore, 1 for CombinedScore
}
TMAligner::~TMAligner(void)
{
	ReleaseTMAlignment();
}

//--------- init ---------//
void TMAligner::InitializeTMAlignment(int maxlen)
{
	TM_DP_sco=new double[maxlen*maxlen];
	TM_ali1=new int[maxlen];
	TM_tmp1=new XYZ[maxlen];
	TM_tmp2=new XYZ[maxlen];
	TM_DP_ali2=new int[maxlen];
	TM_DP_best=new int[maxlen];
	TM_rotmat=new double[12];
	TM_vec_forw1=new XYZ[maxlen];
	TM_vec_forw2=new XYZ[maxlen];
	TM_vec_back1=new XYZ[maxlen];
	TM_vec_back2=new XYZ[maxlen];
	TM_vec_cb1=new XYZ[maxlen];
	TM_vec_cb2=new XYZ[maxlen];
	TM_vec_forw_valid1=new char[maxlen];
	TM_vec_forw_valid2=new char[maxlen];
	TM_vec_back_valid1=new char[maxlen];
	TM_vec_back_valid2=new char[maxlen];
	TM_vec_cb_valid1=new char[maxlen];
	TM_vec_cb_valid2=new char[maxlen];
}
void TMAligner::ReleaseTMAlignment(void)
{
	delete [] TM_DP_sco;
	delete [] TM_ali1;
	delete [] TM_tmp1;
	delete [] TM_tmp2;
	delete [] TM_DP_ali2;
	delete [] TM_DP_best;
	delete [] TM_rotmat;
	delete [] TM_vec_forw1;
	delete [] TM_vec_forw2;
	delete [] TM_vec_back1;
	delete [] TM_vec_back2;
	delete [] TM_vec_cb1;
	delete [] TM_vec_cb2;
	delete [] TM_vec_forw_valid1;
	delete [] TM_vec_forw_valid2;
	delete [] TM_vec_back_valid1;
	delete [] TM_vec_back_valid2;
	delete [] TM_vec_cb_valid1;
	delete [] TM_vec_cb_valid2;
}

//---------- primary functions -------------//
void TMAligner::TM_Align_Init(int moln1,int moln2)
{
	int smaller=moln1<moln2?moln1:moln2;
	Calc_TM_d0(smaller);
	TM_DistCut=d8;
	TM_bound.resize(moln1+1);
	TM_CALC=1;
}
int TMAligner::TM_Align_Get_XYZ(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
	int i;
	int ii,jj;
	int count=0;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]>=0)
		{
			ii=ali2[i];
			jj=i;
#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif	
			TM_tmp1[count]=mol1[ii];
			TM_tmp2[count]=mol2[jj];
			count++;
		}
	}
	return count;
}
int TMAligner::TM_Align_Get_CUT(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,double d8,double *rotmat)
{
	int i;
	int ii,jj;
	double dist2;
	double ori_d8=d8*d8;
	XYZ temp;
	int count=0;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]>=0)
		{
			ii=ali2[i];
			jj=i;
#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif	
			TMs_Cache_Point(mol1,ii,temp,rotmat);
			dist2=mol2[jj].distance_square(temp);
			if(dist2>ori_d8)
			{
				ali2[i]=-1;
				continue;
			}
			TM_tmp1[count]=mol1[ii];
			TM_tmp2[count]=mol2[jj];
			count++;
		}
	}
	return count;
}
double TMAligner::TM_Align_Get_Score_Simp(XYZ *mol1,XYZ *mol2,double *rotmat_,
	int moln1,int moln2,int *ali2)
{
	int i,k;
	double vec1,vec2,vec3;
	double vect;
	int tot_num;
	int range=1;
	double dist2;
	double ori_d=d0*d0;
	double tms;
	int pos;
	XYZ xyz;
	double ws_sco=0.0;
	//calc_score
	TM_Precompute_Vector_Cache(mol1,mol2,moln1,moln2);
	for(i=0;i<moln2;i++)
	{
		pos=ali2[i];
		if(pos<0 || pos>=moln1)continue;
		//single
		if(rotmat_!=0)TMs_Cache_Point(mol1,pos,xyz,rotmat_);
		else xyz=mol1[pos];
		dist2=mol2[i].distance_square(xyz);
		//nearby
		vec1=0.0;
		vec2=0.0;
		//calc
		tot_num=2;
		for(k=1;k<=range;k++)
		{
			vec1+=TM_Vector_Score_Forward_Cached(pos,i,rotmat_,moln1,moln2);
			vec2+=TM_Vector_Score_Backward_Cached(pos,i,rotmat_,moln1,moln2);
		}
		if(TM_cb1!=0 && TM_cb2!=0)
		{
			vec3=TM_Vector_Score_CB_Cached(pos,i,rotmat_,moln1,moln2);
			tot_num=3;
		}
		else vec3=0.0;
		vect=1.0*(vec1+vec2)/range;
		vect=1.0*(vect+vec3)/tot_num;
		tms=1.0/(1.0+dist2/ori_d);
		//total score
		int cur_index=pos*moln2;
		if(TMali_Weight==0) //no weight
		{
			if(TM_Vect_Score==0)ws_sco+=tms;            // original TMScorer
			else if(TM_Vect_Score==1)ws_sco+=vect;      // original Vect_score
			else ws_sco+=vect*tms;  // current Vect_score//
		}
		else                //has weight
		{
			if(TM_Wei_Score==1)ws_sco+=TMali_Weight[cur_index+i];
			else
			{
				if(TM_Vect_Score==0)ws_sco+=TMali_Weight[cur_index+i]*tms;       // original TMScorer
				else if(TM_Vect_Score==1)ws_sco+=TMali_Weight[cur_index+i]*vect; // original TMScorer
				else
				{
					ws_sco+=TMali_Weight[cur_index+i]*vect*tms; // current Vect_score//default:aligner.TMali_Weight = wswei;
				}
			}
		}
	}
	//return
	return ws_sco;
}
//[note]: the following function will take the input of the MatchWei, 
//        whose length should be identical to the match states.
double TMAligner::TM_Align_Get_Score_Simp_MatchWei(XYZ *mol1,XYZ *mol2,double *rotmat_,
	int moln1,int moln2,int *ali2,vector <double> & MatchWei)
{
	int i,k;
	double vec1,vec2,vec3;
	double vect;
	int tot_num;
	int range=1;
	double dist2;
	double ori_d=d0*d0;
	double tms;
	int pos;
	XYZ xyz;
	double ws_sco=0.0;
	//calc_score
	TM_Precompute_Vector_Cache(mol1,mol2,moln1,moln2);
	int lali=0;
	for(i=0;i<moln2;i++)
	{
		pos=ali2[i];
		if(pos<0 || pos>=moln1)continue;
		//single
		if(rotmat_!=0)TMs_Cache_Point(mol1,pos,xyz,rotmat_);
		else xyz=mol1[pos];
		dist2=mol2[i].distance_square(xyz);
		//nearby
		vec1=0.0;
		vec2=0.0;
		//calc
		tot_num=2;
		for(k=1;k<=range;k++)
		{
			vec1+=TM_Vector_Score_Forward_Cached(pos,i,rotmat_,moln1,moln2);
			vec2+=TM_Vector_Score_Backward_Cached(pos,i,rotmat_,moln1,moln2);
		}
		if(TM_cb1!=0 && TM_cb2!=0)
		{
			vec3=TM_Vector_Score_CB_Cached(pos,i,rotmat_,moln1,moln2);
			tot_num=3;
		}
		else vec3=0.0;
		vect=1.0*(vec1+vec2)/range;
		vect=1.0*(vect+vec3)/tot_num;
		tms=1.0/(1.0+dist2/ori_d);
		//total score
		ws_sco+=MatchWei.at(lali)*vect*tms; // current Vect_score
		lali++;
	}
	//return
	return ws_sco;
}

//----------- Test_Score --------------//
//given alignment, return each parts of our score

//====================== minor functions ====================//
//---- vector_score -----//
static int TM_Normalize_Vector(XYZ diff,XYZ &out,double min_len,double max_len)
{
	double len=sqrt(diff.X*diff.X+diff.Y*diff.Y+diff.Z*diff.Z);
	if(len>max_len || len<min_len)return 0;
	out.X=diff.X/len;
	out.Y=diff.Y/len;
	out.Z=diff.Z/len;
	return 1;
}
static void TM_Rotate_Vector(XYZ in,XYZ &out,double *rotmat_)
{
	if(rotmat_==0)
	{
		out=in;
		return;
	}
	out.X=in.X*rotmat_[0]+in.Y*rotmat_[1]+in.Z*rotmat_[2];
	out.Y=in.X*rotmat_[3]+in.Y*rotmat_[4]+in.Z*rotmat_[5];
	out.Z=in.X*rotmat_[6]+in.Y*rotmat_[7]+in.Z*rotmat_[8];
}
static double TM_Vector_Dot_Score(XYZ vec1,XYZ vec2,double *rotmat_)
{
	XYZ rot_vec1;
	double score;
	TM_Rotate_Vector(vec1,rot_vec1,rotmat_);
	score=rot_vec1.X*vec2.X+rot_vec1.Y*vec2.Y+rot_vec1.Z*vec2.Z;
	score=limit(score,-1.0,1.0);
	return (score+1.0)/2;
}
void TMAligner::TM_Precompute_Vector_Cache(XYZ *mol1,XYZ *mol2,int moln1,int moln2)
{
	TM_Precompute_Vector_Cache_One(mol1,TM_cb1,moln1,
		TM_vec_forw1,TM_vec_back1,TM_vec_cb1,
		TM_vec_forw_valid1,TM_vec_back_valid1,TM_vec_cb_valid1);
	TM_Precompute_Vector_Cache_One(mol2,TM_cb2,moln2,
		TM_vec_forw2,TM_vec_back2,TM_vec_cb2,
		TM_vec_forw_valid2,TM_vec_back_valid2,TM_vec_cb_valid2);
}
void TMAligner::TM_Precompute_Vector_Cache_One(XYZ *mol,XYZ *cb,int moln,
	XYZ *forw_vec,XYZ *back_vec,XYZ *cb_vec,
	char *forw_valid,char *back_valid,char *cb_valid)
{
	int i;
	XYZ diff;
	for(i=0;i<moln;i++)
	{
		forw_valid[i]=0;
		back_valid[i]=0;
		cb_valid[i]=0;
		forw_vec[i]=0.0;
		back_vec[i]=0.0;
		cb_vec[i]=0.0;
	}
	for(i=0;i<moln;i++)
	{
		if(i-1>=0)
		{
			diff=mol[i]-mol[i-1];
			forw_valid[i]=(char)TM_Normalize_Vector(diff,forw_vec[i],2.5,4.5);
		}
		if(i+1<moln)
		{
			diff=mol[i]-mol[i+1];
			back_valid[i]=(char)TM_Normalize_Vector(diff,back_vec[i],2.5,4.5);
		}
		if(cb!=0)
		{
			diff=mol[i]-cb[i];
			cb_valid[i]=(char)TM_Normalize_Vector(diff,cb_vec[i],1.0,2.0);
		}
	}
}
double TMAligner::TM_Vector_Score_Forward_Cached(int ii,int jj,double *rotmat_,int moln1,int moln2)
{
	if(ii-1<0||jj-1<0)return 0.5;
	if(ii>=moln1||jj>=moln2)return 0.5;
	if(TM_vec_forw_valid1[ii]==0||TM_vec_forw_valid2[jj]==0)return 0.0;
	return TM_Vector_Dot_Score(TM_vec_forw1[ii],TM_vec_forw2[jj],rotmat_);
}
double TMAligner::TM_Vector_Score_Backward_Cached(int ii,int jj,double *rotmat_,int moln1,int moln2)
{
	if(ii<0||jj<0)return 0.5;
	if(ii+1>=moln1||jj+1>=moln2)return 0.5;
	if(TM_vec_back_valid1[ii]==0||TM_vec_back_valid2[jj]==0)return 0.0;
	return TM_Vector_Dot_Score(TM_vec_back1[ii],TM_vec_back2[jj],rotmat_);
}
double TMAligner::TM_Vector_Score_CB_Cached(int ii,int jj,double *rotmat_,int moln1,int moln2)
{
	if(ii<0||jj<0)return 0.5;
	if(ii>=moln1||jj>=moln2)return 0.5;
	if(TM_vec_cb_valid1[ii]==0||TM_vec_cb_valid2[jj]==0)return 0.0;
	return TM_Vector_Dot_Score(TM_vec_cb1[ii],TM_vec_cb2[jj],rotmat_);
}

double TMAligner::TM_Align_Get_Pair_Score(int i,int j,double distance_square,
	double d0,double *rotmat_,int moln1,int moln2)
{
	double vec1=TM_Vector_Score_Forward_Cached(i,j,rotmat_,moln1,moln2);
	double vec2=TM_Vector_Score_Backward_Cached(i,j,rotmat_,moln1,moln2);
	double vec3=0.0;
	int tot_num=2;
	if(TM_cb1!=0 && TM_cb2!=0)
	{
		vec3=TM_Vector_Score_CB_Cached(i,j,rotmat_,moln1,moln2);
		tot_num=3;
	}
	const double vect=(vec1+vec2+vec3)/tot_num;
	const double tms=1.0/(1.0+distance_square/(d0*d0));
	if(TMali_Weight==0)
	{
		if(TM_Vect_Score==0)return tms;
		if(TM_Vect_Score==1)return vect;
		return vect*tms;
	}
	const double weight=TMali_Weight[i*moln2+j];
	if(TM_Wei_Score==1)return weight;
	if(TM_Vect_Score==0)return weight*tms;
	if(TM_Vect_Score==1)return weight*vect;
	return weight*vect*tms;
}

//simplest
void TMAligner::TM_Align_Get_Matrix(XYZ *mol1,XYZ *mol2,int moln1,int moln2,double d0,double *score)
{
	//calc_score
	TM_Precompute_Vector_Cache(mol1,mol2,moln1,moln2);
	for(int i=0;i<moln1;i++)
	{
		int cur_index=i*moln2;
		for(int j=0;j<moln2;j++)
		{
			const double distance_square=mol1[i].distance_square(mol2[j]);
			score[cur_index+j]=TM_Align_Get_Pair_Score(i,j,distance_square,
				d0,0,moln1,moln2);
		}
	}
}
void TMAligner::TM_Align_Get_Matrix_TMs(XYZ *mol1,XYZ *mol2,int moln1,int moln2,double d0,double *score)
{
	int i,j;
	double dist2;
	double ori_d=d0*d0;
	for(i=0;i<moln1;i++)
	{
		int cur_index=i*moln2;
		for(j=0;j<moln2;j++)
		{
			dist2=mol1[i].distance_square(mol2[j]);
			score[cur_index+j]=1.0/(1.0+dist2/ori_d);  // this is NOT d0 !!!
		}
	}
}

//[note]: the following function will take the input of the MatchWei, 
//        whose length should be identical to the match states.

//input:  original mol1,moln1 and mol2,moln2; 
//        their correspondence set, ali2
//output: TMscore and the score_matrix
double TMAligner::TM_Align_Get_Score(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,int NOTM)
{
	int i;
	int lali;
	double tmscore=0.0;
	//init judge
	int smaller=moln1<moln2?moln1:moln2;
	if(smaller==0)
	{
		for(i=0;i<moln2;i++)ali2[i]=-1;
		return tmscore;
	}
	if(TM_CALC==0)Calc_TM_d0(smaller);
	//get correspondece
	lali=TM_Align_Get_XYZ(mol1,mol2,moln1,moln2,ali2);
	if(lali<=1)return tmscore;
	if(lali<=2)
	{
		kabsch(TM_tmp2,TM_tmp1,lali,finmat);
		tmscore=Calc_TM_Score_Single(TM_tmp1,TM_tmp2,lali,finmat,d0,d8,0);
		tmscore=tmscore/smaller;
	}
	else
	{
		//calculate simple TMscore8
		if(NOTM==1)
		{
			kabsch(TM_tmp2,TM_tmp1,lali,finmat);
			tmscore=Calc_TM_Score_Single(TM_tmp1,TM_tmp2,lali,finmat,d0,d8,0);
			tmscore=tmscore/smaller;
		}
		else
		{
			tmscore=Calc_TM_Score(TM_tmp1,TM_tmp2,lali,d0,d8,1,1);
			tmscore=tmscore/smaller;
		}
	}
	//calculate score matrix
	rot_mol(mol1,TM_tmp1,moln1,finmat);
	if(TM_GAP_TYPE==-1)TM_Align_Get_Matrix_Bound(TM_tmp1,mol2,moln1,moln2,d0,TM_bound,TM_DP_sco); //DynaProg bound//__110720__//
	else
	{
		if(TM_Score_Type==0)TM_Align_Get_Matrix_TMs(TM_tmp1,mol2,moln1,moln2,d0,TM_DP_sco);  //TMscore type
		else TM_Align_Get_Matrix(TM_tmp1,mol2,moln1,moln2,d0,TM_DP_sco);                     //CombinedScore type
	}
	//return
	return tmscore;
}

//----------- TMAligner Functions -------------//
//input:  original mol1,moln1 and mol2,moln2; 
//        their transformation score
//[note]: mol1 and mol2 MUST be superimposed
//output: correspondence path -> ali2

void TMAligner::TM_Align_Get_Ali_coarse_fullDP(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
	int i;
	//init judge
	int smaller=moln1<moln2?moln1:moln2;
	if(smaller==0)
	{
		for(i=0;i<moln2;i++)ali2[i]=-1;
		return;
	}
	if(TM_CALC==0)Calc_TM_d0(smaller);
	//calculate score matrix
	TM_Align_Get_Matrix(mol1,mol2,moln1,moln2,d0,TM_DP_sco);     //CombinedScore type  on coarse stage

	// you can check result at, which show that using CombinedScore type on coarse stage is better than TMscore type
	// /home/ziting/Code/StructureAlignment/EAlign/test_results_fixed/full_alignment_EAlign_Dual_Staged_Evolution_Geometry_CEFragment_GeoFullDP_iteration10_SFP_strategy4_adjust_CLEBLOSUM_weight_2_Evoonly/EAlign_Dual_Staged_Evolution_Geometry_CEFragment_GeoFullDP_iteration10_SFP_strategy4_adjust_CLEBLOSUM_weight_2_Evoonly_alignment_benchmark.csv
	// &
	// /home/ziting/Code/StructureAlignment/EAlign/test_results_fixed/full_alignment_EAlign_Dual_Staged_Evolution_Geometry_CEFragment_CombinedScoreFullDP_iteration10_SFP_strategy4_adjust_CLEBLOSUM_weight_2_Evoonly/EAlign_Dual_Staged_Evolution_Geometry_CEFragment_CombinedScoreFullDP_iteration10_SFP_strategy4_adjust_CLEBLOSUM_weight_2_Evoonly_alignment_benchmark.csv

	//calculate_path
	// TM_Align_Dyna_Prog(moln1,moln2,TM_DP_sco,ali2,0);
	//--- run DynaProg ----//
	int align_len;
	double ali_sco;
	
	align_len=Normal_Align_Dyna_Prog_II(moln1,moln2,TM_DP_sco,0,0,TM_dynaprog,ali_sco);
	
	//--- extract alignment ---//
	int ii,jj;
	int totnum=(int)TM_dynaprog.size();
	memset(ali2,-1,moln2*sizeof(int));
	for(i=0;i<totnum;i++)
	{
		ii=TM_dynaprog[i].first;
		jj=TM_dynaprog[i].second;
		if(ii>0&&jj>0)ali2[jj-1]=ii-1;
	}
}

//------ calc TMScorer -------// return real TMscore
double TMAligner::TM_Align_TM_Score(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,
	int norm_len,double norm_d0,double &rmsd,int &lali,double *MAXSCO)
{
	int i;
	double tmscore8;
	double tmscore=0.0;
	rmsd=0;
	lali=0;
	if(MAXSCO!=NULL)for(int wwi=0;wwi<8;wwi++)MAXSCO[wwi]=0;
	//init judge
	int smaller=moln1<moln2?moln1:moln2;
	if(smaller==0)
	{
		for(i=0;i<moln2;i++)ali2[i]=-1;
		return tmscore;
	}
	if(TM_CALC==0)Calc_TM_d0(smaller);
	//get correspondece
	lali=TM_Align_Get_XYZ(mol1,mol2,moln1,moln2,ali2);
	if(lali<=1)return tmscore;
	if(lali<=2)
	{
		rmsd=kabsch(TM_tmp2,TM_tmp1,lali,finmat);
		if(rmsd>0.0)rmsd=1.0*sqrt(rmsd);
		tmscore=Calc_TM_Score_Single(TM_tmp1,TM_tmp2,lali,finmat,norm_d0,d8,0,MAXSCO);
		tmscore=tmscore/norm_len;
		return tmscore;
	}
	//calculate TMscore8
	if(TM_DIST_CUT==1)
	{
		tmscore8=Calc_TM_Score(TM_tmp1,TM_tmp2,lali,d0,d8,0,1);  //calculate TM8
		tmscore8=tmscore8/smaller;
		//remove dis>d8 in normal TM-score calculation for final report----->
		if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*moln1);
		lali=TM_Align_Get_CUT(mol1,mol2,moln1,moln2,ali2,TM_DistCut,finmat);
		if(lali<=1)return tmscore;
		if(lali<=2)
		{
			rmsd=kabsch(TM_tmp2,TM_tmp1,lali,finmat);
			if(rmsd>0.0)rmsd=1.0*sqrt(rmsd);
			tmscore=Calc_TM_Score_Single(TM_tmp1,TM_tmp2,lali,finmat,norm_d0,d8,0,MAXSCO);
			tmscore=tmscore/norm_len;
			return tmscore;
		}
	}
	//calculate TMscore
	tmscore=Calc_TM_Score(TM_tmp1,TM_tmp2,lali,norm_d0,d8,0,0,MAXSCO);
	tmscore=tmscore/norm_len;
	//return
	rmsd=kabsch(TM_tmp2,TM_tmp1,lali,0);
	if(rmsd>0.0)rmsd=1.0*sqrt(rmsd);
	return tmscore;
}
double TMAligner::TM_Align_TM_Score(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,
	double &rmsd,int &lali,double *MAXSCO)
{
	int smaller=moln1<moln2?moln1:moln2;
	return TM_Align_TM_Score(mol1,mol2,moln1,moln2,ali2,smaller,d0,rmsd,lali,MAXSCO);
}
//return simple TMscore (may a little bit smaller than real TMscore, but more efficient!)
double TMAligner::TM_Align_TM_Score_Simp(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,
	int norm_len,double norm_d0,double &rmsd,int &lali,double *MAXSCO)
{
	int i;
	double tmscore=0.0;
	rmsd=0;
	lali=0;
	if(MAXSCO!=NULL)for(int wwi=0;wwi<8;wwi++)MAXSCO[wwi]=0;
	//init judge
	int smaller=moln1<moln2?moln1:moln2;
	if(smaller==0)
	{
		for(i=0;i<moln2;i++)ali2[i]=-1;
		return tmscore;
	}
	if(TM_CALC==0)Calc_TM_d0(smaller);
	//get correspondece
	lali=TM_Align_Get_XYZ(mol1,mol2,moln1,moln2,ali2);
	if(lali<=1)return tmscore;
	if(lali<=2)
	{
		rmsd=kabsch(TM_tmp2,TM_tmp1,lali,finmat);
		if(rmsd>0.0)rmsd=1.0*sqrt(rmsd);
		tmscore=Calc_TM_Score_Single(TM_tmp1,TM_tmp2,lali,finmat,norm_d0,d8,0,MAXSCO);
		tmscore=tmscore/norm_len;
		return tmscore;
	}
	//calculate TMscore8
	if(TM_DIST_CUT==1)
	{
		//remove dis>d8 in normal TM-score calculation for final report----->
		if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*moln1);
		lali=TM_Align_Get_CUT(mol1,mol2,moln1,moln2,ali2,TM_DistCut,finmat);
		if(lali<=1)return tmscore;
		if(lali<=2)
		{
			rmsd=kabsch(TM_tmp2,TM_tmp1,lali,finmat);
			if(rmsd>0.0)rmsd=1.0*sqrt(rmsd);
			tmscore=Calc_TM_Score_Single(TM_tmp1,TM_tmp2,lali,finmat,norm_d0,d8,0,MAXSCO);
			tmscore=tmscore/norm_len;
			return tmscore;
		}
	}
	//calculate TMscore
	tmscore=Calc_TM_Score(TM_tmp1,TM_tmp2,lali,norm_d0,d8,1,1,MAXSCO);
	tmscore=tmscore/norm_len;
	//return
	rmsd=kabsch(TM_tmp2,TM_tmp1,lali,0);
	if(rmsd>0.0)rmsd=1.0*sqrt(rmsd);
	return tmscore;
}
double TMAligner::TM_Align_TM_Score_Simp(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2,
	double &rmsd,int &lali,double *MAXSCO)
{
	int smaller=moln1<moln2?moln1:moln2;
	return TM_Align_TM_Score_Simp(mol1,mol2,moln1,moln2,ali2,smaller,d0,rmsd,lali,MAXSCO);
}
//return simplest TMscore (just use the raw rotmat)


//================ bound functions =================//__110730__//
void TMAligner::TM_Align_Get_Bound(int moln1,int moln2,vector<pair<int,int> > &align,
	vector<pair<int,int> > &bound,int neib)
{
	int k;
	int ii,jj;
	int pre_ii,pre_jj;
	int size=(int)align.size();
	//[1]get real alignment
	bound[0].first=0;
	bound[0].second=0;
	pre_jj=0;
	for(k=0;k<size;k++)
	{
		ii=align[k].first;
		jj=align[k].second;
		if(ii>0&&jj>0)
		{
			bound[ii].first=jj;  //positive
			bound[ii].second=jj;
			pre_ii=ii+1;
			pre_jj=jj+1;
		}
		else
		{
			if(ii>0)
			{
				bound[ii].first=pre_jj;
				bound[ii].second=pre_jj;
			}
			else
			{
				bound[abs(ii)].second++;
			}
		}
	}
	//[2]get horizontal expand
	for(k=0;k<=moln1;k++)
	{
		ii=bound[k].first;
		jj=bound[k].second;
		if(ii-neib<0)ii=0;
		else ii=ii-neib;
		if(jj+neib>moln2)jj=moln2;
		else jj=jj+neib;
		bound[k].first=ii;
		bound[k].second=jj;
	}
	//[3]get vertical expand
	int i;
	int first=1;
	int wscurr,wsprev;
	int wscurr_start=1;
	int wsprev_start=1;
	int ww1,ww2;
	wsprev=0;
	wscurr=0;
	for(k=0;k<size;k++)
	{
		ii=align[k].first;
		jj=align[k].second;
		if(ii>0&&jj>0)
		{
			if(first==1)
			{
				wscurr=ii;
				wscurr_start=jj;
				first=0;
				//record
				for(i=wsprev-neib;i<=wscurr+neib;i++)
				{
					if(i<1)continue;
					if(i>moln1)break;
					ww1=wsprev_start-neib;
					if(ww1<1)ww1=1;
					ww2=wscurr_start+neib;
					if(ww2>moln2)ww2=moln2;
					if(ww1<bound[i].first)bound[i].first=ww1;
					if(ww2>bound[i].second)bound[i].second=ww2;
				}
			}
		}
		else
		{
			if(first==0)
			{
				wsprev=abs(align[k-1].first)+1;
				wsprev_start=abs(align[k-1].second)+1;
				first=1;
			}
		}
	}
	//final record
	if(wsprev==moln1+1)
	{
		wscurr_start=moln2;
		for(i=wsprev-neib-1;i<=moln1;i++)
		{
			if(i<1)continue;
			if(i>moln1)break;
			ww1=wsprev_start-neib;
			if(ww1<1)ww1=1;
			ww2=wscurr_start+neib;
			if(ww2>moln2)ww2=moln2;
			if(ww1<bound[i].first)bound[i].first=ww1;
			if(ww2>bound[i].second)bound[i].second=ww2;
		}
	}
}
void TMAligner::TM_Align_Get_Ali2_Bound(int moln1,int moln2,int *ali2,vector<pair<int,int> > &bound,int neib)
{
	//[1]from ali2 to ali1
	int i;
	int ii,jj;
	memset(TM_ali1,-1,sizeof(int)*moln1);
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]>=0)
		{
			ii=ali2[i];
			jj=i;
#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif	
			TM_ali1[ii]=jj;
		}
	}
	//[2]from ali1 to alignment
	int j;
	int wlen;
	int pre_ii=0;
	int pre_jj=0;
	int count=0;
	for(i=1;i<=moln1;i++)
	{
		ii=i;
		jj=TM_ali1[i-1];  //ali1 starts from 0, correspondence also from 0
		if(jj==-1)
		{
			continue;
		}
		else
		{
			jj++;
			//previous_path
			wlen=ii-pre_ii;
			for(j=1;j<wlen;j++)
			{
				pre_ii++;
				DP_align1[count]=pre_ii;
				DP_align2[count]=-pre_jj;
				count++;
			}
			wlen=jj-pre_jj;
			for(j=1;j<wlen;j++)
			{
				pre_jj++;
				DP_align1[count]=-pre_ii;
				DP_align2[count]=pre_jj;
				count++;
			}
			//current_path
			DP_align1[count]=ii;
			DP_align2[count]=jj;
			count++;
			//update
			pre_ii=ii;
			pre_jj=jj;
		}
	}
	pre_ii++;
	for(i=pre_ii;i<=moln1;i++)
	{
		DP_align1[count]=i;
		DP_align2[count]=-pre_jj;
		count++;
	}
	pre_jj++;
	for(i=pre_jj;i<=moln2;i++)
	{
		DP_align1[count]=-moln1;
		DP_align2[count]=i;
		count++;
	}
	//[3] from alignment to bound
	vector<pair<int,int> > alignment;
	alignment.resize(count);
	for(i=0;i<count;i++)
	{
		alignment[i].first=DP_align1[i];
		alignment[i].second=DP_align2[i];
	}
	TM_Align_Get_Bound(moln1,moln2,alignment,bound,neib);
}
void TMAligner::TM_Align_Get_Matrix_Bound(XYZ *mol1,XYZ *mol2,int moln1,int moln2,double d0,
	vector<pair<int,int> > &bound,double *score)
{
	int i,j,k;
	int range=1;
	int start,end;
	double dist2;
	double vec1,vec2;
	double vec3;
	double vect;
	double tms;
	int tot_num;
	double ori_d=d0*d0;

	//calc_score
	TM_Precompute_Vector_Cache(mol1,mol2,moln1,moln2);
	for(i=0;i<moln1;i++)
	{
		start=bound[i+1].first;
		end=bound[i+1].second;
		if(start==0)start=1;
		start--;
		int cur_index=i*moln2;
		for(j=start;j<end;j++)
		{
			vec1=0.0;
			vec2=0.0;
			//calc
			tot_num=2;
			for(k=1;k<=range;k++)
			{
				vec1+=TM_Vector_Score_Forward_Cached(i,j,0,moln1,moln2);
				vec2+=TM_Vector_Score_Backward_Cached(i,j,0,moln1,moln2);
			}
			if(TM_cb1!=0 && TM_cb2!=0)
			{
				vec3=TM_Vector_Score_CB_Cached(i,j,0,moln1,moln2);
				tot_num=3;
			}
			else vec3=0.0;
			vect=1.0*(vec1+vec2)/range;
			vect=1.0*(vect+vec3)/tot_num;
			dist2=mol1[i].distance_square(mol2[j]);
			tms=1.0/(1.0+dist2/ori_d);
			//total_score
			if(TMali_Weight==0) //no weight
			{
				if(TM_Vect_Score==0)score[cur_index+j]=tms;        // original TMScorer
				else if(TM_Vect_Score==1)score[cur_index+j]=vect;  // original Vect_score
				else score[cur_index+j]=vect*tms; // current Vect_score
			}
			else                //has weight
			{
				if(TM_Wei_Score==1)score[cur_index+j]=TMali_Weight[cur_index+j];
				else
				{
					if(TM_Vect_Score==0)score[cur_index+j]=TMali_Weight[cur_index+j]*tms;  // original TMScorer
					else if(TM_Vect_Score==1)score[cur_index+j]=TMali_Weight[cur_index+j]*vect;  // original TMScorer
					else score[cur_index+j]=TMali_Weight[cur_index+j]*vect*tms; // current Vect_score
				}
			}
		}
	}
}

/* ******************************************************************* */
/*     Dynamic programming for alignment. */
/*     Input: n1,n2,score(i,j), and gap_open */
/*     Output: invmap(j) */

/*     Please note this subroutine is not a correct implementation of */
/*     the N-W dynamic programming because the score tracks back only */
/*     one layer of the matrix. This code was exploited in TM-align */
/*     because it is about 1.5 times faster than a complete N-W code */
/*     and does not influence much the final structure alignment result. */
/* ******************************************************************* */
/* Subroutine */  /* mol2 is fixed!! superimpose mol1 onto mol2 */
int TMAligner::TM_Align_Dyna_Prog(int n1,int n2,double *score,int *ali2,
	double gapopen,double gapextend,int DP_Type)
{
	//--- run DynaProg ----//
	int align_len;
	double ali_sco;
	if(DP_Type==0)        //normal DynaProg (non-extend)
	{
		align_len=Normal_Align_Dyna_Prog_II(n1,n2,score,gapopen,gapextend,TM_dynaprog,ali_sco);
	}
	else if(DP_Type==1)   //advance DynaProg (four-state)
	{
		align_len=Advance_Align_Dyna_Prog_II(n1,n2,score,gapopen,gapextend,TM_dynaprog,ali_sco);
	}
	else                  //faster DynaProg (consider bound) //__110730__//
	{
		align_len=Normal_Align_Dyna_Prog_Fast(n1,n2,score,gapopen,gapextend,TM_bound,TM_dynaprog,ali_sco);
		TM_Align_Get_Bound(n1,n2,TM_dynaprog,TM_bound,TM_bound_neib);
	}

	//--- extract alignment ---//
	int i;
	int ii,jj;
	int totnum=(int)TM_dynaprog.size();
	memset(ali2,-1,n2*sizeof(int));
	for(i=0;i<totnum;i++)
	{
		ii=TM_dynaprog[i].first;
		jj=TM_dynaprog[i].second;
		if(ii>0&&jj>0)ali2[jj-1]=ii-1;
	}
	//--- final return ----//
	return align_len;
} /* dp_ */

//----------- TM_Align: Main_Function ----------//
//input:  mol1,mol2,moln1,moln2
//        initiail ali_path
//        iteration number
//output: best TMscore and final ali_path
double TMAligner::Calc_TM_Align(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_in,int *ali2_out,
	int INI_SKIP,int ITER_NUM,int SKIPorNOT)
{
	int id;
	//get_initial_matrix
	double TM_best=0.0;
	double TM_old=0.0;
	double TM_cur;
	double diff;
	memcpy(TM_DP_best,ali2_in,moln2*sizeof(int));
	if(INI_SKIP==-1)goto dp_start;
	else
	{
		if(TM_GAP_TYPE==-1)TM_Align_Get_Ali2_Bound(moln1,moln2,ali2_in,TM_bound,TM_bound_neib);//? ali2_in 建立一个沿 alignment path 的窄搜索区域 TM_bound
		TM_best=TM_Align_Get_Score(mol1,mol2,moln1,moln2,ali2_in,TM_NOTM);//mol1*mol2矩阵,TM_DP_sco建立
		if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*moln1);
		if(TM_Score_Type==1)TM_best=TM_Align_Get_Score_Simp(mol1,mol2,finmat,moln1,moln2,ali2_in);  //CombinedScore type
		if(TM_best==0.0)goto dp_end;
		memcpy(TM_rotmat,finmat,12*sizeof(double));
	}
	if(INI_SKIP==1)goto dp_end;//不会调用 DP,留下TM_DP_sco:当前姿态下的下一轮 DP 矩阵,finmat:当前姿态,TM_rotmat:当前姿态副本,TM_bound:可选的 DP band
	//start iteration
dp_start:
	for (id = 1; id <= ITER_NUM; ++id)
	{
		TM_Align_Dyna_Prog(moln1,moln2,TM_DP_sco,TM_DP_ali2,
			0.0,0.0,TM_GAP_TYPE);
		TM_cur=TM_Align_Get_Score(mol1,mol2,moln1,moln2,
			TM_DP_ali2,TM_NOTM);
		if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*moln1);
		if(TM_Score_Type==1)
			TM_cur=TM_Align_Get_Score_Simp(mol1,mol2,finmat,
				moln1,moln2,TM_DP_ali2);
		if(TM_cur>TM_best)
		{
			memcpy(TM_DP_best,TM_DP_ali2,moln2*sizeof(int));
			TM_best=TM_cur;
			memcpy(TM_rotmat,finmat,12*sizeof(double));
		}
		if(SKIPorNOT==1)
		{
			if (id > 1)
			{
				diff = fabs(TM_cur-TM_old);
				if (diff < 1e-6f) goto wsout;
			}
			TM_old = TM_cur;
		}
wsout:
		;
	}
	//final
dp_end:
	memcpy(ali2_out,TM_DP_best,moln2*sizeof(int));
	return TM_best;
}

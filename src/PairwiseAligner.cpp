#include "PairwiseAligner.h"


//------ constructor -------//
PairwiseAligner::PairwiseAligner(int num,int CLESUM)
:TMAligner(num),SequenceScorer(CLESUM),
AlignmentGeometry(num),FragmentChainer(num),
SeedAligner(num),TMRefiner(num)
{
	Pairwise_Kill_Gap=0;   //default: don't Kill Gaps !!
	Kill_Frag=0;     //default: don't Kill Frag !!
	Kill_Frag_Para_0=25;  //-> min_small_len, 25
	Kill_Frag_Para_1=25;  //-> gap_len, 25
	Kill_Frag_Para_2=15;  //-> nongap_len, 15
}
PairwiseAligner::~PairwiseAligner(void)
{
}


//================= Alignment_Gap_Process ============//__130830__//
//[note]: gap_len=25, nongap_len=15


//---------- post_refine ----------//__110830__//
//to kill gaps!!
int PairwiseAligner::IN_Ali_To_Cor(int *Fragment_Cor,int thres,
						int moln1,int moln2,int *ali1,int *ali2)
{
	int i,k;
	int num;
	int ii,jj;
	int count;
	int isFirst;
	int isLast;
	int type;
	int head1,head2;
	int index;

	//init
	num=0;
	Fragment_Cor[0]=0;
	head1=-1;
	head2=-1;
	isLast=0;
	isFirst=1;
	type=0;
	count=INT_MIN_NUM;
	ii=-1;
	jj=-1;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]==-1) //purely blank
		{
			if(isFirst==0)
			{
				if(count>=thres) 
				{
					Fragment_Cor[0]++;
					index=Fragment_Cor[0];
					Fragment_Cor[4*index+0] = type;
					Fragment_Cor[4*index+1] = head1;
					Fragment_Cor[4*index+2] = head2;
					Fragment_Cor[4*index+3] = count;
					num+=count;
				}
				else
				{
					for(k=0;k<count;k++)
					{
						ali1[head1+k]=-1;
						ali2[head2+k]=-1;
					}
				}
				count=0;
				isFirst=1;
			}
			continue;
		}

		if(isFirst==1)
		{
ws_init:
			isFirst=0;
			ii=ali2[i];
			type=1; // >0 mode
			jj=i;
			count=1;
			head1=ii;
			head2=jj;
#ifdef DEBUG
	overrange_debug_test(ii,moln1);
	overrange_debug_test(jj,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif
			continue;
		}
		if(i==jj+1&&ali2[i]==ii+1)
		{
			ii=ali2[i];
			jj=i;
			count++;
			continue;
		}

ws_end:
		if(count>=thres) 
		{
			Fragment_Cor[0]++;
			index=Fragment_Cor[0];
			Fragment_Cor[4*index+0] = type;
			Fragment_Cor[4*index+1] = head1;
			Fragment_Cor[4*index+2] = head2;
			Fragment_Cor[4*index+3] = count;
			num+=count;
		}
		else
		{
			for(k=0;k<count;k++)
			{
				ali1[head1+k]=-1;
				ali2[head2+k]=-1;
			}
		}

		if(isLast==1)goto end;
		else goto ws_init;
	}

	if(count==INT_MIN_NUM)goto end;
	isLast=1;
	goto ws_end;
end:
	return num;
}
void PairwiseAligner::IN_Elongation(double *rotmat_,double thres,int cutoff,int *ali1,int *ali2,
						 XYZ *mol1,XYZ *mol2,int moln1,int moln2)
{
	int i;
	int ii,jj;
	int index1,index2;
	double Dou_Thres=thres*thres;
	int success;
	XYZ temp1,temp2,temp;
	int cut1,cut2,cut;

	//---- ws_record ---//
	int *wsali1=new int[moln1];
	int *wsali2=new int[moln2];
	memset(wsali1,0,moln1*sizeof(int));
	memset(wsali2,0,moln2*sizeof(int));
	for(i=0;i<moln2;i++)
	{
		ii=ali2[i];
		jj=i;
		if(ii<0)continue;
		wsali1[ii]=1;
		wsali2[jj]=1;
	}

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
			if(jj>=0 && ii>=0) break;
			else if(jj>=0 && ii<0)
			{
				//check
				if(wsali2[jj]==1)break;
				//calc
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
			else if(jj<0 && ii>=0)
			{
				//check
				if(wsali1[ii]==1)break;
				//calc
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
				if(mol2[index2].distance_square(temp) < Dou_Thres || cut>cutoff)
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
			if(jj>=0 && ii>=0) break;
			else if(jj>=0 && ii<0)
			{
				//check
				if(wsali2[jj]==1)break;
				//calc
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
			else if(jj<0 && ii>=0)
			{
				//check
				if(wsali1[ii]==1)break;
				//calc
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

	//delete
	delete [] wsali1;
	delete [] wsali2;
}
void PairwiseAligner::Pairwise_Align_Refine(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_in,int *ali2_out,double *rotmat_)
{
	int i;
	int pos;
	int *wsali1=new int[moln1];
	int *wsali2=new int[moln2];
	//init
	memset(wsali1,-1,moln1*sizeof(int));
	memset(wsali2,-1,moln2*sizeof(int));
	for(i=0;i<moln2;i++)
	{
		pos=ali2_in[i];
		if(pos!=-1)
		{
			wsali1[pos]=i;
			wsali2[i]=pos;
		}
	}
	//kill gaps
	IN_Ali_To_Cor(Fragment_Cor,pipeline_MIN_LEN,moln1,moln2,wsali1,wsali2);
	//refine iteration
	IN_Elongation(rotmat_,FIN_CUT,ANG_CUT,wsali1,wsali2,mol1,mol2,moln1,moln2);
	//evaluate
	for(i=0;i<moln2;i++)ali2_out[i]=-1;
	for(i=0;i<moln2;i++)
	{
		pos=wsali2[i];
		if(pos!=-1)ali2_out[i]=pos;
	}
	//delete
	delete [] wsali1;
	delete [] wsali2;

}

//==================== Gap_Kill ================//
//-------- Kill_Gap_Total ----------//
double PairwiseAligner::Fragment_Make_Score(Align_Record & align_record,
								  XYZ *mol1,XYZ *mol2,int moln1,int moln2,
								  double *rotmat,int *ali2)
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
	for(i=0;i<12;i++)align_record.rotmat[i]=rotmat[i];
	//calculate
	double d_0=d0;
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
		Ali_Cache_Point(mol1,ii,temp,rotmat);
		dist2=mol2[jj].distance_square(temp);
		sco+=1.0/(1.0+dist2/ori_d);
		rms+=dist2;
		count++;
	}
	//evaluate
	rmsd=sqrt(1.0*rms/count);
	lali=count;
	int TM_smaller=moln1<moln2?moln1:moln2;
	tms=1.0*sco/TM_smaller;
	//final
	align_record.RMSD=rmsd;
	align_record.lali=lali;
	align_record.TMsco=tms;
	align_record.main_sco=tms;
	return tms;
}
void PairwiseAligner::Fragment_Single_Match_Kill(XYZ *mol1,XYZ *mol2,int moln1,int moln2,vector <Align_Record> &tot)
{
	int i,k;
	int size=(int)tot.size();
	vector <int> refined_alignment(moln2);
	double rotmat[12];
	double rmsd;
	int lali;
	for(i=0;i<size;i++)
	{
		// Remove short AFPs and elongate the retained AFPs under the input pose.
		for(k=0;k<moln2;k++)refined_alignment[k]=tot[i].alignment[k];
		for(k=0;k<12;k++)rotmat[k]=tot[i].rotmat[k];
		if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
		if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*moln1);
		Pairwise_Align_Refine(mol1,mol2,moln1,moln2,&refined_alignment[0],
			&refined_alignment[0],rotmat);

		int TM_DIST_CUT_=TM_DIST_CUT;
		// Establish the pose used by the configured distance cutoff without
		// constructing an unused full pair-score matrix.
		TM_DIST_CUT=0;
		TM_Align_TM_Score_Simp(mol1,mol2,moln1,moln2,
			&refined_alignment[0],rmsd,lali);
		TM_DIST_CUT=TM_DIST_CUT_;
		TM_Align_TM_Score_Simp(mol1,mol2,moln1,moln2,
			&refined_alignment[0],rmsd,lali);

		// Refit a pose consistent with the final cleaned path.
		TM_DIST_CUT=0;
		TM_Align_TM_Score_Simp(mol1,mol2,moln1,moln2,
			&refined_alignment[0],rmsd,lali);
		TM_DIST_CUT=TM_DIST_CUT_;

		// Store the consistent path, pose and TM-score, then rank by CombinedScore.
		if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
		if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*moln1);
		Fragment_Make_Score(tot[i],mol1,mol2,moln1,moln2,finmat,
			&refined_alignment[0]);
		tot[i].main_sco=TM_Align_Get_Score_Simp(
			mol1,mol2,finmat,moln1,moln2,&refined_alignment[0]);
	}
}

//============================= the following is the main function for TMAligner =======================//
void PairwiseAligner::Pairwise_Get_Initial_CLE(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
	//assign SSE
	int *ws_cle1=new int[moln1+moln2];
	int *ws_cle2=new int[moln1+moln2];
	CLE_transform(IN_CLE1,ws_cle1);
	CLE_transform(IN_CLE2,ws_cle2);
	//get alignment path
	int i,j;
	for(i=0;i<moln1;i++)
	{
		int cur_index=i*moln2;
		for(j=0;j<moln2;j++)
		{
			TM_DP_sco[cur_index+j]=Universal_Calc(i,j,1,ws_cle1,ws_cle2,moln1,moln2,1);
		}
	}
	//apply dynamic programming
	double gap_open=-100.0;
	TM_Align_Dyna_Prog(moln1,moln2,TM_DP_sco,ali2,gap_open);
	//delete
	delete [] ws_cle1;
	delete [] ws_cle2;
}

//-------- EnvironmentScorer -------//(obselete)

//================= Pairwise_Align_Series ================//__110720__//
void PairwiseAligner::Alimeng_Add(vector <Align_Record> &align_tot,int moln2,double tms)
{
	Align_Record tmp_ali;
	memcpy(ali2,TM_Main_Ali2,moln2*sizeof(int));
	memcpy(pipeline_rotmat,TM_rotmat,12*sizeof(double));
	Pipeline_Make_Score(tmp_ali);
	tmp_ali.main_sco=tms;
	//init
	int j;
	double comp_check;
	int success=1; //defualt:1
	int cursize=(int)align_tot.size();
	for(j=0;j<cursize;j++)
	{
		comp_check=Alg_Comp_Simp(tmp_ali.alignment,align_tot[j].alignment,Comp_Check_Vari);
		if(comp_check>Comp_Check)
		{
			success=0;
			break;
		}
	}
	if(success==1)align_tot.push_back(tmp_ali);
}
double PairwiseAligner::Pairwise_Align_Total(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_,
									int norm_len,double norm_d0,double *MAXSCO)
{
	int k;
	double TM_CUR;
	double TM_BEST=0.0;
	//total_init
	AliT_Rec=0;  //default:NOT
	AliT_Max=0;  //default:NO
	memset(ali2_,-1,moln2*sizeof(int));
	memset(TM_Main_Ali2,-1,moln2*sizeof(int));

	//============ FragmentAligner ============//
	//FragmentAligner
//ws_init0:
	Pairwise_align_tot.clear();
	TM_CUR=Fragment_Main_Func(Pairwise_align_tot);
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		for(k=0;k<moln2;k++)ali2_[k]=pipeline_MAX.alignment[k];
		for(k=0;k<12;k++)TM_FINMAT[k]=pipeline_MAX.rotmat[k];
		AliT_Max=0;
	}
	//CLE DynaProg
	Pairwise_Get_Initial_CLE(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=8;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init5;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=8;
	}
	Alimeng_Add(Pairwise_align_tot,moln2,TM_CUR);

	//--- get_initial5 ---//local structure superposition(neo)
ws_init5:
	TM_Get_Initial5(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=5;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init6;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=5;
	}
	Alimeng_Add(Pairwise_align_tot,moln2,TM_CUR);

	//--- get_initial6 ---//local structure superposition(old)
ws_init6:
	if(TM_ADDITION==1)
	{
		TM_Get_Initial6(mol1,mol2,moln1,moln2,TM_Main_Ali2);
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=6;
		}
		if(TM_CUR<=TM_BEST*0.2)goto ws_init2;  //check bad...
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=6;
		}
		Alimeng_Add(Pairwise_align_tot,moln2,TM_CUR);
	}

	//--- get_initial2 ---//SSE
ws_init2:
	TM_Get_Initial2(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=2;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init1;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=2;
	}
	Alimeng_Add(Pairwise_align_tot,moln2,TM_CUR);

	//--- get_initial1 ---//gapless threading
ws_init1:
	TM_Get_Initial1(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=1;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init4;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=1;
	}
	Alimeng_Add(Pairwise_align_tot,moln2,TM_CUR);

	//--- get_initial1.5 ---//gapless threading advance
	if(TM_ADDITION==1 && AliT_Rec==1)
	{
		EqualArray(TM_Main_Ali2,TM_Main_AliT,moln2);
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=7;
		}
		if(TM_CUR<=TM_BEST*0.2)goto ws_init4;  //check bad...
		TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
		if(TM_CUR>TM_BEST)
		{
			TM_BEST=TM_CUR;
			memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
			memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
			AliT_Max=7;
		}
		Alimeng_Add(Pairwise_align_tot,moln2,TM_CUR);
	}

	//--- get_initial4 ---//
ws_init4:
	TM_Get_Initial4(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=4;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_init3;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=4;
	}
	Alimeng_Add(Pairwise_align_tot,moln2,TM_CUR);

	//--- get_initial3 ---//best+SSE
ws_init3:
	memcpy(TM_Main_Ali2,ali2_,moln2*sizeof(int));  //assign the current best path!!
	TM_Get_Initial3(mol1,mol2,moln1,moln2,TM_Main_Ali2);
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,1);  //just calc score
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=3;
	}
	if(TM_CUR<=TM_BEST*0.2)goto ws_end;  //check bad...
	TM_CUR=Calc_TM_Align(mol1,mol2,moln1,moln2,TM_Main_Ali2,TM_Main_Ali2,-1); //run DynaProg
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		memcpy(ali2_,TM_Main_Ali2,moln2*sizeof(int));
		memcpy(TM_FINMAT,TM_rotmat,12*sizeof(double));
		AliT_Max=3;
	}
	Alimeng_Add(Pairwise_align_tot,moln2,TM_CUR);

	//========= terminal =========//
ws_end:
	TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
	if(Pairwise_Kill_Gap==1)Fragment_Single_Match_Kill(mol1,mol2,moln1,moln2,Pairwise_align_tot);
	TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
	stable_sort(Pairwise_align_tot.rbegin(),Pairwise_align_tot.rend());   //descending sort
	for(k=0;k<moln2;k++)ali2_[k]=Pairwise_align_tot[0].alignment[k];
	EqualArray(TM_FINMAT,finmat,12);
	TM_BEST=Pairwise_align_tot[0].TMsco; //?TM-score?
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*moln1);

	//---- final add ----//__2016.05.30__//
	{
		//kill
		vector <Align_Record> Pairwise_align_tot_;
		for(int kk=0;kk<(int)Pairwise_align_tot.size();kk++)
		{
			for(int i=0;i<moln2;i++)TM_Main_Ali2[i]=Pairwise_align_tot[kk].alignment[i];
			for(int i=0;i<12;i++)TM_rotmat[i]=Pairwise_align_tot[kk].rotmat[i];
			TM_CUR=Pairwise_align_tot[kk].main_sco;
			Alimeng_Add(Pairwise_align_tot_,moln2,TM_CUR);
		}
		//add
		Pairwise_align_tot.clear();
		for(int kk=0;kk<(int)Pairwise_align_tot_.size();kk++)
		{
			Pairwise_align_tot.push_back(Pairwise_align_tot_[kk]);
		}
	}

	//return
	return TM_BEST;
}

//=========== Fragment_Normal ==========//__110715__//
double PairwiseAligner::Pairwise_Align_Normal(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2_,
									int norm_len,double norm_d0,double *MAXSCO)
{
	int k;
	double TM_CUR;
	double TM_BEST=0.0;

	//total_init
	AliT_Rec=0;  //default:NOT
	AliT_Max=0;  //default:NO
	memset(ali2_,-1,moln2*sizeof(int));
	memset(TM_Main_Ali2,-1,moln2*sizeof(int));

	//============ FragmentAligner ============//
	//FragmentAligner
//ws_init0:
	Pairwise_align_tot.clear();
	TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
	TM_CUR=Fragment_Main_Func(Pairwise_align_tot);//先用 CLE fragment 搜索多个初始刚体姿态，再在每个姿态中用 Low SFP、Kabsch 和 DP 构造候选 alignment，最后只对少量优质候选做 TM-score refinement。
	TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
	if(TM_CUR>TM_BEST)
	{
		TM_BEST=TM_CUR;
		for(k=0;k<moln2;k++)ali2_[k]=pipeline_MAX.alignment[k];
		for(k=0;k<12;k++)TM_FINMAT[k]=pipeline_MAX.rotmat[k];
		AliT_Max=0;
	}

	//---- ws_debug ---//
	int retval=(int)TM_CUR;
	int tm_origin=TM_ADDITION;
	if(retval==-1)TM_ADDITION=1;

ws_end:
	//---- ws_debug ---//
	if(retval==-1)TM_ADDITION=tm_origin;
	if(Pairwise_align_tot.size()==0)return TM_BEST;
	//---- ws_real_emd -//
	TM_GAP_TYPE=-1;  //faster DynaProg_bound //__110720__//
	if(Pairwise_Kill_Gap==1)Fragment_Single_Match_Kill(mol1,mol2,moln1,moln2,Pairwise_align_tot);//default: kill the gap
	TM_GAP_TYPE=0;   //reset to normal DynaProg //__110720__//
	stable_sort(Pairwise_align_tot.rbegin(),Pairwise_align_tot.rend());   //descending sort
	for(k=0;k<moln2;k++)ali2_[k]=Pairwise_align_tot[0].alignment[k];
	//EqualArray(TM_FINMAT,finmat,12);//?finmat 是最近一次 TM-score 计算留下的全局矩阵，未必属于排序后的 Pairwise_align_tot[0]?
	EqualArray(TM_FINMAT, Pairwise_align_tot[0].rotmat.data(), 12);
	TM_BEST=Pairwise_align_tot[0].TMsco;
	if(ALI_CACHE==1)memset(Ali_cache,0,sizeof(int)*moln1);
	if(TM_CACHE==1)memset(TMs_cache,0,sizeof(int)*moln1);

	//---- final add ----//__2016.05.30__//
	{
		//kill
		vector <Align_Record> Pairwise_align_tot_;
		for(int kk=0;kk<(int)Pairwise_align_tot.size();kk++)
		{
			for(int i=0;i<moln2;i++)TM_Main_Ali2[i]=Pairwise_align_tot[kk].alignment[i];
			for(int i=0;i<12;i++)TM_rotmat[i]=Pairwise_align_tot[kk].rotmat[i];
			TM_CUR=Pairwise_align_tot[kk].main_sco;
			Alimeng_Add(Pairwise_align_tot_,moln2,TM_CUR);
		}
		//add
		Pairwise_align_tot.clear();
		for(int kk=0;kk<(int)Pairwise_align_tot_.size();kk++)
		{
			Pairwise_align_tot.push_back(Pairwise_align_tot_[kk]);
		}
	}

	//return
	return TM_BEST;
}


//=========== FragmentAligner With Alignment ==========//


//=========== Fragment_Lite ==========//

//=========== Fragment_Fast ==========//


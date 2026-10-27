#include "AlignmentPipeline.h"

//------ constructor -------//
AlignmentPipeline::AlignmentPipeline(void)
{
	//[macros]
	//global
	HEAD_Chk=1;          // default 1 (do head check)
	TAIL_Chk=1;          // default 1 (do tail check)
	COMP_Chk=1;          // default 1 (do comparison check)
	FAST_Chk=0;          // default 0 (DO NOT use fast check)
	//partial
	ZoomOrNot=1;         // default 1 (do ZoomIn add)
	RefineOrNot=1;       // default 1 (do Refinement)
	NonLinear=1;         // default 1 (do Kill Nonlinear)
	FinalOrNot=1;        // default 1 (do Shrink or Elong)
	TMSorNOT=1;          // default 1 (apply TMscore dynamic programming)
	//[parameters]
	pipeline_MIN_LEN=4;      // default 4    (final cut length)
	ZOOM_ITER=3;         // default 3    (Zoom-In number)
	REFINE_ITER=10;      // default 10   (Optimization number)
	FIN_CUT=5.0;         // default [5.0]  => final distance cutoff
	INI_CUT=10.0;        // default [10.0] => initial distance cutoff
	ZM_TopK=20;          // default 20
	ZM_TopJ=PROT_MAX_NUM;        // default PROT_MAX_NUM (all)
}
AlignmentPipeline::~AlignmentPipeline(void)
{
}

//==================== main process ====================//-> part ZoomIn
//[input]
//m1 must be initially superimposed onto m2,
//SFP_L -> Similar_Fragment_Pair low-valued list
//[output]
//ali2 -> the final returned alignment
double AlignmentPipeline::Pipeline_Part(XYZ *m1,XYZ *m2,int n1,int n2,double *rotmat_,vector <SFP_Record> &SFP_Low,
								 Align_Record & align_record)
{
	int smaller=n1<n2?n1:n2;
	double ret_sco;
	//judge
	if(smaller==0)return -1.0;
	//init
	Pipeline_Input(m1,m2,n1,n2);
	Other_Fragment_Init();
	if(rotmat_==NULL)
	{
		for(int i=0;i<12;i++)pipeline_rotmat[i]=0.0;
		pipeline_rotmat[0]=1.0;
		pipeline_rotmat[4]=1.0;
		pipeline_rotmat[8]=1.0;
	}
	else EqualArray(pipeline_rotmat,rotmat_,12);
	//ZoomIn stragety
	int ret_val=0;
	if(ZoomOrNot==1)ret_val=ZoomIn_Add(ZOOM_ITER,INI_CUT,FIN_CUT,SFP_Low,pipeline_rotmat); //step1 -> ZoomIn: 利用 SFP 逐级收紧空间阈值，更新刚体变换，生成初始 alignment
	if(ret_val==-1)return -1;
	if(RefineOrNot==1)ret_val=Refinement(REFINE_ITER,FIN_CUT,pipeline_rotmat);             //step2 -> Refinement: alignment → Kabsch → 局部候选 DP → 新 alignment，交替优化&iteration<=10
	if(ret_val==-1)return -1;
	if(NonLinear==1)ret_val=Kill_NonLinear(pipeline_rotmat);                               //step3 -> Kill Nonlinear:删除顺序交叉、回折的对齐片段
	if(ret_val==-1)return -1;
	if(FinalOrNot==1)ret_val=Final_Stage(pipeline_rotmat);                                 //step4 -> Final Stage:延长、删除差匹配、重新 Kabsch
	if(ret_val==-1)return -1;
	if(TMSorNOT==1)ret_val=Final_Refine(pipeline_rotmat);                                  //step5 -> TMscore: default TMSorNot=0
	if(ret_val==-1)return -1;
	//terminal
	ret_sco=Pipeline_Make_Score(align_record);//align_record.alignment[i]=ali2[i];align_record.rotmat[i]=pipeline_rotmat[i];
	return ret_sco;
}

//==================== main process ====================//-> full ZoomIn
//[input]
//SFP_H -> Similar_Fragment_Pair high-valued list
//SFP_L -> Similar_Fragment_Pair low-valued list
//[output]
//ali2 -> the final returned alignment
double AlignmentPipeline::Pipeline_Full(XYZ *m1,XYZ *m2,int n1,int n2,vector <SFP_Record> &SFP_High,vector <SFP_Record> &SFP_Low,
								 vector <Align_Record> & align_records)
{
	int i;
	int retv;
	int index;
	int totnum;
	double ret_sco;
	//judge
	Align_Record align_record;
	int smaller=n1<n2?n1:n2;
	if(smaller==0)return -1.0;
	//init
	Pipeline_Input(m1,m2,n1,n2);
	Other_Fragment_Init();
	totnum=Select_Best_Pivot(1,INI_CUT,INI_CUT,SFP_High,ZM_TopJ); 
	Global_Pipeline_Init();
	//start
	int count=0;
	align_records.clear();
	pipeline_MAX.main_sco=0.0;
	for(i=0;i<totnum;i++)
	{
		//break_check
		retv=Pipeline_Break_Check(i,count,totnum);         //Check1 -> whether break or not -> if(count>=ZM_TopK)return 0!最多接受 TopK=20 个 alignment
		if(retv==0)break;
		//head_check
		index=Pipeline_Get_Index(i);
		if(HEAD_Chk==1)                                   //Check2 -> SFP_H check
		{
			retv=Pipeline_Head_Check(SFP_High[index]);
			if(retv==0)continue;
		}
		//Fragment_Part
		Pipeline_Initial_Rotmat(index,pipeline_rotmat);
		ret_sco=Pipeline_Part(m1,m2,n1,n2,pipeline_rotmat,SFP_Low,align_record);
		if(ret_sco<0.0)continue;
		//tail_check
		if(TAIL_Chk==1)                                   //Check3 -> alignment check
		{
			retv=Pipeline_Tail_Check(align_record);
			if(retv==0)continue;
		}
		//comparison_check
		if(COMP_Chk==1)                                   //Check4 -> comparison check
		{
			retv=Pipeline_Comparison_Check(align_record,align_records); //约 90% 的 residue-pair 重合则跳过,但是这里COMP_Chk = 0;
			if(retv==0)continue;
		}
		//terminal
		Pipeline_Make_Score(align_record);
		Pipeline_Terminal_Process(align_record);
		align_records.push_back(align_record);// TopK <= 20
		count++;
		//record
		if(align_record.main_sco>pipeline_MAX.main_sco)pipeline_MAX=align_record;
	}
	if(count==0)return -1.0;
	stable_sort(align_records.rbegin(),align_records.rend());   //descending sort
	return pipeline_MAX.main_sco;
}

#include "FragmentChainer.h"


//--------------constructor--------------------//
FragmentChainer::FragmentChainer(int num)
:AlignmentGeometry(num)
{
	Ali_Fragment_Maximal=num;
	// the following part is KILL_NonLinear
	Fragment_NOL_CUT=1;
	Fragment_GAP_CUT_1=4;               // minlen
	Fragment_GAP_CUT_2=INT_MAX_NUM;     // maxlen
	//create
	Fragment_Create_Linear(Ali_Fragment_Maximal);
}
FragmentChainer::~FragmentChainer(void)
{
	Fragment_Delete_Linear();
}

//---------- create & delete ---------//
void FragmentChainer::Fragment_Create_Linear(int IN_MAX_DIM)
{
	Fragment_temp=new int[IN_MAX_DIM];
	Fragment_indx=new int[IN_MAX_DIM];
	dp_cur_cor=new int[IN_MAX_DIM*5];
	dp_ori_cor=new int[IN_MAX_DIM*5];
}
void FragmentChainer::Fragment_Delete_Linear(void)
{
	delete [] Fragment_temp;
	delete [] Fragment_indx;
	delete [] dp_cur_cor;
	delete [] dp_ori_cor;
}


//Fragment_Cor (data_structure)//
//  Fragment_Cor[0*4+0] -> totnum
//  Fragment_Cor[i*4..] -> the i-th fragment
//  Fragment_Cor[i*4+0] -> RMSD of the frag + type ('+' or normal)
//  Fragment_Cor[i*4+1] -> start position in mol1 (ii)
//  Fragment_Cor[i*4+2] -> start position in mol2 (jj)
//  Fragment_Cor[i*4+3] -> fragment length (winlen)
//----------------//---------------------------------------------------------------------Chapter_II--PART_I
//--- kabsch  ---//
//--------------//
//[1]Fragment_Cor//(ws_modification) -> mol2 is fixed!!

//----------------// -> mol2 is fixed!!
//update_corespond (return len*[MAX_DIST^2 - RMSD^2] of each fragment)
//--------------------------------------------------//
void FragmentChainer::update_corespond(int *Fragment_Cor,double *rotmat_,
							   XYZ *mol1,XYZ *mol2,int moln1,int moln2)
{
	int i;
	int ii,jj,len;
	double SQURE_DIST=MAX_DIST*MAX_DIST;
	double value;
	int totnum=Fragment_Cor[0];
	for(i=1;i<=totnum;i++)
	{
		ii=Fragment_Cor[i*4+1];
		jj=Fragment_Cor[i*4+2];
		len=abs(Fragment_Cor[i*4+3]);
		value=calc_frag_dist(ii,jj,len,rotmat_,mol1,mol2,moln1,moln2);
		Fragment_Cor[i*4+0]=(int)(1.0*len*(SQURE_DIST-value));
	}
}

//----------------// -> mol2 is fixed!!
//update_corespond_II (return RMSD of each fragment)
//--------------------------------------------------//


//----------------// -> mol2 is fixed!!
//Fragment_Cor_Add (We provide two method: ->) 
// [0]check_conflict
// [1]First_Come_Never_Gone
//--------------------------------------------------//
void FragmentChainer::Fragment_Cor_Add(int *Fragment_Cor,double *rotmat_,int method,
						  XYZ *mol1,XYZ *mol2,int moln1,int moln2,
						  int *ali1,int *ali2)  

{
	int i,k,ii,jj;
	int totnum;
	int len;
	int index1,index2;
	double cur_score,ori_score;
	XYZ temp;

	//init
	Ali_Init(moln1,moln2,ali1,ali2);

    //Fragment_Cor->ali
	totnum=Fragment_Cor[0];
	for(i=1;i<=totnum;i++)
	{
		ii=Fragment_Cor[i*4+1];
		jj=Fragment_Cor[i*4+2];
		len=Fragment_Cor[i*4+3];
#ifdef DEBUG
	overrange_debug_test(ii+len-1,moln1);
	overrange_debug_test(jj+len-1,moln2);
	lessrange_debug_test(ii,0);
	lessrange_debug_test(jj,0);
#endif
		for(k=0;k<len;k++)
		{
			//check conflict
			if(ali1[ii+k]==-1 && ali2[jj+k]==-1)
			{
				ali1[ii+k]=jj+k;
				ali2[jj+k]=ii+k;
			}
			else
			{
				//method==1 (First_Come_Never_Gone)
				//method!=1 (check conflict)
				if(method==1)continue;
				//[0]calc ori_dist
				Ali_Cache_Point(mol1,ii+k,temp,rotmat_);
				ori_score=mol2[jj+k].distance_square(temp);
				index1=ali2[jj+k];
				index2=ali1[ii+k];
				//[1]
				if(index1==-1 && index2!=-1)
				{
					Ali_Cache_Point(mol1,ii+k,temp,rotmat_);
					cur_score=mol2[index2].distance_square(temp);
					if(ori_score<cur_score)
					{
						ali2[index2]=-1;
						ali1[ii+k]=jj+k;
						ali2[jj+k]=ii+k;
					}
				}
				//[2]
				if(index1!=-1 && index2==-1)
				{
					Ali_Cache_Point(mol1,index1,temp,rotmat_);
					cur_score=mol2[jj+k].distance_square(temp);
					if(ori_score<cur_score)
					{
						ali1[index1]=-1;
						ali1[ii+k]=jj+k;
						ali2[jj+k]=ii+k;
					}
				}
				//[3]
				if(index1!=-1 && index2!=-1)
				{
					Ali_Cache_Point(mol1,index1,temp,rotmat_);
					cur_score=mol2[index2].distance_square(temp);
					if(ori_score<cur_score)
					{
						ali1[index1]=-1;
						ali2[index2]=-1;
						ali1[ii+k]=jj+k;
						ali2[jj+k]=ii+k;
					}
				}
			}
		}//end of FOR(k)
	}//end of FOR(i)
}

//-------------//--------------------------------------------------------------------Chapter_I--PART_VI
//--Ali_To_Cor-//
//-----------//
int FragmentChainer::Ali_To_Cor(int *Fragment_Cor,int thres,
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

//------------------------------------//========================================================//
//-- FragKill_Dynamic_Programming  --//
//----------------------------------//
//================================== Kill NonLinear ======================//__080501__//
//------Kill_Irrelevant_Cor (ori)
//ori_cor[0][0] totnum
//ori_cor[0][1] realnum
//ori_cor[k][0] score
//ori_cor[k][1] head
//ori_cor[k][2] len
//ori_cor[k][3] winlen  <+,->
//ori_cor[k][4] correspondence[ii]  // because only mol1 will be nonlinear
//cutoff -> gap cutoff between two fragment(combine)
void FragmentChainer::Fragment_Combine_Cor(int *Fragment_Cor,int cutoff)  // mol1 bigger than mol2
{
	int i,j;
	int ii,jj,kk,ll;
	int len,add,winlen;
	int score;
	int cur;
	int div1,div2,div;
	int totnum;
	int count;
	int head;
	int isLast;
	totnum=Fragment_Cor[0];
	if(totnum<=0)return;

	//create dp_tmp_cor
	isLast=0;
	count=1;
	dp_ori_cor[count*5+1]=1;
	winlen=Fragment_Cor[1*4+3];
	for(i=1;i<=totnum-1;i++)
	{
		//cur_AFP
		ii=Fragment_Cor[i*4+1];
		jj=Fragment_Cor[i*4+2];
		add=Fragment_Cor[i*4+3];
		//nxt_AFP
		kk=Fragment_Cor[(i+1)*4+1];
		ll=Fragment_Cor[(i+1)*4+2];
		len=Fragment_Cor[(i+1)*4+3];
		//div
		div1=kk-(ii+add);
		div2=ll-(jj+add);
		div=div1>div2?div1:div2;
		if(div1<0||div>cutoff)
		{
last:
			dp_ori_cor[count*5+2]=i-dp_ori_cor[count*5+1]+1;
			dp_ori_cor[count*5+3]=winlen;
			//score
			score=0;
			cur=dp_ori_cor[count*5+1];
			for(j=0;j<dp_ori_cor[count*5+2];j++)score+=Fragment_Cor[(cur+j)*4+0];
			dp_ori_cor[count*5+0]=score;
			head=dp_ori_cor[count*5+1];
			dp_ori_cor[count*5+4]=Fragment_Cor[head*4+1];
			if(isLast==1)goto end;
			//count
			count++;
			dp_ori_cor[count*5+1]=i+1;
			winlen=Fragment_Cor[(i+1)*4+3];
		}
		else winlen+=len;
	}
	isLast=1;
	goto last;
end:
	dp_ori_cor[0*5+0]=count;
	dp_ori_cor[0*5+1]=totnum;
}

//------Kill_Irrelevant_Cor (cur)
//cur_cor[0][0] totnum
//cur_cor[0][1] totlen
//cur_cor[k][0] score
//cur_cor[k][1] correspondence[i]
//cur_cor[k][2] correspondence[j]
//cur_cor[k][3] winlen
//cur_cor[k][4] ori_cores
void FragmentChainer::Fragment_Create_Cur_Cor(int *Fragment_Cor)
{
	int i;
	int totnum;
	int count;
	int totlen;
	int ii,jj;
	int h2,t2;

	//create normal part
	totnum=dp_ori_cor[0*5+0];
	count=0;
	totlen=0;
	for(i=1;i<=totnum;i++)
	{
		if(dp_ori_cor[i*5+3]<0)continue;
		count++;
		dp_cur_cor[count*5+0]=dp_ori_cor[i*5+0];
		dp_cur_cor[count*5+3]=dp_ori_cor[i*5+3];
		dp_cur_cor[count*5+4]=i;
		totlen+=dp_ori_cor[i*5+3];
		Fragment_temp[count]=dp_ori_cor[i*5+4];
	}
	dp_cur_cor[0*5+0]=count;
	dp_cur_cor[0*5+1]=totlen;
	if(count==0)return;
	//create correspondence part
	fast_sort.fast_sort_1up(Fragment_temp,Fragment_indx,count,1);
	for(i=0;i<count;i++)
	{
		ii=i+1;
		jj=Fragment_indx[i];
		dp_cur_cor[ii*5+2]=jj;
		dp_cur_cor[jj*5+1]=ii;
	}
	//create cut_len_ii
	ii=dp_cur_cor[1*5+4];
	jj=dp_cur_cor[count*5+4];
	ii=dp_ori_cor[ii*5+1];
	jj=dp_ori_cor[jj*5+1]+dp_ori_cor[jj*5+2]-1;
	h2=Fragment_Cor[ii*4+2];
	t2=Fragment_Cor[jj*4+2]+Fragment_Cor[jj*4+3];

//	Fragment_GAP_CUT_2=(t2-h2)/2<100?(t2-h2)/2:100;
//	Fragment_GAP_CUT_2=INT_MAX_NUM;
}

//cutoff -> gap cutoff between two fragment(divide)
void FragmentChainer::Fragment_Linear_Process(int cutoff)
{
	//----starting dynamic programming ----//
	int i,j;
	int totnum;
	int jj,kk,ll;
	int score;
	int totscore;
	int wsmax,maxnum;
	int cur;
	int cutlen,curlen;

	//init
	jj=-1;
	score=INT_MIN_NUM;
	Fragment_max_sco=INT_MIN_NUM;   
	Fragment_max_i=-1;
	Fragment_max_j=-1;
	//create
	totnum=dp_cur_cor[0*5+0];
	Fragment_matrix.resize((totnum+1)*(totnum+1));
	Fragment_prerec.resize((totnum+1)*(totnum+1));
	//process
	for(i=0;i<totnum+1;i++) 
	{  	
		if(i!=0)
		{
			score=dp_cur_cor[i*5+0];  // range [0,500]
			jj=dp_cur_cor[i*5+1];
		}
		for(j=0;j<totnum+1;j++) 
		{
			//init score
			if(i==0||j==0)
			{
				Fragment_matrix[i*(totnum+1)+j]=0;
				Fragment_prerec[i*(totnum+1)+j]=0;
				continue;
			}
			cur=dp_cur_cor[j*5+2];			
			if(jj==j)totscore=score;
			else totscore=-1*(dp_cur_cor[i*5+3]+dp_cur_cor[cur*5+3]);
			//calculate the Fragment_matrix_value
			if(Fragment_matrix[(i-1)*(totnum+1)+j] - dp_cur_cor[i*5+3] >Fragment_matrix[i*(totnum+1)+j-1] - dp_cur_cor[cur*5+3])
			{
				wsmax=Fragment_matrix[(i-1)*(totnum+1)+j] - dp_cur_cor[i*5+3];
				maxnum=1;
			}
			else
			{
				wsmax=Fragment_matrix[i*(totnum+1)+j-1] - dp_cur_cor[cur*5+3];
				maxnum=2;
			}
			if(Fragment_matrix[(i-1)*(totnum+1)+j-1]+totscore >=wsmax)
			{
				wsmax=Fragment_matrix[(i-1)*(totnum+1)+j-1]+totscore;
				maxnum=0;
			}
			if(wsmax<=0)
			{
				Fragment_prerec[i*(totnum+1)+j]=0;
				Fragment_matrix[i*(totnum+1)+j]=0;
				continue;
			}
			//assign the pre_value
			if(maxnum==0)
			{
				kk=i-1;
				ll=j-1;
				if(totscore<0)curlen=-1*totscore;
				else curlen=0;
			}
			else if(maxnum==1)
			{
				kk=i-1;
				ll=j;
				curlen=dp_cur_cor[i*5+3];
			}
			else
			{
				kk=i;
				ll=j-1;
				curlen=dp_cur_cor[cur*5+3];
			}
			cutlen=Fragment_prerec[kk*(totnum+1)+ll]/10;
			if(cutlen+curlen>cutoff)
			{
				Fragment_prerec[i*(totnum+1)+j]=0;
				Fragment_matrix[i*(totnum+1)+j]=0;
				continue;
			}
			else
			{
				Fragment_prerec[i*(totnum+1)+j]=(cutlen+curlen)*10+maxnum;
				Fragment_matrix[i*(totnum+1)+j]=wsmax;
			}
			//record max_value
			if(wsmax>Fragment_max_sco)
			{
				Fragment_max_sco=wsmax;
				Fragment_max_i=i;
				Fragment_max_j=j;
			}
		}//end of j
	}//end of i
}

//cutoff -> the fragment length cutoff
int FragmentChainer::Fragment_Trace_Back(int *Fragment_Cor,int *Fragment_Cor_temp,int cutoff)
{
	int i,j,k,l;
	int prenum;
	int totnum;
	int ii,jj,kk;
	int count;
	int num;
	int totlen;
	int head;
	int len;
	int wsinx;

	//init
	totnum=dp_ori_cor[0*5+0];
	for(i=0;i<=totnum;i++)Fragment_temp[i]=0;
	//trace_back
	count=Fragment_Cor_temp[0*4+0];
	num=Fragment_Cor_temp[0*4+1];
	totlen=0;
	for(i=Fragment_max_i,j=Fragment_max_j;i>0||j>0;) 
	{
		if( Fragment_matrix[i*(totnum+1)+j] <= 0) break;
		prenum=Fragment_prerec[i*(totnum+1)+j]%10; 
		if(prenum==0) 
		{
			ii=i;
			jj=j;
			kk=dp_cur_cor[i*5+1];
			i--;
			j--;
			//check whether valid			
			if(jj!=kk)continue;
			//add to Fragment_Cor
			wsinx=dp_cur_cor[ii*5+4];
			head=dp_ori_cor[wsinx*5+1];
			len=dp_ori_cor[wsinx*5+2];
			for(k=0;k<len;k++)
			{
				num++;				
				for(l=0;l<4;l++)Fragment_Cor_temp[count*4+l]=Fragment_Cor[(k+head)*4+l];
				count--;
			}
			if(dp_ori_cor[wsinx*5+3]<0)fprintf(stderr,"Impossible[2]!!!\n");
			totlen+=dp_ori_cor[wsinx*5+3];
			if(Fragment_temp[wsinx]==1)fprintf(stderr,"Impossible[3]!!!\n");
			Fragment_temp[wsinx]=1;
		}
		else if(prenum==1)i--;
		else j--;
	}
	//check length
	if(totlen<cutoff)return 0;
	else
	{
		Fragment_Cor_temp[0*4+0]=count;
		Fragment_Cor_temp[0*4+1]=num;
		for(i=1;i<=totnum;i++)if(Fragment_temp[i]==1)dp_ori_cor[i*5+3]*=-1;  //used up
		return 1;
	}
}



//-------------- main_usage ------------//
void FragmentChainer::Fragment_Kill_NonLinear(int *ori_cor,int *bak_cor,double *rot_mat)
{
	int i,l;
	int ret_val;
	int count;
	int cur;

	//init_correspondence
	Ali_To_Cor(ori_cor,1,moln1,moln2,ali1,ali2);//alignment 转为连续 AFP
	update_corespond(ori_cor,rot_mat,mol1,mol2,moln1,moln2);//为每个 AFP 计算质量分数.这个分数随后用于决定 crossing 片段之间保留哪一组:片段越长、当前 pose 下距离越小的组
	Fragment_Combine_Cor(ori_cor,Fragment_GAP_CUT_1);//合并相邻 AFP:间隔不超过4
	if(ori_cor[0]==0)return;
	//kill_nonlinear
	bak_cor[0]=ori_cor[0];
	bak_cor[1]=0;
	for(i=0;i<Fragment_NOL_CUT;i++)
	{
		//check length
		Fragment_Create_Cur_Cor(ori_cor);
		if(dp_cur_cor[0*5+1]<Fragment_GAP_CUT_1)break;
		//start kill
		Fragment_Linear_Process(Fragment_GAP_CUT_2);
		ret_val=Fragment_Trace_Back(ori_cor,bak_cor,Fragment_GAP_CUT_1);//Traceback 并重建 alignment
		if(ret_val==0)break;
	}
	//insert_back
	count=0;
	cur=bak_cor[0];
	int ii,jj;
	int pre_ii,pre_jj;
	int totcount1,totcount2;
	//meghod[1]
	pre_ii=0;
	pre_jj=0;
	totcount1=0;
	fast_sort.fast_sort_2up(bak_cor,bak_cor[1],4,cur+1,1);
	for(i=1;i<=bak_cor[1];i++)
	{
		ii=bak_cor[(i+cur)*4+1];
		jj=bak_cor[(i+cur)*4+2];
		if(ii<pre_ii || jj<pre_jj)continue;
		pre_ii=ii;
		pre_jj=jj;
		totcount1+=bak_cor[(i+cur)*4+3];
	}
	//meghod[2]
	pre_ii=0;
	pre_jj=0;
	totcount2=0;
	fast_sort.fast_sort_2up(bak_cor,bak_cor[1],4,cur+1,2);
	for(i=1;i<=bak_cor[1];i++)
	{
		ii=bak_cor[(i+cur)*4+1];
		jj=bak_cor[(i+cur)*4+2];
		if(ii<pre_ii || jj<pre_jj)continue;
		pre_ii=ii;
		pre_jj=jj;
		totcount2+=bak_cor[(i+cur)*4+3];
	}
	//final method
	if(totcount1>totcount2)fast_sort.fast_sort_2up(bak_cor,bak_cor[1],4,cur+1,1);
	else fast_sort.fast_sort_2up(bak_cor,bak_cor[1],4,cur+1,2);
	pre_ii=0;
	pre_jj=0;
	for(i=1;i<=bak_cor[1];i++)
	{
		//final check
		ii=bak_cor[(i+cur)*4+1];
		jj=bak_cor[(i+cur)*4+2];
		if(ii<pre_ii || jj<pre_jj)continue;
		//final add
		count++;
		for(l=0;l<4;l++)ori_cor[count*4+l]=bak_cor[(i+cur)*4+l];
		pre_ii=ii;
		pre_jj=jj;
	}
	ori_cor[0]=count;
	//update_correspondence
	Fragment_Cor_Add(ori_cor,rot_mat,1,mol1,mol2,moln1,moln2,ali1,ali2);
}

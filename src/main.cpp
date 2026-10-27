#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

#include "EAlignConfig.h"
#include "SequenceScorer.h"
#include "PairwiseAligner.h"
#include "BetaCarbonBuilder.h"
#include "Foldseek3Di.h"
#include "Kabsch.h"
#include "StructureLoader.h"
#include "StructureWriter.h"

using namespace std;

StructureWriter mol_out;

int Ori_BLOSUM_62[21][21]={
{  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -5 },  //A
{ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -5 },  //R
{ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3, -5 },  //N
{ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3, -5 },  //D
{  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -5 },  //C
{ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2, -5 },  //Q
{ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2, -5 },  //E
{  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -5 },  //G
{ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3, -5 },  //H
{ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -5 },  //I
{ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -5 },  //L
{ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2, -5 },  //K
{ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -5 },  //M
{ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -5 },  //F
{ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -5 },  //P
{  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2, -5 },  //S
{  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -5 },  //T
{ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -5 },  //W
{ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -5 },  //Y
{  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -5 },  //V
{ -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5 }}; //Z
// A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   Z

//BLOSUM_Mapping//--------------ARNDCQEGHILKMFPSTWYVZ
int Blo_AA_Map_WS[21]=
{ 0,19, 4, 3, 6, 13,7, 8, 9, 17,11,10,12,2, 18,14,5, 1, 15,16,20};
//A  V  C  D  E  F  G  H  I  W  K  L  M  N  Y  P  Q  R   S  T  Z
//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17  18 19 20
//Ori_Mapping//-----------------AVCDEFGHIWKLMNYPQRSTZ
int Ori_AA_Map_WS[26]=
{ 0,20,2,3,4,5,6,7,8,20,10,11,12,13,20,15,16,17,18,19,20, 1, 9,20,14,20};
// A B C D E F G H I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
// 0 1 2 3 4 5 6 7 8  9 10 11 12 14 14 15 16 17 18 19 20 21 22 23 24 25

//------ calculate -------//
int BLOSUM62_Calc(char a,char b)
{
	int ii,jj;
	if(a<'A' || a>'Z')a='Z';
	ii=Blo_AA_Map_WS[Ori_AA_Map_WS[a-'A']];
	if(b<'A' || b>'Z')b='Z';
	jj=Blo_AA_Map_WS[Ori_AA_Map_WS[b-'A']];
	return Ori_BLOSUM_62[ii][jj];
}

//-------------- Single_Output -------------//
//-> protein structure -------//
char *TM_AMI1,*TM_AMI2;
char *TM_CLE1,*TM_CLE2;
int *TM_3DI1,*TM_3DI2;
XYZ *TM_MOL1,*TM_MOL2,*TM_MOL_TMP;
PDB_Residue *TM_PDB1,*TM_PDB2,*TM_PDB_TMP;
int TM_MOLN1,TM_MOLN2;
//-> protein alignment ------//
int *TM_ALIGNMENT;
int *TM_ALI1;
double *TM_ROTMAT;
//-> CB_additional ----//__110830__//
char *tmp_ami;
XYZ *tmp_mcb;
int *ws_reco;
XYZ *TM_MCB1,*TM_MCB2;
//-> weighted matrix ---//
int *wwa1;
int *wwa2;
double *wswei;

namespace
{
void InvertRigidTransform(const double *forward, double *inverse)
{
	inverse[0] = forward[0];
	inverse[1] = forward[3];
	inverse[2] = forward[6];
	inverse[3] = forward[1];
	inverse[4] = forward[4];
	inverse[5] = forward[7];
	inverse[6] = forward[2];
	inverse[7] = forward[5];
	inverse[8] = forward[8];
	inverse[9] = -(inverse[0] * forward[9] + inverse[1] * forward[10] +
		inverse[2] * forward[11]);
	inverse[10] = -(inverse[3] * forward[9] + inverse[4] * forward[10] +
		inverse[5] * forward[11]);
	inverse[11] = -(inverse[6] * forward[9] + inverse[7] * forward[10] +
		inverse[8] * forward[11]);
}

void InvertAlignment(const int *alignment, int input_length, int output_length,
	vector<int> &inverse)
{
	inverse.assign(output_length, -1);
	for (int second = 0; second < input_length; ++second) {
		const int first = alignment[second];
		if (first >= 0 && first < output_length) {
			inverse[first] = second;
		}
	}
}

void RestoreAlignRecord(Align_Record &record, int original_length1,
	int original_length2)
{
	if (record.alignment.empty()) {
		return;
	}
	vector<int> original_alignment;
	const int alignment_length = min(
		static_cast<int>(record.alignment.size()), original_length1);
	InvertAlignment(record.alignment.data(), alignment_length,
		original_length2, original_alignment);
	record.alignment.swap(original_alignment);
	if (record.rotmat.size() >= 12) {
		double original_rotmat[12];
		InvertRigidTransform(record.rotmat.data(), original_rotmat);
		copy(original_rotmat, original_rotmat + 12, record.rotmat.begin());
	}
}

// Present one canonical orientation to the alignment core: the shorter
// structure is always query/mol1 and alignment[target] stores query indices.
// The original command-line orientation is restored at the output boundary.
class CanonicalPairView
{
public:
	explicit CanonicalPairView(PairwiseAligner &aligner)
		: aligner_(aligner), swapped_(TM_MOLN1 > TM_MOLN2), restored_(false),
		  original_length1_(TM_MOLN1), original_length2_(TM_MOLN2)
	{
		if (swapped_) {
			SwapInputs();
		}
	}

	~CanonicalPairView()
	{
		Restore();
	}

	void RestoreResults(int *alignment)
	{
		if (!swapped_) {
			return;
		}
		vector<int> original_alignment;
		InvertAlignment(alignment, original_length1_, original_length2_,
			original_alignment);
		copy(original_alignment.begin(), original_alignment.end(), alignment);
		for (size_t i = 0; i < aligner_.Pairwise_align_tot.size(); ++i) {
			RestoreAlignRecord(aligner_.Pairwise_align_tot[i], original_length1_,
				original_length2_);
		}
		RestoreAlignRecord(aligner_.pipeline_MAX, original_length1_,
			original_length2_);
	}

	void Restore()
	{
		if (restored_) {
			return;
		}
		if (swapped_) {
			SwapInputs();
			aligner_.Fragment_Input_Func(TM_CLE1, TM_CLE2, TM_AMI1, TM_AMI2,
				TM_MOL1, TM_MOL2, TM_MOLN1, TM_MOLN2);
			aligner_.Fragment_3Di_Input(TM_3DI1, TM_3DI2);
			aligner_.TM_Align_Init(TM_MOLN1, TM_MOLN2);
		}
		restored_ = true;
	}

private:
	CanonicalPairView(const CanonicalPairView &);
	CanonicalPairView &operator=(const CanonicalPairView &);

	void SwapInputs()
	{
		swap(TM_AMI1, TM_AMI2);
		swap(TM_CLE1, TM_CLE2);
		swap(TM_3DI1, TM_3DI2);
		swap(TM_MOL1, TM_MOL2);
		swap(TM_PDB1, TM_PDB2);
		swap(TM_MCB1, TM_MCB2);
		swap(wwa1, wwa2);
		swap(TM_MOLN1, TM_MOLN2);
		swap(aligner_.mas1, aligner_.mas2);
		swap(aligner_.TM_cb1, aligner_.TM_cb2);
	}

	PairwiseAligner &aligner_;
	bool swapped_;
	bool restored_;
	int original_length1_;
	int original_length2_;
};
}

//-------- Create -------//
void TM_Align_Init_WS(int maxlen1,int maxlen2)
{
	//-- protein structure --//
	//-> local
	TM_AMI1=new char[maxlen1+1];
	TM_AMI2=new char[maxlen2+1];
	TM_CLE1=new char[maxlen1+1];
	TM_CLE2=new char[maxlen2+1];
	TM_3DI1=new int[maxlen1];
	TM_3DI2=new int[maxlen2];
	for(int i=0;i<maxlen1;i++)TM_3DI1[i]=Foldseek3Di::INVALID_STATE;
	for(int i=0;i<maxlen2;i++)TM_3DI2[i]=Foldseek3Di::INVALID_STATE;
	//-> global
	TM_MOL1=new XYZ[maxlen1];
	TM_MOL2=new XYZ[maxlen2];
	TM_MOL_TMP=new XYZ[maxlen1+maxlen2];
	TM_PDB1=new PDB_Residue[maxlen1];
	TM_PDB2=new PDB_Residue[maxlen2];
	TM_PDB_TMP=new PDB_Residue[maxlen1+maxlen2];
	//-- protein alignment --//
	TM_ALIGNMENT=new int[maxlen1+maxlen2];
	TM_ALI1=new int[maxlen1+maxlen2];
	TM_ROTMAT=new double[12];
	//---- CB_additional ----//__110830__//
	tmp_ami=new char[maxlen1+maxlen2+1];
	tmp_mcb=new XYZ[maxlen1+maxlen2];
	ws_reco=new int[maxlen1+maxlen2];
	TM_MCB1=new XYZ[maxlen1];
	TM_MCB2=new XYZ[maxlen2];
	//---- others ----//
	wwa1=new int[maxlen1];
	wwa2=new int[maxlen2];
	wswei=new double[maxlen1*maxlen2];
}
void TM_Align_dele_WS(void)
{
	//---- protein structure ---//
	//-> local 
	delete [] TM_AMI1;
	delete [] TM_AMI2;
	delete [] TM_CLE1;
	delete [] TM_CLE2;
	delete [] TM_3DI1;
	delete [] TM_3DI2;
	//-> global
	delete [] TM_MOL1;
	delete [] TM_MOL2;
	delete [] TM_MOL_TMP;
	delete [] TM_PDB1;
	delete [] TM_PDB2;
	delete [] TM_PDB_TMP;
	//---- protein alignment ---//
	delete [] TM_ALIGNMENT;
	delete [] TM_ALI1;
	delete [] TM_ROTMAT;
	//---- CB_additional ----//__110830__//
	delete [] tmp_ami;
	delete [] tmp_mcb;
	delete [] ws_reco;
	delete [] TM_MCB1;
	delete [] TM_MCB2;
	//---- others ----//
	delete [] wwa1;
	delete [] wwa2;
	delete [] wswei;
}

//=========== Fasta_Output related ===========//
void FASTA_Output(FILE *fp,string &nam1,string &nam2,char *ami1,char *ami2,int *ali1,int *ali2,int moln1,int moln2)
{
	//ali2->ali1
	int i;
	int ii,jj;
	for(i=0;i<moln1;i++)ali1[i]=-1;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]>=0)
		{
			ii=ali2[i];
			jj=i;
			ali1[ii]=jj;
		}
	}

	//ali1->alignment
	//init
	vector<pair<int,int> > alignment;
	alignment.clear();
	//start
	int j;
	int wlen;
	int pre_ii=0;
	int pre_jj=0;
	for(i=1;i<=moln1;i++)
	{
		ii=i;
		jj=ali1[i-1];  //ali1 starts from 0, correspondence also from 0
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
				alignment.push_back (pair<int,int>(pre_ii, -pre_jj)); //Ix
			}
			wlen=jj-pre_jj;
			for(j=1;j<wlen;j++)
			{
				pre_jj++;
				alignment.push_back (pair<int,int>(-pre_ii, pre_jj)); //Iy
			}
			//current_path
			alignment.push_back (pair<int,int>(ii, jj)); //Match
			//update
			pre_ii=ii;
			pre_jj=jj;
		}
	}
	//termi
	pre_ii++;
	for(i=pre_ii;i<=moln1;i++)alignment.push_back (pair<int,int>(i, -pre_jj)); //Ix
	pre_jj++;
	for(i=pre_jj;i<=moln2;i++)alignment.push_back (pair<int,int>(-moln1, i));  //Iy

	//alignment->output
	//output
	char c;
	int size=(int)alignment.size();
	fprintf(fp,">%s\n",nam1.c_str());
	for(i=0;i<size;i++)
	{
		ii=alignment[i].first;
		if(ii<=0)fprintf(fp,"-");
		else
		{
			c=ami1[ii-1];
			fprintf(fp,"%c",c);
//			if(c=='X'||c=='Z')fprintf(fp,".");
//			else fprintf(fp,"%c",c);
		}
	}
	fprintf(fp,"\n");
	fprintf(fp,">%s\n",nam2.c_str());
	for(i=0;i<size;i++)
	{
		jj=alignment[i].second;
		if(jj<=0)fprintf(fp,"-");
		else
		{
			c=ami2[jj-1];
			fprintf(fp,"%c",c);
//			if(c=='X'||c=='Z')fprintf(fp,".");
//			else fprintf(fp,"%c",c);
		}
	}
	fprintf(fp,"\n");
}
void FASTA_Output_More(string &ws_output_tot,string &nam1_,string &nam2_,
		char *ami1_,char *ami2_,char *cle1_,char *cle2_,
		const int *three_di1,const int *three_di2,XYZ *mol1,XYZ *mol2,
		int *ali1,int *ali2,int moln1,int moln2,double d0,int &seqid,
		int &blosum,int &three_di_score,vector <double> &match_wei,
		vector <double> &rmsd_out,int CLEorNOT,int threscut=80)
{
	ws_output_tot="";
	//ali2->ali1
	int i;
	int ii,jj;
	for(i=0;i<moln1;i++)ali1[i]=-1;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]>=0)
		{
			ii=ali2[i];
			jj=i;
			ali1[ii]=jj;
		}
	}

	//ali1->alignment
	//init
	vector<pair<int,int> > alignment;
	alignment.clear();
	//start
	int j;
	int wlen;
	int pre_ii=0;
	int pre_jj=0;
	for(i=1;i<=moln1;i++)
	{
		ii=i;
		jj=ali1[i-1];  //ali1 starts from 0, correspondence also from 0
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
				alignment.push_back (pair<int,int>(pre_ii, -pre_jj)); //Ix
			}
			wlen=jj-pre_jj;
			for(j=1;j<wlen;j++)
			{
				pre_jj++;
				alignment.push_back (pair<int,int>(-pre_ii, pre_jj)); //Iy
			}
			//current_path
			alignment.push_back (pair<int,int>(ii, jj)); //Match
			//update
			pre_ii=ii;
			pre_jj=jj;
		}
	}
	//termi
	pre_ii++;
	for(i=pre_ii;i<=moln1;i++)alignment.push_back (pair<int,int>(i, -pre_jj)); //Ix
	pre_jj++;
	for(i=pre_jj;i<=moln2;i++)alignment.push_back (pair<int,int>(-moln1, i));  //Iy


	//-------- get start and end ---------//
	int start1,start2;
	int end1,end2;
	int wsstart=-1;
	int wsend=-1;
	int x=0;
	int y=0;
	rmsd_out.clear();
	int size=(int)alignment.size();
	{
		//omit head gap
		for(i=0;i<size;i++)
		{
			x = alignment[i].first;
			y = alignment[i].second;
			wsstart=i;
			if(x>0 && y>0)break;
			rmsd_out.push_back(-1);
		}
		if(x==0 || y==0)return;
		x--;
		y--;
		start1=x;
		start2=y;
		end1=x;
		end2=y;
		//omit tail gap
		for(i=size-1;i>=0;i--)
		{
			x = alignment[i].first;
			y = alignment[i].second;
			wsend=i;
			if(x>0 && y>0)break;
		}
	}

	//---------- calc -----------//
	{
		//get name
		string nam1;
		if(nam1_.length()<14)nam1=nam1_;
		else nam1=nam1_.substr(0,14);
		string nam2;
		if(nam2_.length()<14)nam2=nam2_;
		else nam2=nam2_.substr(0,14);
		//others
		double wscur1,wscur2,wscur;
		int first=1;
		int lali=0;
		char wstmp;
		string cle1,ami1,cle2,ami2;
		string score_out;
		//string clear
		cle1.clear();
		ami1.clear();
		cle2.clear();
		ami2.clear();
		score_out.clear();
		//output related
		char ws_output[300000];
		string ws_output_str;
		ws_output_str.clear();
		//[real process]
		int num=0;
		blosum=0;
		three_di_score=0;
		seqid=0;
		match_wei.clear();
		for(i=wsstart; i <= wsend; i++)
		{
			x = alignment[i].first;
			y = alignment[i].second;
			x--,y--;
			//output
			if(num%threscut==0 && first==0)
			{
				//fragment output
				ws_output_str.clear();
				sprintf(ws_output,"\n");
				ws_output_str+=ws_output;
				if(CLEorNOT==1)
				{
					sprintf(ws_output,"T cle                 %s           \n",cle1.c_str());
					ws_output_str+=ws_output;
				}
				sprintf(ws_output,"T %-14s %4d %s %4d (%d)\n",
					nam1.c_str(),start1+1,ami1.c_str(),end1+1,moln1);
				ws_output_str+=ws_output;
				sprintf(ws_output,"RMSD                  %s           \n",score_out.c_str());
				ws_output_str+=ws_output;
				sprintf(ws_output,"S %-14s %4d %s %4d (%d)\n",
					nam2.c_str(),start2+1,ami2.c_str(),end2+1,moln2);
				ws_output_str+=ws_output;
				if(CLEorNOT==1)
				{
					sprintf(ws_output,"S cle                 %s           \n",cle2.c_str());
					ws_output_str+=ws_output;
				}
				sprintf(ws_output,"\n");
				ws_output_str+=ws_output;
				//string clear
				cle1.clear();
				ami1.clear();
				cle2.clear();
				ami2.clear();
				score_out.clear();
				//position update
				start1=end1+1;
				start2=end2+1;
				//final output
				ws_output_tot+=ws_output_str;
			}
			num++;
			//lali
			if(x>=0 && y>=0)lali++;
			if(x>=0)end1=x;
			if(y>=0)end2=y;
			first=0;
			//record
			//[temp]
			if(x>=0)
			{
				cle1=cle1+cle1_[x];
				ami1=ami1+ami1_[x];
			}
			else
			{
				cle1=cle1+'-';
				ami1=ami1+'-';
			}
			//[targ]
			if(y>=0)
			{
				cle2=cle2+cle2_[y];
				ami2=ami2+ami2_[y];
			}
			else
			{
				cle2=cle2+'-';
				ami2=ami2+'-';
			}
			//[score]
			if(x>=0 && y>=0)
			{
				//calculate local RMSD
				double dist2=mol1[x].distance_square(mol2[y]);
//				double tmscore=1.0/(1.0+dist2/d0/d0);
				double tmscore=sqrt(1.0*dist2);
				int wssco=(int)(tmscore);
				if(wssco<0)wssco=0;
				if(wssco>9)wssco=9;
				wstmp=wssco+'0';
				if(wstmp=='9')wstmp='.';
				score_out=score_out+wstmp;
				rmsd_out.push_back(tmscore);
				// calculate amino-acid and 3Di evidence
				wscur1=BLOSUM62_Calc(ami1_[x],ami2_[y]);
				wscur2=Foldseek3Di::Score(three_di1[x],three_di2[y]);
				blosum+=(int)wscur1;
				three_di_score+=(int)wscur2;
				//calculate sequence identity
				if(ami1_[x]==ami2_[y])
				{
					if(ami1_[x]!='.' && ami1_[x]!='X' && ami1_[x]!='Z')seqid++;
				}
				//add to match_wei
				wscur=Foldseek3Di::CombinedEvidence(wscur1,wscur2);
				match_wei.push_back(wscur);
			}
			else
			{
				score_out=score_out+" ";
				rmsd_out.push_back(-1);
			}
		}

		//rmsd_termi
		for(i=wsend+1;i<size;i++)rmsd_out.push_back(-1);

		//termi_iprocess
	//	if(i%threscut==0 && first==0)
		{
			//fragment output
			ws_output_str.clear();
			sprintf(ws_output,"\n");
			ws_output_str+=ws_output;
			if(CLEorNOT==1)
			{
				sprintf(ws_output,"T cle                 %s           \n",cle1.c_str());
				ws_output_str+=ws_output;
			}
			sprintf(ws_output,"T %-14s %4d %s %4d (%d)\n",
				nam1.c_str(),start1+1,ami1.c_str(),end1+1,moln1);
			ws_output_str+=ws_output;
			sprintf(ws_output,"RMSD                  %s           \n",score_out.c_str());
			ws_output_str+=ws_output;
			sprintf(ws_output,"S %-14s %4d %s %4d (%d)\n",
				nam2.c_str(),start2+1,ami2.c_str(),end2+1,moln2);
			ws_output_str+=ws_output;
			if(CLEorNOT==1)
			{
				sprintf(ws_output,"S cle                 %s           \n",cle2.c_str());
				ws_output_str+=ws_output;
			}
			sprintf(ws_output,"\n");
			ws_output_str+=ws_output;
			//string clear
			cle1.clear();
			ami1.clear();
			cle2.clear();
			ami2.clear();
			score_out.clear();
			//position update
			start1=end1+1;
			start2=end2+1;
			//final output
			ws_output_tot+=ws_output_str;
		}
	}
}

//============= Weight_Matrix_Built ===============//__110815__//
void Linear_Transform(double *in_mat,int moln1,int moln2,double *out_mat)
{
	int i;
	int totnum=moln1*moln2;
	double wscur;
	double wsmax,wsmin;
	wsmax=-99999;
	wsmin=99999;
	for(i=0;i<totnum;i++)
	{
		wscur=in_mat[i];
		if(wscur>wsmax)wsmax=wscur;
		if(wscur<wsmin)wsmin=wscur;
	}
	for(i=0;i<totnum;i++)
	{
		wscur=in_mat[i];
		out_mat[i]=1.0*(wscur-wsmin)/(wsmax-wsmin);
	}
}
void Weight_Matrix_Built_WS(char *ami1,char *cle1,const int *three_di1,
	char *ami2,char *cle2,const int *three_di2,double *wei_out)
{
	SequenceScorer bioinfo_code(2);
	const int w1=(int)strlen(ami1);
	const int w2=(int)strlen(ami2);

	//AMI only
	bioinfo_code.AMI_transform(ami1,wwa1);
	bioinfo_code.AMI_transform(ami2,wwa2);
#if EALIGN_ALIGNMENT_EVIDENCE_MODE == EALIGN_EVIDENCE_AA_AND_CLESUM
	vector<int> cle_code1(w1);
	vector<int> cle_code2(w2);
	bioinfo_code.CLE_transform(cle1,cle_code1.data());
	bioinfo_code.CLE_transform(cle2,cle_code2.data());
	(void)three_di1;
	(void)three_di2;
#else
	(void)cle1;
	(void)cle2;
#endif
	//gen WEI
	int i,j;
	int ii,jj;
	double wscur;
	double wscur1,wscur2;

	//calc
	for(i=0;i<w1;i++)
	{
		for(j=0;j<w2;j++)
		{
			ii=i;
			jj=j;
			wscur1=bioinfo_code.Universal_Calc(ii,jj,1,wwa1,wwa2,w1,w2,3);
#if EALIGN_ALIGNMENT_EVIDENCE_MODE == EALIGN_EVIDENCE_AA_AND_CLESUM
			wscur2=bioinfo_code.Universal_Calc(ii,jj,1,
				cle_code1.data(),cle_code2.data(),w1,w2,1);
#else
			wscur2=Foldseek3Di::Score(three_di1[ii],three_di2[jj]);
#endif

#if EALIGN_ALIGNMENT_EVIDENCE_MODE == EALIGN_EVIDENCE_3DI_ONLY
			wscur=wscur2;
#elif EALIGN_ALIGNMENT_EVIDENCE_MODE == EALIGN_EVIDENCE_AA_ONLY
			wscur=wscur1;
#elif EALIGN_ALIGNMENT_EVIDENCE_MODE == EALIGN_EVIDENCE_AA_AND_CLESUM
			wscur=wscur1+wscur2;
#else
			wscur=Foldseek3Di::CombinedEvidence(wscur1,wscur2);
#endif
			wei_out[i*w2+j]=wscur;

		}
	}
}

//---------------- superimpose_fullatom --------------------//
void Superimpose_FullAtom(Kabsch &kabsch,PDB_Residue *in,int moln,PDB_Residue *out,double *rotmat)
{
	int i,k;
	int num;
	PDB_Residue pdb;
	XYZ xyz;
	int numb;
	double R,tmpr;
	for(i=0;i<moln;i++)
	{
		pdb=in[i];
		//rotate backbone
		num=pdb.get_backbone_totnum();
		for(k=0;k<num;k++)
		{
			if(pdb.get_backbone_part_index(k)==0)continue;
			pdb.get_backbone_atom(k,xyz,numb,R,tmpr);
			kabsch.rot_point(xyz,xyz,rotmat);
			pdb.set_backbone_atom(k,xyz,numb,R,tmpr);
		}
		//rotate sidechain
		num=pdb.get_sidechain_totnum();
		for(k=0;k<num;k++)
		{
			if(pdb.get_sidechain_part_index(k)==0)continue;
			pdb.get_sidechain_atom(k,xyz,numb,R,tmpr);
			kabsch.rot_point(xyz,xyz,rotmat);
			pdb.set_sidechain_atom(k,xyz,numb,R,tmpr);
		}
		out[i]=pdb;
	}
}

//------------- generate CB --------------//
void Generate_CB(BetaCarbonBuilder &confo_beta,int moln,PDB_Residue *pdb,char *ami,char *cle,XYZ *mol,XYZ *mcb)
{
	int k;
	int ws_correct; //default: OK
	PDB_Residue PDB;
	double dist;
	//mol1 process
	ws_correct=1;
	for(k=0;k<moln;k++)
	{
		PDB=pdb[k];
		if(PDB.get_sidechain_part_index(0)==0)
		{
			ws_correct=0;
			ws_reco[k]=0;
		}
		else
		{
			PDB.get_sidechain_atom(0,mcb[k]);
			dist=mcb[k].distance_square(mol[k]);
			if(dist<1.0||dist>4.0)
			{
				ws_correct=0;
				ws_reco[k]=0;
			}
			else ws_reco[k]=1;
		}
	}
	if(ws_correct==0)
	{
		strcpy(tmp_ami,ami);
		for(k=0;k<moln;k++)if(tmp_ami[k]=='G')tmp_ami[k]='A';
		confo_beta.Recon_Beta_21(mol,tmp_mcb,moln,tmp_ami,cle);
		for(k=0;k<moln;k++)if(ws_reco[k]==0)mcb[k]=tmp_mcb[k];
	}
}

void Generate_3Di(PDB_Residue *pdb,XYZ *ca,int moln,int *three_di)
{
	Foldseek3Di::Encode(pdb,ca,moln,three_di);
}

// Run the one alignment mode used by run_full_alignment_benchmark_copy.py.
static int Run_Alignment(PairwiseAligner &aligner, string &name1,
	string &name2, const string &output_base, string &screen_output)
{
	CanonicalPairView pair_view(aligner);
	#if EALIGN_ALIGNMENT_EVIDENCE_MODE == EALIGN_NO_EVIDENCE
		aligner.TMali_Weight = nullptr;  // E_h = V(i,j|R) * g(d_ij(T))
	#else
		Weight_Matrix_Built_WS(
			TM_AMI1, TM_CLE1, TM_3DI1,
			TM_AMI2, TM_CLE2, TM_3DI2, wswei);
		Linear_Transform(wswei, TM_MOLN1, TM_MOLN2, wswei);
		aligner.TMali_Weight = wswei;
	#endif
	aligner.TM_Vect_Score = EALIGN_ALIGNMENT_GEOMETRY_MODE;
	aligner.Pairwise_Kill_Gap = 1;
	aligner.Kill_Frag = 1;

	aligner.Fragment_Input_Func(TM_CLE1, TM_CLE2, TM_AMI1, TM_AMI2,
		TM_MOL1, TM_MOL2, TM_MOLN1, TM_MOLN2);
	aligner.Fragment_3Di_Input(TM_3DI1, TM_3DI2);
	aligner.TM_Align_Init(TM_MOLN1, TM_MOLN2);

	const int norm_len = min(TM_MOLN1, TM_MOLN2);
	const double norm_d0 = aligner.Calc_TM_d0_Simp(norm_len);
	aligner.TM_DIST_CUT = 1;
	const double distance_cutoff = aligner.d8;

	double alignment_score = aligner.Pairwise_Align_Normal(
		TM_MOL1, TM_MOL2, TM_MOLN1, TM_MOLN2, TM_ALIGNMENT,
		norm_len, norm_d0);
	if (alignment_score == 0) {
		alignment_score = aligner.Pairwise_Align_Total(
			TM_MOL1, TM_MOL2, TM_MOLN1, TM_MOLN2, TM_ALIGNMENT,
			norm_len, norm_d0);
	}
	if (alignment_score == 0) {
		pair_view.Restore();
		char line[1024];
		snprintf(line, sizeof(line),
			"%s %s %4d %4d -> %6d %6d %9.2f -> %4d %7.3f %7.3f -> %6.3f %6.3f %6.3f -> %5d %5d %6.3f %.3f\n",
			name1.c_str(), name2.c_str(), TM_MOLN1, TM_MOLN2,
			0, 0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,
			norm_len, distance_cutoff, 0.0);
		screen_output = line;
		return 0;
	}

	pair_view.RestoreResults(TM_ALIGNMENT);
	pair_view.Restore();

	for (int i = 0; i < TM_MOLN2; ++i) {
		TM_ALIGNMENT[i] = aligner.Pairwise_align_tot[0].alignment[i];
	}
	for (int i = 0; i < 12; ++i) {
		TM_ROTMAT[i] = aligner.Pairwise_align_tot[0].rotmat[i];
	}

	Kabsch kabsch;
	kabsch.rot_mol(TM_MOL1, TM_MOL_TMP, TM_MOLN1, TM_ROTMAT);

	int seqid = 0;
	int blosum = 0;
	int three_di_score = 0;
	vector<double> match_wei;
	vector<double> rmsd_out;
	string final_record;
	FASTA_Output_More(final_record, name1, name2,
		TM_AMI1, TM_AMI2, TM_CLE1, TM_CLE2,
		TM_3DI1, TM_3DI2, TM_MOL_TMP, TM_MOL2, TM_ALI1, TM_ALIGNMENT,
		TM_MOLN1, TM_MOLN2, aligner.d0, seqid, blosum, three_di_score,
		match_wei, rmsd_out, 0);

	string filename;
	FILE *fp = 0;

	filename = "./" + output_base + ".fasta";
	fp = fopen(filename.c_str(), "wb");
	if (fp == 0) {
		fprintf(stderr, "ERROR: file %s can't be opened. \n", filename.c_str());
	} else {
		FASTA_Output(fp, name1, name2,
			TM_AMI1, TM_AMI2, TM_ALI1, TM_ALIGNMENT, TM_MOLN1, TM_MOLN2);
		fclose(fp);
	}

	filename = "./" + output_base + ".pdb";
	fp = fopen(filename.c_str(), "wb");
	if (fp == 0) {
		fprintf(stderr, "ERROR: file %s can't be opened. \n", filename.c_str());
	} else {
		Superimpose_FullAtom(kabsch, TM_PDB1, TM_MOLN1, TM_PDB_TMP, TM_ROTMAT);
		mol_out.Output_PDB_III(fp, TM_MOLN1, TM_PDB_TMP, 'A', 1);
		fprintf(fp, "TER                                                                             \n");
		mol_out.Output_PDB_III(fp, TM_MOLN2, TM_PDB2, 'B', 1);
		fprintf(fp, "TER                                                                             \n");
		fprintf(fp, "END                                                                             \n");
		fclose(fp);
	}

	filename = "./" + output_base + ".local";
	fp = fopen(filename.c_str(), "wb");
	if (fp == 0) {
		fprintf(stderr, "ERROR: file %s can't be opened. \n", filename.c_str());
	} else {
		string local_output;
		FASTA_Output_More(local_output, name1, name2,
			TM_AMI1, TM_AMI2, TM_CLE1, TM_CLE2,
			TM_3DI1, TM_3DI2, TM_MOL_TMP, TM_MOL2, TM_ALI1, TM_ALIGNMENT,
			TM_MOLN1, TM_MOLN2, aligner.d0, seqid, blosum, three_di_score,
			match_wei, rmsd_out, 1);
		fprintf(fp, "%s", local_output.c_str());
		fclose(fp);
	}

	double score_parts[8];
	double rms = 0;
	int aligned_length = 0;
	const double tm_score = aligner.TM_Align_TM_Score(
		TM_MOL1, TM_MOL2, TM_MOLN1, TM_MOLN2, TM_ALIGNMENT,
		norm_len, norm_d0, rms, aligned_length, score_parts);
	const double combined_score = aligner.TM_Align_Get_Score_Simp_MatchWei(
		TM_MOL1, TM_MOL2, TM_ROTMAT, TM_MOLN1, TM_MOLN2,
		TM_ALIGNMENT, match_wei);
	const double gdt_ts = (score_parts[1] + score_parts[2] +
		score_parts[3] + score_parts[4]) / (4.0 * norm_len);
	const double gdt_ha = (score_parts[0] + score_parts[1] +
		score_parts[2] + score_parts[3]) / (4.0 * norm_len);
	const double ugdt = 0.25 * (score_parts[1] + score_parts[2] +
		score_parts[3] + score_parts[4]);
	const double maxsub = score_parts[5] / norm_len;

	filename = "./" + output_base + ".score";
	fp = fopen(filename.c_str(), "wb");
	if (fp == 0) {
		fprintf(stderr, "ERROR: file %s can't be opened. \n", filename.c_str());
	} else {
		fprintf(fp, "#name1 name2 len1 len2 -> BLOSUM 3DiScore CombinedScore -> LALI RMSDval TMscore -> MAXSUB GDT_TS GDT_HA -> SeqID nLen dCutoff uGDT\n");
		fprintf(fp, " %s %s %4d %4d -> %6d %6d %9.2f -> %4d %7.3f %7.3f -> %6.3f %6.3f %6.3f -> %5d %5d %6.3f %.3f\n",
			name1.c_str(), name2.c_str(), TM_MOLN1, TM_MOLN2,
			blosum, three_di_score, combined_score, aligned_length, rms, tm_score,
			maxsub, gdt_ts, gdt_ha, seqid, norm_len, distance_cutoff, ugdt);
		fprintf(fp, "#---------------- transformation to superpose 1st structure onto the 2nd ------------------------\n");
		fprintf(fp, "%9.6f %9.6f %9.6f %12.6f\n", TM_ROTMAT[0], TM_ROTMAT[1], TM_ROTMAT[2], TM_ROTMAT[9]);
		fprintf(fp, "%9.6f %9.6f %9.6f %12.6f\n", TM_ROTMAT[3], TM_ROTMAT[4], TM_ROTMAT[5], TM_ROTMAT[10]);
		fprintf(fp, "%9.6f %9.6f %9.6f %12.6f\n", TM_ROTMAT[6], TM_ROTMAT[7], TM_ROTMAT[8], TM_ROTMAT[11]);
		fclose(fp);
	}

	char line[4096];
	string output;
	snprintf(line, sizeof(line), "#------------------------------------------------------#\n");
	output += line;
	snprintf(line, sizeof(line), "#                       EAlign                         #\n");
	output += line;
	snprintf(line, sizeof(line), "#        Pairwise protein structure alignment          #\n");
	output += line;
	snprintf(line, sizeof(line), "#------------------------------------------------------#\n\n");
	output += line;
	snprintf(line, sizeof(line), "1st input protein: %s    length= %4d \n", name1.c_str(), TM_MOLN1);
	output += line;
	snprintf(line, sizeof(line), "2nd input protein: %s    length= %4d \n\n", name2.c_str(), TM_MOLN2);
	output += line;
	snprintf(line, sizeof(line), "#----- Normalization length for TMscore, MAXSUB, GDT_TS, and GDT_HA is %4d\n", norm_len);
	output += line;
	snprintf(line, sizeof(line), "#----- Distance cutoff value is %lf\n", distance_cutoff);
	output += line;
	snprintf(line, sizeof(line), "# BLOSUM 3DiScore CombinedScore SeqID LALI RMSD(A) TMscore MAXSUB GDT_TS GDT_HA uGDT\n");
	output += line;
	snprintf(line, sizeof(line), "  %6d %6d %9.2f %5d %4d %7.3f %7.3f %6.3f %6.3f %6.3f %.3f\n",
		blosum, three_di_score, combined_score, seqid, aligned_length, rms,
		tm_score, maxsub, gdt_ts, gdt_ha, ugdt);
	output += line;
	snprintf(line, sizeof(line), "#----- Transformation to superpose the 1st protein onto the 2nd ---#\n");
	output += line;
	snprintf(line, sizeof(line), "# i          t(i)         u(i,1)         u(i,2)         u(i,3) \n");
	output += line;
	snprintf(line, sizeof(line), "  1      %12.6f   %12.9f   %12.9f   %12.9f\n", TM_ROTMAT[9], TM_ROTMAT[0], TM_ROTMAT[1], TM_ROTMAT[2]);
	output += line;
	snprintf(line, sizeof(line), "  2      %12.6f   %12.9f   %12.9f   %12.9f\n", TM_ROTMAT[10], TM_ROTMAT[3], TM_ROTMAT[4], TM_ROTMAT[5]);
	output += line;
	snprintf(line, sizeof(line), "  3      %12.6f   %12.9f   %12.9f   %12.9f\n\n", TM_ROTMAT[11], TM_ROTMAT[6], TM_ROTMAT[7], TM_ROTMAT[8]);
	output += line;
	output += final_record;
	screen_output = output;
	return 1;
}

static void Usage(const char *program)
{
	fprintf(stderr, "Usage: %s protein_1 protein_2 -p 1 -o output_base\n", program);
}

int main(int argc, char **argv)
{
	if (argc != 7 || strcmp(argv[3], "-p") != 0 || strcmp(argv[4], "1") != 0 ||
		strcmp(argv[5], "-o") != 0 || argv[6][0] == '\0') {
		Usage(argv[0]);
		return 2;
	}

	string pdb1 = argv[1];
	string pdb2 = argv[2];
	const string output_base = argv[6];
	string range = "_";
	StructureLoader input1;
	StructureLoader input2;

	if (input1.XYZ_Input(pdb1, range, 0, TM_MOLN1, 0, 0, 0, 0, 0) != 1) {
		fprintf(stderr, "FILE1 NOT FOUND!!!![%s]\n", pdb1.c_str());
		return -1;
	}
	if (TM_MOLN1 <= 0) {
		fprintf(stderr, "FILE1 IS EMPTY!!!![%s]\n", pdb1.c_str());
		return -1;
	}
	if (input2.XYZ_Input(pdb2, range, 0, TM_MOLN2, 0, 0, 0, 0, 0) != 1) {
		fprintf(stderr, "FILE2 NOT FOUND!!!![%s]\n", pdb2.c_str());
		return -1;
	}
	if (TM_MOLN2 <= 0) {
		fprintf(stderr, "FILE2 IS EMPTY!!!![%s]\n", pdb2.c_str());
		return -1;
	}

	TM_Align_Init_WS(TM_MOLN1, TM_MOLN2);

	const int preload1 = input1.PRE_LOAD;
	const int warning1 = input1.WARNING_out;
	input1.PRE_LOAD = 1;
	input1.WARNING_out = 0;
	input1.XYZ_Input(pdb1, range, 0, TM_MOLN1, TM_MOL1, TM_AMI1, TM_CLE1, 0, TM_PDB1);
	input1.PRE_LOAD = preload1;
	input1.WARNING_out = warning1;

	const int preload2 = input2.PRE_LOAD;
	const int warning2 = input2.WARNING_out;
	input2.PRE_LOAD = 1;
	input2.WARNING_out = 0;
	input2.XYZ_Input(pdb2, range, 0, TM_MOLN2, TM_MOL2, TM_AMI2, TM_CLE2, 0, TM_PDB2);
	input2.PRE_LOAD = preload2;
	input2.WARNING_out = warning2;

	string name1;
	string name2;
	getBaseName(pdb1, name1, '/', '.');
	getBaseName(pdb2, name2, '/', '.');

	const int max_length = max(TM_MOLN1, TM_MOLN2);
	BetaCarbonBuilder confo_beta(max_length);
	Generate_CB(confo_beta, TM_MOLN1, TM_PDB1, TM_AMI1, TM_CLE1, TM_MOL1, TM_MCB1);
	Generate_CB(confo_beta, TM_MOLN2, TM_PDB2, TM_AMI2, TM_CLE2, TM_MOL2, TM_MCB2);
	Generate_3Di(TM_PDB1, TM_MOL1, TM_MOLN1, TM_3DI1);
	Generate_3Di(TM_PDB2, TM_MOL2, TM_MOLN2, TM_3DI2);

	PairwiseAligner aligner(max_length);
	aligner.Geo_Sparse_DP = 1;
	aligner.SparseDP_filter = 0;
	aligner.TM_cb1 = TM_MCB1;
	aligner.TM_cb2 = TM_MCB2;

	string screen_output;
	Run_Alignment(aligner, name1, name2, output_base, screen_output);
	printf("%s", screen_output.c_str());
	TM_Align_dele_WS();
	return 0;
}

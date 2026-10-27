#include "BetaCarbonBuilder.h"

double CB_CLE_Type[17][2]=
{
	{0.9952,-2.0677},
	{1.0543,-2.3374},
	{1.0045,-2.4657},
	{1.0746,-2.5422},
	{1.1064,-2.5927},
	{1.1015,-2.5619},
	{1.0282,-2.5262},
	{0.9489,-1.9581},
	{0.9733,-2.0332},
	{0.9628,-2.0315},
	{1.0555,-2.1765},
	{1.1752,-2.5056},
	{0.9995,-2.1944},
	{1.1792,-2.2233},
	{1.2871,-2.2675},
	{1.0195,-2.1373},
	{1.1023,-2.2232}
};
// first->bend, second->tort

//--------------constructor--------------------//
BetaCarbonBuilder::BetaCarbonBuilder(int num)
{
	InitializeBetaCarbonBuffers(num);
}
BetaCarbonBuilder::~BetaCarbonBuilder(void)
{
	ReleaseBetaCarbonBuffers();
}

//------------ init ------------//
void BetaCarbonBuilder::InitializeBetaCarbonBuffers(int totlen)
{
	//temp_structure
	beta_dist=new double[totlen];
	beta_bend=new double[totlen];
	beta_tort=new double[totlen];
}
void BetaCarbonBuilder::ReleaseBetaCarbonBuffers(void)
{
	//temp_structure
	delete [] beta_dist;
	delete [] beta_bend;
	delete [] beta_tort;
}

//--------------function------------//
//given CA+CB -> return dist,bend,tort
//given CA and dist,bend,tort -> return CB (Zheng & Liu's method)
void BetaCarbonBuilder::BuildFromAngles(XYZ *CA,XYZ *CB,int moln,double *bend,double *tort,double *dist)
{
	int i,j;
	XYZ xyz;
	double beta,radii;
	//process
	for(i=1;i<moln-1;i++)
	{
		//init
		if(dist==NULL)beta=1.53;
		else beta=dist[i];
		if(beta<0.0)
		{
			CB[i]=CA[i];
			continue;
		}
		//calc
		xyz=(CA[i]-CA[i-1]);
		xyz.xyz2double(beta_v1);
		xyz=(CA[i+1]-CA[i]);
		xyz.xyz2double(beta_v2);
		cross(beta_y,beta_v1,beta_v2);
		Universal_Rotation(beta_y,bend[i],beta_rotation);
		Vector_Multiply(beta_x1,beta_rotation,beta_v1);
		Universal_Rotation(beta_v1,tort[i],beta_rotation);
		Vector_Multiply(beta_x2,beta_rotation,beta_x1);
		//evaluate
		radii=sqrt(dot(beta_x2,beta_x2));
		CA[i].xyz2double(beta_ca);
		for(j=0;j<3;j++)beta_cb[j]=beta_ca[j]+beta*beta_x2[j]/radii;
		CB[i].double2xyz(beta_cb);
	}
	//assign head_tail //this might be modified in the future...
	if(moln>1)
	{
		CB[0]=CB[1]-CA[1]+CA[0];
		CB[moln-1]=CB[moln-2]-CA[moln-2]+CA[moln-1];
	}
	else
	{
		CB[0]=CA[0];
		CB[moln-1]=CA[moln-1];
	}
}
//given CA -> return CB (Levitt's method)

//===================== Recon_Beta ==================//
//given CA+CB, and AMI_Dist, update SC
//given CA, return CB (Levitt)
//given CA, and CLE_Type, return CB (Liu & Wang)
void BetaCarbonBuilder::Recon_Beta_21(XYZ *CA,XYZ *CB,int moln,char *ami,char *cle)  //recon CB
{
	int i;
	char code;
	if(moln<=0)return;
	//init
	for(i=0;i<moln;i++)
	{
		//dist
		if(ami[i]=='G')beta_dist[i]=0.0;
		else beta_dist[i]=1.53;
		//bend+tort
		code=cle[i];
		if(code=='R')code='Q'; //use the most random one
		beta_bend[i]=CB_CLE_Type[code-'A'][0];
		beta_tort[i]=CB_CLE_Type[code-'A'][1];
	}	
	//process
	BuildFromAngles(CA,CB,moln,beta_bend,beta_tort,beta_dist);
}

#pragma once
#include "XYZ.h"
#include "GeometryUtils.h"

//====class: BetaCarbonBuilder====//
//=> CA<->CB generating
class BetaCarbonBuilder
{
public:
	BetaCarbonBuilder(int num=PROT_MAX_NUM);
	~BetaCarbonBuilder(void);

//--- variable ---//
public:
	//normal_variable
	double beta_v1[3],beta_v2[3];
	double beta_x1[3],beta_x2[3];
	double beta_y[3],beta_ca[3],beta_cb[3];
	double beta_rotation[3][3];      //temp tormat;
	//temp_variable
	double *beta_dist;
	double *beta_bend;
	double *beta_tort;
	//main_variable

//--- function ---//
public:
	//init
	void InitializeBetaCarbonBuffers(int totlen);
	void ReleaseBetaCarbonBuffers(void);
	//process
	void BuildFromAngles(XYZ *CA,XYZ *CB,int moln,double *bend,double *tort,double *dist);
	void Recon_Beta_21(XYZ *CA,XYZ *CB,int moln,char *ami,char *cle); //given CA, and CLE_Type, return CB (Liu & Wang)
};

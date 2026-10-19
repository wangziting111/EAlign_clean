#include "StructureWriter.h"

//--------------------- Start -----------------//
StructureWriter::StructureWriter(void)
{
	OutPrt=0;       // whether output printf    // (default:no)
}
StructureWriter::~StructureWriter(void)
{
}

//============================ Output_Related ============================//
//---- neo PDB_Residue output ---//__110230__//
//OutType [0: sequential][1: PDB numbering]
//OutMode [-2:CA+CB,-1:CA,0:NCaC+CB,+1:ALL]
//OutGlys [0: GLY 'CB'='CA'][1: normal GLY]
int StructureWriter::Output_PDB_III(FILE *fp,int len,PDB_Residue *mol,char Chain_ID,
	int OutType,int OutMode,int OutGlys,int OutLast) //output full-atom PDB file
{
	int i,k;
	string TER="TER                                                                             ";
	string buf="              ";
	int number;
	char amino;
	const char *dummyaa;
	const char *atomname;
	string pdbind_;
	string pdbind;
	double x,y,z;
	double rfactor,temperature;
	int numb;
	int outtype_atom=1;
	char Ini_Chain;
	char Real_Chain;
	PDB_Residue PDB;
	XYZ xyz;
	int cur_pos;
	char cur_chain;
	char nxt_chain;
	//ws_new//__110430__//
	char output[100];

	//judge
	if(mol==NULL)return -1;
	Ini_Chain=Chain_ID;
	if(Ini_Chain=='!'||Ini_Chain==-1||Ini_Chain=='_')Ini_Chain='_'; //put the chain from pdbind//__110230__//
	//process
	cur_pos=0;
	for(i=0;i<len;i++)
	{
		//init
		PDB=mol[i];
		amino=PDB.get_AA();
		dummyaa=One2Three_III(amino);
		PDB.get_PDB_residue_number(pdbind_);
		pdbind=pdbind_.substr(1,5);
		if(Ini_Chain=='_')Real_Chain=pdbind_[0];
		else Real_Chain=Ini_Chain;
		cur_chain=pdbind_[0];
		//backbone_out
		number=PDB.get_backbone_totnum();  //this value should be 4
		for(k=0;k<number;k++)
		{
			//check
			if(OutMode<0) //only output CA or CB
			{
				if(k!=1)continue;
			}
			//get_backbone
			if(PDB.get_backbone_part_index(k)==0)
			{
				if(OutPrt==1)fprintf(stderr,"BackBone_Bad_Output[chain;%c][pos:%d][atom:%d]\n!!",Real_Chain,i+1,k+1);
				continue;
			}
			atomname=backbone_atom_name_decode(k);
			PDB.get_backbone_atom(k,xyz, numb, rfactor, temperature);
			x=xyz.X;
			y=xyz.Y;
			z=xyz.Z;
			//output
			if(OutType==1) // PDB numbering
			{
				sprintf(output,"ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
					numb,atomname,dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
			}
			else           // sequential numbering
			{
				if(OutType==0)  // partial sequential (atomic number sequential)
				{
					sprintf(output,"ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
				else            // total sequential (atomic + residue sequential)
				{
					sprintf(output,"ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,cur_pos+1,x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
			}
			//real_output//__110430__//
			if(OutLast==1) //-> at column 77, output 'ATOM_NAME' just as column 13
			{
				output[72]=cur_chain;
				output[77]=output[13];
			}
			fprintf(fp,"%s",output);
		}
		if(OutMode==-1)goto termi; //only output CA
		//sidechain_out
		number=PDB.get_sidechain_totnum();
		for(k=0;k<number;k++)
		{
			//check
			if(OutMode<=0) //only output CB
			{
				if(k!=0)continue;
			}
			if(OutGlys==1) //glycine won't output CB
			{
				if(amino=='G')continue;
			}
			//get_sidechain
			if(PDB.get_sidechain_part_index(k)==0)
			{
				if(OutPrt==1)fprintf(stderr,"SideChain_Bad_Output[chain;%c][pos:%d][atom:%d]\n!!",Real_Chain,i+1,k+1);
				continue;
			}
			atomname=sidechain_atom_name_decode(k,amino);
			PDB.get_sidechain_atom(k,xyz, numb, rfactor, temperature);
			x=xyz.X;
			y=xyz.Y;
			z=xyz.Z;
			//output
			if(OutType==1) // PDB numbering
			{
				sprintf(output,"ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
					numb,atomname,dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
			}
			else           // sequential numbering
			{
				if(OutType==0)  // partial sequential (atomic number sequential)
				{
					sprintf(output,"ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
				else            // total sequential (atomic + residue sequential)
				{
					sprintf(output,"ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						outtype_atom,atomname,dummyaa,Real_Chain,cur_pos+1,x,y,z,rfactor,temperature,buf.c_str());
					outtype_atom++;
				}
			}
			//real_output//__110430__//
			if(OutLast==1) //-> at column 77, output 'ATOM_NAME' just as column 13
			{
				output[72]=cur_chain;
				output[77]=output[13];
			}
			fprintf(fp,"%s",output);
		}
termi:
		cur_pos++;
		//terminal_process
		if(Chain_ID=='!'||Chain_ID==-1)
		{
			if(i<len-1)
			{
				PDB_Residue pdb;
				string pdb_ind;
				pdb=mol[i+1];
				pdb.get_PDB_residue_number(pdb_ind);
				nxt_chain=pdb_ind[0];
				if(cur_chain!=nxt_chain)
				{
					fprintf(fp,"%s\n",TER.c_str());
					cur_pos=0;
				}
			}
		}
	}
	//terminal
	return 1;
}


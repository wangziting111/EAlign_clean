#include "StructureInput.h"


//--------------------- Start -----------------//
StructureInput::StructureInput(void)
{
	//path//
	MOL_ROOT="";
	XYZ_FROM_PDB="xyz_from_pdb";
	DSSP_FROM_PDB="dssp_from_pdb";
	//neo_macro//
	OutPrt=0;       // whether output printf    // (default:no)
	TERorNOT=1;     // whether mark STC_TER     // (default:mark)
	BADorNOT=1;     // whether kill bad residue // (default:kill)
	//name_range
	name_range=100; // the range of input structure's name
	//memory limit
	PRE_LOAD=0;
	WARNING_out=1;
	MEMORY_LIMIT=PROT_MAX_NUM;
}
StructureInput::~StructureInput(void)
{
}

//=============================================== Input_Final_Usage ========================================//(including range)
//string_example: "0 1paz.pdb A" ->
//this means we UPLOAD file "1paz.pdb", with range "A"

//range_example: "A:21-179,B.45-@" ->
//this means the range has two parts:
//part one is Chain A from "21" to "179" in PDB numbering
//part two is Chain B from "45" to "END" in SEQUENTIAL numbering
//check the input string contains valid character or not
//invalid 0-9,A-Z,a-z,": . - + , ! _ @"
int StructureInput::check_range_char(char in)
{
	switch(in)	
	{
		case ':':return 1;
		case '.':return 1;
		case '-':return 1;
		case '+':return 1;
		case ',':return 1;
		case '!':return 1;
		case '_':return 1;
		case '@':return 1;
		default:break;
	}	
	if(in>='A'&&in<='Z')return 1;
	if(in>='a'&&in<='z')return 1;
	if(in>='0'&&in<='9')return 1;
	return 0;
}
int StructureInput::check_range_string(string &in)
{
	int i,len;
	int ret_val;
	char test;
	len=(int)in.length();
	for(i=0;i<len;i++)
	{
		test=in[i];
		ret_val=check_range_char(test);
		if(ret_val==0)return 0;
	}
	return 1;
}
int StructureInput::parse_range(string &in,vector <string> &out_map)
{
	int i,j;
	int len;
	int totnum;
	int head,tail;
	int found;
	int tag;
	int ret_val;
	char wwtemp[14];

	//check_range
	ret_val=check_range_string(in);
	if(ret_val==0)return 0;
	

	//get_tot_num
	head=0;
	totnum=-1;   //current position
	found=0;     //range found
	tag=-1;      //default:bad_state
	len=(int)in.length();
	out_map.clear();
	for(i=0;i<len;i++)
	{
		if(in[i]==':') //PDB_Numbering
		{
			//push_back
			wwtemp[0]=in[i-1];
			for(j=1;j<=10;j++)wwtemp[j]=' ';
			wwtemp[j]='1';
			wwtemp[j+1]='0';
			wwtemp[j+2]='\0';
			out_map.push_back((string)wwtemp);
			totnum++;
			//found
			head=i+1;
			found=1;
		}
		else if(in[i]=='.') //SEQ_Numbering
		{
			//push_back
			wwtemp[0]=in[i-1];
			for(j=1;j<=10;j++)wwtemp[j]=' ';
			wwtemp[j]='0';
			wwtemp[j+1]='0';
			wwtemp[j+2]='\0';
			out_map.push_back((string)wwtemp);
			totnum++;
			//found
			head=i+1;
			found=1;
		}
		else if(in[i]=='-')
		{
			if(found==0)return -3;
			if(i==head)continue;
			if(tag!=-1)continue;
			tail=i-head;
			if(tail<=0)return -3;
			if(in[i-1]>='A' && in[i-1]<='Z')
			{
				for(j=0;j<tail;j++)out_map[totnum][5-j]=in[i-1-j];
			}
			else if(in[i-1]=='@')
			{
				out_map[totnum][5]='@';
			}
			else
			{
				for(j=0;j<tail;j++)out_map[totnum][4-j]=in[i-1-j];
			}
			head=i+1;
			found=2;
			if(tag!=-1)return -3;
			tag=0; // head_to_tail mode
		}
		else if(in[i]=='+')
		{
			if(found==0)return -3;
			tail=i-head;
			if(tail<=0)return -3;
			if(in[i-1]>='A' && in[i-1]<='Z')
			{
				for(j=0;j<tail;j++)out_map[totnum][5-j]=in[i-1-j];
			}
			else if(in[i-1]=='@')
			{
				out_map[totnum][5]='@';
			}
			else
			{
				for(j=0;j<tail;j++)out_map[totnum][4-j]=in[i-1-j];
			}
			head=i+1;
			found=3;
			if(tag!=-1)return -3;
			tag=1; // head_plus_len mode
		}
		else if(in[i]==',')
		{
			if(found==2)
			{
				tail=i-head;
				if(tail<=0)return -3;
				if(in[i-1]>='A' && in[i-1]<='Z')
				{
					for(j=0;j<tail;j++)out_map[totnum][10-j]=in[i-1-j];
				}
				else if(in[i-1]=='@')
				{
					out_map[totnum][10]='@';
				}
				else
				{
					for(j=0;j<tail;j++)out_map[totnum][9-j]=in[i-1-j];
				}
				out_map[totnum][12]='0'; // head_to_tail mode
			}
			else if(found==3)
			{
				tail=i-head;
				if(tail<=0)return -3;
				if(in[i-1]>='A' && in[i-1]<='Z')return -3;
				if(in[i-1]=='@')return -3;
				for(j=0;j<tail;j++)out_map[totnum][9-j]=in[i-1-j];
				out_map[totnum][12]='1'; // head_plus_len mode
			}
			else if(found==0)
			{
				//push_back
				wwtemp[0]=in[i-1];
				for(j=1;j<=10;j++)wwtemp[j]=' ';
				wwtemp[j]='1';     //PDB_Numbering
				wwtemp[j+1]='0';   //head_to_tail mode
				wwtemp[j+2]='\0';
				out_map.push_back((string)wwtemp);
				totnum++;
			}
			found=0;
			tag=-1;      //default:bad_state
		}
	}

	//end process
	if(found==2)
	{
		tail=i-head;
		if(tail<=0)return -3;
		if(in[i-1]>='A' && in[i-1]<='Z')
		{
			for(j=0;j<tail;j++)out_map[totnum][10-j]=in[i-1-j];
		}
		else if(in[i-1]=='@')
		{
			out_map[totnum][10]='@';
		}
		else
		{
			for(j=0;j<tail;j++)out_map[totnum][9-j]=in[i-1-j];
		}
		out_map[totnum][12]='0'; // head_to_tail mode
	}
	else if(found==3)
	{
		tail=i-head;
		if(tail<=0)return -3;
		if(in[i-1]>='A' && in[i-1]<='Z')return -3;
		if(in[i-1]=='@')return -3;
		for(j=0;j<tail;j++)out_map[totnum][9-j]=in[i-1-j];
		out_map[totnum][12]='1'; // head_plus_len mode
	}
	else if(found==0)
	{
		//push_back
		wwtemp[0]=in[i-1];
		for(j=1;j<=10;j++)wwtemp[j]=' ';
		wwtemp[j]='1';     //PDB_Numbering
		wwtemp[j+1]='0';   //head_to_tail mode
		wwtemp[j+2]='\0';
		out_map.push_back((string)wwtemp);
		totnum++;
	}
	return 1;
}


//================================ XYZ_Related =================================================================//
//----------- File_Parse ----------//
//------------ Ori_Process ------------//
int StructureInput::Input_XYZ_FULL(string &fn,char chain,int PS,int mode,
	int st,char stt,int ed,char edd,int &totnu,
	XYZ *mol,char *AMI,char *CLE,char *ind) 
{
	ifstream fin;
	string buf,temp;
	char pdbid[9];
	char totn[5];
	int ws_i,i;
	int count;
	int num;
	int wtotnum;
	double x,y,z;
	char ami_c,stc_c,ins_c;
	int index_i;
	int first,last;
	int len,cur;
	int ws_exclaim;


	//--check--//
	ins_c=' ';
	if(totnu<0)return RANGEOVER_ERROR;
	count=totnu;

	fin.open(fn.c_str(), ios::in);
	if (fin.fail()!=0)
	{
		if(OutPrt==1)fprintf(stderr,"XYZ File %s Not Found!!!\n",fn.c_str());
		return FILE_LOAD_ERROR;  // no such file;
	}

	if(chain=='!')ws_exclaim=1;
	else ws_exclaim=0;

	num=0;
	for(;;)
	{
ws_begi:
		if(!getline(fin,buf,'\n'))goto end;  // premature break;
		getline_end(buf,0x0D);
		if(buf[0]=='>')break;
	}
	len=(int)buf.length();

	cur=1;
	for(ws_i=0;;ws_i++)
	{
		if(cur>=len || buf[cur]==' ')break;
		pdbid[ws_i]=buf[cur];
		cur++;
	}
	pdbid[ws_i]='\0';

	for(;;)
	{
		if(buf[cur]!=' ')break;
		cur++;
	}

	for(ws_i=0;;ws_i++)
	{
		if(cur>=len || buf[cur]==' ')break;
		totn[ws_i]=buf[cur];
		cur++;
	}
	totn[ws_i]='\0';

	//check_chain
	if(ws_exclaim==0)
	{
		if(strlen(pdbid)==5) // PDB_XYZ
		{
			if(chain!='_'&&pdbid[4]!=chain)goto ws_begi;
		}
	}

	if(str2int(totn,wtotnum)!=1)
	{
		if(OutPrt==1)fprintf(stderr,"TOTNUM ERROR => XYZ[%s]:NUM[%4s]!\n",pdbid,totn);
		return STR_TRANS_ERROR; //transfer error;
	}

	//assign_head_tail
	if(stt=='@')first=1;
	else first=0;
	if(edd=='@')last=1;
	else last=0;

	if(mode==1) //head_plus_len
	{
		if(ed<=0)last=1;
		else last=0;
	}
	
	for(i=0;i<wtotnum;i++)
	{
		if(!getline(fin,buf,'\n'))return PREMATURE_ERROR; //premature break;
		getline_end(buf,0x0D);
		if(buf.length()!=35)return DATA_LENG_ERROR;      // length error;

		if(buf[0]=='!')    //neglect miss residue
		{
			i--;
			continue;
		}
		if(buf[0]=='x')    //neglect bad residue
		{
			if(BADorNOT==1)
			{
				i--;
				wtotnum--;
				continue;
			}
		}
		if(ws_exclaim==1)goto next;


		//--get_index--//
		temp=buf.substr(2,4);
		if(str2int(temp.c_str(),index_i)!=1)
		{
			if(OutPrt==1)fprintf(stderr,"ERROR => BAD INDEX AT XYZ[%s]:PDB[%s]!\n",pdbid,buf.c_str());
			return STR_TRANS_ERROR; //transfer error;
		}
		ins_c=buf[6];

		//-----first_check-----//
		if(first==0)
		{
			if(PS==1)  // PDB_numbering
			{
				if(index_i==st && ins_c==stt)first=1;
				else continue;
			}
			else       // SEQ_numbering
			{
				if(i==st-1)first=1;
				else continue;
			}
		}

next:
		//=========== Get INDEX map =============//
		if(ind!=NULL)for(ws_i=0;ws_i<6;ws_i++)ind[6*count+ws_i]=buf[1+ws_i];

		//=========== Get AMI and CLE ===========//
		if(AMI!=NULL)
		{
			ami_c=buf[8];
			if(ami_c<'A' || ami_c>'Z')return FILE_FORM_ERROR; //error
			AMI[count]=ami_c;
		}
		if(CLE!=NULL)
		{
			stc_c=buf[9];
			if(stc_c<'A' || stc_c>'R')return FILE_FORM_ERROR; //error
			CLE[count]=stc_c;
		}

		//=========== Get XYZ coordinate ========//
		if(mol!=NULL)
		{
			temp=buf.substr(11,8);  // get PDB_FILE x coordinate (8)
			if(str2dou(temp.c_str(),x)!=1)
			{
				if(OutPrt==1)fprintf(stderr,"ERROR => BAD X_POS AT XYZ[%s]:PDB[%s]!\n",pdbid,buf.c_str());
				return STR_TRANS_ERROR; //transfer error;
			}
			temp=buf.substr(19,8);  // get PDB_FILE y coordinate (8)
			if(str2dou(temp.c_str(),y)!=1)
			{
				if(OutPrt==1)fprintf(stderr,"ERROR => BAD Y_POS AT XYZ[%s]:PDB[%s]!\n",pdbid,buf.c_str());
				return STR_TRANS_ERROR; //transfer error;
			}
			temp=buf.substr(27,8);  // get PDB_FILE z coordinate (8)
			if(str2dou(temp.c_str(),z)!=1)
			{
				if(OutPrt==1)fprintf(stderr,"ERROR => BAD Z_POS AT XYZ[%s]:PDB[%s]!\n",pdbid,buf.c_str());
				return STR_TRANS_ERROR; //transfer error;
			}

			//========== assign coordinate ========//
			mol[count].X=x;
			mol[count].Y=y;
			mol[count].Z=z;
		}


		//-----step_gain-----//
		count++;
		num++;
		if(ws_exclaim==1)continue;

		//-----last_check-----//
		if(last==0)
		{
			if(mode==0) // head_to_tail
			{
				if(PS==1)  // PDB_numbering
				{
					if(index_i==ed && ins_c==edd)break;
				}
				else       // SEQ_numbering
				{
					if(i==ed-1)break;
				}
			}
			else       // head_plus_len
			{
				if(num==ed)break;
			}
		}
	}
	if(ws_exclaim==1)goto ws_begi;

end:
	if(num==0)return RANGEOVER_ERROR;

	if(TERorNOT==1 && CLE!=NULL)
	{
		CLE[totnu]='R';
		CLE[totnu+1]='R';
		CLE[count-1]='R';
	}
	if(CLE!=NULL)CLE[count]='\0';
	if(AMI!=NULL)AMI[count]='\0';	

	totnu=count;
	return 1;
}
//------------ Uni_Process ------------//
int StructureInput::XYZ_Input(string &file,string &range,int form,int &moln,XYZ *mol,char *AMI,char *CLE,char *ind,PDB_Residue *pdb)
{
	int ret_val;
	int i,k;
	int totnum;
	int tag;
	char chain;
	char chain_bak;
	char head[6],tail[6];
	int check;
	int head_tag,tail_tag;
	char head_c,tail_c;
	int PDBorSEQ;

	//range_parse
	ret_val=parse_range(range,out_map);
	if(ret_val!=1)return USER_TYPE_ERROR;  //means range process failed

	//real_process
	totnum=(int)out_map.size();
	moln=0;
	for(i=0;i<totnum;i++)
	{
		//get_init
		chain=out_map[i][0];
		PDBorSEQ=out_map[i][11]-'0';
		tag=out_map[i][12]-'0';
		for(k=0;k<5;k++)
		{
			head[k]=out_map[i][1+k];
			tail[k]=out_map[i][6+k];
		}
		head[k]='\0';
		tail[k]='\0';
		//get_head
		check=0;  //default all ' '
		for(k=0;k<5;k++)if(head[k]==' ')check++;
		if(check==5 || head[4]=='@')
		{
			head_c='@';
			head_tag=-1;			
		}
		else
		{
			head_c=head[4];
			head[4]='\0';
			if(str2int(head,head_tag)!=1)return STR_TRANS_ERROR;
		}
		//get_tail
		check=0;  //default all ' '
		for(k=0;k<5;k++)if(tail[k]==' ')check++;
		if(check==5 || tail[4]=='@')
		{
			tail_c='@';
			tail_tag=-1;
		}
		else
		{
			tail_c=tail[4];
			tail[4]='\0';
			if(str2int(tail,tail_tag)!=1)return STR_TRANS_ERROR;
		}
		//chain_process
		if(chain>='A' && chain<='Z')chain_bak=chain-'A'+'a';
		else if(chain>='a' && chain<='z')chain_bak=chain-'a'+'A';
		else chain_bak=chain;
		//head & tail checking //__110230__//
		if(head_tag!=-1&&tail_tag!=-1)
		{
			if(head_tag>tail_tag)return UPSIDDOWN_ERROR;
		}

		//========== Input_File =============//
		if(form==0)      //PDB mode
		{
			ret_val=Upload_Process(file,chain,chain_bak,PDBorSEQ,tag,
				head_tag,head_c,tail_tag,tail_c,
				moln,mol,AMI,CLE,ind,pdb);
//			if(ret_val!=1)return ret_val; //means PDB mode failed
			if(ret_val!=1) //swith to '@' start
			{
				//----- we need to printf WARNING here -------//start
				if(WARNING_out==1)
				{
					fprintf(stderr,"WARNING !! file %s , the range %s at chain %c that the user designates is not correct. \n EAlign will try to get the intersection for the given range and the current chain. \n",file.c_str(),range.c_str(),chain);
				}
				//----- we need to printf WARNING here -------//end

				int moln_;
				int min_num=9999999;
				int min_cat=-1;
				int PRE_LOAD_=PRE_LOAD;
				PRE_LOAD=1;
				//head
				moln_=0;
				ret_val=Upload_Process(file,chain,chain_bak,PDBorSEQ,tag,
					-1,'@',tail_tag,tail_c,moln_,0,0,0,0,0);
				if(ret_val==1)
				{
					if(moln_<min_num)
					{
						min_num=moln_;
						min_cat=0;  //-> head
					}
				}
				//tail
				moln_=0;
				ret_val=Upload_Process(file,chain,chain_bak,PDBorSEQ,tag,
					head_tag,head_c,-1,'@',moln_,0,0,0,0,0);
				if(ret_val==1)
				{
					if(moln_<min_num)
					{
						min_num=moln_;
						min_cat=1;  //-> tail
					}
				}
				//final judge
				if(min_cat!=-1)
				{
					if(min_cat==0)
					{
						Upload_Process(file,chain,chain_bak,PDBorSEQ,tag,
							-1,'@',tail_tag,tail_c,moln,mol,AMI,CLE,ind,pdb);
					}
					else
					{
						Upload_Process(file,chain,chain_bak,PDBorSEQ,tag,
							head_tag,head_c,-1,'@',moln,mol,AMI,CLE,ind,pdb);
					}
				}
				else
				{
					Upload_Process(file,chain,chain_bak,PDBorSEQ,tag,
						-1,'@',-1,'@',moln,mol,AMI,CLE,ind,pdb);
				}
				//final
				PRE_LOAD=PRE_LOAD_;
			}
		}
		else if(form==1) //XYZ mode
		{
			ret_val=Input_XYZ_FULL(file,chain,PDBorSEQ,tag,
				head_tag,head_c,tail_tag,tail_c,
				moln,mol,AMI,CLE,ind);
			if(ret_val!=1 && chain_bak!=chain)ret_val=Input_XYZ_FULL(file,chain_bak,PDBorSEQ,tag,
				head_tag,head_c,tail_tag,tail_c,
				moln,mol,AMI,CLE,ind);
			if(ret_val!=1)return ret_val; //means XYZ mode failed
		}
		else             //DATABASE mode
		{
			ret_val=Database_Process_XYZ_II(chain,chain_bak,PDBorSEQ,tag,
				head_tag,head_c,tail_tag,tail_c,
				moln,mol,AMI,CLE,ind);
			if(ret_val!=1)return ret_val; //means DATABASE mode failed
		}
	}
	return 1;
}
//------------ Single_Process ------------//
//------------ Batch_Process  ------------//

//======================================== DSSP_Related ==================================================//
//----------- File_Parse ----------//
//------------ Ori_Process ------------//
//------------ Uni_Process ------------//
//------------ Single_Process ------------//
//------------ Batch_Process  ------------//

//============================== Virtual Function ========================//
int StructureInput::Upload_Process(string &file,char chain,char chain_bak,int PDBorSEQ,int tag,
							int head_tag,int head_c,int tail_tag,int tail_c,
							int &moln,XYZ *mol,char *AMI,char *CLE,char *ind,PDB_Residue *pdb)
{
	return UNSUPPORT_ERROR;
}

int StructureInput::Database_Process_XYZ_I(string &name,int wstag)
{
	return UNSUPPORT_ERROR;
}

int StructureInput::Database_Process_XYZ_II(char chain,char chain_bak,int PDBorSEQ,int tag,
									 int head_tag,int head_c,int tail_tag,int tail_c,int &moln,
									 XYZ *mol,char *AMI,char *CLE,char *ind)
{
	return UNSUPPORT_ERROR;
}
int StructureInput::Database_Process_DSSP_I(string &name,int wstag)
{
	return UNSUPPORT_ERROR;
}
int StructureInput::Database_Process_DSSP_II(char chain,char chain_bak,int PDBorSEQ,int tag,
									   int head_tag,int head_c,int tail_tag,int tail_c,int &moln,
									   XYZ *mol,char *AMI,char *CLE,char *ind,
									   int *acc,int *beta,double *tco,double *kappa,
									   double *alpha,double *phi,double *psi)
{
	return UNSUPPORT_ERROR;
}

#include "Foldseek3Di.h"
#include "Foldseek3DiModel.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace
{
	const double PI=3.14159265358979323846;
	const double CA_CB_DISTANCE=1.5336;
	const int FEATURE_COUNT=10;
	const int EMBEDDING_SIZE=2;

	struct Vec3
	{
		double x,y,z;
		Vec3(void):x(0.0),y(0.0),z(0.0) {}
		Vec3(double x_,double y_,double z_):x(x_),y(y_),z(z_) {}
	};

	Vec3 ToVec3(const XYZ &value)
	{
		return Vec3(value.X,value.Y,value.Z);
	}

	Vec3 Add(const Vec3 &a,const Vec3 &b)
	{
		return Vec3(a.x+b.x,a.y+b.y,a.z+b.z);
	}

	Vec3 Subtract(const Vec3 &a,const Vec3 &b)
	{
		return Vec3(a.x-b.x,a.y-b.y,a.z-b.z);
	}

	Vec3 Scale(const Vec3 &a,double factor)
	{
		return Vec3(a.x*factor,a.y*factor,a.z*factor);
	}

	double Dot(const Vec3 &a,const Vec3 &b)
	{
		return a.x*b.x+a.y*b.y+a.z*b.z;
	}

	Vec3 Cross(const Vec3 &a,const Vec3 &b)
	{
		return Vec3(a.y*b.z-a.z*b.y,
			a.z*b.x-a.x*b.z,
			a.x*b.y-a.y*b.x);
	}

	bool Normalize(const Vec3 &input,Vec3 &output)
	{
		const double squared=Dot(input,input);
		if(!std::isfinite(squared) || squared<=1.0e-16)return false;
		output=Scale(input,1.0/std::sqrt(squared));
		return true;
	}

	double Distance(const Vec3 &a,const Vec3 &b)
	{
		const Vec3 delta=Subtract(a,b);
		return std::sqrt(Dot(delta,delta));
	}

	bool IsFinite(const Vec3 &value)
	{
		return std::isfinite(value.x) && std::isfinite(value.y) &&
			std::isfinite(value.z);
	}

	bool ApproximateCBeta(const Vec3 &ca,const Vec3 &n,const Vec3 &c,
		Vec3 &cb)
	{
		Vec3 ca_to_c,ca_to_n;
		if(!Normalize(Subtract(c,ca),ca_to_c) ||
			!Normalize(Subtract(n,ca),ca_to_n))return false;

		const Vec3 b1=Add(ca_to_n,Scale(ca_to_c,1.0/3.0));
		Vec3 u1,u2;
		if(!Normalize(b1,u1) || !Normalize(Cross(ca_to_c,b1),u2))return false;
		const Vec3 tetrahedral=Subtract(Scale(u1,-0.5),
			Scale(u2,std::sqrt(3.0)/2.0));
		const Vec3 direction=Add(Scale(ca_to_c,-1.0/3.0),
			Scale(tetrahedral,std::sqrt(8.0)/3.0));
		cb=Add(ca,Scale(direction,CA_CB_DISTANCE));
		return IsFinite(cb);
	}

	bool VirtualCenter(const Vec3 &ca,const Vec3 &cb,const Vec3 &n,
		Vec3 &center)
	{
		Vec3 axis;
		const Vec3 ca_to_cb=Subtract(cb,ca);
		if(!Normalize(Cross(ca_to_cb,Subtract(n,ca)),axis))return false;
		const double alpha=270.0/180.0*PI;
		const Vec3 rotated=Add(
			Add(Scale(ca_to_cb,std::cos(alpha)),
				Scale(Cross(axis,ca_to_cb),std::sin(alpha))),
			Scale(axis,Dot(axis,ca_to_cb)*(1.0-std::cos(alpha))));
		center=Add(ca,Scale(rotated,2.0));
		return IsFinite(center);
	}

	struct DenseLayer
	{
		int rows,columns,activation;
		std::vector<float> weights;
		std::vector<float> biases;
	};

	class TinyEncoder
	{
	public:
		TinyEncoder(void)
		{
			size_t offset=0;
			const unsigned int layer_count=ReadUnsigned(offset);
			if(layer_count!=3)throw std::runtime_error("unexpected Foldseek 3Di encoder");
			layers.resize(layer_count);
			for(unsigned int layer=0;layer<layer_count;layer++)
			{
				if(ReadUnsigned(offset)!=1)
					throw std::runtime_error("3Di encoder contains a non-dense layer");
				DenseLayer &dense=layers[layer];
				dense.rows=(int)ReadUnsigned(offset);
				dense.columns=(int)ReadUnsigned(offset);
				const int bias_count=(int)ReadUnsigned(offset);
				dense.weights.resize(dense.rows*dense.columns);
				dense.biases.resize(bias_count);
				for(size_t i=0;i<dense.weights.size();i++)
					dense.weights[i]=ReadFloat(offset);
				for(size_t i=0;i<dense.biases.size();i++)
					dense.biases[i]=ReadFloat(offset);
				dense.activation=(int)ReadUnsigned(offset);
				if(bias_count!=dense.columns ||
					(dense.activation!=1 && dense.activation!=2))
					throw std::runtime_error("invalid Foldseek 3Di dense layer");
			}
			if(offset!=foldseek_3di_encoder_model_len)
				throw std::runtime_error("trailing data in Foldseek 3Di encoder");
		}

		void Apply(const double features[FEATURE_COUNT],
			double embedding[EMBEDDING_SIZE]) const
		{
			std::vector<float> input(FEATURE_COUNT);
			for(int i=0;i<FEATURE_COUNT;i++)input[i]=(float)features[i];
			for(size_t layer=0;layer<layers.size();layer++)
			{
				const DenseLayer &dense=layers[layer];
				std::vector<float> output(dense.columns,0.0f);
				for(int column=0;column<dense.columns;column++)
				{
					float value=dense.biases[column];
					for(int row=0;row<dense.rows;row++)
						value+=input[row]*dense.weights[row*dense.columns+column];
					if(dense.activation==2 && value<0.0f)value=0.0f;
					output[column]=value;
				}
				input.swap(output);
			}
			embedding[0]=input[0];
			embedding[1]=input[1];
		}

	private:
		std::vector<DenseLayer> layers;

		unsigned int ReadUnsigned(size_t &offset) const
		{
			if(offset+4>foldseek_3di_encoder_model_len)
				throw std::runtime_error("truncated Foldseek 3Di encoder");
			unsigned int value=0;
			for(int byte=0;byte<4;byte++)
				value|=((unsigned int)foldseek_3di_encoder_model[offset+byte])<<(8*byte);
			offset+=4;
			return value;
		}

		float ReadFloat(size_t &offset) const
		{
			const unsigned int bits=ReadUnsigned(offset);
			float value;
			std::memcpy(&value,&bits,sizeof(value));
			return value;
		}
	};

	const double CENTROIDS[Foldseek3Di::STATE_COUNT][EMBEDDING_SIZE]={
		{-1.0729,-0.3600},{-0.1356,-1.8914},{0.4948,-0.4205},
		{-0.9874,0.8128},{-1.6621,-0.4259},{2.1394,0.0486},
		{1.5558,-0.1503},{2.9179,1.1437},{-2.8814,0.9956},
		{-1.1400,-2.0068},{3.2025,1.7356},{1.7769,-1.3037},
		{0.6901,-1.2554},{-1.1061,-1.3397},{2.1495,-0.8030},
		{2.3060,-1.4988},{2.5522,0.6046},{0.7786,-2.1660},
		{-2.3030,0.3813},{1.0290,0.8772}
	};

	const int SUBSTITUTION_MATRIX[Foldseek3Di::STATE_COUNT][Foldseek3Di::STATE_COUNT]={
		{ 6,-3, 1, 2, 3,-2,-2,-7,-3,-3,-10,-5,-1, 1,-4,-7,-5,-6, 0,-2},
		{-3, 6,-2,-8,-5,-4,-4,-12,-13, 1,-14, 0, 0, 1,-1, 0,-8, 1,-7,-9},
		{ 1,-2, 4,-3, 0, 1, 1,-3,-5,-4,-5,-2, 1,-1,-1,-4,-2,-3,-2,-2},
		{ 2,-8,-3, 9,-2,-7,-4,-12,-10,-7,-17,-8,-6,-3,-8,-10,-10,-13,-6,-3},
		{ 3,-5, 0,-2, 7,-3,-3,-5, 1,-3,-9,-5,-2, 2,-5,-8,-3,-7, 4,-4},
		{-2,-4, 1,-7,-3, 6, 3, 0,-7,-7,-1,-2,-2,-4, 3,-3, 4,-6,-4,-2},
		{-2,-4, 1,-4,-3, 3, 6,-4,-7,-6,-6, 0,-1,-3, 1,-3,-1,-5,-5, 3},
		{-7,-12,-3,-12,-5, 0,-4, 8,-5,-11, 7,-7,-6,-6,-3,-9, 6,-12,-5,-8},
		{-3,-13,-5,-10, 1,-7,-7,-5, 9,-11,-8,-12,-6,-5,-9,-14,-5,-15, 5,-8},
		{-3, 1,-4,-7,-3,-7,-6,-11,-11, 6,-16,-3,-2, 2,-4,-4,-9, 0,-8,-9},
		{-10,-14,-5,-17,-9,-1,-6, 7,-8,-16,10,-9,-9,-10,-5,-10, 3,-16,-6,-9},
		{-5, 0,-2,-8,-5,-2, 0,-7,-12,-3,-9, 7, 0,-2, 2, 3,-4, 0,-8,-5},
		{-1, 0, 1,-6,-2,-2,-1,-6,-6,-2,-9, 0, 4, 0, 0,-2,-4, 0,-4,-5},
		{ 1, 1,-1,-3, 2,-4,-3,-6,-5, 2,-10,-2, 0, 5,-2,-4,-5,-1,-2,-5},
		{-4,-1,-1,-8,-5, 3, 1,-3,-9,-4,-5, 2, 0,-2, 6, 2, 0,-1,-6,-3},
		{-7, 0,-4,-10,-8,-3,-3,-9,-14,-4,-10, 3,-2,-4, 2, 6,-6, 0,-11,-9},
		{-5,-8,-2,-10,-3, 4,-1, 6,-5,-9, 3,-4,-4,-5, 0,-6, 8,-9,-5,-5},
		{-6, 1,-3,-13,-7,-6,-5,-12,-15, 0,-16, 0, 0,-1,-1, 0,-9, 3,-10,-11},
		{ 0,-7,-2,-6, 4,-4,-5,-5, 5,-8,-6,-8,-4,-2,-6,-11,-5,-10, 8,-6},
		{-2,-9,-2,-3,-4,-2, 3,-8,-8,-9,-9,-5,-5,-5,-3,-9,-5,-11,-6, 9}
	};

	bool CalculateFeatures(const std::vector<Vec3> &ca,int i,int partner,
		double features[FEATURE_COUNT])
	{
		Vec3 u1,u2,u3,u4,u5;
		if(!Normalize(Subtract(ca[i],ca[i-1]),u1) ||
			!Normalize(Subtract(ca[i+1],ca[i]),u2) ||
			!Normalize(Subtract(ca[partner],ca[partner-1]),u3) ||
			!Normalize(Subtract(ca[partner+1],ca[partner]),u4) ||
			!Normalize(Subtract(ca[partner],ca[i]),u5))return false;

		features[0]=Dot(u1,u2);
		features[1]=Dot(u3,u4);
		features[2]=Dot(u1,u5);
		features[3]=Dot(u3,u5);
		features[4]=Dot(u1,u4);
		features[5]=Dot(u2,u3);
		features[6]=Dot(u1,u3);
		features[7]=Distance(ca[i],ca[partner]);
		const int separation=partner-i;
		features[8]=(double)std::max(-4,std::min(4,separation));
		features[9]=std::copysign(std::log(std::abs(separation)+1.0),
			(double)separation);
		return true;
	}
}

void Foldseek3Di::Encode(const PDB_Residue *residues,const XYZ *ca_input,
	int length,int *states)
{
	if(states==0)return;
	for(int i=0;i<length;i++)states[i]=INVALID_STATE;
	if(residues==0 || ca_input==0 || length<3)return;

	std::vector<Vec3> ca(length),virtual_centers(length);
	std::vector<bool> valid(length,false);
	for(int i=0;i<length;i++)
	{
		ca[i]=ToVec3(ca_input[i]);
		XYZ n_xyz,c_xyz,cb_xyz;
		if(!IsFinite(ca[i]) ||
			residues[i].get_backbone_atom("N  ",n_xyz)!=0 ||
			residues[i].get_backbone_atom("C  ",c_xyz)!=0)continue;
		const Vec3 n=ToVec3(n_xyz);
		const Vec3 c=ToVec3(c_xyz);
		Vec3 cb;
		bool has_cb=residues[i].get_sidechain_atom("CB ",cb_xyz)==0;
		if(has_cb)
		{
			cb=ToVec3(cb_xyz);
			const Vec3 delta=Subtract(cb,ca[i]);
			const double squared=Dot(delta,delta);
			// EAlign stores CA as a placeholder CB for glycine/missing CB.
			has_cb=std::isfinite(squared) && squared>=1.0 && squared<=4.0;
		}
		if(!has_cb && !ApproximateCBeta(ca[i],n,c,cb))continue;
		if(!IsFinite(n) || !IsFinite(c) || !IsFinite(cb) ||
			!VirtualCenter(ca[i],cb,n,virtual_centers[i]))continue;
		valid[i]=true;
	}

	std::vector<int> partners(length,-1);
	for(int i=1;i<length-1;i++)
	{
		if(!valid[i])continue;
		double nearest_distance=HUGE_VAL;
		for(int j=1;j<length-1;j++)
		{
			if(i==j || !valid[j])continue;
			const double distance=Distance(virtual_centers[i],virtual_centers[j]);
			if(distance<nearest_distance)
			{
				nearest_distance=distance;
				partners[i]=j;
			}
		}
	}

	static const TinyEncoder encoder;
	for(int i=1;i<length-1;i++)
	{
		const int partner=partners[i];
		if(partner<=0 || partner>=length-1 ||
			!valid[i-1] || !valid[i] || !valid[i+1] ||
			!valid[partner-1] || !valid[partner] || !valid[partner+1])continue;
		double features[FEATURE_COUNT];
		if(!CalculateFeatures(ca,i,partner,features))continue;
		double embedding[EMBEDDING_SIZE];
		encoder.Apply(features,embedding);
		double nearest_centroid=HUGE_VAL;
		for(int state=0;state<STATE_COUNT;state++)
		{
			const double dx=embedding[0]-CENTROIDS[state][0];
			const double dy=embedding[1]-CENTROIDS[state][1];
			const double squared=dx*dx+dy*dy;
			if(squared<nearest_centroid)
			{
				nearest_centroid=squared;
				states[i]=state;
			}
		}
	}
}

int Foldseek3Di::Score(int state1,int state2)
{
	if(state1<0 || state1>=STATE_COUNT || state2<0 || state2>=STATE_COUNT)
		return 0;
	return SUBSTITUTION_MATRIX[state1][state2];
}

double Foldseek3Di::CombinedEvidence(double amino_acid_score,
	double three_di_score)
{
	
	return amino_acid_score+three_di_score;
}

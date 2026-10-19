#include "SpatialAlignmentProposal.h"

#include <stdexcept>

namespace
{
	void InsertNeighbor(const SpatialNeighbor &candidate,int max_neighbors,
		SpatialNeighbor *neighbors,int &count)
	{
		int insert_at=count;
		while(insert_at>0)
		{
			const SpatialNeighbor &previous=neighbors[insert_at-1];
			if(previous.distance_square<candidate.distance_square)break;
			if(previous.distance_square==candidate.distance_square &&
				previous.index<candidate.index)break;
			insert_at--;
		}
		if(insert_at>=max_neighbors)return;
		const int last=count<max_neighbors?count:max_neighbors-1;
		for(int i=last;i>insert_at;i--)neighbors[i]=neighbors[i-1];
		neighbors[insert_at]=candidate;
		if(count<max_neighbors)count++;
	}
}

SpatialAlignmentProposal::SpatialAlignmentProposal(int max_candidates)
:max_candidates_(max_candidates),row_count_(0)
{
	if(max_candidates_<=0)
		throw std::invalid_argument("SpatialAlignmentProposal candidate count must be positive");
}

bool SpatialAlignmentProposal::BuildFromCoordinates(const XYZ *mol1,int moln1,
	const XYZ *mol2,int moln2,const double *rotmat,double radius)
{
	row_count_=0;
	counts_.clear();
	candidates_.clear();
	if(mol1==NULL || mol2==NULL || rotmat==NULL || moln1<=0 || moln2<=0 ||
		radius<=0.0)return false;
	row_count_=moln2;
	counts_.assign(row_count_,0);
	candidates_.resize(row_count_*max_candidates_);
	const double radius_square=radius*radius;
	for(int j=0;j<moln2;j++)
	{
		const double dx=mol2[j].X-rotmat[9];
		const double dy=mol2[j].Y-rotmat[10];
		const double dz=mol2[j].Z-rotmat[11];
		XYZ query;
		query.X=dx*rotmat[0]+dy*rotmat[3]+dz*rotmat[6];
		query.Y=dx*rotmat[1]+dy*rotmat[4]+dz*rotmat[7];
		query.Z=dx*rotmat[2]+dy*rotmat[5]+dz*rotmat[8];
		SpatialNeighbor *row=&candidates_[j*max_candidates_];
		int count=0;
		for(int i=0;i<moln1;i++)
		{
			const double distance_square=query.distance_square(mol1[i]);
			if(distance_square>radius_square)continue;
			const SpatialNeighbor candidate={i,distance_square};
			InsertNeighbor(candidate,max_candidates_,row,count);//Top6 Ali3 构建
		}
		counts_[j]=count;
	}
	return true;
}



double SpatialAlignmentProposal::TraceOrderPreservingAlignment(double d0,int moln1,
	std::vector <int> &ali2) const
{
	if(d0<=0.0 || moln1<=0)
	{
		ali2.assign(row_count_,-1);
		return 0.0;
	}
	const double d0_square=d0*d0;
	std::vector <double> candidate_scores(candidates_.size(),0.0);
	for(int j=0;j<row_count_;j++)
	{
		for(int k=0;k<counts_[j];k++)
		{
			const int node=j*max_candidates_+k;
			candidate_scores[node]=1.0/
				(1.0+candidates_[node].distance_square/d0_square);
		}
	}
	return TraceOrderPreservingAlignment(candidate_scores,moln1,ali2);
}
//基于稀疏候选集的顺序保持对齐（Order-Preserving Alignment）动态规划
//在 Top-6 稀疏候选 residue pairs 上求一个“最大权重、保持序列顺序”的 alignment
double SpatialAlignmentProposal::TraceOrderPreservingAlignment(
	const std::vector <double> &candidate_scores,int moln1,
	std::vector <int> &ali2) const
{
	ali2.assign(row_count_,-1);
	if(moln1<=0)return 0.0;
	if(candidate_scores.size()!=candidates_.size())
		throw std::invalid_argument("Sparse candidate score count does not match proposal");
	const int node_count=row_count_*max_candidates_;
	std::vector <double> node_scores(node_count,0.0);//以当前 candidate (i,j) 结尾的最佳 alignment score
	std::vector <int> parents(node_count,-1);//当前节点最佳路径中的前驱节点，用于 traceback
	std::vector <double> best_scores_by_index(moln1,0.0);//已经处理过的 mol2 行中，以 mol1[i] 结尾的最高路径分数
	std::vector <int> best_nodes_by_index(moln1,-1);//产生上述最高分的 candidate node
	std::vector <int> row_nodes(max_candidates_,-1);

	for(int j=0;j<row_count_;j++)
	{
		const int count=counts_[j];// 当前行实际拥有的候选数量（≤ max_candidates_）
		for(int k=0;k<count;k++)
		{
			const int node=j*max_candidates_+k;
			const SpatialNeighbor &candidate=candidates_[node];
			double previous_score=0.0;
			int previous_node=-1;
			for(int i=0;i<candidate.index;i++)//当前候选匹配 : mol2[j] ↔ mol1[candidate.index]
			{
				if(best_scores_by_index[i]>previous_score ||
					(best_scores_by_index[i]==previous_score &&
					best_nodes_by_index[i]>=0 &&
					(previous_node<0 || best_nodes_by_index[i]<previous_node)))
				{
					previous_score=best_scores_by_index[i];
					previous_node=best_nodes_by_index[i];
				}
			}
			node_scores[node]=previous_score+candidate_scores[node];//\[DP(i,j)=w(i,j)+\max_{\substack{i'<i\\j'<j}}DP(i',j')\]
			parents[node]=previous_node;
			row_nodes[k]=node;
		}

		// Delay updates until the complete mol2 row is scored.
		for(int k=0;k<count;k++)
		{
			const int node=row_nodes[k];
			const int index=candidates_[node].index;
			if(node_scores[node]>best_scores_by_index[index] ||
				(node_scores[node]==best_scores_by_index[index] &&
				(best_nodes_by_index[index]<0 ||
				node<best_nodes_by_index[index])))
			{
				best_scores_by_index[index]=node_scores[node];
				best_nodes_by_index[index]=node;
			}
		}
	}

	double best_score=0.0;
	int node=-1;
	for(int i=0;i<moln1;i++)
	{
		if(best_scores_by_index[i]>best_score ||
			(best_scores_by_index[i]==best_score &&
			best_nodes_by_index[i]>=0 &&
			(node<0 || best_nodes_by_index[i]<node)))
		{
			best_score=best_scores_by_index[i];
			node=best_nodes_by_index[i];//这里的 node 就是最优 alignment path 的最后一个 residue pair
		}
	}
	while(node>=0)//traceback
	{
		const int mol2_index=node/max_candidates_;//j = node / K
		ali2[mol2_index]=candidates_[node].index;// i = candidates[node].index,ali2[j] = i
		node=parents[node];
	}
	return best_score;
}

int SpatialAlignmentProposal::CandidateCapacity(void) const
{
	return max_candidates_;
}

int SpatialAlignmentProposal::CandidateCount(int mol2_index) const
{
	if(mol2_index<0 || mol2_index>=row_count_)
		throw std::out_of_range("SpatialAlignmentProposal row index");
	return counts_[mol2_index];
}

const SpatialNeighbor &SpatialAlignmentProposal::Candidate(int mol2_index,
	int candidate_index) const
{
	if(mol2_index<0 || mol2_index>=row_count_ || candidate_index<0 ||
		candidate_index>=counts_[mol2_index])
		throw std::out_of_range("SpatialAlignmentProposal candidate index");
	return candidates_[mol2_index*max_candidates_+candidate_index];
}

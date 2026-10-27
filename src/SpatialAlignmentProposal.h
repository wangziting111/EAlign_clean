#pragma once

#include <vector>
#include "XYZ.h"

struct SpatialNeighbor
{
	int index;// mol1 residue index
	double distance_square;// 与 mol2[j] 的空间距离平方
};

class SpatialAlignmentProposal
{
public:
	explicit SpatialAlignmentProposal(int max_candidates=6);
	bool BuildFromCoordinates(const XYZ *mol1,int moln1,const XYZ *mol2,int moln2,
		const double *rotmat,double radius);
	double TraceOrderPreservingAlignment(double d0,int moln1,
		std::vector <int> &ali2) const;
	double TraceOrderPreservingAlignment(const std::vector <double> &candidate_scores,
		int moln1,std::vector <int> &ali2) const;
	int CandidateCapacity(void) const;
	int CandidateCount(int mol2_index) const;
	const SpatialNeighbor &Candidate(int mol2_index,int candidate_index) const;

private:
	int max_candidates_;
	int row_count_;
	std::vector <int> counts_;
	std::vector <SpatialNeighbor> candidates_;
};

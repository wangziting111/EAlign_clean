#pragma once

#include "PDB_Residue.h"

// Minimal, self-contained Foldseek 3Di encoder used by EAlign.
// Geometry, encoder parameters and substitution scores follow Foldseek
// commit 718d4217 (GPL-3.0). States are numeric [0, 19]; state 2 is
// Foldseek's invalid/coil fallback.
class Foldseek3Di
{
public:
	static const int STATE_COUNT=20;
	static const int INVALID_STATE=2;

	static void Encode(const PDB_Residue *residues,const XYZ *ca,
		int length,int *states);
	static int Score(int state1,int state2);
	static double CombinedEvidence(double amino_acid_score,double three_di_score);
};

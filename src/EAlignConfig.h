#pragma once

// Evolutionary evidence W(i,j) used by the order-preserving alignment stage.
#define EALIGN_EVIDENCE_AA_AND_3DI 0
#define EALIGN_EVIDENCE_3DI_ONLY 1
#define EALIGN_EVIDENCE_AA_ONLY 2
#define EALIGN_EVIDENCE_AA_AND_CLESUM 3
#define EALIGN_NO_EVIDENCE 4

#ifndef EALIGN_ALIGNMENT_EVIDENCE_MODE
#define EALIGN_ALIGNMENT_EVIDENCE_MODE EALIGN_EVIDENCE_AA_AND_3DI
#endif

#if EALIGN_ALIGNMENT_EVIDENCE_MODE < EALIGN_EVIDENCE_AA_AND_3DI || \
	EALIGN_ALIGNMENT_EVIDENCE_MODE > EALIGN_NO_EVIDENCE
#error "EALIGN_ALIGNMENT_EVIDENCE_MODE must be 0, 1, 2, 3 or 4"
#endif

// Geometry multiplying W(i,j): -1 = V*g, 0 = g only, 1 = V only.
#ifndef EALIGN_ALIGNMENT_GEOMETRY_MODE
#define EALIGN_ALIGNMENT_GEOMETRY_MODE -1
#endif

#if EALIGN_ALIGNMENT_GEOMETRY_MODE < -1 || EALIGN_ALIGNMENT_GEOMETRY_MODE > 1
#error "EALIGN_ALIGNMENT_GEOMETRY_MODE must be -1, 0, or 1"
#endif

// Candidate-level DP refinement after the initial alignment hypotheses.
#ifndef EALIGN_SECOND_REFINEMENT
#define EALIGN_SECOND_REFINEMENT 1
#endif

#if EALIGN_SECOND_REFINEMENT != 0 && EALIGN_SECOND_REFINEMENT != 1
#error "EALIGN_SECOND_REFINEMENT must be 0 or 1"
#endif

// Restrict that second refinement to a band around its input alignment.
#ifndef EALIGN_BOUNDED_DP_REFINEMENT
#define EALIGN_BOUNDED_DP_REFINEMENT 1
#endif

#if EALIGN_BOUNDED_DP_REFINEMENT != 0 && EALIGN_BOUNDED_DP_REFINEMENT != 1
#error "EALIGN_BOUNDED_DP_REFINEMENT must be 0 or 1"
#endif

CXX = g++

ALI3_POSE_REFINEMENT_MODE ?= 0
EVOLUTIONARY_EVIDENCE_MODE ?= 6
EALIGN_ALIGNMENT_EVIDENCE_MODE ?= 0
EALIGN_ALIGNMENT_GEOMETRY_MODE ?= -1
EALIGN_SECOND_REFINEMENT ?= 1
EALIGN_BOUNDED_DP_REFINEMENT ?= 1

CPPFLAGS = -DALI3_POSE_REFINEMENT_MODE=$(ALI3_POSE_REFINEMENT_MODE) \
    -DEVOLUTIONARY_EVIDENCE_MODE=$(EVOLUTIONARY_EVIDENCE_MODE) \
    -DEALIGN_ALIGNMENT_EVIDENCE_MODE=$(EALIGN_ALIGNMENT_EVIDENCE_MODE) \
    -DEALIGN_ALIGNMENT_GEOMETRY_MODE=$(EALIGN_ALIGNMENT_GEOMETRY_MODE) \
    -DEALIGN_SECOND_REFINEMENT=$(EALIGN_SECOND_REFINEMENT) \
    -DEALIGN_BOUNDED_DP_REFINEMENT=$(EALIGN_BOUNDED_DP_REFINEMENT)
CXXFLAGS = -O3

TARGET = binary/ealign

SOURCES := \
    StringUtils.cpp GeometryUtils.cpp ConformationLetters.cpp main.cpp \
    StructureLoader.cpp StructureInput.cpp StructureWriter.cpp PdbUtils.cpp \
    ResidueGeometry.cpp PDB_Residue.cpp PDB_Chain.cpp PDB_File.cpp \
    PdbFeatureExtractor.cpp XYZ.cpp AlignmentGeometry.cpp FragmentChainer.cpp AlignmentRefiner.cpp Kabsch.cpp \
    SequenceScorer.cpp AlignmentPipeline.cpp FragmentAligner.cpp SeedAligner.cpp \
    PairwiseAligner.cpp DynamicProgramming.cpp TMScorer.cpp TMAligner.cpp \
    TMRefiner.cpp EnvironmentScorer.cpp BetaCarbonBuilder.cpp \
    SpatialAlignmentProposal.cpp Foldseek3Di.cpp

OBJECTS := $(SOURCES:.cpp=.o)

.PHONY: all clean clean-objects

all: $(TARGET)

$(TARGET): $(OBJECTS) | binary
	$(CXX) -o $@ $(OBJECTS)
	rm -f $(OBJECTS)

binary:
	mkdir -p $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean-objects:
	rm -f $(OBJECTS)

clean:
	rm -f $(TARGET) $(OBJECTS)

all: AnalyzeSinglePionNTuples AnalyzeLepTauNTuples

AnalyzeSinglePionNTuples: AnalyzeSinglePionNTuples.cc
	g++ AnalyzeSinglePionNTuples.cc -o AnalyzeSinglePionNTuples `root-config --cflags --libs` -lTMVA -lGenVector -O3

AnalyzeLepTauNTuples: AnalyzeLepTauNTuples.cc
	g++ AnalyzeLepTauNTuples.cc -o AnalyzeLepTauNTuples `root-config --cflags --libs` -lTMVA -lGenVector -O3


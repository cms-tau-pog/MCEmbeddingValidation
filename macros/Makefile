all: AnalyzeSinglePionNTuples AnalyzeLepTauNTuples AnalyzeInputMuMuNTuples

AnalyzeInputMuMuNTuples: AnalyzeInputMuMuNTuples.cc
	g++ AnalyzeInputMuMuNTuples.cc -o AnalyzeInputMuMuNTuples `root-config --cflags --libs` -lTMVA -lGenVector -O3

AnalyzeSinglePionNTuples: AnalyzeSinglePionNTuples.cc
	g++ AnalyzeSinglePionNTuples.cc -o AnalyzeSinglePionNTuples `root-config --cflags --libs` -lTMVA -lGenVector -O3

AnalyzeLepTauNTuples: AnalyzeLepTauNTuples.cc
	g++ AnalyzeLepTauNTuples.cc -o AnalyzeLepTauNTuples `root-config --cflags --libs` -lTMVA -lGenVector -O3


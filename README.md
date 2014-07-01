MCEmbeddingValidation
=====================

maps to Validation/MCEmbedding, package for validation on TauAnalysis/MCEmbeddingTools

This package runs mostly out of the box in CMSSW_5_3_17. To run the l+tau validation
code, the electron MVA ID needs to be available:

git cms-addpkg EgammaAnalysis/ElectronTools
cd EgammaAnalysis/ElectronTools/data/
cat download.url | xargs wget      

Optionally, the MVA MET package can be installed to study the MVA PF MET in embedded
events.

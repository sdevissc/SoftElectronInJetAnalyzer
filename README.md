Instructions to get the plots of the different variables used in the TMVA
=========================================================================

cd SoftElecInJet/test
echo name_of_the_rootfile_to_be_used > input.txt
root -l -b
.L RecoToGen.C++
SLPlotter a("input.txt","tag",cond) where tag is just a tag, and cond is the number of the condition as defined in SLPlotter::Condition(int cond,EventPurity *Evt)
a.initialize()
a.getDetailedPlot(int log,TString& process)) , log=1 if you want log plots, process in just a tag
a.getPurityPlots(int log,TString& process), log=1 if you want log plots, process in just a tag

Instructions to train the TMVA
==============================
edit scrit_categorize.sh
./script_categorize.sh

Instructions to use the weight files from TMVA
==============================================

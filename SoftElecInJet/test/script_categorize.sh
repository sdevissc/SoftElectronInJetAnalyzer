#!/bin/bash
source /afs/cern.ch/project/eos/installation/0.1.0-22d/etc/setup.sh
rootfiles=("RecoToGenFill_30_80.root")
ptrange=("pt>2 \&\& pt<5" "pt>5 \&\& pt<10" "pt>10")
etarange=("abs(eta)<0.8" "abs(eta)>0.8 \&\& abs(eta)<1.4" "abs(eta)>1.4")



#################################################################################

indfile=0
echo "------------------------------------------------------------------"
echo "|Creating the weight xml files for the soft electron MVA training|"
echo "------------------------------------------------------------------"
echo "We are in `pwd`"
echo "The different categories are "
indpt_1=0
for pt in "${ptrange[@]}"
	do
	indeta_1=0;
	for eta in "${etarange[@]}"
	do
		echo "BDT_${indpt_1}_${indeta_1}: ${pt} \&\& ${eta}"
		let "indeta_1 += 1";
	done
	let "indpt_1 += 1"
done

echo "And the estimation is done from the following rootfiles"
for file in "${rootfiles[@]}"
do
	echo "${file}"
done
echo "are you happy with this choices? (Y/N)"
read answer
if [[ $answer == "N" || $answer == "n" ]]; then
       	echo "Exiting..."
elif [[ $answer == "Y" || $answer == "y" ]]; then
	eval `scram runtime -sh`
	indfile=0
	for file in "${rootfiles[@]}"
	do
		indpt=0;
		for pt in "${ptrange[@]}"
		do
			indeta=0;
			for eta in "${etarange[@]}"
			do
				echo "->${file}"
				echo "->${pt}"
				echo "->${eta}"
				string="${pt} \&\& ${eta}"
				echo "-->$string"
				sed "s/YYY/${string}/g"  TrainElectronMVA_BDTG_automated.C > temp_${indpt}_${indeta}_${indfile}
				sed "s/XXX/${indpt}_${indeta}_${indfile}/g"  temp_${indpt}_${indeta}_${indfile} > temp_out_${indpt}_${indeta}_${indfile}
				sed "s/ZZZ/${file}/g"  temp_out_${indpt}_${indeta}_${indfile} > temp_out2_${indpt}_${indeta}_${indfile}
				rm temp_${indpt}_${indeta}_${indfile}
				rm temp_out_${indpt}_${indeta}_${indfile}	
				mv temp_out2_${indpt}_${indeta}_${indfile} TrainElectronMVA_BDTG_category.C
				root -l -b <<-EOF
					.L TrainElectronMVA_BDTG_category.C++
					TrainElectronMVA_BTDG_NonTrigV1_Cat1_2012()
				EOF
				let "indeta += 1";
			done	
			let "indpt += 1"
		done
		let "indfile += 1"
	done
#	echo "-----------------------------------------------------------------------------------"
#	echo "| Now getting the plots with uncertainties for each variables used in the training|"
#        echo "-----------------------------------------------------------------------------------"		
#	echo "#include \"TFile.h\"" > macroForPlottingTMVA_Var.C
#        echo "#include \"TTree.h\"" >> macroForPlottingTMVA_Var.C
#        echo "#include \"TBranch.h\"" >> macroForPlottingTMVA_Var.C
#        echo "#include \"TH1F.h\"" >> macroForPlottingTMVA_Var.C
#        echo "#include \"TCanvas.h\"" >> macroForPlottingTMVA_Var.C
#        echo "#include \"TROOT.h\"" >> macroForPlottingTMVA_Var.C
#        echo "#include <iostream>" >> macroForPlottingTMVA_Var.C
#        echo "void test(){" >> macroForPlottingTMVA_Var.C
#	echo "TFile *f[${#ptrange[@]}][${#etarange[@]}][${#rootfiles[@]}];" >> macroForPlottingTMVA_Var.C
#	echo "TTree *t[${#ptrange[@]}][${#etarange[@]}][${#rootfiles[@]}];" >> macroForPlottingTMVA_Var.C
#	echo "TIter *it[${#ptrange[@]}][${#etarange[@]}][${#rootfiles[@]}];" >> macroForPlottingTMVA_Var.C
#	indpt_2=0
#	for i in "${ptrange[@]}"
#        do
#		indeta_2=0	
#                for j in "${etarange[@]}"
#                do	
#			indfile_2=0	
#			for k in "${rootfiles[@]}"
#                        do
#				echo "f[${indpt_2}][${indeta_2}][${indfile_2}] =  new TFile(\"SoftElectronMVA__BDTG_${indpt_2}_${indeta_2}_${indfile_2}.root\");" >> macroForPlottingTMVA_Var.C
#				echo "t[${indpt_2}][${indeta_2}][${indfile_2}] =  (TTree*)f[${indpt_2}][${indeta_2}][${indfile_2}]->Get(\"TrainTree\");" >> macroForPlottingTMVA_Var.C
#				echo "it[${indpt_2}][${indeta_2}][${indfile_2}] = new TIter(t[${indpt_2}][${indeta_2}][${indfile_2}]->GetListOfBranches());" >> macroForPlottingTMVA_Var.C
#				let "indfile_2 += 1"
#			done
#			let "indeta_2 += 1"
#		done
#		let "indpt_2 += 1"
#        done
#	echo "TH1F *h[${#ptrange[@]}][${#etarange[@]}][${#rootfiles[@]}][30][2];" >> macroForPlottingTMVA_Var.C
#	indpt_3=0
#	echo "int indexplot=0;" >> macroForPlottingTMVA_Var.C
#	echo "TCanvas c1;" >> macroForPlottingTMVA_Var.C
#	for i in "${ptrange[@]}"
#        do
#		indeta_3=0
#                for j in "${etarange[@]}"
#                do
#			echo "c1.Print(\"hsimple_${indpt_3}_${indeta_3}.pdf[\");" >> macroForPlottingTMVA_Var.C
#			indfile_3=0
#			for k in "${rootfiles[@]}"
#        		do
#				echo "indexplot=0;" >> macroForPlottingTMVA_Var.C
#                		echo "while (TBranch *br = static_cast<TBranch*>(it[${indpt_3}][${indeta_3}][${indfile_3}]->Next())) {" >> macroForPlottingTMVA_Var.C
#                		echo "	t[${indpt_3}][${indeta_3}][${indfile_3}]->Draw(Form(\"%s>>h_${indpt_3}_${indeta_3}_${indfile_3}_%d_%d\",br->GetName(),indexplot,0),\"classID==0\",\"\");" >> macroForPlottingTMVA_Var.C
#				echo "  t[${indpt_3}][${indeta_3}][${indfile_3}]->Draw(Form(\"%s>>h_${indpt_3}_${indeta_3}_${indfile_3}_%d_%d\",br->GetName(),indexplot,1),\"classID==1\",\"same\");" >> macroForPlottingTMVA_Var.C
#				echo "  h[${indpt_3}][${indeta_3}][${indfile_3}][indexplot][0]= (TH1F*)gROOT->FindObjectAny(Form(\"h_${indpt_3}_${indeta_3}_${indfile_3}_%d_%d\",indexplot,0));">> macroForPlottingTMVA_Var.C
#				echo "  h[${indpt_3}][${indeta_3}][${indfile_3}][indexplot][1]= (TH1F*)gROOT->FindObjectAny(Form(\"h_${indpt_3}_${indeta_3}_${indfile_3}_%d_%d\",indexplot,1));">> macroForPlottingTMVA_Var.C
#				echo "  h[${indpt_3}][${indeta_3}][${indfile_3}][indexplot][0]->SetFillColor(2);" >> macroForPlottingTMVA_Var.C
#				echo "  h[${indpt_3}][${indeta_3}][${indfile_3}][indexplot][0]->SetFillStyle(3004);" >> macroForPlottingTMVA_Var.C
#				echo "  h[${indpt_3}][${indeta_3}][${indfile_3}][indexplot][1]->SetFillColor(3);" >> macroForPlottingTMVA_Var.C
#                                echo "  h[${indpt_3}][${indeta_3}][${indfile_3}][indexplot][1]->SetFillStyle(3005);" >> macroForPlottingTMVA_Var.C
#				echo "  h[${indpt_3}][${indeta_3}][${indfile_3}][indexplot][0]->Draw();" >> macroForPlottingTMVA_Var.C
#                                echo "  h[${indpt_3}][${indeta_3}][${indfile_3}][indexplot][1]->Draw(\"same\");" >> macroForPlottingTMVA_Var.C
#               			echo "	c1.Print(\"hsimple_${indpt_3}_${indeta_3}.pdf\");" >> macroForPlottingTMVA_Var.C
#				echo "indexplot++;" >> macroForPlottingTMVA_Var.C
#        			echo "}" >> macroForPlottingTMVA_Var.C
#				let "indfile_3 += 1"
#			done
#			echo "c1.Print(\"hsimple_${indpt_3}_${indeta_3}.pdf]\");" >> macroForPlottingTMVA_Var.C
#			let "indeta_3 += 1"
#		done
#		let "indpt_3 += 1"
#	done	
#
#        echo "}" >> macroForPlottingTMVA_Var.C
#	rm TrainElectronMVA_BDTG_category*
#	root -l -b <<-EOF
#		.L macroForPlottingTMVA_Var.C++
#		test()
#	EOF
fi

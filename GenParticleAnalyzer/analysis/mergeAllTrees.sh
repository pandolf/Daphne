if [ $# -lt 1 ]
then
  echo "Please provide production name and pass it as first argument: ./mergeAllTrees.sh [productionName]"
else
  ./genTreeMerger $1 QCD_Pt_15to30
fi

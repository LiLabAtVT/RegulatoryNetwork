

# to run this script, cd to the base directory of this respository, 
# then use command:
#
# $ sh ./script/runFilterExpMat.sh
#
for m in `ls ./data/ExpressionFromArray/TimeSeries/*.csv`; do
i=$(basename $m)
echo $i; 
j=${i%.csv}_topVar
echo $j;
Rscript ./script/FilterExprMat.R \
        ./data/ExpressionFromArray/TimeSeries/$i \
        ./data/all_interactions_uniq.csv \
        ./data/Ath_TF_list \
        ./data/ExpressionFromArray/TimeSeries/VarFilter/$j
done;

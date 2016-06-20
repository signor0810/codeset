#for R in 70
#do

R=70
d=160


######################################

MainFile="$1"
cd $MainFile

##Enter d and R
cat > Init <<endofdata
d=$d;
R=$R;
endofdata
wait
cat Init ../ModCountConformationSource.py> CountConformation.py 
wait 

##Loop over samples
for ((i=1; i<=10  ; i++)); do
cd S$i
    python ~/0_SourceCode_WithVersionControl/pizza-2Jul14/src/pizza.py -f ../CountConformation.py  
    wait

cd ..
done

rm Init
cd /userdata3/chien/1_SentJobCenter/1_condorJobs/lmp/OptDNAseparation2
#####################################
#wait
#done

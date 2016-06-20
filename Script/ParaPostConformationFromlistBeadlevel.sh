######################################
R=70
d=160

#########################
MainFile="$1"
cd $MainFile

cat > Init <<endofdata
d=$d;
R=$R;
endofdata
wait
cat Init ../BeadLevelCountConformation6Small.py> CountConformationBeadlevel.py 
wait 

for ((i=1; i<=10  ; i++)); do
cd S$i
    python ~/0_SourceCode_WithVersionControl/pizza-2Jul14/src/pizza.py -f ../CountConformationBeadlevel.py  
    wait

cd ..
done

rm Init
cd /userdata3/chien/1_SentJobCenter/1_condorJobs/lmp/OptDNAseparation2
#####################################




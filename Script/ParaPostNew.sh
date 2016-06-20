CurDirec="/mnt/Data/PostMH"
#pizzaDirec="/home/chien/codeset/pizza-9Oct15/src" -> set in .bash
pythonDirec="/home/chien/codeset/WorkVersion"

######################################
R=70
d=160

#########################
MainFile="$1"
cd $MainFile
for ((i=1; i<=10  ; i++)); do
cd S$i
    #python ~/0_SourceCode_WithVersionControl/pizza-2Jul14/src/pizza.py -f ../../PostSubChainAnl.py $d $R  
    python $pythonDirec/PostSinglePost.py $d $R  
    python $pythonDirec/PostSubChainAnl.py $d $R  

    wait
cd ..
done



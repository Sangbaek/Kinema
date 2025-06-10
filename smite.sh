#
#smite.sh
#Feed inputs to kinema 
#
#!/bin/bash

#Default Inputs
Ee=1000;
the=0.2;
Eep=400;

DEST="$HOME/BLAST/Polite/dat" ;
OUTFILE="E1000Gamma60percent.dat" ;

if [ -f $DEST/$OUTFILE ]; then
    echo "Renew $DEST/$OUTFILE" 
    rm -f $DEST/$OUTFILE
fi

for (( i=1 ; i<=101 ; i+=1 )) ;
    do 
	te=`echo $i | awk '{f=$1/100; print f}'`
        kinema -e $Ee -t $te -f $Eep -v -L 2 >> $DEST/$OUTFILE ;
    done;
        
echo "Output file : $DEST/$OUTFILE"


#
#dkinema.sh
#Feed inputs to dkinema 
#
#!/bin/bash

#Default Inputs
Ee=950;
Q2=0.127;
W=1232;

DEST="$HOME/BLAST/Errors" ;
#OUTFILE="dkinema_E950_Q0127_W_BLAST.dat" ;
OUTFILE="dkinema_E950_Q0127_W_BLAST.dat" ;

if [ -f $DEST/$OUTFILE ]; then
    echo "Renew $DEST/$OUTFILE" 
    rm -f $DEST/$OUTFILE
fi


for (( W=1100 ; W<=1400 ; W+=25 )) ;
    do 
	dkinema -e $Ee -Q $Q2 -W $W -L >> $DEST/$OUTFILE ;
    done;
	
echo "Output file : $DEST/$OUTFILE"











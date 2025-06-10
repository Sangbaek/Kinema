#
#dkinema.sh
#Feed inputs to dkinema 
#
#!/bin/bash

#Default Inputs
Ee=500;
TH=4;
W=1232;

DEST="$HOME/BLAST/Errors" ;
OUTFILE="dkinema.E$Ee.TH$TH.W.SMITE.dat" ;

if [ -f $DEST/$OUTFILE ]; then
    echo "Renew $DEST/$OUTFILE" 
    rm -f $DEST/$OUTFILE
fi

for (( W=1100 ; W<=1400 ; W+=25 )) ;
    do 
	dkinema -e $Ee -t $TH -W $W -L -D >> $DEST/$OUTFILE ;
    done;
	
echo "Output file : $DEST/$OUTFILE"











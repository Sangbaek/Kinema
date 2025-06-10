#
#kinema.sh
#Feed inputs to kinema 
#
#!/bin/bash

#Default Inputs
Ee=1100.;
Q2=0.2;
W=1232.;

#DEST="$HOME/NtoDelta/theory/MAID2001/Unpol" ;
#DEST="$HOME/NtoDelta/theory/MAID2001/oldData" ;
#OUTFILE="WScan_E950_Q0127_kinematics_pi0.dat" ;
DEST= ".";
OUTFILE="hamza.dat";

#if [ -f $DEST/$OUTFILE ]; then
#    echo "Renew $DEST/$OUTFILE" 
#    rm -f $DEST/$OUTFILE
#fi

if [ -f hamza.out ]; then
    echo "Renew hamza.out" 
    rm -f hamza.out
fi


#for (( W=1100 ; W<=1340 ; W+=10 )) ;
#for (( W=1230 ; W<=1234 ; W+=2 )) ;
#    do 
#	kinema -e $Ee -Q $Q2 -W $W -T 14.8 -P -L >> $DEST/$OUTFILE ;
	#./kinema -e $Ee -Q $Q2 -W $W -T 14.8 ;
	#./kinema -e $Ee -Q $Q2 -W $W -T 14.8 -P -L 1 >> hamza.out;
	./dkinema -e $Ee -Q $Q2 -W $W -T 14.8 -P -L 0 | tee -a hamza.out;
#    done;
	
#echo "Output file : hamza.out"











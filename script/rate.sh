#!/bin/bash
#rate.sh

DEST=.
E=5.75
tgt=Pb
#tgt=12C
#tgt=Sn 
opt1=normal
#opt1=cornel
#opt1=sergey
opt2=Rate 
#opt2=
dt=0.05

OPTION=
OPTION2=
mode=xs

if [ $opt2 == Rate ] ; then
   OPTION2="-w -D 7 -d $dt"
   mode=rate
   opt1=${opt1}_7days
fi


#========================================================================
#                         CrossSection
#========================================================================
CrossSection(){

OUTFILE=${mode}_${tgt}_${E}GeV_$opt1.dat

if [ -f $DEST/$OUTFILE ]; then
    echo "Renew $DEST/$OUTFILE"
    rm -f $DEST/$OUTFILE
fi
touch $DEST/$OUTFILE

IMAX=`echo $dt | gawk '{print 4/$1}'`

for (( i=0; i<30 ; i++))
  do
    theta=`echo $i | gawk '{print $1/100}'`
    echo -e -n "calculating E=$E GeV, theta_pi=$theta \r"
    gkinema -e $E -t $theta -T $tgt -p $OPTION | PhProd $OPTION2 >> $DEST/$OUTFILE
  done

for (( i=4; i<50 ; i++))
  do
    theta=`echo $i | gawk '{print $1/10}'`
    echo -e -n "calculating E=$E GeV, theta_pi=$theta \r"
    gkinema -e $E -t $theta -T $tgt -p $OPTION | PhProd $OPTION2 >> $DEST/$OUTFILE
  done

}


SimpleLoop(){

for (( i=0; i<$IMAX ; i++))
  do
    theta=`echo $i $dt | gawk '{print $1*$2}'`
    echo -e -n "calculating E=$E GeV, theta_pi=$theta \r"
    gkinema -e $E -t $theta -T $tgt -p $OPTION | PhProd $OPTION2 >> $DEST/$OUTFILE
  done

}



CrossSection;


#!/bin/bash
#RunPlanLowQ.sh

getGeom(){

    if [ $OPTION == "-s" ] ; then
       OPTION=
    fi

    line 1 $KINFILE | specGeom -x $OPTION 
    line 2 $KINFILE | specGeom -x $OPTION 

}


getXS(){    

    line 4 $KINFILE  | specGeom -x -p 1 | XSisolator $OPTION -1 9 -2 8 -3 10 -4 5
    line 7 $KINFILE  | specGeom -x -p 1 | XSisolator $OPTION -1 18 -2 15 -3 12 -4 5
    line 10 $KINFILE  | specGeom -x -p 1 | XSisolator $OPTION -4 2 

}

WriteFile(){    

    OFILE=dat/Sigma.Q0200.W1221.dat
    line 4 $KINFILE | specGeom -x -p 1 | XSisolator -D -1  7 -2 6 -3 7 -4 4 >  $OFILE
    line 7 $KINFILE | specGeom -x -p 1 | XSisolator -D -1 10 -2 6 -3 7 -4 4 >> $OFILE
    echo -e -n "Output:$OFILE\n"

    if [ $BONUS -eq 1 ] ; then
       echo "Not Yet Ready"
    fi

}

DoFitData(){

    DIR=$HOME/MAMI/Fit
    Q2=Q0200
    th1=53
    th2=20
    rv=rg
    
    FILE=$DIR/FitData.$Q2.t$th1.t$th2.$rv.dat
    OFILE=dat/Fit.$Q2.t$th1.t$th2.$rv.dat
    FitData -d -Q -f $FILE > $OFILE
    echo -e "OutFile=$OFILE"
    
}

#############################################################################
#                                     Help                                  #
#############################################################################
help(){
    echo    " " 
    echo    "RunPlanLowQ.sh [-t] [-f <file] [-g] [-x] [-y] [-w] [-b] [-G]" 
    echo    "    : execute spectrometer geometry/cross section calculations"  
    echo    "      input kinematics file. [def]:kinema.LowQ.dat"
    echo    " " 
    echo -e "  -h --help \t give this help" 
    echo -e "  -t        \t Latex table ouput"
    echo -e "  -f <file> \t input kinematics file [def] kinema.LowQ.dat"
    echo -e "  -g        \t calculate spectrometer geometry"
    echo -e "  -G        \t get spectrometer geometry for runtime estimate"
    echo -e "  -x        \t calculate cross section"
    echo -e "  -y        \t calculate yield"
    echo -e "  -w        \t write cross section to file"
    echo -e "  -b        \t write cross section to file including bonus data"
    echo    " " 
    exit;
}

#############################################################################
#                                   Main                                    #
#############################################################################

KINFILE=kinema.Q0200.dat
DATADIR=dat
YLDFILE=$DATADIR/RunPlanLowQ.Yield.dat
XSFILE=$DATADIR/RunPlanLowQ.XS.dat
GMFILE=$DATADIR/RunPlanLowQ.Geom.dat
OPTION=-s
BONUS=0

while test $# -ne 0; do
  case "$1" in
  -t) OPTION=-t ;;
  -g) getGeom ;;
  -G)  OPTION="-p 2"; getGeom ;;
  -x) getXS ;;
  -y) OPTION=-Y; getXS ;;
  -b) BONUS=1 ; WriteFile;;
  -F) DoFitData;;
  -w) WriteFile ;;
  -f) shift; KINFILE="$1" ;;
  -h | --help) help ;;
  *)  echo "Error: Invarid Option"
      help;;
  esac
  shift
done


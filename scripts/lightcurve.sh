#!/bin/bash

################################################################################
#                                                                              #
# EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System                   #
#                                                                              #
# Automatic lightcurve generation of the detected variable sources             #
#                                                                              #
# Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu                         #
# Inés Pastor Marazuela (2018) - ines.pastor.marazuela@gmail.com              #
#                                                                              #
################################################################################

# bash script generating lightcurves and computing the chi-square and
# Kolmogorov-Smirnov probability of constancy for sources detected by
# Variabilitectron as being variable

# bash $SCRIPTS/lightcurve.sh <path_obs> <path_scripts> <instument> <id> <DL> <TW> <GTR> <BoxSize> <output_log>

# Style functions
################################################################################

title1(){
  message=$1
	i=0; x='===='
	while [[ i -lt ${#message} ]]; do x='='$x; ((++i)); done
  echo -e "\n\t  $message \n\t$x\n"
}

title2(){
  message="$1 $observation"
	i=0; x='----'
	while [[ i -lt ${#message} ]]; do x='-'$x; ((++i)); done
  echo -e "\n\t  $message \n\t$x"
}

title3(){
  message="$1 $observation"
  echo -e "\n # $message"
}


################################################################################
#                                                                              #
# Preliminaries                                                                #
#                                                                              #
################################################################################

#read arguments
start=`date +%s`

# Printing help
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]] ; then
	echo -e "\
	Parameters to use :\n\
	@path       : full path to the observation\n\
  	@scripts    : full path to the scripts folder\n\
	@inst	    : type of detector\n\
	@id         : id number of the detected source within the observation\n\
	@DL         : Detection level used for the variable sources detection\n\
	@TW         : Time window used for the variable sources detection\n\
	@GTR        : Good Time Ratio of acceptability for a time window\n\
	@BS         : Box size in pixels used for the variable sources detection\n\
	@output_log : full path to the document storing the information of the detection\n\
	"
	exit

else

path="$1"
scripts="$2"
inst="$3"
id="$4"
DL="$5"
TW="$6"
GTR="$7"
BS="$8"
output_log="$9"
echo $path
folder=${path:0:-11}
observation=${path: -10}

title1 "Lightcurve Obs. $observation Src. $id"

# Selecting the files and paths

#sums=/mnt/xmmcat/3XMM_data/SumSas_files_4Webcat
#fbks=/mnt/data/Ines/data/fbktsr_dr5

sums=$path
fbks=$path
clean_file=$path/${inst}_clean.fits
gti_file=$path/${inst}_gti.fits
img_file=$path/${inst}_image.fits
nosrc_file=$path/${inst}_sourceless.fits
path_out=$path/lcurve_${TW}_${inst}_${id}

if [ ! -d $path_out ]; then mkdir $path_out; fi
cd $path

# Setting SAS tools
export SAS_ODF=$path
export SAS_CCF=$path/ccf.cif
export HEADAS=/home/mgupta/Heasoft/heasoft-6.29/x86_64-pc-linux-gnu-libc2.31
. $HEADAS/headas-init.sh
. /home/mgupta/sas/xmmsas_20210317_1624/setsas.sh

# Summary file
if [ ! -f $path/*SUM.ASC ]; then 
  cp $sums/*$observation*SUM.ASC $path
fi

sum_file=$(ls $path/*SUM.ASC)

# FBKTSR

fbk_file=$(ls $fbks/*$observation*$inst*FBKTSR*)

# cifbuild
if [ ! -f $path/ccf.cif ]; then cifbuild; fi

################################################################################
#                                                                              #
# Source selection                                                             #
#                                                                              #
################################################################################

title2 "Preliminaries"

cd $path_out

if [ ! -f $img_file ]; then
  evselect table=$clean_file imagebinning=binSize imageset=$img_file withimageset=yes xcolumn=X ycolumn=Y ximagebinsize=80 yimagebinsize=80 -V 0
fi

###
# Reading data from the detected_variable_sources file
###

data=$(cat $path/${DL}_${TW}_${BS}_${GTR}_${inst}/ds9_variable_sources.reg | grep 'text="'$id'"')

###
# Defining source position
###

RAd=$(echo $data | awk '{print $2}' | sed "s/.$//")
DEC=$(echo $data | awk '{print $3}' | sed "s/.$//")
R=$(echo $data | awk '{print $4}' | sed "s/.$//")
srcR=$(echo "$R * 64" | bc)

srcXY=$(ecoordconv imageset=$img_file coordtype=eqpos x=$RAd y=$DEC | tee /dev/tty|grep 'X: Y:'|sed 's/X: Y: //'|sed 's/ /,/g'|sed 's/,//')

# Correcting source and background position

srcexp=$(eregionanalyse imageset=$img_file srcexp="(X,Y) in CIRCLE($srcXY,$srcR)" backval=0.1|tee /dev/tty|grep 'SASCIRCLE'|sed 's/SASCIRCLE: //g')
srcR=$(echo $srcexp | sed "s/(X,Y) in CIRCLE([0-9]*.[0-9]*,[0-9]*.[0-9]*,//" | sed "s/)//")
# arcsec
srcRas=$(bc <<< "scale=2; $srcR * 0.05")

###
# Source name and coordinates
###

# Right ascension in hours
RA=$(bc <<< "scale=5; $RAd / 15")

if [[ $RA == .* ]]; then RA=0$RA; fi
h=$(echo $RA | sed "s/.[0-9]*$//g")
m=$(bc <<< "scale=5; ($RA - $h) * 60" | sed "s/.[0-9]*$//g")
s=$(bc <<< "scale=2; (($RA - $h) * 60 - $m) * 60" | sed "s/.[0-9]*$//")
if [ ${#h} == 1 ]; then h=0$h ; elif [ ${#h} == 0 ] ; then h=00 ; fi
if [ ${#m} == 1 ]; then m=0$m ; elif [ ${#m} == 0 ] ; then m=00 ; fi
if [ ${#s} == 1 ]; then s=0$s ; elif [ ${#s} == 0 ] ; then s=00; fi

# Declination

# DEC=$(echo $srccoord | awk '{print $2}')
if [[ $DEC == -* ]]; then DEC=${DEC#"-"}; link="-"; else link="_"; fi
if [[ $DEC == .* ]]; then DEC=0$DEC; fi
dg=$(echo $DEC | sed "s/.[0-9]*$//g")
am=$(bc <<< "scale=5; ($DEC - $dg) * 60" | sed "s/.[0-9]*$//g")
as=$(bc <<< "scale=2; (($DEC - $dg) * 60 - $am) * 60" | sed "s/.[0-9]*$//")
if [ ${#dg} == 1 ]; then dg=0$dg ; elif [ ${#dg} == 0 ] ; then dg=00 ; fi
if [ ${#am} == 1 ]; then am=0$am ; elif [ ${#am} == 0 ] ; then am=00 ; fi
if [ ${#as} == 1 ]; then as=0$as ; elif [ ${#as} == 0 ] ; then as=00; fi
if [[ $link == "-" ]]; then DEC=$link$DEC; fi

sleep 1
# Source namedg$am$as
src=J$h$m$s$link$dg$am$as
echo -e "\n\t$src\n"

###
# Background extraction region
###

if [ ! -f $nosrc_file ]; then
evselect table=$clean_file withfilteredset=Y filteredset=$nosrc_file destruct=Y keepfilteroutput=T expression="region($fbk_file:REGION,X,Y)" -V 0
fi

bgdXY=$(ebkgreg withsrclist=no withcoords=yes imageset=$img_file x=$RAd y=$DEC r=$srcRas coordtype=EQPOS | grep 'X,Y Sky Coord.' | head -1 | awk '{print$5$6}')
bgdexp="(X,Y) in CIRCLE($bgdXY,$srcR)"

sleep 1
echo -e "\nExtracting data obs. $observation source $src with the following coordinates: \n  Source     : $srcexp\n  Background : $bgdexp"


################################################################################
#                                                                              #
# Lightcurve generation                                                        #
#                                                                              #
################################################################################

title2 "Extracting lightcurve"
if [ ! -f $path/${inst}_gti.wi ]; then gti2xronwin -i $path/${inst}_gti.fits -o $path/${inst}_gti.wi; fi

title3 "Frame time"
frmtime=$(hexdump -e '80/1 "%_p" "\n"' $clean_file | grep -m 1 FRMTIME | awk '{print $3}')
frmtime=$(bc <<< "$frmtime * 0.001")
frmtime=$(bc <<< "scale=3; ($frmtime + 0.001)/1")
if [[ $frmtime == .* ]]; then frmtime=0$frmtime; fi

sleep 1
echo -e "\nFrame time = $frmtime"

title3 "evselect $frmtime s"
evselect table=$clean_file energycolumn=PI expression="$srcexp" withrateset=yes rateset=$path_out/${src}_lc_${frmtime}_src.lc timebinsize=$frmtime maketimecolumn=yes makeratecolumn=yes -V 0
evselect table=$nosrc_file energycolumn=PI expression="$bgdexp" withrateset=yes rateset=$path_out/${src}_lc_${frmtime}_bgd.lc timebinsize=$frmtime maketimecolumn=yes makeratecolumn=yes -V 0

title3 "evselect $TW s"
evselect table=$clean_file energycolumn=PI expression="$srcexp" withrateset=yes rateset=$path_out/${src}_lc_${TW}_src.lc timebinsize=$TW maketimecolumn=yes makeratecolumn=yes -V 0
evselect table=$nosrc_file energycolumn=PI expression="$bgdexp" withrateset=yes rateset=$path_out/${src}_lc_${TW}_bgd.lc timebinsize=$TW maketimecolumn=yes makeratecolumn=yes -V 0

title3 "epiclccorr"
epiclccorr srctslist=$path_out/${src}_lc_${frmtime}_src.lc eventlist=$clean_file outset=$path_out/${src}_lccorr_${frmtime}.lc bkgtslist=$path_out/${src}_lc_${frmtime}_bgd.lc withbkgset=yes applyabsolutecorrections=yes -V 0


sleep 1

title3 "lcstats"
lcstats cfile1="$path_out/${src}_lccorr_${frmtime}.lc" window=$path/${inst}_gti.wi dtnb=$frmtime nbint=1000000 tchat=2 logname="$path_out/${src}_xronos.log"
P=$(lcstats cfile1="$path_out/${src}_lccorr_${frmtime}.lc" window=$path/${inst}_gti.wi dtnb=$frmtime nbint=1000000 tchat=2 logname="$path_out/${src}_xronos.log" | grep "Prob of constancy")

sleep 5
P_chisq=$(echo $P | sed "s/Chi-Square Prob of constancy. //" | sed "s/ (0 means.*//")
P_KS=$(echo $P | sed "s/.*Kolm.-Smir. Prob of constancy //" |  sed "s/ (0 means.*//")

echo -e "Probabilities of constancy : \n\tP_chisq = $P_chisq\n\tP_KS    = $P_KS"

title3 "lcurve"

python3 $scripts/lcurve.py -path $folder -obs $observation -inst $inst -name $src -tw $TW -ft $frmtime -mode medium -pcs $P_chisq -pks $P_KS -n $id

end=`date +%s`
runtime=$((end-start))

###
# Writing output to files, ending script
###

echo -e > $path_out/${src}_region.txt "Source     = $srcexp\nBackground = $bgdexp\nTotal time = $runtime"
echo -e >> $output_log "$observation $id $src $inst $DL $TW $P_chisq $P_KS"

echo -e " # Total time obs. $observation : $runtime seconds"
echo -e "\nObservation $observation ended\nTotal time = $runtime seconds"
fi

#!/bin/sh
# This script adds a record dimension to the raw output from tht FiT event
# tracking code and then concatenates all files together in time
# Concatenation is done in groups to get over a number of files limit

while getopts ":d:v:s:e:" opt; do
        case $opt in
                d) dataV="$OPTARG"
                ;;
                v) versionV="$OPTARG"
                ;;
                s) startyrV="$OPTARG"
                ;;
                e) endyrV="$OPTARG"
                ;;
                \?) echo "Invalid option -$OPTARG"
                ;;
        esac
done

cd /home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/${dataV}_output/${versionV}${startyrV}/raw/

filestart=${dataV}_${versionV}_${startyrV}-${endyrV}
fileend='_4Dobjects.nc'
nyears=$((endyrV-startyrV + 1))
nstart=0
echo $filestart
nend=$((nyears * 366 * 8))
ngroups=$((nend/10000 + 1))
filei=$nstart

echo $nend
echo $ngroups

while [ $filei -le $nend ] 
do
	echo $filei
	ncecat -O -u time ${filestart}_`printf %06d $filei`${fileend} ${filestart}_`printf %06d $filei`${fileend}
	ncks -O --mk_rec_dmn time ${filestart}_`printf %06d $filei`${fileend} ${filestart}_`printf %06d $filei`${fileend}
	echo ncks -O --mk_rec_dmn time ${filestart}_`printf %06d $filei`${fileend} ${filestart}_`printf %06d $filei`${fileend}

	#next line can be used if you accidentally run this twice on some files and get two time dimensions!
	#ncwa -O -a time -d time,0,0 TRMMtest_4th_`printf %05d $filei`_4Dobjects.nc TRMMtest_out1_`printf %05d $filei`_4Dobjects.nc

	((filei=filei+1))
done


groupnum=0

while [ $groupnum -le $ngroups ]
do
        echo ncrcat -h ${filestart}_`printf %02d $groupnum`0000${fileend} group`printf %02d $groupnum`_ts_${filestart}${fileend} 

        ncrcat -O -h ${filestart}_`printf %02d $groupnum`*${fileend} group`printf %02d $groupnum`_ts_${filestart}${fileend}
        ((groupnum=groupnum+1))
done

echo ncrcat group*_ts_${filestart}${fileend} ts_${filestart}${fileend}
ncrcat -O group*_ts_${filestart}${fileend} ts_${filestart}${fileend}
rm -f group*_ts_${filestart}${fileend}




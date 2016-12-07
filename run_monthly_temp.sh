#!/bin/sh

# set defaults:
email=rhwhite@uw.edu	# the email address you want failure messages to be sent to
runThresh=0
runCodeFit=0
runProc=0
runAna=0
runRegrid=0
runFig=0
minGBs=-1

while getopts "t:c:p:a:r:f:d:v:s:e:n:l:m:" opt; do
        case $opt in
		t) runThresh="$OPTARG"
		;;
		c) runCodeFit="$OPTARG"
		;;
		p) runProc="$OPTARG"
		;;
		a) runAna="$OPTARG"
		;;
		r) runRegrid="$OPTARG"
		;;
      	f) runFig="$OPTARG"
		;;
		d) dataV="$OPTARG"
        ;;
        v) versionV="$OPTARG"
        ;;
        s) startyrV="$OPTARG"
        ;;
        e) endyrV="$OPTARG"
        ;;
        n) numfilesV="$OPTARG"
        ;;
		l) numlatsV="$OPTARG"
		;;
		m) minGBs="$OPTARG"
		;;
                \?) echo "Invalid option -$OPTARG"
                ;;
        esac
done

echo 'I have inputs of'
echo 'runThresholding: ' $runThresh
echo 'run Fit Code: ' $runCodeFit
echo 'run Post-processing: ' $runProc
echo 'run Analysis: ' $runAna
echo 'plot Figures: ' $runFig

echo 'Data: ' $dataV
echo 'Version: ' $versionV
echo 'startyr: ' $startyrV
echo 'endyr: ' $endyrV
echo 'numlats/lons: ' $numlatsV
echo 'min gridboxes: ' $minGBs

runThresh(){
	python thresholding.py --Data $dataV --Version $versionV --startyr $startyrV --endyr $endyrV
}

runCodeFit(){
	cd /home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/
	./link_files.sh -d $dataV -v $versionV -s $startyrV -e $endyrV

	mkdir -p ${dataV}_output/${versionV}${startyrV}/raw/
	echo running with namelist namelist_${dataV}_${versionV}_${startyrV}-${endyrV}.dat
	./FiT_Object_analysis_basic_with_NetCDF4.exex namelist_${dataV}_${versionV}_${startyrV}-${endyrV}.dat

	./add_record_concat.sh -d $dataV -v $versionV -s $startyrV -e $endyrV -n $numfilesV 
}

runProc(){
	cd /home/disk/eos4/rachel/git/Python/EventTracking/
	echo python process_full_output_arg.py --Data ${dataV} --Version ${versionV}
	python process_full_output_arg.py --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV}
}

runAna() {

    echo python map_time_speed_monthly_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --minGB ${minGBs}
    python map_time_speed_monthly_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --minGB ${minGBs} || return 1

    echo python map_time_speed_monthly_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -1000 -30 -10 -6 -3 3 6 10 30 --tbound2 -30 -10 -6 -3 3 6 10 30 1000 --unit ms-1 --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --minGB ${minGBs}
    python map_time_speed_monthly_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -1000 -30 -10 -6 -3 3 6 10 30 --tbound2 -30 -10 -6 -3 3 6 10 30 1000 --unit ms-1 --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --minGB ${minGBs} || return 1

}

runRegrid(){
	echo python regrid_time_speed_ann_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV --minGB ${minGBs}
	python regrid_time_speed_ann_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV --minGB ${minGBs}	

	echo python regrid_time_speed_ann_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -1000 -30 -10 -6 -3 3 6 10 30 --tbound2 -30 -10 -6 -3 3 6 10 30 1000 --unit ms-1 --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV --minGB ${minGBs}
	python regrid_time_speed_ann_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -1000 -30 -10 -6 -3 3 6 10 30 --tbound2 -30 -10 -6 -3 3 6 10 30 1000 --unit ms-1 --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV --minGB ${minGBs}
}


runFig(){
	echo python paperplot_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 --tbound2 1 2 5 100 --unit day --Data ${dataV} --Version ${versionV} --anstartyr ${startyrV} --anendyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV --minGB ${minGBs}
	python paperplot_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 --tbound2 1 2 5 100 --unit day --Data ${dataV} --Version ${versionV} --anstartyr ${startyrV} --anendyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV --minGB ${minGBs}

	python paperplot_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -30 -10 -6 -3 3 6 10 --tbound2 -10 -6 -3 3 6 10 30 --unit ms-1 --Data ${dataV} --Version ${versionV} --anstartyr ${startyrV} --anendyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV --minGB ${minGBs}
}

report(){
	[ "$1" != "ok" ] && echo "Failure type $1 with ${dataV} ${versionV} " | mailx -s "$1" $email
	exit 0
}

if [ $runThresh -eq 1 ]; then
	runThresh || report thresholdingerror
fi
if [ $runCodeFit -eq 1 ]; then
	runCodeFit || report fiterror
fi
if [ $runProc -eq 1 ]; then
	runProc || report processingerror
fi
if [ $runAna -eq 1 ]; then		# if running analysis code
	runAna || report analysiserror	# run code, or email to say there was an error
fi
if [ $runRegrid -eq 1 ]; then
	runRegrid || report regriderror
fi
if [ $runFig -eq 1 ]; then		# if running plot figures code
	runFig || report plottingerror
fi




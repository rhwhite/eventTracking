#!/bin/sh
while getopts ":d:v:s:e:n:l:" opt; do
        case $opt in
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
                \?) echo "Invalid option -$OPTARG"
                ;;
        esac
done

echo 'I have inputs of'
echo 'Data: ' $dataV
echo 'Version: ' $versionV
echo 'startyr: ' $startyrV
echo 'endyr: ' $endyrV
echo 'numlats/lons: ' $numlatsV

#echo python regrid_time_speed_ann_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV

#python regrid_time_speed_ann_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV

#python regrid_time_speed_ann_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -1000 -30 -10 -6 -3 3 6 10 30 --tbound2 -30 -10 -6 -3 3 6 10 30 1000 --unit ms-1 --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV

python paperplot_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 --tbound2 1 2 5 100 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV

python paperplot_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -30 -10 -6 -3 3 6 10 --tbound2 -10 -6 -3 3 6 10 30 --unit ms-1 --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV


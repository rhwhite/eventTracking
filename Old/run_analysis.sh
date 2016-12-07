#!/bin/sh
while getopts ":d:v:s:e:n:l:m:" opt; do
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
		m) minGBs="$OPTARG"
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

cd /home/disk/eos4/rachel/git/Python/EventTracking/
echo python process_full_output_arg.py --Data ${dataV} --Version ${versionV}
python process_full_output_arg.py --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV}

echo python extract_time_speed_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --minGB ${minGBs}
python extract_time_speed_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --minGB ${minGBs}

echo python extract_time_speed_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -1000 -30 -10 -6 -3 3 6 10 30 --tbound2 -30 -10 -6 -3 3 6 10 30 1000 --unit ms-1 --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --minGB ${minGBs}
python extract_time_speed_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -1000 -30 -10 -6 -3 3 6 10 30 --tbound2 -30 -10 -6 -3 3 6 10 30 1000 --unit ms-1 --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --minGB ${minGBs}

echo python map_time_speed_annual_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --minGB ${minGBs}
python map_time_speed_annual_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --minGB ${minGBs}

echo python map_time_speed_annual_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -1000 -30 -10 -6 -3 3 6 10 30 --tbound2 -30 -10 -6 -3 3 6 10 30 1000 --unit ms-1 --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --minGB ${minGBs}
python map_time_speed_annual_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -1000 -30 -10 -6 -3 3 6 10 30 --tbound2 -30 -10 -6 -3 3 6 10 30 1000 --unit ms-1 --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --minGB ${minGBs}

echo python regrid_time_speed_ann_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV
python regrid_time_speed_ann_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV

echo python regrid_time_speed_ann_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -1000 -30 -10 -6 -3 3 6 10 30 --tbound2 -30 -10 -6 -3 3 6 10 30 1000 --unit ms-1 --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV
python regrid_time_speed_ann_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -1000 -30 -10 -6 -3 3 6 10 30 --tbound2 -30 -10 -6 -3 3 6 10 30 1000 --unit ms-1 --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV

echo python paperplot_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV
python paperplot_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV} --sumlons $numlatsV --sumlats $numlatsV



#!/bin/sh

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

echo 'I have inputs of'
echo 'Data: ' $dataV
echo 'Version: ' $versionV
echo 'startyr: ' $startyrV
echo 'endyr: ' $endyrV

echo python process_full_output_arg.py --Data ${dataV} --Version ${versionV}
python process_full_output_arg.py --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV}

#echo python map_time_speed_annual_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV}
#python map_time_speed_annual_arg.py --splittype day --speedtspan 4 --tbound1 0 1 2 5 1 --tbound2 1 2 5 100 5 --unit day --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV}

#echo python map_time_speed_annual_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -1000 -30 -10 -6 -3 3 6 10 30 --tbound2 -30 -10 -6 -3 3 6 10 30 1000 --unit ms-1 --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV}
#python map_time_speed_annual_arg.py --splittype maxspeed --speedtspan 4 --tbound1 -1000 -30 -10 -6 -3 3 6 10 30 --tbound2 -30 -10 -6 -3 3 6 10 30 1000 --unit ms-1 --Data ${dataV} --Version ${versionV} --startyr ${startyrV} --endyr ${endyrV}





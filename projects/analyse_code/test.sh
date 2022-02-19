#!/bin/bash

calculate_path=/home/chaugiang/Dropbox/POLARIS-/projects/globule/example/dust
polaris_path=/home/chaugiang/Dropbox/POLARIS-/projects/Bok_globule

runfile () 
{
	grain_size=$4
	f_high_J=$5
	zoom_factor=$6
	mag_type=$1
	file_name=$2
	number=$3

	# file data directory
	file_in_data=$polaris_path/$grain_size/internal_alignment/$f_high_J/with_ISRF/no_larmor/$zoom_factor/$mag_type/data/$file_name

	echo $file_in_data

	# direction of polaris-plot
	file_out_data=$calculate_path/data/$file_name

	file_check=$calculate_path/data/$file_name
	
	#check if file alrealdy exist, if yes, replace
	if test -f "$file_check"; then
		rm $file_check
	fi

	# copy and paste fits file to polaris-plot directory
	cp $file_in_data $file_out_data

	# plot figure using polaris-plot
	polaris-plot globule example dust map $number 1 1 7 --title
}

copyfigure () 
{
	file_name=$1
	grain_size=$2
	f_high_J=$3
	zoom_factor=$4

	file_in_figure=$calculate_path/plots/$file_name
	file_out_figure=$polaris_path/figure/$grain_size/$f_high_J/tau/$zoom_factor/$file_name

	cp $file_in_figure $file_out_figure
}
       
changename ()
{
	grain_size=$1
	size=$2
	f_high_J=$3
	state=$4
	zoom_factor=$5
	zoom=$6
	mag_type=$7

	cd $polaris_path/figure/$grain_size/$f_high_J/tau/$zoom_factor 

	mv emission_map_detector_nr0001.pdf $state"_"$mag_type"_30um_"$size"_"$zoom.pdf
	mv emission_map_detector_nr0002.pdf $state"_"$mag_type"_53um_"$size"_"$zoom.pdf
	mv emission_map_detector_nr0003.pdf $state"_"$mag_type"_89um_"$size"_"$zoom.pdf
	mv emission_map_detector_nr0004.pdf $state"_"$mag_type"_100um_"$size"_"$zoom.pdf
	mv emission_map_detector_nr0005.pdf $state"_"$mag_type"_250um_"$size"_"$zoom.pdf
	mv emission_map_detector_nr0006.pdf $state"_"$mag_type"_450um_"$size"_"$zoom.pdf
}

runscript ()
{
	grain_size=$1
	size=$2
	f_high_J=$3
	state=$4
	zoom_factor=$5
	zoom=$6
	
	
	direct=$polaris_path/figure/$grain_size/$f_high_J/tau/$zoom_factor
	if [ -d "$direct" ]; then
		echo "File exist"
		rm -rf  $direct
	fi

	mkdir $direct

	for j in "para_grain" "Ncl_100" "Ncl_10000"
	do
		number=0
		for i in "polaris_detector_nr0001.fits" "polaris_detector_nr0002.fits" "polaris_detector_nr0003.fits" "polaris_detector_nr0004.fits" "polaris_detector_nr0005.fits" "polaris_detector_nr0006.fits"
		do
			let number+=1
			runfile $j $i $number $grain_size $f_high_J $zoom_factor
		done 

		for i in "emission_map_detector_nr0001.pdf" "emission_map_detector_nr0002.pdf" "emission_map_detector_nr0003.pdf" "emission_map_detector_nr0004.pdf" "emission_map_detector_nr0005.pdf" "emission_map_detector_nr0006.pdf"
		do
			copyfigure $i $grain_size $f_high_J $zoom_factor
		done

		changename $grain_size $size $f_high_J $state $zoom_factor $zoom $j
	done
}

bigrun()
{
	grain_size="amax_"$1"um"
	size="a"$1

	for i in "4000" "1000" "500"
	do
		zoom_factor="zoom_in_"$i"au"
		zoom="z"$i
	
		echo $grain_size, $zoom_factor, $zoom

		runscript $grain_size $size "only_low_J" "low" $zoom_factor $zoom
		runscript $grain_size $size "only_high_J" "high" $zoom_factor $zoom
		runscript $grain_size $size "low+high_J" "lh" $zoom_factor $zoom
	done
}

bigrun 100
bigrun 10

#!bin/sh

for file in ./*.in #./gouge_rupture_inputfiles/*.in
do
	fname="$(basename -- ${file})"
	#https://unix.stackexchange.com/questions/52313/how-to-get-execution-time-of-a-script-effectively
	
	# skip if the simulation case has been already done
	casename=$( echo "$fname" | cut -c15- )
	# echo ${casename%.*}
	if [ -d ${casename%.*}-DataFiles ]; then
  		echo ${casename%.*}-DataFiles already exists. skipping.
  		continue
	fi

	st=`date +%s.%N`
	echo start running ${fname} at `date`
	exec="./rupture_gougepatch_linear_coulomb_friction_law_selfhealing_Gaussnuc"
	mpirun -np 8 --oversubscribe ${exec} ${fname} # --oversubscribe for Mac
	et=`date +%s.%N`
	runtime=$( echo "$et - $st" | bc -l )
	echo case ${fname} runtime ${runtime}s 
done

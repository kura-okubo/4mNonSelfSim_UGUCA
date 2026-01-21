#!/bin/sh
cd ..
ctest
cd simulation_revision_coshselfheal

for file in ./*.in
do
    fname="$(basename -- "${file}")"

    # skip if the simulation case has been already done
    casename=$( echo "$fname" | cut -c15- )
    outdir="${casename%.*}-DataFiles"
    if [ -d "${outdir}" ]; then
        echo "${outdir} already exists. skipping ${fname}."
        continue
    fi

    echo "start running ${fname} at $(date)"

    # ruptype の判定
    if echo "$fname" | grep -q "ruptype=pulse"; then
        # --- ruptype=pulse → selfhealing ---
        exec="./rupture_gougepatch_linear_coulomb_friction_law_selfhealing_Gaussnuc"
        echo "  → ruptype=pulse detected: running ${exec}"
        st=$(date +%s.%N)
        mpirun -np 8 --oversubscribe "${exec}" "${fname}"
        et=$(date +%s.%N)
        runtime=$(echo "$et - $st" | bc -l)
        echo "    -> runtime ${runtime}s"

    elif echo "$fname" | grep -q "ruptype=smoothselfheal"; then
        # --- ruptype=smoothselfheal → coshselfhealing ---
        exec="./rupture_gougepatch_linear_coulomb_friction_law_coshselfhealing_Gaussnuc"
        echo "  → ruptype=smoothselfheal detected: running ${exec}"
        st=$(date +%s.%N)
        mpirun -np 8 --oversubscribe "${exec}" "${fname}"
        et=$(date +%s.%N)
        runtime=$(echo "$et - $st" | bc -l)
        echo "    -> runtime ${runtime}s"

    else
        echo "  [WARN] unknown ruptype in ${fname} – skipping."
    fi

done
###
# @file   CMakeLists.txt
#
# @author David S. Kammer <dkammer@ethz.ch>
# @author Gabriele Albertini <ga288@cornell.edu>
# @author Chun-Yu Ke <ck659@cornell.edu>
#
# @date creation: Fri Feb 5 2021
# @date last modification: Fri Feb 5 2021
#
# @brief  TODO
#
#
# Copyright (C) 2021 ETH Zurich (David S. Kammer)
#
# This file is part of uguca.
#
# uguca is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# uguca is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with uguca.  If not, see <https://www.gnu.org/licenses/>.
####

# gouge patch simulations
add_simulation(rupture_gougepatch_linear_coulomb_friction_law_selfhealing_Gaussnuc rupture_gougepatch_linear_coulomb_friction_law_selfhealing_Gaussnuc.cc)

# copy all input file into build directroy
# ref: https://stackoverflow.com/questions/2183324/looping-over-files-matching-wildcard-in-cmake
file(GLOB files "../simulations_main_casestudy_hpc_paramstudymaster/gouge_rupture_inputfiles*/*.in")
foreach(file ${files})
    get_filename_component(fname "${file}" NAME)
    configure_file(${file} ${fname} COPYONLY)
endforeach()


#configure_file(qsubrunner_casestudy qsubrunner_casestudy COPYONLY)
#configure_file(run_rupgougepatch_casestudy_hpc.sh run_rupgougepatch_casestudy_hpc.sh COPYONLY)
configure_file(qsubrunner_casestudy_paramstudy qsubrunner_casestudy_paramstudy COPYONLY)
configure_file(submit_rupgougepatch_casestudy_hpc_paramstudy.sh submit_rupgougepatch_casestudy_hpc_paramstudy.sh COPYONLY)
#configure_file(run_rupgougepatch_casestudy.sh run_rupgougepatch_casestudy.sh COPYONLY)
#configure_file(plot_checkresult_pythonscript.py plot_checkresult_pythonscript.py COPYONLY)


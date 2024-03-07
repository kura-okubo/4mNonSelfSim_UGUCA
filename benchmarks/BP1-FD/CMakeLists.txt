###
# @file   CMakeLists.txt
#
# @author Chun-Yu Ke <ck659@cornell.edu>
#
# @date creation: Wed Mar 7 2024
# @date last modification: Wed Mar 7 2024
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
add_simulation(BP1-FD BP1-FD.cc)
configure_file(BP1-FD.sh BP1-FD.sh COPYONLY)
configure_file(result_summary.py result_summary.py COPYONLY)
configure_file(generate_output.py generate_output.py COPYONLY)

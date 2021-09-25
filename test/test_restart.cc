/**
 * @file   test_restart.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 *
 * @date creation: Sat Sept 25 2021
 * @date last modification: Sat Sept 25 2021
 *
 * @brief  TODO
 *
 *
 * Copyright (C) 2021 ETH Zurich (David S. Kammer)
 *
 * This file is part of uguca.
 *
 * uguca is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * uguca is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with uguca.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "uca_restart.hh"

#include <iostream>
/*
#include <stdio.h>
#include <random>
#include <sstream>
#include <string>
#include <vector>
*/

using namespace uguca;


int main(){
  std::cout << "start test: test_restart" << std::endl;

  Restart restart_dump("test_restart",".");

  Restart restart_load("test_restart",".");
  
  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}

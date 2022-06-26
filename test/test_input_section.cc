/**
 * @file   test_input_section.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 *
 * @date creation: Sun Jun 26 2022
 * @date last modification: Sun Jun 26 2022
 *
 * @brief  test input section class
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

#include "uca_input_section.hh"

#include <iostream>

using namespace uguca;

int main() {
  std::cout << "start test: input_section" << std::endl;

  InputSection is;

  std::string sname = "toto";
  std::string svalue = "mimi";
  
  std::string dname = "tata";
  double dvalue = 33.3;

  std::string bname = "roro";
  bool bvalue = false;
  
  std::cout << "check 'insert' method" << std::endl;
  is.insert(sname,svalue);
  is.insert(dname,std::to_string(dvalue));
  if (bvalue)
    is.insert(bname,"True");
  else
    is.insert(bname,"False");
  std::cout << "'insert' method works -> success" << std::endl;
    
  std::cout << "check 'has' method" << std::endl;
  if (!is.has(dname)) {
    std::cerr << "should have " << dname << " as key" << std::endl;
    return 1; // failure
  }
  if (is.has(sname+"2")) {
    std::cerr << "should not have this key" << std::endl;
    return 1; // failure
  }
  std::cout << "'has' method works -> success" << std::endl;

  std::cout << "check 'get' method" << std::endl;

  // check get for string
  auto s = is.get<std::string>(sname);
  std::cout << s << std::endl;
  if (s != svalue) {
    std::cerr << "found: " << s << " but should have: " << svalue << std::endl;
    return 1; // failure
  }
  
  // check get for double
  auto d = is.get<double>(dname);
  std::cout << d << std::endl;
  if (std::abs((d-dvalue)/dvalue) > 1e-6) {
    std::cerr << "found: " << d << " but should have: " << dvalue << std::endl;
    return 1; // failure
  }

  // check get for int
  auto i = is.get<int>(dname);
  std::cout << i << std::endl;
  if (i != int(dvalue)) {
    std::cerr << "found: " << i << " but should have: " << int(dvalue) << std::endl;
    return 1; // failure
  }
  
  // check get for unsigned int
  auto ui = is.get<unsigned int>(dname);
  std::cout << ui << std::endl;
  if (ui != (unsigned int)(dvalue)) {
    std::cerr << "found: " << ui << " but should have: "
	      << (unsigned int)(dvalue) << std::endl;
    return 1; // failure
  }

  // check get for boolean
  auto b = is.get<bool>(bname);
  std::cout << b << std::endl;
  if (b != bvalue) {
    std::cerr << "found: " << b << " but should have: "
	      << bvalue << std::endl;
    return 1; // failure
  }
  
  // check that get fails for negative value and unsigned int
  std::string name2 = "titi";
  double value2 = -55.2;
  is.insert(name2,std::to_string(value2));
  bool caught_exception = true;
  try {
    auto nui = is.get<unsigned int>(name2);
    std::cout << nui << std::endl;
    caught_exception = false;
  }
  catch (std::runtime_error &e) {
    std::cout << "caught exception -> success" << std::endl;
  }
  if (!caught_exception) {
    std::cerr << "should have caught exception for negative value -> unsigned int" << std::endl;
    return 1; // failure
  }
  std::cout << "'get' method works -> success" << std::endl;
  
  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}

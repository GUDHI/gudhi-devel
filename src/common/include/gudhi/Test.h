/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#ifndef TEST_H_
#define TEST_H_

#include <list>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>


#define TEST(a) std::cout << "TEST: " << (a) << std::endl
#define TESTMSG(a, b) std::cout << "TEST: " << a << b << std::endl
#define TESTVALUE(a) std::cout << "TEST: " <<  #a << ": " << a << std::endl

/**
 * Class to perform test
 */

class Test {
 private:
  std::string name;
  bool (*test)();

  std::string separation() const {
    return "+++++++++++++++++++++++++++++++++++++++++++++++++\n";
  }

  std::string print_between_plus(std::string& s) const {
    std::stringstream res;
    res << "+++++++++++++++++" << s << "+++++++++++++++++\n";
    return res.str();
  }

 public:
  Test(std::string name_, bool (*test_)()) {
    name = name_;
    test = test_;
  }

  bool run() {
    std::cout << print_between_plus(name);
    return test();
  }

  std::string getName() {
    return name;
  }
};

class Tests {
 private:
  std::list<Test> tests;

 public:
  void add(std::string name_, bool (*test_)()) {
    Test test(name_, test_);
    tests.push_back(test);
  }

  bool run() {
    bool tests_succesful(true);
    std::vector<bool> res;
    for (Test test : tests) {
      res.push_back(test.run());
    }
    std::cout << "\n\n results of tests : " << std::endl;
    int i = 0;
    for (Test t : tests) {
      std::cout << "Test " << i << " \"" << t.getName() << "\"  -->  ";
      if (res[i++]) {
        std::cout << "OK" << std::endl;
      } else {
        std::cout << "Fail" << std::endl;
        tests_succesful = false;
        break;
      }
    }
    return tests_succesful;
  }
};

#endif  // TEST_H_

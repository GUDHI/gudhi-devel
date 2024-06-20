/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>

#include <gudhi/Fields/Z2_field.h>
#include <gudhi/Fields/Z2_field_operators.h>
#include <gudhi/Fields/Zp_field.h>
#include <gudhi/Fields/Zp_field_shared.h>
#include <gudhi/Fields/Zp_field_operators.h>

using Gudhi::persistence_fields::Z2_field_element;
using Gudhi::persistence_fields::Z2_field_operators;
using Gudhi::persistence_fields::Zp_field_element;
using Gudhi::persistence_fields::Shared_Zp_field_element;
using Gudhi::persistence_fields::Zp_field_operators;

template<class Field_element>
void field_element_example(int ini){
  Field_element e = ini;
  std::cout << "Element e initialized with " << ini << "\n";
  std::cout << "Characteristic: " << e.get_characteristic() << "\n";
  std::cout << "e: " << e << "\n";
  std::cout << "e + 3: " << (e + 3) << "\n";
  std::cout << "3 + e: " << (3 + e) << "\n";
  e += 3;
  std::cout << "e += 3: " << e << "\n";
  int t = 3;
  t += e;   //standard integer addition
  std::cout << "3 += e: " << t << "    //standard integer addition\n";
  std::cout << "e - 4: " << (e - 4) << "\n";
  std::cout << "4 - e: " << (4 - e) << "\n";
  e -= 4;
  std::cout << "e -= 4: " << e << "\n";
  t = 4;
  t -= e;   //standard integer addition
  std::cout << "4 -= e: " << t << "    //standard integer substraction\n";
  std::cout << "e * 6: " << (e * 6) << "\n";
  std::cout << "6 * e: " << (6 * e) << "\n";
  e *= 6;
  std::cout << "e *= 6: " << e << "\n";
  t = 6;
  t *= e;   //standard integer addition
  std::cout << "6 *= e: " << t << "    //standard integer multiplication\n";
  std::cout << "Inverse: " << e.get_inverse() << "\n";
}

template<class Field_operator>
void field_operator_example(const Field_operator& op){
  std::cout << "Characteristic: " << op.get_characteristic() << "\n";
  std::cout << "2 + 3: " << op.add(2u, 3u) << "\n";
  std::cout << "3 + 2: " << op.add(3u, 2u) << "\n";
  std::cout << "10 - 4: " << op.substract(10u, 4u) << "\n";
  std::cout << "4 - 10: " << op.substract(4u, 10u) << "\n";
  std::cout << "3 * 6: " << op.multiply(3u, 6u) << "\n";
  std::cout << "6 * 3: " << op.multiply(6u, 3u) << "\n";
  std::cout << "Value of 7: " << op.get_value(7u) << "\n";
  std::cout << "Inverse of 7: " << op.get_inverse(7u) << "\n";
}

int main() {
  std::cout << "=== Example for Z2 field elements ===\n\n";
  field_element_example<Z2_field_element>(2);

  std::cout << "\n=== Example for Z3 field elements ===\n\n";
  field_element_example<Zp_field_element<3> >(5);

  std::cout << "\n=== Example for Z5 field elements ===\n\n";
  Shared_Zp_field_element<>::initialize(5);
  field_element_example<Shared_Zp_field_element<> >(4);

  std::cout << "\n+++ Example for Z2 field operator +++\n\n";
  Z2_field_operators z2op;
  field_operator_example(z2op);

  std::cout << "\n+++ Example for Z3 field operator +++\n\n";
  Zp_field_operators zpop;
  zpop.set_characteristic(3);
  field_operator_example(zpop);

  std::cout << "\n+++ Example for Z5 field operator +++\n\n";
  zpop.set_characteristic(5);
  field_operator_example(zpop);
}

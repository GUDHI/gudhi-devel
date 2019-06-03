/* 
 *    Author: Alejandro Luque cerpa
 */

#include <gudhi/Persistence_intervals.h>
#include <gudhi/Persistent_entropy.h>


#include <iostream>
#include <utility>
#include <vector>

using namespace std;

using Persistence_intervals = Gudhi::Persistence_representations::Persistence_intervals;


int main(int argc, char** argv) {
  if (argc != 2 ) {
    std::cout << "To run this program, please provide the name of a file with a list of files with persistence diagrams \n";
    return 1;
}
  // results vector ignoring infinity bars and not normalized
  std::vector<double> res = Gudhi::Persistent_entropy::persistent_entropy_from_files_list(argv[1], false, -1, false);
  // results vector changing infinity bar to 50 (example value)
  std::vector<double> res2 = Gudhi::Persistent_entropy::persistent_entropy_from_files_list(argv[1], true, 50, true);
  // results vector normalized
  std::vector<double> res3 = Gudhi::Persistent_entropy::persistent_entropy_from_files_list(argv[1], false, -1, true);
  
  // display
  for(std::size_t i=0; i<res.size(); i++)
	std::cout << "The result for diagram " << i << " (standard) is " << res[i] << endl; 
  std::cout << "-------------------------------------------------------------------" << endl;
  for(std::size_t i=0; i<res2.size(); i++)
	std::cout << "The result for diagram " << i << " (infinity_value = 50) is " << res2[i] << endl;
  std::cout << "-------------------------------------------------------------------" << endl;
  for(std::size_t i=0; i<res3.size(); i++)
	std::cout << "The result for diagram " << i << " (standard normalized) is " << res3[i] << endl;
  std::cout << "-------------------------------------------------------------------" << endl;
  return 0;
}

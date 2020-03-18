#include <gudhi/Points_off_io.h>
#include <gudhi/pick_n_random_points.h>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Random.h>

#include <boost/program_options.hpp>

using Point = CGAL::Epick_d< CGAL::Dynamic_dimension_tag>::Point_d;
using Vector_of_points = std::vector<Point>;

const double PI  = 3.141592653589793238463;
#define _USE_MATH_DEFINES

class PointSetGen {
  private:
  	double unirand(){return (double) rand()/(double) RAND_MAX;}
  public:
    void program_options(int argc, char * const argv[]
                         , std::size_t & number_of_points
                         , double & steps
                         , double & end_thresold
                         , int 	  & repetetions
                         , char   & manifold
                         , int 	  & dimension
                         , int    & dim_max
                         , std::string & in_file_name
                         , std::string & out_file_name
                         ) {
      namespace po = boost::program_options;

      po::options_description visible("Allowed options", 100);
      visible.add_options()
        ("help,h", "produce help message")
        ("number,n", po::value<std::size_t>(&number_of_points)->default_value(0),
           "Number of generated point_vector.")
      	("steps,s", po::value<double>(&steps)->default_value(0.1),
       		"Steps of the threshold")
        ("end_thresold,e", po::value<double>(&end_thresold)->default_value(1),
      		"Final threshold for rips complex.")
          
        ("repetetions,r", po::value<int>(&repetetions)->default_value(1),
        	"Num of repetetions of the experiments.")
      	("manifold,m", po::value<char>(&manifold)->default_value('s'),
       		"Type of manifold")
           
        ("dimensions,D", po::value<int>(&dimension)->default_value(2),
         "Dimension of the manifold.")

         ("dim_max,k ", po::value<int>(&dim_max)->default_value(2),
         "Maximum allowed dimension of the Rips complex.")

        ("input_file_name,i", po::value<std::string>(&in_file_name),
         "The input file.")
        ("out_file_name,o", po::value<std::string>(&out_file_name),
         "The output file.");

      po::options_description all;
      all.add(visible);

      po::variables_map vm;
      po::store(po::command_line_parser(argc, argv).
                options(all).run(), vm);
      po::notify(vm);

      if (vm.count("help")) {
        std::cout << std::endl;
        std::cout << "Computes rips complexes of different threshold values, to 'end_thresold', with priodic steps of 'steps' from a n random uniform point_vector on a selected manifold, . \n";
        std::cout << "Strongly collapses all the rips complexes and output the results in out_file. \n";
        std::cout << "The experiments are repeted 'repete' num of times for each threshold value. \n";
        std::cout << "type -m for manifold options, 's' for uni sphere, 'b' for unit ball, 'f' for file. \n";
        std::cout << "type -i 'filename' for Input file option for exported point sample. \n";
        std::cout << std::endl << std::endl;

        std::cout << "Usage: " << argv[0] << " [options]" << std::endl << std::endl;
        std::cout << visible << std::endl;
        std::abort();
      }
    }

   
    void generate_points_sphere(Vector_of_points& W, int nbP, int dim, double radius) {
    	CGAL::Random_points_on_sphere_d<Point> rp(dim+1, radius);
    	for (int i = 0; i < nbP; i++)
    		W.push_back(*rp++);
    }
    // void generate_fibonaci_grid_sphere(Vector_of_points& W, int nbP, int dim, double radius)
    // {

    // }
    
    void generate_grid_2sphere(Vector_of_points& W, int nbP, int r )
    {
    	std::vector<double> coords;
    	int Ncount = 0;
    	double p,v; //the angles phi and psi
    	int M_p;
		
  		double a = (4*PI*pow(r,2))/nbP;
  		double d = sqrt(a);
  		
  		int M_v = PI/d;
  		
  		double d_v = PI/M_v; 
  		double d_p = a/d_v;
  		
  		for( int m = 0; m < M_v ; m++) {
  			v = (PI*(m + 0.5))/M_v;
  			M_p = ((2*PI*sin(v))/d_p);
  			for(int n = 0;  n < M_p ; n++) {
  				p = (2*PI*n)/M_p;
  				coords = {r*sin(v)*cos(p), r*sin(v)*sin(p), r*cos(v)};
  				W.push_back(Point(coords));
  				Ncount += 1;
  			}
  		}
    }

    
    /*
      Generates point sets on <nbSpheres> spheres wedged at origin, sphere can have different radii from <init_radius> with steps of <multiplier_step>
      Number of points on the sphere can also be different from <init_nbP> for the smallest sphere and then multiplied by <multiplier_step>
    */
    void generate_points_wedged_sphere(Vector_of_points& W, int init_nbP, int dim, double init_radius, int multiplier_step, int nbSpheres)  {
      double radius = init_radius;
      int nbP = init_nbP;
      std::vector<double> translation;
      for(int d = 0; d< dim; d++) {
        translation.push_back(0);
      }
      for(int s = 0; s < nbSpheres; s++) {
        CGAL::Random_points_on_sphere_d<Point> rp(dim+1, radius); 
        for (int i = 0; i < nbP; i++) {
          W.push_back(add_point(*rp++, translation, dim)); 
        }
        nbP = nbP*multiplier_step;
        radius = radius*multiplier_step;
        translation.at(dim-1) = (radius - init_radius);  
      }
    }

    void generate_points_concentric_sphere(Vector_of_points& W, int init_nbP, int dim, int init_radius, int multiplier_step, int nbSpheres) {
      double radius = init_radius;
      int nbP = init_nbP;
     
      for(int s = 0; s < nbSpheres; s++) {
        CGAL::Random_points_on_sphere_d<Point> rp(dim+1, radius); 
        for (int i = 0; i < nbP; i++) {
          W.push_back(*rp++); 
        }
        nbP = nbP*(pow(multiplier_step,2));
        radius = radius*multiplier_step;
      }
      
    }
    void generate_points_2annulus(Vector_of_points& W, int nbP, double r_min, double r_max) {
    	double rho, theta;
        double x, y;
	    std::vector<double> coords;
	    double r_min_sq = pow(r_min,2);
	    double r_max_sq = pow(r_max,2);
	
	    srand(time(NULL));
	    for (int i=0; i<nbP; i++) {
	        rho = sqrt((r_max_sq - r_min_sq)*unirand() + r_min_sq );
	        theta = 2.*M_PI*unirand();
	        x = rho*cos(theta);
	        y = rho*sin(theta);

	        coords = {x,y};
	        W.push_back(Point(coords));
	    }
	}
    void generate_points_spherical_shell(Vector_of_points& W, int nbP, double r_min, double r_max) {
    	double rho, phi, theta;
        double x, y, z;
	    std::vector<double> coords;
	    double r_min_cube = pow(r_min,3);
	    double r_max_cube = pow(r_max,3);
	
	    srand(time(NULL));
	    for (int i=0; i<nbP; i++) {
	        rho = cbrt((r_max_cube - r_min_cube)*unirand() + r_min_cube );
	        phi = 2.*M_PI*unirand();
	        theta = acos(1. - 2.*unirand());
	        	 
	        x = rho*sin(theta)*cos(phi);
	        y = rho*sin(theta)*sin(phi);
	        z = rho*cos(theta);

	        coords = {x,y,z};
	        W.push_back(Point(coords));
	    }
	}   

    void generate_points_ball(Vector_of_points& W, int nbP, int dim, double radius) {
    	CGAL::Random_points_in_ball_d<Point> rp(dim, radius); 
    	for (int i = 0; i < nbP; i++)
    		W.push_back(*rp++);
    }

    void generate_points_cube(Vector_of_points& W, int nbP, int dim) {
    	CGAL::Random_points_in_cube_d<Point> rp(dim, 6);
    	for (int i = 0; i < nbP; i++)
    		W.push_back(*rp++);
    }

    void add_point_vectors(Vector_of_points& V, Vector_of_points& U, int nbP, int dim) { // Adds two point vectors of the same size (nbP), by Modifying the first one, V = V+W.
    	for (int i = 0; i < nbP; i++)
    	{
    		V[i] = add_point(V[i], U[i], dim);
    	}
    }

    //returns x = x+y;
    Point add_point(const Point & x, const Point & y, int dim) {
      std::vector<double> coords;
      for(int i =0; i< dim; i++)
        coords.push_back(x[i]+y[i]); 
      return Point(coords);
    }
    void print_point(const Point & x) { 
      std::cout<< "(";
      for(auto & p : x){
        std::cout<< p << ", " ;
      }
      std::cout<< ")" << std::endl;  
    }
    void output_points(Vector_of_points & W, std::string outFile) {
      std::ofstream myfile (outFile, std::ios::app);
        if (myfile.is_open()) {    
          myfile << "OFF"  << "  " << W.size() << " " << W.at(0).size() << std::endl;
          for(auto & v : W){
            for(auto & x : v){
                myfile<< x << " " ;
            }
            myfile << std::endl;
          }
          myfile << "# Tower updated for the additional subcomplex.\n";
          myfile.close();
        }
        else {
            std::cerr << "Unable to open file";
            exit(-1) ;
        } 
    }

    Point noise_point(double noise_param, int dim, double radius) {
	    std::vector<double> noise;
	    for(int d = 0; d< dim; d++){
	    	if(d % 2)
	    		noise.push_back(-noise_param*radius);
	    	else
	    		noise.push_back(noise_param*radius);
	    }
	    return Point(noise);
    }
    //add noise to the points in W.
    void add_noise(Vector_of_points& W, int nbP, int dim, double radius, double noise_param) {
    	Point noise = noise_point(noise_param, dim, radius);
    	for(Vector_of_points::iterator it = W.begin(); it != W.end(); it++ ) {
    		*it = add_point(*it,noise,dim);
    	}
    }
    
    PointSetGen(){}
    ~PointSetGen(){}
};
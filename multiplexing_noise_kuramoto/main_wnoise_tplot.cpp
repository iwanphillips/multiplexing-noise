#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <random>
#include <boost/random.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

typedef boost::numeric::ublas::vector<double> vectord;
typedef boost::numeric::ublas::matrix<double> matrixd;
using namespace boost::numeric::odeint;
using namespace std;

#include "observer_tplot.h"
#include "system_s.h"
#include "noise.h"
#include "ou_noise.h"
#include "stepper.h"

int main(int argc, char **argv)
{
    string file = string(argv[1]);

    boost::property_tree::ptree desc;
    boost::property_tree::json_parser::read_json(file, desc);
    
    const int   n         = 200;
    const int   r_one     = 70;
    float       om_loc    = 0;
    double      g         = atof(argv[5]);
    const float a         = atof(argv[8]);
    string      path_mid  = desc.get<string>("path");
    string      model     = desc.get<string>("model");
    //float       alpha_    = desc.get<float>("PSalpha");
    
    const double sim_time = 199;
    
    cout << "n: "         << n        << '\n' << "r: "         << r_one    << '\n' <<
            "omega loc: " << om_loc   << '\n' << "g: "         << g        << '\n' <<
            "Alpha: "     << a        << '\n' << "Sim time: "  << sim_time << '\n' <<
            "path mid: "  << path_mid  << '\n';
     
    clock_t tStart = clock();
    
    //double tau = alpha_;
    double fac = atof(argv[3]);
    double om_add = atof(argv[4]);
    
    float  l_one = atof(argv[6]);
    double l_two = atof(argv[7]);
    
    double dt = 0.01; 
    matrixd x( n , 2 , 0. );
    matrixd x_init( n , 2 , 0. );
    
    double tau = atof(argv[9]);
    int seed2 = atoi(argv[10]);
    double D = atof(argv[11]);
    
    string path_start  = "/";
    string path_middle = "chim_wnoise_ltwo" + string(argv[7]) + "_seed" + string(argv[10]) + string(argv[2]) + "r70_D" + string(argv[11]) + "_aphase" + string(argv[8]) + "_g" + string(argv[5]) + "_lone" + string(argv[6]);
    string path_end    = "om" + string(argv[4]) + ".txt";
    
    string path =  path_start + path_middle + "_200s_" + path_end;
    string path2 = path_start + path_middle + "_400s_" + path_end;
    string path3 = path_start + path_middle + "_600s_" + path_end;
    string path4 = path_start + path_middle + "_800s_" + path_end;

    ofstream data_out(path); ofstream data_out2(path2); ofstream data_out3(path3); ofstream data_out4(path4);

    
    cout << "total path = " << path << '\n';
    cout << "Starting integration..." << endl;

    size_t D_count = 0;
            
    static uniform_real_distribution<double> u_dist(-0.5,0.5);
    default_random_engine generator(seed2);
            
    for( size_t i = 0 ; i < n ; ++i )
    {
    	double pos = i*2*M_PI/(n-1) - M_PI;
        double r1 = u_dist(generator); double r2 = u_dist(generator);
	size_t j = (i + 100) % 200;
        double pos2 = j*2*M_PI/(n-1) - M_PI;

        x(j,0) = 6*r1*exp(-0.76*pos*pos); x(i,1) = 6*r2*exp(-0.76*pos*pos);
        //x(j,0) = 6*r1*exp(-0.76*pos2*pos2); x(i,1) = 6*r2*exp(-0.76*pos*pos);
	//x(j,0) = 0; x(i,1) = 6*r2*exp(-0.76*pos*pos);
    }
                
    kuramoto_det det( n , l_one , om_loc , g , r_one , a , fac , om_add );
    kuramoto_stoch stoch;
                
    observer_s< kuramoto_det , kuramoto_stoch > obs(  det, stoch, x, data_out,  sim_time );
    observer_s< kuramoto_det , kuramoto_stoch > obs2( det, stoch, x, data_out2, sim_time );
    observer_s< kuramoto_det , kuramoto_stoch > obs3( det, stoch, x, data_out3, sim_time );
    observer_s< kuramoto_det , kuramoto_stoch > obs4( det, stoch, x, data_out4, sim_time );

                
    obs.set_params( l_two , D ); obs2.set_params( l_two , D );
    obs3.set_params( l_two , D ); obs4.set_params( l_two , D );

    obs.reset(); obs2.reset(); obs3.reset(); obs4.reset();

                
    det.set_params(l_two); stoch.set_params(D*l_two);
                
    stoch_RK srk( 1 );
                
    integrate_const( srk, make_pair( det , stoch ), x, 0.0, sim_time, dt, boost::ref( obs ));
    integrate_const( srk, make_pair( det , stoch ), x, 0.0, sim_time, dt, boost::ref( obs2 ));
    integrate_const( srk, make_pair( det , stoch ), x, 0.0, sim_time, dt, boost::ref( obs3 ));
    integrate_const( srk, make_pair( det , stoch ), x, 0.0, sim_time, dt, boost::ref( obs4 ));

    data_out.close();
    
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
    
}

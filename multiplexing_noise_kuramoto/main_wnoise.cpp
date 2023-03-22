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

#include "observer.h"
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
    
    const double sim_time = 500;
    
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
    
    matrixd x( n , 2 , 0. );
    matrixd x_init( n , 2 , 0. );
    
    double tau = atof(argv[9]);
    
    int seed2 = atoi(argv[10]);
    double D = atof(argv[11]);
    double dt = atof(argv[12]);
    int nseed = atoi(argv[13]);
    
    string path_start  = "/" + string(argv[2]) + "/";
    string path_middle = "/1sy1ch_wnoise_ltwo" + string(argv[7]) + "_seed" + string(argv[10]) + "r70_D" + string(argv[11]) + "_aphase" + string(argv[8]) + "_g" + string(argv[5]) + "_lone" + string(argv[6]) + "_dt" + string(argv[12]) + "_nseed" + string(argv[13]);
    string path_end    = "om" + string(argv[4]) + ".txt";
    
    string path = path_start + "500s" + path_middle + "_500s_" + path_end;
    string path2 = path_start + "1000s" + path_middle + "_1000s_" + path_end;
    ofstream data_out(path); ofstream data_out2(path2);
    
    cout << "total path = " << path << '\n';
    cout << "Starting integration..." << endl;
    
    size_t D_count = 0;
    //for( double D  = 0.0 ; D < 10*Dstep + 0.001 ; D += Dstep ){
        for( int seed3  = seed2 ; seed3 < seed2 + nseed*100 + 0.001 ; seed3 += 100 ){
            
            static uniform_real_distribution<double> u_dist(-0.5,0.5);
            default_random_engine generator(seed3);
            
            for( size_t i = 0 ; i < n ; ++i )
            {
                double pos = i*2*M_PI/(n-1) - M_PI;
                double r1 = u_dist(generator); double r2 = u_dist(generator);
                size_t j = (i + 100) % 200;
                double pos2 = j*2*M_PI/(n-1) - M_PI;

		//x(j,0) = 6*r1*exp(-0.76*pos*pos); x(i,1) = 6*r2*exp(-0.76*pos*pos);
                x(j,0) = 6*r1*exp(-0.76*pos2*pos2); x(i,1) = 6*r2*exp(-0.76*pos*pos);
		//x(j,0) = 0; x(i,1) = 6*r2*exp(-0.76*pos*pos);
            }
                
            kuramoto_det det( n , l_one , om_loc , g , r_one , a , fac , om_add );
            kuramoto_stoch stoch;
                
            observer_s< kuramoto_det , kuramoto_stoch > obs(  det, stoch, x, data_out,  sim_time, dt, r_one );
            observer_s< kuramoto_det , kuramoto_stoch > obs2( det, stoch, x, data_out2, sim_time, dt, r_one );
                
            obs.set_params( l_two , D , D_count ); obs2.set_params( l_two , D , D_count );
            obs.reset(); obs2.reset();
                
            det.set_params(l_two); stoch.set_params(D*l_two);
                
            stoch_RK srk( 1 );
                
            integrate_const( srk, make_pair( det , stoch ), x, 0.0, sim_time, dt, boost::ref( obs ));
            data_out  << D << '\t' << tau << '\t' << fac << '\t' << l_two << '\t' << seed3 << '\t' << obs.get_R_inter() << '\t' << obs.get_delta() << '\t' <<  obs.get_om() << '\t' << obs.R().first << '\t' << obs.R().second << '\n';
            integrate_const( srk, make_pair( det , stoch ), x, 0.0, sim_time, dt, boost::ref( obs2 ));
            data_out2 << D << '\t' << tau << '\t' << fac << '\t' << l_two << '\t' << seed3 << '\t' << obs2.get_R_inter() << '\t' << obs2.get_delta() << '\t' << obs2.get_om() << '\t' << obs2.R().first << '\t' << obs2.R().second << '\n';
                
            ++D_count;
        }

    data_out.close();
    
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
    
}

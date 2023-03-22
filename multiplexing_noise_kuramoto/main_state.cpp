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

#include "noise.h"
#include "system.h"
#include "observer.h"
#include "ou_noise.h"

#include "system_s.h"
#include "f_stepper.h"
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
    
    //double dt = 0.01; //0.003; // 0.002;
    matrixd x( n , 2 , 0. );
    matrixd x_init( n , 2 , 0. );
    
    double tau = atof(argv[9]);
    double gamma = 1/tau;
    double alpha_ = tau;
    
    //double Dstep = atof(argv[9]);
    int seed2 = atoi(argv[10]);
    double D = atof(argv[11]);
    double dt = atof(argv[12]);
    int nseed = atoi(argv[13]);
    
    /*int steps = (float) sim_time / dt;
    int n_pts = 16384;
    
    while(steps > n_pts){n_pts = 2*n_pts;}

    float ha = 5 / dt;
    float Q = ha / (2*pow(6.283185, alpha_));
    float Q_d = Q; //Q / pow(dt, 1-alpha_);
        
    int c = 5;
    int max = c*n_pts;
    while(max < 9000000){max += steps; c += 1;}
    vector<float> X;

    int order = 1;
    for( int i=0 ; i<c ; ++i ){
        float X_i[order*n_pts];
        memset(X_i, 0.0, sizeof X_i);
        long utime; utime=(long)time(NULL); long seed=utime;
        f_alpha(n_pts, X_i, Q_d, alpha_, &seed);
        vector<double> Xi (X_i, X_i + sizeof(X_i) / sizeof(int) );
        X.insert(X.end(), Xi.begin(), Xi.end());
    }*/

    long utime; utime=(long)time(NULL); long seed=utime;
    //vector<float> X = ou_euler( gamma, 0, 2*gamma, 0.0, dt, 2000000 , seed );
    vector<float> X = ou_euler( gamma, 0, sqrt(2*gamma), 0.0, dt, 2000000 , seed );
    
    string path_start  = "/" + string(argv[2]) + "/";
    string path_middle = "/chim_ou_ltwo" + string(argv[7]) + "_a" + string(argv[9]) + "_seed" + string(argv[10]) + "r70_D" + string(argv[11]) + "_aphase" + string(argv[8]) + "_g" + string(argv[5]) + "_lone" + string(argv[6]) + "_dt" + string(argv[12]) + "_nseed" + string(argv[13]);
    string path_end    = "om" + string(argv[4]) + ".txt";
    
    string path = path_start + "500s" + path_middle + "_500s_" + path_end;
    string path2 = path_start + "1000s" + path_middle + "_1000s_" + path_end;
    ofstream data_out(path); ofstream data_out2(path2);
    
    cout << "total path = " << path << '\n';
    cout << "Starting integration..." << endl;
    
    size_t D_count = 0;

        for( int seed3  = seed2 ; seed3 < seed2 + nseed*100 + 0.001 ; seed3 += 100 ){
            
            static uniform_real_distribution<double> u_dist(-0.5,0.5);
            default_random_engine generator(seed3);
            
            for( size_t i = 0 ; i < n ; ++i )
            {
                double pos = i*2*M_PI/(n-1) - M_PI;
                double r1 = u_dist(generator); double r2 = u_dist(generator);
                //size_t j = (i + 100) % 200;
                //x(i,0) = 6*r1*exp(-0.003*pos*pos); x(i,1) = 6*r2*exp(-0.003*pos*pos);
                x(i,0) = 0; x(i,1) = 6*r2*exp(-0.76*pos*pos);
                //x(i,0) = 6.2831853*r1; x(i,1) = 6*r2*exp(-0.76*pos*pos);
            }
            
        if(tau < 0.001){
                
            kuramoto_det det( n , l_one , om_loc , g , r_one , a , fac , om_add );
            kuramoto_stoch stoch;
                
            observer_s< kuramoto_det , kuramoto_stoch > obs(  det, stoch, x, data_out,  sim_time, dt, r_one );
            observer_s< kuramoto_det , kuramoto_stoch > obs2( det, stoch, x, data_out2, sim_time, dt, r_one );

                
            obs.set_params( l_two , D , D_count ); obs2.set_params( l_two , D , D_count );
            obs.reset(); obs2.reset();
                
            det.set_params(l_two); stoch.set_params(D*l_two);
                
            stoch_RK srk( 1 );
            //stochastic_euler se( X );
            //fbm_milstein mil( X );
                
            integrate_const( srk, make_pair( det , stoch ), x, 0.0, sim_time, dt, boost::ref( obs ));
            integrate_const( srk, make_pair( det , stoch ), x, 0.0, sim_time, dt, boost::ref( obs2 ));
                
            data_out  << D << '\t' << tau << '\t' << fac << '\t' << l_two << '\t' << seed3 << '\t' << obs.get_R_inter() << '\t' << obs.get_delta() << '\t' <<  obs.get_om() << '\t' << obs.R().first << '\t' << obs.R().second << '\n';
            data_out2 << D << '\t' << tau << '\t' << fac << '\t' << l_two << '\t' << seed3 << '\t' << obs2.get_R_inter() << '\t' << obs2.get_delta() << '\t' << obs2.get_om() << '\t' << obs2.R().first << '\t' << obs2.R().second << '\n';
                
        } else {
                
            k_ring system( X, n, l_one, l_two, om_loc, g, r_one, a, fac, om_add );
                
            observer< k_ring > obs( system, x, data_out,  sim_time, r_one, l_two);
            observer< k_ring > obs2(system, x, data_out2, sim_time, r_one, l_two);
                
            obs.set_params( l_two , D , D_count ); obs2.set_params( l_two , D , D_count );
            obs.reset(); obs2.reset();
                
            system.set_params( l_two , D*l_two );
            
            runge_kutta4< matrixd > rk;
                
            integrate_const(rk, system, x, 0.0, sim_time, dt, boost::ref( obs ));
            integrate_const(rk, system, x, 0.0, sim_time, dt, boost::ref( obs2 ));
                
            data_out  << D << '\t' << tau << '\t' << fac << '\t' << l_two << '\t' << seed3 << '\t' << obs.get_R_inter() << '\t' << obs.get_delta() << '\t' <<  obs.get_om() << '\t' << obs.R().first << '\t' << obs.R().second << '\n';
            data_out2 << D << '\t' << tau << '\t' << fac << '\t' << l_two << '\t' << seed3 << '\t' << obs2.get_R_inter() << '\t' << obs2.get_delta() << '\t' << obs2.get_om() << '\t' << obs2.R().first << '\t' << obs2.R().second << '\n';
        }
            
        ++D_count;
    }

    data_out.close();
    
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
    
}

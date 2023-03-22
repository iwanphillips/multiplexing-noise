#ifndef stepper_h
#define stepper_h

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>

#include <vector>
#include <iostream>
#include <boost/random.hpp>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

double r8_normal_01 ( int *seed );
double r8_uniform_01 ( int *seed );
void timestamp ( );


class stoch_RK
{
    double a21 =   0.66667754298442; double a31 =   0.63493935027993; double a32 =   0.00342761715422; double a41 = - 2.32428921184321;
    double a42 =   2.69723745129487; double a43 =   0.29093673271592; double a51 =   0.25001351164789; double a52 =   0.67428574806272;
    double a53 = - 0.00831795169360; double a54 =   0.08401868181222;

    double q1 = 3.99956364361748; double q2 = 1.64524970733585; double q3 = 1.59330355118722; double q4 = 0.26330006501868;
    
    int m_seed; double m_q; vector<float> m_X;
    
public:

    typedef boost::numeric::ublas::matrix<double> state_type;
    typedef boost::numeric::ublas::matrix<double> deriv_type;
    typedef double value_type;
    typedef double time_type;
    typedef unsigned short order_type;

    typedef boost::numeric::odeint::stepper_tag stepper_category;

    static order_type order( void ) { return 4; }
    
    stoch_RK( double q ) : m_q( q ) {}

    template< class System >
    void do_step( System system , state_type &x , time_type t , time_type dt ) const
    {
        int h = t / dt;
        
        double t1 = t;
        matrixd x1 = x;
        state_type k1( x.size1() , 2 , 0. );
        deriv_type det1( x.size1() , 2 , 0. ) , stoch1( x.size1() , 2 , 0. )  ;
        system.first( x1 , det1 , t1 );
        system.second( x1 , stoch1 , t1 );
        
        int utime; utime=(int)time(NULL); int seed=utime;
        double w1; double w2; double w3; double w4;
        
        //double w1 = r8_normal_01 ( &seed ) * sqrt ( q1 * m_q / dt );
        for( size_t i=0 ; i<x.size1() ; ++i ){
            for( size_t j=0 ; j<x.size2() ; ++j ){
                w1 = r8_normal_01 ( &seed ) * sqrt ( q1 * m_q / dt );
                k1(i,j) = dt * det1(i,j) + dt * stoch1(i,j) * w1;
            }
        }
        
        double t2 = t1 + a21 * dt;
        matrixd x2 = x1 + a21 * k1;
        state_type k2( x.size1() , 2 , 0. );
        deriv_type det2( x.size1() , 2 , 0. ) , stoch2( x.size1() , 2 , 0. ) ;
        system.first( x2 , det2 , t2 );
        system.second( x2 , stoch2 , t2 );
        
        //double w2 = r8_normal_01 ( &seed ) * sqrt ( q2 * m_q / dt );
        for( size_t i=0 ; i<x.size1() ; ++i ){
            for( size_t j=0 ; j<x.size2() ; ++j ){
                w2 = r8_normal_01 ( &seed ) * sqrt ( q2 * m_q / dt );
                k2(i,j) = dt * det2(i,j) + dt * stoch2(i,j) * w2;
            }
        }
        
        double t3 = t1 + a31 * dt  + a32 * dt;
        matrixd x3 = x1 + a31 * k1 + a32 * k2;
        state_type k3( x.size1() , x.size2() , 0. );
        deriv_type det3( x.size1() , 2 , 0. ) , stoch3( x.size1() , 2 , 0. ) ;
        system.first( x3 , det3 , t3 );
        system.second( x3 , stoch3 , t3 );
        
        //double w3 = r8_normal_01 ( &seed ) * sqrt ( q3 * m_q / dt );
        for( size_t i=0 ; i<x.size1() ; ++i ){
            for( size_t j=0 ; j<x.size2() ; ++j ){
                w3 = r8_normal_01 ( &seed ) * sqrt ( q3 * m_q / dt );
                k3(i,j) = dt * det3(i,j) + dt * stoch3(i,j) * w3;
            }
        }
        
        double t4 = t1 + a41 * dt  + a42 * dt  + a43 * dt;
        matrixd x4 = x1 + a41 * k1 + a42 * k2 + a43 * k3;
        state_type k4( x.size1() , x.size2() , 0. );
        deriv_type det4( x.size1() , 2 , 0. ) , stoch4( x.size1() , 2 , 0. );
        system.first( x4 , det4 , t4 );
        system.second( x4 , stoch4 , t4 );
        
        //double w4 = r8_normal_01 ( &seed ) * sqrt ( q4 * m_q / dt );
        for( size_t i=0 ; i<x.size1() ; ++i ){
            for( size_t j=0 ; j<x.size2() ; ++j ){
                w4 = r8_normal_01 ( &seed ) * sqrt ( q4 * m_q / dt );
                k4(i,j) = dt * det4(i,j) + dt * stoch4(i,j) * w4;
            }
        }

                x = x1 + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4;
    }
};

//stoch_RK( vector<float> X , int seed , double q ) : m_X ( X ) , m_seed( seed ) , m_q( q ) {}

/*
double r8_normal_01 ( int *seed )
{
# define R8_PI 3.141592653589793

  double r1;
  double r2;
  static int seed2 = 0;
  static int seed3 = 0;
  static int used = 0;
  double v1;
  static double v2 = 0.0;

  if ( ( used % 2 ) == 1 )
  {
    if ( *seed != seed2 )
    {
      used = 0;
      seed2 = 0;
      seed3 = 0;
      v2 = 0.0;
    }
  }

  if ( ( used % 2 ) == 0 )
  {
    r1 = r8_uniform_01 ( seed );

    if ( r1 == 0.0 )
    {
      cerr << "\n";
      cerr << "R8_NORMAL_01 - Fatal error!\n";
      cerr << "  R8_UNIFORM_01 returned a value of 0.\n";
      exit ( 1 );
    }

    seed2 = *seed;
    r2 = r8_uniform_01 ( seed );
    seed3 = *seed;
    *seed = seed2;

    v1 = sqrt ( - 2.0 * log ( r1 ) ) * cos ( 2.0 * R8_PI * r2 );
    v2 = sqrt ( - 2.0 * log ( r1 ) ) * sin ( 2.0 * R8_PI * r2 );
  }

  else
  {
    v1 = v2;
    *seed = seed3;
  }

  used = used + 1;

  return v1;
# undef R8_PI
}


double r8_uniform_01 ( int *seed )
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}


void timestamp ( )
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}

*/
 
#endif /* stepper_h */



#ifndef system_s_h
#define system_s_h

typedef boost::numeric::ublas::matrix<double> state_type;

struct kuramoto_det
{
    int     m_N; double  m_lamb1; double  m_sig12; double m_D;
    vectord m_om1; vectord m_om2;  int m_range; float m_alpha;
    double m_fac; double m_om_add;
    
    kuramoto_det( int N, double lamb1, float om_loc, double g,
                  int range, float alpha, double fac, double om_add ) :
    m_N ( N ) , m_lamb1 ( lamb1 ) , m_om1 ( N , 0.0 ) , m_om2 ( N , 0.0 ) ,
    m_range ( range ) , m_alpha ( alpha ) , m_fac ( fac ) , m_om_add ( om_add )
    {
        create_frequencies( om_loc, g );
    }
    
    void create_frequencies( float om_loc, double g )
    {
        boost::mt19937 rng;
        boost::cauchy_distribution<> cauchy( om_loc , g );
        boost::variate_generator< boost::mt19937&, boost::cauchy_distribution<> > gen( rng , cauchy );
        generate( m_om1.begin() , m_om1.end() , gen );
        generate( m_om2.begin() , m_om2.end() , gen );
    }
    
    void set_params( double lamb2 ) { m_sig12 = lamb2; }
    
    void operator() ( const state_type &x , state_type &dxdt , const double t )
    {
        state_type loc_coup( m_N , 2 , 0. );
        
        for( int i=0 ; i< x.size1() ; ++i ){
            for( int k=0; k< x.size1() ; ++k ){
                
                //loc_coup(i,0) += sin( x(i,0) - x(k,0) + m_alpha );
                
                float dist = abs(k-i);
                dist = abs(dist - round(dist/( (float) m_N ) ) * m_N);
                if(dist <= m_range){ // && dist > 0
                    loc_coup(i,0) += sin( x(i,0) - x(k,0) + m_alpha );
                    loc_coup(i,1) += sin( x(i,1) - x(k,1) + m_alpha );
                }
            }
            //dxdt(i,0) = m_om1(i) - (m_lamb1/(2*m_range))*loc_coup(i,0) + m_sig12*sin( x(i,1) - x(i,0) );
            //dxdt(i,1) = m_om2(i) + m_om_add - (m_fac*m_lamb1/(2*m_range))*loc_coup(i,1) - m_sig12*sin( x(i,1) - x(i,0) );
            
            dxdt(i,0) = 0 - (m_lamb1/(2*m_range))*loc_coup(i,0) + m_sig12*sin( x(i,1) - x(i,0) );
            dxdt(i,1) = 0 + m_om_add - (m_fac*m_lamb1/(2*m_range))*loc_coup(i,1) - m_sig12*sin( x(i,1) - x(i,0) );
        }
    }
};

struct kuramoto_stoch
{
    double m_D;
    
    void set_params( double D ) { m_D = D; }
    
    void operator()( const state_type &x , state_type &dxdt , const double t  )
    {
        for( int i=0 ; i< x.size1() ; ++i ){
            dxdt(i,0) = m_D*sin( x(i,1) - x(i,0) );
            dxdt(i,1) = m_D*sin( x(i,0) - x(i,1) );
        }
    }
    
    void order2( const state_type &x , state_type &dxdt , const double t  )
    {
        state_type loc_coup( x.size1() , 2 , 0. );
        
        for( int i=0 ; i< x.size1() ; ++i ){
            
            for( int k=0; k< x.size1() ; ++k ){
                loc_coup(i,0) += sin( x(k,0) - x(k,1) );
                loc_coup(i,1) += sin( x(k,1) - x(k,0) );
            }
            dxdt(i,0) = 2*m_D*m_D*loc_coup(i,0)*cos( x(i,1) - x(i,0) );
            dxdt(i,1) = 2*m_D*m_D*loc_coup(i,1)*cos( x(i,0) - x(i,1) );
        }
    }
}; // Why does to milstein method not work for the deterministic case when the
// second order has nothing to do with the deterministic case.

#endif /* system_s_h */

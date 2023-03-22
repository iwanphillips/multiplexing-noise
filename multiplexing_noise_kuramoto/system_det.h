#ifndef system_h
#define system_h

struct k_ring
{
    int     m_N; double  m_lamb1; double  m_lamb2;
    vectord m_om1; vectord m_om2; int m_range;
    float   m_alpha;
    double m_fac; double m_om_add;
    
    k_ring( int N, double lamb1, double lamb2, float om_loc, double g,
            int range, float alpha, double fac, double om_add ) :
    m_N ( N ) , m_lamb1 ( lamb1 ) , m_lamb2 ( lamb2 ) , m_om1 ( N , 0.0 ) ,
    m_om2 ( N , 0.0 )  , m_range ( range ) , m_alpha ( alpha ) , m_fac ( fac ) ,
    m_om_add ( om_add )
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
    
    void set_params( double lamb2 ) { m_lamb2 = lamb2; }
    
    void operator() (const matrixd &x, matrixd &dxdt, const double t )
    {
        matrixd loc_coup( m_N , 2 , 0. );
        
        for( int i=0 ; i<m_N ; ++i ){
            for( int k=0; k<m_N ; ++k ){
                
                float dist = abs(k-i);
                dist = abs(dist - round(dist/( (float) m_N ) ) * m_N);
                if(dist <= m_range && dist > 0){
                    
                    loc_coup(i,0) += sin( x(i,0) - x(k,0) + m_alpha );
                    loc_coup(i,1) += sin( x(i,1) - x(k,1) + m_alpha );
                }
            }
            dxdt(i,0) = 0 - (m_lamb1/(2*m_range))*loc_coup(i,0) + m_lamb2*sin( x(i,1) - x(i,0) );
            dxdt(i,1) = 0 + m_om_add - (m_fac*m_lamb1/(2*m_range))*loc_coup(i,1) - m_lamb2*sin( x(i,1) - x(i,0) );
        }
    }
};

#endif /* system_h */





//double a = i*2*M_PI/(m_N-1) - M_PI;
//double b = k*2*M_PI/(m_N-1) - M_PI;
//m_loc_coup(i,0) += (1/(2*M_PI))*(1 + 0.995*cos(a-b))*sin( x(i,0) - x(k,0) + m_alpha );
//m_loc_coup(i,1) += (1/(2*M_PI))*(1 + 0.995*cos(a-b))*sin( x(i,1) - x(k,1) + m_alpha );

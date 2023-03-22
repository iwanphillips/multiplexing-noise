#ifndef observer_h
#define observer_h

pair< double , double > calc_mean_field( const vectord &x , int range = 0 , int i = 1 )
{
    size_t n = x.size();
    double cos_sum = 0.0 , sin_sum = 0.0;
     
    for( int j=0 ; j<n ; ++j){
        float dist = abs(j-i);
        dist = abs(dist - round(dist/( (float) n ) ) * n);
             
        if(dist <= range && dist > 0){
        cos_sum += cos( x[j] );
        sin_sum += sin( x[j] );
        }
    }
    cos_sum /= double( 2 * (float) range );
    sin_sum /= double( 2 * (float) range );
         
     
    double K = sqrt( cos_sum * cos_sum + sin_sum * sin_sum );
    double Theta = atan2( abs(sin_sum) , abs(cos_sum) );
  
    return make_pair( K , Theta );
}

pair< double , double > calc_mean_alltoall( const vectord &x )
{
    size_t n = x.size();
    double cos_sum = 0.0 , sin_sum = 0.0;
     
    for( int j=0 ; j<n ; ++j){
        cos_sum += cos( x[j] );
        sin_sum += sin( x[j] );
    }
    
    cos_sum /= double( n );
    sin_sum /= double( n );
         
    double K = sqrt( cos_sum * cos_sum + sin_sum * sin_sum );
    double Theta = atan2( abs(sin_sum) , abs(cos_sum) );
  
    return make_pair( K , Theta );
}

double sync_error( const matrixd &x )
{
    size_t n = x.size1();
    double psi_sum = 0.0;
    
    for( size_t i=0 ; i<n ; ++i )
    {
        double psi = fmod(abs(x(i,1)-x(i,0)), 2*M_PI);
        if(psi > M_PI){psi = 2*M_PI - psi;}
        psi_sum += psi * psi;
    }
    
    psi_sum /= n;

    return psi_sum;
}

double om ( const matrixd &dxdt )
{
    size_t n = dxdt.size1();
    double psi_sum = 0.0;
    
    for( size_t i=0 ; i<n ; ++i )  // I may need modulus or similar here
    {
        double psi = dxdt(i,1)-dxdt(i,0);
        psi_sum += psi * psi;
    }
    
    psi_sum /= n;

    return psi_sum;
}

template<class system>
struct observer
{
    ostream&        m_outfile;
    const double    m_sim_time;
    string          m_output;
    
    double m_om;
    double m_delta;
    double m_R_int;
    size_t m_count;
    int m_lap_count;
    double m_lamb2;
    double m_D;
    size_t m_D_count;
    int m_r_one;
    
    double m_R1;
    double m_R2;
    
    matrixd m_dxdt;
    system m_odefun;
    
    observer( system odefun, matrixd x, ostream &out, const float &sim_time, int r_one, double lamb2 = 0.0 ) :
        m_odefun(odefun) , m_dxdt(x) , m_om( 0.0 ) , m_count( 0 ) , m_lap_count( 0 ) ,
        m_delta( 0 ) , m_R1( 0.0 ), m_R2( 0.0 ), m_R_int( 0.0 ), m_outfile( out ) , m_sim_time( sim_time ) ,
        m_r_one( r_one ) , m_lamb2( lamb2 ) { }
    
    void set_params( double lamb2 , double D , double D_count ) { m_lamb2 = lamb2; m_D = D; m_D_count = D_count; }
    
    void operator()( matrixd &x, double t)
    {
        m_odefun( x, m_dxdt, t );
        
        /*if(m_count % 10 == 0 ){              // for animations
            for(int i=0; i < x.size1(); ++i){
                m_outfile << t << '\t' << m_dxdt(i,1) - m_dxdt(i,0) << '\n';
            }
        }
        ++m_count; */
            
        /*if(m_count % 100 == 0 ){              // for tplot
            for(int i=0; i < x.size1(); ++i){
                m_outfile << m_D << '\t' << m_lamb2 << '\t' << t << '\t' << i
                << '\t' << x(i,1) << '\t' << x(i,0) << '\t' << m_dxdt(i,1) <<
                '\t' << m_dxdt(i,0) << '\n'; 
            }
        }
        ++m_count;*/
        
        if( t > m_sim_time-30 ){
            ++m_count;
            
            double mean = sync_error( x );
            double der_mean = om ( m_dxdt );
            
            m_delta += mean;
            m_om += der_mean;
            
            double R_inter = 0;
            double cos_sum; double sin_sum;
            for(int i=0; i < x.size1(); ++i){
                cos_sum = cos( x(i,0) ) + cos( x(i,1) );
                sin_sum = sin( x(i,0) ) + sin( x(i,1) );
                cos_sum /= 2; sin_sum /= 2;
                R_inter += sqrt( cos_sum * cos_sum + sin_sum * sin_sum );
            }
            /*for(int i=0; i < x.size1(); ++i){
                
                double cos_sum = 0; double sin_sum = 0;
                for(int j=0; j < x.size1(); ++j){
                    cos_sum += cos( x(i,0) ) + cos( x(j,1) );
                    sin_sum += sin( x(i,0) ) + sin( x(j,1) );
                    cos_sum /= 2*double( x.size1() ); sin_sum /= 2*double( x.size1() );
                }
                R_inter = sqrt( cos_sum * cos_sum + sin_sum * sin_sum );
            }*/
            R_inter /= double( x.size1() );
            m_R_int += R_inter;
            
            double r1 = 0; double r2 = 0;
            
            boost::numeric::ublas::matrix_column<matrixd> x1 (x,0);
            boost::numeric::ublas::matrix_column<matrixd> x2 (x,1);
            
            for(int i=0; i < x.size1(); ++i){
                
                r1 += calc_mean_field( x1 , m_r_one , i ).first;
                r2 += calc_mean_field( x2 , m_r_one , i ).first;
            }
            
            m_R1 += r1 / x.size1();
            m_R2 += r2 / x.size1();
            
            //m_R1 += calc_mean_alltoall( x1 ).first;
            //m_R2 += calc_mean_alltoall( x2 ).first;
                
            /*if(t > m_sim_time-0.01 && m_D_count % 5 == 0 ){
                for(int i=0; i < x.size1(); ++i){
                    
                    boost::numeric::ublas::matrix_column<matrixd> x1 (x,0);
                    boost::numeric::ublas::matrix_column<matrixd> x2 (x,1);
                    //pair< double , double > mean1 = calc_mean_field( x1 , m_r_one , i );
                    //pair< double , double > mean2 = calc_mean_field( x2 , m_r_one , i );
                    pair< double , double > mean1 = calc_mean_alltoall( x1 );
                    pair< double , double > mean2 = calc_mean_alltoall( x2 );
                    
                    m_outfile << m_D << '\t' << m_lamb2 << '\t' << t << '\t' << i << '\t'
                    << x(i,1) << '\t' << x(i,0) << '\t' << m_dxdt(i,1) <<'\t' <<
                    m_dxdt(i,0) << '\t' << mean1.first << '\t' << mean2.first << '\n';
                }
            }*/
        }
    }
    
    double get_om( void ) const { return ( m_count != 0 ) ? sqrt( m_om / double( m_count ) ) : 0.0 ; }
    
    double get_delta( void ) const { return ( m_count != 0 ) ? sqrt( m_delta / double( m_count ) ) : 0.0 ; }
    
    double get_R_inter( void ) const { return ( m_R_int / double( m_count ) ); }
    
    pair< double , double > R( void ){ return make_pair( ( m_R1 / double( m_count ) ) , ( m_R2 / double( m_count ) ) ); }
    
    double get_lap( void ) const { return m_lap_count; }
    
    void reset( void ) { m_om = 0.0; m_count = 0; m_lap_count = 0; m_delta = 0; m_R1 = 0; m_R2 = 0; m_R_int = 0; }
};



template<class det, class stoch>
struct observer_s
{
    ostream&        m_outfile;
    const double    m_sim_time;
    string          m_output;
    
    double m_om;
    double m_delta;
    double m_R_int;
    size_t m_count;
    int m_lap_count;
    double m_lamb2;
    double m_D;
    size_t m_D_count;
    double m_dt;
    int m_r_one;
    
    double m_R1;
    double m_R2;

    matrixd m_dxdt_d;
    matrixd m_dxdt_s;
    det m_odedet;
    stoch m_odestoch;
    
    observer_s( det odedet, stoch odestoch, matrixd x, ostream &out, const float &sim_time, double dt, int r_one ) :
        m_dxdt_d(x) , m_dxdt_s(x) , m_om( 0.0 ) , m_count( 0 ) , m_lap_count( 0 ) , m_odedet(odedet) ,
        m_R1( 0.0 ), m_R2( 0.0 ), m_odestoch(odestoch) , m_delta( 0 ) , m_R_int( 0.0 ), m_outfile( out ) ,
        m_sim_time( sim_time ) , m_dt( dt ) , m_r_one( r_one ) { }
    
    void set_params( double lamb2 , double D , double D_count ) { m_lamb2 = lamb2; m_D = D; m_D_count = D_count; }
    
    void operator()( matrixd &x, double t )
    {
        m_odedet( x, m_dxdt_d, t );
        m_odestoch( x, m_dxdt_s, t );
        matrixd dxdt = m_dxdt_d + m_dxdt_s;
        
        if(t > m_sim_time-30){
            ++m_count;
            
            double mean = sync_error( x );
            double der_mean = om ( dxdt );
            
            m_delta += mean;
            m_om += der_mean;
            
            double R_inter = 0;
            double cos_sum; double sin_sum;
            
            for(int i=0; i < x.size1(); ++i){
                cos_sum = cos( x(i,0) ) + cos( x(i,1) );
                sin_sum = sin( x(i,0) ) + sin( x(i,1) );
                cos_sum /= 2; sin_sum /= 2;
                R_inter += sqrt( cos_sum * cos_sum + sin_sum * sin_sum );
            }

            R_inter /= double( x.size1() );
            m_R_int += R_inter;
            
            double r1 = 0; double r2 = 0;
            
            boost::numeric::ublas::matrix_column<matrixd> x1 (x,0);
            boost::numeric::ublas::matrix_column<matrixd> x2 (x,1);
            
            for(int i=0; i < x.size1(); ++i){
                
                r1 += calc_mean_field( x1 , m_r_one , i ).first;
                r2 += calc_mean_field( x2 , m_r_one , i ).first;
            }
            
            m_R1 += r1 / x.size1();
            m_R2 += r2 / x.size1();
            
        }
    }
    
    double get_om( void ) const { return ( m_count != 0 ) ? sqrt( m_om / double( m_count ) ) : 0.0 ; }
    
    double get_delta( void ) const { return ( m_count != 0 ) ? sqrt( m_delta / double( m_count ) ) : 0.0 ; }
    
    double get_R_inter( void ) const { return ( m_R_int / double( m_count ) ); }
    
    pair< double , double > R( void ){ return make_pair( ( m_R1 / double( m_count ) ) , ( m_R2 / double( m_count ) ) ); }
    
    void reset( void ) { m_om = 0.0; m_count = 0; m_lap_count = 0; m_delta = 0; m_R1 = 0; m_R2 = 0; m_R_int = 0; }
};


#endif /* observer_h */

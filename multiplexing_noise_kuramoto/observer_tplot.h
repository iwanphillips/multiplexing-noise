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
    
    psi_sum /= n; // sqrt( psi_sum );

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

template<class system> // not necassary
struct observer
{
    ostream&        m_outfile;
    const double    m_sim_time;
    string          m_output;
    
    double m_om;
    double m_delta;
    size_t m_count;
    int m_lap_count;
    double m_lamb2;
    double m_D;
    
    double m_R1;
    double m_R2;
    
    matrixd m_dxdt;
    system m_odefun;
    
    observer( system odefun, matrixd x, ostream &out, const float &sim_time, double lamb2 = 0.0 ) :
        m_odefun(odefun) , m_dxdt(x) , m_delta( 0 ) , m_om( 0.0 ) , m_count( 0 ) , m_lap_count( 0 ) ,
        m_R1( 0.0 ), m_R2( 0.0 ), m_outfile( out ) , m_sim_time( sim_time ) , m_lamb2( lamb2 ) { }
    
    void set_params( double lamb2 , double D ) { m_lamb2 = lamb2; m_D = D; }
    
    void operator()( matrixd &x, double t)
    {
        m_odefun( x, m_dxdt, t );
        
        /*if(m_count % 10 == 0 ){              // for animations
            for(int i=0; i < x.size1(); ++i){
                m_outfile << t << '\t' << m_dxdt(i,1) - m_dxdt(i,0) << '\n';
            }
        }
        ++m_count; */
            
        if(m_count % 100 == 0 ){
            for(int i=0; i < x.size1(); ++i){
                m_outfile << m_D << '\t' << m_lamb2 << '\t' << t << '\t' << i
                << '\t' << x(i,1) << '\t' << x(i,0) << '\t' << m_dxdt(i,1) <<
                '\t' << m_dxdt(i,0) << '\n';
            }
        }
        ++m_count;
    }
    
    double get_om( void ) const { return ( m_count != 0 ) ? sqrt( m_om / double( m_count ) ) : 0.0 ; }
    
    double get_delta( void ) const { return ( m_count != 0 ) ? sqrt( m_delta / double( m_count ) ) : 0.0 ; }
    
    pair< double , double > R( void ){ return make_pair( ( m_R1 / double( m_count ) ) , ( m_R2 / double( m_count ) ) ); }
    
    double get_lap( void ) const { return m_lap_count; }
    
    void reset( void ) { m_om = 0.0; m_count = 0; m_lap_count = 0; m_R1 = 0; m_R2 = 0; }
};



template<class det, class stoch>
struct observer_s
{
    ostream&        m_outfile;
    const double    m_sim_time;
    string          m_output;
    
    double m_lamb2;
    double m_D;
    size_t m_count;

    matrixd m_dxdt_d;
    matrixd m_dxdt_s;
    det m_odedet;
    stoch m_odestoch;
    
    observer_s( det odedet, stoch odestoch, matrixd x, ostream &out, const float &sim_time ) :
        m_dxdt_d(x) , m_dxdt_s(x) , m_odedet(odedet) , m_count( 0 ) ,
        m_odestoch(odestoch) , m_outfile( out ) ,
        m_sim_time( sim_time ) { }
    
    void set_params( double lamb2 , double D ) { m_lamb2 = lamb2; m_D = D; }
    
    void operator()( matrixd &x, double t )
    {
        m_odedet( x, m_dxdt_d, t );
        m_odestoch( x, m_dxdt_s, t );
        matrixd dxdt = m_dxdt_d + m_dxdt_s;

	if(m_count % 100 == 0 ){
            for(int i=0; i < x.size1(); ++i){
                m_outfile << m_D << '\t' << m_lamb2 << '\t' << t << '\t' << i
                << '\t' << x(i,1) << '\t' << x(i,0) << '\t' << dxdt(i,1) <<
                '\t' << dxdt(i,0) << '\n';
            }
        }
        ++m_count;

    }

    void reset( void ) { m_count = 0; }

};




#endif /* observer_h */

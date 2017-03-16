#ifndef _ESTIMATOR_HEADER_
#define _ESTIMATOR_HEADER_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <typeinfo>

#include "Particle.h"
#include "Material.h"
#include "Reaction.h"

// base estimator class
class estimator {
private:
    std::string estimator_name;
protected:
    unsigned long long nhist;
public:
    estimator( std::string label ) : estimator_name(label) {};
    ~estimator() {};
    
    virtual std::string name() final { return estimator_name; };
    
    virtual void score( particle*, double, double ) = 0;
    template< typename T >
    // particle, distance travelled, macro cap xs, cell name
    void score( particle*, double,  double,  T ) { assert(false); };
    void score( particle*, double,  double, std::shared_ptr< material > )  {};
    
    virtual void endHistory()       = 0;
    virtual void report( int )      = 0;
};

// derived class for simple estimators like current or scalar flux
// future estimators could get spectra or flux on a mesh
class single_valued_estimator : public estimator {
private:
    
protected:
    double tally_hist, tally_sum, tally_squared;
    double theMacroCapXs;
public:
    using estimator::score;
    
    single_valued_estimator(std::string label ) : estimator(label) {
        nhist              = 0;
        tally_hist         = 0.0;
        tally_sum          = 0.0;
        tally_squared      = 0.0;
    };
    ~single_valued_estimator() {};
    
    virtual void endHistory()    final {
        nhist++;
        tally_sum     += tally_hist;
        tally_squared += tally_hist * tally_hist;
        tally_hist     = 0.0;
    }
    
    virtual void score( particle*, double, double) = 0;
    virtual void report(int) = 0;
    
};

// current crossing a surface-----------------------------------------------------
class surface_current_estimator : public single_valued_estimator {
private:
    
public:
    surface_current_estimator( std::string label ) : single_valued_estimator(label) {};
    ~surface_current_estimator() {};
    
    void score( particle*, double, double );
    
    void report(int) {
        double mean = tally_sum / nhist;
        double var  = ( tally_squared / nhist - mean*mean ) / nhist;
        std::cout << name() << "   " << mean << "   " << std::sqrt( var ) / mean << std::endl;
    };
    
};

// volume-averaged scalar flux in a cell--------------------------------------------
class cell_pathLengthFlux_estimator : public single_valued_estimator {
private:
    double volume; //we don't use the volume in this hw
public:
    cell_pathLengthFlux_estimator( std::string label , double vol ):
    single_valued_estimator(label) , volume(vol) {};
    ~cell_pathLengthFlux_estimator() {};
    
    void score( particle* , double, double );
    
    void report(int numMove) {
        
        //for all cells that have estimator
        double mean = tally_sum / nhist;
        double var  = ( tally_squared / nhist - mean*mean ) / nhist;
        
        //output the results
        std::cout <<"   Capture rate : " << mean  << "   " << std::sqrt( var ) / mean << std::endl;
        std::cout <<"   FOM          : " << 1.0/(numMove * pow( (std::sqrt( var ) / mean ) ,2) )<< std::endl;
        //std::cout <<"Number of tracks: " << numMove << std::endl;
        
    };
};

//counting_estimator---------------------------------------------------------------
class counting_estimator : public estimator {
private:
    int count_hist;
    std::vector< double > tally;
public:
    counting_estimator( std::string label ) : estimator(label) { count_hist = 0; };
    ~counting_estimator() {};
    
    void score( particle*, double, double);
    void endHistory();
    void report(int);
};

#endif

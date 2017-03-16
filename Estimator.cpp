#include <cmath>
#include <iostream>

#include "Estimator.h"
#include "Material.h"
#include "Particle.h"

//surface_current_estimator------------------------------------------------------
void surface_current_estimator::score( particle* p, double null, double macroCapXs) {
    tally_hist += p->wgt();
}

//counting_estimator--------------------------------------------------------------
void counting_estimator::score( particle* p, double null, double macroCapXs) {
    count_hist++;
}

void counting_estimator::endHistory() {
    if ( tally.size() < count_hist + 1 ) { tally.resize( count_hist + 1, 0.0 ); }
    tally[ count_hist ] += 1.0;
    nhist++;
    count_hist = 0;
}

void counting_estimator::report(int) {
    std::cout << name() << std::endl;
    double s1 = 0.0, s2 = 0.0;
    for ( int i = 0 ; i < tally.size() ; i++ ) {
        double p = tally[i] / nhist;
        std::cout << i << " " << p << "   " << std::sqrt( p * ( 1.0 - p ) / nhist ) / p <<  std::endl;
        s1 += p * i;
        s2 += p * i * i;
    }
    std::cout << "   mean = " << s1 << std::endl;
    std::cout << "   var  = " << s2 - s1*s1 << std::endl;
}

//cell_pathLengthFlux_estimator---------------------------------------------------
void cell_pathLengthFlux_estimator::score( particle* p , double path_length, double macroCapXs){
    
    //for all cells that have path length flux estimator, get the capture rate
    tally_hist += p->wgt()* path_length * macroCapXs;
    
}


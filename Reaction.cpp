#include "Reaction.h"
#include "Particle.h"
#include "Distribution.h"

void  capture_reaction::sample( particle* p, std::stack< particle >* bank ) {
    // kill the particle and leave the bank unmodified
    p->kill();
}

void  scatter_reaction::sample( particle* p, std::stack< particle >* bank ) {
    // scatter the particle and leave the bank unmodified
    double mu0 = scatter_dist->sample();
    p->scatter( mu0 );
}

void  fission_reaction::sample( particle* p, std::stack< particle >* bank ) {
    
    int n = multiplicity_dist->sample();
    if ( n <= 0 ) {
        p->kill();
    }
    else {
        // bank all but last particle (skips if n = 1)
        for ( int i = 0 ; i < (n - 1) ; i++ ) {
            particle q( p->pos(), isotropic->sample() ); //created with weight = 1
            q.recordCell( p->cellPointer() );
            q.adjustWeight(p->wgt()); //but this new particle should have the same weight as its parent
            bank->push( q );
        }
        // set working particle to last one
        particle q( p->pos(), isotropic->sample() );
        q.recordCell( p->cellPointer() );
        q.adjustWeight(p->wgt());
        *p = q;
    }
}

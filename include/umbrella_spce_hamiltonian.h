/* 
 * File:   umbrella_spce_hamiltonian.h
 * Author: BENSON J MA
 *
 * Created on February 1, 2012, 3:56 PM
 */

#ifndef UMBRELLA_SPCE_HAMILTONIAN_H
#define	UMBRELLA_SPCE_HAMILTONIAN_H

class UmbrellaSPCEHamiltonian : public SPCEHamiltonian {
private:
    double &WINDOW_LOWER_BOUND, &WINDOW_UPPER_BOUND;
    bool OVER_THE_WINDOW;
    
    virtual double energy_between_two_ions(int i, int j, WHICH_TYPE typ);

public:
    UmbrellaSPCEHamiltonian(WaterSystem &sys);
    ~UmbrellaSPCEHamiltonian();
    virtual void initialize_calculations();
    virtual double total_energy_difference();
    virtual void undo_calculations();
};

#endif	/* UMBRELLA_SPCE_HAMILTONIAN_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //
#define POW2(X) ((X)*(X))
#define POW4(X) ((X)*(X)*(X)*(X))
#ifndef XTRACK_PulsedLaser_H
#define XTRACK_PulsedLaser_H

/*gpufun*/
void PulsedLaser_track_local_particle(PulsedLaserData el, LocalParticle* part0){

    //The algorithm is partially from https://anaconda.org/petrenko/psi_beam_vs_laser

    double nx  = PulsedLaserData_get_laser_direction_nx(el);
    double ny  = PulsedLaserData_get_laser_direction_ny(el);
    double nz  = PulsedLaserData_get_laser_direction_nz(el);
    
    double laser_x = PulsedLaserData_get_laser_x(el);
    double laser_y = PulsedLaserData_get_laser_y(el);
    double laser_z = PulsedLaserData_get_laser_z(el);
    double w0 = PulsedLaserData_get_laser_waist_radius(el);

    double laser_energy = PulsedLaserData_get_laser_energy(el);
    double laser_waist_shift = PulsedLaserData_get_laser_waist_shift(el);

    double laser_sigma_t = PulsedLaserData_get_laser_duration_sigma(el);
    //double laser_sigma_w = 1/laser_sigma_t; // rad/sec -- assuming Fourier-limited pulse

    double laser_wavelength = PulsedLaserData_get_laser_wavelength(el); // Hz
    
    double ion_excited_lifetime = PulsedLaserData_get_ion_excited_lifetime(el); // sec
    double ion_excitation_energy = PulsedLaserData_get_ion_excitation_energy(el); // eV
       
    // constants for the Map_of_Excitation vs OmegaRabi and Detuning:
    int64_t N_OmegaRabiTau_values = PulsedLaserData_get_N_OmegaRabiTau_values(el);
    int64_t N_DeltaDetuningTau_values = PulsedLaserData_get_N_DeltaDetuningTau_values(el);
    double  OmegaRabiTau_max = PulsedLaserData_get_OmegaRabiTau_max(el);
    double  DeltaDetuningTau_max = PulsedLaserData_get_DeltaDetuningTau_max(el);
    double  dOmegaRabiTau = OmegaRabiTau_max/(N_OmegaRabiTau_values-1.0);
    double  dDeltaDetuningTau = DeltaDetuningTau_max/(N_DeltaDetuningTau_values-1.0);
    
    double laser_Rayleigh_length = PI*w0*w0/laser_wavelength;
    
    double p0c = LocalParticle_get_p0c(part0); // eV
    double m0  = LocalParticle_get_mass0(part0); // eV/c^2
    double hbar = 1.054571817e-34; // J*sec
    
    double OmegaTransition = ion_excitation_energy*QELEM/hbar; // rad/sec
    
    // Maximum laser intensity (at the focal point)
    double I0 = sqrt(2/PI)*(laser_energy/laser_sigma_t)/(PI*w0*w0); // W/m^2
    
    //start_per_particle_block (part0->part)
    
        double state = LocalParticle_get_state(part);
        double delta = LocalParticle_get_delta(part);
        double z     = LocalParticle_get_zeta(part);
        double x     = LocalParticle_get_x(part);
        double y     = LocalParticle_get_y(part);
        double px    = LocalParticle_get_px(part);
        double py    = LocalParticle_get_py(part);
    
        double pc = p0c*(1.0+delta); // eV
        double gamma = sqrt(1.0 + pc*pc/(m0*m0));
        double beta  = sqrt(1.0 - 1.0/(gamma*gamma));
        double beta_x  = px*p0c/m0/gamma;
        double beta_y  = py*p0c/m0/gamma;
        double beta_z  = sqrt(beta*beta - beta_x*beta_x -beta_y*beta_y);

        double vx  = C_LIGHT*beta_x; // m/sec
        double vy  = C_LIGHT*beta_y; // m/sec
        double vz  = C_LIGHT*beta_z; // m/sec
            
        // Collision of ion with the laser pulse:
        // The position of the laser beam center is rl=rl0+ct*n. We can find the moment
        // when a particle with a position r=r0+vt collides with the laser as the moment
        // when r−rl is perpendicular to n. Then (r−rl,n)=0, which yields the equation
        // (r0,n)+(v,n)t−(rl0,n)−ct(n,n)=0. Hence
        // tcol=(r0−rl0,n)/[c−(v,n)]

        double tcol = ( (x-laser_x)*nx + (y-laser_y)*ny + (z-laser_z)*nz ) / (C_LIGHT - (vx*nx+vy*ny+vz*nz)); // sec
        
	    double xcol = x + vx*tcol; // m
	    double ycol = y + vy*tcol; // m
	    double zcol = z + vz*tcol; // m
        
        // r^2 to the laser center = |r-rl| at the moment tcol:
        double r2 = (\
                POW2(xcol - (laser_x+C_LIGHT*nx*tcol)) + \
                POW2(ycol - (laser_y+C_LIGHT*ny*tcol)) + \
                POW2(zcol - (laser_z+C_LIGHT*nz*tcol)) \
             ); // m

        double Z_to_laser_focus = laser_waist_shift - tcol*C_LIGHT; // m
        
        // Laser beam size at the point of collision:
        double w = w0*sqrt(1.0+POW2(Z_to_laser_focus/laser_Rayleigh_length));
        // Max. laser intensity experienced by the ion (in the ion rest frame):
        double I = 4.0*gamma*gamma * I0*(w0/w)*(w0/w)*exp(-2.0*r2/(w*w)); // W/m^2
        double OmegaRabi = \
            (hbar*C_LIGHT/(ion_excitation_energy*QELEM)) * \
            sqrt(I*2*PI/(ion_excitation_energy*QELEM*ion_excited_lifetime)); // rad/sec
        
        double OmegaRabiTau = OmegaRabi*laser_sigma_t/(2.0*gamma); // in the ion rest frame
        
        // Detuning from the ion transition resonance in the ion rest frame:
        double cos_theta = -(nx*vx + ny*vy + nz*vz)/(beta*C_LIGHT);
        double laser_omega_ion_frame = (2.0*PI*C_LIGHT/laser_wavelength)*(1.0+beta*cos_theta)*gamma;
        double DeltaDetuningTau = fabs(
            (OmegaTransition - laser_omega_ion_frame)*laser_sigma_t/(2.0*gamma)
        );
        
        //only apply laser cooling to particles that have not been lost. In Xsuite, this means positive state
        if (state > 0)
            {    
 
        if (DeltaDetuningTau < DeltaDetuningTau_max && OmegaRabiTau > OmegaRabiTau_max)
            {
                //printf("high field case: DeltaDetuningTau = %f\n", DeltaDetuningTau);
                // printf("high field case: OmegaRabiTau = %e\n", OmegaRabiTau);
                // In case of a very high laser field:
                LocalParticle_set_state(part, 2); // Excited particle
                double rnd = (float)rand()/(float)(RAND_MAX);
                LocalParticle_add_to_energy(part,-ion_excitation_energy*2.0*rnd*2.0*gamma, 0); // eV
            }
    
        else if (DeltaDetuningTau < DeltaDetuningTau_max &&
            OmegaRabiTau > dOmegaRabiTau/10.0)
            {
            // N_OmegaRabiTau_values  N_DeltaDetuningTau_values
            //   OmegaRabiTau_max       DeltaDetuningTau_max
            int64_t row = (int)floor(OmegaRabiTau/dOmegaRabiTau);
            int64_t col = (int)floor(DeltaDetuningTau/dDeltaDetuningTau);
            int64_t idx = row*N_DeltaDetuningTau_values + col;
            
            double excitation_probability = PulsedLaserData_get_Map_of_Excitation(el, idx);
                                        
            double rnd = (float)rand()/(float)(RAND_MAX); //between 0 and 1
            if ( rnd < excitation_probability )
                {
                LocalParticle_set_state(part, 2); // Excited particle
                // photon recoil (from emitted photon!):
                double rnd = (float)rand()/(float)(RAND_MAX);

                // If particle is excited, reduce its energy by, on average, the excitation energy with Lorentz boost
                // 2.0*rnd ensures that the average energy lost is the excitation energy
                // 2.0*gamma is the Lorentz boost
                LocalParticle_add_to_energy(part,-ion_excitation_energy*2.0*rnd*2.0*gamma, 0); // eV

                // //Also do transverse emission. This doesnt have a Lorentz boost
                // double rnd_x = ((double)rand() / (RAND_MAX)) * 2.0 - 1.0;//between -1 and +1
                // LocalParticle_add_to_px(part,rnd_x*ion_excitation_energy/p0c);
                // double rnd_y = ((double)rand() / (RAND_MAX)) * 2.0 - 1.0;//between -1 and +1
                // LocalParticle_add_to_py(part,rnd_y*ion_excitation_energy/p0c);
                }	
                else
                {
                LocalParticle_set_state(part, 1); // Still particle
                }

            }
            }
            
	//end_per_particle_block
    
}

#endif

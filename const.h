#ifndef CONST_H
#define CONST_H

/* Avogadro's Number [/mol] */
const double AVOGADRO     = 6.022045e23;

/* Electronic Charge [Coulomb] */
const double COULOMB      = 1.602189e-19;

/* Plank Constant [MeV*cm] */
const double hc           = 197.327053;

/* Fine Structure Constant */
const double ALPHA        = 1/137.036;

/* Particle Mass */
const double MASS_ATOMIC_UNIT = 931.49432;
const double MASS_PROTON      = 938.2723128; 
const double MASS_NEUTRON     = 939.5656; 
const double MASS_PIPLUS      = 139.56995; 
const double MASS_PI0         = 134.9764;       
const double MASS_ELECTRON    = 0.51099906;       
const double MASS_PHOTON      = 0.0;       

/* Proton RMS radius */
const double r_proton     = 0.862 ; // fm 
//[ref.] G.G.Simon et al. NPA333(1980)380

/* PiPlus Decay Const */
const double PIPLUS_DECAY = 26.03e-9 * 3.0e8 ; // unit = [m]

/* Unit Conversion Factors */
const double rad2deg      = 180/M_PI;
const double deg2rad      = M_PI/180;
const double r2d          = rad2deg;
const double d2r          = deg2rad;
const double cm2mcb       = 1e30;       // [cm^2] -> [mcb] 
const double mcb2cm       = 1e-30;      // [mcb]  -> [cm^2]
const double fm2mb        = 10;         // [fm^2] -> [mb]
const double fm2mcb       = 1e4;        // [fm^2] -> [mcb]
const double mcb2mb       = 1e-3;       // [mcb]  -> [mb]
const double M2G          = 1e-3;       // [MeV]  -> [GeV]
const double G2M          = 1e3;        // [GeV]  -> [MeV]

#endif /* CONST_H */


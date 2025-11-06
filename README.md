# Kinema

## How to install

```
git clone https://github.com/Sangbaek/Kinema.git
cd Kinema
make
```

## How to use


If we want to reproduce the Kin-1 of Spring 2026 plan for VCS experiment,
which is the following.
$$
Q^2 = .05 \mathrm{GeV}^2/c^2\\
p_{e'} = 1055.28 \mathrm{MeV}/c\\
\theta_{e'} = 10.48^{\circ}\\
p_{p'} = 493.93 \mathrm{MeV}/c\\
\theta_{p'} = 56.01^{\circ}\\	
\theta_{gg} = 110^{\circ}\\	
\phi = 0
$$

Then, first, run the kinema with only electron side, i.e.,

```
./kinema -P -e 1419 -Q 0.05 -W 1230 -G
```

The last flag `-G` is for the VCS. For the N-Delta, simply do not use `-G`.

This will print out the angle of virtual photon to be 26.7275 degrees, like the following.

```
	 Incident Energy  = 1419 [MeV]
	 Scattering Angle = 10.4843 [deg.]
	 Final Energy     = 1055.28 [MeV]

 w        = 363.724 [MeV]
 W        = 1230 [MeV]
 q        = 426.961 [MeV/c]     2.16372[fm^-1]
 theta_q  = 26.7275 [deg.]
 q^2(4vec)= -0.05 [GeV/c]^2     -1.28409 [fm^-2]
 Q^2      = -0.05 [GeV/c]^2     1.28409 [fm^-2]
 ```

Then, consider the SHMS range 13.6 to 56.3 degrees.

When $$\theta_{p'} = 56.01^{\circ}$$, $$theta_{pq}= 56.01-26.7275^{\circ}\sim29.3^{\circ}$$.

Then,

`./kinema -P -e 1419 -Q 0.05 -W 1230 -T 29.3 -G`

will print out 

```
H    938.27  122.02  493.83   29.30000   34.60  257.13   70.031     4.70
H      0.00  241.70  241.70  -90.87616  257.13  257.13  109.969     0.88
```

and $$p_{p'} = 493.83$$ MeV/$$c$$ well agrees with the kinematics setting within the numerical uncertainties.

Take 109.969 for the $\theta_{gg}.
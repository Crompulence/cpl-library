/*

Overview of the CPL library LAMMPS interface
=============================================

LAMMPS works through a series of user defined fixes. 
These are objects the user can define by inheretiting from a base class called fix with a pre-definied structure.
These are automatically registered with LAMMPS and instantiated (created when LAMMPS runs).
The use can then defined a range of "hook"* functions with names like post_setup or pre_force, which can be used to inject code where you want (after the setup or before the force calcultation respectivly) in the overall LAMMPS solver. 
The addition code is included in LAMMPS by writing user add on packages which can be included when LAMMPS is complied. 
Provided LAMMPS has been built with the extra package, they can then be used, often by keywords to switch these on from the user input. 

For CPL library, the package is USER_CPL and the code inside includes fix_cpl_init, fix_cpl_force and CPLSocketLAMMPS.
The socket is the bulk of the code which handles getting information from LAMMPS, fix_cpl_force applies the constraint force and fix_cpl_init is the top level function to handle communication and call the CPLSocketLAMMPS/fix_cpl_force routines.


* The concept of a hook here is a programming term for a pre-defined place in a software package when you can stick in your own code. LAMMPS is designed to be almost entirly built up of various user packages so provides hooks to get at almost every possible point in the solver.


Section on fix_cpl_init
============================
For CPL library, the fix we use is fix_cpl_init, which handles setup of the CPL library, using code in the "setup" hook to set up a coupled topology, extracting all required parameters from LAMMPS and sending these via MPI (to coupled code).
The required averages are performed in the post_force hook, then information is exchanged between the codes and the constraint forces applied.

The coupling fix, fix_cpl_init, is turned by by adding a command of this form to LAMMPS input:

    fix ID group-ID cpl/init region all forcetype X sendtype Y

The first two arguments are standard LAMMPS:
ID = user-assigned name for the fix
group-ID = ID of the group of atoms to apply the fix to
The next part, cpl/init, is the "style", i.e. the fixes' name. 
"region all" specifies that the fix is applied to the entire box in space. 
Note that in the granular case, this is clearly consistent (as CFD overlaps the whole domain).
For domain decompositional coupling, region should still be all, as we specify the extents of coupling through the COUPLER.in file. 
The remaining words are the args which are currently as follows:

forcetype X -- This allows you to specify which constraint force system to use. 
    The options include 
        1) test -- A simple test for debugging (sends cell indices)
        2) Flekkoy -- stress based coupling
        3) Velocity -- applies (U_particle - U_CFD) which is the most important part of the constraint of Nie, Chen, E and Robbins (2004)
        4) Drag -- base drag class for granular type drag forces
        5) Di_Felice -- an example granular drag class which applies the Di Felice drag correlation (untested!!)
        6+) This has been design so the user can add anything they want here easily. This process is described below. 

sendtype Y -- Specifiy which data is sent to the CFD solver. 
    This Y can be 
        1) velocity --4 values including 3 velocity sums and number of particles
        2) gran -- which sends voidratio and force 
        3) granfull - is designed for SediFOAM and sends velocity, force, sum of drag coefficients and voidratio.

Optional: bndryavg Z -- the location of the bottom cell in domain-decompositional coupling, Z can be either below, above or midplane.


Section on sendtypes
====================

Section on forcetype and fix_cpl_force 
======================================

fix_cpl_force
--------------------

The force class have a pre-determined interface, set by an abstract base class, to allow them to always be called from fix_cpl_force in the same way.
Think of an interface like a contract, you guarentee that you will always take the same arguments in and return the same things.
In this case, we have an interface which takes in position, velocity, acceleration, mass, radius (sigma) and interaction strength (epsilon) and works out summed up arrays needed for the particular force type (in the pre_force function) and returns the actual force (in the get_force function).
As each type of coupling force you could want to apply always has the same form:

1) pre_force : Get some stuff from the particle system (e.g. add up current void fraction, get latest velocity from CFD solver) 
2) get_force : Use the pre_force information, with the CFD solver information, to get the force on each particle. 

These means we can make use of c++ polymorphism, where we choose which type of force we want based on the forcetype input argment.
The required force type object is instantiated using a "factory" which takes the user input and returns fxyz, a pointer to the right object: 
*/

//Force factory
if (fxyzType.compare("Flekkoy") == 0) {
    fxyz.reset(new CPLForceFlekkoy(9, cfdBuf->shape(1), 
                                      cfdBuf->shape(2), 
                                      cfdBuf->shape(3)));
    fxyz->calc_preforce = true;
} else if (fxyzType.compare("test") == 0) {
    fxyz.reset(new CPLForceTest(3, cfdBuf->shape(1), 
                                   cfdBuf->shape(2), 
                                   cfdBuf->shape(3)));
} else if (fxyzType.compare("Velocity") == 0) {
    fxyz.reset(new CPLForceVelocity(3, cfdBuf->shape(1), 
                                       cfdBuf->shape(2), 
                                       cfdBuf->shape(3)));
} else if (fxyzType.compare("Drag") == 0) {
    fxyz.reset(new CPLForceDrag(9, cfdBuf->shape(1), 
                                   cfdBuf->shape(2), 
                                   cfdBuf->shape(3),
                                   false)); 
    fxyz->calc_preforce = true;

} else if (fxyzType.compare("Di_Felice") == 0) {
    fxyz.reset(new CPLForceGranular(9, cfdBuf->shape(1), 
                                       cfdBuf->shape(2), 
                                       cfdBuf->shape(3))); 

} else {
    std::string cmd("CPLForce type ");
    cmd += fxyzType + " not defined";
    throw std::runtime_error(cmd);
}

//Once we have this pointer, it can then be used by looping over all particles in LAMMPS (nlocal here)

    //Pre-force calculation, get quantities from discrete system needed to apply force
	for (int i = 0; i < nlocal; ++i)
	{
   		if (mask[i] & groupbit)
    	{
	        //Get local molecule data
	        mi = rmass[i];
	        radi = radius[i];
	        for (int n=0; n<3; n++){
	            xi[n]=x[i][n]; 
	            vi[n]=v[i][n]; 
	            ai[n]=f[i][n];
	        }

	        // Sum all the weights for each cell.
	        fxyz->pre_force(xi, vi, ai, mi, radi, pot);

    	}
    }


//This pre-force has added all the needed  things for the force you want, so we can simply get the force now:

   // Calculate force and apply
    for (int i = 0; i < nlocal; ++i)
    {
        if (mask[i] & groupbit)
        {

            //Get local molecule data
            mi = rmass[i];
            radi = radius[i];
            for (int n=0; n<3; n++){
                xi[n]=x[i][n]; 
                vi[n]=v[i][n]; 
                ai[n]=f[i][n];
            }

            //Get force from object
            fi = fxyz->get_force(xi, vi, ai, mi, radi, pot);

        }
    }

/*
Section on ndArray
=====================

Coupler library is build on the concept of exchanging 4D uniform grids of data. As a result, the first requirement is to setup a 4D array of data based on the number of cells on a given processes. 

The most useful form of CPL_ndArray is a 4D array, that is a dimensionality (i.e. 3 for a vector, 9 for a stress tensor, 100 for a timeseries) and then the three physical dimensions in space, ordered as nd, icell, jcell, kcell. CPL_ndArray has a range of helper functions which mirror numpy or Fortran style arrays.


Section on CPL_field
=====================

CPL_field extends the ndarray to build up a physical field. This is an encapsulated relationship, the field contains an array of field data, extending it to include domain extents, functions to interpolate or average values in space. So far, there is no inheritance used with field as we can get the required flexitbility from dynamic dispatch (e.g. add_to_array can include a radius, in which case the fraction inside a cell will be calculated).



Section on CPL_force
=====================

CPL_force combines multiple field classes, including the retrieved values from coupling and accumulated MD/DEM results before returning the required force. Heavy use is made of inheritance here, with a common interface of a pre_force and a get_force function. The different types of coupling force are then created using a factory method and extension of this is as simple as inhereting from the most relevant force class and adapting to your needs.
A rough schematic of the inheretence diagram is included below (not names have been shortened): 

                         CPL_vel
Abstract Base Class     / 
   |          _________/__CPL_test
   |         /         \ 
   |        /           CPL_flekkoy
   v       /
CPL_force /                       __CPL_Di_Felice
          \__ CPL_drag__ CPL_gran/
                      \          \__CPL_Tang
                       \            
                        CPL_dragtest

Tutorial on Drag Forces
=====================

The process of designing an appropriate drag force is given as an example here. 
Note we have cut lots of boilerplate code so this will not work as expected.
The user should refer to the examples in CPL_force.cpp 

We will outline the process of designing the Ergun drag force,

F = 150.0*e*nu*rho/(pow(e*D, 2.0)) + 1.75*rho()*Ur/(e*D)

where we need to work out both porosity e and the fluid velocity Ur = U_MD-U_CFD 
from the particle simulation. We use an existing class as a basis, 
choose the one that is closest to your current case: 

1) Choose a CPL_force type to inheret, probably the CPLForceGranular in this case.
To define an inheretemce frp, CPLForceGranular class, we use 
class CPLForceErgun : public CPLForceGranular 

*/

class CPLForceErgun : public CPLForceGranular {

public:

    //Constructors
    CPLForceErgun(CPL::ndArray<double> field);
    CPLForceErgun(int nd, int icell, int jcell, int kcell);

    //Pre force collection and get force calculation
    // position, velocity, acceleration, mass, radius, interaction
    void pre_force(double r[], double v[], double a[], 
                   double m, double s, double e);
    std::vector<double> get_force(double r[], double v[], double a[], 
                                  double m, double s, double e);

    bool calc_preforce = true;

    CPL::ndArray<double> eSums;

private:

    void initialisesums(CPL::ndArray<double> f);
    void resetsums();

};

/*
2) Now we have specified the header, we then need to override any routines which are unique to CPLForceErgun in 
the CPL_force.cpp file.
But first, we need to define the constructors of our new class, these take in either the size of the
grid array (local to a processes) or an existing array.
These functions essentially creates the field which is used to get data from the CFD (in the parent).
This is referred internally in functions of CPLForceErgun as fieldptr and the array data can be obtained.
The child automatically calls the parent constructor and then we explicitly call initialisesums to setup 
various fields which will be populated pre_force.
The logic here is the CPL array expected from the CFD code can be used and all other fields created 
with a consistent size.
The density and viscosity (and anything else) can be added in here.
*/

//Constructor using cells
CPLForceErgun::CPLForceErgun(int nd, int icells, 
                             int jcells, int kcells, 
                             double rho_, double mu_) : CPLForceGranular(nd, icells, jcells, kcells){
    rho = rho_;
    mu = mu_;
    initialisesums(fieldptr->get_array());
}

//Constructor of datatype
CPLForceErgun::CPLForceErgun(CPL::ndArray<double> arrayin, 
                             double rho_, double mu_) : CPLForceDrag(arrayin){
    rho = rho_;
    mu = mu_;
    initialisesums(arrayin);
}

/*
3) Next, we need  to develop pre force and get_force for the Ergun equation.
As Ergun needs porosity, this must be collected pre-force.
The field classes provide an elegent way to get these (currently not used in most of cpl_force as this is new).
We simply add to an array based on molecular position r, the value we want to, here a count for 
no. of molecules for nSums (1) and the sphere volume for eSums.
*/


//Pre force collection of sums (should this come from LAMMPS fix chunk/atom bin/3d)
void CPLForceErgun::pre_force(double r[], double v[], double a[], double m, double s, double e) {
    
    nSums.add_to_array(r, 1.0);
    double vol = (4./3.)*M_PI*pow(s,3);
    eSums.add_to_array(r, vol);
}

/*
If you wanted to use partial overlap calculation (fraction of sphere in box) then simply add the
 radius of the shere in as a second argument,
    eSums.add_to_array(r, s, vol);

4) We need to work out the value of force defined as a 3 vector.
*/

//Get force using sums collected in pre force
std::vector<double> CPLForceErgun::get_force(double r[], double v[], double a[], double m, double s, double e) {

    std::vector<double> f(3), Ui(3);

    //Porosity e is cell volume - sum in volume
    double e = 1.0 - eSums.get_array_value(r)/eSums.dV;
    double D = 2.0*s;

    //Should use std::vector<double> Ui(3) = field.interpolate(r);
    CPL::ndArray<int> indices = {0,1,2};
    Ui = fieldptr->get_array_value(indices, r);
    for (int i=0; i<3; i++){
        double Ur = Ui[i]-v[i];
        f[i] = 150.0*e*nu*rho/(pow(e*D, 2.0)) + 1.75*rho*Ur/(e*D);
        FSums(i, cell[0], cell[1], cell[2]) += f[i];
    }
    return f;

}

/*
5) Finally, we need to add out new force type into the fix_cpl_force
*/

//Force factory
if (fxyzType.compare("Flekkoy") == 0) {
    fxyz.reset(new CPLForceFlekkoy(9, cfdBuf->shape(1), 
                                      cfdBuf->shape(2), 
                                      cfdBuf->shape(3)));
...
...
...
} else if (fxyzType.compare("Ergun") == 0) {
    fxyz.reset(new CPLForceErgun(9, cfdBuf->shape(1), 
                                   cfdBuf->shape(2), 
                                   cfdBuf->shape(3),
                                   rho, mu));


/*
and we can turn on by adding Ergun to the input line,

    fix ID group-ID cpl/init region all forcetype Ergun sendtype Granfull

*/

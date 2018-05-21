# Basic functionality test

Basic tests which ensure all the CPL library utilities for defining ndArrays, creating fields 
(with overlap and interpolation) and applying drag forces work as expected.
This is important as the major complexity of coupled simulation is in these details, which
have been generalised and abstracted away so they can be tested in isolation.

This uses google tests.

##The tests are as follows

### Tests for CPL::ndArray
 - Setup and array size 
 - Set elements
 
### Tests for CPL::fields

 - Test setting and changing of fields
 - Test field constructor with optional inputs and setting with names
 - More examples of field setting and use of pointers to array
 - Check overlap calculation (uses [overlap library](https://github.com/Crompulence/cpl-library/tree/master/src/utils/overlap) and [Eigen](http://eigen.tuxfamily.org))
 - Check getters for field properties

### Test for CPL::forces
 - Test constructors for base class
 - Check getting and setting internal field (whihc is recv from coupling)
 - Get cell method to bin molecules
 - Get force method
 - Test constructor for internal fields and name retrival
 - Test for velocity class with pre force and get force
 
### Test for CPL::forces Flekkoy
 - Test for Flekkoy style force
 - Check gweight calculation 
 - Check Flekkoy pre foce 
 - Test with varying domain sizes
 
### Test for CPL::forces Drag
 - Test constructors for drag class
 - Test constructor with optional inputs and setting with names
 - Test summing volume and forces in cells
 - Test summing and volumes for random and extensive generated values
 - Test overlap calculation used in force
 - Test Drag force with specified drag coefficient
 
### Test for CPL::forces Granular
 - Test granular constructor
 - Test force calculation
 
### Test for interpolation

 - Test mapping of cell to node
 - Test interpolation within cells
 - Test force calculation
  

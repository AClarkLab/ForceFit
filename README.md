# ForceFit

The *ForceFit* program package has been developed for fitting classical force field parameters based upon a force matching algorithm to quantum mechanical gradients of configurations that span the potential energy surface of the system. The program, which runs under UNIX and is written in C++, is an easy-to-use, nonproprietary platform that enables gradient fitting of a wide variety of functional force field forms to quantum mechanical information obtained from Gaussian and NWChem. All aspects of the fitting process are run from a graphical user interface, from the parsing of quantum mechanical data, assembling of a potential energy surface database, setting the force field, and variables to be optimized, choosing a molecular mechanics code for comparison to the reference data, and finally, the initiation of a least squares minimization algorithm. Furthermore, the code is based on a modular templated code design that enables the facile addition of new functionality to the program.

*ForceFit* is offered under a LGPL license, wherein the user may not redistribute the code, but modifications may be made and if sent back to Prof. Clark, will be incorporated into the next version of the code.

Any results obtained with *ForceFit* should refer to the following publication:

Waldher, B.; Kuta, J.; Chen, S.; Henson, N.; Clark, A. E. "*ForceFit*: A Code to Fit Classical Force Fields to Quantum Mechanical Potential Energy Surfaces", J. Comp. Chem. **2010**, *31*, 2307-2316.

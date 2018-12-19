# About
ArgusInterface is a thin layer between cloth dynamics codes and  so-bogus
Coulomb friction solver.
It also provide a way to export friction problems, and a corresponding
loader application for offline testing of the solver algorithms.

# Dependencies

Building this software requires CMake and a C++98 compiler, and has
only been tested on GNU/Linux with recent versions of `g++`.

Other dependencies are:

  - [Eigen][1] and [So-bogus][2], which are already included in this archive
  - [boost\_serialization][3]

So-bogus is provided as a git submodule of argus' repository, so no further installation should be required. (Maybe `git submodule update`).

# Compiling 
Standard CMake procedure:

	> mkdir build
	> cd build
	> cmake ../interface
	> make
	> [make install]

Install path may be specified with the `cmake -DCMAKE_INSTALL_PREFIX=[...]` switch.

# Algorithms

ArgusInterface currently provides 4 algorithms for solving the friction problems:

 - NodalContact corresponds to the algorithm described in our SCA 2015 poster, 
   that is, a primal SOCQP minimization inside a Cadoux fixed-point loop. 
   This only works for mesh-cloth contacts that coincide with the cloth vertices.
 - AlternatingMinimization correspond to a variant of ADMM implemented in bogus as DualAMA:
   Minimization of (.5 z M^(-1) z ) under the constraint that r is in the friction cone and z = H'r.
   Its main advantage is to require only multiplications by the stiffness matrix M, 
   while being able to handle arbitrary contact configurations. 
   However, precise tuning of the step sizes is necessary to obtain satifying covergence. 
   There is probably a lot of room for improvement.
 - DualProjectedGradient corresponds to a dual SOCQP minimization inside a Cadoux fixed-point loop.
   The dual (Delassus) operator W = H M^(-1) H' is not explicitly assembled, 
   which means that a linear solve is performed at each iteration of the projected gradient
   (or projected gradient descent) algorithm. 
   Several variants are implemented, e.g. with Nesterov acceleration. 
 - DualGaussSeidel corresponds to a standard Gauss-Seidel algorithm on the dual formulation,
   as e.g. in our 2011 Sig Asia paper. The (dense) W matrix is computed, which is severely expensive.


  [1]: http://eigen.tuxfamily.org     "Eigen, template library for linear algebra"
  [2]: http://gdaviet.fr/code/bogus   "So-bogus, Coulomb friction solver"
  [3]: http://www.boost.org/doc/libs/release/libs/serialization/ "Boost serialization library"

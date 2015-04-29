! E. Higham
! ME4 Final Year Project
! Pre-Processor Variable Declaration

! ================================================================================================ !
!                                   Physical Mesh Structure                                        !
! ================================================================================================ !

! Uncomment if you DO NOT wish to use grid refinement
#define xGridRefinement = 1
#define yGridRefinement = 1

! ================================================================================================ !
!                                  Poisson's Equation Solver                                       !
! ================================================================================================ !

! Uncomment if you want a 2nd order accurate Poisson's equation solver (remaining code will still be order 4)
! --> Note: using a 2nd order accurate solver will be computationally cheaper, but if refinement is used, the 
! grid stretching will bring about larger error than if 4th order is used.
!#define Order2 = 1

! Comment if you DO NOT wish to use MultiGrid methods of solution
#define MultiGrid = 1

! IF MultiGrid IS ENABLED, Uncomment the number of grid levels that are required (DEFAULT: 2 Grids with MultiGrid)
! --> Note, nGrids must be uncommented successively
#define nGrids3 = 1
#define nGrids4 = 1
#define nGrids5 = 1
!#define nGrids6 = 1

! ================================================================================================ !
!                                  			Simulation Type                                        !
! ================================================================================================ !

!#define SimulationVars VortexRebound
#define SimulationVars CavityFlow
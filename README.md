    !###################################################################################################################
    ! HDG (SEM & DG) Hybrod Solver with a regular mesh for the 1D elastic wave equation 
    ! NOTE: This program is yet to be tested against analytical solutions.
    ! 
    ! Language: Fortran 90, with parralel impelementation using OpenMP API
    ! 
    ! Sources:  Igel 2017, Leveque 2002, Hesthaven & Warburton 2008
    ! 
    ! The code used for arbitrary GLL points and weights was created by M.I.T departement of engineering. Link is hereafter
    ! https://geodynamics.org/cig/doxygen/release/specfem3d/    | file name: gll_library.f90
    ! 
    ! This is part of the Numerical Modelling Workshop.
    ! 
    ! Supervisor: Pr. Emmanual Chaljub
    ! Author    : Mus Benziane
    !
    ! Input file: [Example] [Last BC options are yet to be impelented...]
    ! testing              ! Model name prefix (model names: prefix_vp, prefix_rho)
    ! 6                    ! Polynomial order
    ! 200                  ! Number of elements SEM
    ! 200                  ! Number of elements DG
    ! 30.                  ! Element size
    ! 20.                  ! Wavelet's peak frequency
    ! 0.0002               ! Time step
    ! 20000                ! Number of time steps
    ! 150                  ! Source location (element number)
    ! 4                    ! Source location (gll point)
    ! 50                   ! Snapshot interval
    ! 1                    ! [1/2/3/4] 1: Free surface, 2: Rigid wall, 3: Periodic, 4: Sponge 
    ! 0                    ! Boundary condition only on the left side [not for periodic BC | on the right side for absorbing]
    ! 20                   ! Sponge layer width in gll points
    ! .65                  ! Att constant for sponge layer
    ! 0                    ! Use intial conditions instead of point source, and produce analytical solution
    ! 100                  ! Width of Gaussian used as intial condition


    ! -> Model files in C-Style binary floats [doubles]: Vs, Rho files are needed.
    !                                                  : For simple models, use create1Dmodel_files.f90
    ! 
    ! -> Outputs are created in OUTPUT/ if OUTPUT/ is not created by the user, the program will not handle it.
    !    Output files in OUTPUT directory:
    !
    !####################################################################################################################
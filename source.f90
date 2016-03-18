! Charles McEachern

! Fall 2015

! Note: This document wraps at column 80. 

! Note: This code uses camelCase for legibility. Fortran is caps-insensitive. 

! #############################################################################
! ################################################################### Main Loop
! #############################################################################

! To encourage the compiler to optimize it as much as possible, the main loop
! does not bring in any other modules. All arrays and helper functions used
! within the loop module, are defined within the loop module. On the other 
! hand, the loop module is used by several other modules -- the fields module
! initializes and outputs the field arrays, for example. 

module loop
  implicit none

  ! ===========================================================================
  ! ======================================================== Physical Constants
  ! ===========================================================================

  ! Most of these constants aren't used in the main loop, but it seems best to
  ! have them all in the same place. 
  double precision, parameter :: pi = 3.14159265358979
  ! Speed of light squared in (Mm/s)**2
  double precision, parameter :: cc = 299.792458**2
  ! Magnetic and electric constants in nH/m and mF/m. 
  double precision, parameter :: mu0 = 1256.63706
  double precision, parameter :: eps0 = 1/(cc*mu0)
  ! Proton and electron masses in g, and elementary charge in MC
  double precision, parameter :: mp = 1.67262e-24
  double precision, parameter :: me = 9.10938e-28
  double precision, parameter :: qe = 1.60218e-25
  ! Radius from Earth's center to the surface and to the ionosphere, in Mm. 
  double precision, parameter :: RE = 6.378388
  double precision            :: RI

  ! ===========================================================================
  ! ============================================================= Grid Indexing
  ! ===========================================================================

  ! Number of field lines and number of points per field line respectively. 
  integer                            :: n1, n3
  ! First and last indeces for easy handling of the edges. 
  integer, dimension(0:1)            :: hh, ii, kk, zz
  ! Azimuthal modenumber. 
  integer                            :: azm

  ! ===========================================================================
  ! ====================================================================== Time
  ! ===========================================================================

  ! Time step, time between outputs, maximum time, current time. Several of 
  ! these are user parameters. The time step is not -- it's computed based on 
  ! the Alfven speed. 
  double precision                              :: dt, dtout, tmax, t = 0

  ! ===========================================================================
  ! ==================================================================== Fields
  ! ===========================================================================

  ! Electric and magnetic fields in the bulk of the simulation. Also 
  ! field-aligned current. 
  double complex, dimension(:,:), allocatable   :: B1, B2, B3, E1, E2, E3, j3
  ! Curl components of the bulk electric and magnetic fields. The J concedes
  ! that these are actually off by a factor of the Jacobian, which gets bundled
  ! with the coefficients. 
  double complex, dimension(:,:), allocatable   :: JCsup1, JCsup2, JCsup3,    &
                                                   JFsup1, JFsup2, JFsup3
  ! Magnetic fields and scalar magnetic potential at the atmospheric boundary. 
  double complex, dimension(:,:), allocatable   :: B1E, B1I, B2E, B2I, PsiE,  &
                                                   PsiI
  ! Driving fields for compressional and current driving respectively. 
  double precision, dimension(:), allocatable   :: B3drive
  double precision, dimension(:,:), allocatable :: j2drive

  ! ===========================================================================
  ! ============================================================== Coefficients
  ! ===========================================================================

  ! For updating bulk magnetic fields. 
  double precision, dimension(:,:), allocatable   :: B1_JCsup1, B1_JCsup3,    &
                                                     B2_JCsup2, B3_JCsup1,    &
                                                     B3_JCsup3
  ! For updating bulk electric fields and field-aligned current. 
  double precision, dimension(:,:), allocatable   :: E1_E1, E1_E2, E1_E3,     &
                                                     E1_JFsup1, E1_JFsup2,    &
                                                     E1_j2drive, E1_JFsup3,   &
                                                     E1_j3, E2_E1, E2_E2,     &
                                                     E2_E3, E2_JFsup1,        &
                                                     E2_JFsup2, E2_j2drive,   &
                                                     E3_E3, E3_j3, E3_JFsup1, &
                                                     E3_JFsup3, j3_j3, j3_E3
  ! For computing the scalar magnetic potential. 
  double precision, dimension(:,:,:), allocatable :: PsiI_B1, PsiI_B3
  ! We also get ground signatures and harmonic weights. 
  double precision, dimension(:,:), allocatable   :: alpha_Br, PsiE_Br
  double complex, dimension(:, :), allocatable    :: alpha
  ! For computing edge electric field values from the scalar magnetic potential. 
  double precision, dimension(:,:), allocatable   :: E1_B1I, E1_B2I, E2_B1I,  &
                                                     E2_B2I

  ! ===========================================================================
  ! ====================== Linearized Differentiation and Interpolation Weights
  ! ===========================================================================

  ! For computing derivatives with respect to usup1 and usup3. 
  double precision, allocatable, dimension(:,:) :: d1w, d3w
  ! Each field is defined only on a single parity in each dimension. When an 
  ! off-parity field value is needed, it's interpolated with these weights. 
  double precision, allocatable, dimension(:,:) :: i1w, i3w

  ! ===========================================================================
  ! ======================================================== Driving Parameters
  ! ===========================================================================

  ! Waveform index, angular frequency, characteristic time. 
  integer                                         :: idrive
  double precision                                :: wdrive, tdrive
  ! If driving with a spectrum, we also need frequencies and phase offsets. 
  integer, parameter                              :: nspectrum = 20
  double precision, dimension(1:nspectrum)        :: wdrives, pdrives

  contains

  ! ===========================================================================
  ! =========================================== Time Component of Driving Field
  ! ===========================================================================

  ! The driving is broken down into a static array (either B3drive or j2drive)
  ! which indicates the relative strangth of driving at different locations, 
  ! and a time-dependent scale factor, which is recomputed each time step. The
  ! scale factor waveform depends on the input parameter idrive. 
  double precision function getDriveScale(t0)
    double precision, intent(in) :: t0
    ! The shape of the driving waveform in time is set by idrive. 
    select case (idrive)
      ! Sine wave. 
      case(1)
        getDriveScale = sin(wdrive*t0)
      ! Pulse. 
      case(2)
        if(t0 .le. tdrive) then
          getDriveScale = 0.5*( 1 - cos(2*pi*t0/tdrive) )
        else
          getDriveScale = 0
        endif
      ! Ramp to constant. 
      case(3)
        if(t0 .le. tdrive) then
          getDriveScale = t0/tdrive
        else
          getDriveScale = 1
        endif
      ! Ramp to a spectrum. 
      case(4)
        ! Sum over the frequencies and phase offsets computed earlier. 
        getDriveScale = sum( sin( wdrives(:)*t0 + pdrives(:) ) )
        ! If we're still ramping, scale it down proportionally. 
        if (t0 .le. tdrive) getDriveScale = getDriveScale*t0/tdrive
      ! Wave packet. 
      case(5)
        if(t0 .le. tdrive) then
          getDriveScale = 0.5*( 1 - cos(2*pi*t0/tdrive) )*cos(wdrive*t0)
        else
          getDriveScale = 0
        endif
      ! Damped oscillator. 
      case(6)
        getDriveScale = sin(wdrive*t0)*exp(-1.*t0/tdrive)
    end select
  end function getDriveScale

  ! ===========================================================================
  ! ======================== Differentiation and Interpolation Helper Functions
  ! ===========================================================================

  ! These functions are actually used during the main loop. We include them in
  ! the loop module in hopes that the compiler will inline them. 

  ! ---------------------------------------------------------------------------
  ! ----------------------------------------------- Differentiation (Azimuthal)
  ! ---------------------------------------------------------------------------

  ! There are no boundaries to worry about in the azimuthal direction. 
  double complex function d2(field, i, k)
    double complex, intent(in), dimension(0:, 0:) :: field
    integer, intent(in)                           :: i, k
    d2 = (0., 1.)*azm*field(i, k)
  end function d2

  ! ---------------------------------------------------------------------------
  ! --------------------------- Differentiation (Dirichlet Boundary Conditions)
  ! ---------------------------------------------------------------------------

  ! With a Dirichlet boundary condition, we specify the value of the function
  ! at the edge. We choose to have the fields go to zero. 

  ! Alternatively, we can try specifying that the second derivative goes to
  ! zero. That's probably... similar? 

  double complex function d1D(field, i, k)
    double complex, intent(in), dimension(0:, 0:) :: field
    integer, intent(in)                           :: i, k
    if (i .eq. 0) then
      d1D = d1w(0, 1)*field(1, k)
!      d1D = d1w(0, 0)*field(1, k) + d1w(0, 1)*field(3, k)
    else if (i .eq. n1) then
      d1D = d1w(n1, 0)*field(n1-1, k)
!      d1D = d1w(n1, 0)*field(n1-3, k) + d1w(n1, 1)*field(n1-1, k)
    else
      d1D = d1w(i, 0)*field(i-1, k) + d1w(i, 1)*field(i+1, k)
    end if
  end function d1D

  double complex function d3D(field, i, k)
    double complex, intent(in), dimension(0:, 0:) :: field
    integer, intent(in)                           :: i, k
    if (k .eq. 0) then
      d3D = d3w(0, 1)*field(i, 1)
!      d3D = d3w(0, 0)*field(i, 1) + d3w(0, 1)*field(i, 3)
    else if (k .eq. n3) then
      d3D = d3w(n3, 0)*field(i, n3-1)
!      d3D = d3w(n3, 0)*field(i, n3-3) + d3w(n3, 1)*field(i, n3-1)
    else
      d3D = d3w(k, 0)*field(i, k-1) + d3w(k, 1)*field(i, k+1)
    end if
  end function d3D

  ! ---------------------------------------------------------------------------
  ! ----------------------------- Differentiation (Neumann Boundary Conditions)
  ! ---------------------------------------------------------------------------

  ! Under Neumann boundary conditions, we specify the value of the derivative. 
  ! We choose to have field derivatives go to zero. 

  double complex function d1N(field, i, k)
    double complex, intent(in), dimension(0:, 0:) :: field
    integer, intent(in)                           :: i, k
    if (i .eq. 0) then
      d1N = 0
    else if (i .eq. n1) then
      d1N = 0
    else
      d1N = d1w(i, 0)*field(i-1, k) + d1w(i, 1)*field(i+1, k)
    end if
  end function d1N

  double complex function d3N(field, i, k)
    double complex, intent(in), dimension(0:, 0:) :: field
    integer, intent(in)                           :: i, k
    if (k .eq. 0) then
      d3N = 0
    else if (k .eq. n3) then
      d3N = 0
    else
      d3N = d3w(k, 0)*field(i, k-1) + d3w(k, 1)*field(i, k+1)
    end if
  end function d3N

  ! ---------------------------------------------------------------------------
  ! ----------------------------- Interpolation (Dirichlet Boundary Conditions)
  ! ---------------------------------------------------------------------------

  ! With a Dirichlet boundary condition, we specify the value of the function 
  ! at the edge. We choose to have the fields go to zero. 

  ! Alternatively, we can try specifying that the second derivative goes to
  ! zero. That's probably... similar? 

  double complex function i1D(field, i, k)
    double complex, intent(in), dimension(0:, 0:) :: field
    integer, intent(in)                           :: i, k
    if (i .eq. 0) then
      i1D = 0
!      i1D = i1w(0, 0)*field(1, k) + i1w(0, 1)*field(3, k)
    else if (i .eq. n1) then
      i1D = 0
!      i1D = i1w(n1, 0)*field(n1-3, k) + i1w(n1, 1)*field(n1-1, k)
    else
      i1D = i1w(i, 0)*field(i-1, k) + i1w(i, 1)*field(i+1, k)
    end if
  end function i1D

  double complex function i3D(field, i, k)
    double complex, intent(in), dimension(0:, 0:) :: field
    integer, intent(in)                           :: i, k
    if (k .eq. 0) then
      i3D = 0
!      i3D = i3w(0, 0)*field(i, 1) + i3w(0, 1)*field(i, 3)
    else if (k .eq. n3) then
      i3D = 0
!      i3D = i3w(n3, 0)*field(i, n3-3) + i3w(n3, 1)*field(i, n3-1)
    else
      i3D = i3w(k, 0)*field(i, k-1) + i3w(k, 1)*field(i, k+1)
    end if
  end function i3D

  ! ---------------------------------------------------------------------------
  ! ------------------------------- Interpolation (Neumann Boundary Conditions)
  ! ---------------------------------------------------------------------------

  ! Under Neumann boundary conditions, we specify the value of the derivative. 
  ! We choose to have field derivatives go to zero. That means that the value
  ! of a function at the edge is the same as the value of that function just
  ! inside the edge. 

  double complex function i1N(field, i, k)
    double complex, intent(in), dimension(0:, 0:) :: field
    integer, intent(in)                           :: i, k
    if (i .eq. 0) then
      i1N = field(1, k)
    else if (i .eq. n1) then
      i1N = field(n1-1, k)
    else
      i1N = i1w(i, 0)*field(i-1, k) + i1w(i, 1)*field(i+1, k)
    end if
  end function i1N

  double complex function i3N(field, i, k)
    double complex, intent(in), dimension(0:, 0:) :: field
    integer, intent(in)                           :: i, k
    if (k .eq. 0) then
      i3N = field(i, 1)
    else if (k .eq. n3) then
      i3N = field(i, n3-1)
    else
      i3N = i3w(k, 0)*field(i, k-1) + i3w(k, 1)*field(i, k+1)
    end if
  end function i3N

  ! ===========================================================================
  ! ========================================================== Main Loop Driver
  ! ===========================================================================

  ! This subroutine uses the coefficients to advance the fields through an
  ! amount of time dtOut. It then returns to the main program for error
  ! checking and output. The computation is parallelized using OpenMP. 
  subroutine advanceFields()
    use omp_lib
    double precision :: tLoop, driveScale
    integer          :: h, i, ip, k
    !$omp parallel private(h, i, ip, k)
    tLoop = 0
    do while (tloop < dtOut)

      ! -----------------------------------------------------------------------
      ! ---------------------------------------------------- Get Driving Factor
      ! -----------------------------------------------------------------------

      !$omp single
      driveScale = getDriveScale(t + tLoop)
      !$omp end single

      ! -----------------------------------------------------------------------
      ! ------------------------------------------------------------- Compute C
      ! -----------------------------------------------------------------------

      ! Compute JCsup1 on odd i, odd k.
      !$omp do
      do k=1,n3,2
        do i=1,n1,2
          JCsup1(i, k) = d2(E3, i, k) - d3D(E2, i, k)
        end do
      end do
      !$omp end do
      ! Compute JCsup2 on even i, odd k.
      !$omp do
      do k=1,n3,2
        do i=0,n1,2
          JCsup2(i, k) = d3D(E1, i, k) - d1D(E3, i, k)
        end do
      end do
      !$omp end do
      ! Compute JCsup3 on even i, even k. 
      !$omp do
      do k=0,n3,2
        do i=0,n1,2
          JCsup3(i, k) = d1D(E2, i, k) - d2(E1, i, k)
        end do
      end do
      !$omp end do

      ! -----------------------------------------------------------------------
      ! ------------------------------------------------------------- Advance B
      ! -----------------------------------------------------------------------

      ! Interpolate C3 to odd i. 
      !$omp do
      do k=0,n3,2
        do i=1,n1,2
          JCsup3(i, k) = i1N(JCsup3, i, k)
        end do
      end do
      !$omp end do
      ! Interpolate C3 to odd k, then advance B1 on odd i, odd k. 
      !$omp do
      do k=1,n3,2
        do i=1,n1,2
          JCsup3(i, k) = i3N(JCsup3, i, k)
          B1(i, k) = B1(i, k) + B1_JCsup1(i, k)*JCsup1(i, k) +                &
                                B1_JCsup3(i, k)*JCsup3(i, k)
        end do
      end do
      !$omp end do
      ! Advance B2 on even i, odd k.
      !$omp do
      do k=1,n3,2
        do i=0,n1,2
          B2(i, k) = B2(i, k) + B2_JCsup2(i, k)*JCsup2(i, k)
        end do
      end do
      !$omp end do
      ! Interpolate C1 to even i. 
      !$omp do
      do k=1,n3,2
        do i=0,n1,2
          JCsup1(i, k) = i1N(JCsup1, i, k)
        end do
      end do
      !$omp end do
      ! Interpolate C1 to even k, then advance B3 on even i, even k. 
      !$omp do
      do k=0,n3,2
        do i=0,n1,2
          JCsup1(i, k) = i3N(JCsup1, i, k)
          B3(i, k) = B3(i, k) + B3_JCsup1(i, k)*JCsup1(i, k) +                &
                                B3_JCsup3(i, k)*JCsup3(i, k)
        end do
      end do
      !$omp end do

      ! -----------------------------------------------------------------------
      ! ---------------------------------------------------- Apply Driving to B
      ! -----------------------------------------------------------------------

      !$omp single
      B3(n1, 0:n3:2) = B3drive(0:n3:2)*driveScale
      !$omp end single

!      ! -----------------------------------------------------------------------
!      ! ---------------------------------------- Apply Boundary Conditions to B
!      ! -----------------------------------------------------------------------

!      ! Magnetic fields use Neumann boundary conditions: the derivative is zero
!      ! at the boundary. We set this explicitly on the inner and outer
!      ! boundaries. Ionospheric boundaries are handled using Psi. 

!      ! B1 is defined at odd i, so it isn't explicitly valued at the edge. 

!      ! B2 is defined on even i, odd k. 
!      !$omp single
!      B2(0, 1:n3:2) = B2(2, 1:n3:2)
!      B2(n1, 1:n3:2) = B2(n1-2, 1:n3:2)
!      !$omp end single

!      ! B3 is defined on even i, even k. We don't touch the outer boundary,
!      ! since that's where driving is delivered. 
!      !$omp single
!      B3(0, 0:n3:2) = B3(2, 0:n3:2)
!      !$omp end single

      ! -----------------------------------------------------------------------
      ! ------------------------------------------------------------- Advance j
      ! -----------------------------------------------------------------------

      ! Like the magnetic field components, the current is offset from the
      ! electric field by half a time step. That means we want to update the
      ! current, then update the electric fields based on the new current
      ! values, just like we do with the magnetic fields. 

      ! Advance j3 on odd i, odd k. 
      !$omp do
      do k=1,n3,2
        do i=1,n1,2
          j3(i, k) = j3_j3(i, k)*j3(i, k) + j3_E3(i, k)*E3(i, k)
        end do
      end do
      !$omp end do

      ! -----------------------------------------------------------------------
      ! ------------------------------------------------------------- Compute F
      ! -----------------------------------------------------------------------

      ! Compute JFsup1 on even i, even k.
      !$omp do
      do k=0,n3,2
        do i=0,n1,2
          JFsup1(i, k) = d2(B3, i, k) - d3N(B2, i, k)
        end do
      end do
      !$omp end do
      ! Compute JFsup2 on odd i, even k. 
      !$omp do
      do k=0,n3,2
        do i=1,n1,2
          JFsup2(i, k) = d3N(B1, i, k) - d1N(B3, i, k)
        end do
      end do
      !$omp end do
      ! Compute JFsup3 on odd i, odd k. 
      !$omp do
      do k=1,n3,2
        do i=1,n1,2
          JFsup3(i, k) = d1N(B2, i, k) - d2(B1, i, k)
        end do
      end do
      !$omp end do

      ! -----------------------------------------------------------------------
      ! ------------------------------------------------------------- Advance E
      ! -----------------------------------------------------------------------

      ! There's a lot of interdependence between electric field components --
      ! E1 and E2 both depend on each other, for instance. Luckily, they are
      ! all defined on different parities. As long as we do all of our
      ! interpolation first, we're safe. When we update E1, we'll still have
      ! old values of E1 on off-parity points, which is what we need for E2. 

      ! Interpolate E1 and JFsup1 to odd i, even k (for E2). Interpolate JFsup1
      ! to odd i, odd k (for E3). 
      !$omp do
      do k=0,n3,2
        do i=1,n1,2
          E1(i, k) = i1D(E1, i, k)
          JFsup1(i, k) = i1D(JFsup1, i, k)
        end do
      end do
      !$omp end do
      !$omp do
      do k=1,n3,2
        do i=1,n1,2
          JFsup1(i, k) = i3D(JFsup1, i, k)
        end do
      end do
      !$omp end do
      ! Interpolate E2 and JFsup2 to even i, even k (for E1). 
      !$omp do
      do k=0,n3,2
        do i=0,n1,2
          E2(i, k) = i1D(E2, i, k)
          JFsup2(i, k) = i1D(JFsup2, i, k)
        end do
      end do
      !$omp end do
      ! Interpolate E3, JFsup3, and j3 to even i, even k (for E1). Along the
      ! way, interpolate E3 at odd i, even k (for E2). 
      !$omp do
      do k=0,n3,2
        do i=1,n1,2
          E3(i, k) = i3D(E3, i, k)
          JFsup3(i, k) = i3D(JFsup3, i, k)
          j3(i, k) = i3D(j3, i, k)
        end do
      end do
      !$omp end do
      !$omp do
      do k=0,n3,2
        do i=0,n1,2
          E3(i, k) = i1D(E3, i, k)
          JFsup3(i, k) = i1D(JFsup3, i, k)
          j3(i, k) = i1D(j3, i, k)
        end do
      end do
      !$omp end do
      ! Advance E1 on even i, even k. 
      !$omp do
      do k=0,n3,2
        do i=0,n1,2
          ! Note that E1 depends on tons of different things. This is because 
          ! there are contributions from both Esup1 (in the xhat direction) and
          ! E3 (in the zhat direction). 
          E1(i, k) = E1_E1(i, k)*E1(i, k) + E1_JFsup1(i, k)*JFsup1(i, k) +    &
                     E1_E2(i, k)*E2(i, k) + E1_JFsup2(i, k)*JFsup2(i, k) +    &
                     E1_E3(i, k)*E3(i, k) + E1_JFsup3(i, k)*JFsup3(i, k) +    &
                     E1_j3(i, k)*j3(i, k)
        end do
      end do
      !$omp end do
      ! Advance E2 on odd i, even k.
      !$omp do
      do k=0,n3,2
        do i=1,n1,2
          E2(i, k) = E2_E1(i, k)*E1(i, k) + E2_JFsup1(i, k)*JFsup1(i, k) +    &
                     E2_E2(i, k)*E2(i, k) + E2_JFsup2(i, k)*JFsup2(i, k) +    &
                     E2_E3(i, k)*E3(i, k)
        end do
      end do
      !$omp end do
      ! Advance E3 on odd i, odd k. 
      !$omp do
      do k=1,n3,2
        do i=1,n1,2
          E3(i, k) =                        E3_JFsup1(i, k)*JFsup1(i, k) +    &
                     E3_E3(i, k)*E3(i, k) + E3_JFsup3(i, k)*JFsup3(i, k) +    &
                     E3_j3(i, k)*j3(i, k)
        end do
      end do
      !$omp end do

      ! -----------------------------------------------------------------------
      ! ---------------------------------------------------- Apply Driving to E
      ! -----------------------------------------------------------------------

      ! It would be marginally more efficient to handle this above, but it's 
      ! more legible here. 

      ! Apply current driving to E1 on even i, even k. 
      !$omp do
      do k=0,n3,2
        do i=0,n1,2
          E1(i, k) = E1(i, k) + E1_j2drive(i, k)*j2drive(i, k)*driveScale
        end do
      end do
      !$omp end do
      ! Apply current driving to E2 on odd i, even k.
      !$omp do
      do k=0,n3,2
        do i=1,n1,2
          E2(i, k) = E2(i, k) + E2_j2drive(i, k)*j2drive(i, k)*driveScale
        end do
      end do
      !$omp end do

      ! -----------------------------------------------------------------------
      ! ----------------------------------------------------------- Compute Psi
      ! -----------------------------------------------------------------------

      ! Interpolate B1 and B3 to both i at the edges. 
      !$omp do
      do i=1,n1,2
        do h=0,1
          B1(i, h*n3)  = i3N(B1, i, h*n3)
          B3(i, h*n3)  = i1N(B3, i, h*n3)
        end do
      end do
      !$omp end do
      !$omp do
      do i=0,n1,2
        do h=0,1
          B1(i, h*n3)  = i1N(B1, i, h*n3)
        end do
      end do
      !$omp end do
      ! Compute PsiI at both i at the edges. 
      !$omp do
      do i=0,n1
        do h=0,1
          PsiI(i, h) = sum( PsiI_B1(i, :, h)*B1(:, h*n3) +                    &
                            PsiI_B3(i, :, h)*B3(:, h*n3) )
        end do
      end do
      !$omp end do

      ! -----------------------------------------------------------------------
      ! --------------------------------------------------- Advance Edge Fields
      ! -----------------------------------------------------------------------

      ! Interpolate B2 to both i at the edges. 
      !$omp do
      do i=0,n1,2
        do h=0,1
          B2(i, h*n3)  = i3N(B2, i, h*n3)
        end do
      end do
      !$omp end do
      !$omp do
      do i=1,n1,2
        do h=0,1
          B2(i, h*n3)  = i1N(B2, i, h*n3)
        end do
      end do
      !$omp end do
      ! Compute B1I and B2I, the fields below the ionospheric current sheet, 
      ! from PsiI on all i. 
      !$omp do
      do i=0,n1
        do h=0,1
          B1I(i, h) = d1D(PsiI, i, h)
          B2I(i, h) = d2(PsiI, i, h)
        end do
      end do
      !$omp end do
      ! Based on the jump over the current sheet, advance E1 on even i. 
      !$omp do
      do i=0,n1,2
        do h=0,1
          E1(i, h*n3) = E1_B1I(i, h)*( B1(i, h*n3) - B1I(i, h) ) +            &
                        E1_B2I(i, h)*( B2(i, h*n3) - B2I(i, h) )
        end do
      end do
      !$omp end do
      ! Based on the jump over the current sheet, advance E2 on odd i. 
      !$omp do
      do i=1,n1,2
        do h=0,1
          E2(i, h*n3) = E2_B1I(i, h)*( B1(i, h*n3) - B1I(i, h) ) +            &
                        E2_B2I(i, h)*( B2(i, h*n3) - B2I(i, h) )
        end do
      end do
      !$omp end do

!      ! -----------------------------------------------------------------------
!      ! ---------------------------------------- Apply Boundary Conditions to E
!      ! -----------------------------------------------------------------------

!      ! Electric fields use Dirichlet boundary conditions: they go to zero at
!      ! the boundary. We set this explicitly on the inner and outer boundaries.
!      ! Ionospheric boundaries are handled using Psi. 

!      ! E1 is defined on even i, even k. 
!      !$omp single
!      E1(0, 0:n3:2) = E1(2, 0:n3:2)
!      E1(n1, 0:n3:2) = E1(n1-2, 0:n3:2)
!      !$omp end single

!      ! E2 is defined on odd i, so we don't have to worry about it. 

!      ! E3 is defined on odd i, so we don't have to worry about it. 

      ! -----------------------------------------------------------------------
      ! ------------------------------------------------------ End of Time Loop
      ! -----------------------------------------------------------------------

      !$omp single
      tLoop = tLoop + dt
      !$omp end single

    end do
    !$omp end parallel
    ! Apply the elapsed time to the global time variable. 
    t = t + dtOut
  end subroutine advanceFields

  ! ===========================================================================
  ! =================================================== End of Main Loop Module
  ! ===========================================================================

end module loop

! #############################################################################
! ######################################################### Input/Output Module
! #############################################################################

! The IO module also exists at the top level, bringing in no other modules. 
! It's used for accessing files, and provides convenient output handling. 

module io
  implicit none
  ! Location of the file with input parameters. 
  character(128) :: paramfile = 'params.in'
  ! Location of the ionospheric profiles. Change the 0 to select a model. 
  character(21)  :: ionosfile = '../models/ionpar0.dat'
  ! Format for writing out eigenvectors for inspection. 
  character(11)  :: evecfile = 'evec000.out'

  contains

  ! ===========================================================================
  ! ===================================================== Count Lines in a File
  ! ===========================================================================

  integer function lines(filename)
    character(len=*), intent(in) :: filename
    character(128)               :: line
    ! Offset by -1 so we can iterate from 0 to lines. 
    lines = -1
    open(unit=99, file=filename(1:len_trim(filename)), action='read')
    do 
      read(99,*,end=10) line
      lines = lines + 1
    end do
    10 close(99)
  end function lines

  ! ===========================================================================
  ! =================================================== Read Column from a File
  ! ===========================================================================

  function readColumn(filename, column)
    character(len=*), intent(in)                :: filename
    integer, intent(in)                         :: column
    integer                                     :: line, nLines
    double precision, dimension(0:column)       :: lineData
    double precision, allocatable, dimension(:) :: readColumn
    nLines = lines(filename)
    allocate( readColumn(0:nLines) )
    open(unit=99, file=filename, action='read')
    ! Grab columnNum values, store the last one, and move on to the next line. 
    do line=0,nLines
      read(99,*,end=20) lineData
      readColumn(line) = lineData(column)
    end do
    20 close(99)
  end function readColumn

  ! ===========================================================================
  ! ======================================================= Get Parameter Value
  ! ===========================================================================

  double precision function readParam(varname)
    character(len=*), intent(in) :: varname
    character(128)               :: label, eq
    open(unit=99, file=paramfile(1:len_trim(paramfile)), action='read')
    do
      ! Read parameter label and value. If it matches varname, return it. 
      read(99,*,end=30) label, eq, readParam
      if ( varname .eq. label(1:len_trim(label)) ) then
        close(99)
        return
      end if
    end do
    30 close(99)
    ! If we don't find the parameter in the file, we use the default value. 
    readParam = defaultParam(varname)
  end function readParam

  ! ===========================================================================
  ! ================================================== Default Parameter Values
  ! ===========================================================================

  double precision function defaultParam(varname)
    character(len=*), intent(in) :: varname
    ! Geometric parameters.
    if (varname .eq. 'n1'  ) defaultParam = 150        ! Number of field lines.
    if (varname .eq. 'n3'  ) defaultParam = 350        ! Grid points per line. 
    if (varname .eq. 'lmin') defaultParam = 2          ! Innermost L value. 
    if (varname .eq. 'lmax') defaultParam = 10         ! Outermost L value. 
    if (varname .eq. 'sfac') defaultParam = 1.03       ! Geometric spacing
                                                       ! factor along outermost
                                                       ! field line. 
    if (varname .eq. 'zi'  ) defaultParam = 100        ! Ionosphere height, km.
    if (varname .eq. 'azm' ) defaultParam = 0          ! Azimuthal modenumber. 
    if (varname .eq. 'modes'  ) defaultParam = 0.25    ! Harmonics to keep. 
    ! Time parameters.
    if (varname .eq. 'tmax'  ) defaultParam = 10.      ! Simulation time, s. 
    if (varname .eq. 'dtout' ) defaultParam = 1.       ! Output period, s. 
    if (varname .eq. 'cour'  ) defaultParam = 0.1      ! Courant condition. 
    ! Parallel physics handling parameters. 
    if (varname .eq. 'epsfac'  ) defaultParam = -1.    ! Boris factor for eps0.
                                                       ! Negative is automatic.
    if (varname .eq. 'inertia'  ) defaultParam = -1.   ! Positive includes
                                                       ! electron inertia. 
    if (varname .eq. 'fudge' ) defaultParam = 0.2      ! Factor to stabilize
                                                       ! plasma oscillations. 
    ! Physical parameter profiles.
    if (varname .eq. 'model') defaultParam = 1         ! Ionospheric profile. 
                                                       ! 1 is active dayside. 
                                                       ! 2 is quiet dayside. 
                                                       ! 3 is active nightside. 
                                                       ! 4 is quiet nightside. 
                                                       ! 5 is QUIET flank. 
                                                       ! 6 is ACTIVE flank. 
    if (varname .eq. 'naz'  ) defaultParam = 10        ! Density at auroral
                                                       ! base, cm^-3. 
    if (varname .eq. 'haz'  ) defaultParam = 1         ! Scale height (in RE)
                                                       ! of auroral density. 
    if (varname .eq. 'nps'  ) defaultParam = 1e4       ! Number density, cm^-3,
                                                       ! at plasmasphere base. 
    if (varname .eq. 'lps'  ) defaultParam = 1.0857    ! Plasmasphere scale L.
    if (varname .eq. 'lpp'  ) defaultParam = 4         ! Plasmapause location, 
                                                       ! RE at equator. 
    if (varname .eq. 'dlpp' ) defaultParam = 0.1       ! Plasmapause thickness,  
                                                       ! RE at equator. 
    ! Debugging factors, set to 0 to turn things off. 
    if (varname .eq. 'i1fac' ) defaultParam = 1.       ! Interpolation in i. 
    if (varname .eq. 'i3fac' ) defaultParam = 1.       ! Interpolation in k. 
    if (varname .eq. 'sighfac' ) defaultParam = 1.     ! Hall conductivity. 
    if (varname .eq. 'sigpfac' ) defaultParam = 1.     ! Pedersen conductivity.
    ! Drive parameters.
    if (varname .eq. 'idrive'   ) defaultParam = 1     ! Waveform index. 
    if (varname .eq. 'bdrive'   ) defaultParam = 0.    ! Strength of driving B
                                                       ! field, nT. 
    if (varname .eq. 'jdrive'   ) defaultParam = 0.    ! Strength of driving
                                                       ! current, uA/m^2. 
    if (varname .eq. 'fdrive'   ) defaultParam = 0.015 ! Frequency, Hz. 
    if (varname .eq. 'tdrive'   ) defaultParam = 60.   ! Ramp/wave packet
                                                       ! duration, s. 
    if (varname .eq. 'latdrive' ) defaultParam = 5.    ! Latitude, degrees. 
    if (varname .eq. 'dlatdrive') defaultParam = 5.    ! Spread in latitude. 
    if (varname .eq. 'ldrive'   ) defaultParam = 5.    ! L shell for driving
                                                       ! current. 
    if (varname .eq. 'dldrive'  ) defaultParam = 0.5   ! Spread in L shell. 
    ! Integrated atmospheric conductivities. Negative means automatic. 
    if (varname .eq. 'sig0atm') defaultParam = -1      ! Parallel. 
    if (varname .eq. 'sighatm') defaultParam = -1      ! Hall. 
    if (varname .eq. 'sigpatm') defaultParam = -1      ! Pedersen
  end function defaultParam

  ! ===========================================================================
  ! ===================================================== Write Parameter Value
  ! ===========================================================================

  ! While most of the output is in the form of arrays (in DAT files), we also
  ! want to be able to report snapshots of the parameter scales we're working
  ! with. These just go to stdout. 
  subroutine writeParam(key, val, units)
    character(len=*), intent(in)           :: key
    double precision, intent(in)           :: val
    character(len=*), intent(in), optional :: units
    ! Output is lined up nicely in columns, optionally with units.  
    if ( present(units) ) then
      write(*, '(a30, a, es10.2, 2a)') key, ' = ', val, ' ', units
    else
      write(*, '(a30, a, es10.2)') key, ' = ', val
    end if
  end subroutine writeParam

  ! ===========================================================================
  ! ======================================================= Write Complex Array
  ! ===========================================================================

  subroutine writeComplexArray(filename, arr, onlyheader, onlydata, noskip, nt)
    character(len=*), intent(in)   :: filename
    double complex, dimension(:,:) :: arr
    logical, intent(in), optional  :: onlyheader, onlydata, noskip
    integer, intent(in), optional  :: nt
    integer                        :: nx, nz, stride
    ! By default, we print out every other grid point to reduce data size. 
    if ( .not. present(noskip) ) then
      nx = 1 + size(arr, 1)/2
      nz = 1 + size(arr, 2)/2
      stride = 2
    ! If requested, we can print out the entire array instead. 
    else
      nx = size(arr, 1)
      nz = size(arr, 2)
      stride = 1
    end if
    open(unit=99, file=filename, action='write', access='append')
    ! Unless asked not to, print out the header. 
    if ( .not. present(onlydata) ) then
      ! Dimensions of the (sliced) 2d array, optionally including time steps. 
      if (present(nt)) then
        write(99,*) nx, nz, nt
      else
        write(99,*) nx, nz
      end if
    end if
    ! Unless asked not to, print the array. 
    if (.not. present(onlyheader)) then
      ! If there are exactly two columns, always print them both. 
      if (size(arr, 2) .eq. 2) then
        write(99,*) arr(::stride, :)
      else
        write(99,*) arr(::stride, ::stride)
      end if
    end if
  end subroutine writeComplexArray

  ! ===========================================================================
  ! ========================================================== Write Real Array
  ! ===========================================================================

  subroutine writeRealArray(filename, arr, onlyheader, onlydata, noskip, nt)
    character(len=*), intent(in)     :: filename
    double precision, dimension(:,:) :: arr
    logical, intent(in), optional    :: onlyheader, onlydata, noskip
    integer, intent(in), optional    :: nt
    integer                          :: nx, nz, stride
    ! By default, we print out every other grid point to reduce data size. 
    if ( .not. present(noskip) ) then
      nx = 1 + size(arr, 1)/2
      nz = 1 + size(arr, 2)/2
      stride = 2
    ! If requested, we can print out the entire array instead. 
    else
      nx = size(arr, 1)
      nz = size(arr, 2)
      stride = 1
    end if
    open(unit=99, file=filename, action='write', access='append')
    ! Unless asked not to, print out the header. 
    if ( .not. present(onlydata) ) then
      ! Dimensions of the (sliced) 2d array, optionally including time steps. 
      if (present(nt)) then
        write(99,*) nx, nz, nt
      else
        write(99,*) nx, nz
      end if
    end if
    ! Unless asked not to, print the array. 
    if (.not. present(onlyheader)) then
      ! If there are exactly two columns, always print them both. 
      if (size(arr, 2) .eq. 2) then
        write(99,*) arr(::stride, :)
      else
        write(99,*) arr(::stride, ::stride)
      end if
    end if
  end subroutine writeRealArray

  ! ===========================================================================
  ! ======================================================== Write Real Columns
  ! ===========================================================================

  ! Writes up to ten columns. If integers are given, they are cast as doubles. 
  subroutine writeRealColumns(filename, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9)
    character(len=*), intent(in)                          :: filename
    double precision, dimension(0:), intent(in)           :: c0
    double precision, dimension(0:), intent(in), optional :: c1, c2, c3, c4,  &
                                                             c5, c6, c7, c8, c9
    integer                                               :: line, nLines
    open(unit=99, file=filename, action='write', access='append')
    ! Print the number of lines. 
    write(99, *) size(c0)
    ! Assume all columns are of the same length. 
    nLines = size(c0) - 1
    do line=0,nLines
      write(99, '(es12.4)', advance='no') c0(line)
      if ( present(c1) ) write(99, '(es12.4)', advance='no') c1(line)
      if ( present(c2) ) write(99, '(es12.4)', advance='no') c2(line)
      if ( present(c3) ) write(99, '(es12.4)', advance='no') c3(line)
      if ( present(c4) ) write(99, '(es12.4)', advance='no') c4(line)
      if ( present(c5) ) write(99, '(es12.4)', advance='no') c5(line)
      if ( present(c6) ) write(99, '(es12.4)', advance='no') c6(line)
      if ( present(c7) ) write(99, '(es12.4)', advance='no') c7(line)
      if ( present(c8) ) write(99, '(es12.4)', advance='no') c8(line)
      if ( present(c9) ) write(99, '(es12.4)', advance='no') c9(line)
      write(99,*) ''
    end do
    close(99)
  end subroutine writeRealColumns

  ! ===========================================================================
  ! ========================================================== Write an Integer
  ! ===========================================================================

  subroutine writeInteger(filename, val)
    character(len=*), intent(in) :: filename
    integer                      :: val
    open(unit=99, file=filename, action='write', access='append')
    write(99, *) val
    close(99)
  end subroutine writeInteger

  ! ===========================================================================
  ! ============================================================== Write a Real
  ! ===========================================================================

  subroutine writeReal(filename, val)
    character(len=*), intent(in) :: filename
    double precision             :: val
    open(unit=99, file=filename, action='write', access='append')
    write(99, *) val
    close(99)
  end subroutine writeReal

  ! ===========================================================================
  ! ========================================================== End of IO Module
  ! ===========================================================================

end module io

! #############################################################################
! #################################################################### Geometry
! #############################################################################

! The geometry module keeps track of the grid. That includes grid positions, of
! course, as well as nonorthogonal coordinate, differentiation and
! interpolation weight setup, and mapping between covariant, contravariant, and
! physical bases. The solution to Laplace's equation in the ionospheric shell
! is also handled by the geometry module. 

module geometry
  use io
  use loop
  implicit none

  ! ===========================================================================
  ! ============================================================== Grid Spacing
  ! ===========================================================================

  ! Spherical coordinates r and theta. Note that the nonorthogonal coordinates
  ! usup1 and usup3 are not actually stored as arrays; they are computed on the
  ! fly using the functions below. This saves on memory, and doesn't
  ! significantly affect run time since the functions are not invoked during
  ! the main loop. The same is the case for metric tensor components, etc. 
  double precision, allocatable, dimension(:,:) :: r, q

  ! ===========================================================================
  ! ============================= Numerical Harmonics for the Atmospheric Shell
  ! ===========================================================================

  ! The eigenvectors Y are analogous to spherical harmonics (and look a lot
  ! like them, too). They are computed based on the grid spacing and
  ! derivative formulation. 
  double precision, allocatable, dimension(:,:) :: Y, Yinv
  ! Each Y has a corresponding eigenvalue, nu. 
  double precision, allocatable, dimension(:)   :: nu
  ! We can choose to get rid of high-order harmonics, which may exhibit
  ! unstably-large gradients. 
  integer                                       :: nModes

  contains

  ! ===========================================================================
  ! ============================================================ Geometry Setup
  ! ===========================================================================

  ! ---------------------------------------------------------------------------
  ! ------------------------------------------ Establish Grid Size and Indexing
  ! ---------------------------------------------------------------------------

  subroutine indexSetup()
    ! Get the grid size from the parameters file. 
    n1 = readParam('n1')
    n3 = readParam('n3')
    ! Make sure that n1 and n3 are even. 
    if (mod(n1, 2) .ne. 0) n1 = n1 + 1
    if (mod(n3, 2) .ne. 0) n3 = n3 + 1
    ! Define little helpers to make edge handling more convenient. 
    hh = [0, 1]
    ii = [0, n1]
    kk = [0, n3]
    zz = [0, 0]
  end subroutine indexSetup

  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------- Set Grid Along Equator
  ! ---------------------------------------------------------------------------

  ! Field lines are spaced uniformly in usup1 space. This allows us to
  ! determine r and q at the ionosphere and at the equator. Once those are set
  ! up, it'll be safe for other subroutines to call the function usup1(). 
  subroutine usup1Setup()
    double precision :: usup1min, usup1max, usup1temp
    integer          :: i
    ! Grab ionosphere height in km, and turn it into ionospheric radius in Mm. 
    RI = RE + readParam('zi') / 1000.
    ! Field lines are spaced evenly in usup1. This the ionosphere and equator. 
    usup1min = -RI/( RE*readParam('lmin') )
    usup1max = -RI/( RE*readParam('lmax') )
    ! Allocate arrays for spatial coordinates. 
    allocate( r(0:n1, 0:n3), q(0:n1, 0:n3) )
    do i=0,n1
      usup1temp = usup1min + i*(usup1max - usup1min)/n1
      q(i, 0) = acos( sqrt(usup1temp + 1) )
      r(i, 0) = RI
      q(i, n3/2) = pi/2
      r(i, n3/2) = -RI/usup1temp
    end do
  end subroutine usup1Setup

  ! ---------------------------------------------------------------------------
  ! --------------------------------------- Set Grid Along Outermost Field Line
  ! ---------------------------------------------------------------------------

  ! This helper function computes the distance along the outermost field line
  ! from angle qI (theta ionosphere) to qp (theta prime), by integrating
  ! ds = RE*L*sqrt(1 + 3*cosq^2) sinq dq. 
  double precision elemental function s(qp)
    double precision, intent(in) :: qp
    double precision             :: qI
    ! Note that we can't have s depend explicitly on the q array, since s is
    ! used to figure out q. Instead, we compute the colatitude based on L,
    ! which is already set. 
    qI = acos( sqrt( 1 - RI/r(n1, n3/2) ) )
    s = ( sqrt(3.)*asinh( sqrt(3.)*cos(qI) )                                  &
        + 3*cos(qI)*sqrt( 1+3*cos(qI)**2 )                                    &
        - sqrt(3.)*asinh( sqrt(3.)*cos(qp) )                                  &
        - 3*cos(qp)*sqrt( 1+3*cos(qp)**2 ) )*r(n1, n3/2)/6
  end function s

  ! Spacing in usup3 is determined along the outermost field line. 
  subroutine usup3Setup()
    double precision, dimension(0:n3)   :: qMin, qMax, sFinal, sTry
    double precision                    :: sFac
    integer                             :: k, m

    ! Factor by which grid spacing increases along the outermost field line.
    ! This is unexpectedly crucial in getting stability... model 3 in
    ! particular is prone to wigglies when the spacing is too large/small. 
!    sFac = readParam('sfac')
    if (readParam('model') .eq. 3) then
      sFac = 1.03
      write(*,*) 'Found model 3. Setting sfac = ', sfac
    else
      sFac = 1.025
      write(*,*) 'Found model 1, 2, or 4. Setting sfac = ', sfac
    end if

    ! Spacing along the outermost line increases geometrically to the equator,
    ! then is reflected to the southern hemisphere. 
    do k=0,n3/2
      sFinal(k) = s( q(n1, n3/2) )*( 1 - sFac**k )/( 1 - sFac**(n3/2) )
      sFinal(n3 - k) = 2*s( q(n1, n3/2) ) - sFinal(k)
    end do
    ! Placement is determined via bisection search. Fifty iterations is plenty.
    qMin = 0
    qMax = pi
    do m=0,50
      q(n1, :) = 0.5*(qMin + qMax)
      sTry = s( q(n1, :) )
      ! Compare this step's s to the goal. 
      where ( sTry .lt. sFinal )
        ! Where sTry is too small, increase qMin. 
        qMin = q(n1, :)
      elsewhere
        ! Where sTry is too large, decrease qMax. 
        qMax = q(n1, :)
      end where
    end do
    ! After establishing q along the outermost line, use usup1 to get r. 
    r(n1, :) = -RI*sin( q(n1, :) )**2 / maxval( usup1() )
  end subroutine usup3Setup

  ! ---------------------------------------------------------------------------
  ! --------------------------------------- Map from usup1 and usup3 to r and q
  ! ---------------------------------------------------------------------------

  ! Points are placed on the remaining field lines in accordance with the usup1
  ! values given by the equatorial spacing, and the usup3 values given by the
  ! outermost field line. 
  subroutine rqSetup()
    double precision, dimension(0:n1, 0:n3) :: qMin, qMax, usup1Final,        &
                                               usup3Final, usup3Try
    integer                                 :: m
    ! There's no easy mapping from a position in usup1-usup3 space to a
    ! position in r-q space, so we again place points via bisection search. 
    usup1Final = usup1()
    usup3Final = usup3()
    qMin = 0
    qMax = pi
    do m=0,50
      q = 0.5*(qMin + qMax)
      ! Assuming the above q, compute r along the field line, then usup3. 
      r = -RI*sin(q)**2 / usup1Final
      usup3Try = (RI/r)**2 * cos(q) / sqrt(1 + usup1Final)
      ! Compare this step's usup3 to the goal. 
      where ( usup3Try .lt. usup3Final )
        ! Where usup3Try is too large, increase qmax. 
        qMax = q
      elsewhere
        ! Where usup3Try is too large, increase qmin. 
        qMin = q
      end where
    end do
  end subroutine rqSetup

  ! ---------------------------------------------------------------------------
  ! ------------------------------------- Compute Linearized Derivative Weights
  ! ---------------------------------------------------------------------------

  subroutine dwSetup()
    ! Allocate arrays for linearized differential weights. 
    allocate( d1w(0:n1,0:1), d3w(0:n3,0:1) ) 
    ! Weights are determined by the local rate of change of usup1. 
    d1w(:, 0) = -1/( eoshift(usup1_(), 1) - eoshift(usup1_(), -1) )
    ! The derivative at i==0 is the same as the derivative at i==2, since both
    ! are based on the two nearest odd grid points. 
    d1w(0, 0) = d1w(2, 0)
    d1w(n1, 0) = d1w(n1-2, 0)
    d1w(:, 1) = -d1w(:, 0)
    ! Derivatives in the azimuthal direction are based on the modenumber. 
    azm = readParam('azm')
    ! Derivatives in usup3 are the same as those in usup1. 
    d3w(:, 0) = -1/( eoshift(usup3_(), 1) - eoshift(usup3_(), -1) )
    d3w(0, 0) = d3w(2, 0)
    d3w(n3, 0) = d3w(n3-2, 0)
    d3w(:, 1) = -d3w(:, 0)
  end subroutine dwSetup

  ! ---------------------------------------------------------------------------
  ! ---------------------------------- Compute Linearized Interpolation Weights
  ! ---------------------------------------------------------------------------

  function iw(usup)
    double precision, allocatable, dimension(:,:) :: iw
    double precision, dimension(0:), intent(in)   :: usup
    integer                                       :: i, n
    n = size(usup) - 1
    allocate( iw(0:n, 0:1) )
    ! f(0) is computed from f(1) and f(3)
    iw(0, 0) =  ( usup(3) - usup(0) )/( usup(3) - usup(1) )
    ! f(i) is computed from f(i-1) and f(i+1)
    do i=1,size(usup)-2
      iw(i, 0) = ( usup(i+1) - usup(i) )/( usup(i+1) - usup(i-1) )
    end do
    ! f(n) is computed from f(n-3) and f(n-1)
    iw(n, 0) = -( usup(n) - usup(n-1) )/( usup(n-1) - usup(n-3) )
    ! The above and below weights must sum to 1 to be correct to zeroth order. 
    iw(:, 1) = 1. - iw(:, 0)
  end function iw

  subroutine iwSetup()
    ! Allocate arrays for linearized interpolation weights. 
    allocate( i1w(0:n1,0:1), i3w(0:n3,0:1) ) 
    ! Interpolation isn't as well-behaved as differentiation at the boundary. 
    i1w = iw( usup1_() )*readParam('i1fac')
    i3w = iw( usup3_() )*readParam('i3fac')
  end subroutine iwSetup

  ! ---------------------------------------------------------------------------
  ! ------------------------------ Solve for Harmonics in the Atmospheric Shell
  ! ---------------------------------------------------------------------------

  function laplaceMatrix()
    complex, dimension(0:n1, 0:n1)    :: laplaceMatrix
    double precision, dimension(0:n1) :: u1, coeff0, coeff1, coeff2
    double precision                  :: du
    integer                           :: i
    u1 = usup1_()
    du = abs( u1(1) - u1(0) )
    coeff0 = azm**2 / u1
    coeff1 = (2 + 3*u1) / du
    coeff2 = 4*u1*(1 + u1) / du**2
    laplaceMatrix = 0
    do i=0,n1
      if (i .ne. 0) laplaceMatrix(i-1, i) = coeff2(i) - coeff1(i)
      laplaceMatrix(i, i) = -coeff0(i) - 2*coeff2(i)
      if (i .ne. n1) laplaceMatrix(i+1, i) = coeff2(i) + coeff1(i)
    end do
  end function laplaceMatrix

  subroutine harmonicSetup()
    ! Intel matrix library. 
    use lapack95
    interface
      ! Eigenvector solver. 
      subroutine geev(mat, val, vl, vr, info)
        implicit none
        complex, dimension(:,:), intent(inout)         :: mat
        complex, dimension(:), intent(out)             :: val
        complex, dimension(:,:), intent(out), optional :: vl,vr
        integer, intent(out), optional                 :: info
      end subroutine geev
      ! Decomposition of matrix. 
      subroutine getrf(eve, ipiv, info)
        implicit none
        double precision, dimension(:,:), intent(inout) :: eve
        integer, dimension(:), intent(out)              :: ipiv
        integer, intent(out)                            :: info
      end subroutine getrf
      ! Inversion from decomposed matrix. 
      subroutine getri(eve, ipiv, info)
        implicit none
        double precision, dimension(:,:), intent(inout) :: eve
        integer, dimension(:), intent(in)               :: ipiv
        integer, intent(out)                            :: info
      end subroutine getri
    end interface
    ! The geev subroutine expects single precision complex values. 
    complex, dimension(0:n1, 0:n1)    :: matrix, evecs
    complex, dimension(0:n1)          :: evals
    double precision, dimension(0:n1) :: Yswap
    double precision                  :: nuswap
    integer, dimension(0:n1)          :: ipiv
    integer                           :: info, inu, inup
    ! Allocate the arrays we'll use to hold our harmonics. 
    allocate( Y(0:n1, 0:n1), Yinv(0:n1, 0:n1), nu(0:n1) )
    ! Grab the matrix of Laplace's equation and solve for its eigenvectors. 
    matrix = laplaceMatrix()
    call geev(matrix, evals, vr=evecs)
    ! The eigenvalues are nu*(nu+1). Get nu, the values we actually care about.
    nu = real( 0.5*(sqrt(1 + 4*evals) - 1) )
    ! The eigenvectors are our numerical harmonics. They should not have an
    ! imaginary component. 
    Y = real(evecs)
    ! Sort harmonics by nu value. 
    do inu=0,n1
      do inup=inu+1,n1
        ! If we have a higher eigenvalue later in the list, do nothing. 
        if ( abs( nu(inu) ) < abs( nu(inup) ) ) cycle
        ! Otherwise, swap this eigenvector with that one. 
        Yswap = Y(:, inu)
        Y(:, inu) = Y(:, inup)
        Y(:, inup) = Yswap
        ! Also swap their eigenvalues. 
        nuswap = nu(inu)
        nu(inu) = nu(inup)
        nu(inup) = nuswap
      end do
    end do
    ! Invert the harmonics. This guarantees orthonormality. 
    Yinv = Y
    call getrf(Yinv, ipiv, info)
    call getri(Yinv, ipiv, info)
    ! High order harmonics exhibit very large gradients. We can leave them out
    ! if we're worried about stability. The modes parameter can specify the
    ! number of modes to use (if greater than 1) or the fraction of modes to
    ! use (if less than or equal to 1). 
    if (readParam('modes') .gt. 1) then
      nModes = int( readParam('modes') )
    else
      nModes = int( readParam('modes')*n1 )
    end if
    Y(:, nModes:n1) = 0
    Yinv(nModes:n1, :) = 0
  end subroutine harmonicSetup

  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------- Write Geometric Arrays
  ! ---------------------------------------------------------------------------

  subroutine writeGeometry()
    integer :: inu
    call writeRealArray('q.dat', q)
    ! Convert distances from Mm to RE for easier plotting. 
    call writeRealArray('r.dat', r/RE)
    ! Write out the eigenvalues and eigenvectors at the ionosphere. 
    call writeRealColumns( 'evals.dat', nu(0:nModes-1) )
    call writeRealArray('evecs.dat', Yinv(0:nModes-1, :), noskip=.true.)
  end subroutine writeGeometry

  ! ---------------------------------------------------------------------------
  ! --------------------------------------------------------- Grid Setup Driver
  ! ---------------------------------------------------------------------------

  ! Establish the placement of grid points. The grid is well-behaved in
  ! usup1-usup3 space, so we start by establishing usup1 and usup3. From there,
  ! we compute each grid point's position in spherical coordinates. To minimize
  ! memory use, only arrays used during the main loop are stored; values used
  ! only during setup are computed on the fly. 
  subroutine geometrySetup()
    ! Determine grid size and allocate arrays. 
    call indexSetup()
    ! Place grid points along the equator to determine usup1. 
    call usup1Setup()
    ! Place grid points along the outermost field line to determine usup3. 
    call usup3Setup()
    ! Based on usup1 and usup3, determine r and q everywhere.
    call rqSetup()
    ! Compute linearized differentiation and interpolation weights. 
    call dwSetup()
    call iwSetup()
    ! Solve for the harmonics in the atmospheric shell. 
    call harmonicSetup()
    ! Write geometric arrays to file. 
    call writeGeometry()
    ! Also grab the total simulation time and output interval. 
    tmax = readParam('tmax')
    dtout = readParam('dtout')
  end subroutine geometrySetup

  ! ===========================================================================
  ! ======================================================== Geometric "Arrays"
  ! ===========================================================================

  ! ---------------------------------------------------------------------------
  ! -------------------------------------------------------- McIlwain Parameter
  ! ---------------------------------------------------------------------------

  function L_()
    double precision, dimension(0:n1) :: L_
    L_ = r(:, 0)/( RE*sin( q(:, 0) )**2 )
  end function L_

  function L()
    double precision, dimension(0:n1, 0:n3) :: L
    L = r/( RE*sin(q)**2 )
  end function L

  ! ---------------------------------------------------------------------------
  ! ------------------------------------------------------ Invariant Colatitude
  ! ---------------------------------------------------------------------------

  function q0()
    double precision, dimension(0:n1, 0:n3) :: q0
    q0 = spread( q(:, 0), 2, n3+1)
  end function q0

  ! ---------------------------------------------------------------------------
  ! ------------------------------------------ Nonorthogonal Dipole Coordinates
  ! ---------------------------------------------------------------------------

  function usup1_()
    double precision, dimension(0:n1) :: usup1_
    usup1_ = -RI/r(:, n3/2)
  end function usup1_

  function usup1()
    double precision, dimension(0:n1, 0:n3) :: usup1
    usup1 = spread(usup1_(), 2, n3+1)
  end function usup1

  function usup3_()
    double precision, dimension(0:n3) :: usup3_
    usup3_ = ( RI/r(n1, :) )**2 * cos( q(n1, :) )/cos( q(n1, 0) )
  end function usup3_

  function usup3()
    double precision, dimension(0:n1, 0:n3) :: usup3
    usup3 = transpose( spread(usup3_(), 2, n1 + 1) )
  end function usup3

  ! ---------------------------------------------------------------------------
  ! ----------------------------------------------- Covariant Metric Components
  ! ---------------------------------------------------------------------------

  function g11()
    double precision, dimension(0:n1, 0:n3) :: g11
    g11 = ( (1 - RI/r)**2 + ( 0.5*( 1 + 3*cos( q0() )**2 )/tan(q) )**2 )*r**4 &
          / ( RI**2 * ( 1 + 3*cos(q)**2 )**2 * cos( q0() )**4 )
  end function g11

  function g13()
    double precision, dimension(0:n1, 0:n3) :: g13
    g13 = r**4 * 0.5*cos(q)/( RI**2 * cos( q0() ) * ( 1 + 3*cos(q)**2 ) )
  end function g13

  function g22()
    double precision, dimension(0:n1, 0:n3) :: g22
    g22 = ( r*sin(q) )**2
  end function g22

  function g31()
    double precision, dimension(0:n1, 0:n3) :: g31
    g31 = r**4 * 0.5*cos(q)/( RI**2 * cos( q0() ) * ( 1 + 3*cos(q)**2 ) )
  end function g31

  function g33()
    double precision, dimension(0:n1, 0:n3) :: g33
    g33 = r**6 * cos( q0() )**2 / ( RI**4 * ( 1 + 3*cos(q)**2 ) )
  end function g33

  ! ---------------------------------------------------------------------------
  ! ------------------------------------------- Contravariant Metric Components
  ! ---------------------------------------------------------------------------

  function gsup11()
    double precision, dimension(0:n1, 0:n3) :: gsup11
    gsup11 = RI**2 * sin(q)**2 * ( 1 + 3*cos(q)**2 ) / r**4
  end function gsup11

  function gsup13()
    double precision, dimension(0:n1, 0:n3) :: gsup13
    gsup13 = -0.5*RI**3 * sin( q0() )**2 * cos(q)*( 1 + 3*cos(q)**2 )         &
              /( cos( q0() )**3 * r**5 )
  end function gsup13

  function gsup22()
    double precision, dimension(0:n1, 0:n3) :: gsup22
    gsup22 = 1/( r*sin(q) )**2
  end function gsup22

  function gsup31()
    double precision, dimension(0:n1, 0:n3) :: gsup31
    gsup31 = -0.5*RI**3 * sin( q0() )**2 * cos(q)*( 1 + 3*cos(q)**2 )         &
             /( cos( q0() )**3 * r**5 )
  end function gsup31

  function gsup33()
    double precision, dimension(0:n1, 0:n3) :: gsup33
    gsup33 = RI**4 * ( ( 0.5*cos(q)*( 1 + 3*cos( q0() )**2 ) )**2             &
                       + ( sin(q)*(1 - RI/r) )**2 ) / ( r*cos( q0() ) )**6
  end function gsup33

  ! ---------------------------------------------------------------------------
  ! ------------------------------------------------------ Jacobian Determinant
  ! ---------------------------------------------------------------------------

  function J()
    double precision, dimension(0:n1, 0:n3) :: J
    J = r**6 * cos( q0() ) / ( RI**3 * ( 1+3*cos(q)**2 ) )
  end function J

  ! ---------------------------------------------------------------------------
  ! -------------------- Scale Factors for Mapping to Field-Aligned Coordinates
  ! ---------------------------------------------------------------------------

  function h1()
    double precision, dimension(0:n1, 0:n3) :: h1
    h1 = 1/sqrt( gsup11() )
  end function h1

  function h2()
    double precision, dimension(0:n1, 0:n3) :: h2
    h2 = 1/sqrt( gsup22() )
  end function h2

  function h3()
    double precision, dimension(0:n1, 0:n3) :: h3
    h3 = sqrt( g33() )
  end function h3

  ! ---------------------------------------------------------------------------
  ! ------------------------ Scale Factors for Mapping to Spherical Coordinates
  ! ---------------------------------------------------------------------------

  ! We always give hr at RI because Br .eq. 0 at RE by construction. 
  function hr()
    double precision, dimension(0:n1, 0:1) :: hr
    hr =  -2*RI*cos( q(:, zz) )**3 / ( cos( q(:, kk) ) * ( 1+3*cos( q(:, zz) )**2 ) )
  end function hr

  ! By default, compute the value at the ionosphere. 
  function hf(R0)
    double precision, dimension(0:n1, 0:1) :: hf
    double precision, optional, intent(in) :: R0
    if ( present(R0) ) then
      hf = R0*sin( q(:, kk) )
    else
      hf = RI*sin( q(:, kk) )
    end if
  end function hf

  ! By default, compute the value at the ionosphere. 
  function hq(R0)
    double precision, dimension(0:n1, 0:1) :: hq
    double precision, optional, intent(in) :: R0
    if ( present(R0) ) then
      hq = -R0**2 / ( RI * sin( 2*q(:,kk) ) )
    else
      hq = -RI / sin( 2*q(:,kk) )
    end if
  end function hq

  ! ===========================================================================
  ! ==================================================== End of Geometry Module
  ! ===========================================================================

end module geometry

! #############################################################################
! ######################################################### Ionospheric Profile
! #############################################################################

! The ionos module is a wrapper around the io module (which we might say is
! itself a wrapped around Fortran's read() and write() functions). It allows
! for convenient access to the input files describing ionospheric conductivity,
! etc. Values are mapped onto our grid spacing, as well as scaled to our units. 

module ionos
  use geometry
  use io
  implicit none

  ! ===========================================================================
  ! ======================================================== Parameter Profiles
  ! ===========================================================================

  ! Parallel, Hall, and Pedersen conductivity in mS/m. 
  double precision, allocatable, dimension(:,:) :: sig0, sigH, sigP
  ! Electric constants in mF/m. We use epsPara to apply the Boris
  ! approximation to the parallel electric constant. 
  double precision                              :: epsPara
  double precision, allocatable, dimension(:,:) :: epsPerp
  ! Integrated atmospheric conductivities in kS. 
  double precision                              :: intSig0, intSigH, intSigP

  contains

  ! ===========================================================================
  ! ========================================================== Ionosphere Setup
  ! ===========================================================================

  ! ---------------------------------------------------------------------------
  ! -------------------------------------------------------- File Access Helper
  ! ---------------------------------------------------------------------------

  function getIonos(varname)
    double precision, allocatable, dimension(:) :: getIonos
    character(len=*), intent(in)                :: varname
    allocate( getIonos( 0:lines(ionosfile) ) )
    select case (varname)
      case ('epsp')
        ! Electric constant in mF/m. Convert from units of eps0. 
        getIonos = eps0*readColumn(ionosfile, 1)
      case ('sigP')
        ! Hall conductivity in mS/m. Convert from cgs. 
        getIonos = 4*pi*eps0*readColumn(ionosfile, 2)
      case ('sigH')
        ! Pedersen conductivity in mS/m. Convert from cgs. 
        getIonos = 4*pi*eps0*readColumn(ionosfile, 3)
      case ('eta')
        ! Diffusivity in Mm**2/s. Convert from km**2/s. 
        getIonos = 1e-6*readColumn(ionosfile, 4)
      case ('sig0')
        ! Parallel conductivity in mS/m. Convert from diffusivity in km**2/s. 
        getIonos = 1e6/( mu0*readColumn(ionosfile, 4) )
      case default
        ! Radius in Mm. Convert from altitude in km. 
        getIonos = RE + 1e-3*readColumn(ionosfile, 0)
    end select
  end function getIonos

  ! ---------------------------------------------------------------------------
  ! -------------------------------------- Integrate Atmospheric Conductivities
  ! ---------------------------------------------------------------------------

  double precision function integrate(f, x, xmin, xmax)
    double precision, dimension(0:), intent(in) :: f, x
    double precision, intent(in)                :: xmin, xmax
    integer                                     :: i
    integrate = 0.
    ! For simplicity, we assume f = ( f(i) + f(i+1) )/2 from x(i) to x(i+1).
    do i=0,size(x)-2
      if (x(i) .le. xmin) cycle
      if (x(i) .ge. xmax) return
      integrate = integrate + ( min( xmax, x(i+1) ) - max( xmin, x(i) ) )     &
                              *( f(i) + f(i+1) )/2
    end do
  end function integrate

  subroutine atmSetup()
    ! Convert integrated atmospheric conductivities from Mho to kS. 
    intSig0 = 0.001*readParam('sig0atm')
    intSigH = 0.001*readParam('sighatm')
    intSigP = 0.001*readParam('sigpatm')
    ! If negative values are given, integrate the values ourselves. 
    if ( intSig0 .le. 0 ) intSig0 = integrate(getIonos('sig0'), getIonos('r'), RE, RI )
    if ( intSigH .le. 0 ) intSigH = integrate(getIonos('sigH'), getIonos('r'), RE, RI )
    if ( intSigP .le. 0 ) intSigP = integrate(getIonos('sigP'), getIonos('r'), RE, RI )

    call writeParam('Sigma0, integrated RE to RI', 1000*intSig0, 'Mho')
    call writeParam('SigmaP, integrated RE to RI', 1000*intSigP, 'Mho')
    call writeParam('SigmaH, integrated RE to RI', 1000*intSigH, 'Mho')

    call writeParam('Sigma0, integrated everywhere', 1000*integrate(getIonos('sig0'), getIonos('r'), RE, maxval( getIonos('r') ) ), 'Mho')
    call writeParam('SigmaP, integrated everywhere', 1000*integrate(getIonos('sigP'), getIonos('r'), RE, maxval( getIonos('r') ) ), 'Mho')
    call writeParam('SigmaH, integrated everywhere', 1000*integrate(getIonos('sigH'), getIonos('r'), RE, maxval( getIonos('r') ) ), 'Mho')

  end subroutine atmSetup

  ! ----------------------------------------------------------------------------------------------
  ! --------------------------------------------------------- Map Ionospheric Profiles to the Grid
  ! ----------------------------------------------------------------------------------------------

  function mapIonos(varname, power)
    double precision, dimension(0:n1, 0:n3)     :: mapIonos, w0, w1
    character(len=*), intent(in)                :: varname
    integer, intent(in), optional               :: power
    integer                                     :: line, nLines
    double precision, allocatable, dimension(:) :: rIon, vIon
    ! Grab the ionos values for radius and whichever var we'relooking at. 
    nLines = lines(ionosfile)
    allocate( rIon(0:nLines), vIon(0:nLines) )
    rIon = getIonos('r')
    vIon = abs( getIonos(varname) )
    ! Above the ionospheric profile, the values follow a power law or just go to 0. Well, not
    ! quite zero, since we don't want to deal with infinity when we invert. 
    mapIonos = 1e-20
    if ( present(power) ) mapIonos = vIon(nLines) * ( r/rIon(nLines) )**power
    ! For each line of the ionos profile, update the grid points that fall at that radius. 
    do line=1,nLines
      where ( r .gt. rIon(line-1) .and. r .le. rIon(line) )
        ! Values between rIon radii are interpolated linearly. 
        w0 = rIon(line) - r
        w1 = r - rIon(line-1)
        mapIonos = ( w0*vIon(line-1) + w1*vIon(line) )/(w0 + w1)
      end where
    end do
  end function mapIonos

  subroutine profileSetup()
    allocate( epsPerp(0:n1, 0:n3), sig0(0:n1, 0:n3), sigH(0:n1, 0:n3), sigP(0:n1, 0:n3) )
    ! The perpendicular electric constant comes from (isotropic) tabulated values, then we add a
    ! term for the plasmaspheric density contribution. The parallel electric constant is handled
    ! with the time step. 
    epsPerp = mapIonos('epsp', power=5) + mp*mmm()*n()/BB()
    ! The mapIonos function drops profiles to zero at large distances, but field-aligned
    ! conuctivity should actually go to infinity. We handle that by letting eta go to zero, then
    ! taking its inverse. 
    sig0 = 1/( mu0*mapIonos('eta') )
    sigH = readParam('sighfac')*mapIonos('sigH')
    sigP = readParam('sigpfac')*mapIonos('sigP')
  end subroutine profileSetup

  ! ----------------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------- Compute Time Step
  ! ----------------------------------------------------------------------------------------------

  ! This routine computes the time step according to the Alfven speed and grid spacing, the time
  ! step based on parallel plasma motions, adjusts them based on any Boris corrections, and
  ! chooses a time step for the experiment. 
  subroutine dtSetup()
    ! One over the zone crossing times in each direction. 
    double precision                       :: oodtx, oodty, oodtz
    ! The Alfven time step is constrained by zone crossing times. The inertial time step is
    ! constrained by the plasma frequency. The compressional time step is constrained by
    ! speed-of-light crossing time. 
    double precision                       :: dtAlfven, dtInertial, dtCompressional
    ! We can use a Boris correction to increase eps0 above its vacuum value. This increases the
    ! plasma frequency, allowing us to stabilize parallel and inertial terms at a reasonable time
    ! step. This is safe if the plasma frequency remains much higher than the waves we're driving.
    double precision                       :: boris, borisMax
    double precision                       :: c, wdrive
    double precision, dimension(0:n1,0:n3) :: wp, np
    double precision, dimension(0:n1,0:n3) :: dx, dy, dz
    ! Compute the grid spacing. We'll need it when we get to crossing times. 
    dx = 0.5*h1()*( eoshift(usup1(), 1, dim=1) - eoshift(usup1(), -1, dim=1) )
    dy = h2()/azm
    dz = -0.5*h3()*( eoshift(usup3(), 1, dim=2) - eoshift(usup3(), -1, dim=2) )
    ! We used eoshift, so the edges are garbage. Set them to "infinity." 
    dx(ii,:) = 1e10
    dz(:,kk) = 1e10
    ! Let's print out some characteristic scales for our grid. 
    call writeParam('Ionospheric Min dx', 1000*minval( dx(:, 0) ), 'km')
    call writeParam('Ionospheric Max dx', 1000*maxval( dx(1:n1-1, 0) ), 'km')
    call writeParam('Equatorial Min dx', 1000*minval( dx(:, n3/2) ), 'km')
    call writeParam('Equatorial Max dx', 1000*maxval( dx(1:n1-1, n3/2) ), 'km')
    call writeParam('Outermost Min dz', 1000*minval( dz(n1, :) ), 'km')
    call writeParam('Outermost Max dz', 1000*maxval( dz(n1, 1:n3-1) ), 'km')
    call writeParam('Innermost Min dz', 1000*minval( dz(0, :) ), 'km')
    call writeParam('Innermost Max dz', 1000*maxval( dz(0, 1:n3-1) ), 'km')
    ! Our time step is constrained by the zone-crossing time of an Alfven wave. 
    oodtx = maxval(vA()/dx)
    oodty = maxval(vA()/dy)
    oodtz = maxval(vA()/dz)
    call writeParam('Min dx/vA', 1/oodtx, 's')
    call writeParam('Min dy/vA', 1/oodty, 's')
    call writeParam('Min dz/vA', 1/oodtz, 's')
    ! Get diagonal zone crossing -- still working in one over time. 
    if (azm .eq. 0) then
      dtAlfven = readParam('cour')/ ( sqrt(2.)*sqrt( oodtx**2 + oodtz**2 ) )
    else
      dtAlfven = readParam('cour')/ ( sqrt(3.)*sqrt( oodtx**2 + oodty**2 + oodtz**2 ) )
    end if
    call writeParam('Alfven dt', dtAlfven, 's')
    ! Without electron inertial effects, the Alfven speed sets the time step, and we don't change
    ! the parallel dielectric constant. This means that there are no parallel electric fields. 
    if ( readParam('inertia') .le. 0) then
      dt = dtAlfven
      epsPara = eps0
    ! If we're using electron inertial effects, we also need to consider
    ! compressional crossing time and the plasma frequency as time step
    ! constraints. We lower these with a Boris correction: we raise the
    ! electric constant, which lowers the plasma frequency and speed of light,
    ! while keeping the electron inertial length unchanged.  
    else
      ! Our Boris factor is constrained by wdrive**2/wp**2 << 1 and wdrive*np/wp**2 << 1. Let's keep
      ! four orders of magnitude, just to be safe. Note that wp and np are the plasma frequency and
      ! the parallel collision frequency, respectively. 
      wp = sqrt( n()*qe**2 / (me*eps0) )
      np = n()*qe**2 / (me*sig0)
      wdrive = 2*pi*readParam('fdrive')
      ! Note that wp**2 scales as 1/eps0. 
      call writeParam( 'Max w**2/wp**2', maxval(wdrive/wp)**2 )
      call writeParam( 'Max w*nu/wp**2', maxval(wdrive*np/(wp)**2) )
      borisMax = 1e-4 / max( maxval( wdrive**2 / wp**2 ), maxval( wdrive*np / wp**2 ) )
      call writeParam('Maximum Boris Factor', borisMax)
      ! If we're given a Boris factor, use that. 
      if (readParam('epsfac') .gt. 0) then
        boris = readParam('epsfac')
      ! Otherwise, use the biggest one we can get away with. 
      else
        boris = borisMax
      end if
      call writeParam('Boris Factor', boris)
      ! Recompute the plasma frequency with the Boris factor. Also compute adjusted speed of light. 
      epsPara = eps0*boris
      wp = sqrt( n()*qe**2 / (me*epsPara) )
      c = sqrt( 1/(mu0*epsPara) )
      ! Confirm that we're not worried about stability. 
      call writeParam( 'Max Boris-Adjusted w**2/wp**2', maxval(wdrive/wp)**2 )
      call writeParam( 'Max Boris-Adjusted w*nu/wp**2', maxval(wdrive*np/(wp)**2) )
      ! Now figure out the time step. We are constrained by Alfven zone crossing time, compressional
      ! (speed of light) zone crossing time, and the plasma frequency. Let's compute the grid spacing
      ! everywhere. 
      call writeParam('Min Boris-Adjusted 1/wp', minval(1/wp), 's')
      ! The inertial effects are particularly sensitive to the size of the time step. We use a fudge
      ! factor (in addition to the Courant condition) to ensure stability. 
      dtInertial = readParam('cour')*readParam('fudge')*minval(1/wp)
      call writeParam('Inertial dt', dtInertial, 's')
      ! The zone crossing time of the compressional mode constrains our time step. 
      oodtx = maxval(c/dx)
      oodty = maxval(c/dy)
      oodtz = maxval(c/dz)
      call writeParam('Min Boris-Adjusted dx/c', 1/oodtx, 's')
      call writeParam('Min Boris-Adjusted dy/c', 1/oodty, 's')
      call writeParam('Min Boris-Adjusted dz/c', 1/oodtz, 's')
      ! Get diagonal zone crossing -- still working in one over time. 
      if (azm .eq. 0) then
        dtCompressional = readParam('cour')/ ( sqrt(2.)*sqrt( oodtx**2 + oodtz**2 ) )
      else
        dtCompressional = readParam('cour')/ ( sqrt(3.)*sqrt( oodtx**2 + oodty**2 + oodtz**2 ) )
      end if
      call writeParam('Compressional dt', dtCompressional, 's')
      ! The time step is based on the Alfven crossing time, the compressional crossing time, or
      ! the plasma timescale, whichever is strictest. 
      dt = min(dtAlfven, dtInertial, dtCompressional)
    end if
    ! Write out the time step. Let's also get a sense for the relative cost of the grid setup. 
    call writeParam('dt', dt, 's')
    call writeParam('Floating Point Operations', n1*n3/dt, 'per second of time simulated')
  end subroutine dtSetup

  ! ----------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------- Write Ionospheric Profiles
  ! ----------------------------------------------------------------------------------------------

  subroutine writeIonos()
    double precision, dimension(0:n1, 0:n3) :: dz, v
    ! Conductivity profiles are printed in S/m (Mho/m). 
    call writeRealArray('sig0.dat', 1000*sig0)
    call writeRealArray('sigH.dat', 1000*sigH)
    call writeRealArray('sigP.dat', 1000*sigP)
    ! Perpendicular electric constant in units of eps0. 
    call writeRealArray('epsPerp.dat', epsPerp/eps0)
    ! Number density in cm^-3. 
    call writeRealArray( 'n.dat', 1e-24*n() )
    ! We don't print out profiles for the plasma frequency and Alfven speed, since those can be
    ! easily computed from the parameters we do print out. Instead, we give snapshots. 
    call writeParam('Alfven Speed Min', 1000*minval( vA() ), 'km/s')
    call writeParam('Alfven Speed Max', 1000*maxval( vA() ), 'km/s')

    call writeParam('Electron Inertial Length Min', 1000*minval( sqrt( me / ( mu0*n()*qe**2 ) ) ), 'km')
    call writeParam('Electron Inertial Length Max', 1000*maxval( sqrt( me / ( mu0*n()*qe**2 ) ) ), 'km')

    call writeParam('Plasma Frequency Min', minval( sqrt( n()*qe**2 / (me*eps0) ) ), 'rad/s')
    call writeParam('Plasma Frequency Max', maxval( sqrt( n()*qe**2 / (me*eps0) ) ), 'rad/s')

    call writeParam('Adjusted Plasma Frequency Min', minval( sqrt( n()*qe**2 / (me*epsPara) ) ), 'rad/s')
    call writeParam('Adjusted Plasma Frequency Max', maxval( sqrt( n()*qe**2 / (me*epsPara) ) ), 'rad/s')

    call writeParam('epsPerp/sigP Min', minval(epsPerp/sigP), 's')
    call writeParam('epsPerp/sigH Min', minval(epsPerp/sigH), 's')
    call writeParam('eps0/sig0    Max', maxval(eps0/sig0), 's')

    call writeParam('vA*epsPerp/sigP Min', 1000*minval(vA()*epsPerp/sigP), 'km')
    call writeParam('vA*epsPerp/sigH Min', 1000*minval(vA()*epsPerp/sigH), 'km')

    call writeParam('epsPara/sig0 Max', maxval(epsPerp/sig0), 's')
    call writeParam('Collision Frequency Max', maxval( n()*qe**2 / (me*sig0) ), 'Hz')
    call writeParam('Adjusted Speed of Light', 1000*sqrt( 1./(mu0*epsPara) ), 'km/s')
!    ! Alfven speed in km/s. 
!    call writeRealArray('vA.out', 1000*vA())
!    ! Plasma frequencies for electrons and protons in radians/s. 
!    call writeRealArray( 'wpe.out', sqrt( n()*qe**2 / (me*epsPara) ) )
!    call writeRealArray( 'wpp.out', sqrt( n()*qe**2 / (mp*epsPara) ) )
!    ! Gyrofrequencies for electrons and protons in radians/s. 
!    call writeRealArray('wce.out', qe*sqrt( BB() )/me)
!    call writeRealArray('wcp.out', qe*sqrt( BB() )/mp)
    ! Bounce time (factor of 2 for there and back) in seconds. Compute the distance along the
    ! field line between each pair of points, divide by the average Alfven speed, and sum. 
    dz = abs( eoshift(usup3()*h3(), 1, dim=2) - usup3()*h3() )
    v = 0.5*( eoshift(vA(), 1, dim=2) + vA() )
    call writeRealColumns( 'L.dat', L_() )
    call writeRealColumns( 'tBounce.dat', 2*sum( dz(:, 0:n3-1)/v(:, 0:n3-1), dim=2 ) )
  end subroutine writeIonos

  ! ----------------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------- Ionosphere Setup Driver
  ! ----------------------------------------------------------------------------------------------

  subroutine ionosSetup()
    ! Alter the ionosfile path to get the model we want. 
    write(ionosfile(len(ionosfile) - 4:len(ionosfile) - 4),'(i1)') int( readParam('model') )
    ! Integrate atmospheric conductivities. 
    call atmSetup()
    ! Compute grid-resolved electric constant and conductivities. 
    call profileSetup()
    ! Figure out the time step based on the now-known Alfven speed. 
    call dtSetup()
    ! Write out ionos profiles lined up with our grid. 
    call writeIonos()
  end subroutine ionosSetup

  ! ==============================================================================================
  ! =================================================================== Parameter Profile "Arrays"
  ! ==============================================================================================

  ! ----------------------------------------------------------------------------------------------
  ! --------------------------------------------------------- Zeroth-Order Field Strength, Squared
  ! ----------------------------------------------------------------------------------------------

  function BB()
    double precision, dimension(0:n1, 0:n3) :: BB
    ! At the equator, Earth's magnetic field is 31.1e3 nT. We scale it as an ideal dipole. 
    BB = (31.1e3)**2 * (RE/r)**6 * abs( 1+3*cos(q)**2 )
  end function BB

  ! ----------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------- Number Density
  ! ----------------------------------------------------------------------------------------------

  function n()
    double precision, dimension(0:n1, 0:n3) :: n
    ! This profile depends on several parameters, which we don't use for anything else. As long as
    ! we only call n() a few times, we are fine not holding on to them. 
    double precision                        :: haz, naz, nps, Lps, Lpp, dLpp
    ! Convert number density from cm**-3 to Mm**-3 at the base of the auroral zone, plasmasphere.
    naz = 1e24 * readParam('naz')
    nps = 1e24 * readParam('nps')
    ! Convert auroral zone scale height from RE to Mm. 
    haz = RE * readParam('haz')
    ! Plasmasphere scale L. 
    Lps = readParam('lps')
    ! Plasmapause location and width, in L. 
    Lpp = readParam('lpp') ! plasmapause mcilwain parameter
    dLpp = readParam('dlpp') ! plasmapause spread in lpp
    ! The number density is the sum of the auroral and plasmaspheric profiles. 
    n = naz*haz/r + 0.5*nps*exp( -L()/Lps )*( 1-tanh( (L()-Lpp)/dLpp ) )
  end function n

  ! ----------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------- Mean Molecular Mass
  ! ----------------------------------------------------------------------------------------------

  ! Right at the atmospheric boundary, there's a ton of oxygen, but at high altitude it's all
  ! protons. This profile is eyeballed from Lysak [2013], which credits Kelley [1989]. This factor
  ! has practically no effect on the perpendicular electric constant (and thus the Alfven speed). 
  ! Near the atmospheric boundary, the magnetic field is small enough that the density term has 
  ! little effect. And the exponential is tiny by the time it gets to IAR altitudes. 
  function mmm()
    double precision, dimension(0:n1, 0:n3) :: mmm
    ! At the atmosphere, the mean molecular mass is 28 amu. It decays to 1 amu by about 1300km. 
    mmm = 1 + 27*exp( -(r - RI)/0.7 )
  end function mmm

  ! ----------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------- Alfven Speed
  ! ----------------------------------------------------------------------------------------------

  function vA()
    double precision, dimension(0:n1, 0:n3) :: vA
    ! This formulation of the Alfven speed is relativity-safe. 
    vA = sqrt( cc / (epsPerp/eps0) )
  end function vA

  ! ==============================================================================================
  ! ============================================================ End of Ionospheric Profile Module
  ! ==============================================================================================

end module ionos

! ################################################################################################
! ################################################################################### Coefficients
! ################################################################################################

! Extensive precomputation of coefficients makes Tuna's main loop more legible, while also
! increasing its speed. We know that every time step we'll update B1 from JCsup1 by applying a
! factor of -g11*dt/J. Multiplying those together here lets us accomplish that in one floating
! point operation during the main loop instead of four. The coefficients all follow Lysak's notes;
! however, they do not follow his code exactly. Lysak takes derivatives of covariant field 
! components, uses them to update contravariant components, then maps back to the covariant basis
! using the metric tensor. Tuna, instead, wraps the metric tensor components into the coefficients
! and stores only covariant fields explicitly. The coefficients module handles the assembly of
! these coefficients from the ionospheric profiles and geometric factors. 

module coefficients
  use loop
  use geometry
  use ionos
  implicit none

  contains

  ! ==============================================================================================
  ! ========================================================================= Coefficient Assembly
  ! ==============================================================================================

  ! ----------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------- Assemble Bulk Coefficients
  ! ----------------------------------------------------------------------------------------------

  subroutine bulkCoefficientSetup()
    ! We define a handful of intermediate factors for convenience. 
    double precision, dimension(0:n1,0:n3) :: sine, sinf, cose, cosf, expe, expf, gg12, gg21
    sine = sin(        sigh * dt / epsPerp )
    sinf = sin(  0.5 * sigh * dt / epsPerp )
    cose = cos(        sigh * dt / epsPerp )
    cosf = cos(  0.5 * sigh * dt / epsPerp )
    expe = exp(       -sigp * dt / epsPerp )
    expf = exp( -0.5 * sigp * dt / epsPerp )
    ! B1 coefficients.
    B1_JCsup1 = -g11()*dt/J()
    B1_JCsup3 = -g13()*dt/J()
    ! B2 coefficients.
    B2_JCsup2 = -g22()*dt/J()
    ! B3 coefficients.
    B3_JCsup1 = -g31()*dt/J()
    B3_JCsup3 = -g33()*dt/J()
    ! Only use E3 and j3 if we're worrying about electron inertial effects. 
    if (readParam('inertia') .gt. 0) then
      ! j3 coefficients (solved with integrating factors).
      j3_j3 = exp( -dt*n()*qe**2 / (me*sig0) )
      j3_E3 = dt*n()*qe**2/me * exp( -dt*n()*qe**2 / (2*me*sig0) )
      ! E3 coefficients (advanced directly). 
      E3_E3 = 1
      E3_JFsup1 = dt*g31() / ( mu0*epsPara*J() )
      E3_JFsup3 = dt*g33() / ( mu0*epsPara*J() )
      E3_j3 = -dt/epsPara
    else
      ! j3 coefficients (without electron inertia, we don't have j3). 
      j3_j3 = 0
      j3_E3 = 0
      ! E3 coefficients (solved with integrating factors like E1 and E2). 
      E3_E3 = exp(-sig0*dt/epsPara)
      E3_JFsup1 = dt * exp( -sig0*dt/(2*epsPara) ) * g31() / ( mu0*epsPara*J() )
      E3_JFsup3 = dt * exp( -sig0*dt/(2*epsPara) ) * g33() / ( mu0*epsPara*J() )
      E3_j3 = 0
    end if
    ! E1 coefficients. Note that some E1 coefficients are defined in terms of E3 coefficients. 
    ! This is because we have an equation for updating Esup1, and an equation for updating E3, and
    ! we have to match those together to get our expressions for updating E1. 
    E1_E1 =                                   cose*expe
    E1_E2 =       sqrt( gsup22() / gsup11() )*sine*expe
    E1_E3 =           ( gsup13() / gsup11() )*expe*cose                                          &
                    - ( gsup13() / gsup11() )*E3_E3
    E1_JFsup1 =       1/( J()*gsup11() )     *cosf*expf*dt*vA()**2                               &
                    - ( gsup13() / gsup11() )*E3_JFsup1
    E1_JFsup2 =   sqrt( 1 / g33() )          *sinf*expf*dt*vA()**2
    E1_JFsup3 =      -( gsup13() / gsup11() )*E3_JFsup3
    E1_j2drive =  sqrt( gsup22() / gsup11() )*sinf*expf*dt/epsPerp
    E1_j3 =          -( gsup13() / gsup11() )*E3_j3
    ! E2 coefficients. 
    E2_E1 =      -sqrt( gsup11() / gsup22() )*sine*expe
    E2_E2 =                                   cose*expe
    E2_E3 = -gsup13() / sqrt( gsup11()*gsup22() )*expe*sine
    E2_JFsup1 =  -sqrt( 1 / g33() )          *sinf*expf*dt*vA()**2
    E2_JFsup2 =       1/( J()*gsup22() )     *cosf*expf*dt*vA()**2
    E2_j2drive =                              cosf*expf*dt/epsPerp
  end subroutine bulkCoefficientSetup

  ! ----------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------- Assemble Edge Coefficients
  ! ----------------------------------------------------------------------------------------------

  subroutine edgeCoefficientSetup()
    integer                                :: i, ip
    ! A bunch of temporary variables for legibility. 
    double precision, dimension(0:n1)      :: npo, tnpo, nuface, nufaci
    double precision, dimension(0:n1,0:n3) :: scratch
    double precision, dimension(0:n1, 0:1) :: Br_B1, Br_B3
    double complex, dimension(0:n1, 0:n1)  :: PsiI_Br
    double precision, dimension(0:n1, 0:1) :: cosa, sina, intSig11, intSig12, intSig21,          &
                                              intSig22, intSigzz
    ! For legibility, coefficients for Psi are assembled in two steps. First, we compute the
    ! coefficients to get from B1 and B3 to Br at the ionosphere. Then, we compute coefficients to
    ! get from Br to Psi. We slam those together into a single set of coefficients, then, since we
    ! don't store Br (proportional to Bsup3) in memory during runtime. 
    scratch = gsup31()
    Br_B1 = hr()*scratch(:, kk)
    scratch = gsup33()
    Br_B3 = hr()*scratch(:, kk)
    ! Going from Br to Psi is a bit more complicated, so we use some further intermediates. 
    npo = nu + 1
    tnpo = 2*nu + 1
    nuface = (RI/nu)*(tnpo/npo)*(RE/RI)**nu / ( 1 - (RE/RI)**tnpo )
    nufaci = (RI/nu)*( 1 + (nu/npo)*(RE/RI)**tnpo ) / ( 1 - (RE/RI)**tnpo )
    do i=0,n1
      do ip=0,n1
        PsiE_Br(i, ip) = sum( Y(i, :)*Yinv(:, ip)*nuface(:) )
        PsiI_Br(i, ip) = sum( Y(i, :)*Yinv(:, ip)*nufaci(:) )
      end do
    end do
    ! We want to be able to output the harmonic weights, alpha_nu. They are used as:
    ! Psi = Sum_nu alpha_nu Y_nu ( r^nu + nu/(nu+1) RE^(2nu+1) r^(-nu-1) )
    ! This comes from the definition of Psi in terms of harmonics:
    ! Psi = Sum_nu Y_nu (  alpha_nu r^nu + beta_nu r^(-nu-1) )
    ! Combined with the fact that Br = dPsi/dr is zero at RE. 
    allocate( alpha(0:nModes-1, 0:1), alpha_Br(0:nModes-1, 0:n1) )
    ! The form for alpha_Br comes from working out the expression for Psi in terms of Br, just
    ! like we worked it out for PsI_Br. Here, though, we only have to do half the work, since we
    ! stop at the harmonic weights (instead of summing them up again to get Psi). 
    do i=0,n1
      alpha_Br(:, i) = Yinv(i, 0:nModes-1)                                                       &
                     / ( nu(0:nModes-1)*RI**(nu(0:nModes-1)-1)                                   &
                         *(1 - (RE/RI)**tnpo(0:nModes-1) ) )
    end do
    ! Note that we don't have to worry about the second step for alpha or PsiE. We don't compute
    ! those during runtime -- only when we output fields -- so we can rely on the (relatively
    ! slow) Br() helper function. 
      do i=0,n1
        do ip=0,n1
          PsiI_B1(i, ip, hh) = PsiI_Br(i, ip)*Br_B1(ip, hh)
          PsiI_B3(i, ip, hh) = PsiI_Br(i, ip)*Br_B3(ip, hh)
        end do
      end do
    ! assemble coefficients for computing boundary e from bi
    cosa = -2*cos( q(:, kk) )/sqrt( 1 + 3*cos( q(:, kk) )**2 )
    sina = sqrt( 1 - cosa**2 )
    ! Compute the integrated conductivity tensor for getting boundary electric fields from the
    ! magnetic field jump over the ionospheric current sheet. 
    intSigzz = intSig0*cosa**2 + intSigP*sina**2
    intSig11 = ( hf()/hq() ) * intSig0*intSigP/intSigzz
    intSig12 = -intSig0*intSigH*cosa/intSigzz
    intSig21 = -intSig12
    intSig22 = ( hq()/hf() ) *( intSigP + sina**2*intSigH**2 / intSigzz )
    ! Assemble coefficients for computing edge electric fields from edge magnetic fields. 
    E1_B1I = -intSig12/( mu0*( intSig11*intSig22 - intSig12*intSig21 ) )
    E1_B2I = -intSig22/( mu0*( intSig11*intSig22 - intSig12*intSig21 ) )
    E2_B1I =  intSig11/( mu0*( intSig11*intSig22 - intSig12*intSig21 ) )
    E2_B2I =  intSig21/( mu0*( intSig11*intSig22 - intSig12*intSig21 ) )
  end subroutine edgeCoefficientSetup

  ! ----------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------- Coefficient Setup Driver
  ! ----------------------------------------------------------------------------------------------

  subroutine coefficientSetup()

    allocate( B1_JCsup1(0:n1, 0:n3), B1_JCsup3(0:n1, 0:n3), B2_JCsup2(0:n1, 0:n3),               &
              B3_JCsup1(0:n1, 0:n3), B3_JCsup3(0:n1, 0:n3), E1_E1(0:n1, 0:n3),                   &
              E1_E2(0:n1, 0:n3), E1_E3(0:n1, 0:n3), E1_JFsup1(0:n1, 0:n3),                       &
              E1_JFsup2(0:n1, 0:n3), E1_JFsup3(0:n1, 0:n3), E1_j2drive(0:n1, 0:n3),              &
              E1_j3(0:n1, 0:n3), E2_E1(0:n1, 0:n3), E2_E2(0:n1, 0:n3), E2_E3(0:n1, 0:n3),        &
              E2_JFsup1(0:n1, 0:n3), E2_JFsup2(0:n1, 0:n3), E2_j2drive(0:n1, 0:n3),              &
              E3_E3(0:n1, 0:n3), E3_j3(0:n1, 0:n3), E3_JFsup1(0:n1, 0:n3),                       &
              E3_JFsup3(0:n1, 0:n3), j3_j3(0:n1, 0:n3), j3_E3(0:n1, 0:n3) )

  allocate( E1_B1I(0:n1, 0:1), E1_B2I(0:n1, 0:1), E2_B1I(0:n1, 0:1), E2_B2I(0:n1, 0:1),          &
            PsiE_Br(0:n1, 0:n1), PsiI_B1(0:n1, 0:n1, 0:1), PsiI_B3(0:n1, 0:n1, 0:1) )

    call bulkCoefficientSetup()
    call edgeCoefficientSetup()
  end subroutine coefficientSetup

  ! ==============================================================================================
  ! =================================================================== End of Coefficients Module
  ! ==============================================================================================

end module coefficients

! ################################################################################################
! ######################################################################################### Fields
! ################################################################################################

! The fields module assembles the fields which are advanced in the main loop. It also includes helper functions for mapping between covariant fields (which we store) and physical fields (which we output). Subroutines for interpolating and printing field values are also stored in the fields module. 

module fields
  use loop
  use coefficients
  use geometry
  use ionos
  implicit none

  contains

  ! ==============================================================================================
  ! ================================================================================== Field Setup
  ! ==============================================================================================

  ! ----------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------- Write Field Data Headers
  ! ----------------------------------------------------------------------------------------------

  subroutine writeFieldHeaders()
    integer :: nt
    nt = int(tmax/dtout)
    ! Time stamps. 
    call writeInteger('t.dat', nt)
    call writeInteger('driveScale.dat', nt)
    ! Bulk magnetic fields mapped to physical coordinates in nT. 
    call writeComplexArray('Bx.dat', Bx(), onlyheader=.true., nt=nt)
    call writeComplexArray('By.dat', By(), onlyheader=.true., nt=nt)
    call writeComplexArray('Bz.dat', Bz(), onlyheader=.true., nt=nt)
    ! Bulk electric fields mapped to physical coordinates in mV/m. 
    call writeComplexArray('Ex.dat', Ex(), onlyheader=.true., nt=nt)
    call writeComplexArray('Ey.dat', Ey(), onlyheader=.true., nt=nt)
    call writeComplexArray('Ez.dat', Ez(), onlyheader=.true., nt=nt)
    ! Parallel and driving current mapped to physical coordinates in uA/m^2. 
    call writeComplexArray('Jz.dat', jz(), onlyheader=.true., nt=nt)
    call writeRealArray( 'JyDrive.dat', j2drive/h2() )
    ! Scalar magnetic potential in... nT*Mm? This may not be a physically meaningful quantity. 
    call writeComplexArray('PsiE.dat', PsiE, onlyheader=.true., nt=nt)
    call writeComplexArray('PsiI.dat', PsiI, onlyheader=.true., nt=nt)
    ! Edge magnetic fields in nT. 
    call writeComplexArray('BfE.dat', Bf(RE), onlyheader=.true., nt=nt)
    call writeComplexArray('BfI.dat', Bf(RI), onlyheader=.true., nt=nt)
    call writeComplexArray('BqE.dat', Bq(RE), onlyheader=.true., nt=nt)
    call writeComplexArray('BqI.dat', Bq(RI), onlyheader=.true., nt=nt)
    call writeComplexArray('Br.dat', Br(), onlyheader=.true., nt=nt)
    ! Harmonic weights. One per harmonic, per hemisphere, per time step. 
    call writeComplexArray('eweights.dat', alpha, onlyheader=.true., noskip=.true., nt=nt)
  end subroutine writeFieldHeaders

  ! ----------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------- Driving Field Setup
  ! ----------------------------------------------------------------------------------------------

  subroutine driveSetup()
    double precision                       :: qdrive, dqdrive, Ldrive, dLdrive, Bdrive, jdrive
    double precision, dimension(0:n1,0:n3) :: scratch
    integer                                :: i
    ! Let's zero all driving that would be delivered inside the ionosphere. That's just asking for
    ! instabilities. 
    double precision                       :: rmin
    ! Grab driving parameters: waveform index, frequency, characteristic timescale. 
    idrive = readParam('idrive')
    wdrive = 2*pi*readParam('fdrive')
    tdrive = readParam('tdrive')
    ! Get the latitude and spread in latitude, and convert from degrees to radians. 
    qdrive = (pi/180)*( 90 - readParam('latdrive') )
    dqdrive = (pi/180)*readParam('dlatdrive')
    ! Current driving also needs to be delivered at a radial distance. 
    Ldrive = readParam('ldrive')
    dLdrive = readParam('dldrive')
    ! Get the magnitude of the current and compressional driving. One of these should be zero. 
    Bdrive = readParam('bdrive')
    jdrive = readParam('jdrive')
    ! Compressional driving is delivered at the outer boundary. Map from Bz to B3. 
    scratch = h3()
    B3drive = Bdrive*scratch(n1, :)*exp( -( ( q(n1, :) - qdrive ) / dqdrive )**2 )
    ! Current driving is delivered through the electric field. Like the compressional driving, 
    ! it's gaussian in latitude, and it's also gaussian radial distribution. 
    j2drive = jdrive*exp( -0.5*( (q - qdrive)/dqdrive )**2 )*                                    &
              exp( -0.5*( (L() - Ldrive)/dLdrive )**2 )/( h2()*gsup22() )
    ! Find the maximum radius in the ionospheric profile, and the radius corresponding to the
    ! inner boundary. Don't allow any driving lower than there. 
    rmin = max( readParam('lmin')*RE, maxval( getIonos('r') ) )
    where (r(n1, :) .le. rmin) B3drive = 0.
    where (r .le. rmin) j2drive = 0.
!    write(*,*) 'max j2drive = ', maxval(j2drive)
!    write(*,*) 'max B3drive = ', maxval(B3drive)
    ! If we're driving with a spectrum, set up an ensemble of frequencies and phase offsets. 
    if (idrive == 4) then
      call random_seed()
      do i=1,nspectrum
        wdrives(i) = i*wdrive/nspectrum
        call random_number( pdrives(i) )
      end do
    end if
  end subroutine driveSetup

  ! ----------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------- Field Setup Driver
  ! ----------------------------------------------------------------------------------------------

  subroutine fieldSetup
    ! Allocate arrays to hold the fields. 
    allocate( B1(0:n1, 0:n3), B2(0:n1, 0:n3), B3(0:n1, 0:n3), E1(0:n1, 0:n3), E2(0:n1, 0:n3),    &
              E3(0:n1, 0:n3), JCsup1(0:n1, 0:n3), JCsup2(0:n1, 0:n3), JCsup3(0:n1, 0:n3),        &
              JFsup1(0:n1, 0:n3), JFsup2(0:n1, 0:n3), JFsup3(0:n1, 0:n3), B1E(0:n1, 0:1),        &
              B1I(0:n1, 0:1), B2E(0:n1, 0:1), B2I(0:n1, 0:1), PsiE(0:n1, 0:1), PsiI(0:n1, 0:1),  &
              B3drive(0:n3), j2drive(0:n1, 0:n3), j3(0:n1, 0:n3) )

    ! Set up the driving fields. 
    call driveSetup()
    ! Initialize bulk fields to zero. Most compilers do this on allocation, but it's not 
    ! guaranteed. Note that we don't have to initialize edge fields or curl components, since they
    ! are computed from scratch each step (rather than changing incrementally). 
    B1 = 0
    B2 = 0
    B3 = 0
    E1 = 0
    E2 = 0
    E3 = 0
    j3 = 0
    ! Write headers onto the field files in preparation of getting data. 
    call writeFieldHeaders()
  end subroutine fieldSetup

  ! ==============================================================================================
  ! ========================================================================= Runtime Field Output
  ! ==============================================================================================

  ! ----------------------------------------------------------------------------------------------
  ! -------------------------------------------------------- Interpolate a Field at Even Locations
  ! ----------------------------------------------------------------------------------------------

  ! Get field values (in place) at even i and even k for plotting. Use Neumann boundary conditions
  ! for all fields here, since we want to be able to slice at the edge for plotting. 
  subroutine interpolateEvens(field, i0, k0)
    double complex, intent(inout), dimension(0:n1, 0:n3) :: field
    integer, intent(in)                                  :: i0, k0
    integer                                              :: i, k
    ! If the field isn't defined on even i, interpolate in i. 
    if (i0 .ne. 0) then
      do k=k0,n3,2
        do i=0,n1,2
          field(i, k) = i1N(field, i, k)
        end do
      end do
    end if
    ! If the field isn't defined on even k, interpolate in k. 
    if (k0 .ne. 0) then
      do k=0,n3,2
        do i=0,n1,2
          field(i, k) = i3N(field, i, k)
        end do
      end do
    end if
  end subroutine interpolateEvens

  ! ----------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------- Interpolate a Field Everywhere
  ! ----------------------------------------------------------------------------------------------

  ! This isn't used for plotting, but rather for debugging. Sometimes we want to be able to look
  ! at any slice of any field. Interpolation happens in place. 
  subroutine interpolateEverywhere(field, i0, k0)
    double complex, intent(inout), dimension(0:n1, 0:n3) :: field
    integer, intent(in)                                  :: i0, k0
    integer                                              :: i, k
    ! Note: flip(0) .eq. 1 and flip(1) .eq. 0. 
    integer, dimension(0:1)                       :: flip = [1, 0]
    ! Interpolate to have both parities in i. Note that abs(i0-1) swaps between 1 and 0. 
    do k=k0,n3,2
      do i=flip(i0),n1,2
        field(i, k) = i1N(field, i, k)
      end do
    end do
    ! Interpolate to have both parities in k. 
    do k=flip(k0),n3,2
      do i=0,n1
        field(i, k) = i3N(field, i, k)
      end do
    end do
  end subroutine interpolateEverywhere

  ! ----------------------------------------------------------------------------------------------
  ! ----------------------------------------------------- Interpolate all Fields at Even Locations
  ! ----------------------------------------------------------------------------------------------

  ! Get field values (in place) at even i and even k for plotting. Use Neumann boundary conditions
  ! for all fields here, since we want to be able to slice at the edge for plotting. 
  subroutine interpolateEvenFields
    call interpolateEvens(B1, 1, 1)
    call interpolateEvens(JCsup1, 1, 1)
    call interpolateEvens(B2, 0, 1)
    call interpolateEvens(JCsup2, 0, 1)
    call interpolateEvens(B3, 0, 0)
    call interpolateEvens(JCsup3, 0, 0)
    call interpolateEvens(E1, 0, 0)
    call interpolateEvens(JFsup1, 0, 0)
    call interpolateEvens(E2, 1, 0)
    call interpolateEvens(JFsup2, 1, 0)
    call interpolateEvens(E3, 1, 1)
    call interpolateEvens(JFsup3, 1, 1)
    call interpolateEvens(j3, 1, 1)
  end subroutine interpolateEvenFields

  ! ----------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------ Interpolate all Fields Everywhere
  ! ----------------------------------------------------------------------------------------------

  subroutine interpolateEverything()
    call interpolateEverywhere(B1, 1, 1)
    call interpolateEverywhere(JCsup1, 1, 1)
    call interpolateEverywhere(B2, 0, 1)
    call interpolateEverywhere(JCsup2, 0, 1)
    call interpolateEverywhere(B3, 0, 0)
    call interpolateEverywhere(JCsup3, 0, 0)
    call interpolateEverywhere(E1, 0, 0)
    call interpolateEverywhere(JFsup1, 0, 0)
    call interpolateEverywhere(E2, 1, 0)
    call interpolateEverywhere(JFsup2, 1, 0)
    call interpolateEverywhere(E3, 1, 1)
    call interpolateEverywhere(JFsup3, 1, 1)
    call interpolateEverywhere(j3, 1, 1)
  end subroutine interpolateEverything

  ! ----------------------------------------------------------------------------------------------
  ! ----------------------------------------------------------------- Compute Fields at the Ground
  ! ----------------------------------------------------------------------------------------------

  subroutine computeGroundFields()
    integer                              :: h, i, k, m
    double complex, dimension(0:n1, 0:2) :: BrTemp
    ! Grab the radial magnetic field at RI. 
    BrTemp = Br()
    ! Compute the scalar magnetic potential at RE from BR at RI. 
    do i=0,n1
      do h=0,1
        PsiE(i, h) = sum( PsiE_Br(i, :)*BrTemp(:, h) )
      end do
    end do
    ! Take derivatives of PsiE to get B1 and B2 at RE. 
    do i=0,n1,2
      do h=0,1
        B1E(i, h) = d1D(PsiE, i, h)
        B2E(i, h) = d2(PsiE, i, h)
      end do
    end do
    ! Let's also compute the weights of the eigenvectors (hamonics). 
    do m=0,nModes-1
      do h=0,1
        alpha(m, h) = sum( alpha_Br(m, :)*BrTemp(:, h) )
      end do
    end do
  end subroutine computeGroundFields

  ! ----------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------- Write Fields to File
  ! ----------------------------------------------------------------------------------------------

  subroutine writeFields()
    ! Time stamp. 
    call writeReal('t.dat', t)
    call writeReal( 'driveScale.dat', getDriveScale(t) )
    ! Each field is defined on its own parity in i and k, but all fields are printed out at even i
    ! and even k. That means some interpolation is necessary before output. 
    call interpolateEvenFields()
    ! Fields at the ground are not used during computation, but we have access to them because the
    ! ground is the lower boundary to our Laplace's equation. We just have to go get them. 
    call computeGroundFields()
    ! Bulk magnetic fields mapped to physical coordinates in nT. 
    call writeComplexArray('Bx.dat', Bx(), onlydata=.true.)
    call writeComplexArray('By.dat', By(), onlydata=.true.)
    call writeComplexArray('Bz.dat', Bz(), onlydata=.true.)
    ! Bulk electric fields mapped to physical coordinates in mV/m. 
    call writeComplexArray('Ex.dat', Ex(), onlydata=.true.)
    call writeComplexArray('Ey.dat', Ey(), onlydata=.true.)
    call writeComplexArray('Ez.dat', Ez(), onlydata=.true.)
    ! Parallel current mapped to physical coordinates in uA/m^2. 
    call writeComplexArray('Jz.dat', jz(), onlydata=.true.)
    ! Scalar magnetic potential in... nT*Mm? This may not be a physically meaningful quantity. 
    call writeComplexArray('PsiE.dat', PsiE, onlydata=.true.)
    call writeComplexArray('PsiI.dat', PsiI, onlydata=.true.)
    ! Edge magnetic fields in nT. 
    call writeComplexArray('BfE.dat', Bf(RE), onlydata=.true.)
    call writeComplexArray('BfI.dat', Bf(RI), onlydata=.true.)
    call writeComplexArray('BqE.dat', Bq(RE), onlydata=.true.)
    call writeComplexArray('BqI.dat', Bq(RI), onlydata=.true.)
    call writeComplexArray('Br.dat', Br(), onlydata=.true.)
    ! Harmonic weights. One per harmonic, per hemisphere, per time step. 
    call writeComplexArray('eweights.dat', alpha, onlydata=.true., noskip=.true.)
  end subroutine writeFields

  ! ==============================================================================================
  ! ====================================================================== Physical Field "Arrays"
  ! ==============================================================================================

  ! ----------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------- Bulk Magnetic Fields
  ! ----------------------------------------------------------------------------------------------

  function Bx()
    double complex, dimension(0:n1, 0:n3) :: Bx
    Bx = (gsup11()*B1 + gsup13()*B3)*h1()
  end function Bx

  function By()
    double complex, dimension(0:n1, 0:n3) :: By
    By = gsup22()*B2*h2()
  end function By

  function Bz()
    double complex, dimension(0:n1, 0:n3) :: Bz
    Bz = B3/h3()
  end function Bz

  ! ----------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------- Bulk Electric Fields
  ! ----------------------------------------------------------------------------------------------

  function Ex()
    double complex, dimension(0:n1, 0:n3) :: Ex
    Ex = (gsup11()*E1 + gsup13()*E3)*h1()
  end function Ex

  function Ey()
    double complex, dimension(0:n1, 0:n3) :: Ey
    Ey = gsup22()*E2*h2()
  end function Ey

  function Ez()
    double complex, dimension(0:n1, 0:n3) :: Ez
    Ez = E3/h3()
  end function Ez

  ! ----------------------------------------------------------------------------------------------
  ! ---------------------------------------------------------- Bulk Magnetic Field Curl Components
  ! ----------------------------------------------------------------------------------------------

  function Fx()
    double complex, dimension(0:n1, 0:n3) :: Fx
    Fx = h1()*JFsup1/J()
  end function Fx

  function Fy()
    double complex, dimension(0:n1, 0:n3) :: Fy
    Fy = h2()*JFsup2/J()
  end function Fy

  function Fz()
    double complex, dimension(0:n1, 0:n3) :: Fz
    Fz = (g31()*JFsup1 + g33()*JFsup3)/( h3()*J() )
  end function Fz

  ! ----------------------------------------------------------------------------------------------
  ! ---------------------------------------------------------- Bulk Electric Field Curl Components
  ! ----------------------------------------------------------------------------------------------

  function Cx()
    double complex, dimension(0:n1, 0:n3) :: Cx
    Cx = h1()*JCsup1/J()
  end function Cx

  function Cy()
    double complex, dimension(0:n1, 0:n3) :: Cy
    Cy = h2()*JCsup2/J()
  end function Cy

  function Cz()
    double complex, dimension(0:n1, 0:n3) :: Cz
    Cz = (g31()*JCsup1 + g33()*JCsup3)/( h3()*J() )
  end function Cz

  ! ----------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------- Edge Magnetic Fields
  ! ----------------------------------------------------------------------------------------------

  ! Radial magnetic field at the ionosphere. 
  function Br()
    double complex, dimension(0:n1, 0:1)  :: Br
    double complex, dimension(0:n1, 0:n3) :: Bsup3
    Bsup3 = gsup31()*B1 + gsup33()*B3
    Br = Bsup3(:, kk)*hr()
  end function Br

  ! Azimuthal magnetic field at the ionosphere or at the ground. 
  function Bf(R0)
    double precision, optional, intent(in) :: R0
    double complex, dimension(0:n1, 0:1)   :: Bf
    if (present(R0) .and. R0 .eq. RE) then
      Bf = B2E/hf(RE)
    else
      Bf = B2I/hf()
    end if
  end function Bf

  ! Theta-aligned magnetic field at the ionosphere or at the ground. 
  function Bq(R0)
    double precision, optional, intent(in) :: R0
    double complex, dimension(0:n1, 0:1)   :: Bq
    if (present(R0) .and. R0 .eq. RE) then
      Bq = B1E/hq(RE)
    else
      Bq = B1I/hq()
    end if
  end function Bq

  ! ----------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------ Field-Aligned Current
  ! ----------------------------------------------------------------------------------------------

  function jz()
    double complex, dimension(0:n1, 0:n3) :: jz
    jz = j3/h3()
  end function jz

  ! ==============================================================================================
  ! ========================================================================= End of Fields Module
  ! ==============================================================================================

end module fields

! ################################################################################################
! ########################################################################################## Debug
! ################################################################################################

! The debug module allows us to check for NaNs and suspiciously large field values between calls
! to the main loop. 

module debug
  use loop
  use coefficients
  use fields
  use geometry
  implicit none

  contains

  ! ==============================================================================================
  ! ================================================================================ Debug Helpers
  ! ==============================================================================================

  ! ----------------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------- Check an Array for NaNs
  ! ----------------------------------------------------------------------------------------------

  logical function nanArray(arr, arrname, silent)
    double precision, dimension(0:, 0:), intent(in) :: arr
    character(len=*), intent(in)                    :: arrname
    logical, optional, intent(in)                   :: silent
    integer                                         :: i, k
    ! Figure out if there are any NaNs in the given array. 
    nanArray = any( isnan(arr) )
    ! If there are no NaNs, return that information. 
    if (nanArray .eq. .false.) return
    ! Also return if we're supposed to be running silently. 
    if ( present(silent) ) then
      if (silent .eq. .true.) return
    end if
    ! Otherwise, go through the array and report the location of all NaNs. 
    do i=0,size(arr, 1) - 1
      do k=0,size(arr, 2) - 1
        if ( isnan( arr(i, k) ) ) write(*,'(2a6, a10, a6, 2i6)') '', arrname, 'NaN', 'at', i, k
      end do
    end do
  end function nanArray

  ! ----------------------------------------------------------------------------------------------
  ! ------------------------------------------------ Check an Array for Conspicuously Large Values
  ! ----------------------------------------------------------------------------------------------

  logical function maxArray(arr, arrname, thresh, silent)
    double precision, dimension(0:, 0:), intent(in)     :: arr
    character(len=*), intent(in)                        :: arrname
    double precision, intent(in), optional              :: thresh
    logical, optional, intent(in)                       :: silent
    integer                                         :: i, k
    ! Check the maximum value against the threshhold (default 1e6). 
    if ( present(thresh) ) then
      maxArray = (maxval( abs(arr) ) .gt. thresh)
    else
      maxArray = (maxval( abs(arr) ) .gt. 1.e6)
    end if
    ! If running silently, return without saying anything. 
    if ( present(silent) ) then
      if (silent .eq. .true.) return
    end if
    ! Otherwise, report the maximum value in the array. 
    write(*,'(2a6, es10.2, a6, 2i6)') 'Max', arrname, maxval( abs(arr) ), 'at', i, k
  end function maxArray

  ! ==============================================================================================
  ! ================================================================================ Debug Drivers
  ! ==============================================================================================

  ! ----------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------ Check Fields for NaNs
  ! ----------------------------------------------------------------------------------------------

  logical function nanFields(silent)
    logical, optional, intent(in) :: silent
    logical                       :: s
    ! By default, we do ask the max checker to print out the values for each field. 
    if ( present(silent) ) then
      s = silent
    else
      s = .false.
    end if
    ! Check the covariant field components for NaNs, to avoid unnecessarily convoluting variables.
    nanFields =                nanArray(abs(B1), 'B1', silent=s)
    nanFields = nanFields .or. nanArray(abs(B2), 'B2', silent=s)
    nanFields = nanFields .or. nanArray(abs(B3), 'B3', silent=s)
    nanFields = nanFields .or. nanArray(abs(E1), 'E1', silent=s)
    nanFields = nanFields .or. nanArray(abs(E2), 'E2', silent=s)
    nanFields = nanFields .or. nanArray(abs(E3), 'E3', silent=s)
    nanFields = nanFields .or. nanArray(abs(B1I), 'B1I', silent=s)
    nanFields = nanFields .or. nanArray(abs(B2I), 'B2I', silent=s)
    ! If any NaNs are found, we'll return .true. to stop the run. 
  end function nanFields

  ! ----------------------------------------------------------------------------------------------
  ! -------------------------------------------------- Check Fields for Conspicuously Large Values
  ! ----------------------------------------------------------------------------------------------

  logical function maxFields(silent)
    logical, optional, intent(in) :: silent
    logical                       :: s
    ! By default, we do ask the max checker to print out the values for each field. 
    if ( present(silent) ) then
      s = silent
    else
      s = .false.
    end if
    ! Maximum values are checked on the physical fields, which can sometimes vary from the
    ! covariant fields by a few orders of magnitude. 
    maxFields =                maxArray(abs( Bx() ), 'Bx', silent=s)
    maxFields = maxFields .or. maxArray(abs( By() ), 'By', silent=s)
    maxFields = maxFields .or. maxArray(abs( Bz() ), 'Bz', silent=s)
    maxFields = maxFields .or. maxArray(abs( Ex() ), 'Ex', silent=s)
    maxFields = maxFields .or. maxArray(abs( Ey() ), 'Ey', silent=s)
    maxFields = maxFields .or. maxArray(abs( Ez() ), 'Ez', silent=s)
    maxFields = maxFields .or. maxArray(abs( Br() ), 'Br', silent=s)
    maxFields = maxFields .or. maxArray(abs( Bf() ), 'Bf', silent=s)
    maxFields = maxFields .or. maxArray(abs( Bq() ), 'Bq', silent=s)
    ! If any too-large value is found, return .true. to stop the run. 
  end function maxFields

  ! ==============================================================================================
  ! =========================================================================== Peek at the Fields
  ! ==============================================================================================

  ! Print out what the fields look like at X, Z.  
  subroutine peekFields(X, Z)
    double precision, intent(in)          :: X, Z
    integer, dimension(0:1)               :: ik
    integer                               :: i, k
    double complex, dimension(0:n1, 0:n3) :: arrx, arry, arrz
    double complex, dimension(0:n1, 0:1)  :: arrr, arrf, arrq
    ! Figure out the grid point closest to X, Z. 
    ik = minloc( sqrt( (r*sin(q) - X)**2 + (r*cos(q) - Z)**2 ) )
    ! Insist that both i and k be even. That's where all fields are defined (for output). 
    i = 2*(ik(0)/2)
    k = 2*(ik(1)/2)
    ! Include the peek location in the printout. 
    write(*,'(a, f5.2, a, f5.2, a)') 'Peeking at X, Z = ', X/RE, ', ', Z/RE, ' RE: '
    ! Use arrays to grab the physical fields, then print off the value at X, Z. 
    arrx = Bx()
    arry = By()
    arrz = Bz()
    write(*,'(a14, 3es9.1)') 'Bx By Bz', abs( arrx(i, k) ), abs( arry(i, k) ), abs( arrz(i, k) )
    arrx = Ex()
    arry = Ey()
    arrz = Ez()
    write(*,'(a14, 3es9.1)') 'Ex Ey Ez', abs( arrx(i, k) ), abs( arry(i, k) ), abs( arrz(i, k) )
    arrx = Cx()
    arry = Cy()
    arrz = Cz()
    write(*,'(a14, 3es9.1)') 'Cx Cy Cz', abs( arrx(i, k) ), abs( arry(i, k) ), abs( arrz(i, k) )
    arrx = Fx()
    arry = Fy()
    arrz = Fz()
    write(*,'(a14, 3es9.1)') 'Fx Fy Fz', abs( arrx(i, k) ), abs( arry(i, k) ), abs( arrz(i, k) )
    ! Also show the edge fields at the ionospheric end of that same field line. 
    arrr = Br()
    arrf = Bf()
    arrq = Bq()
    write(*,'(a14, 3es9.1)') 'Br Bf Bq', abs( arrr(i, 0) ), abs( arrf(i, 0) ), abs( arrq(i, 0) )
  end subroutine peekFields


  ! ==============================================================================================
  ! ===================================================================== Peek at the Coefficients
  ! ==============================================================================================

  ! Print out what the coefficients look like at X, Z.  
  subroutine peekCoefficients(X, Z)
    double precision, intent(in)          :: X, Z
    integer, dimension(0:1)               :: ik
    integer                               :: i, k
    double complex, dimension(0:n1, 0:n3) :: arr1, arr2, arr3
    ! Figure out the grid point closest to X, Z. 
    ik = minloc( sqrt( (r*sin(q) - X)**2 + (r*cos(q) - Z)**2 ) )
    i = ik(0)
    k = ik(1)
    ! Include the peek location in the printout. 
    write(*,'(a, f5.2, a, f5.2, a)') 'Peeking at X, Z = ', r(i, k)*sin( q(i, k) )/RE, ', ',      &
                                     r(i, k)*cos( q(i, k) )/RE, ': '
    write(*,'(a)') 'Directions are messy, but everything is scaled to a normal basis. '
    arr1 = E1_E1*sqrt( g11()/g11() )
    arr2 = E1_E2*sqrt( g22()/g11() )
    arr3 = E1_E3*sqrt( g33()/g11() )
    write(*,'(a10, 3es9.1)') 'E1_E*', abs( arr1(i, k) ), abs( arr2(i, k) ), abs( arr3(i, k) )
    arr1 = E1_JFsup1*sqrt( gsup11()/g11() )*J()
    arr2 = E1_JFsup2*sqrt( gsup22()/g11() )*J()
    arr3 = E1_JFsup3*sqrt( gsup33()/g11() )*J()
    write(*,'(a10, 3es9.1)') 'E1_F*', abs( arr1(i, k) ), abs( arr2(i, k) ), abs( arr3(i, k) )
    arr1 = 0
    arr2 = E1_j2drive*sqrt( g22()/g11() )
    arr3 = E1_j3*sqrt( g33()/g11() )
    write(*,'(a10, 3es9.1)') 'E1_j*', abs( arr1(i, k) ), abs( arr2(i, k) ), abs( arr3(i, k) )
    arr1 = E2_E1*sqrt( g11()/g22() )
    arr2 = E2_E2*sqrt( g22()/g22() )
    arr3 = E2_E3*sqrt( g33()/g22() )
    write(*,'(a10, 3es9.1)') 'E2_E*', abs( arr1(i, k) ), abs( arr2(i, k) ), abs( arr3(i, k) )
    arr1 = E2_JFsup1*sqrt( gsup11()/g22() )*J()
    arr2 = E2_JFsup2*sqrt( gsup22()/g22() )*J()
    arr3 = 0
    write(*,'(a10, 3es9.1)') 'E2_F*', abs( arr1(i, k) ), abs( arr2(i, k) ), abs( arr3(i, k) )
    arr1 = 0
    arr2 = E2_j2drive*sqrt( g22()/g22() )
    arr3 = 0
    write(*,'(a10, 3es9.1)') 'E2_j*', abs( arr1(i, k) ), abs( arr2(i, k) ), abs( arr3(i, k) )
    arr1 = 0
    arr2 = 0
    arr3 = E3_E3*sqrt( g33()/g33() )
    write(*,'(a10, 3es9.1)') 'E3_E*', abs( arr1(i, k) ), abs( arr2(i, k) ), abs( arr3(i, k) )
    arr1 = E3_JFsup1*sqrt( gsup11()/g33() )*J()
    arr2 = 0
    arr3 = E3_JFsup3*sqrt( gsup33()/g33() )*J()
    write(*,'(a10, 3es9.1)') 'E3_F*', abs( arr1(i, k) ), abs( arr2(i, k) ), abs( arr3(i, k) )
    arr1 = 0
    arr2 = 0
    arr3 = E3_j3
    write(*,'(a10, 3es9.1)') 'E3_j*', abs( arr1(i, k) ), abs( arr2(i, k) ), abs( arr3(i, k) )
  end subroutine peekCoefficients

  ! ==============================================================================================
  ! ========================================================================== End of Debug Module
  ! ==============================================================================================

end module debug

! ################################################################################################
! ################################################################################# Program Driver
! ################################################################################################

! All of the actual work is delegated out of the main program. All that happens here is a series
! of high-level subroutine calls. 

program tuna
  use coefficients
  use debug
  use fields
  use geometry
  use io
  use loop
  ! Time iterator. Fields are advanced in increments of dtOut. Each step, we call advanceFields()
  ! to evolve the fields by dtOut seconds, then we write the field values to file, then we check
  ! the field values for NaNs or conspicuously large values. 
  integer :: step
  ! Set the grid spacing. Based on the grid spacing, set linearized interpolation and
  ! differentiation weights. Compute numerical harmonics for Laplace's equation in the atmospheric
  ! shell. The geometry module also gives us access to metric tensor components and other factors
  ! for mapping between covariant, contravariant, and physical field components. 
  call geometrySetup()
  ! Read in ionospheric parameter profiles and map them to the grid. Most notably, the ionos
  ! module gives us access to ionospheric conductivities, the perpendicular electric constant, the
  ! time step, and the Alfven speed. 
  call ionosSetup()
  ! Initialize field arrays. The fields module does very little setup, actually, since fields all
  ! start at zero. The module also contains helper functions to cast the fields (which are stored
  ! in terms of their covariant coordinates) into physical coordinates, as well as to write those
  ! values to file. 
  call fieldSetup()
  ! Assemble coefficients. Tuna's speed is largely a product of its extensive precomputation of
  ! coefficients. Geometric scale factors and ionospheric parameter values are combined here, so
  ! that in the main loop we can advance fields in as few floating point operations as possible. 
  call coefficientSetup()





  ! Let's look at coefficients related to E1 at the equatorial corner. That's (i, k) = (0, 0). 

!  write(*,*) 'Looking at coefficients at r = ', r(0, 0), ', q = ', q(0, 0)
!  write(*,*) 'E1_E1 = ', E1_E1(0, 0)
!  write(*,*) 'E1_E2 = ', E1_E2(0, 0)
!  write(*,*) 'E1_E3 = ', E1_E3(0, 0)
!  write(*,*) 'E1_F1 = ', E1_JFsup1(0, 0)
!  write(*,*) 'E1_F2 = ', E1_JFsup2(0, 0)
!  write(*,*) 'E1_F3 = ', E1_JFsup3(0, 0)
!  write(*,*) 'E1_j3 = ', E1_j3(0, 0)

!  call peekCoefficients(0.*RE, 0.1*RE)
!  stop

!  E1_E1 = sqrt( n()*qe**2 / (me*epsPara) ) / 2*pi

!  write(*,'(a20, es9.1)') 'min w / wp = ', readParam('fdrive') / maxval( E1_E1 )
!  write(*,'(a20, es9.1)') 'max w / wp = ', readParam('fdrive') / minval( E1_E1 )

!  E2_E2 = n()*qe**2 / (me*sig0)

!  write(*,'(a20, es9.1)') 'min w / nu = ', readParam('fdrive') / maxval( E2_E2 )
!  write(*,'(a20, es9.1)') 'max w / nu = ', readParam('fdrive') / minval( E2_E2 )

!  write(*,'(a20, es9.1)') 'min w nu / wp wp = ', readParam('fdrive') / maxval( E1_E1**2 / E2_E2 )
!  write(*,'(a20, es9.1)') 'max w nu / wp wp = ', readParam('fdrive') / minval( E1_E1**2 / E2_E2 )

!  write(*,*) ''
!  write(*,*) 'Minval sig0 / eps0 ', minval( sig0/eps0 ), ' Hz'

!  write(*,*) ''
!  write(*,*) 'Maxval cc eps0^2 / sig0^2 ', maxval( cc*eps0**2/sig0**2 ), ' Mm^2'

!  write(*,*) ''
!  write(*,*) 'Minval m omega sigH / n e e ', minval( me*(0.02)*sigH/(n()*qe*qe) ), ' (dimensionless)'
!  write(*,*) 'Maxval m omega sigH / n e e ', maxval( me*(0.02)*sigH/(n()*qe*qe) ), ' (dimensionless)'

!  write(*,*) ''
!  write(*,*) 'Minval m omega sigP / n e e ', minval( me*(0.02)*sigP/(n()*qe*qe) ), ' (dimensionless)'
!  write(*,*) 'Maxval m omega sigP / n e e ', maxval( me*(0.02)*sigP/(n()*qe*qe) ), ' (dimensionless)'

!  write(*,*) ''
!  write(*,*) 'Minval epsPerp / eps0 ', minval( epsPerp / eps0 ), ' (dimensionless)'
!  write(*,*) 'Maxval epsPerp / eps0 ', maxval( epsPerp / eps0 ), ' (dimensionless)'

!  write(*,*) 'Or... 1/', 1e6/maxval( sqrt( n()*qe**2 / (2*pi*me*eps0) ) ), 'us'
!  write(*,'(a, es10.2, a)') 'n e e / m sig0 ... maxval ', maxval( n()*qe**2 / (me*sig0) ), 'Hz'
!  write(*,'(a, es10.2, a)') '       1/                 ', 1/maxval( n()*qe**2 / (me*sig0) ), 's'
!  write(*,'(a, es10.2, a)') '       *dt                ', dt*maxval( n()*qe**2 / (me*sig0) )
!  write(*,*) ''
!  write(*,'(a, es10.2, a)') 'n e e dt / m ... maxval ', dt*maxval( n()*qe**2 / me ), 'uA/m^2 per mV/m'
!  stop

  ! The subroutine advanceFields() moves each field forward in time by dtOut seconds. This usually
  ! means iterating through several thousand time steps. Then it returns to here, so that we can
  ! print out the field values and check them for errors. 
  do step=1,int(tMax/dtOut)
    ! Actual field advancement happens here. Almost all of the execution time is spent in calls to
    ! this subroutine. The main loop is parallelized using OpenMP. 
    call advanceFields()
    ! Print out a notification each time we return to main. 
!    write(*,'(a6, f6.2, a2)') 't = ', t, 's'

!    call peekFields(0.75*RE, 0.75*RE)

    ! Write field values (and the time stamp) to file. 
    call writeFields()
    ! Check fields for NaNs or conspicuously large values. If it seems that the run has become
    ! unstable, stop execution. 
    if ( maxFields(silent=.true.) ) stop 'Field value too large! '
    if ( nanFields(silent=.true.) ) stop 'NaN in field! '
  end do
  write(*,'(a)') 'Finished loop. '


end program tuna




!=============================================================================80
!>
!! Crossterm Sinusoidal MS
!!
!! This manufactured solution is three sinusoidals in each direction and
!! crossterm ones as well. The equation takes the form:
!!   f = a1
!!     + a2 *sin( a3 *pi*x  /l   + a4 *pi )
!!     + a5 *sin( a6 *pi*y  /l   + a7 *pi )
!!     + a8 *sin( a9 *pi*z  /l   + a10*pi )
!!     + a11*sin( a12*pi*x*y/l^2 + a13*pi )
!!     + a14*sin( a15*pi*y*z/l^2 + a16*pi )
!!     + a17*sin( a18*pi*x*z/l^2 + a19*pi )
!!
!! There are 5 such equations, one for each primitive variable. If a lower
!! dimension problem is being solved then the appropriate terms will be set to
!! zero.
!!
!! Inputs: for each line the 19 coefficients and then the reference length.
!!
!! Works for all forms of MMS
!<
!================================ constructor ================================80
!
! Creates an instance of the Crossterm Sinusoidal MS
!
!=============================================================================80
  function constructor( neq, coef_in, length_in )

    use set_constants,    only : zero, half, one, two, three, four, five,      &
                                 third, fourth, eighth, fifth, sixth, onep5
    use project_inputs,   only : dimen
    use reference_inputs, only : l_ref, l_ref_grid, rho_ref, a_ref
    use exact_inputs,     only : mms_number
    use message,          only : error_message

    integer,                                       intent(in) :: neq
    real(dp), dimension(cts_ncoef, neq), optional, intent(in) :: coef_in
    real(dp), dimension(neq),            optional, intent(in) :: length_in
    type(cts_ms)                                              :: constructor

    real(dp), dimension(cts_ninputs,neq) :: inputs

    continue

    allocate( constructor%length( neq ) )
    allocate( constructor%coef( cts_ncoef, neq ) )

    ! Set common variables
    call constructor%initialize_super( cts_manufactured, neq, cts_ninputs,     &
                                       neq, cts_name )

    ! Parameters are given
    if ( present(coef_in) .and. present(length_in) ) then

      constructor%coef   = coef_in
      constructor%length = length_in

    ! Default paramters
    else if (mms_number == -1) then

      !if (neq/=5) then
      !  err = error_message( routine_name,                                     &
      !                       "Default CTS MS only supports 5 equations!" )
      !end if

      constructor%length = l_ref

      ! Density
      constructor%coef(:,1) = [ one,                                           &
                                0.15_dp, two,   third,                         &
                               -0.1_dp,  one,  -fifth,                         &
                                0.2_dp,  half,  sixth,                         &
                               -0.05_dp, one,   fourth,                        &
                                0.25_dp, half, -third,                         &
                                0.15_dp, two,   eighth                         &
                              ]

      ! X-Velocity
      constructor%coef(:,2) = [ 40._dp,                                        &
                                2.5_dp, two*third,    half,                    &
                               -1.5_dp, onep5,       -eighth,                  &
                                2._dp, five*third,   third,                    &
                               -3._dp, onep5,        eighth,                   &
                                1.25_dp, five*eighth, -fourth,                 &
                                2.5_dp, half,         fourth                   &
                              ]

      ! Y-Velocity
      constructor%coef(:,3) = [ 40._dp,                                        &
                               -3.75_dp, five*third,  -half,                   &
                                2._dp, three,        zero,                     &
                                3._dp, four,         zero,                     &
                                3._dp, onep5,        eighth,                   &
                               -1.25_dp, five*eighth, -fourth,                 &
                               -3._dp, three*eighth, fourth                    &
                              ]

      ! Z-Velocity
      constructor%coef(:,4) = [ 45._dp,                                        &
                               -1.5_dp, half,          zero,                   &
                                4._dp, onep5,        -third,                   &
                                2.25_dp, three*fourth, -fourth,                &
                                1.5_dp, one,           fourth,                 &
                               -2.5_dp, half,         -third,                  &
                                2._dp, onep5,         third                    &
                              ]

      ! Pressure
      constructor%coef(:,5) = [ 100000._dp,                                    &
                                -20000._dp, half,          half,               &
                                 50000._dp, one,          -half,               &
                                 30000._dp, third,         zero,               &
                                -10000._dp, three*fourth, -fourth,             &
                                 20000._dp, half,         -fourth,             &
                                -30000._dp, onep5,         eighth              &
                              ]

      if (neq==6) then
        ! nutilde (these are non-dimensional)
        constructor%coef(:,6) = [ 100.0_dp,                                    &
                                  15.0_dp, four*third, eighth,                 &
                                 -10.0_dp, half,      -fifth,                  &
                                  20.0_dp, one,        sixth,                  &
                                 -5.0_dp,  onep5,     -half,                   &
                                  25.0_dp, two*third,  third,                  &
                                  15.0_dp, one,       -fourth                  &
                                ]
      end if

      if (neq==7) then
        ! TKE (these are non-dimensional)
        constructor%coef(:,6) = [ 1.0e-7_dp,                                   &
                                  1.5e-8_dp, four*third, eighth,               &
                                 -1.0e-8_dp, half,      -fifth,                &
                                  2.0e-8_dp, one,        sixth,                &
                                 -5.0e-9_dp, onep5,     -half,                 &
                                  2.5e-8_dp, two*third,  third,                &
                                  1.5e-8_dp, one,       -fourth                &
                                ]
        ! omega (these are non-dimensional)
        constructor%coef(:,7) = [ 1.0e-5_dp,                                   &
                                  2.5e-7_dp, four*third, eighth,               &
                                 -1.4e-7_dp, half,      -fifth,                &
                                  2.2e-7_dp, one,        sixth,                &
                                 -8.0e-8_dp, onep5,     -half,                 &
                                  3.5e-7_dp, two*third,  third,                &
                                  3.0e-7_dp, one,       -fourth                &
                                ]
      end if

    ! Parameters from file
    else

      inputs = constructor%read_mms_config()

      constructor%coef   = inputs( 1:cts_ncoef, : )
      constructor%length = inputs( cts_ninputs, : )


    end if

    ! Clear out unused dimensions
    if (dimen < 3) then
      constructor%coef([8,14,17],:) = zero
      constructor%coef(:,4)         = zero
    end if

    if (dimen < 2) then
      constructor%coef([5,11],:) = zero
      constructor%coef(:,3)      = zero
    end if

    ! Normalize coefficients to make MMS nondimensional
    constructor%length           = constructor%length                          &
                                 * l_ref_grid / l_ref
    constructor%coef(1,1)        = constructor%coef(1,1)                       &
                                 / rho_ref
    constructor%coef(1,2:4)      = constructor%coef(1,2:4)                     &
                                 / a_ref
    constructor%coef(1,5)        = constructor%coef(1,5)                       &
                                 / ( rho_ref * a_ref**2 )
    constructor%coef(2:17:3,1)   = constructor%coef(2:17:3,1)                  &
                                 / rho_ref
    constructor%coef(2:17:3,2:4) = constructor%coef(2:17:3,2:4)                &
                                 / a_ref
    constructor%coef(2:17:3,5)   = constructor%coef(2:17:3,5)                  &
                                 / ( rho_ref * a_ref**2 )

  end function constructor


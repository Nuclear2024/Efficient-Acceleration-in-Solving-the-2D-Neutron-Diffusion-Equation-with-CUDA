    subroutine solve_system(A, b, phi)
    
      real(8) :: phi(3999)
      real(8) :: M(3999,3999)

!     .. Parameters ..
      INTEGER          N, NRHS
      PARAMETER        ( N = 3999, NRHS = 1 )
      INTEGER          LDA, LDB
      PARAMETER        ( LDA = N, LDB = N )
!
!     .. Local Scalars ..
      INTEGER          INFO
!
!     .. Local Arrays ..
      INTEGER          IPIV( N )
      DOUBLE PRECISION A( LDA, N ), B(LDB )
!
!     .. External Subroutines ..
      EXTERNAL         DGESV
!
!     Solve the equations A*X = B.
      M = A
      CALL DGESV( N, NRHS, M, LDA, IPIV, B, LDB, INFO )
      phi = B
!     Check for the exact singularity.
!
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The diagonal element of the triangular factor of A,'
         WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
         WRITE(*,*)'A is singular; the solution could not be computed.'
         STOP
      END IF

      end subroutine solve_system
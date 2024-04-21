program neutron_diffusion

    implicit none

    ! Parâmetros do problema
    integer, parameter :: N = 4000
    real(8), parameter :: L = 4.0d2
    real(8), parameter :: D = 1.320
    real(8), parameter :: Sigma_a = 0.041
    real(8), parameter :: nu_Sigma_f = 0.045
    real(8), parameter :: k_eff = 1.0d0
    real(8) :: dx
    real(8) :: x(N+1)
    real(8), dimension(N-1, N-1) :: A
    real(8), dimension(N-1) :: b
    real(8) :: phi(N-1)
    real(8) :: k_eff_old, k_eff_new
    integer :: i
    
	integer :: Hora_Inicial            , Minuto_Inicial          , Segundo_Inicial          , &
			   Hora_Final              , Minuto_Final            , Segundo_Final            , &
			   Centesimos_Inicial      , Centesimos_Final
    

!    Abertura do arquivo para gravar o tempo de execução
!    ===================================================

     open(unit=1, file='TempoCalculo.dat', status='unknown')
	
     call gettim (Hora_Inicial, Minuto_Inicial, Segundo_Inicial, Centesimos_Inicial )
     
     
    ! Calcular dx e x
    dx = L / real(N, 8)
    
    ! Preencher o vetor x usando um loop DO
    do i = 0, N
        x(i+1) = -L/2.0d0 + real(i) * dx
    end do

    ! Inicializar matriz A e vetor b
    A = 0.0d0
    b = 0.0d0

    do i = 1, N-1
        ! Coeficiente na diagonal principal
        A(i, i) = 2.0d0 * D / dx**2 + Sigma_a

        ! Coeficientes nas diagonais secundárias
        if (i > 1) then
            A(i, i-1) = -D / dx**2
        end if
        if (i < N-1) then
            A(i, i+1) = -D / dx**2
        end if

        ! Termo fonte
        b(i) = nu_Sigma_f / k_eff
    end do

    ! Iterar até a convergência do k_eff
    k_eff_old = k_eff

    phi = 0.0D0

    do
        ! Resolver o sistema linear
        call solve_system(A, b, phi)

        ! Calcular o novo k_eff
        k_eff_new = nu_Sigma_f * sum(phi(1:N-1)) * dx / (k_eff * sum(Sigma_a * phi(1:N-1)) * dx)

        ! Verificar a convergência
        if (abs(k_eff_new - k_eff_old) < 1.0d-6) then
            exit
        else
            k_eff_old = k_eff_new
            b = (nu_Sigma_f / k_eff_new) * phi(1:N-1)

        end if
    end do

!    Finaliza e grava o tempo de cálculo
!    ===================================
     
     call gettim (Hora_Final, Minuto_Final, Segundo_Final, Centesimos_Final )

     write (1,*) Hora_Inicial, Minuto_Inicial, Segundo_Inicial, Centesimos_Inicial
	 write (1,*) '                                                               '
     write (1,*) Hora_Final, Minuto_Final, Segundo_Final, Centesimos_Final

     close(unit=1)
     
     open(unit=2, file='Saida.dat', status='unknown')
          
        write (2,*) k_eff_new
        write (2,*) ''
        write (2,*) phi
        
    close(unit=2)
        
        
        
end program neutron_diffusion

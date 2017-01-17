
  !============== Variable Module ====================================
  MODULE my_vrbls
    IMPLICIT NONE

    double precision, parameter :: eps    = 1.d-14            ! very small number
    double precision, parameter :: tol    = 0.15d0            ! tolerance for Cor

    integer :: Lx
    double precision :: Beta
    integer :: NBlck=1024
    integer :: MXNCPU = 100
    integer :: iblck,TotSamp,Seed
    double precision :: JJ,Q

    character(8 ), parameter :: ident = 'hs_sqa 0'            ! identifier
    character(12), parameter :: file1 = 'ext_hs_sqa0.dat'        ! datafile
    character(100) :: file3 
    character(100) :: file4
    !-----------------------------------------------------------------

    !-- Observables --------------------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    integer, parameter :: NObs_b = 23                         ! #basic observables
    integer, parameter :: NObs_c = 3                          ! #composite observables
    integer, parameter :: NObs   = NObs_b+NObs_c              ! Total # observables
    !-----------------------------------------------------------------

    !-- Statistics ---------------------------------------------------
    ! THIS IS PROJECT-INDEPENDENT 
    double precision, allocatable :: Quan(:,:), Err(:,:)       ! 1st--#quan.  2nd--#block
    double precision, allocatable :: Ave(:), Dev(:), Cor(:)           ! average, error bars, and correlation of observables
    !-----------------------------------------------------------------
  END MODULE my_vrbls

  !=====Main routine for bond percolation on square lattice ==========
  PROGRAM main
    use my_vrbls
    implicit none
    integer :: whichone,stat,NB,ind
    character(100) :: filename
    integer :: i,j,flag

    allocate(Quan(1:NObs, 1:MXNCPU), Err(1:NObs, 1:MXNCPU))
    allocate(Ave(1:NObs), Dev(1:NObs), Cor(1:NObs))

    open(8,file="hs_sqa0.dat") 
    ind=0
    flag=1
    do while(ind<=MXNCPU)
        read(8, *,iostat=stat);
        if(stat/=0) exit

        ind=ind+1
	read(8,40) filename, Lx, beta, JJ, Q, Seed, TotSamp, NB, NBlck
	40 format(a8,i6,3f10.6,i12,i8,i8,i8)

	do j = 1, Nobs
	  read(8,41) i, Quan(j,ind), Err(j,ind), Cor(j)
	  41 format(i4,2f25.15,f12.5)
	enddo
    enddo
    close(8)

    if(ind>=1) then
        call stat_alan(ind)
        call write2file()
    endif

  CONTAINS


  SUBROUTINE stat_alan(nfile)
    implicit none
    integer          :: i, j, k, k0, nfile
    logical          :: prt
    double precision :: devn, devp, nor
    double precision, allocatable :: Aux(:)

    nor  = 1.d0/(nfile*1.d0)

    ! -- calculate average -------------------------------------------
    do j = 1, NObs
      Ave(j) = nor*Sum(Quan(j,1:nfile))
    enddo

    DO j = 1, NObs
      Dev(j) = 0.d0
      do k = 1,  nfile
        devn   = Quan(j,k)-Ave(j)
        Dev(j) = Dev(j)+devn*devn
      enddo 
      Dev(j)   = Dev(j)*nor
      Dev(j)   = dsqrt(Dev(j)/(nfile-1.d0))
    ENDDO 

    return
  END SUBROUTINE stat_alan

  SUBROUTINE write2file 
    implicit none
    integer       :: i, j, k, Nwri

    open (8,file=file1, access='append') 
    write(8, *);          write(6,*)

    write(8,40) ident, Lx, beta, JJ, Q, Seed, TotSamp, NBlck, iblck
    40 format(a8,i6,3f10.6,i12,i8,i8,i8)

    do j = 1, Nobs
      write(8,41) j, Ave(j), Dev(j), Cor(j)
      41 format(i4,2f25.15,f12.5)
    enddo
    close(8)
    return
  END SUBROUTINE write2file

END PROGRAM main
!=====================================================================

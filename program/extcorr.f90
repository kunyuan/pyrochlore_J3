
  !============== Variable Module ====================================
  MODULE my_vrbls
    IMPLICIT NONE

    double precision, parameter :: eps    = 1.d-14            ! very small number
    double precision, parameter :: tol    = 0.15d0            ! tolerance for Cor
    double precision, parameter :: Pi = 3.141592653

    integer :: Lx, Dim=3, Vol, num
    double precision :: Beta
    double precision :: JJ,Q

    double precision,allocatable :: Momentum(:,:)
    double precision,allocatable :: StaticCorr(:,:,:), FreqCorr(:,:,:)
    double precision,allocatable :: WeightStaticCorr(:), WeightFreqCorr(:,:)
    double precision,allocatable :: AveStaticCorr(:,:), AveFreqCorr(:,:)
    double precision,allocatable :: DevStaticCorr(:,:), DevFreqCorr(:,:)

    character(8 ), parameter :: ident = 'hs_sqa 0'            ! identifier
    character(12), parameter :: file1 = 'hs_sqa0.dat'        ! datafile
    character(100) :: file3 
    character(100),allocatable :: corrfile(:), freqfile(:, :)
    !-----------------------------------------------------------------

    !-- Observables --------------------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    integer, parameter :: MxOmega = 128
    integer, parameter :: Nk = 2
    integer, parameter :: NSub = 4
    !-----------------------------------------------------------------
  END MODULE my_vrbls

  !=====Main routine for bond percolation on square lattice ==========
  PROGRAM main
    use my_vrbls
    implicit none
    integer :: whichone,stat,NB,ind
    double precision :: k1, k2, k3, scorr, fcorr, err, omega
    character(100) :: filename, head
    integer :: i,j,k, flag, numstat, numfreq(1:10),  weight
    double precision :: freqweight(1:10), statweight
    double precision :: nor
    character(8) :: string

    !print *, 'Outputfile,MinBlck'
    read  *,  num,Lx,Beta
    Vol = Lx*Lx*Lx*NSub

    allocate(corrfile(1:num))
    allocate(freqfile(1:Nk, 1:num))

    do whichone = 1, num
      if(whichone<=10) then
          write(corrfile(whichone),"(a14,i1,a4)") "static_corr_1_",whichone-1,".txt"
          write(freqfile(1, whichone),"(a10,i1,a4)") "corr_k1_1_",whichone-1,".txt"
          write(freqfile(2, whichone),"(a10,i1,a4)") "corr_k2_1_",whichone-1,".txt"
      else
        write(corrfile(whichone),"(a14,i2,a4)")   "static_corr_1_",whichone-1,".txt"
        write(freqfile(1, whichone),"(a10,i2,a4)")   "corr_k1_1_",whichone-1,".txt"
        write(freqfile(2, whichone),"(a10,i2,a4)")   "corr_k2_1_",whichone-1,".txt"
      endif
    enddo

    allocate(StaticCorr(1:NSub, 1:Vol, 1:num))
    allocate(FreqCorr(1:Nk, 1:MxOmega+1, 1:num))
    allocate(WeightStaticCorr(1:num))
    allocate(WeightFreqCorr(1:Nk, 1:num))
    allocate(AveStaticCorr(1:NSub, 1:Vol))
    allocate(AveFreqCorr(1:Nk, 1:MxOmega+1))
    allocate(DevStaticCorr(1:NSub, 1:Vol))
    allocate(DevFreqCorr(1:Nk, 1:MxOmega+1))
    StaticCorr(:,:,:) = 0.d0
    FreqCorr(:,:,:) = 0.d0

    allocate(Momentum(1:Nk, 1:Dim))
    Momentum(1, :) = (/0.d0, 0.d0, 0.d0/)
    Momentum(2, :) = (/0.d0, 0.d0, 2.d0*Pi/)

    numstat = 0
    numfreq(:) = 0
    statweight=0.d0
    freqweight(:) = 0.d0

    LoopFile: do i = 1, num
	do j = 1, Nk
	  open (8,file=freqfile(j, i)) 
	  read(8, *, iostat=stat) head, weight
	  if(stat/=0) cycle LoopFile
	  WeightFreqCorr(j, numfreq(j)+1) = weight*1.d0
	  read(8, 41) head, k1, k2, k3
	  41 format(a8,3f10.6)

	  do k = 1, MxOmega+1
	    read(8,*) omega, fcorr, err
	    FreqCorr(j, k, numfreq(j)+1) = fcorr
	  enddo
	  read(8, *)
	  numfreq(j) = numfreq(j)+1
	  close(8)
	enddo

      open(9, file=corrfile(i))
      read(9, *, iostat=stat) head, weight
      if(stat/=0) exit
      WeightStaticCorr(numstat+1) = weight*1.d0
      read(9,42) head
      42 format(a20)

      do k = 1, NSub
        do j = 1, Vol
          read(9,43) scorr, string
          43 format(f25.15, a3)
	  StaticCorr(k, j, numstat+1) = scorr
        enddo
        read(9, *)
      enddo
      numstat = numstat + 1
      close(9)
    enddo LoopFile

    do j = 1, Nk
	if(numfreq(j)/=numfreq(1)) then
	    print *, numfreq(j), numfreq(1), "num of files not equal!"
	    exit
	endif
    enddo

    print *, "num of valid static correlation files:", numstat
    print *, "num of valid k-frequency correlation files:", numfreq(1)

    call stat_alan2d(StaticCorr(:,:,:), WeightStaticCorr(:), NSub, Vol, numstat, AveStaticCorr(:,:), DevStaticCorr(:,:))
    do i = 1, Nk
	call stat_alan(FreqCorr(i,:,:), WeightFreqCorr(i,:),  MxOmega+1, numfreq(i), AveFreqCorr(i,:), DevFreqCorr(i,:))
    enddo

    call write2file()

  CONTAINS

  SUBROUTINE stat_alan(Data1d, Weight, d1, num, Ave, Dev)
    implicit none
    integer          :: i, j, k, k0, d1,  num
    double precision :: devn, nor
    double precision :: Data1d(1:d1, 1:num), Ave(1:d1), Dev(1:d1)
    double precision :: Weight(1:num)

    nor  = 1.d0/sum(Weight(1:num))

    ! -- calculate average -------------------------------------------
    do i = 1, d1
	Ave(i) = nor*Sum(Data1d(i, 1:num)*Weight(1:num))
    enddo

    ! -- calculate error and t=1 correlation for basics obs.-------- 
    Dev(:) = 0.d0
    DO i = 1, d1
	do k = 1,  num
	  devn   = Data1d(i, k)-Ave(i)
	  Dev(i) = Dev(i)+devn*devn*Weight(k)
	enddo 
	Dev(i)   = Dev(i)*nor
	Dev(i)   = dsqrt(Dev(i)/(num-1.d0))
    ENDDO 
    return
  END SUBROUTINE stat_alan

  SUBROUTINE stat_alan2d(Data2d, Weight, d1, d2, num, Ave, Dev)
    implicit none
    integer          :: i, j, k, k0, d1, d2, num
    double precision :: devn, nor
    double precision :: Data2d(1:d1, 1:d2, 1:num), Ave(1:d1, 1:d2), Dev(1:d1, 1:d2)
    double precision :: Weight(1:num)

    nor  = 1.d0/sum(Weight(1:num))

    ! -- calculate average -------------------------------------------
    do i = 1, d1
	do j = 1, d2
	    Ave(i, j) = nor*Sum(Data2d(i, j, 1:num)*Weight(1:num))
	enddo
    enddo

    ! -- calculate error and t=1 correlation for basics obs.-------- 
    Dev(:, :) = 0.d0
    DO i = 1, d1
	do j = 1, d2
	    do k = 1,  num
	      devn   = Data2d(i,j, k)-Ave(i, j)
	      Dev(i, j) = Dev(i, j)+devn*devn*Weight(k)
	    enddo 
	    Dev(i, j)   = Dev(i, j)*nor
	    Dev(i, j)   = dsqrt(Dev(i, j)/(num-1.d0))
	enddo
    ENDDO 
    return
  END SUBROUTINE stat_alan2d


  SUBROUTINE write2file 
    implicit none
    integer       :: i, j, k, Nwri
    open(9,file="static_corr.txt") 
    write(9, *) "{ 'Correlations': [ ["
    do j = 1, Vol
      write(9,47) AveStaticCorr(1, j), ', '
      47 format(f25.15, a2)
    enddo
    write(9, *) "], ["
    do j = 1, Vol
      write(9,47) AveStaticCorr(2, j), ', '
    enddo
    write(9, *) "], ["
    do j = 1, Vol
      write(9,47) AveStaticCorr(3, j), ', '
    enddo
    write(9, *) "], ["
    do j = 1, Vol
      write(9,47) AveStaticCorr(4, j), ', '
    enddo
    write(9, *) "]]} "
    close(9)

    open (10,file="corr_k1.txt") 
    write(10, *) "##k=", Momentum(1,:)
    do j = 1, MxOmega+1
      write(10,*) 2.d0*Pi*(j-1.d0)/Beta, AveFreqCorr(1, j), DevFreqCorr(1, j)
    enddo
    write(10, *) 
    close(10)

    open (11,file="corr_k2.txt") 
    write(11, *) "##k=", Momentum(2, :)
    do j = 1, MxOmega+1
      write(11,*) 2.d0*Pi*(j-1.d0)/Beta, AveFreqCorr(2, j), DevFreqCorr(2, j)
    enddo
    write(11, *) 
    close(11)
    return
  END SUBROUTINE write2file

END PROGRAM main
!=====================================================================





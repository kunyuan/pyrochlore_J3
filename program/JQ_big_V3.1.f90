  !******************************************************************
  ! Bond percolation on the square lattice with built-in error bars
  ! using the blocking technique. 

  ! Error bars are calculated using the blocking technique. 
  ! Let ' \chi^2 = \sum_{t=1}^{T} <(O_t -<O>)^2> ' be the squared chi for
  ! 'T' blocks of observable 'O'. Assuming each block of data is independent
  ! of each other, the error bar of 'O' is given by 'Dev = \sqrt{\chi^2/T(T-1)}.

  ! Reliabity of the obtained errors is monitored by t=1 correlation,
  ! for which tolerance is set by variable 'tol' (default: tol=0.15d0).

  ! To calculate the error bar of composite quantity like Binder ratios 'Q' and 
  ! specific heat 'C', an auxilliary variable is defined as derivative of
  ! 'Q' or 'C'. For instance,
  ! (1). 'Q = <O_2>/<O_1>, the auxillary variable is 
  !      AQ = (1/<O_1>^2) (<O_1> Q_2 - <Q_2> Q_1), which is again a time series
  ! (2). 'C = V(<O^2>-<O>^2),
  !      AC = V(O^2 - 2<O> O).
  ! The error bar of 'Q' ('C') is taken that of 'AQ' ('AC').

  ! Results are written into a special file 'dat.***' if the number of
  ! blocks is less than 125 or correlation is too big. Data in each 
  ! block will be also printed out in this case.

  ! Default number of extensive simulation is 'NBlck=1024'.

  ! For test purpose, for which huge amount of information will be 
  ! printed out, 'NBlck' should be set smaller but >2.

  ! Dynamical behavior is not studied.

  ! 'my_vrbls.f90', 'carlo.f90', 'monte.f90', and 'measure.f90'
  ! need to be modified for new projects.

  !  Look for 'PROJECT-DEPENDENT'.

  !  Author: Youjin Deng
  !  Date  : Jan. 25th, 2011.
  !*******************************************************************


  !============== Variable Module ====================================
  MODULE my_vrbls
    IMPLICIT NONE

    !-- common parameters and variables ------------------------------
    ! THIS IS ALMOST PROJECT-INDEPENDENT 
    double precision, parameter :: tm32   = 1.d0/(2.d0**32.d0)
    double precision, parameter :: eps    = 1.d-14            ! very small number
    double precision, parameter :: tol    = 0.15d0            ! tolerance for Cor
    logical                     :: prt                        ! flag for write2file
    integer,          parameter :: Mxint  = 2147483647        ! maximum integer
    integer,          parameter :: Mnint  =-2147483647        ! minimum integer

    integer,          parameter :: MxBlck =8192              ! maximum number of blocks for statistics
    integer,          parameter :: MnBlck = 1                 ! minimum number of blocks

    integer            :: NBlck                               ! # blocks
    integer            :: Nswee                               ! # MC sweeps between measurements
    integer            :: Ntoss                               ! # MC sweeps for equilibration
    integer            :: Nsamp                               ! # samples in unit 'NBlck'
    integer            :: TotSamp                             
    integer            :: IsLoad

    integer            :: Vol                                 ! actual volumn
    double precision   :: Norm                                ! normalization factor (1/Vol)

    double precision   :: rmuv,radv                           ! auxillary variables to draw a random vertex
    double precision   :: rmue,rade                           ! auxillary variables to draw a random neighboring edge
    !-----------------------------------------------------------------

    !-- parameters and variables -------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    integer, parameter :: D   = 2                             ! dimensionality
    integer, parameter :: nnb = 4                             ! # neighboring edges
    integer, parameter :: MxL = 192                        ! maximum lattice size 
    integer,parameter  :: MaxKink=512
    integer,parameter  :: NSave=1
    integer, parameter :: MxV = MxL*MxL                       ! maximum volumn
    integer            :: Lx, Ly                              ! actual linear lattice size 

    character(8 ), parameter :: ident = 'hs_sqa 0'            ! identifier
    character(11), parameter :: file1 = 'hs_sqa0.dat'       ! datafile
    character(100) :: file2 
    character(100) :: file3 
    character(100) :: file4 
    !-----------------------------------------------------------------

    !-- Lattice and state --------------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    integer, dimension(nnb,MxV) :: Ngs                        ! neighbouring-vertex table
    integer(1), dimension(nnb)  :: dx,dy                      ! auxillary variables to measure wrapping probability
    !-----------------------------------------------------------------

    !-- Observables --------------------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    integer, parameter :: NObs_b = 25                        ! #basic observables
    integer, parameter :: NObs_c = 12                         ! #composite observables
    integer, parameter :: NObs   = NObs_b+NObs_c              ! Total # observables
    !-----------------------------------------------------------------

    !-- Statistics ---------------------------------------------------
    ! THIS IS PROJECT-INDEPENDENT 
    double precision, allocatable :: Quan(:)                  ! Measured quantities
    double precision, allocatable :: Obs(:,:)                 ! 1st--#quan.  2nd--#block
    double precision, allocatable :: Ave(:), Dev(:), Cor(:)   ! average, error bars, and correlation of observables
    !-----------------------------------------------------------------

    !-- Random-number generator---------------------------------------
    ! THIS IS PROJECT-INDEPENDENT 
    integer, parameter           :: mult=32781
    integer, parameter           :: mod2=2796203, mul2=125
    integer, parameter           :: len1=9689,    ifd1=471
    integer, parameter           :: len2=127,     ifd2=30
    integer, dimension(1:len1)   :: inxt1
    integer, dimension(1:len2)   :: inxt2
    integer, dimension(1:len1)   :: ir1
    integer, dimension(1:len2)   :: ir2
    integer                      :: ipnt1, ipnf1
    integer                      :: ipnt2, ipnf2
    integer, parameter           :: mxrn = 10000
    integer, dimension(1:mxrn)   :: irn(mxrn)

    integer                      :: Seed                      ! random-number seed
    integer                      :: nrannr                    ! random-number counter
    !-----------------------------------------------------------------

    !-- time-checking variables --------------------------------------
    ! THIS IS PROJECT-INDEPENDENT 
    character( 8)         :: date
    character(10)         :: time
    character(5 )         :: zone
    integer, dimension(8) :: tval
    double precision      :: t_prev, t_curr, t_elap
    integer               :: h_prev, h_curr 
    double precision      :: t_simu, t_mcmc, t_meas
    !-----------------------------------------------------------------
    
    !-- data structs -------------------------------------------------

    double precision       :: Comp                             ! common variable in probability
    double precision       :: JJ                                ! common variable 
    double precision       :: Q
    double precision       :: JJQ
    double precision       :: beta                                ! common variable 
    double precision       :: Number
    double precision       :: Const
    integer                :: MxKink
    integer                :: rank
                                               
    double precision      :: KinkTime(MxV,MaxKink)      !Time point of every kink
    integer               :: SegmentNum(MxV)             !Track the number of segments
    integer               :: SegmentState(MxV,MaxKink)  !Track the state of segments    
    integer               :: TailPoint(MxV)
    
    
    integer               :: PrevSite(MxV,MaxKink)
    integer               :: NextSite(MxV,MaxKink)
    integer               :: NeighSite(MxV,MaxKink)        
    integer               :: NeighVertexNum(MxV,MaxKink)
    integer               :: NearestSite(MxV,MaxKink,nnb) 
    double precision      :: KinkWeight(MxV,MaxKink)
    integer               :: KinkWeightChangeVertex(MaxKink)
    integer               :: KinkWeightChangeSite(MaxKink)
    double precision      :: KinkWeightChange(MaxKink)
    integer               :: ListIndex=1
    
    integer               :: PairedSite(MxV,MaxKink)
    integer               :: PairedVertexNum(MxV,MaxKink)
    integer               :: PairedKinkNum
    
    double precision      :: IraTime
    integer               :: IraVertex
    double precision      :: MashaTime
    integer               :: MashaVertex
    integer :: IraSite,MashaSite
    
    !--- auxillary data for measuring ----------------------------
    double precision     :: W0
    double precision     :: W1
    double precision     :: TotalW
    double precision     :: EnergyCheck
    integer    :: dWx,dWy,dWt
    integer    :: Wx,Wy,Wt
    integer    :: WindX,WindY,WindT
    integer    :: Direction
    integer    :: State
    character(13)       :: RunName(13)
    double precision    :: RunNum(10)
    double precision    :: RunNum0(10)
    double precision  :: P1=1.0
    double precision  :: P2=1.0
    double precision  :: P3=1.0
    double precision  :: P4=2.0
    double precision  :: P5=2.0
    double precision  :: P6=1.0
  END MODULE my_vrbls

  !=====Main routine for bond percolation on square lattice ==========
  PROGRAM main
    use my_vrbls
    implicit none
!!!!! MPI WRAPPING !!!!!
  include 'mpif.h'
  integer :: comm_size, ierr, outproc
  double precision :: step_beta,P
  integer :: itoss,isamp,iblck,stat

  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, comm_size, ierr )

  if(rank>9) then
    write(file2,2345) rank
  else
    write(file2,2343) rank
  endif

2343 format('collect',i1,'.dat')
2345 format('collect',i2,'.dat')

  if(rank>9) then
    write(file3,2349) rank
  else
    write(file3,2347) rank
  endif

2347 format('c',i1,'.dat')
2349 format('c',i2,'.dat')

  if(rank>9) then
    write(file4,2353) rank
  else
    write(file4,2351) rank
  endif

2351 format('status',i1,'.dat')
2353 format('status',i2,'.dat')

!--these will erase previous outputs---
    !open(1,file=hs_sqa_file)
    !write(1,*)'   '
    !close(1)
!-----------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!  MPI !!!!!

  open(1,file="in00")
  read(1,*) Lx, beta, JJ,Q, Ntoss,Nsamp,Nswee,Seed,NBlck,step_beta,IsLoad
  close(1)
  beta=beta+step_beta*rank
  Seed=Seed+rank
  MxKink=MaxKink-2
!!!!!!!!!!!!!!!!!!!!! 
    JJ=JJ/4.0
    Q=Q/16.0
    JJQ=JJ+2*Q

    P=P1+P2+P3+P4+P5+P6
    P1=P1/P
    P2=P1+P2/P
    P3=P2+P3/P
    P4=P3+P4/P
    P5=P4+P5/P

    TotSamp = Nsamp*NBlck/1024
    RunNum(:)=0
    W0=0.0;W1=0.0
    RunName(1)="Delete worm:"
    RunName(2)="Move   worm:"
    RunName(3)="Create kink:"
    RunName(4)="Delete kink:"
    RunName(5)="Create Dkink:"
    RunName(6)="Delete Dkink:"
    RunName(7)="Create worm:"
    ListIndex=1

    open(18,file=file4,access="append")
       write(18,*) "Initialization..."
    close(18)

    !--- Initialization ----------------------------------------------
    call t_elapse(0)         ! '0' for the 1st time.
    call carlo
    call t_elapse(1)         ! '1' for initialization

    open(18,file=file4,access="append")
       write(18,*) "Initialization done!"
    close(18)

    call winding_number()
    Wx=WindX
    Wy=WindY
    Wt=WindT

    EnergyCheck=potential_energy()

    call check_config()

    open(18,file=file4,access="append")
       write(18,*) EnergyCheck
       write(18,*) Wx,Wy,Wt
       write(18,*) "Begin Simulation"
    close(18)

    !--- Equilibration -----------------------------------------------
    do iblck = 1, NBlck
        open(18,file=file4,access="append")
           write(18,*) iblck,"toss monte"
        close(18)
        call check_config()
        do itoss = 1, Ntoss
          call monte
        enddo
    enddo
    call t_elapse(2)         ! '2' for equilibration

    call check_config()
    call saveconfig()

  !open (8,file=hs_sqa_file , access='append') 
  !!--- Simulation --------------------------------------------------
    !write(8,40) rank, Lx, beta, JJ*4.0, Q*16.0, Seed, Nsamp*NBlck/1024
!40  format(i6,i6,3f10.6,2i8)    
    !--- Simulation --------------------------------------------------
    do iblck = 1, NBlck
      write(6,*) rank,iblck
      open(18,file=file4,access="append")
         write(18,*) iblck,"do monte",Number
      close(18)
      do isamp = 1, Nsamp
        call monte
        call measure
        call coll_data(iblck)
      enddo
      call norm_Nsamp(iblck)
      call t_elapse(-1)      ! '-1' just for trace the time
      if(iblck-int(iblck/NSave)*NSave==0) then
          call check_config()
          call midwrite2file(iblck)
          call saveconfig()
          open(18,file=file4,access="append")
             write(18,*) iblck,"save configuration"
          close(18)
      endif
    enddo
    call t_elapse(3)         ! '3' for markov-chain, including simulation and measurement

    !--- Statistics --------------------------------------------------
    call stat_alan
    call t_elapse(3)         ! '3' for markov-chain, including simulation and measurement
    call write2file
    call saveconfig
    call t_elapse(4)         ! '4' total time 

!!!!!! MPI  ending !!!!!!
  call MPI_finalize(ierr); STOP
!!!!!!!!!!!!!!!!

  CONTAINS

  !*******************************************************************
  !            Beginning of PROJECT-DEPENDENT part 
  !*******************************************************************
  !==============Initialization ======================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE carlo
    implicit none
    integer  :: i,j, k, x, y, xp, xm, ym, yp
    integer  :: kb,temp
    
    Ly=Lx
    Vol=Lx*Ly
    Comp=1.0/beta
    Number=0
    !Comp=1/(beta/<n>)^2/<n>=<n>/beta^2 
    !<n> is the average kink number per site
    ! <n>=beta*2Vol*<SxSx+SySy>/Vol~beta
    
    if(Lx<=2) then
        write(6,*) 'L must bigger than 2'
    endif

    if(Lx-Lx/2*2==1) then
        write(6,*) 'Give me even L,please!'
    endif

    if(Lx>MxL) then
      write(6,*) 'L < MxL?';                              stop
    endif

    !if(Ntoss==0) then
      !write(6,*) 'to compare simulation and measurement &
      !&  time, please set Ntoss/=0';                      stop
    !endif

    if((NBlck>MxBlck).or.(NBlck<MnBlck)) then
      write(6,*) 'MnBlck <= NBlck <=MxBlck?';             stop
    endif

    !if((NBlck>200).and.(NBlck/=MxBlck)) then
      !write(6,*) '"NBlck>200" is supposed for extensive &
      !& simulation. "NBlck=MxBlk" is suggested!';         stop
    !endif

    write(6,40) Lx
    40 format(' Simulation is on L=',i6,2x,'square lattice.')

    write(6,45) Beta
    45 format(' Simulation is on Beta=',f12.8,2x,'Time Line.')

    write(6,41) JJ
    41 format(' Coupling J is',f12.8)

    write(6,42) Ntoss*NBlck
    42 format(' Will throw away    ',i10,2x,'steps ')

    write(6,43) Nsamp*NBlck
    43 format(' Will simulate      ',i10,2x,'steps ')

    write(6,44) Nswee      
    44 format(' A step contains    ',i10,2x,'sweeps')

    !-- define lattice -----------------------------------------------
    lat_def: do k = 1,Vol
      y  = (k-1)/Lx+1;           x = k-(y-1)*Lx

      xm =  -1 ;      if(x== 1) xm = Lx-1
      xp =  +1 ;      if(x==Lx) xp = 1-Lx

      ym =  -Lx;      if(y== 1) ym = Vol-Lx
      yp =  +Lx;      if(y==Ly) yp = Lx-Vol

      Ngs(1,k) = k+xp;    Ngs(4,k) = k+xm
      Ngs(2,k) = k+yp;    Ngs(3,k) = k+ym

      !            2
      !            |
      !       4---   ---1
      !            |
      !            3

    enddo lat_def

    !-- define time-line list ----------------------------------------
    SegmentState(:,:)=0
    do k=1,Vol
        SegmentNum(k)=1;
        SegmentState(k,1)=1
        SegmentState(k,2)=-1
        KinkTime(k,1)=0.0
        KinkTime(k,2)=beta
        NeighSite(k,1)=0
        NeighSite(k,2)=0
        PrevSite(k,1)=0
        NextSite(k,1)=2
        NearestSite(k,1,:)=1
        PrevSite(k,2)=1
        NextSite(k,2)=0
        NearestSite(k,2,:)=2
        do i=3,MaxKink-1
            NextSite(k,i)=i+1
            PrevSite(k,i)=i-1
            NearestSite(k,i,:)=0
        enddo
        PrevSite(k,MaxKink)=MaxKink-1
        NextSite(k,MaxKink)=0
        NearestSite(k,MaxKink,:)=0
        !TailPoint points to the first empty node
        TailPoint(k)=3
        PairedSite(k,:)=-1 !Set it to -1,then we know it is not pairedkink
        PairedVertexNum(k,:)=0
    enddo
    if(IsLoad==1) then
        call readconfig()
    endif

    !-- set random-number generator-------------------------------------
    call set_RNG
    
    !-- measurement initialization -----------------------------------
    Norm = 1.d0/Vol/beta
    allocate(Quan(1:NObs_b))
    allocate(Obs(1:NObs, 1:NBlck));       Obs = 0.d0
    allocate(Ave(1:NObs), Dev(1:NObs), Cor(1:NObs))
    return
  END SUBROUTINE carlo
  !===================================================================


  !==============Simulation ==========================================
  !! THIS IS PROJECT-DEPENDENT 
    SUBROUTINE monte
        implicit none
        integer :: i,j
        logical :: flag
        integer :: Site
        double precision :: x

        do i=1,Nswee
            dWx=0;dWy=0;dWt=0
            Number=Number+1
            call create_worm(flag)
            if(flag) then
                do while(flag)
                    Number=Number+1
                    call choose_Ira()
                    x=rn()
                    if(x<P1) then
                        call move_worm()
                    else if(x<P2) then
                       call create_kink()
                    else if(x<P3) then
                       call delete_kink()
                    else if(x<P4) then
                       call create_paired_kink()
                    else if(x<P5) then
                       call delete_paired_kink()
                    else
                       call annihilate_worm(flag)
                    endif
                enddo
            endif
            Wx=Wx+(dWx/Lx)
            Wy=Wy+(dWy/Lx)
            Wt=Wt+dWt
        enddo
        return
    END SUBROUTINE monte

    logical FUNCTION is_in_same_segment(IraVertex,IraSite,MashaVertex,MashaSite)
        implicit none
        integer ::     IraVertex
        integer ::     MashaVertex
        integer ::  IraSite,MashaSite
        if(IraVertex/=MashaVertex) then
            is_in_same_segment=.false.
        else if(IraSite/=NextSite(MashaVertex,MashaSite).and.IraSite/=PrevSite(MashaVertex,MashaSite)) then
            is_in_same_segment=.false.
        else
            is_in_same_segment=.true.
        endif
        return
    END FUNCTION

    SUBROUTINE choose_Ira()
        implicit none
        integer :: vertex,site
        double precision :: time
        if(rn()<0.5) then
            Vertex=IraVertex
            IraVertex=MashaVertex
            MashaVertex=vertex

            site=IraSite
            IraSite=MashaSite
            MashaSite=site

            time=IraTime
            IraTime=MashaTime
            MashaTime=time

            direction=-direction
        endif
        return
    END SUBROUTINE
    
    SUBROUTINE create_worm(flag)
        !Choose an Interval to create a worm with Ira and Masha
        implicit none
        integer  :: Site1,Site2,Vertex,temp,i
        double precision :: prevtime,nexttime,itime,mtime,time,Tmin,Tmax
        double precision :: Factor,total_E,TempTime
        logical :: flag
        Vertex=int(rn()*Vol)+1
        !temp=int(SegmentNum(Vertex)*rn())
        !Site1=1
        !do i=1,temp
            !Site1=NextSite(Vertex,Site1)
        !enddo
        Site1=floor(rn()*MxKink)+2
        do while(SegmentState(Vertex,Site1)==0)
            Site1=floor(rn()*MxKink)+2
        enddo
        if(Site1==2)Site1=1
        Site2=NextSite(Vertex,Site1)
        Tmin=KinkTime(Vertex,Site1)
        Tmax=KinkTime(Vertex,Site2)
        itime=rn()*(Tmax-Tmin)+Tmin
        if(itime==Tmin) return
        mtime=rn()*(Tmax-Tmin)+Tmin
        if(mtime==Tmin)return
        if(itime>mtime) then
            TempTime=itime;itime=mtime;mtime=TempTime
        endif
        RunNum0(7)=RunNum0(7)+1
        total_E=single_energy_kinkweight(Vertex,Site1,itime,mtime,Factor,0)
        Factor=Comp*Factor*SegmentNum(vertex)*(Tmax-Tmin)**2
        if((Factor>=1.0 .and. total_E<=0.0) .or. rn()<Factor*exp(-total_E)) then              
            RunNum(7)=RunNum(7)+1
            EnergyCheck=EnergyCheck+total_E

            if(SegmentState(Vertex,Site1)==-1) then
                direction=-1
            else
                direction=1
            endif

            IraVertex=Vertex
            MashaVertex=Vertex
            IraTime=itime
            MashaTime=mtime
            IraSite=insert_site(Vertex,Site1,IraTime)
            MashaSite=insert_site(Vertex,IraSite,MashaTime)

            call change_kink_weight()

            flag=.true.
        else
            flag=.false.
        endif
        return
    END SUBROUTINE create_worm
    
    SUBROUTINE annihilate_worm(flag)
        implicit none
        double precision :: Factor,total_E,Tmin,Tmax
        integer :: EndSite
        logical :: flag
        if(.not. is_in_same_segment(IraVertex,IraSite,MashaVertex,MashaSite)) then
            return
        endif
        if(IraTime<MashaTime) then
            EndSite=NextSite(MashaVertex,MashaSite)
            Tmin=KinkTime(IraVertex,PrevSite(IraVertex,IraSite))
            Tmax=KinkTime(IraVertex,EndSite)
            total_E=single_energy_kinkweight(IraVertex,IraSite,IraTime,MashaTime,Factor,0)
        else
            EndSite=NextSite(IraVertex,IraSite)
            Tmin=KinkTime(IraVertex,PrevSite(MashaVertex,MashaSite))
            Tmax=KinkTime(IraVertex,EndSite)
            total_E=single_energy_kinkweight(MashaVertex,MashaSite,MashaTime,IraTime,Factor,0)
        endif
        Factor=Factor/Comp/(SegmentNum(IraVertex)-2)/(Tmax-Tmin)**2

        RunNum0(1)=RunNum0(1)+1
        if((Factor>=1.0 .and. total_E<=0.0) .or. rn()<Factor*exp(-total_E)) then              

            RunNum(1)=RunNum(1)+1
            EnergyCheck=EnergyCheck+total_E

            call remove_site(IraVertex,IraSite)
            call remove_site(MashaVertex,MashaSite)

            call change_kink_weight()
            flag=.false.
        else
            flag=.true.
        endif
        return    
    END SUBROUTINE annihilate_worm
    
    SUBROUTINE move_worm()
        implicit none
        integer :: Site1,Site2,indicator,TempSite
        double precision :: Time1,Time2,tarTime,delta_Time,OldIraTime
        double precision :: total_E,factor,ran
        indicator=1
        Site1=PrevSite(IraVertex,IraSite)
        Site2=NextSite(IraVertex,IraSite)
        if(Site1==1) then
            TempSite=PrevSite(IraVertex,2)
            Time1=KinkTime(IraVertex,TempSite)
            Time2=KinkTime(IraVertex,Site2)
            delta_Time=Time2+beta-Time1
            tarTime=rn()*delta_Time
            if(tarTime==0.d0)return
            tarTime=tarTime+Time1-beta
            if(tarTime<0) then 
                indicator=2
                tarTime=tarTime+beta
                total_E=single_energy_kinkweight(IraVertex,1,0.0d0,IraTime,Factor,0)+ &
                    & single_energy_kinkweight(IraVertex,TempSite,tarTime,beta,Factor,1)
            endif
        else if(Site2==2) then
            TempSite=NextSite(IraVertex,1)
            Time1=KinkTime(IraVertex,Site1)
            Time2=KinkTime(IraVertex,TempSite)
            delta_Time=Time2+beta-Time1
            tarTime=rn()*delta_Time
            if(tarTime==0.d0)return
            tarTime=tarTime+Time1
            if(tarTime>beta) then
                indicator=3
                tarTime=tarTime-beta
                total_E=single_energy_kinkweight(IraVertex,IraSite,IraTime,beta,Factor,0) &
                    & +single_energy_kinkweight(IraVertex,1,0.d0,tarTime,Factor,1)
            end if
        else
            Time1=KinkTime(IraVertex,Site1)
            Time2=KinkTime(IraVertex,Site2)
            delta_Time=Time2-Time1
            tarTime=rn()*delta_Time
            if(tarTime==0.d0)return
            tarTime=tarTime+Time1
        endif
        if(indicator==1) then
            if(tarTime<IraTime) then
                total_E=single_energy_kinkweight(IraVertex,Site1,TarTime,IraTime,Factor,0)
                indicator=0
             else
                total_E=single_energy_kinkweight(IraVertex,IraSite,IraTime,TarTime,Factor,0)
             endif
        endif
        RunNum0(2)=RunNum0(2)+1
        if((Factor>=1.0 .and. total_E<=0.0) .or. rn()<Factor*exp(-total_E)) then              

            RunNum(2)=RunNum(2)+1
            EnergyCheck=EnergyCheck+total_E
            
            OldIraTime=IraTime
            IraTime=tarTime
            if(indicator==0) then
                !Ira moves in the same interval ,the simplest case
                call move_site(IraVertex,IraSite,IraTime,1)
            else if(indicator==1) then
                call move_site(IraVertex,IraSite,IraTime,0)
            else if(indicator==2)  then
                !Ira moves from the first interval to the last
                dWt=dWt-direction

                SegmentState(IraVertex,1)=SegmentState(IraVertex,IraSite)
                SegmentState(IraVertex,2)=-SegmentState(IraVertex,1)
                call remove_site(IraVertex,IraSite)
                IraSite=insert_site(IraVertex,TempSite,IraTime)
            else
                !Ira moves from the last interval to the first
                dWt=dWt+direction

                SegmentState(IraVertex,1)=-SegmentState(IraVertex,IraSite)
                SegmentState(IraVertex,2)=-SegmentState(IraVertex,1)
                call remove_site(IraVertex,IraSite)
                IraSite=insert_site(IraVertex,1,IraTime)
            endif
            call change_kink_weight()
        endif
    END SUBROUTINE move_worm
    
    SUBROUTINE create_kink()
        implicit none
        integer :: TarVertex,OldIraVertex,OldIraSite
        integer :: Site0,Site1,Site2
        double precision    :: TargetTime,NewKinkTime
        integer :: IsPast,num
        double precision    :: Factor,total_E,KinkFactor
        double precision    :: weight

        num=int(rn()*4)+1
        TarVertex=Ngs(num,IraVertex)
        Site0=NearestSite(IraVertex,IraSite,num)
        if(SegmentState(IraVertex,IraSite)==SegmentState(TarVertex,Site0)) then
            IsPast=1
            Site2=Site0
            Site1=PrevSite(IraVertex,IraSite)
        else
            IsPast=0
            Site2=NextSite(TarVertex,Site0)
            Site1=NextSite(IraVertex,IraSite)
        endif
        !Site0,Site2 on TarVertex,Site1 on IraVertex
        TargetTime=compare(KinkTime(IraVertex,Site1),KinkTime(TarVertex,Site2),IsPast)
        NewKinkTime=rn()*(IraTime-TargetTime)+TargetTime
        if(NewKinkTime==TargetTime)return

        if(IsPast==1) then
            total_E=double_energy_kinkweight(IraVertex,Site1,TarVertex,Site2,num,NewKinkTime,IraTime,KinkFactor)
            weight=JJQ-two_plaquette_weight(IraVertex,Site1,TarVertex,Site2,num,NewKinkTime)
        else
            total_E=double_energy_kinkweight(IraVertex,IraSite,TarVertex,Site0,num,IraTime,NewKinkTime,KinkFactor)
            weight=JJQ-two_plaquette_weight(IraVertex,IraSite,TarVertex,Site0,num,NewKinkTime)
        endif
        Factor=4.0*abs(IraTime-TargetTime)*weight*KinkFactor

        RunNum0(3)=RunNum0(3)+1
        if((Factor>=1.0 .and. total_E<=0.0) .or. rn()<Factor*exp(-total_E)) then              

            RunNum(3)=RunNum(3)+1
            dWx=dWx+direction*go_where(num,1)
            dWy=dWy+direction*go_where(num,-1)
            EnergyCheck=EnergyCheck+total_E

            OldIraVertex=IraVertex
            OldIraSite=IraSite

            IraVertex=TarVertex
            !Must move old Ira To NewKinkTime first,then NearestSite of New Ira
            !will point to Ira
            call move_site(OldIraVertex,OldIraSite,NewKinkTime,IsPast)
            !Move IraSite to new kink
            if(IsPast==1) then
                Site2=insert_site(TarVertex,Site0,NewKinkTime)
                !Notice we use Site0=Site1 here
                IraSite=insert_site(TarVertex,Site2,IraTime)
            else
                IraSite=insert_site(TarVertex,Site0,IraTime)
                !Notice we use Site0=PrevSite(TarVertex,Site1) here
                Site2=insert_site(TarVertex,IraSite,NewKinkTime)
            endif
            NeighVertexNum(OldIraVertex,OldIraSite)=num
            NeighVertexNum(IraVertex,Site2)=5-num
            !Set the neighbor vertex of the new kink
            NeighSite(OldIraVertex,OldIraSite)=Site2
            NeighSite(IraVertex,Site2)=OldIraSite
            !Set the neighbor site of the new kink
            KinkWeight(OldIraVertex,OldIraSite)=weight
            KinkWeight(IraVertex,Site2)=weight

            !IraTime will remain the same here.
            call change_kink_weight()
        endif
        return
    END SUBROUTINE create_kink

    SUBROUTINE delete_kink()
        implicit none
        integer :: IsPast,i,num,OldIraVertex,OldIraSite
        integer :: TargetSite,TarVertex,NSite,Site0
        double precision    :: TargetTime,OldKinkTime
        double precision    :: Factor,total_E,weight,ran,KinkFactor
        IsPast=int(rn()*2)
        TargetSite=nearest_site(IraVertex,IraSite,IsPast)
        if(TargetSite==1 .or. TargetSite==2) return
        if(IraVertex==MashaVertex .and. TargetSite==MashaSite) return
        if(PairedSite(IraVertex,TargetSite)/=-1) return
        num=NeighVertexNum(IraVertex,TargetSite)
        NSite=NeighSite(IraVertex,TargetSite)
        Site0=NearestSite(IraVertex,IraSite,num)
        TarVertex=Ngs(num,IraVertex)
        OldKinkTime=KinkTime(IraVertex,TargetSite)
        if(IsPast==1) then
            if(NSite/=Site0)return
            TargetTime=compare(KinkTime(IraVertex,PrevSite(IraVertex,TargetSite)),   &
                &  KinkTime(TarVertex,PrevSite(TarVertex,NSite)),1)
            !total_E=double_energy(IraVertex,TargetSite,TarVertex,NSite,num,OldKinkTime,IraTime)
            total_E=double_energy_kinkweight(IraVertex,TargetSite,TarVertex,NSite,num,OldKinkTime,IraTime,KinkFactor)
            Site0=NextSite(TarVertex,NSite)
            !KinkFactor=kink_factor_double(IraVertex,IraSite,TarVertex,Site0,OldKinkTime,IraTime,num)
        else
            if(NSite/=NextSite(TarVertex,Site0))return
            TargetTime=compare(KinkTime(IraVertex,NextSite(IraVertex,TargetSite)),   &
                &  KinkTime(TarVertex,NextSite(TarVertex,NSite)),0)
            !total_E=double_energy(IraVertex,IraSite,TarVertex,Site0,num,IraTime,OldKinkTime)
            total_E=double_energy_kinkweight(IraVertex,IraSite,TarVertex,Site0,num,IraTime,OldKinkTime,KinkFactor)
            !KinkFactor=kink_factor_double(IraVertex,TargetSite,TarVertex,NSite,IraTime,OldKinkTime,num)
        endif

        weight=KinkWeight(IraVertex,TargetSite)
        Factor=0.25/abs(IraTime-TargetTime)/weight*KinkFactor

        RunNum0(4)=RunNum0(4)+1
        ran=rn()
        if((total_E<=0.0 .and. ran<Factor) .or. ran<Factor*exp(-total_E)) then              

            RunNum(4)=RunNum(4)+1
            dWx=dWx+direction*go_where(num,1)
            dWy=dWy+direction*go_where(num,-1)
            EnergyCheck=EnergyCheck+total_E

            OldIraVertex=IraVertex
            OldIraSite=IraSite
            IraVertex=TarVertex
            IraSite=NeighSite(OldIraVertex,TargetSite)
            call remove_site(OldIraVertex,OldIraSite) 
            call remove_site(OldIraVertex,TargetSite)
            call move_site(IraVertex,IraSite,IraTime,1-IsPast)

            call change_kink_weight()
        endif
        return
    END SUBROUTINE delete_kink

    SUBROUTINE create_paired_kink()
        implicit none
        integer :: OldIraVertex,OldIraSite
        integer :: Vertex1,Vertex2,Vertex3,Vertex4
        integer :: Site0,Site1,Site2,Site3,Site4
        double precision    :: TargetTime,NewPairedTime
        integer :: IsPast,num,NeighNum
        double precision    :: Factor,total_E,weight,ran,KinkFactor
        !1: Vertex with Ira
        !2: Vertex sharing kink with Ira
        !3: NeighborVertex of Vertex1
        !4: neighborVertex of Vertex2
        !     1------kink----2
        !   (Ira)            |
        !     |              |
        !     |              |
        !     |              |
        !     3-----kink-----4

        num=int(rn()*4)+1
        IsPast=int(rn()*2)
        Vertex3=Ngs(num,IraVertex)
        Site3=NearestSite(IraVertex,IraSite,num)
        if(IsPast==0) Site3=NextSite(Vertex3,Site3)
        if(PairedSite(Vertex3,Site3)/=-1) return
        if(Site3==1 .or. Site3==2) return
        if(Vertex3==MashaVertex .and. Site3==MashaSite) return
        NeighNum=NeighVertexNum(Vertex3,Site3)
        if(NeighNum==num .or. NeighNum==5-num) return
        Vertex1=IraVertex
        Vertex2=Ngs(NeighNum,Vertex1)
        Site0=NearestSite(Vertex1,IraSite,NeighNum)
        Vertex4=Ngs(NeighNum,Vertex3)
        Site4=NeighSite(Vertex3,Site3)

        if(IsPast==1 .and. SegmentState(Vertex1,IraSite)/=SegmentState(Vertex2,Site0)) return
        if(IsPast==0 .and. SegmentState(Vertex1,IraSite)==SegmentState(Vertex2,Site0)) return

        NewPairedTime=KinkTime(Vertex3,Site3)

        if(IsPast==1) then
            Site2=Site0
            Site1=PrevSite(Vertex1,IraSite)
            if(KinkTime(Vertex1,Site1)>NewPairedTime .or. KinkTime(Vertex2,Site2)>NewPairedTime)return
            total_E=double_energy_kinkweight(IraVertex,Site1,Vertex2,Site2,NeighNum,NewPairedTime,IraTime,KinkFactor)
        else
            Site2=NextSite(Vertex2,Site0)
            Site1=NextSite(Vertex1,IraSite)
            if(KinkTime(Vertex1,Site1)<NewPairedTime .or. KinkTime(Vertex2,Site2)<NewPairedTime)return
            total_E=double_energy_kinkweight(IraVertex,IraSite,Vertex2,Site0,NeighNum,IraTime,NewPairedTime,KinkFactor)
        endif
        Factor=8.0*Q/KinkWeight(Vertex3,Site3)*KinkFactor
        !(1/4)*(1/2)*K1*Exp(-Ei)*P(i-->f)=(1/2)*K2*Exp(-Ef)*P(f-->i)
        RunNum0(5)=RunNum0(5)+1
        ran=rn()
        if((total_E<=0.0 .and. ran<Factor) .or. ran<Factor*exp(-total_E)) then              

            RunNum(5)=RunNum(5)+1
            dWx=dWx+direction*go_where(NeighNum,1)
            dWy=dWy+direction*go_where(NeighNum,-1)
            EnergyCheck=EnergyCheck+total_E

            Site1=IraSite
            IraVertex=Vertex2

            call move_site(Vertex1,Site1,NewPairedTime,IsPast)
            !Move IraSite to new kink
            if(IsPast==1) then
                Site2=insert_site(Vertex2,Site0,NewPairedTime)
                !Notice we use Site0=Site1 here
                IraSite=insert_site(Vertex2,Site2,IraTime)
            else
                IraSite=insert_site(Vertex2,Site0,IraTime)
                !Notice we use Site0=PrevSite(TarVertex,Site1) here
                Site2=insert_site(Vertex2,IraSite,NewPairedTime)
            endif
            NeighVertexNum(Vertex1,Site1)=NeighNum
            NeighVertexNum(Vertex2,Site2)=5-NeighNum
            !Set the neighbor vertex of the new kink
            NeighSite(Vertex1,Site1)=Site2
            NeighSite(Vertex2,Site2)=Site1
            !Set the neighbor site of the new kink
            PairedKinkNum=PairedKinkNum+1
            PairedSite(Vertex1,Site1)=Site3
            PairedSite(Vertex2,Site2)=Site4
            PairedSite(Vertex3,Site3)=Site1
            PairedSite(Vertex4,Site4)=Site2
            PairedVertexNum(Vertex1,Site1)=num
            PairedVertexNum(Vertex2,Site2)=num
            PairedVertexNum(Vertex3,Site3)=5-num
            PairedVertexNum(Vertex4,Site4)=5-num
            call change_kink_weight()
        endif
    END SUBROUTINE create_paired_kink

    SUBROUTINE delete_paired_kink()
        implicit none
        integer :: IsPast,num,OldIraVertex,OldIraSite,NeighNum,IsChange
        integer :: Vertex1,Vertex2,Vertex3,Vertex4
        integer :: Site1,Site2,Site3,Site4,NSite
        double precision    :: PairedTime,TargetTime
        double precision    :: Factor,total_E,weight,KinkFactor
        !1: Vertex sharing kink with Ira
        !2: Vertex with Ira
        !3: NeighborVertex of Vertex1
        !4: neighborVertex of Vertex2
        !     1------kink----2
        !     |            (Ira)
        !     |              |
        !     |              |
        !     |              |
        !     3-----kink-----4
        IsPast=int(rn()*2)
        Site2=nearest_site(IraVertex,IraSite,IsPast)
        if(PairedSite(IraVertex,Site2)==-1) return
        NeighNum=NeighVertexNum(IraVertex,Site2)
        Vertex2=IraVertex
        Vertex1=Ngs(NeighNum,Vertex2)
        Site1=NeighSite(Vertex2,Site2)
        num=PairedVertexNum(Vertex2,Site2)
        Vertex3=Ngs(num,Vertex1)
        Vertex4=Ngs(num,Vertex2)
        Site3=PairedSite(Vertex1,Site1)
        Site4=PairedSite(Vertex2,Site2)
        PairedTime=KinkTime(Vertex2,Site2)

        if(rn()<0.5 .and. SegmentState(Vertex1,Site1)==SegmentState(Vertex4,Site4)) then
            NSite=NearestSite(Vertex2,IraSite,num)
            if(IsPast==1) then
                if(Site4/=NSite)return
                if(KinkTime(Vertex3,NextSite(Vertex3,Site3))<IraTime)return
                total_E=double_energy_kinkweight(Vertex2,Site2,Vertex4,Site4,num,PairedTime,IraTime,KinkFactor)
            else
                if(Site4/=NextSite(Vertex4,NSite))return
                if(KinkTime(Vertex3,PrevSite(Vertex3,Site3))>IraTime)return
                total_E=double_energy_kinkweight(Vertex2,IraSite,Vertex4,NSite,num,IraTime,PairedTime,KinkFactor)
            endif
            IsChange=1
            weight=JJQ-two_plaquette_weight(Vertex3,Site3,Vertex1,Site1,5-num,PairedTime)
        else 
            NSite=NearestSite(Vertex2,IraSite,NeighNum)
            if(IsPast==1) then
                if(Site1/=NSite)return
                if(KinkTime(Vertex3,NextSite(Vertex3,Site3))<IraTime)return
                total_E=double_energy_kinkweight(Vertex2,Site2,Vertex1,Site1,NeighNum,PairedTime,IraTime,KinkFactor)
            else
                if(Site1/=NextSite(Vertex1,NSite))return
                if(KinkTime(Vertex3,PrevSite(Vertex3,Site3))>IraTime)return
                total_E=double_energy_kinkweight(Vertex2,IraSite,Vertex1,NSite,NeighNum,IraTime,PairedTime,KinkFactor)
            endif
            IsChange=0
            weight=JJQ-two_plaquette_weight(Vertex3,Site3,Vertex4,Site4,5-NeighNum,PairedTime)
        endif

        !!(1/4)*(1/2)*K1*Exp(-Ei)*P(i-->f)=(1/2)*K2*Exp(-Ef)*P(f-->i)
        Factor=weight/8.0/Q*KinkFactor

        RunNum0(6)=RunNum0(6)+1
        if((Factor>=1.0 .and. total_E<=0.0) .or. rn()<Factor*exp(-total_E)) then              

            RunNum(6)=RunNum(6)+1
            EnergyCheck=EnergyCheck+total_E

            OldIraVertex=IraVertex
            OldIraSite=IraSite

            PairedKinkNum=PairedKinkNum-1
            PairedSite(Vertex1,Site1)=-1
            PairedSite(Vertex2,Site2)=-1
            PairedSite(Vertex3,Site3)=-1
            PairedSite(Vertex4,Site4)=-1

            call remove_site(OldIraVertex,OldIraSite) 
            call remove_site(OldIraVertex,Site2)
            if(IsChange) then
                dWx=dWx+direction*go_where(num,1)
                dWy=dWy+direction*go_where(num,-1)

                IraVertex=Vertex4
                IraSite=Site4
                call move_site(Vertex4,Site4,IraTime,1-IsPast)
                NeighVertexNum(Vertex1,Site1)=num
                NeighSite(Vertex1,Site1)=Site3
                NeighVertexNum(Vertex3,Site3)=5-num
                NeighSite(Vertex3,Site3)=Site1
                kinkweight(Vertex1,Site1)=weight
                kinkweight(Vertex3,Site3)=weight
            else
                dWx=dWx+direction*go_where(NeighNum,1)
                dWy=dWy+direction*go_where(NeighNum,-1)

                IraVertex=Vertex1
                IraSite=Site1
                call move_site(Vertex1,Site1,IraTime,1-IsPast)
                kinkweight(Vertex3,Site3)=weight
                kinkweight(Vertex4,Site4)=weight
            endif
            call change_kink_weight()
        endif
        return
    END SUBROUTINE delete_paired_kink

    integer FUNCTION go_where(num,IsX)
        implicit none
        integer :: num,IsX
        if(IsX==1) then
            if(num==1) then
                go_where=1
            else if(num==4) then
                go_where=-1
            endif
        else
            if(num==2) then
                go_where=1
            else if(num==3) then
                go_where=-1
            endif
        endif
    end FUNCTION

    integer FUNCTION insert_site(CurrentVertex,CurrentSite,NewSiteTime)
        implicit none
        integer,intent(in):: CurrentSite,CurrentVertex
        integer :: NewSite,CNSite,i,Temp,NVertex
        double precision :: NewSiteTime
        if(SegmentNum(CurrentVertex)>=MxKink) then
            open(18,file=file4,access="append")
            write(18,*) "ERROR!  Too much kinks!"
            close(18)
            call MPI_finalize(ierr); STOP
        endif

        NewSite=TailPoint(CurrentVertex)
        TailPoint(CurrentVertex)=NextSite(CurrentVertex,NewSite)

        PrevSite(CurrentVertex,NewSite)=CurrentSite
        CNSite=NextSite(CurrentVertex,CurrentSite)
        NextSite(CurrentVertex,NewSite)=CNSite
        KinkTime(CurrentVertex,NewSite)=NewSiteTime

        PrevSite(CurrentVertex,CNSite)=NewSite
        NextSite(CurrentVertex,CurrentSite)=NewSite

        SegmentNum(CurrentVertex)=SegmentNum(CurrentVertex)+1
        do i=1,4
            Temp=NearestSite(CurrentVertex,CNSite,i)
            NVertex=Ngs(i,CurrentVertex)
            do while(KinkTime(NVertex,Temp)>NewSiteTime)
                Temp=PrevSite(NVertex,Temp)
            enddo
            NearestSite(CurrentVertex,NewSite,i)=Temp
        enddo

        call change_nearest(CurrentVertex,CNSite,CurrentSite,NewSite)
        SegmentState(CurrentVertex,NewSite)=-SegmentState(CurrentVertex,CurrentSite);
        insert_site=NewSite
        return
    END FUNCTION
    
    SUBROUTINE remove_site(CurrentVertex,CurrentSite)
        implicit none
        integer :: CurrentVertex,CurrentSite,BeginSite,EndSite

        BeginSite=PrevSite(CurrentVertex,CurrentSite)
        EndSite=NextSite(CurrentVertex,CurrentSite)

        NextSite(CurrentVertex,BeginSite)=EndSite 
        PrevSite(CurrentVertex,EndSite)=BeginSite

        SegmentState(CurrentVertex,CurrentSite)=0

        call change_nearest(CurrentVertex,EndSite,CurrentSite,BeginSite)

        NextSite(CurrentVertex,CurrentSite)=TailPoint(CurrentVertex)
        TailPoint(CurrentVertex)=CurrentSite
        SegmentNum(CurrentVertex)=SegmentNum(CurrentVertex)-1
        return
    END SUBROUTINE remove_site

    SUBROUTINE move_site(CurrentVertex,CurrentSite,NewSiteTime,IsPast)
        implicit none
        integer :: CurrentVertex,CurrentSite,IsPast
        double precision :: NewSiteTime
        integer :: i,Temp,NVertex,PSite
        if(IsPast==1) then
            do i=1,4
                NVertex=Ngs(i,CurrentVertex)
                Temp=NearestSite(CurrentVertex,CurrentSite,i)
                do while(KinkTime(NVertex,Temp)>NewSiteTime)
                    NearestSite(NVertex,Temp,5-i)=CurrentSite
                    Temp=PrevSite(NVertex,Temp)
                enddo
                if(KinkTime(NVertex,Temp)==NewSiteTime)NearestSite(NVertex,Temp,5-i)=CurrentSite
                NearestSite(CurrentVertex,CurrentSite,i)=Temp
            enddo
        else
            PSite=PrevSite(CurrentVertex,CurrentSite)
            do i=1,4
                NVertex=Ngs(i,CurrentVertex)
                Temp=NearestSite(CurrentVertex,CurrentSite,i)
                if(NearestSite(NVertex,Temp,5-i)/=CurrentSite) then
                    Temp=NextSite(NVertex,Temp)
                endif
                do while(KinkTime(NVertex,Temp)<NewSiteTime)
                    NearestSite(NVertex,Temp,5-i)=PSite
                    Temp=NextSite(NVertex,Temp)
                enddo
                if(KinkTime(NVertex,Temp)==NewSiteTime) then
                    NearestSite(CurrentVertex,CurrentSite,i)=Temp
                else
                    NearestSite(CurrentVertex,CurrentSite,i)=PrevSite(NVertex,Temp)
                endif
            enddo
        endif
        KinkTime(CurrentVertex,CurrentSite)=NewSiteTime
        return
    end SUBROUTINE move_site

    double precision FUNCTION compare(Time1,Time2,IsPast)
        implicit none
        double precision :: Time1,Time2
        integer :: IsPast
        if(IsPast==1) then
            if(Time1>Time2) then
                compare=Time1
            else
                compare=Time2
            endif            
        else 
            if(Time1<Time2) then
                compare=Time1
            else
                compare=Time2
            endif
        endif
        return
    END FUNCTION compare

    SUBROUTINE change_nearest(Vertex,EndSite,OldSite,NewSite)
        implicit none
        integer ::Vertex,OldSite,NewSite,EndSite,NVertex,tarSite
        integer ::i
        double precision ::Time
        Time=KinkTime(Vertex,NewSite)
        do i=1,4
            NVertex=Ngs(5-i,Vertex)
            tarSite=NearestSite(Vertex,EndSite,5-i)
            if(NearestSite(NVertex,tarSite,i)==EndSite) then
                tarSite=PrevSite(NVertex,tarSite)
            endif
            do while(KinkTime(NVertex,tarSite)>=Time .and. &
                &  NearestSite(NVertex,tarSite,i)==OldSite .and. tarSite/=1)
                NearestSite(NVertex,tarSite,i)=NewSite
                tarSite=PrevSite(NVertex,tarSite)
            enddo
        enddo
        return
    END SUBROUTINE

    integer FUNCTION nearest_site(CurrentVertex,CurrentSite,IsPast)
        implicit none
        integer :: CurrentVertex,CurrentSite
        integer ::IsPast
        if(IsPast/=1) then
            nearest_site=NextSite(CurrentVertex,CurrentSite)
        else
            nearest_site=PrevSite(CurrentVertex,CurrentSite)
        endif 
        return       
    END FUNCTION

    SUBROUTINE change_kink_weight()
        implicit none
        integer :: i
        integer :: Vertex,Site
        double precision :: weight,Time
        do i=1,ListIndex-1
            Vertex=KinkWeightChangeVertex(i)
            Site=KinkWeightChangeSite(i)
            KinkWeight(Vertex,Site)=KinkWeightChange(i)
            KinkWeight(Ngs(NeighVertexNum(Vertex,Site),Vertex),NeighSite(Vertex,Site))=KinkWeightChange(i)
        enddo
    end SUBROUTINE change_kink_weight

    double precision FUNCTION two_plaquette_weight(Vertex1,Site1,Vertex2,Site2,num,Time)
    !       k------(1)-------m
    !       |       |        |
    !       |       |        |
    !       l------(2)-------n
    !if (k,l) or (m,n) is also a kink, two_plaquette_weight will not change!!!
        implicit none
        integer :: Vertex1,Vertex2,num,Site1,Site2
        double precision :: Time
        integer :: k,l,m,n
        integer :: neigh1,neigh2
        if(num==2 .or. num==3) then
            neigh1=1;neigh2=4
        else
            neigh1=2;neigh2=3
        endif
        k=get_state(Ngs(neigh1,Vertex1),NearestSite(Vertex1,Site1,neigh1),Time)
        m=get_state(Ngs(neigh2,Vertex1),NearestSite(Vertex1,Site1,neigh2),Time)
        l=get_state(Ngs(neigh1,Vertex2),NearestSite(Vertex2,Site2,neigh1),Time)
        n=get_state(Ngs(neigh2,Vertex2),NearestSite(Vertex2,Site2,neigh2),Time)
        two_plaquette_weight=Q*(k*l+m*n)
        return
    end FUNCTION

    integer FUNCTION get_state(Vertex,BeginSite,Time)
        implicit none
        integer :: Vertex,BeginSite,TarSite
        double precision :: Time
        TarSite=BeginSite
        do while(KinkTime(Vertex,TarSite)<Time)
            TarSite=NextSite(Vertex,TarSite)
        enddo
        get_state=-SegmentState(Vertex,TarSite)
        return
    end FUNCTION get_state

    double precision FUNCTION single_energy_kinkweight(CurVertex,BeginSite,BeginTime,EndTime,KinkFactor,IsAdd) 
        implicit none
        integer :: CurVertex,BeginSite,State,IsAdd
        integer :: NSite,NVertex,Vertex(3),Site(3),num,neigh1,neigh2
        double precision :: Energy,BeginTime,EndTime,KinkFactor
        Energy=0.0
        if(.not. IsAdd) KinkFactor=1.0
        State=SegmentState(CurVertex,BeginSite)
        neigh1=1;neigh2=4
        num=2
        Vertex(1)=Ngs(num,CurVertex)
        Site(1)=NearestSite(CurVertex,BeginSite,num)
        Vertex(2)=Ngs(neigh1,CurVertex);Site(2)=NearestSite(CurVertex,BeginSite,neigh1)
        Vertex(3)=Ngs(num,Vertex(2));Site(3)=NearestSite(Vertex(2),Site(2),num)
        Energy=Energy+plaquette_one_bond_kinkweight(State,Vertex,Site,5-neigh1,5-num,BeginTime,EndTime,KinkFactor,IsAdd)

        Vertex(2)=Ngs(neigh2,CurVertex);Site(2)=NearestSite(CurVertex,BeginSite,neigh2)
        Vertex(3)=Ngs(num,Vertex(2));Site(3)=NearestSite(Vertex(2),Site(2),num)
        Energy=Energy+plaquette_one_bond_kinkweight(State,Vertex,Site,5-neigh2,5-num,BeginTime,EndTime,KinkFactor,1)
        
        num=3
        Vertex(1)=Ngs(num,CurVertex)
        Site(1)=NearestSite(CurVertex,BeginSite,num)
        Vertex(2)=Ngs(neigh1,CurVertex);Site(2)=NearestSite(CurVertex,BeginSite,neigh1)
        Vertex(3)=Ngs(num,Vertex(2));Site(3)=NearestSite(Vertex(2),Site(2),num)
        Energy=Energy+plaquette_one_bond_kinkweight(State,Vertex,Site,5-neigh1,5-num,BeginTime,EndTime,KinkFactor,1)

        Vertex(2)=Ngs(neigh2,CurVertex);Site(2)=NearestSite(CurVertex,BeginSite,neigh2)
        Vertex(3)=Ngs(num,Vertex(2));Site(3)=NearestSite(Vertex(2),Site(2),num)
        Energy=Energy+plaquette_one_bond_kinkweight(State,Vertex,Site,5-neigh2,5-num,BeginTime,EndTime,KinkFactor,1)

        single_energy_kinkweight=-2.0*State*Energy
    end FUNCTION single_energy_kinkweight

    double precision FUNCTION double_energy_kinkweight(Vertex1,Site1,Vertex2,Site2,NeighNum,BeginTime,EndTime,KinkFactor)
        implicit none
        integer :: Vertex1,Vertex2,NeighNum,Site1,Site2,num
        double precision :: BeginTime,EndTime,Energy1,Energy2,KinkFactor
        integer :: neigh1,neigh2,Vertex(3),Site(3),State
        Energy1=0.0;Energy2=0.0;KinkFactor=1.0
        if(NeighNum==2 .or. NeighNum==3) then
            neigh1=1;neigh2=4
        else
            neigh1=2;neigh2=3
        endif
        num=5-NeighNum
        State=SegmentState(Vertex1,Site1)
        Vertex(1)=Ngs(num,Vertex1)
        Site(1)=NearestSite(Vertex1,Site1,num)
        Vertex(2)=Ngs(neigh1,Vertex1);Site(2)=NearestSite(Vertex1,Site1,neigh1)
        Vertex(3)=Ngs(num,Vertex(2));Site(3)=NearestSite(Vertex(2),Site(2),num)
        Energy1=Energy1+plaquette_two_bond_kinkweight(State,Vertex,Site,5-neigh1,5-num,BeginTime,EndTime,KinkFactor,0)

        Vertex(2)=Ngs(neigh2,Vertex1);Site(2)=NearestSite(Vertex1,Site1,neigh2)
        Vertex(3)=Ngs(num,Vertex(2));Site(3)=NearestSite(Vertex(2),Site(2),num)
        Energy1=Energy1+plaquette_two_bond_kinkweight(State,Vertex,Site,5-neigh2,5-num,BeginTime,EndTime,KinkFactor,1)

        num=NeighNum
        State=SegmentState(Vertex2,Site2)
        Vertex(1)=Ngs(num,Vertex2)
        Site(1)=NearestSite(Vertex2,Site2,num)
        Vertex(2)=Ngs(neigh1,Vertex2);Site(2)=NearestSite(Vertex2,Site2,neigh1)
        Vertex(3)=Ngs(num,Vertex(2));Site(3)=NearestSite(Vertex(2),Site(2),num)
        Energy2=Energy2+plaquette_two_bond_kinkweight(State,Vertex,Site,5-neigh1,5-num,BeginTime,EndTime,KinkFactor,1)

        Vertex(2)=Ngs(neigh2,Vertex2);Site(2)=NearestSite(Vertex2,Site2,neigh2)
        Vertex(3)=Ngs(num,Vertex(2));Site(3)=NearestSite(Vertex(2),Site(2),num)
        Energy2=Energy2+plaquette_two_bond_kinkweight(State,Vertex,Site,5-neigh2,5-num,BeginTime,EndTime,KinkFactor,1)

        double_energy_kinkweight=-2.0*(SegmentState(Vertex1,Site1)*Energy1+SegmentState(Vertex2,Site2)*Energy2)
        return
    end FUNCTION double_energy_kinkweight

    double precision FUNCTION delta_E(NVertex,NSite,BeginTime,EndTime)
        integer :: NVertex,NSite,TarSite
        double precision :: BeginTime,EndTime,TarTime,PrevTime
        double precision :: Energy
        Energy=0
        if(NSite==2) then
            TarSite=2
        else
            TarSite=NextSite(NVertex,NSite)
            do while(KinkTime(NVertex,TarSite)<BeginTime .and. TarSite/=2)
                TarSite=NextSite(NVertex,TarSite)
            enddo
        endif
        TarTime=KinkTime(NVertex,TarSite)
        PrevTime=BeginTime
        do while(TarTime<EndTime .and. TarSite/=2)
            Energy=Energy-SegmentState(NVertex,TarSite)*(TarTime-PrevTime)
            PrevTime=TarTime
            TarSite=NextSite(NVertex,TarSite)
            TarTime=KinkTime(NVertex,TarSite)
        enddo
        Energy=Energy-SegmentState(NVertex,TarSite) &
            &  *(EndTime-PrevTime)
        delta_E=JJQ*Energy
        !For ferromagnetic case, there will be a minus sign here
    END FUNCTION delta_E

    double precision FUNCTION two_plaquette_energy(CurrentVertex,num,BeginSite,BeginTime,EndTime)
        implicit none
        integer :: CurrentVertex,num,BeginSite
        double precision :: BeginTime,EndTime,Energy
        integer :: neigh1,neigh2,Vertex(3),Site(3)
        Energy=0.0
        Vertex(1)=Ngs(num,CurrentVertex)
        Site(1)=NearestSite(CurrentVertex,BeginSite,num)
        if(num==2 .or. num==3) then
            neigh1=1;neigh2=4
        else
            neigh1=2;neigh2=3
        endif
        Vertex(2)=Ngs(neigh1,CurrentVertex);Site(2)=NearestSite(CurrentVertex,BeginSite,neigh1)
        Vertex(3)=Ngs(num,Vertex(2));Site(3)=NearestSite(Vertex(2),Site(2),num)
        Energy=Energy+plaquette_energy(Vertex,Site,BeginTime,EndTime)

        Vertex(2)=Ngs(neigh2,CurrentVertex);Site(2)=NearestSite(CurrentVertex,BeginSite,neigh2)
        Vertex(3)=Ngs(num,Vertex(2));Site(3)=NearestSite(Vertex(2),Site(2),num)
        Energy=Energy+plaquette_energy(Vertex,Site,BeginTime,EndTime)
        two_plaquette_energy=Energy
        return
    end FUNCTION two_plaquette_energy

    double precision FUNCTION plaquette_energy(Vertex,Site,BeginTime,EndTime)
        implicit none
        integer :: Vertex(3),i,j
        integer :: Site(3),TarSite(3)
        double precision :: BeginTime,EndTime,TarTime,PrevTime
        double precision :: Energy
        Energy=0.0
        do i=1,3
            if(Site(i)==2) then
                TarSite(i)=2
            else
                TarSite(i)=NextSite(Vertex(i),Site(i))
                do while(KinkTime(Vertex(i),TarSite(i))<BeginTime .and. TarSite(i)/=2)
                    TarSite(i)=NextSite(Vertex(i),TarSite(i))
                enddo
            endif
        enddo
        TarTime=KinkTime(Vertex(1),TarSite(1));j=1
        do i=2,3
            if(TarTime>KinkTime(Vertex(i),TarSite(i))) then
                TarTime=KinkTime(Vertex(i),TarSite(i))
                j=i
            endif
        enddo
        PrevTime=BeginTime
        do while(TarTime<EndTime .and. TarSite(j)/=2)
            Energy=Energy-SegmentState(Vertex(1),TarSite(1))* &
               & SegmentState(Vertex(2),TarSite(2))*          &
               & SegmentState(Vertex(3),TarSite(3))*(TarTime-PrevTime)
            PrevTime=TarTime
            TarSite(j)=NextSite(Vertex(j),TarSite(j))
            TarTime=KinkTime(Vertex(1),TarSite(1));j=1
            do i=2,3
                if(TarTime>KinkTime(Vertex(i),TarSite(i))) then
                    TarTime=KinkTime(Vertex(i),TarSite(i))
                    j=i
                endif
            enddo
        enddo
        Energy=Energy-SegmentState(Vertex(1),TarSite(1))* &
           & SegmentState(Vertex(2),TarSite(2))*          &
           & SegmentState(Vertex(3),TarSite(3))*(EndTime-PrevTime)
        plaquette_energy=-2*Q*Energy
        return
    end FUNCTION plaquette_energy

    double precision FUNCTION plaquette_one_bond_kinkweight(State,Vertex,Site,num1,num2,BeginTime,EndTime,KinkFactor,IsAdd)
        implicit none
        integer :: Vertex(3),i,j,num1,num2,State,num
        integer :: Site(3),TarSite(3),IsAdd
        double precision :: BeginTime,EndTime,TarTime,PrevTime
        double precision :: Energy,KinkFactor,weight,final_weight
        Energy=0.0
        if(.not. IsAdd) ListIndex=1 
        do i=1,3
            TarSite(i)=NextSite(Vertex(i),Site(i))
            do while(KinkTime(Vertex(i),TarSite(i))<BeginTime .and. TarSite(i)/=2)
                TarSite(i)=NextSite(Vertex(i),TarSite(i))
            enddo
        enddo
        TarTime=KinkTime(Vertex(3),TarSite(3));j=3
        do i=1,2
            if(TarTime>KinkTime(Vertex(i),TarSite(i))) then
                TarTime=KinkTime(Vertex(i),TarSite(i))
                j=i
            endif
        enddo
        PrevTime=BeginTime
        do while(TarTime<EndTime .and. TarSite(j)/=2)
            Energy=Energy-(JJQ/2.0*(SegmentState(Vertex(1),TarSite(1))+SegmentState(Vertex(2),TarSite(2)))-2.0*Q*  &
               & SegmentState(Vertex(1),TarSite(1))*          &
               & SegmentState(Vertex(2),TarSite(2))*          &
               & SegmentState(Vertex(3),TarSite(3)))*(TarTime-PrevTime)
            if(j==3 .and. PairedSite(Vertex(3),TarSite(3))==-1) then
                if(Vertex(3)/=MashaVertex .or. TarSite(3)/=MashaSite) then
                    num=NeighVertexNum(Vertex(3),TarSite(3))
                    if(num==num1) then
                        weight=KinkWeight(Vertex(3),TarSite(3))
                        final_weight=weight-2*Q*State*SegmentState(Vertex(2),TarSite(2))
                        KinkFactor=KinkFactor*final_weight/weight
                        KinkWeightChangeVertex(ListIndex)=Vertex(3)
                        KinkWeightChangeSite(ListIndex)=TarSite(3)
                        KinkWeightChange(ListIndex)=final_weight
                        ListIndex=ListIndex+1
                    elseif(num==num2) then
                        weight=KinkWeight(Vertex(3),TarSite(3))
                        final_weight=weight-2*Q*State*SegmentState(Vertex(1),TarSite(1))
                        KinkFactor=KinkFactor*final_weight/weight
                        KinkWeightChangeVertex(ListIndex)=Vertex(3)
                        KinkWeightChangeSite(ListIndex)=TarSite(3)
                        KinkWeightChange(ListIndex)=final_weight
                        ListIndex=ListIndex+1
                    endif
                endif
            endif
            PrevTime=TarTime
            TarSite(j)=NextSite(Vertex(j),TarSite(j))
            TarTime=KinkTime(Vertex(3),TarSite(3));j=3
            do i=1,2
                if(TarTime>KinkTime(Vertex(i),TarSite(i))) then
                    TarTime=KinkTime(Vertex(i),TarSite(i))
                    j=i
                endif
            enddo
        enddo
        Energy=Energy-(JJQ/2.0*(SegmentState(Vertex(1),TarSite(1))+SegmentState(Vertex(2),TarSite(2)))-2.0*Q*  &
           & SegmentState(Vertex(1),TarSite(1))*          &
           & SegmentState(Vertex(2),TarSite(2))*          &
           & SegmentState(Vertex(3),TarSite(3)))*(EndTime-PrevTime)
        plaquette_one_bond_kinkweight=Energy
        return
    end FUNCTION plaquette_one_bond_kinkweight

    double precision FUNCTION plaquette_two_bond_kinkweight(State,Vertex,Site,num1,num2,BeginTime,EndTime,KinkFactor,IsAdd)
        implicit none
        integer :: Vertex(3),i,j,num1,num2,num
        integer :: Site(3),TarSite(3),IsAdd,State
        double precision :: BeginTime,EndTime,TarTime,PrevTime
        double precision :: Energy,KinkFactor,weight,final_weight
        Energy=0.0
        if(.not. IsAdd) ListIndex=1 
        do i=1,3
            TarSite(i)=NextSite(Vertex(i),Site(i))
            do while(KinkTime(Vertex(i),TarSite(i))<BeginTime .and. TarSite(i)/=2)
                TarSite(i)=NextSite(Vertex(i),TarSite(i))
            enddo
        enddo
        TarTime=KinkTime(Vertex(3),TarSite(3));j=3
        do i=1,2
            if(TarTime>KinkTime(Vertex(i),TarSite(i))) then
                TarTime=KinkTime(Vertex(i),TarSite(i))
                j=i
            endif
        enddo
        PrevTime=BeginTime
        do while(TarTime<EndTime .and. TarSite(j)/=2)
            Energy=Energy-(JJQ*(SegmentState(Vertex(1),TarSite(1))/2.0+SegmentState(Vertex(2),TarSite(2)))-2.0*Q*  &
               & SegmentState(Vertex(1),TarSite(1))*          &
               & SegmentState(Vertex(2),TarSite(2))*          &
               & SegmentState(Vertex(3),TarSite(3)))*(TarTime-PrevTime)
            if(j==3 .and. PairedSite(Vertex(3),TarSite(3))==-1) then
                if(Vertex(3)/=MashaVertex .or. TarSite(3)/=MashaSite) then
                    num=NeighVertexNum(Vertex(3),TarSite(3))
                    if(num==num1) then
                        weight=KinkWeight(Vertex(3),TarSite(3))
                        final_weight=weight-2*Q*State*SegmentState(Vertex(2),TarSite(2))
                        KinkFactor=KinkFactor*final_weight/weight
                        KinkWeightChangeVertex(ListIndex)=Vertex(3)
                        KinkWeightChangeSite(ListIndex)=TarSite(3)
                        KinkWeightChange(ListIndex)=final_weight
                        ListIndex=ListIndex+1
                    elseif(num==num2) then
                        weight=KinkWeight(Vertex(3),TarSite(3))
                        final_weight=weight-2*Q*State*SegmentState(Vertex(1),TarSite(1))
                        KinkFactor=KinkFactor*final_weight/weight
                        KinkWeightChangeVertex(ListIndex)=Vertex(3)
                        KinkWeightChangeSite(ListIndex)=TarSite(3)
                        KinkWeightChange(ListIndex)=final_weight
                        ListIndex=ListIndex+1
                    endif
                endif
            endif
            PrevTime=TarTime
            TarSite(j)=NextSite(Vertex(j),TarSite(j))
            TarTime=KinkTime(Vertex(3),TarSite(3));j=3
            do i=1,2
                if(TarTime>KinkTime(Vertex(i),TarSite(i))) then
                    TarTime=KinkTime(Vertex(i),TarSite(i))
                    j=i
                endif
            enddo
        enddo
        Energy=Energy-(JJQ*(SegmentState(Vertex(1),TarSite(1))/2.0+SegmentState(Vertex(2),TarSite(2)))-2.0*Q*  &
           & SegmentState(Vertex(1),TarSite(1))*          &
           & SegmentState(Vertex(2),TarSite(2))*          &
           & SegmentState(Vertex(3),TarSite(3)))*(EndTime-PrevTime)
        plaquette_two_bond_kinkweight=Energy
        return
    end FUNCTION plaquette_two_bond_kinkweight
  !===================================================================


  !==============Measurement =========================================
  SUBROUTINE measure
    implicit none
    integer :: i
    double precision :: Vbs_x,Vbs_y,TotalKinks

    !--Observables definition ----------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    if(Wt==0 .and. Wx==0 .and. Wy==0) then
        W0=1;W1=0
    else
        W0=0;W1=1
    endif
    Quan( 1)=EnergyCheck
    Quan( 2)=PairedKinkNum/Beta
    Quan( 3)=kink_number(TotalKinks)/Beta
    Quan( 4)=(Wx)**2
    Quan( 5)=(Wy)**2
    Quan( 6)=(Wt)**2
    Quan( 7)=Quan(4)+Quan(5)+Quan(6)
    Quan( 8)=W1
    Quan( 9)=W0
    Quan(10)=SM()**2
    !call vbs(Vbs_x,Vbs_y)
    !Quan(11)=Vbs_x**2+Vbs_y**2
    !call vbs1(Vbs_x,Vbs_y)
    !Quan(12)=Vbs_x**2+Vbs_y**2
    call SU2_vbs(Vbs_x,Vbs_y,Beta)
    Quan(11)=(Vbs_x+Vbs_y)/2.0
    call SU2_vbs(Vbs_x,Vbs_y,Beta/Lx)
    Quan(12)=(Vbs_x+Vbs_y)/2.0
    if(W0==1) then
        Quan(13)=Quan(10)
        Quan(14)=0.0
        Quan(15)=Quan(11)
        Quan(16)=0.0
    else
        Quan(13)=0.0
        Quan(14)=Quan(10)
        Quan(15)=0.0
        Quan(16)=Quan(11)
    endif
    call vbs1(Vbs_x,Vbs_y)
    Quan(17)=(Vbs_x**2+Vbs_y**2)/2.0
    call vbs2(Vbs_x,Vbs_y)
    Quan(18)=(Vbs_x**2+Vbs_y**2)/2.0
    Quan(19)=Quan(10)**2
    Quan(20)=Quan(11)**2
    Quan(21)=Quan(12)**2
    Quan(22)=potential_J()
    Quan(23)=Quan(7)*Quan(22)
    Quan(24)=Quan(10)
    Quan(25)=Susceptibility_total()
  END SUBROUTINE measure

  double precision FUNCTION potential_energy()
    implicit none
    double precision :: BeginTime,EndTime,temp1,temp2
    integer :: i,j,BeginSite,NVertex,NSite,EndSite
    integer :: Vertex(3),Site(3)
    double precision :: Energy1,Energy2
    Energy1=0.0;Energy2=0.0
    do i=1,Vol
        BeginSite=1
        EndSite=NextSite(i,1)
        temp1=0.0
        temp2=0.0
        Vertex(1)=Ngs(1,i)
        Vertex(2)=Ngs(2,i)
        Vertex(3)=Ngs(1,Vertex(2))
        Site(:)=1
        do while(BeginSite/=2)
            BeginTime=KinkTime(i,BeginSite)
            EndTime=KinkTime(i,EndSite)
            do j=1,2
                NSite=NearestSite(i,BeginSite,j)
                NVertex=Ngs(j,i)
                temp1=temp1+SegmentState(i,BeginSite)*delta_E(NVertex,NSite,BeginTime,EndTime)
            enddo
            temp2=temp2+SegmentState(i,BeginSite)*plaquette_energy(Vertex,Site,BeginTime,EndTime)
            Site(1)=NearestSite(i,EndSite,1)
            Site(2)=NearestSite(i,EndSite,2)
            Site(3)=NearestSite(Vertex(2),Site(2),1)
            BeginSite=EndSite
            EndSite=NextSite(i,EndSite)
        enddo
        Energy1=Energy1+temp1
        Energy2=Energy2+temp2
    enddo
    potential_energy=Energy1+Energy2
    !potential_energy=Energy1
  end FUNCTION

  double precision FUNCTION potential_J()
    implicit none
    double precision :: BeginTime,EndTime,temp1,temp2
    integer :: i,j,BeginSite,NVertex,NSite,EndSite
    integer :: Vertex(3),Site(3)
    double precision :: Energy1,Energy2
    Energy1=0.0;Energy2=0.0
    do i=1,Vol
        Energy1=Energy1+SegmentState(i,1)*(SegmentState(Ngs(1,i),1)+SegmentState(Ngs(2,i),1))
    enddo
    potential_J=Energy1
    !potential_J=JJ*Energy1/Beta
  end FUNCTION

  double precision FUNCTION kink_number(TotalKinks)
    implicit none
    integer :: i
    double precision :: TotalKinks
    kink_number=0.0
    TotalKinks=0.0
    do i=1,Vol
        kink_number=kink_number+SegmentNum(i)
    enddo
    TotalKinks=(kink_number-Vol)/2
    kink_number=(kink_number-Vol-4*PairedKinkNum)/2
  end FUNCTION kink_number

  SUBROUTINE winding_number()
      implicit none
      integer dx,dy,dt
      integer Vertex,i
      integer Site
      dx=0;dy=0;dt=0
      do Vertex=1,Vol
          dt=dt+SegmentState(Vertex,1)
      enddo
      do i=1,Lx
          Vertex=i
          Site=NextSite(Vertex,1)
          do while(Site/=2)
              if(NeighVertexNum(Vertex,Site)==2) then
                  dy=dy-SegmentState(Vertex,Site)
              endif
              Site=NextSite(Vertex,Site)
          enddo
      enddo
      do i=1,Ly
          Vertex=(i-1)*Lx+1
          Site=NextSite(Vertex,1)
          do while(Site/=2)
              if(NeighVertexNum(Vertex,Site)==1) then
                  dx=dx-SegmentState(Vertex,Site)
              endif
              Site=NextSite(Vertex,Site)
          enddo
      enddo
      WindX=dx
      WindY=dy
      WindT=dt/2
  end SUBROUTINE winding_number

  SUBROUTINE vbs(Vbs_x,Vbs_y)
      implicit none
      double precision :: Vbs_x,Vbs_y
      integer Vertex,Site,X,Y
      integer KinkNum_x,KinkNum_y
      Vbs_x=0.0;Vbs_y=0.0
      do Y=0,Ly-1
          do X=1,Lx
              Vertex=Y*Lx+X
              Site=1
              KinkNum_x=0;KinkNum_y=0
              do while(Site/=2)
                  if(NeighVertexNum(Vertex,Site)==2) then
                      KinkNum_y=KinkNum_y+1
                  endif
                  if(NeighVertexNum(Vertex,Site)==1) then
                      KinkNum_x=KinkNum_x+1
                  endif
                  Site=NextSite(Vertex,Site)
              enddo
              Vbs_x=Vbs_x+KinkNum_x*(-1)**X
              Vbs_y=Vbs_y+KinkNum_y*(-1)**(Y+1)
          enddo
      enddo
      Vbs_x=Vbs_x/Vol/beta
      Vbs_y=Vbs_y/Vol/beta
  end SUBROUTINE vbs

  double precision FUNCTION SM()
    implicit none
    double precision :: M
    integer :: Vertex,X,Y
    M=0.0
    do Y=0,Ly-1
        do X=1,Lx
            Vertex=Y*Lx+X
            M=M+SegmentState(Vertex,1)*(-1)**(X+Y+1)
        enddo
    enddo
    SM=M/Vol/2.0
  end FUNCTION SM

  SUBROUTINE vbs1(Vbs_x,Vbs_y)
      implicit none
      double precision :: Vbs_x,Vbs_y
      integer Vertex,X,Y
      Vbs_x=0.0;Vbs_y=0.0
      do Y=0,Ly-1
          do X=1,Lx
              Vertex=Y*Lx+X
              Vbs_x=Vbs_x+SegmentState(Vertex,1)*SegmentState(Ngs(1,Vertex),1)*(-1)**X
              Vbs_y=Vbs_y+SegmentState(Vertex,1)*SegmentState(Ngs(2,Vertex),1)*(-1)**(Y-1)
          enddo
      enddo
      Vbs_x=Vbs_x/Vol/4.0
      Vbs_y=Vbs_y/Vol/4.0
  end SUBROUTINE

  SUBROUTINE vbs2(Vbs_x,Vbs_y)
      implicit none
      double precision :: Vbs_x,Vbs_y,temp1,temp2,BeginTime,EndTime
      integer Vertex,X,Y,NSite,NVertex,BeginSite,EndSite
      Vbs_x=0.0;Vbs_y=0.0
      do Y=0,Ly-1
          do X=1,Lx
                Vertex=Y*Lx+X
                temp1=0.0
                temp2=0.0
                BeginSite=1
                EndSite=NextSite(Vertex,1)
                do while(BeginSite/=2)
                    BeginTime=KinkTime(Vertex,BeginSite)
                    EndTime=KinkTime(Vertex,EndSite)
                    NSite=NearestSite(Vertex,BeginSite,1)
                    NVertex=Ngs(1,Vertex)
                    temp1=temp1+SegmentState(Vertex,BeginSite)*delta_E(NVertex,NSite,BeginTime,EndTime)
                    NSite=NearestSite(Vertex,BeginSite,2)
                    NVertex=Ngs(2,Vertex)
                    temp2=temp2+SegmentState(Vertex,BeginSite)*delta_E(NVertex,NSite,BeginTime,EndTime)
                    BeginSite=EndSite
                    EndSite=NextSite(Vertex,EndSite)
                enddo
                Vbs_x=Vbs_x+temp1/JJQ*(-1)**X
                Vbs_y=Vbs_y+temp2/JJQ*(-1)**(Y-1)
          enddo
      enddo
      Vbs_x=Vbs_x/Vol/beta/4.0
      Vbs_y=Vbs_y/Vol/beta/4.0
  end SUBROUTINE 

  SUBROUTINE SU2_vbs(Vbs_x,Vbs_y,T)
      implicit none
      double precision :: Vbs_x,Vbs_y,temp1,temp2,temp3,temp4,BeginTime,EndTime,T
      integer Vertex,X,Y,NSite,NVertex,BeginSite,EndSite
      Vbs_x=0.0;Vbs_y=0.0
      temp3=0.0;temp4=0.0
      do Y=0,Ly-1
          do X=1,Lx
                Vertex=Y*Lx+X
                temp1=0.0
                temp2=0.0
                BeginSite=1
                EndSite=NextSite(Vertex,1)
                BeginTime=KinkTime(Vertex,BeginSite)
                EndTime=KinkTime(Vertex,EndSite)
                do while(BeginTime<T)
                    if(EndTime>=T)EndTime=T
                    NSite=NearestSite(Vertex,BeginSite,1)
                    NVertex=Ngs(1,Vertex)
                    temp1=temp1+SegmentState(Vertex,BeginSite)*delta_E(NVertex,NSite,BeginTime,EndTime)/JJQ
                    NSite=NearestSite(Vertex,BeginSite,2)
                    NVertex=Ngs(2,Vertex)
                    temp2=temp2+SegmentState(Vertex,BeginSite)*delta_E(NVertex,NSite,BeginTime,EndTime)/JJQ
                    if(BeginSite/=1 .and. BeginSite/=2 .and. PairedSite(Vertex,BeginSite)==-1) then
                        if(NeighVertexNum(Vertex,BeginSite)==1) then
                            temp1=temp1-1/KinkWeight(Vertex,BeginSite)
                            temp3=temp3+1/KinkWeight(Vertex,BeginSite)**2
                        elseif(NeighVertexNum(Vertex,BeginSite)==2) then
                            temp2=temp2-1/KinkWeight(Vertex,BeginSite)
                            temp4=temp4+1/KinkWeight(Vertex,BeginSite)**2
                        endif
                    endif
                    BeginSite=EndSite
                    EndSite=NextSite(Vertex,EndSite)
                    BeginTime=KinkTime(Vertex,BeginSite)
                    EndTime=KinkTime(Vertex,EndSite)
                enddo
                Vbs_x=Vbs_x+temp1*(-1)**X
                Vbs_y=Vbs_y+temp2*(-1)**(Y-1)
          enddo
      enddo
      Vbs_x=(Vbs_x**2-temp3)/Vol/Vol/T/T/16.0
      Vbs_y=(Vbs_y**2-temp4)/Vol/Vol/T/T/16.0
  end SUBROUTINE 
  
  double precision FUNCTION Susceptibility_total()
    implicit none
    double precision :: M,Sz,BeginTime,EndTime
    integer :: Vertex,X,Y
    integer :: BeginSite,EndSite
    M=0.0
    do Y=0,Ly-1
        do X=1,Lx
            Sz=0.0
            Vertex=Y*Lx+X
            BeginSite=1
            EndSite=NextSite(Vertex,BeginSite)
            do while(BeginSite/=2)
                BeginTime=KinkTime(Vertex,BeginSite)
                EndTime=KinkTime(Vertex,EndSite)
                Sz=Sz+SegmentState(Vertex,BeginSite)*(EndTime-BeginTime)
                BeginSite=EndSite
                EndSite=NextSite(Vertex,EndSite)
            enddo
            M=M+Sz*(-1)**(X+Y+1)
        enddo
    enddo
    Susceptibility_total=M**2/beta/Vol/4.0
  end FUNCTION Susceptibility_total

  !==============Calculate composite observables =====================
  !! THIS IS PROJECT-INDEPENDENT 
  !! call in 'stat_alan'
  SUBROUTINE cal_Obs_comp
    implicit none
    integer :: k
    double precision  :: nor

    !-- calculate the average ----------------------------------------
    Ave(NObs_b+1:NObs) = 0.d0

    !-- Q1=<C1>^2/<C1^2> ---------------------------------------------
    Ave(NObs_b+1)=(Ave(4)+Ave(5))/Ave(6)/2
    Ave(NObs_b+2)=Ave(8)/Ave(9)
    Ave(NObs_b+3)=Ave(13)/Ave(14)
    Ave(NObs_b+4)=Ave(15)/Ave(16)
    Ave(NObs_b+5)=Ave(17)/Ave(18)
    Ave(NObs_b+6)=5*(1-Ave(19)/Ave(10)**2/3.0)/2.0
    Ave(NObs_b+7)=1-Ave(20)/Ave(11)**2/3.0
    Ave(NObs_b+8)=1-Ave(21)/Ave(12)**2/3.0
    Ave(Nobs_b+9)=Ave(23)-Ave(7)*Ave(22)
    Ave(Nobs_b+10)=Ave(25)/Ave(24)
    Ave(Nobs_b+11)=Ave(11)/Ave(12)
    Ave(Nobs_b+12)=Ave(18)/Ave(17)

    !-- Obs(j,k) series ----------------------------------------------
    do k=1,NBlck
        if(abs(Obs(6,k))>eps) then
            Obs(NObs_b+1,k)=(Obs(4,k)+Obs(5,k))/Obs(6,k)/2
        else
            Obs(NObs_b+1,k)=0.0
        endif
    enddo
    Obs(NObs_b+2,:)=Obs(8,:)/Obs(9,:)
    Obs(NObs_b+3,:)=Obs(13,:)/Obs(14,:)
    Obs(NObs_b+4,:)=Obs(15,:)/Obs(16,:)
    Obs(NObs_b+5,:)=Obs(17,:)/Obs(18,:)
    Obs(NObs_b+6,:)=5*(1-Obs(19,:)/Obs(10,:)**2/3.0)/2.0
    Obs(NObs_b+7,:)=1-Obs(20,:)/Obs(11,:)**2/3.0
    Obs(NObs_b+8,:)=1-Obs(21,:)/Obs(12,:)**2/3.0
    Obs(NObs_b+9,:)=Obs(23,:)-Obs(7,:)*Obs(22,:)
    Obs(NObs_b+10,:)=Obs(25,:)/Obs(24,:)
    Obs(NObs_b+11,:)=Obs(11,:)/Obs(12,:)
    Obs(NObs_b+12,:)=Obs(18,:)/Obs(17,:)

   return
  END SUBROUTINE cal_Obs_comp
  !===================================================================
  !*******************************************************************
  !            End of PROJECT-DEPENDENT part 
  !*******************************************************************


  !*******************************************************************
  !            Beginning of PROJECT-INDEPENDENT part 
  !*******************************************************************
  !============== Write to files =====================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE write2file 
    implicit none
    integer       :: i, j, k, Nwri

    open (8,file=file1, access='append') 
    write(8, *);          write(6,*)

    write(8,40) ident, Lx, beta, JJ*4.0, Q*16.0, Seed, rank, TotSamp, NBlck
    40 format(a8,i6,3f10.6,i12,i8,i8,i8)

    do j = 1, Nobs
      write(8,41) j, Ave(j), Dev(j), Cor(j)
      41 format(i4,2f25.15,f12.5)
      write(6,41) j, Ave(j), Dev(j), Cor(j)
    enddo

    do j=1, 7
      write(8,43) RunName(j),RunNum(j)/Number,RunNum(j)/RunNum0(j),RunNum0(j)/Number*7.0
      write(6,43) RunName(j),RunNum(j)/Number,RunNum(j)/RunNum0(j),RunNum0(j)/Number*7.0
    enddo
    43 format(a13,3f12.5)

    close(8)

    if(prt) return
    Nwri = NObs_b;                   if(Nwri>5) Nwri = 7
    write(6,*)
    do k = 1, NBlck
       write(6,42) k,(Obs(j,k),j=1,Nwri)
       42 format(i6,7f16.8) 
    end do

    return
  END SUBROUTINE write2file

  SUBROUTINE midwrite2file(iblck)
    implicit none
    integer :: iblck,i, j, start
    open (8,file=file2, access='append') 
    start=iblck-NSave+1
    do i=start,iblck
        write(8,*)
        write(8,50) ident, Lx, beta, JJ*4.0, Q*16.0, Seed,rank, NSamp, i
        50 format(a8,i6,3f10.6,i12,i8,i8,i8)
        do j = 1, NObs_b
          write(8,51) j,Obs(j,i)
          51 format(i4,f25.15)
        enddo
    enddo
    close(8)
  END SUBROUTINE midwrite2file
  !===================================================================

  SUBROUTINE readconfig()
    implicit none
    logical  :: IsPaired
    integer  :: i,j,k,Vertex,Site,num,State,Lxx,SegNum
    integer  :: CSite
    integer  :: PSite,NSite,NPSite,PVertexNum,NVertexNum,PVertex,NVertex,NPVertex
    integer  :: BeginSite
    double precision :: Time,CTime,BB,JJJ,QQ,Factor
    write(6,*) "read configuration"
    open(9,file=file3,form="binary")
    read(9) Lxx
    read(9) BB
    read(9) JJJ
    read(9) QQ 
    if(Lxx/=Lx .or. BB/=Beta .or. JJJ/=JJ .or. QQ/=Q) then
        open(18,file=file4,access="append")
          write(18,*) "Warning!"
          write(18,*) "Parameters are different in the old configuration!"
        close(18)
    endif
    Factor=Beta/BB
    do i=1,Vol
      read(9) SegmentState(i,1)
      SegmentState(i,2)=-SegmentState(i,1)
    enddo
    PairedKinkNum=0
    do i=1,Vol
       read(9) num
       BeginSite=1
       do j=1,num
           read(9) IsPaired
           read(9) NVertexNum
           read(9) Time
           Time=Time*Factor
           if(Time>Beta) Time=Beta-eps
           CSite=search_and_insert(i,Time)
           NVertex=Ngs(NVertexNum,i)
           NSite=search_and_insert(NVertex,Time)
           NeighVertexNum(i,CSite)=NVertexNum
           NeighVertexNum(NVertex,NSite)=5-NVertexNum
           NeighSite(i,CSite)=NSite
           NeighSite(NVertex,NSite)=CSite
           if(IsPaired==0) then
               PairedSite(i,CSite)=-1
               PairedSite(NVertex,NSite)=-1
           else
               if(NVertexNum==1) then
                   PVertexNum=2
               else
                   PVertexNum=1
               endif
               PVertex=Ngs(PVertexNum,i)
               PSite=search_and_insert(PVertex,Time)
               NPVertex=Ngs(NVertexNum,PVertex)
               NPSite=search_and_insert(NPVertex,Time)
               NeighVertexNum(PVertex,PSite)=NVertexNum
               NeighVertexNum(NPVertex,NPSite)=5-NVertexNum
               NeighSite(PVertex,PSite)=NPSite
               NeighSite(NPVertex,NPSite)=PSite
               PairedVertexNum(i,CSite)=PVertexNum
               PairedVertexNum(NVertex,NSite)=PVertexNum
               PairedVertexNum(PVertex,PSite)=5-PVertexNum
               PairedVertexNum(NPVertex,NPSite)=5-PVertexNum
               PairedSite(i,CSite)=PSite
               PairedSite(PVertex,PSite)=CSite
               PairedSite(NVertex,NSite)=NPSite
               PairedSite(NPVertex,NPSite)=NSite
               PairedKinkNum=PairedKinkNum+1
           endif
       enddo
    enddo
    close(9)
    do Vertex=1,Vol
        Site=NextSite(Vertex,1)
        do while(Site/=2)
            SegmentState(Vertex,Site)=-SegmentState(Vertex,PrevSite(Vertex,Site))
            Site=NextSite(Vertex,Site)
        enddo
    enddo
    do Vertex=1,Vol
        SegNum=1
        Site=NextSite(Vertex,1)
        do while(Site/=2)
            if(PairedSite(Vertex,Site)==-1) then
                Time=KinkTime(Vertex,Site)
                num=NeighVertexNum(Vertex,Site)
                KinkWeight(Vertex,Site)=JJQ-two_plaquette_weight(Vertex,Site,Ngs(num,Vertex),NeighSite(Vertex,Site),num,Time)
            endif
            SegNum=SegNum+1
            Site=NextSite(Vertex,Site)
        enddo
        SegmentNum(Vertex)=SegNum
    enddo
    write(*,*) "read configuration complete"
  end SUBROUTINE 

  SUBROUTINE saveconfig
    implicit none
    integer  :: i,j,k,num
    write(*,*) "Save configuration"
    open(9,file=file3,form="binary")
    write(9) Lx
    write(9) Beta
    write(9) JJ
    write(9) Q 
    do i=1,Vol
      write(9) SegmentState(i,1)
    enddo
    do i=1,Vol
      num=0
      j=NextSite(i,1)
      do while(j/=2)
          if(PairedSite(i,j)==-1) then
              if(NeighVertexNum(i,j)==1 .or. NeighVertexNum(i,j)==2) then
                  num=num+1
              endif
          else
              if(NeighVertexNum(i,j)==1 .and. PairedVertexNum(i,j)==2 .or. &
            &   NeighVertexNum(i,j)==2 .and. PairedVertexNum(i,j)==1) then
                  num=num+1
              endif
          endif
          j=NextSite(i,j)
      enddo
      write(9) num
      j=NextSite(i,1)
      do while(j/=2)
        if(PairedSite(i,j)==-1) then
            if(NeighVertexNum(i,j)==1 .or. NeighVertexNum(i,j)==2) then
                write(9) 0
                write(9) NeighVertexNum(i,j)
                write(9) KinkTime(i,j)
            endif
        else
            if(NeighVertexNum(i,j)==1 .and. PairedVertexNum(i,j)==2 .or. &
          &   NeighVertexNum(i,j)==2 .and. PairedVertexNum(i,j)==1) then
                write(9) 1
                write(9) NeighVertexNum(i,j)
                write(9) KinkTime(i,j)
            endif
        endif
        j=NextSite(i,j)
      enddo
    enddo
    close(9)
    write(*,*) "Save configuration complete"
  end SUBROUTINE

  integer FUNCTION search_and_insert(Vertex,Time)
      implicit none
      integer :: Vertex,Site
      double precision :: Time
      Site=NextSite(Vertex,1)
      do while(Site/=2 .and. KinkTime(Vertex,Site)<Time)
          Site=NextSite(Vertex,Site)
      enddo
      Site=PrevSite(Vertex,Site)
      search_and_insert=insert_site(Vertex,Site,Time)
      return
  end FUNCTION
  
  SUBROUTINE check_config()
    implicit none
    integer :: i,Site,Vertex,num,pnum,NVertex,KNum,PKNum
    double precision :: Energy
    open(18,file=file4,access="append")
    write(18,*) "checking..."
    state=0
    PKNum=0
    do Vertex=1,Vol
        KNum=0
        Site=NextSite(Vertex,1)
        do while(Site/=2)
            KNum=KNum+1
            if(PairedSite(Vertex,Site)/=-1)then
                PKNum=PKNum+1
            endif
            if(KinkTime(Vertex,Site)>=Beta)then
                write(18,*) Number, Vertex,Site
                state=error_print("Kink Time is Beta!")
            endif
            if(KinkTime(Vertex,Site)<=KinkTime(Vertex,PrevSite(Vertex,Site)))then
                write(18,*) Number, Vertex,Site
                write(18,*) KinkTime(Vertex,Site),KinkTime(Vertex,PrevSite(Vertex,Site))
                state=error_print("Kink Time order is wrong!")
            endif
            if(Vertex/=MashaVertex .or. Site/=MashaSite) then
                if(Vertex/=IraVertex .or. Site/=IraSite) then
                   num=NeighVertexNum(Vertex,Site)
                   if(PairedSite(Vertex,Site)==-1) then
                       if(num==1 .or. num==2) then
                           call check_kink_weight(Vertex,Site,num)
                       endif
                       pnum=0
                   else
                       pnum=PairedVertexNum(Vertex,Site)
                   endif
                   call check_kink_order(Vertex,Site,num,pnum)
                endif
            endif
            Site=NextSite(Vertex,Site)
            if(SegmentState(Vertex,Site)==0) then
                write(18,*) Number,Vertex,Site
                state=error_print("State=0!")
            endif
            if(SegmentState(Vertex,Site)/=-SegmentState(Vertex,PrevSite(Vertex,Site))) then
                write(18,*) Number,Vertex,Site,PrevSite(Vertex,Site)
                write(18,*) SegmentState(Vertex,Site),SegmentState(Vertex,PrevSite(Vertex,Site))
                state=error_print("State dismatch!")
            endif
        enddo
        if(KNum/=SegmentNum(Vertex))then
            write(*,*) "Segment Number wrong!",KNum,SegmentNum(Vertex),Vertex
        endif
    enddo
    if(PKNum/=PairedKinkNum)then
        write(*,*) "Paired Kink Number wrong!",PKNum,PairedKinkNum
    endif
    Energy=potential_energy()
    if(abs(EnergyCheck-Energy)/Energy>1e-10) then
        write(18,*) number,EnergyCheck,Energy
        write(18,*) "Error! Energy is wrong!"
    endif
    EnergyCheck=Energy
    write(18,*) "check done"
    close(18)
    if(state/=0) then
        call saveconfig()
        call MPI_finalize(ierr); STOP
    endif
  end SUBROUTINE

  SUBROUTINE check_kink_order(Vertex,Site,num,pnum)
    implicit none
    integer :: Site,Vertex,num,i,NVertex,pnum
    double precision  :: Time
    do i=1,4
        NVertex=Ngs(i,Vertex)
        if(KinkTime(NVertex,NearestSite(Vertex,Site,i))>KinkTime(Vertex,Site))then
            write(18,*) number,Vertex,Site
            write(18,*) NearestSite(Vertex,Site,i)
            write(18,*) KinkTime(NVertex,NearestSite(Vertex,Site,i)),KinkTime(Vertex,Site)
            state=error_print("NearestSite is wrong!")
        endif
        if(i==num .and. NearestSite(Vertex,Site,num)/=NeighSite(Vertex,Site)) then
            write(18,*) number,Vertex,Site
            write(18,*) NeighSite(Vertex,Site),NearestSite(Vertex,Site,num)
            write(18,*) KinkTime(Vertex,Site),KinkTime(NVertex,NeighSite(Vertex,Site)),  &
                &  KinkTime(NVertex,NearestSite(Vertex,Site,num))
            state=error_print("Neigh Sites are not at the same time!")
        endif
        if(i==pnum .and. NearestSite(Vertex,Site,pnum)/=PairedSite(Vertex,Site)) then
            write(18,*) number,Vertex,Site
            write(18,*) PairedSite(Vertex,Site),NearestSite(Vertex,Site,pnum)
            write(18,*) KinkTime(Vertex,Site),KinkTime(NVertex,PairedSite(Vertex,Site)),  &
                &  KinkTime(NVertex,NearestSite(Vertex,Site,pnum))
            state=error_print("Paired Sites are not at the same time!")
        endif
    enddo
    return
  end SUBROUTINE

  SUBROUTINE check_kink_weight(Vertex,Site,num)
    implicit none
    integer :: Site,Vertex,num,NVertex,NSite
    double precision  :: weight,Time
    Time=KinkTime(Vertex,Site)
    NVertex=Ngs(num,Vertex)
    NSite=NeighSite(Vertex,Site)
    weight=JJQ-two_plaquette_weight(Vertex,Site,NVertex,NSite,num,Time)
    if(abs(kinkweight(Vertex,Site)-weight)>1e-8) then
        write(18,*) number,Vertex,Site,kinkweight(Vertex,Site),weight
        state=error_print("Kink Weight is wrong!")
    endif
    if(abs(KinkWeight(Vertex,Site)-KinkWeight(NVertex,NSite))>1e-8) then
        write(18,*) number,Vertex,Site,kinkweight(Vertex,Site),KinkWeight(NVertex,NSite)
        state=error_print("Kink Weight dismatch!")
    endif
  end SUBROUTINE

  integer FUNCTION error_print(str)
      implicit none
      character(*)  :: str
      write(18,*) "ERROR!"
      write(18,*) str
      error_print=1
  end FUNCTION

  !==============Collect data ========================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE coll_data(iblck)
    implicit none
    integer, intent(in) :: iblck
    integer             :: j
    do j = 1, NObs_b
      Obs(j,iblck) = Obs(j,iblck)+ Quan(j)
    enddo
  END SUBROUTINE coll_data 
  !===================================================================

  !==============Normalize by Nsamp ==================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE norm_Nsamp(iblck)
    implicit none
    integer, intent(in) :: iblck
    integer             :: j
    double precision    :: nor
    nor = 1.d0/(Nsamp*1.d0)
    do j = 1, NObs_b
      Obs(j,iblck) = Obs(j,iblck)*nor
    enddo
  END SUBROUTINE norm_Nsamp
  !===================================================================

  !==============Statistics ==========================================
  !! THIS IS PROJECT-INDEPENDENT 
  !! In a way learned from Alan
  SUBROUTINE stat_alan 
    implicit none
    integer          :: i, j, k, k0
    double precision :: devn, devp, nor
    double precision, allocatable :: Aux(:)

    nor  = 1.d0/(NBlck*1.d0)

    ! -- calculate average -------------------------------------------
    do j = 1, NObs_b
      Ave(j) = nor*Sum(Obs(j,1:NBlck))
    enddo

    Coarsen: do
      prt = .true.
      ! -- calculate error and t=1 correlation for basics obs.--------
      DO j = 1, NObs_b
        devp = 0.d0;  Cor(j) = 0.d0;  Dev(j) = 0.d0
        do k = 1,  NBlck
          devn   = Obs(j,k)-Ave(j)
          Dev(j) = Dev(j)+devn*devn
          Cor(j) = Cor(j)+devn*devp
          devp   = devn
        enddo 
        Dev(j)   = Dev(j)*nor;        Cor(j) = Cor(j)*nor
        if(Dev(j)>eps)                Cor(j) = Cor(j)/Dev(j)
        Dev(j)   = dsqrt(Dev(j)/(NBlck-1.d0))
        if(dabs(Cor(j))>tol) prt = .false.
      ENDDO 

      IF(prt)                         EXIT Coarsen 
      IF(NBlck<64)    THEN
        prt = .false.;                EXIT Coarsen 
      ENDIF

      ! -- coarsen blocking ------------------------------------------
      NBlck = NBlck/2;    nor=nor*2.d0
      DO j = 1, NObs
        k0 = 1
        do k   = 1, NBlck
          Obs(j,k) = (Obs(j,k0)+Obs(j,k0+1))*0.5d0
          k0 = k0 +2
        enddo 
      ENDDO 
    enddo Coarsen 

    ! -- define auxillary variables and average of composite obs.-----
    call cal_Obs_comp

    ! -- calculate error and t=1 correlation for composite obs.-----
    do j = 1+NObs_b, NObs
      devp = 0.d0;  Cor(j) = 0.d0;  Dev(j) = 0.d0
      DO k = 1,  NBlck
        devn   = Obs(j,k)-Ave(j)  
        Dev(j) = Dev(j)+devn*devn
        Cor(j) = Cor(j)+devn*devp
        devp   = devn
      ENDDO
      Dev(j)   = Dev(j)*nor;        Cor(j) = Cor(j)*nor
      IF(Dev(j)>eps)                Cor(j) = Cor(j)/Dev(j)
      Dev(j)   = dsqrt(Dev(j)/(NBlck-1.d0))
    enddo
    return
  END SUBROUTINE stat_alan
  !===================================================================

  !===============Shift register random number generator =============
  !  very long period sequential version
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE set_RNG
    implicit none
    integer :: i_r,k_r,k1_r
    integer :: iseed

    nrannr = mxrn
    iseed  = iabs(Seed)+1
    k_r    = 3**18+2*iseed
    k1_r   = 1313131*iseed
    k1_r   = k1_r-(k1_r/mod2)*mod2

    do i_r = 1, len1
      k_r  = k_r *mult
      k1_r = k1_r*mul2
      k1_r = k1_r-(k1_r/mod2)*mod2
      ir1(i_r) = k_r+k1_r*8193
    enddo

    do i_r = 1, len2
      k_r  = k_r *mult
      k1_r = k1_r*mul2
      k1_r = k1_r-(k1_r/mod2)*mod2
      ir2(i_r) = k_r+k1_r*4099
    enddo

    do i_r = 1, len1
      inxt1(i_r) = i_r+1
    enddo
    inxt1(len1) = 1
    ipnt1 = 1
    ipnf1 = ifd1+1

    do i_r = 1, len2
      inxt2(i_r) = i_r+1
    enddo
    inxt2(len2) = 1
    ipnt2 = 1
    ipnf2 = ifd2 + 1
    return
  END SUBROUTINE set_RNG 
  !===================================================================

 !===============Calculate next random number =======================
  !! THIS IS ALMOST PROJECT-INDEPENDENT 
  double precision function rn()
  !integer function rn()
    implicit none
    integer   :: i_r, l_r, k_r
    nrannr = nrannr +1
    if(nrannr>=mxrn) then
      nrannr = 1
      do i_r= 1, mxrn
        l_r = ieor(ir1(ipnt1),ir1(ipnf1))
        k_r = ieor(ir2(ipnt2),ir2(ipnf2))
        irn(i_r) = ieor(k_r,l_r)
        ir1(ipnt1)=l_r
        ipnt1 = inxt1(ipnt1)
        ipnf1 = inxt1(ipnf1)
        ir2(ipnt2) = k_r
        ipnt2 = inxt2(ipnt2)
        ipnf2 = inxt2(ipnf2)
      enddo
    endif 
    !rn = irn(nrannr)
    rn = irn(nrannr)*tm32+0.5d0
  end function rn
  !===================================================================

  !==============Trace elapsed time ==================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE t_elapse(t0)
    implicit none
    integer, intent(in) :: t0   ! '0' for the 1st time
    double precision    :: dt
    
    call date_and_time(date, time, zone, tval)
    t_curr = tval(5)*3600.d0+tval(6)*60.d0+tval(7)+tval(8)*0.001d0  ! seconds 
    h_curr = tval(5)

    if(t0==0) then      ! first time
      t_elap = 0.d0;                          dt = 0.d0
    else
      dt = t_curr-t_prev;   if(h_curr<h_prev) dt = dt+24*3600.d0    ! across midnight
      t_elap = t_elap+dt
    endif
    t_prev = t_curr;      h_prev = h_curr 

    select case(t0)
    case(1)            ! initialization time
      write(6,40) dt
      40 format(/'        set up time:',f16.7,2x,'s')
    case(2)            ! equilibration
      write(6,41) dt
      41 format( ' equilibration time:',f16.7,2x,'s.')
      write(6,47) real(Nsamp)/Ntoss*dt/60
      47 format( ' I need more:',f16.7,2x,'minutes.')
      t_simu = dt/(Ntoss*NBlck*1.d0);
      t_mcmc = t_elap
    case(3)            ! for the whole mcmc chain
      t_mcmc = t_elap-t_mcmc
      t_mcmc = t_mcmc/number
      !t_mcmc = t_mcmc/(Nsamp*NBlck*1.d0)
      t_meas = t_mcmc-t_simu
      write(6,42) t_meas
      42 format( '   measurement time:',f16.7,2x,'s./step ')
      write(6,43) t_mcmc
      43 format( '  markov-chain time:',f16.7,2x,'s./step ')
      write(6,44) (t_meas/t_mcmc)
      44 format( ' rela. measure time:',f16.7)
    case(4)
      write(6,45) dt
      45 format(/' stat.-write   time:',f16.7,2x,'s.')
      write(6,46) t_elap/60.d0
      46 format( ' program ends at   :',f16.7,2x,'minutes')
    end select
    return
  END SUBROUTINE t_elapse
  !===================================================================
  !*******************************************************************
  !            End of PROJECT-INDEPENDENT part 
  !*******************************************************************

END PROGRAM main
!=====================================================================

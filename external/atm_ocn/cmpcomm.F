      MODULE CMP_COMM

      implicit none

! MPI variables
      include 'mpif.h'
 
      integer Coupler_id /0/   ! this is Coupler's id, used to address
                               ! Coupler. This is a default value,
                               ! possibly to be redefined later
!
!     Make Coupler's id 0 if it is active (i.e. communnicating with
! the Component.) Otherwise, make it a negative integer; in this case,
! the Component is standalone.
!

      integer ibuffer_size
      parameter (ibuffer_size=10)
      integer Coupler_rank,my_id,COMM_local,  &
      component_master_rank_global,process_rank_global,  &
      component_master_rank_local /-1/ ,process_rank_local,    &
      component_nprocs,FlexLev,ibuffer(ibuffer_size),nprocs_global  

      integer kind_REAL,kind_INTEGER,MPI_kind_REAL,   &
      kind_alt_REAL,MPI_kind_alt_REAL
      parameter (kind_REAL=8,kind_INTEGER=4)
      parameter (kind_alt_REAL=12-kind_REAL)
!       kind_INTEGER must be number of bytes equal to number of bytes
!     implied by MPI_INTEGER MPI constant; all integers sent/received
!     are of this kind. No value other than 4 is anticipated as of now
!       kind_REAL is type of real data to communicate. The corresponding
!     MPI data type variable MPI_kind_REAL is assigned in CMP_INIT.
!       kind_alt_REAL is alternative type of real data to communicate. 
!     The corresponding MPI data type variable MPI_kind_alt_REAL is
!     assigned in CMP_INIT. (It is used in subroutines CMP_alt_SEND
!     and CMP_alt_RECV,)

      integer,parameter:: max_component_nprocs=10000
      integer,parameter:: max_nparmax=20
      integer,parameter:: component_master_rank_local_dummy=0
      integer nparmax /1/
      integer,dimension(1:max_nparmax):: CCR,CCR_local
      logical cmp_domain,comm_process
      data cmp_domain /.false./
      data comm_process /.false./

      save

      END MODULE CMP_COMM
!
!***********************************************************************
!
      SUBROUTINE CMP_INIT(id,flex)
!                         in  in
!
!     This subroutine must be called by every Component right upon
!     calling MPI_INIT. It assigns a value to the Component communicator
!     COMM_local (which is a global variable in module CMP), to be 
!     thereafter used by the Component in place of
!     MPI_COMM_WORLD wherever it is used by the Component's
!     standalone version. Besides, it stores the Component's id,
!     the process's ranks, and the "flexibility level" (flex) requested
!     by the Component in glob. variables. (The latter parameter affects
!     the mode of communications; for its description, see CMP_SEND and
!     CMP_RECV.) Finally, it starts handshaking with Coupler, receiving
!     the unique (global, i.e. in MPI_COMM_WORLD) Coupler process 
!     rank Coupler_rank from Coupler
                                        ! ibuffer may include additional
                                        ! info to be received
!
      USE CMP_COMM

      implicit none

      integer id,flex

      integer ierr,color,key,status(MPI_STATUS_SIZE),tag,dummy
      character*10 s
      logical izd
!

!        Determine if MPI is initialized, if not initialize
      call MPI_INITIALIZED(izd,ierr)
      if (.not.izd) call MPI_INIT(ierr)

!        Determine MPI send/receive types according to prescribed
!        types for arrays to be communicated
      if (kind_REAL.eq.8) then
        MPI_kind_REAL=MPI_REAL8
        MPI_kind_alt_REAL=MPI_REAL4
      else if (kind_REAL.eq.4) then
        MPI_kind_REAL=MPI_REAL4
        MPI_kind_alt_REAL=MPI_REAL8
      else
        write(s,'(i0)') kind_REAL
        call GLOB_ABORT(1,      &
        'CMP_INIT: illegal value of kind_REAL='//s,1)
      end if
      if (kind_INTEGER.ne.4) then
        write(s,'(i0)') kind_INTEGER
        call GLOB_ABORT(1,      &
        'CMP_INIT: illegal value of kind_INTEGER='//s,1)
      end if

!        Store the Component's id
!
      my_id=id

!        Store the Component's "flexibility level"
!
      IF (flex.ne.1) THEN ! in version with comm_process "flexibility
                          ! level" other than 1 forbidden
        call GLOB_ABORT(flex-1,'CMP_INIT: flex.ne.1, ABORTED',1)
      END IF

      FlexLev=flex

!        Assign a value to the Component communicator
!        COMM_local, to be thereafter used by the Component in place of
!        MPI_COMM_WORLD wherever it is used by the Component's
!        standalone version
!
      color=id
      key=1
!           print*,'CMP_INIT: to call MPI_COMM_SPLIT, color=',color
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,COMM_local,ierr)
      call GLOB_ABORT(ierr,'CMP_INIT: error in MPI_COMM_SPLIT',1)

!        Store the process's global and local ranks
!
!           print*,'CMP_INIT: to call MPI_COMM_RANK for global rank'
      call MPI_COMM_RANK(MPI_COMM_WORLD,process_rank_global,ierr)
      call GLOB_ABORT(ierr,     &
      'CMP_INIT: error in MPI_COMM_RANK(MPI_COMM_WORLD...)',1)
!           print*,'CMP_INIT: to call MPI_COMM_RANK for local rank'
      call MPI_COMM_RANK(COMM_local,process_rank_local,ierr)
      call GLOB_ABORT(ierr,     &
      'CMP_INIT: error in MPI_COMM_RANK(COMM_local...)',1)

!        Store component_nprocs - component's number of processes;
!        calculate global number number of processes;
!        determine whether it is standalone mode and if it is, make
!        Coupler's id negative and return
!
      call MPI_COMM_SIZE(COMM_local,component_nprocs,ierr)

      if (component_nprocs .gt. max_component_nprocs) then
        print*,'CMP_INIT: component_nprocs: ',component_nprocs
        call GLOB_ABORT(ierr,     &
        'CMP_INIT: component_nprocs > max_component_nprocs',1)
      end if

      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs_global,ierr)
      if (component_nprocs.eq.nprocs_global) then
        if(process_rank_local.eq.0) print*,'CMP_INIT: standalone mode'
        Coupler_id=-1
        RETURN
      end if

!        Start handshaking with Coupler (all processes):
!        receive the unique (global, i.e. in MPI_COMM_WORLD) Coupler 
!        process rank Coupler_rank from Coupler
!
      tag=Coupler_id+23456
!           print*,'CMP_INIT: to call MPI_RECV'
      call MPI_RECV(ibuffer,ibuffer_size,MPI_INTEGER,MPI_ANY_SOURCE,tag,   &
      MPI_COMM_WORLD,status,ierr)
      call GLOB_ABORT(ierr,'CMP_INIT: error in MPI_RECV',1)
      Coupler_rank=ibuffer(2)
      if (ibuffer(1).ne.Coupler_id) then
        print*,'CMP_INIT: stopped, rcvd ibuffer(1) value is not C id: ',   &
        ibuffer(1)
        CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
      end if
      if (ibuffer(3).ne.ibuffer_size) then
        print*,'CMP_INIT: stopped, rcvd ibuffer(3) value ',ibuffer(3),   &
        ' is not ibuffer_size=',ibuffer_size
        CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
      end if

!        Inform Coupler that this components exists and is active
!
      call MPI_GATHER(id,1,MPI_INTEGER,dummy,1,MPI_INTEGER,     &
      Coupler_rank,MPI_COMM_WORLD,ierr)

!
!     print*,
!    >'CMP_INIT: ranks: process local, global, Coupler; Coupler_id: ',
!    >process_rank_local,process_rank_global,Coupler_rank,Coupler_id

      return
      END
!
!***********************************************************************
!
      SUBROUTINE CMP_MASTER_RANK(n)
!                               in
! This subroutine is to be called if
!   (1) CMP_INTRO is called with actual argument corresponding to
!     -nparmax, rather than to component_master_rank_local; and
!   (2) component_master_rank_local_dummy
!     is not a suitable value for component_master_rank_local.
! (If called at all) this procedure must be called before CMP_INTRO .

      USE CMP_COMM

      implicit none
 
      integer n
!

      component_master_rank_local=n

      return
      END
!
!***********************************************************************
!
      SUBROUTINE CMP_INTRO(nn)
!                          in
! The actual argument corresponding to nn must be:
! -- if processing in the component is serial, the master rank of the
!    component (hence former name of this dummy - master_rank_local);
! -- if processing in the component is parallel, minus the number of the
!    parallel groups of processes
!
! This routine must be called by all Component's processes (which must
! all know the value of the actual argument corresponding to nn)

      USE CMP_COMM

      implicit none
 
      integer nn

      integer ierr,ibuf(3),color,key,tag
!

!     print*,'CMP_INTRO: entered ',nn,process_rank_local
!    >,Coupler_rank

      if (nn.lt.0) then
        if (component_master_rank_local.lt.0) then
          component_master_rank_local=component_master_rank_local_dummy
        end if ! otherwise, component_master_rank_local has already
               ! been assigned
        nparmax=-nn
      else
        component_master_rank_local=nn
      end if

      if (Coupler_id.lt.0) return    !   <- standalone mode

!        If this process is the Component's master process,
!        complete handshaking with Coupler:
!        "register", i.e. send Component master process global rank 
!        to Coupler. Also, send the requested "flexibility level".
!        (Sending Component's id (in ibuf(1)) is for double-check only.)
!
      if (process_rank_local .eq. component_master_rank_local) then
        component_master_rank_global=process_rank_global
        ibuf(1)=my_id  ! redundant, sent for control only
        if (nn.ge.0) then
          ibuf(2)=process_rank_global
        else
          ibuf(2)=-nparmax
        end if
        ibuf(3)=FlexLev
        tag=my_id+54321
            print*,'CMP_INTRO: to call MPI_SEND ',process_rank_local, &
            process_rank_global
        call MPI_SEND(ibuf,3,MPI_INTEGER,Coupler_rank,tag,   &
        MPI_COMM_WORLD,ierr)
        if (ierr.ne.0) then
          print*,'CMP_INTRO: error in MPI_SEND, process ',    &
          process_rank_global
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
        end if
      end if
!           print*,'CMP_INTRO: returning ',process_rank_local,
!    >      process_rank_global,Coupler_rank
      return
      END
!
!***********************************************************************
!
      SUBROUTINE CMP_INTRO_CCR(CCR_loc_d)
!
      USE CMP_COMM

      implicit none

      integer CCR_loc_d(max_nparmax)
 
      integer i,k,tag,ierr
!
      CCR_local(1:nparmax)=CCR_loc_d(1:nparmax)
      component_master_rank_local=CCR_local(1)
         print*,'CMP_INTRO_CCR: CCR_local: ',CCR_local(1:nparmax)
!
      if (Coupler_id.lt.0) return    !   <- standalone mode
!

!  Here it is assumed that CCR_local(k) is the least process rank
!  in the group k of processes:
!x    k=1
!x    do i=0,component_nprocs-1
!x      if (k.lt.nparmax .and. i.ge.CCR_local(k+1)) then
!x        k=k+1
!x      end if
!x      master(i)=CCR_local(k)
!x    end do
!
      do k=1,nparmax
        if (process_rank_local.eq.CCR_local(k)) then
          CCR(k)=process_rank_global
        end if
        call MPI_BCAST(CCR(k),1,MPI_INTEGER, &
        CCR_local(k),COMM_local,ierr)
        call GLOB_ABORT(ierr,'CMP_INTRO_CCR: error in MPI_BCAST',1)
      end do
 
      if (process_rank_local.eq.component_master_rank_local) then
        tag=my_id+54321
        call MPI_SEND(CCR,nparmax,MPI_INTEGER,Coupler_rank,tag,   &
        MPI_COMM_WORLD,ierr)
      end if

      return
      END
!
!***********************************************************************
!
      SUBROUTINE CMP_SEND(F,N)
!
      USE CMP_COMM

      implicit none
 
      integer N,ierr,tag
      real(kind=kind_REAL) F(N)
!
      if (Coupler_id.lt.0) return    !   <- standalone mode

!           call CMP_DBG_CR(6,'CMP_SEND: entered')

      IF (process_rank_local.ne.component_master_rank_local .AND. &
          .not. comm_process) THEN
        if (FlexLev.eq.0) then
!         With "flexibility level" FlexLev=0, only Component master 
!         process is supposed to call this subroutine.
          print '("CMP_SEND: process_rank_local=",i4)',  &
          process_rank_local
          print '("and component_master_rank_local=",i4," differ:")', &
          component_master_rank_local 
          print '("*** STOPPED ***")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
        else if (FlexLev.eq.1) then
!         With "flexibility level" FlexLev=1, any Component process is 
!         allowed to call this subroutine but only the Component
!         master process can actually send data (so the
!         others just make a dummy call), as the Coupler process only 
!         receives data from the Component master process.
          return
        else if (FlexLev.ne.2 .and. FlexLev.ne.3) then
          print '("*** CMP_SEND: illegal value of FlexLev",i9)',FlexLev
          print '("*** STOPPED")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
        end if
!         With "flexibility level" FlexLev=2 or FlexLev=3, any 
!         Component process is allowed to actually send data.
!         [In this case, the Coupler process (in CPL_RECV) receives 
!         from MPI_ANY_SOURCE rather than component_master_rank_global,
!         and it is only identification by  tag  which enables Coupler
!         to receive the data from the right source.]
!         But in any case only one Component process may actually be
!         engaged in a particular exchange of data with Coupler.
      END IF

      tag=my_id

      call MPI_SEND(F,N,MPI_kind_REAL,Coupler_rank,tag,   &
      MPI_COMM_WORLD,ierr)
      call GLOB_ABORT(ierr,'CMP_SEND: error in MPI_SEND',1)

!           call CMP_DBG_CR(6,'CMP_SEND: exiting')
      return
      END
!
!***********************************************************************
!
      SUBROUTINE CMP_alt_SEND(F,N)
!
      USE CMP_COMM

      implicit none
 
      integer N,ierr,tag
      real(kind=kind_alt_REAL) F(N)
!
      if (Coupler_id.lt.0) return    !   <- standalone mode

!           call CMP_DBG_CR(6,'CMP_alt_SEND: entered')

      IF (process_rank_local.ne.component_master_rank_local .AND. &
          .not. comm_process) THEN
        if (FlexLev.eq.0) then
!         With "flexibility level" FlexLev=0, only Component master 
!         process is supposed to call this subroutine.
          print '("CMP_alt_SEND: process_rank_local=",i4)', &
          process_rank_local
          print '("and component_master_rank_local=",i4," differ:")', &
          component_master_rank_local
          print '("*** STOPPED ***")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
        else if (FlexLev.eq.1) then
!         With "flexibility level" FlexLev=1, any Component process is 
!         allowed to call this subroutine but only the Component
!         master process can actually send data (so the
!         others just make a dummy call), as the Coupler process only 
!         receives data from the Component master process.
          return
        else if (FlexLev.ne.2 .and. FlexLev.ne.3) then
          print '("*** CMP_alt_SEND: illegal value of FlexLev",i9)',FlexLev
          print '("*** STOPPED")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
        end if
!         With "flexibility level" FlexLev=2 or FlexLev=3, any 
!         Component process is allowed to actually send data.
!         [In this case, the Coupler process (in CPL_RECV) receives 
!         from MPI_ANY_SOURCE rather than component_master_rank_global,
!         and it is only identification by  tag  which enables Coupler
!         to receive the data from the right source.]
!         But in any case only one Component process may actually be
!         engaged in a particular exchange of data with Coupler.
      END IF

      tag=my_id

      call MPI_SEND(F,N,MPI_kind_alt_REAL,Coupler_rank,tag,   &
      MPI_COMM_WORLD,ierr)
      call GLOB_ABORT(ierr,'CMP_SEND: error in MPI_SEND',1)

!           call CMP_DBG_CR(6,'CMP_SEND: exiting')
      return
      END
!
!***********************************************************************
!
      SUBROUTINE CMP_gnr_SEND(F,N,MPI_DATATYPE)
!
      USE CMP_COMM

      implicit none
 
      integer N,MPI_DATATYPE
      integer F(1)

      integer ierr,tag
!

      if (Coupler_id.lt.0) return    !   <- standalone mode

!           call CMP_DBG_CR(6,'CMP_alt_SEND: entered')

      IF (process_rank_local.ne.component_master_rank_local .AND. &
          .not. comm_process) THEN
        if (FlexLev.eq.0) then
!         With "flexibility level" FlexLev=0, only Component master 
!         process is supposed to call this subroutine.
          print '("CMP_gnr_SEND: process_rank_local=",i4)', &
          process_rank_local
          print '("and component_master_rank_local=",i4," differ:")', &
          component_master_rank_local
          print '("*** STOPPED ***")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
        else if (FlexLev.eq.1) then
!         With "flexibility level" FlexLev=1, any Component process is 
!         allowed to call this subroutine but only the Component
!         master process can actually send data (so the
!         others just make a dummy call), as the Coupler process only 
!         receives data from the Component master process.
          return
        else if (FlexLev.ne.2 .and. FlexLev.ne.3) then
          print '("*** CMP_gnr_SEND: illegal value of FlexLev",i9)',FlexLev
          print '("*** STOPPED")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
        end if
!         With "flexibility level" FlexLev=2 or FlexLev=3, any 
!         Component process is allowed to actually send data.
!         [In this case, the Coupler process (in CPL_RECV) receives 
!         from MPI_ANY_SOURCE rather than component_master_rank_global,
!         and it is only identification by  tag  which enables Coupler
!         to receive the data from the right source.]
!         But in any case only one Component process may actually be
!         engaged in a particular exchange of data with Coupler.
      END IF

      tag=my_id

      call MPI_SEND(F,N,MPI_DATATYPE,Coupler_rank,tag,     &
      MPI_COMM_WORLD,ierr)
      call GLOB_ABORT(ierr,'CMP_SEND: error in MPI_SEND',1)

!           call CMP_DBG_CR(6,'CMP_SEND: exiting')
      return
      END
!
!***********************************************************************
!
      SUBROUTINE CMP_INTEGER_SEND(F,N)
!
      USE CMP_COMM

      implicit none
 
      integer N,ierr,tag
      integer F(N)
!
      if (Coupler_id.lt.0) return    !   <- standalone mode

!           print*,'CMP_INTEGER_SEND: entered with N=',N,' F=',F,
!    >      '; my_id=',my_id,'Coupler_rank=',Coupler_rank

      IF (process_rank_local.ne.component_master_rank_local .AND. &
          .not. comm_process) THEN
        if (FlexLev.eq.0) then
!         With "flexibility level" FlexLev=0, only Component master 
!         process is supposed to call this subroutine.
          print '("CMP_INTEGER_SEND: process_rank_local=",i4)', &
          process_rank_local
          print '("and component_master_rank_local=",i4," differ:")', &
          component_master_rank_local
          print '("*** STOPPED ***")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
        else if (FlexLev.eq.1) then
!         With "flexibility level" FlexLev=1, any Component process is 
!         allowed to call this subroutine but only the Component
!         master process can actually send data (so the
!         others just make a dummy call), as the Coupler process only 
!         receives data from the Component master process.
          return
        else if (FlexLev.ne.2 .and. FlexLev.ne.3) then
          print '("CMP_INTEGER_SEND: illegal value of FlexLev",i9)',FlexLev
          print '("*** STOPPED")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
        end if
!         With "flexibility level" FlexLev=2 or FlexLev=3, any 
!         Component process is allowed to actually send data.
!         [In this case, the Coupler process (in CPL_RECV) receives 
!         from MPI_ANY_SOURCE rather than component_master_rank_global,
!         and it is only identification by  tag  which enables Coupler
!         to receive the data from the right source.]
!         But in any case only one Component process may actually be
!         engaged in a particular exchange of data with Coupler.
      END IF

      tag=my_id
            print*,'CMP_INTEGER_SEND: to call MPI_SEND; F=',     &
            F,' N=',N,' Coupler_rank=',Coupler_rank,' tag=',tag
      call MPI_SEND(F,N,MPI_INTEGER,Coupler_rank,tag,    &
      MPI_COMM_WORLD,ierr)
      call GLOB_ABORT(ierr,'CMP_INTEGER_SEND: error in MPI_SEND',1)
            print*,'CMP_INTEGER_SEND: to return'

      return
      END
!
!***********************************************************************
!
      SUBROUTINE CMP_RECV(F,N)
!
      USE CMP_COMM

      implicit none
 
      integer N,ierr,tag,ibuf(3),status(MPI_STATUS_SIZE)
      real(kind=kind_REAL) F(N)
!
      if (Coupler_id.lt.0) return    !   <- standalone mode

!           call CMP_DBG_CR(6,'CMP_RECV: entered')

      IF (process_rank_local.ne.component_master_rank_local .AND. &
          .not. comm_process) THEN

        if (FlexLev.eq.0) then

!         With "flexibility level" FlexLev=0, only Component master 
!         process is supposed to call this subroutine.

          print '("CMP_RECV: process_rank_local=",i4)', &
          process_rank_local
          print '("and component_master_rank_local=",i4," differ:")', &
          component_master_rank_local
          print '("*** STOPPED ***")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)

        else if (FlexLev.eq.1 .or. FlexLev.eq.2) then

!         With "flexibility level" FlexLev=1 or FlexLev=2, any 
!         Component process is allowed to call this subroutine but 
!         only the Component master process is supposed to actually 
!         receive data (so the others just make a dummy call), as
!         the Coupler process only sends data to the Component master
!         process.

          return

        else if (FlexLev.eq.3) then

!         With "flexibility level" FlexLev=3, any Component process
!         may actually receive data.
!         [In this case, the Coupler process (in CPL_SEND) first
!         receives the Component process global rank 
!         (process_rank_global) from this subroutine, the source being
!         MPI_ANY_SOURCE, so it is only identification by  tag  which 
!         enables Coupler to receive process_rank_global from the right
!         source. Upon the receipt, the Coupler process (in CPL_SEND)
!         sends the data to this Component process, rather than to 
!         the Component master process as is the case with lower 
!         "flexibility levels".]
!         But in any case only one Component process may actually be
!         engaged in a particular exchange of data with Coupler.

          ibuf(1)=my_id
          ibuf(2)=process_rank_global
          tag=my_id
          call MPI_SEND(ibuf,2,MPI_INTEGER,Coupler_rank,tag,   &
          MPI_COMM_WORLD,ierr)
          call GLOB_ABORT(ierr,'CMP_RECV: error in MPI_SEND',1)

        else

          print '("*** CMP_RECV: illegal value of FlexLev",i9)',FlexLev
          print '("*** STOPPED")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)

        end if

      END IF

      tag=my_id
      call MPI_RECV(F,N,MPI_kind_REAL,Coupler_rank,tag,  &
      MPI_COMM_WORLD,status,ierr)
      call GLOB_ABORT(ierr,'CMP_RECV: error in MPI_RECV',1)

!           call CMP_DBG_CR(6,'CMP_RECV: exiting')

      return
      END
!
!***********************************************************************
!
      SUBROUTINE CMP_alt_RECV(F,N)
!
      USE CMP_COMM

      implicit none
 
      integer N,ierr,tag,ibuf(3),status(MPI_STATUS_SIZE)
      real(kind=kind_alt_REAL) F(N)
!
      if (Coupler_id.lt.0) return    !   <- standalone mode

!           call CMP_DBG_CR(6,'CMP_alt_RECV: entered')

      IF (process_rank_local.ne.component_master_rank_local .AND. &
          .not. comm_process) THEN

        if (FlexLev.eq.0) then

!         With "flexibility level" FlexLev=0, only Component master 
!         process is supposed to call this subroutine.

          print '("CMP_alt_RECV: process_rank_local=",i4)',&
          process_rank_local
          print '("and component_master_rank_local=",i4," differ:")', &
          component_master_rank_local
          print '("*** STOPPED ***")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)

        else if (FlexLev.eq.1 .or. FlexLev.eq.2) then

!         With "flexibility level" FlexLev=1 or FlexLev=2, any 
!         Component process is allowed to call this subroutine but 
!         only the Component master process is supposed to actually 
!         receive data (so the others just make a dummy call), as
!         the Coupler process only sends data to the Component master
!         process.

          return

        else if (FlexLev.eq.3) then

!         With "flexibility level" FlexLev=3, any Component process
!         may actually receive data.
!         [In this case, the Coupler process (in CPL_SEND) first
!         receives the Component process global rank 
!         (process_rank_global) from this subroutine, the source being
!         MPI_ANY_SOURCE, so it is only identification by  tag  which 
!         enables Coupler to receive process_rank_global from the right
!         source. Upon the receipt, the Coupler process (in CPL_SEND)
!         sends the data to this Component process, rather than to 
!         the Component master process as is the case with lower 
!         "flexibility levels".]
!         But in any case only one Component process may actually be
!         engaged in a particular exchange of data with Coupler.

          ibuf(1)=my_id
          ibuf(2)=process_rank_global
          tag=my_id
          call MPI_SEND(ibuf,2,MPI_INTEGER,Coupler_rank,tag, &
          MPI_COMM_WORLD,ierr)
          call GLOB_ABORT(ierr,'CMP_alt_RECV: error in MPI_SEND',1)

        else

          print '("CMP_alt_RECV: illegal value of FlexLev",i9)',FlexLev
          print '("*** STOPPED")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)

        end if

      END IF

      tag=my_id
      call MPI_RECV(F,N,MPI_kind_alt_REAL,Coupler_rank,tag,  &
      MPI_COMM_WORLD,status,ierr)
      call GLOB_ABORT(ierr,'CMP_alt_RECV: error in MPI_RECV',1)

!           call CMP_DBG_CR(6,'CMP_alt_RECV: exiting')

      return
      END
!
!***********************************************************************
!
      SUBROUTINE CMP_gnr_RECV(F,N,MPI_DATATYPE)
!
      USE CMP_COMM

      implicit none
 
      integer N,MPI_DATATYPE
      integer F(1)

      integer ierr,tag,ibuf(3),status(MPI_STATUS_SIZE)
!

      if (Coupler_id.lt.0) return    !   <- standalone mode

!           call CMP_DBG_CR(6,'CMP_gnr_RECV: entered')

      IF (process_rank_local.ne.component_master_rank_local .AND. &
          .not. comm_process) THEN

        if (FlexLev.eq.0) then

!         With "flexibility level" FlexLev=0, only Component master 
!         process is supposed to call this subroutine.

          print '("CMP_gnr_RECV: process_rank_local=",i4)', &
          process_rank_local
          print '("and component_master_rank_local=",i4," differ:")',&
          component_master_rank_local
          print '("*** STOPPED ***")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)

        else if (FlexLev.eq.1 .or. FlexLev.eq.2) then

!         With "flexibility level" FlexLev=1 or FlexLev=2, any 
!         Component process is allowed to call this subroutine but 
!         only the Component master process is supposed to actually 
!         receive data (so the others just make a dummy call), as
!         the Coupler process only sends data to the Component master
!         process.

          return

        else if (FlexLev.eq.3) then

!         With "flexibility level" FlexLev=3, any Component process
!         may actually receive data.
!         [In this case, the Coupler process (in CPL_SEND) first
!         receives the Component process global rank 
!         (process_rank_global) from this subroutine, the source being
!         MPI_ANY_SOURCE, so it is only identification by  tag  which 
!         enables Coupler to receive process_rank_global from the right
!         source. Upon the receipt, the Coupler process (in CPL_SEND)
!         sends the data to this Component process, rather than to 
!         the Component master process as is the case with lower 
!         "flexibility levels".]
!         But in any case only one Component process may actually be
!         engaged in a particular exchange of data with Coupler.

          ibuf(1)=my_id
          ibuf(2)=process_rank_global
          tag=my_id
          call MPI_SEND(ibuf,2,MPI_INTEGER,Coupler_rank,tag,  &
          MPI_COMM_WORLD,ierr)
          call GLOB_ABORT(ierr,'CMP_gnr_RECV: error in MPI_SEND',1)

        else

          print '("CMP_gnr_RECV: illegal value of FlexLev",i9)',FlexLev
          print '("*** STOPPED")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)

        end if

      END IF

      tag=my_id
      call MPI_RECV(F,N,MPI_DATATYPE,Coupler_rank,tag,   &
      MPI_COMM_WORLD,status,ierr)
      call GLOB_ABORT(ierr,'CMP_gnr_RECV: error in MPI_RECV',1)

!           call CMP_DBG_CR(6,'CMP_gnr_RECV: exiting')

      return
      END
!
!***********************************************************************
!
      SUBROUTINE CMP_ANNOUNCE(nunit,s)
!
      USE CMP_COMM

      implicit none

      character*(*) s
 
      integer nunit,ierr
!

      IF (process_rank_local.ne.component_master_rank_local .AND. &
          .not. comm_process) THEN
        if (FlexLev.eq.0) then
!         With "flexibility level" FlexLev=0, only an entitled 
!         process is supposed to call this subroutine.
          print '("CMP_ANNOUNCE: *wrong* process_rank_local=",i4)', &
          process_rank_local
          print '("*** STOPPED ***")'
          CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
        end if
        RETURN
      END IF

      write(nunit,*) trim(s)

      return
      END
!
!***********************************************************************
!
      SUBROUTINE CMP_STDOUT(s)
!
!     USE CMP_COMM, ONLY: Coupler_id,process_rank_global
        ! <- These values may not have the right value by this moment,
        ! as this routine may be called before CMP_INIT  - 02/23/05

      implicit none

      character*(*) s
      integer ios
      character*4 mess
!

! -> For debugging:
      OPEN(12345,   &
      file='/nfsuser/g01/wx20ds/C/cmp.stdout',   &
      form='formatted',status='old',iostat=ios)
      if (ios.eq.0) then
        read(12345,*) mess
        if (mess.eq.'mess') then
!         print*,'CMP_STDOUT: unit 6 left alone, process ',
!    >    process_rank_global
        ! <- process_rank_global may be undefined by this moment, as
        !    this routine may be called before CMP_INIT  - 02/23/05
          RETURN
        end if
        CLOSE(12345)
      end if
! <- for debugging

!     if (Coupler_id.lt.0) RETURN    ! Nothing is to occur if there is
                                     ! no communication with Coupler,
                                     ! i.e. if Component is standalone
        ! <- Coupler_id may not have the right value by this moment,
        ! as this routine may be called before CMP_INIT  - 02/23/05

      if (len_trim(s).eq.0) RETURN

      close(6)
      
      open(6,file=trim(s),form='formatted',status='unknown')

      print*,'CMP_STDOUT: unit 6 closed, reopened as '//trim(s)

      return
      END
!
!***********************************************************************
!
      SUBROUTINE CMP_DBG_CR(nunit,s)
!
!       Debugging routine: mainly, prints Coupler_rank
!
      USE CMP_COMM

      implicit none

      character*(*) s
      integer nunit

      integer ncall/0/,ncallmax/5000/
      save
!

      if (s(5:6).eq.'m:') then
        if (process_rank_local .ne. component_master_rank_local) RETURN
      end if

      if (ncall.ge.ncallmax) RETURN
      ncall=ncall+1

      write(nunit,*)process_rank_global,ncall,Coupler_id,Coupler_rank,s

! The following assumes that Coupler_rank must be =0, comment out if
! this is not the case
      call GLOB_ABORT(Coupler_rank,    &
      'CMP_DBG_CR: Coupler_rank.ne.0, aborting',1)

      return
      END
!
!***********************************************************************
!
      SUBROUTINE CMP_FLUSH(nunit)

      USE CMP_COMM

      implicit none

      integer nunit

      integer i,ierr,rc
!

      do i=0,component_nprocs-1
        call MPI_BARRIER(COMM_local,ierr)
        call GLOB_ABORT(ierr,'CMP_FLUSH: MPI_BARRIER failed, aborting',   &
        rc)
        if (i.eq.process_rank_local) call FLUSH(nunit)
      end do

      return
      END
!
!***********************************************************************
!
      subroutine CMP_FINALIZE(izd,ierr)

      USE CMP_COMM

      implicit none

      logical izd
      integer ierr

      integer ierr1,ierr2
!

      ierr=0
      ierr2=0
      call MPI_INITIALIZED(izd,ierr1)
      if (izd) call MPI_FINALIZE(ierr2)
      ierr=abs(ierr1)+abs(ierr2)

      return
      END

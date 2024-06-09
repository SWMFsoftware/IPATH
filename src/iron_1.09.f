c=======================================================================
c      transport_back.f
c
c      Author: Junxiang Hu
c      Date:   02/2016
c       
c      A backward transport code using Monte Carlo method
c      follows focused transport equation from Ruffolo 1995
c
c
c      V1.07:
c      - Introduced ESP modeling
c      - Simulation time is extended to after shock arrival time
c
c      V1.07.5:
c      - Removed p_big checking
c
c=======================================================================
       program transport_back
       implicit none
       include       'mpif.h'
       include       'common.h'

       
             
       integer       i,     j,     k,     ishl,  ip,    iphi,  itime,
     1               ii, jp
       real*8        p_ini, p0,    v0,    e0,    t0,    dt0,   miu0,
     1               e_p,   gamma_p
       real*8        tnew,  rnew,  pnew,  miunew, vnew, phinew,
     1               t,     r,     p,     miu,   v,     theta, phi,
     2               dt,    ds,    dr,    dp,    dmiu,  sec_psi,
     3               dmiu2, miumo, total, psi,   r_min, dt_perp
     
       real*8        LZ,    p11,   e11

       
       real*8        eps1,         eps2,         gyrof,        gamma0,
     1               waveK,        Dmm,          dDmmdm,       Dmiumiu,
     2               eps3  

       integer       phi_e_index,  arrival_index,       t_num
       real*8        arrival_time,        leave_time,   arrival_time_2,
     1               sys_time
       

       !record control
       integer       jrec (seed_num)
       real*8        prec (seed_num),     mrec (seed_num),
     1               trec (seed_num),     phirec(seed_num),
     2               rrec (seed_num)
       real*8        rtrace_1(seed_num),         rtrace_10(seed_num),
     1               phitrace_1(seed_num),       phitrace_10(seed_num)  
       real*8        fpAU (p_num),        feAU (p_num)

       ! system time variables
       character(8)  :: date
       character(10) :: time
       character(5)  :: zone
       integer,dimension(8) :: values



       ! declare for function
       real*8        fp0,   randpc,       parker_field,        inv_erf,
     1               dierfc, deltabs, Aq, meanfreep,    fe0,   getft
      
       ! declare for parallel calculation
       integer       strlen,       myrank,       nodes
       character*2   rankstr
       character*512 flnam

       integer       iwarmUp,      iseed,        ierr,      if_parallel

       character*2   filecha1, filecha2
       integer       il               ! string length      
       

       ! for reading distr files
       integer       line,         Reason
       real*8        temp_a(5)
       real*8, allocatable ::      dist_param(:,:),     momenta(:,:),
     1                             fort(:,:),           fort2(:,:),
     1                             fort3(:,:)
       real*8, allocatable ::      dist_fp(:,:,:),      shell_time(:),
     1                             dist_fp_hi(:,:,:),   shell_loc(:,:),
     2                             dist_fp_up(:,:,:),   p_big(:,:),
     3                             dist_fp_dn(:,:,:),   p_big_in(:,:)  
     

       integer       prev_lines,   leave_index,  arrival_index_2

       real*8, allocatable ::      bndy_loc_in(:,:),  
     1                             bndy_loc_out(:,:)

       real*8        rsh,  p_big_r,  time_norm,  phi_intsec,  fp_intsec,
     1               rsh_in 
       integer       ip_shell, im_shell, ip_phi, im_phi, ip_phi_in, 
     1               im_phi_in 

       !Kappa_perp related
       real*8        kappa_perp0(p_num),  kappa_perp,   check(p_num),
     1               xtemp,        ytemp,        dX  , ratio_temp,
     2               dX1,          dX3,          B_r, kappa_perp1(p_num)
     3               , ratio_temp1


c---------------   Common Block declaration and value  -----------------

       include       'assign_common.h'

c-----------------------------------------------------------------------
       print*, 'Aq_AU', Aq_AU
       print*, 'C(nu)', 1./2./dsqrt(pi)*gamma(5./6.)/gamma(1./3.)
c       Aq_AU =  1./2./dsqrt(pi)*gamma(5./6.)/gamma(1./3.) /2. /pi


c------   Read in the distribution function within shock complex  ------
       
       open(21, file="example_data/dist_at_shock.dat", 
     1        form="formatted")
       line = 0
       do
           read(21,*,IOSTAT=Reason) temp_a(5)
           IF (Reason < 0 ) EXIT
           line= line+1   
       enddo
       
       allocate ( dist_param(5, line))
       rewind 21

       read(21,*) dist_param
       close(21)

       shell_no = line/p_number/phi_no

       print*, "line, shell_no:", line, shell_no

       allocate ( dist_fp(shell_no, p_number, phi_no))
       allocate ( dist_fp_up(shell_no, p_number, phi_no))
       allocate ( dist_fp_dn(shell_no, p_number, phi_no))

       allocate ( dist_fp_hi(shell_no, p_number, phi_no))
       allocate ( shell_time(shell_no),  shell_loc(shell_no, phi_no))

       do i = 1, shell_no
          shell_time(i) = dist_param(4, phi_no*p_number*(i-1)+1)
          do j = 1, phi_no
             shell_loc(i,j)  = dist_param(5, phi_no*p_number*(i-1)+ 
     1                             p_number*(j-1) + 1)
          enddo   
       enddo

       do i = 1, p_number  ! see the acceleration part for this
          p0 = min_p + (i - 1) * del_p
          p_0(i) = exp(p0)* mp * vo
          energy_0(i) = (dsqrt((mp*co**2.0)**2.0 + (mp*exp(p0)*vo*co)
     1          **2.0) - mp*co**2.0 )/ eo
       enddo
       
       do i = 1, phi_no  ! see the acceleration part for this
          phi_0(i) = int(cme_center + del_phi*(i-1-floor(phi_no/2.)))
          if (phi_0(i) .eq. phi_e ) then
              phi_e_index = i
          endif
       enddo

       ishl = 1;     iphi = 1;       ip = 1

       do i = 1, line         
c          dist_fp(ishl, ip, iphi) = dist_param(2, i)
          dist_fp(ishl, ip, iphi) = dist_param(3, i)
          ip = ip + 1
          if (ip .gt. p_number) then
              ip = 1
              iphi = iphi + 1
              if (iphi .gt. phi_no) then
                 iphi = 1
                 ishl = ishl + 1   
              endif
          endif
       enddo

       open(22, file="example_data/shock_momenta.dat",
     1        form="formatted")
       
       allocate ( momenta(6, phi_no*shell_no))
       allocate ( p_big(shell_no, phi_no))
       read(22,*) momenta
       close(22)
       
       do i = 1, shell_no
          do j = 1, phi_no
              p_big(i,j) = momenta(4, (i-1)*phi_no + j)*1.6d-19/co
                                                      ! due to a strange norm
              ! now p_big is in SI unit
          enddo
       ! decide the earth encounter time
          if (shell_loc(i, phi_e_index) .le. r0/AU ) then
              arrival_index = i
          endif
       enddo

!       print*, (shell_time(10)-time_start)*t_o/3600, 
!     1        momenta(2, 9*phi_no+2)

       arrival_time = ((r0/AU- shell_loc(arrival_index, phi_e_index)) *
     1        shell_time(arrival_index+1) + (shell_loc(arrival_index+1,
     2        phi_e_index) -r0/AU)*shell_time(arrival_index) )
     3        /( shell_loc(arrival_index+1, phi_e_index) - 
     4           shell_loc(arrival_index, phi_e_index)  )


       print*, shell_no, phi_no, p_number

 
       open(24, file="example_data/esc_distr_dn.dat",
     1        form="formatted")
       
       allocate ( fort2(4, phi_no*shell_no*p_number))

       read(24,*) fort2
       close(24)
       ishl = 1;     iphi = 1;       ip = 1
       do i = 1, phi_no*shell_no*p_number         
          dist_fp_dn(ishl, ip, iphi) = fort2(2, i)
          ip = ip + 1
          if (ip .gt. p_number) then
              ip = 1
              iphi = iphi + 1
              if (iphi .gt. phi_no) then
                 iphi = 1
                 ishl = ishl + 1   
              endif
          endif
       enddo
       !calculate downstream p_big: p_big_in

       allocate ( p_big_in(shell_no, phi_no))

       do i = 1, shell_no
          do jp = 1, phi_no
              do j = 1, p_number-1
                 if ( (dist_fp_dn(i, j, jp) .ne. 0.0) .AND.    
     1               (dist_fp_dn(i, j+1, jp) .eq. 0.0 ) ) then
                     p_big_in(i,jp) = p_0(j)
                 endif
              enddo
          enddo
       enddo    



       open(25, file="example_data/all_shell_bndy.dat",
     1        form="formatted")
       
       allocate ( fort3(4, (shell_no*(shell_no+1)/2+shell_no)*phi_no))

       read(25,*) fort3
       close(25)


       allocate ( bndy_loc_out(shell_no, phi_no),
     1            bndy_loc_in(shell_no, phi_no))


       print*, phi_e_index 
       do i = 1, shell_no
          prev_lines = (i*(i-1)/2+i-1)*phi_no
          do j = 1, phi_no
             bndy_loc_out(i, j) = fort3(4, prev_lines + (j -1)*(i +1)
     1                      + i + 1)
              ! outter boundary location at each zeus time
             BNdy_loc_in(i, j) = fort3(4, prev_lines + (j - 1)*(i +1)
     1                      + 1 )
              ! inner boundary location at each zeus time
          enddo
       enddo

       do i = 1, shell_no

c          if (bndy_loc_out(i, phi_e_index) .le. 1.0 ) then
c              arrival_index_2 = i
c          endif
          if (bndy_loc_in(i, phi_e_index) .le. r0/AU ) then
              leave_index = i
          endif

       enddo
c       arrival_time_2 =((1.- bndy_loc_out(arrival_index_2,phi_e_index))*
c     1  shell_time(arrival_index_2+1) + (bndy_loc_out(arrival_index_2+1,
c     2       phi_e_index) -1.)*shell_time(arrival_index_2) )
c     3        /( bndy_loc_out(arrival_index_2+1, phi_e_index) - 
c     4           bndy_loc_out(arrival_index_2, phi_e_index)  )

       leave_time = (  (r0/AU - bndy_loc_in(leave_index, phi_e_index)) * 
     1        shell_time(leave_index + 1) + (bndy_loc_in(leave_index+1,
     2        phi_e_index) - r0/AU) * shell_time(leave_index)   )
     3        /( bndy_loc_in(leave_index+1, phi_e_index) 
     4         - bndy_loc_in(leave_index, phi_e_index) )

c       print*, "arrival_index, arrival_index_2",
c     1               arrival_index, arrival_index_2

       print*, "shock complex arrives at the observer at", arrival_time, 
     1        "and leaves at", leave_time, "(in system time)"

c       print*, arrival_time_2
c       print*, 'leave_index', leave_index
c       print*, shell_loc(arrival_index, phi_e_index),
c     1         bndy_loc_out(arrival_index, phi_e_index)

c       print*, shell_loc(arrival_index_2, phi_e_index),
c     1         bndy_loc_out(arrival_index_2, phi_e_index)  



c       pause

c       print*, shell_time(leave_index), shell_time(leave_index+1),
c     1        leave_time
c       print*, bndy_loc_in(leave_index),bndy_loc_in(leave_index+1)
c       pause

!       print*, arrival_index
!       print*, "arrival time", arrival_time, phi_e_index,
!     1        (arrival_time-time_start)*t_o/3600.
!       pause
       


!        do i = 1, shell_no
!        	print*,"shell_no, time, location(au)",
!      1	   	i,shell_time(i), shell_loc(i, 10)
!        enddo
!        pause

c-----------------------------------------------------------------------       

       
       if (if_parallel .eq. 1) then
         call initnodes(myrank,nodes)
         call getrankstr(rankstr,myrank)

         print*, "myrank and nodes", myrank, nodes
       endif

!
!        !calculate mean free path at 1AU at 1MeV

!      calculate Kappa_perp at 1AU numerically
       do j = 1, p_num
          e0 = min_e * ((max_e/min_e)**(1./(p_num-1)))**(j-1.)
          gamma0 = (e0*eo + mp*co**2.0)/(mp*co**2.0)
          p0 = dsqrt(((e0*eo+ mp*co**2.0)**2.-(mp*co**2.0)**2.0)/co**2.)
c          print*, gamma0-dsqrt(1+ (p0/mp/co)**2.)    !should be 0
          v0 = p0/mp/gamma0
          total = 0.0
         do i = 1 , 18999
           miu = real(i/19000., 8)
           miumo = real((i-1)/19000., 8)
           if (i .eq. 1) then
               total = total+ (1.-miu**2.)**2./Dmiumiu(miu,AU,v0,dDmmdm)
     1              *1./19000./2.0
           else
           total = total+0.5*((1.-miu**2.)**2./Dmiumiu(miu,AU,v0,dDmmdm)
     1       +(1.-miumo**2.)**2./Dmiumiu(miumo,AU,v0,dDmmdm))*1./19000. 
           endif
         enddo 
         total = total * 0.25* v0* v0
         ! now this total is K_parallel
         
c        print*, 'Mean free path at 1AU at 100MeV is ', total,'or', total
c     1        /AU/AU, 'AU^2/s'
         kappa_perp0(j) = total *2.* Cturb_AU*s2d /2.
         kappa_perp1(j) = (1./sqrt(3.)*pi*gamma(5./6.)*v0/2./dsqrt(pi)
     1   /(gamma(1./3.))*Cturb_AU*s2d*l2d)**(2./3.)*total**(1./3.)

         ratio_temp = total*kappa_perp0(j)/v0/v0/(3*l2d**2)
         ratio_temp1 = total*kappa_perp1(j)/v0/v0/(3*l2d**2)

         PRINT*, 'energy(eV):', e0, 
     1        'lambda_parallel(au):', total*3./v0/AU, 
     1        'lambda_perp(au):', Kappa_perp1(j)*3./v0/AU
     1        ,'Kappa_perp(m^2/s)', Kappa_perp1(j)
     1        , mp*v0/eo/parker_field(AU)/0.03/AU


c         PRINT*, 'satisfy first condition?', ratio_temp1
     
c         print*, "lambda_para/perp ratio:",Kappa_perp1(j)/total
         
c         if (ratio_temp .lt. 1) then
c            check(j) = 1
c         else
c            check(j) = 0
c            kappa_perp0(j) = (sqrt(3.)*pi*gamma(5./6.)*v0/2./dsqrt(pi)
c     1       /(gamma(1./3.))*Cturb_AU*0.01*AU )**(2./3.)*total**(1./3.)
c         endif

c       print*, kappa_perp0(j),ratio_temp, kappa_perp1(j), ratio_temp1 

       enddo

c       print*, kappa_perp0(j)/v0*3. /AU, kappa_perp1(j)/v0*3. /AU

c       print*, 'kappa_perp:', kappa_perp0
c       print*, check


c       print*, "kappa_perp is", Kappa_perp0
c       print*, kappa_perp0*total/v0/v0, 3./kl_au/kl_au, 3*(0.01*AU)**2

c        pause 
!
!     11 -> distribution of f in p, 
!     12 -> conditioned dist. f
!     13 -> state variable upon crossing inner boundary {energy, time, momentum, pitch} 
!     14 -> distribution of f in T (T is kinetic energy) 
c-----------------------------------------------------------------------
       
       if (if_parallel .eq. 1) then

          open(11, file ="fp_"//rankstr(1:2), 
     1      form='formatted',  status = 'unknown')
          open(13, file ="RawData_"//rankstr(1:2), 
     1      form='formatted', status = 'unknown')

          if (myRank .eq. 0) then
              open(12, file = 'Time_steps', form = 'formatted',
     1             status = 'unknown')
          endif
c          open(14, file ="fe_"//rankstr(1:2), 
c     1      form='formatted', status = 'unknown')

!      for tracing particles
!          open(15, file ="Trace_1hr_"//rankstr(1:2), 
!     1      form='formatted', status = 'unknown')
!          open(16, file ="Trace_10hr_"//rankstr(1:2), 
!     1      form='formatted', status = 'unknown')

          iseed = 100 + (myRank +1)**2 ! random seed.

       else
          open(11, file ="fp_single" ,
     1      form='formatted',  status = 'unknown')
          open(13, file ="RawData_single", 
     1      form='formatted', status = 'unknown')
c          open(14, file ="fe_single", 
c     1      form='formatted', status = 'unknown')
!          open(15, file ="Trace_1hr_single", 
!     1      form='formatted', status = 'unknown')
!          open(16, file ="Trace_10hr_single", 
!     1      form='formatted', status = 'unknown')
          
          ! choose seed based on the current system time
          call date_and_time(date, time, zone, values)
          iseed = 100 + values(6)*values(7) ! (mm x ss)
       endif
c---- generate iseed based on rank 

       do iwarmUp = 1, 100
          call srandpc(iseed) ! warm up the seed.         
       enddo
       
       if ((if_parallel .eq. 1) .and. (myRank .eq. 0)) then
            write(12,*) "shock complex arrives at the observer at", 
     1        arrival_time,"and leaves at", leave_time,
     2         "(in system time)"
       endif     

       t_num = 40
       do itime = 1, t_num ! in V1.05 we only consider up to 1AU
          !time for the interested output f(p)
          time_fp = time_start + itime *(arrival_time -time_start)/t_num

          if ((if_parallel .eq. 1) .and. (myRank .eq. 0)) then
              write(12,*) itime, time_fp
          endif

       do i = 1 , p_num
          
c         p_ini = min_p + (i - 1) * del_p
c         p0    = Exp(p_ini) * mp * vo     ! momentum in SI
c         e0    = (dsqrt((mp*co**2.0)**2.0 + (p0*co)**2.0) - mp*co**2.0)
c     1           / eo   ! Kinetic energy in EV

         !initial e ranges from 1KeV to 1GeV, uniformly on log scale

          e0 = min_e * ((max_e/min_e)**(1./(p_num-1)))**(i-1.)
          gamma0 = (e0*eo + mp*co**2.0)/(mp*co**2.0)
          p0 = dsqrt(((e0*eo + mp*co**2.0)**2. -(mp*co**2.0)**2.0))/co
          v0    = p0 / mp /gamma0 ! speed in SI
          
c     print*, 'energy',e0, v0, gamma0
          

          t0    = 0.0
          fpAU(i) = 0.0

! decide dt0 based on v0
! We don't do this because dt may end up being too small
! dt has to be larger than the gyro period
! dt0   = 4d-4 * (AU /v0 )         
         

          do j = 1, seed_num  

c             miu0 = 8.0D-1       
 
!             try random miu0
             eps3 = randpc()
             miu0 = (eps3-0.5)*2.0D0
c             print*,"miu0 is: ", miu0
             k = 1              ! k is number of steps. [since dt can vary, this is not fully correct]
             
!     initialization
           tnew    = t0
           rnew    = r0
           pnew    = p0
           miunew  = miu0
           phinew  = phi_e  ! phi is the particle location in phi direction
                            ! we use degree here

           r_min   = r0     ! rmin is to keep track of minimum r to sun

           print*, "itime, p no.,seed no.", itime, i, j             

c--------------  start monte carlo for each seed particle

10         continue         ! use if ... continue to be safe for fortran 77 compilers

c           if ( (k .lt. iloop) .and. (rnew .ge. r_bnd) .and.
c     1         (dabs(tnew) .le. (time_fp - time_start)*t_o)) then

           if ( (k .lt. iloop) .and. ( dabs(tnew) .le. (time_fp -
     1          shell_time(1))*t_o) ) then
                ! check if its within the step limit 

c----------- fill the gap where miu is about 0  ------------------------
              if ( dabs(miunew) .le. 1d-2) then
                 eps3   = randpc()
                 miunew = (eps3 - 0.5)/dabs(eps3 - 0.5)*1d-2
                     ! kick lazy particles out of the 'valley of death'
              endif

              t      = tnew
              r      = rnew
              p      = pnew
              miu    = miunew
              phi    = phinew
              
              if (r .lt. r_min) then
                 r_min = r
              endif


              e_p = (dsqrt((mp*co**2.0)**2.0 + (p*co)**2.0)-mp*co**2.0)
     1           / eo   ! Kinetic energy in EV
              gamma_p = (e_p*eo + mp*co**2.0)/(mp*co**2.0)
              
              v      = p/mp/gamma_p
c              print*, 'v0',v0
              theta  = acos(miu)

             
              dt     = - 30.*(r/AU)**1.5

30            continue      ! continue tag for reducing dt



         ! Updating (r, miu, p) using RK4 method
!              print*, r, miu, p
              call RK4(r, miu, p, dt, rnew, miunew, pnew, vnew)            
              
c              print*, rnew, miunew, pnew
c              pause
              if (miunew .gt. 1) then
                 print*, 'Warning! miu is larger than 1!'
              endif
              

c----------- second part of dmiu due to turbulence scattering ----------
c              eps1   = rand()      ! eps is only [0,1), how to make it
c              eps2   = rand()      ! [0,1] or (0,1) ?

c--- use the following way to generate the random numbers.
              eps1   = randpc()   
              eps2   = randpc()   

c----------- fill the gap where miu is about 0  ------------------------
              if ( dabs(miunew) .le. 1d-2) then
                 eps3   = randpc()
                 miunew = (eps3 - 0.5)/dabs(eps3 - 0.5)*1d-2
                     ! kick particles out of the 'valley of death'
              endif

              Dmm = Dmiumiu(miunew, rnew, vnew, dDmmdm)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~              
c              Dmm = Dmiumiu(miunew, rnew, vnew, dDmmdm)
       !      dDmmdm = (Dmiumiu(miu+0.00002*miu, r, v) -Dmm)/0.00002/miu

c----------- fill the gap where miu is about 0  ------------------------
c              if ( dabs(miu) .le. 0.0002) then
c                 Dmm = Dmiumiu(2d-4, r, v, dDmmdm)
c                 dDmmdm = 0.0
c              endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

              dmiu2  = (eps1 - 0.5)/abs(eps1-0.5) *inv_erf(eps2)
     1           *dsqrt(4.*Dmm*dabs(dt)) + dDmmdm * dt
                                        ! changed to + in v1.07              


c----------- When dmiu is too big, we reduce dt by half ----------------
              if ( dabs(miunew - miu + dmiu2) .ge. 0.1) then
                 dt = dt/2.
                 goto 30
              endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c              print*, '1. phinew, rnew', phinew, rnew/AU
c              Kappa_Perp = Kappa_Perp0(i) *(rnew/AU)**1.2
c              print*,'simple    Kperp:',Kappa_perp,     rnew/AU

              B_r = parker_field(rnew)    ! introduce this parameter to optimize speed
              sec_psi = dsqrt(1.0+(omega*rnew/usw)**2.0)
              psi    = acos(1./sec_psi)   ! this psi is the angle 
                                          ! between r and z


              Kappa_perp = (vnew/v0)**(10./9.)*((parker_field(AU)/
     1             B_r))**(7./9.)*(AU/rnew)**(7./6.)* Kappa_Perp1(i)
              ! changed in V1.05

c              print*,'alternate Kperp:', Kappa_perp
c              PRINT*, 'Vnew, v0', vnew/1000, v0/1000
c              if (k .eq. 300) then
c                     print*, Kappa_perp0(i)
c                     pause
c              endif
              eps1   = randpc()   
              eps2   = randpc()   
              dX1     = inv_erf(eps2) *dsqrt(4.*kappa_Perp*dabs(dt))

              ! added in #1.04
c              dX2    = Kappa_perp*B_r* dlog(B_r) * br_bnd* r_bnd**2.
c     1          *(-2.*sec_psi/(rnew**3.)+(omega/usw)**2./(rnew*sec_psi))
c     2          *sin(psi) * dt  ! dt is negative so dX2 is negative
              
              dX3    = (Kappa_perp -(v/v0)**(10./9.)*(parker_field(AU)/
     1        parker_field(r))**(7./9.)*(AU/r)**(7./6.)* Kappa_perp1(i))
     2        /(rnew-r)/sin(psi) * dt

              dX = (eps1-0.5)/abs(eps1-0.5) *( dX1 + dX3 )
       
c              print*, 'dX, dX1, dX3', dX, dX1, dX3
c              if (k .eq. 100) then
c                pause
c              endif

              phinew = phi + omega/usw*(r - rnew)*180./pi
     1               + omega*dt*180./pi
                       ! added in v1.05
c              print*, omega/usw*(r - rnew)*180./pi, omega*dt*180./pi
c              if (k .eq. 40) then
c              pause
c              endif


              xtemp  = rnew*cos(phinew*pi/180.) - 
     1         dX * sin(phinew*pi/180. - psi)  
              ytemp  = rnew*sin(phinew*pi/180.) + 
     1         dX * cos(phinew*pi/180. - psi)
              

              if (xtemp .gt. 0.) then
                 if (ytemp .ge. 0.) then                 
                     phinew = atan(ytemp/xtemp)*180./pi
                 else
                     phinew = atan(ytemp/xtemp)*180./pi + 360.
                 endif
              else if (xtemp .lt. 0) then
                 phinew = atan(ytemp/xtemp)*180./pi + 180.
              else   ! xtemp = 0
                 if (ytemp .gt. 0.) then
                     phinew = 90.0
                 else
                     phinew = 270.0
                 endif
              endif
             
              rnew   = dsqrt(xtemp**2. + ytemp**2.)

c              print*, phinew
c              print*, phi + omega/usw*(r - rnew)*180./pi
c              pause

c              print*, '3. phinew, rnew', phinew, rnew/AU


              phi_at_inner =  phinew + (rnew-r_bnd)/usw*Omega *180./pi

c              gradKappaPerp = ...
c              ds2 = -gradKappaPerp *dt  


!               print*,'miu, dmiu, dmiu2', miu, dmiu,  dmiu2
!               print*,'dmm, dDmmdm', Dmm, dDmmdm

              tnew   = t + dt  !? do we need to add dt_perp ?
              miunew = miunew + dmiu2

!               print*,'tnew, rnew, pnew, miunew', tnew, rnew, pnew
!      1  , miunew
!               pause

c--------------    when miunew > 1 or < -1 -----------------------------

              if (miunew .gt. 1.) then
                !print*, "warning!!@!"
                miunew = 1. - (miunew - 1.)
              endif

              if (miunew .lt. -1.) then
                miunew = -1. + (-1.-miunew)
              endif 

c              if (miunew* miu .lt. 0) then
c                     print*, "cross happens", miunew, miu
c              endif



c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


              
c              print*,'r, miu = ', rnew, miunew
              
c----- jump out of iteration when particle reaches inner boundary              
c              if (rnew .le. r_inner) then

        
c----- jump out of iteration when particle reaches the shock front              
              ! Find the intersect point between the field line and 
              ! the shock front. Decide the fp at the intersection 

c              phinew = phi_e + omega/usw* (AU - rnew)*180./pi
              ! this expression is different when we introduce Kappa_perp

              if ( phinew .gt. phi_0(phi_no) 
     1        .or. phinew .lt. phi_0(1)) then
                 ! particle is out of phi range, keep going
                 k = k + 1
                 goto 10
              endif


              ! check inner and outter boundary
              if (rnew .lt. r_bnd) then
                 print*, "Particle gets too close."
                 goto 14
              endif

              if (rnew .gt. 5*AU) then
                 print*, "Particle is too far outside.",
     1            "minimum distance is:", r_min/AU,"(AU)"
                 ! this is not necessary, but can speed up the code
                 goto 14
              endif
                   
              time_norm = time_fp + tnew /t_o
       
              ip_shell = 0;		im_shell = 0
              ip_phi	= 0;		im_phi	= 0

              do ii = 1, shell_no
                 if (shell_time(ii) .ge. time_norm) then
                     ip_shell = ii
                     im_shell = ii-1
                     exit
                 endif
              enddo

              if (ip_shell .eq. 1) then
                 print*, "Transport for too long"
                 ! current time is smaller than the first shell time
                 ! this particle is not accelerated by shock
                 goto 14
              endif
  
       
              do ii = 1, phi_no
                 if (phi_0(ii) .ge. phinew) then
                     ip_phi = ii
                     im_phi = ii-1
                     exit
                 endif
              enddo


             ! added this in v1.07 to speed up calculations
             ! maybe not necessary
c              if ((rnew/AU .gt. shell_loc(ip_shell, floor(phi_no/2.)+1))
c     1   .and.(rnew/AU .gt. shell_loc(ip_shell, floor(phi_no/2.)))) then
c                  ! particle is still far upstream, keep going
c                 k = k + 1
c                 goto 10
c              endif
       
c-----------------------------------------------------------------------

              call get_rsh(time_norm, phinew, shell_time, shell_loc, 
     1               bndy_loc_in, rsh, rsh_in, ip_shell, im_shell, 
     2               ip_phi, im_phi)

c              print*, rsh/AU, rsh_in/AU
c              pause
   
c-----------------------------------------------------------------------
c      added in V1.07 to check if the particle is from downstream
c      of the shock when the shock complex passes the Earth
c-----------------------------------------------------------------------
              goto 666
              if ((rnew .ge. rsh_in) .and. (r .lt. rsh_in)) then
                 !encounter shock downstream

              ! record particle as originated from downstream
c                   print*, "rsh_in, rsh", rsh_in/AU, rsh/AU
                   fp_intsec = getft(pnew, tnew, rnew, shell_time, 
     1                      bndy_loc_in, phinew, dist_fp_dn, ip_shell, 
     2                      im_shell, ip_phi, im_phi)
                   fp_intsec = fp_intsec * (rnew/r0)**2.0

                   print*, "origin at downstream of the shock complex"

                   goto 15
              endif

666           continue


c-----------------------------------------------------------------------
        ! criterion changed in #1.04 to make sure particle is from upstream
              if ((rnew .le. rsh) .and. (r .gt. rsh)) then

c                print*, 'rsh/AU, rnew/AU', rsh/AU, rnew/AU
c                print*, 'phinew', phinew
                ! record the particle

c                print*, "rsh_in, rsh", rsh_in/AU, rsh/AU
                fp_intsec = getft(pnew, tnew, rnew, shell_time, 
     1               shell_loc, phinew, dist_fp, ip_shell, im_shell, 
     2               ip_phi, im_phi)

                fp_intsec = fp_intsec * (rnew/r0)**2.0
                
                fp_intsec = fp_intsec*exp(-2.0)
                !  introduce the effect of escape length here

                goto 15

              endif
               
              if (isnan(rnew)) then
              ! don't record particle, skip to next one
              ! this is a very rare case

                     jrec(j)  = 0
                     prec(j)  = 0.
                     mrec(j)  = 0.
                     trec(j)  = 0.
                     phirec(j) = 0.
                     rrec(j)   =0.

                     print*, "Not recorded"
                     print*, "total step is", k, "r is (AU)", r/AU
                     print*, "minimum distance is:", r_min/AU,"(AU)"
                     goto 20
              endif
           
              k = k + 1
              goto 10
           endif


14         continue         ! record particle, but no contirubtion
                            ! only adds to the denominator

           fp_intsec = 0.0  

15         continue


           if (j .eq. 1) then
               fpAU(i) = fp_intsec
           else
               fpAU(i) = (fpAU(i)*(j-1) + fp_intsec)/(j*1.)
           endif    

           if (fp_intsec .ne. 0.0) then

                jrec(j) = 1
                prec(j) = pnew
                mrec(j) = miunew
                trec(j) = tnew
                phirec(j) = phi_at_inner
                rrec(j) = rnew/AU

                print*, "recorded, miu = ", miunew, 'f(p) = ', fp_intsec

c                print*, "total step is", k
c                print*, "current r is:", rnew/AU, "(AU)"           else
           else 
                jrec(j)  = 0
                prec(j)  = 0.
                mrec(j)  = 0.
                trec(j)  = 0.
                phirec(j) = 0.
                rrec(j)   = 0.               
                
                print*, "No impact"
                ! recorded, but only contribute to denominator
           endif

20         continue ! now onto the next test particle
         enddo
         
         do j = 1, seed_num
           !record all the data
           write(13,*) i, jrec(j), prec(j), mrec(j), trec(j)
     1                      , phirec(j), rrec(j)
         enddo 

         write(11,*) i, fpAU(i), e0, fp0(p0)


       enddo  ! for i loop
       enddo  ! for itime loop

       close(11)
       close(13)

       if (if_parallel .eq. 1) then
          if (myRank .eq. 0) then
              close(12)
          endif
          call finishnodes
       endif

       END
       
c=======================================================================
c
c
c                   Additional Function and Subroutines
c
c
c=======================================================================
       
c========== Distribution function at inner boundary  ===================
      function fp0(p)
      real*8         p, fp0
      fp0 = p**(-4.5d0)
      return
      end
       
c========== Distribution function at inner boundary  ===================
      function fe0(e)
      real*8         e, fe0
      fe0 = e**(-2.75d0)
      return
      end
       
c========== Parker Field ===============================================
c---------- total B at r -----------------------------------------------

      function parker_field(r)
      real*8         parker_field, r, sec_psi

      include       'common.h'
	     
      sec_psi = dsqrt(1.0+(omega*r/usw)**2.0)
c      print*, 'sec(psi)=', sec_psi
      parker_field = br_bnd*(r_bnd/r)**2.* sec_psi

    
      return
      end

c================= Inverse Error Function ==============================
c--------------------------------------------------------------------
c--- The following function gives the inverse of error function
c--- define through erf(x) = 2/sqrt(pi) \int_0^x dt exp(-t^2).
c---------------------------------------------------------------------
      FUNCTION inv_erf(kesi)
      real*8         inv_erf,      kesi

      x = 1.
 11   continue
      if (abs( erf(x) - kesi) .gt. 1.E-4) then
         
         if (erf(x) .lt. kesi) then
            x = x * 1.2
         else
            x = x*0.8
         endif
         goto 11
      endif
      inv_erf = x

      return 
      end

c============ inverse of error function in double precision ============
c-----------------------------------------------------------
c     inverse of error function in double precision
c     erf^-1(x) = dierfc(1-x) 
c-----------------------------------------------------------
      FUNCTION dierfc(y)
      implicit real*8 (a-h, k-z)

      parameter (
     &    qa = 9.16461398268964d-01,
     &    qb = 2.31729200323405d-01,
     &    qc = 4.88826640273108d-01,
     &    qd = 1.24610454613712d-01,
     &    q0 = 4.99999303439796d-01,
     &    q1 = 1.16065025341614d-01,
     &    q2 = 1.50689047360223d-01,
     &    q3 = 2.69999308670029d-01,
     &    q4 = -7.28846765585675d-02)
      parameter (
     &    pa = 3.97886080735226000d+00,
     &    pb = 1.20782237635245222d-01,
     &    p0 = 2.44044510593190935d-01,
     &    p1 = 4.34397492331430115d-01,
     &    p2 = 6.86265948274097816d-01,
     &    p3 = 9.56464974744799006d-01,
     &    p4 = 1.16374581931560831d+00,
     &    p5 = 1.21448730779995237d+00,
     &    p6 = 1.05375024970847138d+00,
     &    p7 = 7.13657635868730364d-01,
     &    p8 = 3.16847638520135944d-01,
     &    p9 = 1.47297938331485121d-02,
     &    p10 = -1.05872177941595488d-01,
     &    p11 = -7.43424357241784861d-02)
      parameter (
     &    p12 = 2.20995927012179067d-03,
     &    p13 = 3.46494207789099922d-02,
     &    p14 = 1.42961988697898018d-02,
     &    p15 = -1.18598117047771104d-02,
     &    p16 = -1.12749169332504870d-02,
     &    p17 = 3.39721910367775861d-03,
     &    p18 = 6.85649426074558612d-03,
     &    p19 = -7.71708358954120939d-04,
     &    p20 = -3.51287146129100025d-03,
     &    p21 = 1.05739299623423047d-04,
     &    p22 = 1.12648096188977922d-03)
      z = y
      if (y .gt. 1) z = 2 - y
      w = qa - log(z)
      u = w**0.5
      s = (qc + log(u)) / w
      t = 1 / (u + qb)

      x = u * (1 - s * (0.5d0 + s * qd)) -
     &    ((((q4 * t + q3) * t + q2) * t + q1) * t + q0) * t
      t = pa / (pa + x)
      u = t - 0.5d0
      s = (((((((((p22 * u + p21) * u + p20) * u +
     &    p19) * u + p18) * u + p17) * u + p16) * u +
     &    p15) * u + p14) * u + p13) * u + p12
      s = ((((((((((((s * u + p11) * u + p10) * u +
     &    p9) * u + p8) * u + p7) * u + p6) * u + p5) * u +
     &    p4) * u + p3) * u + p2) * u + p1) * u + p0) * t -
     &    z * exp(x * x - pb)
      x = x + s * (1 + x * s)
      if (y .gt. 1) x = -x

      dierfc = x
      return
      end





c======== Calculate the turbulence level Delta_B^2 at r ================

      real*8  function deltaBs(r)       
      real*8         r, parker_field

      include       'common.h'

      deltaBs = Cturb_AU * parker_field(AU)**2. * (r/AU)**(-3.5)
      
      return
      end



     
c======== Calculate the normalization for Power spectrum ===============
c
c      The purpose for this function is to avoid repeating calculation 
c      for the normalization factor for different K at same R

      real*8  function Aq(r) 
      real*8         r,     total,        beta,  ks,    kl
      real*8         k,     kmo,   parker_field, lambdac
      integer        knum, ii
       
      include       'common.h'
      
      knum = 5000
      beta = 5./3.
      total = 0.0
      ks = ks_au
      kl = kl_au
      lambdac = lambdac_au
      
       do ii = 1, knum
         k   = exp(dlog(ks) + ii*(dlog(kl)-dlog(ks))/knum)
         kmo = exp(dlog(ks) + (ii-1)*(dlog(kl)-dlog(ks))/knum)
        
         total = total + lambdac *((1.+ (k*lambdac)**2.)**(-beta/2.) + 
     1        (1.0 +(kmo*lambdac)**2.)**(-beta/2.))/2. * (k-kmo)              
       enddo

       Aq = 1/ total

!      Aq = (lambdac**(1.-beta) / (1.-beta) * 
!     1        (kl**(1.-beta) - ks**(1.-beta)))**(-1.) 
      !Looks like Aq is going to blow to NaN if r is bigger than 1.806
       
      return
      end






c======== Calculate D_miumiu  ==========================================
      real*8  function Dmiumiu(miu,r,v, dDmmdm)
      implicit none
      real*8         miu, r, v, gyrof, waveK, dDmmdm, dkdmiu, dPkdk
      real*8         gamma0, parker_field, deltabs
c      real*8         sslab, s2d
       
      real*8         powerspec, ks, kl, lambdac, scale_factor
      real*8         Aq,  beta,  k,  total
      integer        knum, ii, jj
      real*8         Amass, Acharge
      
      include       'common.h'
     
      Amass = 56.
      Acharge = 14.

      beta = 5./3.
c      sslab  = 0.2   ! b_slab^2/b^2 put this into common.h later
c      s2d    = 0.8   ! b_2d^2/b^2 
      !  First Calculate the Power Spectrum P(k)

      scale_factor =  parker_field(AU)/parker_field(r) 
      lambdac = lambdac_au * scale_factor 
      ks      = ks_au / scale_factor
      kl      = kl_au / scale_factor

      gamma0 = 1./dsqrt(1.-(v/co)**2.)
      gyrof  = Acharge/Amass * eo * parker_field(r) /mp/gamma0

c      print*,'eo, parker_field(r), gamma0', eo, parker_field(r), gamma0 
      wavek  = gyrof*(v*dabs(miu))**(-1.)
       

      ! power spectrum

c      PowerSpec = Aq_AU* lambdac_au**(1.-beta) * waveK**(-1.*beta) 
c     1               * deltaBs(r)

       PowerSpec = Aq_AU * lambdac*(1.0+ (waveK*lambdac)**2.)
     1            **(-beta/2.) * deltaBs(r) *sslab ! changed in V1.05
                                            ! since it's deltaB_slab   
    

c      print*, 'miu, waveK, gyrof, PowerSpec, B', miu, waveK, gyrof, 
c     1        PowerSpec,    parker_field(r)

       Dmiumiu = pi/4.*(1.-miu**2.)* gyrof * waveK * PowerSpec
     1     / parker_field(r)**2.

c      Dmiumiu = pi/4.*(1.-miu**2)*gyrof**(2.-beta) * Aq_AU* (lambdac_au)
c     1    **(1.-beta) * deltaBs(r) * v **(beta-1.)/ parker_field(r)**2.
c     2    *((dabs(miu))**(beta-1.) + 0.2)

c------------  now calculate dDmiumiu/dmiu -----------------------------
       
      if (miu .ge. 0.) then
         dkdmiu = - gyrof /v/(miu)**2.
      else
         dkdmiu = gyrof /v/(miu)**2. 
      endif

c       dPkdk = - Aq_AU * lambdac **(beta+1.) * deltabs(r) *0.2 * beta* 
c     1        wavek**(beta-1.) *(1. + (waveK*lambdac)**beta)**(-2.) 
       dPkdk = - beta * PowerSpec * lambdac**2. * 
     1         (1.0+ (waveK*lambdac)**2.)**(-1.) * waveK

      dDmmdm = pi/4.*gyrof / parker_field(r)**2.0 * (
     1        -2.*miu * wavek * PowerSpec + (1. - miu**2.) * dkdmiu
     2        * PowerSpec + (1. - miu**2.)* wavek * dPkdk * dkdmiu  )    

      
C      dDmmdm = pi/4.*gyrof**(2.-beta) *Aq_AU * (lambdac_au)
C     1    **(1.-beta) * deltaBs(r) * v **(beta-1.)/ parker_field(r)**2.
C     2    *(-2.*miu*((dabs(miu))**(beta-1.)+0.2)+ (1.-miu**2.)*((beta
C     3     -1.)*(dabs(miu))**(beta-2.)*dabs(miu)/miu))

      return
      end



c======================================================================
c       real*8 function meanfreep(v, r)
       
c       include       'common.h'

c       real*8        v,     r,     beta,  const, Aq,    parker_field

c       beta = 5./3.

c       const  = pi/4.*Aq(r)*lambdac_au**(1.-beta) * Cturb_AU 
c     1        *(eo*parker_field(r)/mp)**(2.-beta)         
c       meanfreep = 0.75* v**(2.-beta) / const *2./(2.-beta)/(4.-beta)

c       return
c       end





c======================================================================

      SUBROUTINE SRANDPC(ISEED)
C
C  This subroutine sets the integer seed to be used with the
C  companion RAND function to the value of ISEED.  A flag is 
C  set to indicate that the sequence of pseudo-random numbers 
C  for the specified seed should start from the beginning.
C
      INTEGER NEXTN
      COMMON  /SEED/JSEED,IFRST, NEXTN
C
      JSEED = ISEED
      IFRST = 0
C
      RETURN
      END

      REAL *8 FUNCTION RANDpc()
      IMPLICIT REAL *8 (A-H, L-Z)      
C
C  This function returns a pseudo-random number for each invocation.
C  It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal 
C  standard number generator whose Pascal code appears in the article:
C
C     Park, Steven K. and Miller, Keith W., "Random Number Generators: 
C     Good Ones are Hard to Find", Communications of the ACM, 
C     October, 1988.
C
      INTEGER MPLIER, MODLUS, MOBYMP, MOMDMP
      PARAMETER (MPLIER=16807,MODLUS=2147483647,MOBYMP=127773,
     +           MOMDMP=2836)
C
      INTEGER HVLUE, LVLUE, TESTV, NEXTN
      COMMON  /SEED/JSEED, IFRST, NEXTN
c      SAVE    NEXTN
C
      IF (IFRST .EQ. 0) THEN
        NEXTN = JSEED
        IFRST = 1
      ENDIF
C
      HVLUE = NEXTN / MOBYMP
      LVLUE = MOD(NEXTN, MOBYMP)
      TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
      IF (TESTV .GT. 0) THEN
        NEXTN = TESTV
      ELSE
        NEXTN = TESTV + MODLUS
      ENDIF
      RANDPC = REAL(NEXTN)/REAL(MODLUS)
C
      RETURN
      END
      BLOCKDATA RANDBD
      COMMON  /SEED/JSEED,IFRST, NEXTN
C
      DATA JSEED,IFRST/123456789,0/
C
      END


c======================================================================
c==== standard R-K methods
c==== http://www.myphysicslab.com/runge_kutta.html
c== xn+1 = xn + h⁄6 (a + 2 b + 2 c + d) 
c== where a = f (tn, xn)
c== b = f (tn + h⁄2, xn + h⁄2 a)
c== c = f (tn + h⁄2, xn + h⁄2 b)
c== d = f (tn + h, xn + h c)
c======================================================================
       Subroutine RK4(r, miu, p, h, rnew, miunew, pnew, vnew)
       implicit none
       real *8       r,     miu,   p,     v,     h,     rnew,  miunew,
     1               pnew,  vnew
       real *8       e_p,   gamma_p,      ds,    dr,    sec_psi,      
     1               LZ
       real *8       RKr1,  RKp1,  RKmiu1,       r1,    p1,    miu1,
     1               RKr2,  RKp2,  RKmiu2,       r2,    p2,    miu2,
     2               RKr3,  RKp3,  RKmiu3,       r3,    p3,    miu3,
     3               RKr4,  RKp4,  RKmiu4

       real *8       parker_field
       
       include       'common.h'

c=== first get (Rkr1, Rkz1)
       
       e_p = (dsqrt((mp*co**2.0)**2.0 + (p*co)**2.0)-mp*co**2.0) /eo
       gamma_p = (e_p*eo + mp*co**2.0)/(mp*co**2.0)
       v      = p/mp/gamma_p              
       sec_psi = dsqrt(1.0+(omega*r/usw)**2.0)
       
       ds     = v * h * miu
       dr     = ds / sec_psi


       LZ     = - parker_field(r) * ds/4.
     1           /(parker_field(r + dr/4.) - parker_field(r))

                

       RKr1 = v * miu / sec_psi
       RKp1 = - p * usw *(sec_psi/2./LZ*(1.-miu**2.) +
     1             1./ sec_psi * (omega/usw)**2.*
     2             r /dsqrt(1.+(omega* r /usw)**2.)* miu**2.)
       RKmiu1  = v/(2.*LZ)*(1. + miu * usw/v*sec_psi - miu*usw*v
     1             /co**2. * sec_psi)*(1.-miu**2.) - usw *(1./sec_psi
     2             * (omega/usw)**2.*r/dsqrt(1.+(omega*r/usw)**2.))
     3             *miu *(1.-miu**2.)

c=== Next get (Rkr2, Rkz2)

       r1 = r + h/2. * RKr1
       p1 = p + h/2. * RKp1
       miu1 = miu + h/2. * RKmiu1

c      Bmag = getBmag(r1/AU, u_sw, Br, Bphi) 
c      drds = Br/Bmag
c      RKmu2 = v * (1-mu1**2)/r1 * drds      
       e_p = (dsqrt((mp*co**2.0)**2.0 + (p1*co)**2.0)-mp*co**2.0) /eo
       gamma_p = (e_p*eo + mp*co**2.0)/(mp*co**2.0)
       v      = p1/mp/gamma_p   
       sec_psi = dsqrt(1.0+(omega*r1/usw)**2.0)
       LZ     = - parker_field(r1) * ds/4.
     1           /(parker_field(r1 + dr/4.) - parker_field(r1))
       
       RKr2 =  v * miu1 / sec_psi
       RKp2 =  - p1 * usw *(sec_psi/2./LZ*(1.-miu1**2.) +
     1             1./ sec_psi * (omega/usw)**2.*
     2             r1 /dsqrt(1.+(omega* r1 /usw)**2.)* miu1**2.)
       Rkmiu2 = v/(2.*LZ)*(1. + miu1 * usw/v*sec_psi - miu1 *usw*v
     1             /co**2. * sec_psi)*(1.-miu1**2.) - usw *(1./sec_psi
     2             * (omega/usw)**2.*r1 /dsqrt(1.+(omega*r1 /usw)**2.))
     3             *miu1 *(1.-miu1**2.)

c=== Next get (Rkr3, Rkz3)

       r2 = r + h/2. * RKr2
       p2 = p + h/2. * RKp2
       miu2 = miu + h/2. * RKmiu2

       e_p = (dsqrt((mp*co**2.0)**2.0 + (p2*co)**2.0)-mp*co**2.0) /eo
       gamma_p = (e_p*eo + mp*co**2.0)/(mp*co**2.0)
       v      = p2/mp/gamma_p   
       sec_psi = dsqrt(1.0+(omega*r2/usw)**2.0)
       LZ     = - parker_field(r2) * ds/4.
     1           /(parker_field(r2 + dr/4.) - parker_field(r2))
       
       RKr3 =  v * miu2 / sec_psi
       RKp3 =  - p2 * usw *(sec_psi/2./LZ*(1.-miu2**2.) +
     1             1./ sec_psi * (omega/usw)**2.*
     2             r2 /dsqrt(1.+(omega* r2 /usw)**2.)* miu2**2.)
       Rkmiu3 = v/(2.*LZ)*(1. + miu2 * usw/v*sec_psi - miu2 *usw*v
     1             /co**2. * sec_psi)*(1.-miu2**2.) - usw *(1./sec_psi
     2             * (omega/usw)**2.*r2 /dsqrt(1.+(omega*r2 /usw)**2.))
     3             *miu2 *(1.-miu2**2.)

c=== Next get (Rkr4, Rkz4)
       r3 = r + h * RKr3
       p3 = p + h * RKp3
       miu3 = miu + h * RKmiu3

       e_p = (dsqrt((mp*co**2.0)**2.0 + (p3*co)**2.0)-mp*co**2.0) /eo
       gamma_p = (e_p*eo + mp*co**2.0)/(mp*co**2.0)
       v      = p3/mp/gamma_p   
       sec_psi = dsqrt(1.0+(omega*r3/usw)**2.0)
       LZ     = - parker_field(r3) * ds/4.
     1           /(parker_field(r3 + dr/4.) - parker_field(r3))
       
       RKr4 =  v * miu3 / sec_psi
       RKp4 =  - p3 * usw *(sec_psi/2./LZ*(1.-miu3**2.) +
     1             1./ sec_psi * (omega/usw)**2.*
     2             r3 /dsqrt(1.+(omega* r3 /usw)**2.)* miu3**2.)
       Rkmiu4 = v/(2.*LZ)*(1. + miu3 * usw/v*sec_psi - miu3 *usw*v
     1             /co**2. * sec_psi)*(1.-miu3**2.) - usw *(1./sec_psi
     2             * (omega/usw)**2.*r3 /dsqrt(1.+(omega*r3 /usw)**2.))
     3             *miu3 *(1.-miu3**2.)

c=================Now get the new (r, mu, p) ===================

c      print*, "=== RKr1, RKr2, RKr3, RKr4 ===", RKr1, RKr2, RKr3, RKr4
c      print*, 'r1,r2,r3', r1,r2,r3
c      print*, 'p1,p2,p3', p1,p2,p3
c      print*, 'm1,m2,m3', miu1,miu2,miu3
c      print*, "=== RKmu1, RKmu2, RKmu3, RKmu4 ===", 
c     X     RKmiu1, RKmiu2, RKmiu3, RKmiu4
c      print*, "=== in RK4 , r, mu, v, p ===", r, miu, v, p


       rnew = r + h/6. * ( RKr1 + 2.*RKr2 + 2.*RKr3 + RKr4)
       pnew = p + h/6. * ( RKp1 + 2.*RKp2 + 2.*RKp3 + RKp4)
       miunew = miu + h/6. * ( RKmiu1 + 2.*RKmiu2 + 2.*RKmiu3 + Rkmiu4)

       e_p = (dsqrt((mp*co**2.0)**2.0 + (pnew*co)**2.0)-mp*co**2.0) /eo
       gamma_p = (e_p*eo + mp*co**2.0)/(mp*co**2.0)
       vnew      = pnew/mp/gamma_p   

c      print*, "v, vnew, Omega", v, vnew, Omega
c      pause
       end



c=======================================================================
c=============  get the shock location at time T  ======================
       Subroutine get_rsh(sys_time, phinew, shell_time, shell_loc, 
     1               bndy_loc_in, rsh, rsh_in, ip_shell, im_shell, 
     2               ip_phi, im_phi)

c      output are rsh, rsh_in

       implicit none
       include       'common.h'
       real*8        sys_time,  pnew, phinew, shell_time(shell_no), 
     1               shell_loc(shell_no, phi_no),
     2               bndy_loc_in(shell_no, phi_no ),
     3               footpoint_phi,       rsh,   rsh_in
     
       real*8        time
       Integer       i, j, k, ip_shell, im_shell, ip_phi, im_phi
       


! 	print*, "shell_time(ip_shell), time, shell_time(im_shell)"
!      1        , shell_time(ip_shell), time, shell_time(im_shell)
!        print*, "ip_phi, im_phi, loc(ip_shell), loc(im-shell)",
!      1          ip_phi, im_phi, shell_loc(ip_shell, ip_phi),
!      2          shell_loc(im_shell, ip_phi), phi_e_at_inner, phi_0(1) 

       time = sys_time

       rsh=((shell_time(ip_shell) -time) *shell_loc(ip_shell, ip_phi)
     1    + (time - shell_time(im_shell)) * shell_loc(im_shell, ip_phi))
     2     /(shell_time(ip_shell) - shell_time(im_shell)) * AU
c          print*, 'rsh = ', rsh/AU, 'AU'

       rsh_in=((shell_time(ip_shell) -time)
     1               *bndy_loc_in(ip_shell, ip_phi)
     1    + (time - shell_time(im_shell)) 
     2               * bndy_loc_in(im_shell, ip_phi))
     2     /(shell_time(ip_shell) - shell_time(im_shell)) * AU

c     	if (rsh .ge. AU ) then      ! this part is removed in #1.04
c     	   rsh = r_inner
c     	else
                     
c     	endif

       ! need average over phi?
       end


c=======================================================================
c=============  get fp at time T  ======================


       real*8 Function getft(pnew, tnew, rnew, shell_time, shell_loc, 
     1       phinew, dist_fp_up, ip_shell, im_shell, ip_phi, im_phi)

       implicit none
       include       'common.h'
       real*8        pnew,  tnew,  rnew, phinew, shell_time(shell_no), 
     1               shell_loc(shell_no, phi_no ),
     2               dist_fp_up(shell_no, p_number, phi_no)
       Integer       ip_shell, im_shell, ip_phi, im_phi, ip_p, im_p 
       real*8        tm, tp, phim, phip, pm, pp, time
       integer       i,     j

       time = time_fp + tnew /t_o   

c-------- debug      ---------------------------------------------------
c       print*, 'rnew, tnew', rnew/AU, tnew  
c       print*, 'ip_shell, im_shell', ip_shell, im_shell
c       print*,  phi_0(phi_no) + omega/usw*(shell_loc(ip_shell,phi_no)*AU
c     1                   - r_bnd)*180./pi  , phi_e_at_inner
c       print*, shell_loc(ip_shell, phi_no-1)
c     1   , shell_loc(im_shell, phi_no-1)
c       print*, 'im_phi, ip_phi', im_phi, ip_phi
c-----------------------------------------------------------------------

       tm     = shell_time(im_shell)
       tp     = shell_time(ip_shell)

c       print*, 'phinew in getft', phinew
       phim   = phi_0(im_phi)
       phip   = phi_0(ip_phi)
       ! check again to see if phinew is in between phim, phip
       if (phinew .gt. phip .and. ip_phi .lt. phi_no) then
           im_phi = im_phi + 1
           ip_phi = ip_phi + 1
       else if (phinew .lt. phim .and. im_phi .gt. 1) then
           im_phi = im_phi - 1
           ip_phi = ip_phi - 1
       endif
       phim   = phi_0(im_phi)
       phip   = phi_0(ip_phi)
c       print*, '2done'
       

       do i = 1, p_number
          if (p_0(i) .ge. pnew) then
             ip_p = i
             im_p = i-1
             goto 37
          endif       
       end do 
37     continue

       pm     = p_0(im_p)
       pp     = p_0(ip_p)
c       print*, '3done'

       getft = (dist_fp_up (im_shell, im_p, im_phi) * (tp - time) *
     1               (pp - pnew) * (phip - phinew) +   
     2          dist_fp_up (im_shell, im_p, ip_phi) * (tp - time) *
     1               (pp - pnew) * (phinew- phim) +
     2          dist_fp_up (im_shell, ip_p, im_phi) * (tp - time) *
     1               (pnew - pm) * (phip - phinew) +
     2          dist_fp_up (ip_shell, im_p, im_phi) * (time - tm) *
     1               (pp - pnew) * (phip - phinew) +
     2          dist_fp_up (im_shell, ip_p, ip_phi) * (tp - time) *
     1               (pnew - pm) * (phinew - phim) +
     2          dist_fp_up (ip_shell, im_p, ip_phi) * (time - tm) *
     1               (pp - pnew) * (phinew - phim) +
     2          dist_fp_up (ip_shell, ip_p, im_phi) * (time - tm) *
     1               (pnew - pm) * (phip - phinew) +
     2          dist_fp_up (ip_shell, ip_p, ip_phi) * (time - tm) *
     1               (pnew - pm) * (phinew - phim) )
     2        / ((pp - pm)*(tp - tm)*(phip - phim))  

c       print*, '4done'
       if (getft .lt. 0.0) then
         print*, "warning! getft<0 !"
         print*,"phim, phip, phinew", phim, phip, phinew
         print*,"getft=", getft
         print*, "tp, tm, time", tp, tm, time
         print*," pp, pm, pnew", pp, pm, pnew
c         print*, dist_fp_up (im_shell, im_p, im_phi)
c         print*, dist_fp_up (im_shell, im_p, ip_phi)
c         print*, dist_fp_up (im_shell, ip_p, im_phi)
c         print*, dist_fp_up (ip_shell, im_p, im_phi)
c         print*, dist_fp_up (ip_shell, ip_p, im_phi)
c         print*, dist_fp_up (im_shell, ip_p, ip_phi)
c         print*, dist_fp_up (ip_shell, im_p, ip_phi)
c         print*, dist_fp_up (ip_shell, ip_p, ip_phi)
c         pause
       endif
       
       return
       end









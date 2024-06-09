       min_e = 1.0D5;       max_e = 1.0D9;       mp = 1.67D-27

       AU = 1.5D11;         eo = 1.6D-19
       pi   = 3.141592653589793116d0
       bo = 2.404D-9;       t_o = 2858068.3;     vo = 52483.25 
       co = 3.0D8;          n_0 = 1.0D6
     
       ! parameter for specific solar wind settings    
       omega = 2.87D-6;     usw = 442727;        br_bnd =476.5*bo
       r0  = 1.0 * AU;      r_bnd = 0.05 * AU

       ! cycle 23: usw = 8.116* vo;  br_bnd = 196.90*bo
       ! cycle 24: usw = 7.666* vo;  br_bnd = 145.85*bo
       ! low_B: usw = 8.116  , br = 152.43

       ! omega --- rotation angular speed of the sun
       ! usw   --- solar wind speed (m/s)
       ! br_bnd--- radial magnetic field at inner bndry (T)
       ! r0    --- observer location (m)
       ! r_bnd --- inner boundary

       ! parameter for turbulence:
       cturb_au = 0.5;             lambdac_au = 0.03*AU
       kl_au = 2.d-7;              ks_au = 1.d-10
       l2d   = 0.003*AU;           sslab =0.2  ! b_slab^2/b^2
       s2d = 0.8     ! b_2d^2/b^2 

       !r_inner = 0.025*AU

       !Zeus parameters
       min_p = 1.0;         max_p = 11.0;        maxshells = 200
       cme_center = 100;         cme_width = 120
       del_phi = 5.0

       del_p   = (max_p - min_p) / (p_number*1.0)
       Aq_AU   = Aq(AU)
       
c       Aq_AU =  1./2./dsqrt(pi)*gamma(5./6.)/gamma(1./3.)

c       !time for the interested output f(p)

       time_start = 0.40
       time_end   = 0.47
       !earth's longitudinal position on the ecliptic plane
       phi_e = 100.0
       ! this is used to mark the field line passing the observer
       
       phi_e_at_inner = phi_e + omega/usw* (r0 - r_bnd)*180./pi
       
       if (phi_e_at_inner .ge. 360.) then
              phi_e_at_inner = phi_e_at_inner - 360.*floor(
     1        phi_e_at_inner/360.)
       endif

       print*, "iiiiiiiiiiiiiii",phi_e_at_inner

       if_parallel = 1

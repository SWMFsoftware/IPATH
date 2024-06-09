c--------------------    COMMON BLOCK DECLARATION    -------------------

       integer       p_num,        seed_num,     iloop
       integer       p_number,     phi_no,       shell_no    

       parameter     (p_num  = 40,  seed_num = 100,  iloop  = 60000)
       parameter     (p_number = 400, phi_no = 25)

       real*8        mp,           AU,           pi,           eo,
     1               bo,           t_o,          vo,           co,
     2               n_0
       real*8        lambdac_au,   kl_au,        ks_au,        Cturb_AU,
     1               Aq_AU,        l2d,          sslab,        s2d
       real*8        omega,        usw,          br_bnd,       r_bnd,
     1               r0 

       real*8        min_e,        max_e,        r_inner

       real*8        max_p,        min_p,        maxshells,    del_p,
     1               cme_center,   cme_width,    del_phi
       real*8        time_fp,      phi_e,        phi_e_at_inner,
     1               time_start,   time_end,     phi_at_inner

       real*8        energy_0(p_number),  phi_0(phi_no), p_0(p_number)
               

       COMMON / NORMAL /
     1               mp,           AU,           eo,           pi,
     2               bo,           t_o,          vo,           co,
     3               n_0,          min_e,        max_e,        r_inner,
     4               time_fp,      phi_e,        phi_e_at_inner,
     5               time_start,   time_end,     phi_at_inner
       COMMON / SWIND /
     1               omega,        usw,          br_bnd,       r_bnd,
     2               r0  
       COMMON / TURBULENCE /
     1               lambdac_au,   kL_au,        ks_au,        Cturb_AU,
     2               Aq_AU,        l2d,          sslab,        s2d                 
       COMMON / ZEUS /
     1               maxshells,    max_p,        min_p,        del_p,
     2               cme_center,   cme_width,    del_phi
       COMMON / DISTR /
     1               energy_0,     phi_0,        p_0,          shell_no

C---------------     END OF COMMON BLOCK DECLARATION     ---------------

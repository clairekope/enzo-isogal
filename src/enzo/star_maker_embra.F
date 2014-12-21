#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c///////////////////// SUBROUTINE STAR_MAKER_EMBRA \\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_maker_embra(nx, ny, nz,
     &                      d, dm, temp, u, v, w, cooltime,
     &                      dt, r, metal, dx, t,
     &                      d1, x1, v1, t1,
     &                      nmax, xstart, ystart, zstart, ibuff, 
     &                      imetal, imethod, 
     &                      odthresh, dthresh, masseff,
     &                      np, xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf, tip,
     &                      cur_np, cur_xp, cur_yp, cur_zp,
     &                      cur_up, cur_vp, cur_wp,
     &                      cur_mp, cur_tcp, cur_metal, cur_tip)

c
c  CREATE EMBRA STAR PARTICLES
c
c  written by: Britton Smith
c  date:       December 2014
c
c  Embra particles are formed at a rate given by the KS relation
c  where the mass that goes into a star particle is:
c  m_star = m_cell * c * dt / tdyn,
c  where c is the star formation efficiency.
c
c  Star particles grow in mass over multiple timesteps until reaching 
c  a final mass of (m_cell * c), at which time they are turned on.
c  When a star is turned on, it emits early stellar feedback as 
c  thermal energy for a period of time until it explodes.
c  All supernova feedback is injected in a single timestep to yield 
c  a given change in temperature, dT.
c
c  INPUTS:
c
c    nx,ny,nz - dimensions of field arrays
c    d        - density field
c    dm       - dark matter field
c    temp     - temperature field
c    u,v,w    - velocity fields
c    cooltime - cooling time in code units
c    dt       - current timestep
c    r        - refinement field (non-zero if zone is further refined)
c    metal    - metal field (converted to a fraction)
c    dx       - zone size (code units)
c    t        - current time (code units)
c    d1,x1,v1,t1 - density, length, velocity, and time units
c    x/y/z start - starting position of grid origin
c    ibuff    - number of buffer zones at each end of grid
c    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS)
c    odthresh - overdensity threshold
c    dthresh  - density threshold in proper CGS
c    masseff  - fraction of gas mass converted into stars
c
c  OUTPUTS:
c
c    np      - number of particles created
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle
c    metalf   - metallicity fraction of particle
c    tip      - particle ignition time
c    nmax     - particle array size specified by calling routine
c    cur_np   - current number of particles in grid
c    cur_xp, cur_yp, cur_zp, 
c    cur_up, cur_vp, cur_wp,
c    cur_mp, cur_tcp, cur_metal, 
c    cur_tip  - position, velocity, mass, creation time, 
c               metal fraction, and ignition time of current particles
c
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, nmax, np, imetal, imethod,
     &      cur_np
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), cooltime(nx,ny,nz), metal(nx,ny,nz)
      R_PREC    dt, dx
      R_PREC    d1, x1, v1, t1
      P_PREC    xstart, ystart, zstart, t
      P_PREC    xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax),
     &     mp(nmax), tdp(nmax), tcp(nmax), 
     &     metalf(nmax), tip(nmax)
      P_PREC    cur_xp(cur_np), cur_yp(cur_np), cur_zp(cur_np)
      R_PREC    cur_up(cur_np), cur_vp(cur_np), cur_wp(cur_np),
     &     cur_mp(cur_np), cur_tcp(cur_np), cur_metal(cur_np),
     &     cur_tip(cur_np)
      R_PREC    odthresh, dthresh, masseff
c
c  Locals:
c
      INTG_PREC  i, j, k, ii, n, newst
      R_PREC   div, tdyn, dtot
      R_PREC   targetmass, newstarmass, newmassthisstar
      R_PREC   starid(nx,ny,nz)
      R_PREC   pi, G
      parameter (pi=pi_val, G=GravConst)
c
c     Check to see if a grid cell contains a star that has not
c     reached its final mass.
c
      do k=1+ibuff,nz-ibuff
         do j=1+ibuff,ny-ibuff
            do i=1+ibuff,nx-ibuff
               starid(i,j,k) = -1
            enddo
         enddo
      enddo

      do n=1,cur_np
c
c     Check to see if this particle is a star.
c
         if (cur_tcp(n) .lt. 0._RKIND) goto 34
c
c     If star has already turned on, skip it.
c
         if (cur_tip(n) .ge. 0._RKIND) goto 34
c
c     Compute index of the cell that the star particle resides in.
c
         i = int((cur_xp(n) - xstart)/dx,IKIND) + 1
         j = int((cur_yp(n) - ystart)/dx,IKIND) + 1
         k = int((cur_zp(n) - zstart)/dx,IKIND) + 1
c
c     If star has reached target mass, turn it on.
c
         targetmass = masseff * d(i,j,k)
         if (cur_mp(n) .ge. targetmass) then
            cur_tip(n) = t
            goto 34
         endif

         starid(i,j,k) = n

 34      continue

      enddo

      ii = np
c
c  For each zone, a star particle is created if answers to all the
c  following questions are affirmative:
c
c    is this the finest level of refinement ?
c    is the density/overdensity greater than a critical density/overdensity ?
c    is the flow convergent ?
c    is the cooling time less than a dynamical time ? 
c
c  Before forming a new star, we check to see if a star already exists in
c  that grid cell that has not reached its target mass.  If so, we put the 
c  stellar mass into that particle.
c
      do k=1+ibuff,nz-ibuff
         do j=1+ibuff,ny-ibuff
            do i=1+ibuff,nx-ibuff
c
c              0) Do not generate more star particles than available
c
               if (ii .eq. nmax) goto 20
c
c              1) is this finest level of refinement?
c
               if (r(i,j,k) .ne. 0._RKIND) goto 10
c
c              2) is density greater than density threshold?
c
               if (d(i,j,k)*d1 .lt. dthresh) goto 10
c
c              2.5) is density greater than overdensity threshold?
c
               if (d(i,j,k) .lt. odthresh) goto 10
c
c              3) is divergence negative?
c                 (the first calculation is face centered for ZEUS, 
c                  the second is cell-centered for PPM)
c
               if (imethod .eq. 2) then
                  div = u(i+1,j  ,k  ) - u(i,j,k)
     &                + v(i  ,j+1,k  ) - v(i,j,k)
     &                + w(i  ,j  ,k+1) - w(i,j,k)
               else
                  div = u(i+1,j  ,k  ) - u(i-1,j  ,k  )
     &                + v(i  ,j+1,k  ) - v(i  ,j-1,k  )
     &                + w(i  ,j  ,k+1) - w(i  ,j  ,k-1)
               endif
               if (div .ge. 0._RKIND) goto 10
c
c              4) t_cool < t_free-fall (if T < 1.1e4 skip this check)
c
               dtot = ( d(i,j,k) + dm(i,j,k) )*d1
               tdyn  = sqrt(3._RKIND*pi/32._RKIND/G/dtot)/t1

               if (tdyn .lt. cooltime(i,j,k) .and. 
     &             temp(i,j,k) .gt. 1.1e4_RKIND) goto 10
c
c              Calculate amount of mass to convert into stars.
c
               targetmass = masseff * d(i,j,k)
c
c              Use at most half the cell mass
c
               newstarmass = min((targetmass * dt / tdyn),
     &              (0.5_RKIND * d(i,j,k)))
c
c              5) Check if a star already exists at this location.
c
               if (starid(i,j,k) .gt. 0) then
                  n = starid(i,j,k)
c
c                 Calculate how much mass to add to un-ignited star.
c
                  newmassthisstar = min(newstarmass,
     &                 (targetmass - cur_mp(n)))
                  newmassthisstar = max(newmassthisstar, 0._RKIND)

c
c                 Add metals to star.
c
                  if (imetal .eq. 1) then
                     cur_metal(n) = (cur_metal(n) * cur_mp(n) +
     &                    metal(i,j,k) * newmassthisstar) /
     &                    (cur_mp(n) + newmassthisstar)
                  endif
c
c                 Add momentum.
c
                  if (imethod .ne. 2) then
                     cur_up(n) = (cur_mp(n) * cur_up(n) +
     &                    newmassthisstar * u(i,j,k)) /
     &                    (cur_mp(n) + newmassthisstar)
                     cur_vp(n) = (cur_mp(n) * cur_vp(n) +
     &                    newmassthisstar * v(i,j,k)) /
     &                    (cur_mp(n) + newmassthisstar)
                     cur_wp(n) = (cur_mp(n) * cur_wp(n) +
     &                    newmassthisstar * w(i,j,k)) /
     &                    (cur_mp(n) + newmassthisstar)
                  else
                     cur_up(n) = (cur_mp(n) * cur_up(n) +
     &                    newmassthisstar * 0.5 *
     &                    (u(i,j,k) + u(i+1,j,k))) /
     &                    (cur_mp(n) + newmassthisstar)
                     cur_vp(n) = (cur_mp(n) * cur_vp(n) +
     &                    newmassthisstar * 0.5 *
     &                    (v(i,j,k) + v(i,j+1,k))) /
     &                    (cur_mp(n) + newmassthisstar)
                     cur_wp(n) = (cur_mp(n) * cur_wp(n) +
     &                    newmassthisstar * 0.5 *
     &                    (w(i,j,k) + w(i,j,k+1))) /
     &                    (cur_mp(n) + newmassthisstar)
                  endif

c
c                 Add mass to star and remove from grid cell 
c                 and total star formation budget.
c
                  cur_mp(n) = cur_mp(n) + newmassthisstar
                  d(i,j,k) = d(i,j,k) - newmassthisstar
                  newstarmass = newstarmass - newmassthisstar
c
c                 If target mass reached, turn this star on.
c
                  if (((targetmass - cur_mp(n)) / targetmass) .le.
     &                 1.e-3_RKIND) then
                     cur_tip(n) = t
                  endif

               endif
c
c              Create new star particles until mass is gone
c
               newst = int(newstarmass / targetmass) + 1
               do n=1,newst

                  if (newstarmass .le. 0._RKIND) goto 10

                  ii = ii + 1
                  newmassthisstar = min(newstarmass, targetmass)
                  mp(ii)  = newmassthisstar
                  tcp(ii) = t
                  tdp(ii) = tdyn
                  xp(ii) = xstart + (REAL(i,RKIND)-0.5_RKIND)*dx
                  yp(ii) = ystart + (REAL(j,RKIND)-0.5_RKIND)*dx
                  zp(ii) = zstart + (REAL(k,RKIND)-0.5_RKIND)*dx

c
c                 If target mass reached, turn this star on.
c
                  if (((targetmass - mp(ii)) / targetmass) .le.
     &                 1.e-3_RKIND) then
                     tip(ii) = t
                  endif

                  if (imethod .eq. 2) then
                     up(ii) = 0.5_RKIND*(u(i,j,k)+u(i+1,j,k))
                     vp(ii) = 0.5_RKIND*(v(i,j,k)+v(i,j+1,k))
                     wp(ii) = 0.5_RKIND*(w(i,j,k)+w(i,j,k+1))
                  else
                     up(ii) = u(i,j,k)
                     vp(ii) = v(i,j,k)
                     wp(ii) = w(i,j,k)
                  endif

c
c                 Set the particle metal fraction
c
                  if (imetal .eq. 1) then
                     metalf(ii) = metal(i,j,k) ! in here metal is a fraction
                  else
                     metalf(ii) = 0._RKIND
                  endif
c
c                 Remove mass from grid
c
                  d(i,j,k) = d(i,j,k) - newmassthisstar
                  newstarmass = newstarmass - newmassthisstar

               enddo

 10            continue

            enddo
         enddo
      enddo
 20   continue
c
      if (ii .ge. nmax) then
         write(6,*) 'star_maker_embra: reached max new particle count'
         ERROR_MESSAGE
      endif
      np = ii

      return
      end
c
c=======================================================================
c////////////////// SUBROUTINE STAR_FEEDBACK_EMBRA \\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback_embra(nx, ny, nz,
     &                      d, te, ge, u, v, w, metal,
     &                      idual, imetal, imethod,
     &                      dt, dx, t,
     &                      d1, x1, v1, t1, temp1, gamma,
     &                      m_eject, yield,
     &                      tdelay, deltatemp, iearly, uvlum,
     &                      npart, xstart, ystart, zstart,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf, tip, type)

c
c  RELEASES STAR PARTICLE ENERGY, MASS AND METALS
c
c  written by: Britton Smith
c  date:       December 2014
c
c  INPUTS:
c
c    nx,ny,nz - dimensions of field arrays
c    d        - density field
c    te,ge    - total energy and gas energy fields
c    u,v,w    - velocity fields
c    metal    - metallicity density field
c    idual    - dual energy flag
c    imetal   - metallicity flag (0 - none, 1 - yes)
c    imethod  - hydro method (0 - PPMDE, 1 - PPMLR, 2 - ZEUS)
c    dt       - current timestep
c    dx       - zone size (code units)
c    t        - current time (code units)
c    d1,x1,v1,t1,temp1 - density, length, velocity, time, and
c                        temperature units
c    gamma    - adiabatic index
c    m_eject  - fraction of stellar mass ejected back to gas
c    yield    - fraction of stellar mass that is converted to metals
c    tdelay   - time delay between star ignition and supernova feedback (in years)
c    deltatemp - temperature change resulting from supernova feedback (in K)
c    iearly   - flag to enable early stellar feedback
c    uvlum    - uv luminosity from early stellar feedback (in erg/s/Msun)
c    npart    - particle array size specified by calling routine
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of particles
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle (-1 if not a star particle)
c    metalf   - star particle metal fraction
c    tip      - ignitino time of particles (when star turns on)
c    type     - particle type
c
c  OUTPUTS:
c    d,u,v,w,te,ge - modified fields
c
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, npart, idual, imetal, imethod,
     &      iearly
      R_PREC    d(nx,ny,nz), te(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    metal(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1, temp1, gamma
      R_PREC    tdelay, deltatemp, uvlum
      P_PREC    xstart, ystart, zstart, t
      P_PREC    xp(npart), yp(npart), zp(npart)
      R_PREC    up(npart), vp(npart), wp(npart)
      R_PREC    mp(npart), tdp(npart), tcp(npart),
     &     metalf(npart), tip(npart)
      INTG_PREC type(npart)
c
c  Locals
c
      INTG_PREC i, j, k, n
      R_PREC    m_eject, yield, dratio, distmass, tdelaycode, 
     &     snenergy, p_age, uv_dt, en_early
      LOGIC_PREC tooyoung
c-----------------------------------------------------------------------
c
c     Convert delay time to code units.
c
      tdelaycode = tdelay * yr / t1
c
c     Calculate feedback energy corresponding to desired temperature change.
c
      snenergy = deltatemp / (gamma - 1._RKIND) / temp1
c
c     Loop over particles
c
      do n=1, npart
         if (tcp(n) .gt. 0 .and. mp(n) .gt. 0 .and. type(n) .eq. 2) then
c
c           Skip particles that have not turned on.
c
            if (tip(n) .lt. 0._RKIND) goto 10
c
c           Calculate particle age (since turning on.)
c
            p_age = t - tip(n)
c
c           Skip if supernova already occurred.
c
            if (p_age .gt. tdelaycode) goto 10
c
c           Check if SN feedback will occur this timestep.
c
            if (p_age + dt .lt. tip(n)) then
               tooyoung = .true.
            else
               tooyoung = .false.
            endif
c
c           Skip if star too young and no early feedback.
c
            if ((iearly .eq. 0) .and. tooyoung) goto 10
c
c           Compute index of the cell that the star particle
c           resides in.
c 
            i = int((xp(n) - xstart)/dx,IKIND) + 1
            j = int((yp(n) - ystart)/dx,IKIND) + 1
            k = int((zp(n) - zstart)/dx,IKIND) + 1
c
c           Check bounds - if star particle is outside of this grid
c           then exit and give a warning.
c
            if (i .lt. 1 .or. i .gt. nx .or. j .lt. 1 .or. j .gt. ny
     &          .or. k .lt. 1 .or. k .gt. nz) then
               write(6,*) 'warning: star particle out of grid',i,j,k
               goto 100
            endif
c
c           Inject thermal energy for early stellar feedback.
c
            if (iearly .eq. 1) then
               uv_dt = min(dt, (tdelaycode - p_age))
               uv_dt = max(uv_dt, 0._RKIND)

               en_early = uvlum * mp(n) / SolarMass *
     &              uv_dt * t1 / d(i,j,k) / v1**2

               te(i,j,k) = te(i,j,k) + en_early
               if (idual .eq. 1) then
                  ge(i,j,k) = ge(i,j,k) + en_early
               endif
            endif
c
            if (tooyoung) goto 10
c
c           Calculate ejected mass from particle.
c
            distmass = mp(n) * m_eject
c
c           Subtract ejected mass from particle (ejection due
c           to winds, supernovae.)
c
            mp(n) = mp(n) * (1.0_RKIND - m_eject)
c
c           Add energy to energy field
c
            te(i,j,k) = te(i,j,k) + snenergy
            if (idual .eq. 1) then
               ge(i,j,k) = ge(i,j,k) + snenergy
            endif
c
c           Metal feedback (note that in this function gas metal is
c             a fraction (rho_metal/rho_gas) rather than a density.
c             The conversion has been done in the handling routine)
c
            if (imetal .eq. 1) then
c
c           "Cen method".  This takes into account gas recycling.
c
               dratio = 1._RKIND / (d(i,j,k) + distmass)
               metal(i,j,k) = 
     &              (metal(i,j,k) * d(i,j,k) + 
     &              mp(n) * 
     &              (yield * (1._RKIND-metalf(n)) + 
     &              m_eject * metalf(n))) * dratio
            endif
c
c           Mass and momentum feedback
c
            u(i,j,k) = u(i,j,k)*d(i,j,k) +
     &           distmass * up(n)
            v(i,j,k) = v(i,j,k)*d(i,j,k) +
     &           distmass * vp(n)
            w(i,j,k) = w(i,j,k)*d(i,j,k) +
     &           distmass * wp(n)
            d(i,j,k) = d(i,j,k) + distmass
            u(i,j,k) = u(i,j,k)/d(i,j,k)
            v(i,j,k) = v(i,j,k)/d(i,j,k)
            w(i,j,k) = w(i,j,k)/d(i,j,k)
c
c           If te is really total energy (and it is unless imethod=2),
c             then just set this value
c
            if (imethod .ne. 2 .and. idual .eq. 1) then
               te(i,j,k) = 0.5_RKIND*(u(i,j,k)**2 + 
     &              v(i,j,k)**2 + w(i,j,k)**2) +
     &              ge(i,j,k)
            endif

 10         continue
         endif

 100     continue

      enddo

      return
      end

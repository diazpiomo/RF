      program langevin_md
      implicit none
      include 'var.inc'

      double precision cpu0

      call cpu_time(cpu0)
      call fecha
      call readdata
      call posiciones
      call constantes
      call pbc
      call configuracion(0)
      call movie(0,0)
      call letrero
      call informa
      call cabeza
      call promedios(0)
      call muestra(0)
      call lista
      call force
      call mdloop
      call promedios(2)
      call muestra(2)
      call configuracion(1)
      call movie(2,0)
      call tiempo(cpu0)
      call fecha

      stop
      end
!======================================================================
      subroutine fecha
!======================================================================
      implicit none

      integer anno, mes, dia, hora, minuto, segundo, milisegundo

      call calendar (anno,mes,dia,hora,minuto,segundo,milisegundo)
      write(*,002)dia,mes,anno,hora,minuto,segundo,milisegundo

002   format(/,t20,'fecha: ',i2,'-',i2.2,'-',i4,/
     :        ,t20,'hora : ',i2,':',i2.2,':',i2,'.',i3,/)

      return
      end
!=====================================================================
      subroutine calendar(anno,mes,dia,hora,minuto,segundo,milisegundo)
!=====================================================================
      implicit none

      integer      anno, mes, dia, hora, minuto, segundo, milisegundo
      integer      values(8)
      character (len=5)  zone
      character (len=8)  date
      character (len=10) time

      call date_and_time (date,time,zone,values)

      anno        = values(1)
      mes         = values(2)
      dia         = values(3)
      hora        = values(5)
      minuto      = values(6)
      segundo     = values(7)
      milisegundo = values(8)

      return
      end
!=====================================================================
      subroutine tiempo(cpu0)
!=====================================================================
      implicit none
      include 'var.inc'

      double precision  fin, fin_s, fin_m, fin_h, cpu0, cpu1

      call cpu_time(cpu1)

      fin   = cpu1 - cpu0
      fin_s = fin
      fin_m = fin/60.0d0
      fin_h = fin/3600.0d0
      write(*,001)'tiempo de la simulacion:'
      write(*,002)fin_h, fin_m, fin_s
      write(*,003)fin/dble(nstep)

001   format(t10,a)
002   format(t10,'horas   ',t20,f13.3,/,
     :       t10,'minutos ',t20,f13.3,/,
     :       t10,'segundos',t20,f13.3)
003   format(t10,'tiempo por paso:',t26,f10.6)

      return
      end
!=====================================================================
      subroutine readdata
!=====================================================================
      implicit none
      include 'var.inc'

      integer i, j

      open (1,file='run.txt',status='unknown')
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*) lattice
       read(1,*) tipo
       read(1,*) nat
       read(1,*) boxx,boxy,boxz
       read(1,*) keypo
       read(1,*) temp
       read(1,*) nstep
       read(1,*) nequil
       read(1,*) tstep
       read(1,*) nprome
       read(1,*) nprint
       read(1,*) nfort
       read(1,*) nmovie
       read(1,*) rcut,rlist
       read(1,*)
       read(1,*) lrdf
       read(1,*) cte_lj
       read(1,*) !para langevin md
       read(1,*) llang
       read(1,*) gama_lin
       read(1,*) gama_ang
       read(1,*) dd
       read(1,*) lam
       read(1,*) !campo
       read(1,*) lcampo
       read(1,*) cte_h
       read(1,*) hx, hy, hz
       !hz = -1.0*hz
       read(1,*) !Reaction field
       read(1,*) eprf 
       read(1,*) !Lees Edwards 
       read(1,*) le
       read(1,*) delx
       read(1,*) gmadot
       read(1,*) 
      close(1)

      call ucase(keypo)
      call ucase(tipo)

      return
      end
!=====================================================================
      subroutine ucase (string)
!=====================================================================
      implicit none

      character(len=*), intent(inout) :: string
      integer :: i                 ! loop index
      integer :: length            ! length of input string

      length = len(string)

      do i=1,length
       if ( lge(string(i:i),'a') .and. lle(string(i:i),'z') ) then
         string(i:i) = achar ( iachar ( string(i:i) ) - 32 )
       end if
      end do

      return
      end
!=====================================================================
      subroutine posiciones
!=====================================================================
      implicit none
      include 'var.inc'

      integer i

      if(lattice)then

       if(nat>mnat) stop 'nat es mayor que mnat ...'
       vol = boxx*boxy*boxz
       rho = dble(nat)/vol
       if(tipo=='CSS')then 
        call cubica_simple
       elseif(tipo=='RND')then
        call azar
       else
        stop 'tipo no implementado ... '
       endif

       call velocidades

      else

       ukin = 0.0
       urot = 0.0
       open(1,file='fort.1',status='old')
        read(1,*) nat,boxx,boxy,boxz
        if(nat>mnat) stop 'nat es mayor que mnat ...'
        do i=1,nat
         read(1,*)rx(i),ry(i),rz(i)
         read(1,*)vx(i),vy(i),vz(i)
         read(1,*)ex(i),ey(i),ez(i)
         read(1,*)ux(i),uy(i),uz(i)
         ukin  = ukin + vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)
         urot  = urot + ux(i)*ux(i) + uy(i)*uy(i) + uz(i)*uz(i)
        end do
       close(1)
       inercia = 2.0d0/5.0d0
       ukin  = 0.50*ukin
       urot  = 0.50*urot*inercia
       tempi = 2.0d0*(ukin+urot)/dble(5*nat-3)
      end if

      call centro_masa

      return
      end
!=====================================================================
      subroutine velocidades
!=====================================================================
      implicit none
      include 'var.inc'

      integer          i
      double precision mean, xi1, xi2, xisq, xi, dot, osq, norm, o
      double precision rangauss, ran_uniform

      do i=1,nat
       vx(i) = rangauss()
       vy(i) = rangauss()
       vz(i) = rangauss()
      end do

!      call momento

      inercia = 2.0d0/5.0d0
      mean = 2.0d0*temp/inercia

      do i=1,nat
       xisq = 1.0
1000   if (xisq>=1.0 ) then
        xi1  = ran_uniform()*2.0d0 - 1.0d0
        xi2  = ran_uniform()*2.0d0 - 1.0d0
        xisq = xi1*xi1 + xi2*xi2
        go to 1000
       endif

       xi    = sqrt(1.0d0-xisq)
       ux(i) = 2.0d0*xi1*xi
       uy(i) = 2.0d0*xi2*xi
       uz(i) = 1.0d0 - 2.0d0*xisq

       dot   = ux(i)*ex(i) + uy(i)*ey(i) + uz(i)*ez(i)
       ux(i) = ux(i) - dot*ex(i)
       uy(i) = uy(i) - dot*ey(i)
       uz(i) = uz(i) - dot*ez(i)

       osq   = ux(i)*ux(i) + uy(i)*uy(i) + uz(i)*uz(i)
       norm  = dsqrt(osq)
       ux(i) = ux(i)/norm
       uy(i) = uy(i)/norm
       uz(i) = uz(i)/norm

       osq   = -mean*dlog (ran_uniform())
       o     = dsqrt ( osq )
       ux(i) = o*ux(i)
       uy(i) = o*uy(i)
       uz(i) = o*uz(i)
      end do

      ukin = 0.0d0
      urot = 0.0d0
      do i=1,nat
       ukin = ukin + vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)
       urot = urot + ux(i)*ux(i) + uy(i)*uy(i) + uz(i)*uz(i)
      end do
      ukin = 0.50d0*ukin
      urot = 0.50d0*urot*inercia

      tempi = 2.0d0*(ukin+urot)/dble(5*nat-3)
      scale = sqrt(temp/tempi)

      do i=1,nat
       vx(i) = scale*vx(i)
       vy(i) = scale*vy(i)
       vz(i) = scale*vz(i)

       ux(i) = scale*ux(i)
       uy(i) = scale*uy(i)
       uz(i) = scale*uz(i)
      end do

      return
      end
!=====================================================================
      subroutine cubica_simple
!=====================================================================
      implicit none
      include 'var.inc'

      integer          i, j, k, num, nplace
      double precision desx, desy, desz, des, ei
      double precision rangauss, dot_pdto

      num  = idint((dble(nat)**(1.0d0/3.0d0)) + 1.0d0)
      desx = boxx/dble(num)
      desy = boxy/dble(num)
      desz = boxz/dble(num)
      des  = min(desx,desy,desz)

      nplace = 1
      do i=1,num
       do j=1,num
        do k=1,num
         if (nplace <= nat) then
          rx(nplace) = dble(i)*des
          ry(nplace) = dble(j)*des
          rz(nplace) = dble(k)*des
          ex(nplace) = rangauss()
          ey(nplace) = rangauss()
          ez(nplace) = rangauss()
          nplace     = nplace + 1
         endif
        enddo
       enddo
      enddo

      do i=1,nat
       ei = dot_pdto(ex(i),ey(i),ez(i),ex(i),ey(i),ez(i))  
       ex(i) = ex(i)/ei
       ey(i) = ey(i)/ei
       ez(i) = ez(i)/ei
      enddo

      return
      end
!=====================================================================
      subroutine azar
!=====================================================================
      implicit none
      include 'var.inc'

      integer          i, j, insert
      double precision dx,dy,dz, ei
      double precision ran_uniform,rangauss,minima,dot_pdto

      insert = 0

      do i=1,nat
100    rx(i)  = boxx*ran_uniform()
       ry(i)  = boxy*ran_uniform()
       rz(i)  = boxz*ran_uniform()
       insert = insert + 1
       do j=1,i-1
        dx = rx(i) - rx(j)
        dy = ry(i) - ry(j)
        dz = rz(i) - rz(j)

        dx = minima(dx,boxx)
        dy = minima(dy,boxy)
        dz = minima(dz,boxz)

        rij = sqrt(dx*dx + dy*dy + dz*dz)
        if (rij <= 1.0d0) goto 100
       enddo
      !write(11,*)i,insert
      insert = 0
      enddo

      do i=1,nat
       ex(i) = rangauss()
       ey(i) = rangauss()
       ez(i) = rangauss()
       ei    = dot_pdto(ex(i),ey(i),ez(i),ex(i),ey(i),ez(i))  
       ex(i) = ex(i)/ei
       ey(i) = ey(i)/ei
       ez(i) = ez(i)/ei
      enddo

      return
      end
!=====================================================================
      subroutine centro_masa
!=====================================================================
      implicit none
      include 'var.inc'

      integer i
      double precision rx_cm,ry_cm,rz_cm

      rx_cm = 0.0d0
      ry_cm = 0.0d0
      rz_cm = 0.0d0

      do i=1,nat
       rx_cm = rx_cm + rx(i)
       ry_cm = ry_cm + ry(i)
       rz_cm = rz_cm + rz(i)
      enddo

      rx_cm=rx_cm/dble(nat)
      ry_cm=ry_cm/dble(nat)
      rz_cm=rz_cm/dble(nat)

      do i=1,nat
       rx(i) = rx(i) + (boxx/2.0d0 - rx_cm)
       ry(i) = ry(i) + (boxy/2.0d0 - ry_cm)
       rz(i) = rz(i) + (boxz/2.0d0 - rz_cm)
      enddo

      return
      end
!=====================================================================
      subroutine configuracion(switch)
!=====================================================================
      implicit none
      include 'var.inc'

      integer i, switch

      if(switch==1) then

       write(*,*)'guardando configuracion ...'
       open(1,file='fort.1',status='unknown')
        write(1,*)nat,boxx,boxy,boxz
        do i=1,nat
         write(1,*)rx(i),ry(i),rz(i)
         write(1,*)vx(i),vy(i),vz(i)
         write(1,*)ex(i),ey(i),ez(i)
         write(1,*)ux(i),uy(i),uz(i)
        end do
       close(1)

       open(97,file='Foto.ray',status='unknown')
        write(97,*)'2,'
        write(97,*)'0.5,0.5,'
        write(97,*)nat,','
        write(97,*)boxx,',',boxy,',',boxz,','
        do i=1,nat
         write(97,*)rx(i),',',ry(i),',',rz(i),','
         write(97,*)ex(i),',',ey(i),',',ez(i),','
        end do
       close(97)

      else

       open(97,file='Foto.ray',status='unknown')
        write(97,*)'2,'
        write(97,*)'0.5,0.5,'
        write(97,*)nat,','
        write(97,*)boxx,',',boxy,',',boxz,','
        do i=1,nat
         write(97,*)rx(i),',',ry(i),',',rz(i),','
         write(97,*)ex(i),',',ey(i),',',ez(i),','
        end do
       close(97)

      endif

      return
      end
!=====================================================================
      subroutine letrero
!=====================================================================
      implicit none
      include 'var.inc'

      write(*,001)
      write(*,002)'langevin molecular dynamics                '
      write(*,002)'    verlet list                            '        
      write(*,001)

001   format(t10,50('_'))
002   format(t15,a)

      return
      end
!=====================================================================
      subroutine constantes
!=====================================================================
      implicit none
      include 'var.inc'

      double precision ircut, ircut3, ircut6, cte, ircutn

      hboxx = 0.50*boxx
      hboxy = 0.50*boxy
      hboxz = 0.50*boxz

      boxxin=1.0d0/boxx
      boxyin=1.0d0/boxy
      boxzin=1.0d0/boxz

      rct   = min(hboxx,hboxy,hboxz)
      mct   = max(boxx,boxy,boxz)

      if(rcut>rct) then
       rcut  = rct - 1.0
       rlist = rct
       write(*,001)
       write(*,002)
       write(*,003)rcut
       write(*,004)rlist
       write(*,005)
      end if

      c1_lin = dexp(-gama_lin*tstep/2.0d0)
      c2_lin = dsqrt( (1.0d0 -c1_lin**2)*temp)

      c1_ang = dexp(-gama_ang*tstep/2.0d0)
      c2_ang = dsqrt( (1.0d0 -c1_ang**2)*temp)

      if(delx>0.50) delx = 0.50d0
      delvx  = gmadot*boxy
      delx   = delx*boxx

      rcutsq  = rcut*rcut
      rlistsq = rlist*rlist

      cte_rf = -(2.0d0*(eprf-1.0d0))/(2.0d0*eprf+1.0d0)
      cte_rf = cte_rf/(rcutsq*rcut)

      tstep2  = tstep*tstep
      tsteph  = 0.5d0*tstep
      tstepd  = 2.0d0*tstep
  
      tercio = 1.0d0/3.0d0
      pi  = 4.0*datan(1.0d0)

      vol = boxx*boxy*boxz
      rho = dble(nat)/vol

      eps0 = 1.0d0
      eps4 = 4.0d0*eps0
      eps48= 12.0d0*eps4

      inercia = 2.0d0/5.0d0
      pass    = 0.0d0

!      ircut  = 1.0d0/rcut
!      ircut3 = ircut*ircut*ircut
!      ircut6 = ircut3*ircut3

      ircut  = 1.0d0/rcut
      ircutn = ircut**cte_lj
      
      cte   = 0.0d0 !8.0d0*pi*rho*ircut3/9.0d0
      etail = 0.0d0 !cte*(ircut6 - 3.0d0)
      ptail = 0.0d0 !2.0d0*cte*rho*(2.0d0*ircut6 - 3.0d0)

      if(keypo=='TR')then 
       ucut = 0.0d0
       fcut = 0.0d0
      elseif(keypo=='TS')then
       ucut = eps4*ircutn
       fcut = 0.0d0
      elseif(keypo=='SF')then
       ucut = eps4*ircutn
       fcut = ucut*cte_lj*ircut
      endif

001   format(t20,'|-------W  A  R  N  I  N  G-----------|')
002   format(t20,'| radio de corte mayor que box/2 !    |')
003   format(t20,'| usando rcut ==>',f6.3,t58,'|')
004   format(t20,'| usando rlist==>',f6.3,t58,'|')
005   format(t20,'|-------W  A  R  N  I  N  G-----------|',/)

      return
      end
!=====================================================================
      subroutine informa
!=====================================================================
      implicit none
      include 'var.inc'

      if(keypo=='TR')then
       write(*,001)'truncate: [u(rij)]'
      elseif(keypo=='TS')then
       write(*,001) 'truncate & shifted:[u(rij)-u(rc)]'
      elseif(keypo=='SF')then
       write(*,001) 'shifted force: [u(rij)-u(rc)+f(rc)*(r-rc)]'
      endif
      if(lattice)then
       if(tipo=='CSS') then
        write(*,007)'posiciones               :    cs (T=',tempi,')'
       elseif(tipo=='RND') then
        write(*,007)'posiciones               :    random (T=',tempi,')'
       endif
      else
       write(*,007)'posiciones                :    file (T=',tempi,')'
      endif
      write(*,004) 'numero de atomos          : ', nat
      write(*,003) 'caja(x,y,z)               : ', boxx,boxy,boxz
      write(*,005) 'densidad                  : ', rho
      write(*,005) 'temperatura               : ', temp
      write(*,005) 'paso de integracion       : ', tstep
      write(*,005) 'radio de corte            : ', rcut
      write(*,005) 'radio de verlet           : ', rlist
      write(*,005) 'cut-off energia           : ', ucut
      write(*,005) 'cut-off force             : ', fcut
      write(*,004) 'numero de pasos           : ', nstep
      write(*,004) 'promedio cada             : ', nprome
      write(*,004) 'salva las posiciones cada : ', nfort
      write(*,005) 'momento de inercia:esfera : ', inercia
      write(*,*)

001   format(t15,a)
002   format(t15,a,4X,a)
003   format(t15,a,3f10.3)
004   format(t15,a,i7)
005   format(t15,a,f10.5)
006   format(t15,a,3(1X,i5))
007   format(t15,a,f8.4,a)

      return
      end
!=====================================================================
      subroutine cabeza
!=====================================================================
      implicit none
      include 'var.inc'

      write(*,0010)
      write(*,0001)
      write(*,0010)

0010  format(t2,105('_')) 
0001  format(t4,'step',t12,'utot',t23,'ukin',t34,'upot',t45,'tempi',
     :       t56,'press',t67,'p1',t78,'p2',t89,'p4',t100,'udip',
     :       t113,'p1m0',t126,'p1m',t139,'p2m')!nvt

      return
      end
!=====================================================================
      double precision function rangauss()
!=====================================================================
      implicit none

      double precision ran_uniform, v1, v2, rsq

 100  v1  = 2.0d0*ran_uniform() - 1.0d0
      v2  = 2.0d0*ran_uniform() - 1.0d0
      rsq = v1*v1 + v2*v2
      if ( (rsq>=1.0d0) .or. (rsq<= 0.0d0) ) goto 100

      rangauss = v1*sqrt(-2.0d0*log(rsq)/rsq)

      return
      end
!=====================================================================
      double precision function ran_uniform()
!=====================================================================
      implicit none

      integer          idum
      double precision ran2
      data idum /1234567890/
      save idum


      ran_uniform = ran2(idum)

      return
      end
!=====================================================================
      double precision function ran2(idum)
!=====================================================================
      implicit none

      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      double precision am, eps, rnmx
      parameter ( im1 = 2147483563, im2 = 2147483399,
     :     am  = 1.0d0/im1,  imm1 = im1-1, ia1 = 40014,
     :     ia2 = 40692, iq1 = 53668, iq2 = 52774, ir1 = 12211,
     :     ir2 = 3791, ntab = 32, ndiv = 1+imm1/ntab,
     :     eps = 1.2d-7, rnmx = 1.0d0-eps )

      integer idum2, j, k, iv(ntab), iy
      save iv, iy, idum2
      data idum2 /123456789/, iv /ntab*0/, iy /0/

      if (idum .le. 0) then
         idum  = max(-idum,1)
         idum2 = idum
         do j = ntab+8,1,-1
            k   = idum/iq1
            idum= ia1*(idum-k*iq1)-k*ir1
            if (idum .lt. 0) idum = idum+im1
            if (j.le.ntab) iv(j) = idum
         enddo
         iy = iv(1)
      endif
      k   = idum/iq1
      idum= ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum =idum+im1
      k     = idum2/iq2
      idum2 = ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2 = idum2+im2
      j    = 1+iy/ndiv
      iy   = iv(j)-idum2
      iv(j)= idum
      if(iy.lt.1) iy = iy+imm1

      ran2 = min(am*iy,rnmx)

      return
      end
!=====================================================================
      subroutine minima_le(dx,dy,dz)
!=====================================================================
      implicit none
      include 'var.inc'

      double precision cory,dx,dy,dz

      cory = anint(dy*boxyin)
      dx   = dx - cory*delx
      dx   = dx - anint(dx*boxxin)*boxx
      dy   = dy - cory*boxy
      dz   = dz - anint(dz*boxzin)*boxz

      return
      end
!=====================================================================
      double precision function minima(ddd,lado)
!=====================================================================
      implicit none

      double precision ddd, lado, hlado

      hlado = 0.50*lado

      if(ddd > hlado) then
       ddd = ddd - lado
      elseif(ddd < -hlado) then
       ddd = ddd + lado
      endif

      minima = ddd

      return
      end
!=====================================================================
      subroutine force
!=====================================================================
      implicit none
      include 'var.inc'

      integer          i, j, k
      double precision dx, dy, dz, fxi, fyi, fzi, ei
      double precision irijn,irijn1, irijn2, ir3, ulj, gparalela
      double precision dot_pdto,minima
      double precision gsumx,gsumy,gsumz

      gsumx=0.0
      gsumy=0.0
      gsumz=0.0

      do i=1,nat
       fx(i) = 0.0d0
       fy(i) = 0.0d0
       fz(i) = 0.0d0

       gx(i) = 0.0d0
       gy(i) = 0.0d0
       gz(i) = 0.0d0
      end do

      do j=1,6
       press(j) = 0.0d0
      enddo

      upot = 0.0d0
      udip = 0.0d0
      umu = 0.0d0
      urc = 0.0d0

      do k=1,npar

       i = nblist1(k)
       j = nblist2(k)

       dx = rx(i) - rx(j)
       dy = ry(i) - ry(j)
       dz = rz(i) - rz(j)

       if(le)then
        call minima_le(dx,dy,dz)
       else
        dx = minima(dx,boxx)
        dy = minima(dy,boxy)
        dz = minima(dz,boxz)
       endif

       r2  = dx*dx + dy*dy + dz*dz

       if(r2<rcutsq)then
        rij  = dsqrt(r2)
        irijn= 1.0d0/rij**cte_lj
        ir3  = 1.0d0/(rij*r2)

        call ctes_gza(i,j,dx,dy,dz)
        gx(i) =  gx(i) - dd*ir3*gza_x - dd*cte_rf*ex(j)
        gy(i) =  gy(i) - dd*ir3*gza_y - dd*cte_rf*ey(j)
        gz(i) =  gz(i) - dd*ir3*gza_z - dd*cte_rf*ez(j)

        call ctes_gza(j,i,-dx,-dy,-dz)
        gx(j) =  gx(j) - dd*ir3*gza_x - dd*cte_rf*ex(i)
        gy(j) =  gy(j) - dd*ir3*gza_y - dd*cte_rf*ey(i)
        gz(j) =  gz(j) - dd*ir3*gza_z - dd*cte_rf*ez(i)

        call ctes_fza(i,j,dx,dy,dz)

        ulj = eps4*(irijn-lam*irijn*cte_lj*0.5)-ucut + (rcut-rij)*fcut 
       
!       urf = cte_rf*dd*cor_dd + cte_rf*dd
        urf = cte_rf*dd*cor_dd 
        uqq = dd*ir3*cor

        udip= udip + uqq + urf
        upot= upot + ulj + uqq + urf

        umu = umu + uqq
        urc = urc + urf
        duij = eps4*(irijn*cte_lj-lam*irijn*cte_lj*0.5)/rij + fcut

        fxi = duij*dx/rij + dd*ir3*fza_x
        fyi = duij*dy/rij + dd*ir3*fza_y
        fzi = duij*dz/rij + dd*ir3*fza_z

        fx(i) = fx(i) + fxi
        fy(i) = fy(i) + fyi
        fz(i) = fz(i) + fzi

        fx(j) = fx(j) - fxi
        fy(j) = fy(j) - fyi
        fz(j) = fz(j) - fzi

        press(1) = press(1) + fxi*dx
        press(2) = press(2) + fxi*dy
        press(3) = press(3) + fxi*dz
        press(4) = press(4) + fyi*dy
        press(5) = press(5) + fyi*dz
        press(6) = press(6) + fzi*dz
       endif
      enddo
      do i=1,nat
         gsumx=gsumx+gx(i)
         gsumy=gsumy+gy(i)
         gsumz=gsumz+gz(i)
      end do
!     write(78,*)gsumx,gsumy,gsumz

      !if(lcampo) call campo_P2
      if(lcampo) call campo

      do i=1,nat
       ei = dsqrt(dot_pdto(ex(i),ey(i),ez(i),  ex(i),ey(i),ez(i)))
       gparalela = dot_pdto(gx(i),gy(i),gz(i), ex(i),ey(i),ez(i))/ei
       gxp(i) = gx(i) - gparalela*ex(i)/ei
       gyp(i) = gy(i) - gparalela*ey(i)/ei
       gzp(i) = gz(i) - gparalela*ez(i)/ei
      end do

      return
      end
!=====================================================================
      subroutine mdloop
!=====================================================================
      implicit none
      include 'var.inc'

!     integer iequil, paso
      integer iequil
      logical update

      update = .false.
      iequil = 0
      paso   = 0

      do while(iequil<=nequil) 
       iequil = iequil + 1
       call checka(update)
       if(update) call lista
       call velocity(0)
       call velocity(1)
       call force
       call velocity(2)
       call velocity(0)
       call order_parameter
       if(mod(iequil,nfort )==0) call configuracion(1)
       if(mod(iequil,nprint)==0) call screen(iequil)
       if(mod(paso,nmovie)==0  ) call movie(1,paso/nmovie)
      enddo

      do while(paso<=nstep)
       paso = paso + 1
       pass = paso
       call checka(update)
       if(update) call lista
       call velocity(0)
       call velocity(1)
       call force
       call velocity(2)
       call velocity(0)
       call order_parameter
       call screen(paso)
       if(mod(paso,nprome)==0) call promedios(1)
       if(mod(paso,nprome)==0) call muestra(1)
       if(mod(paso,nfort )==0) call configuracion(1)
       if(mod(paso,nmovie)==0) call movie(1,paso/nmovie)
      end do

      return
      end
!=====================================================================
      subroutine velocity(switch)
!=====================================================================
      implicit none
      include 'var.inc'

      integer          i,switch
      double precision lambdap,ei
      double precision rangauss,dot_pdto
      double precision sumex,sumey,sumez

      if(switch==0)then

       if(llang)then
        do i=1,nat
         vx(i) = c1_lin*vx(i) + c2_lin*rangauss()
         vy(i) = c1_lin*vy(i) + c2_lin*rangauss()
         vz(i) = c1_lin*vz(i) + c2_lin*rangauss()

         ux(i) = c1_ang*ux(i) + c2_ang*rangauss()
         uy(i) = c1_ang*uy(i) + c2_ang*rangauss()
         uz(i) = c1_ang*uz(i) + c2_ang*rangauss()
        enddo
       endif

      elseif(switch==1)then

       do i=1,nat
        vx(i) = vx(i) + fx(i)*tsteph
        vy(i) = vy(i) + fy(i)*tsteph
        vz(i) = vz(i) + fz(i)*tsteph

        ux(i) = ux(i) + gxp(i)*tsteph/inercia
        uy(i) = uy(i) + gyp(i)*tsteph/inercia
        uz(i) = uz(i) + gzp(i)*tsteph/inercia
       end do

       call fincham

       do i=1,nat
        ei = dsqrt(dot_pdto(ex(i),ey(i),ez(i),ex(i),ey(i),ez(i)))
        ux(i) = ux(i) - lambda(i)*ex(i)*tsteph/ei
        uy(i) = uy(i) - lambda(i)*ey(i)*tsteph/ei
        uz(i) = uz(i) - lambda(i)*ez(i)*tsteph/ei
       end do

       do i=1,nat
        rx(i) = rx(i) + vx(i)*tstep
        ry(i) = ry(i) + vy(i)*tstep
        rz(i) = rz(i) + vz(i)*tstep

        ex(i) = ex(i) + ux(i)*tstep
        ey(i) = ey(i) + uy(i)*tstep
        ez(i) = ez(i) + uz(i)*tstep
       end do

       call pbc
       
!      do i=1,nat
!      sumex = sumex + ex(i)
!      sumey = sumey + ey(i)
!      sumez = sumez + ez(i)
!      end do
!     write(*,*) sumex,sumey,sumez

      elseif(switch==2)then

       ukin = 0.0d0
       urot = 0.0d0

       do i=1,nat
        vx(i) = vx(i) + fx(i)*tsteph
        vy(i) = vy(i) + fy(i)*tsteph
        vz(i) = vz(i) + fz(i)*tsteph
        ukin  = ukin  + vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)

        ei = dsqrt(dot_pdto(ex(i),ey(i),ez(i),ex(i),ey(i),ez(i)))
        lambdap = dot_pdto(ux(i),uy(i),uz(i), ex(i),ey(i),ez(i))/ei
        lambdap = 2.0d0*lambdap/tstep

        ux(i) = ux(i) + ( gxp(i)/inercia - lambdap*ex(i)/ei )*tsteph
        uy(i) = uy(i) + ( gyp(i)/inercia - lambdap*ey(i)/ei )*tsteph
        uz(i) = uz(i) + ( gzp(i)/inercia - lambdap*ez(i)/ei )*tsteph
        urot  = urot  + ux(i)*ux(i) + uy(i)*uy(i) + uz(i)*uz(i)
       end do

       ukin  = 0.50d0*ukin
       urot  = 0.50d0*urot*inercia
       tempi = 2.0d0*(ukin+urot)/dble(5*nat-3)

       do i=1,nat
        press(1) = press(1) + vx(i)*vx(i)
        press(2) = press(2) + vx(i)*vy(i)
        press(3) = press(3) + vx(i)*vz(i)
        press(4) = press(4) + vy(i)*vy(i)
        press(5) = press(5) + vy(i)*vz(i)
        press(6) = press(6) + vz(i)*vz(i)
       enddo

       do i=1,6
        press(i) = press(i)/vol
       enddo

       if(.not.llang)then
        scale = dsqrt(temp/tempi)
        do i=1,nat
         vx(i) = vx(i)*scale
         vy(i) = vy(i)*scale
         vz(i) = vz(i)*scale

         ux(i) = ux(i)*scale
         uy(i) = uy(i)*scale
         uz(i) = uz(i)*scale
        enddo
       endif
!       call momento
      endif

      return
      end
!=====================================================================
      subroutine momento
!=====================================================================
      implicit none
      include 'var.inc'

      integer i
      double precision momx,momy,momz

      momx = 0.0d0
      momy = 0.0d0
      momz = 0.0d0

      do i=1,nat
       momx = momx + vx(i)
       momy = momy + vy(i)
       momz = momz + vz(i)
      end do

      momx = momx/dble(nat)
      momy = momy/dble(nat)
      momz = momz/dble(nat)

      do i=1,nat
       vx(i) = vx(i) - momx
       vy(i) = vy(i) - momy
       vz(i) = vz(i) - momz
      end do

      return
      end
!=====================================================================
      subroutine pbc
!=====================================================================
      implicit none
      include 'var.inc'

      integer i
      double precision cory

      if(le)then
       do i=1,nat
        cory  = anint(ry(i)*boxyin)
        rx(i) = rx(i) - cory*delx
        rx(i) = rx(i) - anint(rx(i)*boxxin)*boxx
        ry(i) = ry(i) - cory*boxy
        rz(i) = rz(i) - anint(rz(i)*boxzin)*boxz
        vx(i) = vx(i)  - cory*delvx
       enddo
      else
       do i=1,nat
        if(rx(i)<0   ) rx(i) = rx(i) + boxx
        if(rx(i)>boxx) rx(i) = rx(i) - boxx
        if(ry(i)<0   ) ry(i) = ry(i) + boxy
        if(ry(i)>boxy) ry(i) = ry(i) - boxy
        if(rz(i)<0   ) rz(i) = rz(i) + boxz
        if(rz(i)>boxz) rz(i) = rz(i) - boxz
       enddo
      endif

      return
      end
!=====================================================================
      subroutine checka(update)
!=====================================================================
      implicit none
      include 'var.inc'

      double precision  dspm, dspmx, dspmy, dspmz
      double precision  dx, dy, dz, minima
      logical update
      integer i

      dspmx = 0.0d0
      dspmy = 0.0d0
      dspmz = 0.0d0

      do i=1,nat
       dx = rx(i) - rxo(i)
       dy = ry(i) - ryo(i)
       dz = rz(i) - rzo(i)

       if(le)then
        call minima_le(dx,dy,dz)
       else
        dx = minima(dx,boxx)
        dy = minima(dy,boxy)
        dz = minima(dz,boxz)
       endif

       ddx(i) = ddx(i) + dx
       ddy(i) = ddy(i) + dy
       ddz(i) = ddz(i) + dz

       dspmx = max(dabs(ddx(i)),dspmx)
       dspmy = max(dabs(ddy(i)),dspmy)
       dspmz = max(dabs(ddz(i)),dspmz)
      end do

      dspm   = sqrt(dspmx*dspmx  + dspmy*dspmy + dspmz*dspmz)
      update = (dspm>(rlist-rcut))

      return
      end
!=====================================================================
      subroutine screen()
!=====================================================================
      implicit none
      include 'var.inc'

!     integer paso

      ukin = ukin/nat
      upot = upot/nat
      udip = udip/nat
      umu = umu/nat  
      urc = urc/nat
      utot = ukin + upot
      presi=(press(1)+press(4)+press(6))/3.0d0 

      if(mod(paso,nprint)==0) then
       write(*, 001)paso,utot,ukin,upot,tempi,presi,p1,p2,p4,
     : p1m0, p1m,p2m,umu,urc,udip
       write(99,001)paso,utot,ukin,upot,tempi,presi,p1,p2,p4,
     : p1m0, p1m,p2m,umu,urc,udip
      endif

001   format(i7,14(1x,f10.5))

      return
      end
!=====================================================================
      subroutine lista
!=====================================================================
      implicit none
      include 'var.inc'

      integer i, j
      double precision dx, dy, dz, rij2, minima

      npar = 0

      do i=1,nat-1
       do j=i+1,nat
        dx = rx(i) - rx(j)
        dy = ry(i) - ry(j)
        dz = rz(i) - rz(j)

        if(le)then
         call minima_le(dx,dy,dz)
        else
         dx = minima(dx,boxx)
         dy = minima(dy,boxy)
         dz = minima(dz,boxz)
        endif

        rij2 = dx*dx + dy*dy + dz*dz

        if(rij2<rlistsq)then
         npar          = npar + 1
         nblist1(npar) = i
         nblist2(npar) = j
        end if
       enddo
      enddo

      do i=1,nat
       ddx(i) = 0.0d0
       ddy(i) = 0.0d0
       ddz(i) = 0.0d0

       rxo(i) = rx(i)
       ryo(i) = ry(i)
       rzo(i) = rz(i)
      end do
 
      if(npar>mlist)stop 'incrementa mlist ...'

      return
      end
!=====================================================================
      subroutine promedios(switch)
!=====================================================================
      implicit none
      include 'var.inc'

      integer          switch, i, nmpro
      double precision autot,aupot,aukin,apresi,atempi,arho,avol,aboxx
      double precision vutot,vupot,vukin,vpresi,vtempi,vrho,vvol,vboxx
      double precision a2utot,a2upot,a2ukin,a2presi,a2tempi,a2rho,a2vol,
     :                 a2boxx
      double precision kt,cv,c4
      double precision mmu2, ydiel, dielec, ydielm, dielma, dielecm,
     :                 errdielm


      if(switch==0)then

       veces   = 0.0d0
       a4upot  = 0.0d0
       do i=1,2
        a_utot(i)  = 0.0d0
        a_upot(i)  = 0.0d0
        a_ukin(i)  = 0.0d0
        a_presi(i) = 0.0d0
        a_tempi(i) = 0.0d0
        a_rho(i)   = 0.0d0
        a_vol(i)   = 0.0d0
        a_boxx(i)  = 0.0d0
       enddo
        diela=0.0d0
        diel2a=0.0d0
        dielmxa=0.0d0
        dielmya=0.0d0
        dielmza=0.0d0
        dielm2a=0.0d0
        dielmsqta=0.0d0

      elseif(switch==1)then

       veces      = veces + 1.0d0
       a_utot(1)  = a_utot(1) + utot
       a_utot(2)  = a_utot(2) + utot*utot
       a_upot(1)  = a_upot(1) + upot
       a_upot(2)  = a_upot(2) + upot*upot
       a4upot     = a4upot    + upot*upot*upot*upot
       a_ukin(1)  = a_ukin(1) + ukin
       a_ukin(2)  = a_ukin(2) + ukin*ukin
       a_presi(1) = a_presi(1)+ presi
       a_presi(2) = a_presi(2)+ presi*presi
       a_tempi(1) = a_tempi(1)+ tempi
       a_tempi(2) = a_tempi(2)+ tempi*tempi
       
       diela   = diela + diel
       diel2a  = diel2a + diel2
       dielmxa = dielmxa + dielmx
       dielmya = dielmya + dielmy
       dielmza = dielmza + dielmz
       dielm2a = dielm2a + dielm2
       dielmsqta = dielmsqta + dielmsqt

       mmu2 =dd
       nmpro=veces
       ydiel   = 4.0d0*pi*mmu2*nat*nat/(3.0d0*vol*temp)
       dielec  = 1.0d0 + ydiel*(diel2a/nmpro-(diela/nmpro)**2)
!      write(6,*)"<y>**2,<y**2>", (diela/nmpro)**2,diel2a/nmpro,ydiel
       write(67,100)dielec,(diela/nmpro)**2,diel2a/nmpro

       ydielm  = 4.0d0*pi*mmu2/(3.0d0*vol*temp)
       dielma  = dielmxa*dielmxa+dielmya*dielmya+dielmza*dielmza
       dielecm = 1.0d0 + ydielm*(dielm2a/nmpro - dielma/nmpro**2)
       errdielm= (dielma/nmpro**2)/(dielm2a/nmpro)

!      write(6,*)"<y>**2,<y**2>", dielma/nmpro**2,dielm2a/nmpro
       write(68,100)dielecm,dielma/nmpro**2,dielm2a/nmpro
 100   format(3(1x,e15.7))

      elseif(switch==2)then
       autot  = a_utot(1)/veces
       a2utot = a_utot(2)/veces
       aupot  = a_upot(1)/veces
       a2upot = a_upot(2)/veces
       aukin  = a_ukin(1)/veces
       a2ukin = a_ukin(2)/veces
       apresi = a_presi(1)/veces
       a2presi= a_presi(2)/veces
       atempi = a_tempi(1)/veces
       a2tempi= a_tempi(2)/veces
       a4upot = a4upot/veces

       vutot  = a2utot - autot*autot
       vupot  = a2upot - aupot*aupot
       vukin  = a2ukin - aukin*aukin
       vpresi = a2presi- apresi*apresi
       vtempi = a2tempi-atempi*atempi

       cv = vupot/atempi**2
       c4 = 1.0d0 - a4upot/(3.0d0*a2upot**2)
        
       mmu2    = dd
       ydiel   = 4.0d0*pi*mmu2*nat*nat/(3.0d0*vol*temp)
       dielec  = 1.0d0 + ydiel*(diel2a/veces - (diela/veces)**2)
       ydielm  = 4.0d0*pi*mmu2/(3.0d0*vol*temp)
       dielma = dielmxa*dielmxa+dielmya*dielmya+dielmza*dielmza
       dielecm = 1.0d0 + ydielm*(dielm2a/veces - dielma/veces**2)
       errdielm= (dielma/veces**2)/(dielm2a/veces)
       write(6,*)"<y>**2,<y**2>", dielma/veces**2,dielm2a/veces
        write(69,*) dielma/veces**2,dielm2a/veces
        write(*,001)
        call informa

        write(*,002)'--- a v e r a g e s ---'
        write(*,003)'ave utot :', autot,  ' +/-', sqrt(vutot)
        write(*,003)'ave ukin :', aukin,  ' +/-', sqrt(vukin)
        write(*,003)'ave upot :', aupot,  ' +/-', sqrt(vupot)
        write(*,003)'ave tempi:', atempi, ' +/-', sqrt(vtempi)
        write(*,003)'ave presi:', apresi, ' +/-', sqrt(vpresi)
        write(*,005)'cv       :', cv
        write(*,005)'c4       :', c4
        write(*,001)
        write(*,003)'dielc_P1 :', dielec 
        write(*,003)'dielc_M  :', dielecm, ' +/-', sqrt(errdielm)
        write(*,003)'M        :', mmu2*dielmsqta/veces/nat, ' +/-', 
     :       mmu2*sqrt(dielm2a/veces/nat - dielmsqta**2/veces**2/nat**2)
        endif



001   format(t10,50('_'))
002   format(t15,a)
003   format(t15,a,f10.5,a,f8.4)
004   format(t15,2(a,es12.3))
005   format(t15,a,es12.3)

      return
      end
!=====================================================================
      subroutine muestra(switch)
!=====================================================================
      implicit none
      include 'var.inc'

      integer          switch,i,j,bin
      double precision din, gder110, gder220
      double precision dx, dy, dz, minima
      double precision volyz, volxy, volxz
      double precision rhozi, rhoyi, rhoxi, gder, dr, rm
      double precision cte, vshell, r, rk, ri, zi
      save             volyz, volxz, volxy, cte, dr


      if(lrdf)then

      if(switch==0)then

       cont = 0.0d0
       dr   = 0.050d0

       cte = 4.0d0*pi/3.0d0

       volyz = dr*boxy*boxz
       volxz = dr*boxx*boxz
       volxy = dr*boxx*boxy

       nbinr = int(mct/dr) + 1
       do i=1,nbinr
        rhox(i)  = 0.0d0
        rhoy(i)  = 0.0d0
        rhoz(i)  = 0.0d0
       enddo
  
        nbing = int(rct/dr) + 1
        do i=1,nbing
         gdr(i) = 0.0d0
         gdr110(i) = 0.0d0 !new line 
         gdr220(i) = 0.0d0 !new line
        end do

      elseif(switch==1)then

       cont = cont + 1.0
     
       if(le)then
        do i=1,nat
         rx(i) = rx(i) + hboxx
         ry(i) = ry(i) + hboxy
         rz(i) = rz(i) + hboxz
        enddo
       endif

       do i=1,nat
        bin       = int(rx(i)/dr)  + 1
        rhox(bin) = rhox(bin)      + 1.0
        bin       = int(ry(i)/dr)  + 1
        rhoy(bin) = rhoy(bin)      + 1.0
        bin       = int(rz(i)/dr)  + 1
        rhoz(bin) = rhoz(bin)      + 1.0
       enddo

       open(3,file='rho.dat',status='unknown')
        write(3,*)'# bin  rho_z,y,x'
        do i=1,nbinr-1
         r  = dr*dble(i)
         ri = r - 0.50*dr
         rhozi=rhoz(i)/(volxy*cont)
         rhoyi=rhoy(i)/(volxz*cont)
         rhoxi=rhox(i)/(volyz*cont)
         write(3,'(4f10.6)')ri,rhozi,rhoyi,rhoxi
        enddo
       close(3)

        do i=1,nat-1
         do j=i+1,nat
          dx = rx(i) - rx(j)
          dy = ry(i) - ry(j)
          dz = rz(i) - rz(j)
          dx = minima(dx,boxx)
          dy = minima(dy,boxy)
          dz = minima(dz,boxz)
          rij= sqrt(dx*dx + dy*dy + dz*dz)
          bin= int(rij/dr) + 1
          din = (ex(i)*ex(j)+ey(i)*ey(j)+ez(i)*ez(j))
          if(bin<=nbing)then 
           gdr(bin) = gdr(bin)  + 2.0d0
           gdr110(bin) = gdr110(bin) + 2.0*din
           gdr220(bin) = gdr220(bin) + 2.0*(3.0d0*din*din -1.0d0)*0.5d0
          endif
         end do
        end do

        open(2,file='rdf.dat',status='unknown')
         write(2,*)'# bin  rdf_total'
         rm = 0.0
         do i=1,nbing
          r      = dr*dble(i)
          ri     = r - 0.50*dr
          vshell = cte*(r**3-rm**3)
          gder   = gdr(i)/(nat*cont)
          gder   = gder/(vshell*rho)
         gder110 = gdr110(i)/(nat*cont)
         gder110 = gder110/(vshell*rho)
         gder220 = gdr220(i)/(nat*cont)
         gder220 = gder220/(vshell*rho)
          write(2,'(4f12.7)')ri,gder, gder110, gder220
          rm = r
         enddo
        close(2)

       if(le)then
        do i=1,nat
         rx(i) = rx(i) - hboxx
         ry(i) = ry(i) - hboxy
         rz(i) = rz(i) - hboxz
        enddo
       endif

      elseif(switch==2)then

       do i=1,nbinr
        rhox(i) = rhox(i)/(cont*volyz)
        rhoy(i) = rhoy(i)/(cont*volxz)
        rhoz(i) = rhoz(i)/(cont*volxy)
       enddo

       open(3,file='rho.dat',status='unknown')
        write(3,*)'# bin  rho_z,y,x'
        do i=1,nbinr-1
         r  = dr*dble(i)
         ri = r - 0.50*dr
         write(3,'(7f10.6)')ri,rhoz(i),rhoy(i),rhox(i)
        enddo
       close(3)

        open(2,file='rdf.dat',status='unknown')
         write(2,*)'# bin  rdf_total'
         rm = 0.0
         do i=1,nbing
          r      = dr*dble(i)
          ri     = r - 0.50*dr
          vshell = cte*(r**3-rm**3)
          gdr(i)= gdr(i)/(cont*nat)
          gdr(i)= gdr(i)/(vshell*rho)
          gdr110(i)= gdr110(i)/(cont*nat)
          gdr110(i)= gdr110(i)/(vshell*rho)
          gdr220(i)= gdr220(i)/(cont*nat)
          gdr220(i)= gdr220(i)/(vshell*rho)
          write(2,'(4f12.7)')ri,gdr(i), gdr110(i), gdr220(i)
          rm = r
         enddo
        close(2)

      endif

      else

       return

      endif

      return
      end
!=====================================================================
      double precision function dot_pdto( v1x,v1y,v1z, v2x,v2y,v2z)
!=====================================================================
      implicit none

      integer i,j
      double precision v1x,v1y,v1z, v2x,v2y,v2z

      dot_pdto = v1x*v2x + v1y*v2y + v1z*v2z

      return
      end
!=====================================================================
      subroutine ctes_gza(i,j,dx,dy,dz)
!=====================================================================
      implicit none
      include 'var.inc'

      double precision dot_pdto, dx, dy, dz
      double precision ei, ej, eirij, ejrij, eiej
      integer          i,j

      ei = dsqrt(dot_pdto(ex(i),ey(i),ez(i),ex(i),ey(i),ez(i)))
      ej = dsqrt(dot_pdto(ex(j),ey(j),ez(j),ex(j),ey(j),ez(j)))
!      eirij= dot_pdto(dx, dy, dz, ex(i), ey(i), ez(i) )/(rij*ei)
      ejrij= dot_pdto(dx,dy,dz,ex(j),ey(j),ez(j))/(rij*ej)
!      eiej = dot_pdto(ex(i),ey(i),ez(i),ex(j),ey(j),ez(j) )/(ei*ej)

*-----gza dpole --------------------------------------------------------
      gza_x = (ex(j)/ej - 3.0d0*ejrij*dx/rij)/ei
      gza_y = (ey(j)/ej - 3.0d0*ejrij*dy/rij)/ei
      gza_z = (ez(j)/ej - 3.0d0*ejrij*dz/rij)/ei
*-----gza reaction field --------------------------------------------------------
!     gza_x = gza_x + cte_rf*ex(j)
!     gza_y = gza_y + cte_rf*ey(j)
!     gza_z = gza_z + cte_rf*ez(j)

      return
      end
!=====================================================================
      subroutine ctes_fza( i, j, dx, dy, dz )
!=====================================================================
      implicit none
      include 'var.inc'

      double precision dot_pdto, dx, dy, dz
      double precision cte_1, cte_2, cte_3
      double precision ei, ej, eirij, ejrij, eiej
      integer          i,j

      ei = dsqrt(dot_pdto(ex(i),ey(i),ez(i),ex(i),ey(i),ez(i)))
      ej = dsqrt(dot_pdto(ex(j),ey(j),ez(j),ex(j),ey(j),ez(j)))
      eirij = dot_pdto(dx,dy,dz, ex(i),ey(i),ez(i))/(rij*ei)
      ejrij = dot_pdto(dx,dy,dz, ex(j),ey(j),ez(j))/(rij*ej)
      eiej  = dot_pdto(ex(i),ey(i),ez(i),ex(j),ey(j),ez(j))/(ei*ej)

      cor    = eiej - 3.0d0*eirij*ejrij
      cor_dd = eiej

*-----fza qpole --------------------------------------------------------
      fza_x=3.0d0*(cor*dx/r2+eirij*(ex(j)/(ej*rij)-ejrij/(rij*rij)*dx)+ 
     :                     ejrij*(ex(i)/(ei*rij)-eirij/(rij*rij)*dx))
      fza_y=3.0d0*(cor*dy/r2+eirij*(ey(j)/(ej*rij)-ejrij/(rij*rij)*dy)+ 
     :                     ejrij*(ey(i)/(ei*rij)-eirij/(rij*rij)*dy))
      fza_z=3.0d0*(cor*dz/r2+eirij*(ez(j)/(ej*rij)-ejrij/(rij*rij)*dz)+ 
     :                     ejrij*(ez(i)/(ei*rij)-eirij/(rij*rij)*dz))
*------fza lj-----------------------------------------------------------
!      duij = 6.0d0*eps4*(2.0d0/rij**13-cte_lj/rij**7) - fcut*eps48

      return
      end
!=====================================================================
      subroutine fincham
!=====================================================================
      implicit none
      include 'var.inc'

      integer i
      double precision util2,lambda1,lambda2,dot_pdto,disc1

      do i=1,nat
       util2 = dot_pdto(ux(i),uy(i),uz(i),ux(i),uy(i),uz(i))
       disc1 = (4.0d0/tstep**2)**2 - (16.0d0/tstep**2)*util2

       if(disc1<0.d0) then
        write(*,*)'discr1 en finch',i,util2
        stop
       end if

       lambda1 = (4.0d0/tstep**2)/2.0d0  + dsqrt(disc1)/2.0d0
       lambda2 = (4.0d0/tstep**2)/2.0d0  - dsqrt(disc1)/2.0d0

       lambda(i) = min(lambda1,lambda2)

       if(lambda2<0.0) then
        write(*,*)'******************************'
        write(*,*)'negativa',lambda1,lambda2
        write(*,*)'******************************'
        stop
       end if
      end do

      return
      end
!=====================================================================
!
!=====================================================================
      subroutine campo
      implicit none
      include 'var.inc'

      integer i
      double precision u_campo, cte0, dot_pdto

      do i=1,nat
       cte0    = dot_pdto(ex(i),ey(i),ez(i),hx,hy,hz)
       u_campo = -sqrt(dd)*cte0
       upot    = upot  + u_campo
       gx(i)   = gx(i) + sqrt(dd)*hx
       gy(i)   = gy(i) + sqrt(dd)*hy
       gz(i)   = gz(i) + sqrt(dd)*hz
      enddo

      return
      end
!=====================================================================
!
!=====================================================================
      subroutine campo_P2
      implicit none
      include 'var.inc'

      integer i
      double precision u_campo, ei, hi, cte0, dot_pdto
      double precision ee0, ee0x, ee0y, omega

!      ee0 = 7.0
!      omega = 8.380

      do i=1,nat
       ei   = dsqrt(dot_pdto(ex(i),ey(i),ez(i),ex(i),ey(i),ez(i)))
       hi   = dsqrt(dot_pdto(hx,hy,hz,hx,hy,hz))
       cte0 = dot_pdto(ex(i),ey(i),ez(i),hx,hy,hz)
       cte0 = cte0/ei
       u_campo = 1.50d0*cte0*cte0 - 0.5d0
       upot  = upot - cte_h*u_campo
       gx(i) = gx(i) + 3.0d0*cte_h*cte0*hx/hi
       gy(i) = gy(i) + 3.0d0*cte_h*cte0*hy/hi
       gz(i) = gz(i) + 3.0d0*cte_h*cte0*hz/hi
       !gz(i) = gz(i) + ee0*cos(omega*tstep*pass)
       !gy(i) = gy(i) + ee0*sin(omega*tstep*pass)
!       gz(i) = gz(i) + 3.0d0*cte_h*cte0*hz/hi
      enddo

      return
      end
!=====================================================================
      subroutine order_parameter
!=====================================================================
      implicit none
      include 'var.inc'

      integer          i
      double precision qq(3,3),eigval(3),eigvec(3,3)
      double precision eigval1,eigval2,eigval3
      double precision qxx,qxy,qxz,qyy,qyz,qzz
      double precision cte,norm,dot_pdto
      double precision px1,py1,pz1
      double precision px2,py2,pz2

      p1 = 0.0d0
      p4 = 0.0d0
c---<>
      px1 = 0.0d0
      py1 = 0.0d0
      pz1 = 0.0d0
c---<>
      px2 = 0.0d0
      py2 = 0.0d0
      pz2 = 0.0d0
c--tensor qq
      qxx = 0.0d0
      qxy = 0.0d0
      qxz = 0.0d0
      qyy = 0.0d0
      qyz = 0.0d0
      qzz = 0.0d0

      do i=1,nat
       px1 = px1 + ex(i)
       py1 = py1 + ey(i)
       pz1 = pz1 + ez(i)

       px2 = px2 + 0.5d0*(3.0d0*ex(i)*ex(i)- 1.0d0)
       py2 = py2 + 0.5d0*(3.0d0*ey(i)*ey(i)- 1.0d0)
       pz2 = pz2 + 0.5d0*(3.0d0*ez(i)*ez(i)- 1.0d0)

       qxx = qxx + ex(i)*ex(i) - 1.0/3.0
       qxy = qxy + ex(i)*ey(i)
       qxz = qxz + ex(i)*ez(i)
       qyy = qyy + ey(i)*ey(i) - 1.0/3.0
       qyz = qyz + ey(i)*ez(i)
       qzz = qzz + ez(i)*ez(i) - 1.0/3.0
      enddo

      px1 = px1/dble(nat)
      py1 = py1/dble(nat)
      pz1 = pz1/dble(nat)

      px2 = px2/dble(nat)
      py2 = py2/dble(nat)
      pz2 = pz2/dble(nat)

      qq(1,1) = qxx/dble(nat)
      qq(1,2) = qxy/dble(nat)
      qq(1,3) = qxz/dble(nat)
      qq(2,1) = qq(1,2)
      qq(2,2) = qyy/dble(nat)
      qq(2,3) = qyz/dble(nat)
      qq(3,1) = qq(1,3)
      qq(3,2) = qq(2,3)
      qq(3,3) = qzz/dble(nat)

      call jacobi(3,qq,eigval,eigvec)

      eigval1 = 1.50d0*eigval(1)
      eigval2 = 1.50d0*eigval(2)
      eigval3 = 1.50d0*eigval(3)

      p2 = eigval3

      director(1) = eigvec(1,3) !dirx
      director(2) = eigvec(2,3) !diry
      director(3) = eigvec(3,3) !dirz
      norm = sqrt(director(1)**2+director(2)**2+director(3)**2)

      director(1) = director(1)/norm
      director(2) = director(2)/norm
      director(3) = director(3)/norm

      p1m0= director(1)*px1 + director(2)*py1 + director(3)*pz1
      p2m = director(1)*px2 + director(2)*py2 + director(3)*pz2
      p1m = abs(p1m0)

      do i=1,nat
       cte = dot_pdto(ex(i),ey(i),ez(i),
     :                director(1),director(2),director(3))
       p4 =  p4 +(35.0d0*cte**4 - 30.0d0*cte**2 + 3.0d0)/8.0d0
       p1  = p1 + cte
      enddo
      p1 = p1/dble(nat)
      p4 = p4/dble(nat)
      p1=dabs(p1)
!      write(98,'(5f12.8)')director(1),director(2),director(3),p1,p2
      
      dielmx=0.0d0
      dielmy=0.0d0
      dielmz=0.0d0

      do i = 1,nat
      dielmx=dielmx+ex(i)
      dielmy=dielmy+ey(i)
      dielmz=dielmz+ez(i)
      end do
      diel=p1
      diel2=diel*diel
      dielm2=dielmx**2+dielmy**2+dielmz**2
      dielmsqt=dsqrt(dielm2)

      return
      end
!=====================================================================
      subroutine jacobi(np,a,d,v)
!=====================================================================
      implicit none
!      include 'var.inc'

      integer i,j,k,ip,iq,n,np,nrot,maxrot
      real*8 sm,tresh,ss,c,t,theta,tau,h,gg,p
      real*8 a(np,np),d(np),v(np,np),b(np),z(np)
c
c     setup and initialization
c
      maxrot = 100
      nrot = 0
      n = np
      do ip = 1, n
         do iq = 1, n
            v(ip,iq) = 0.0d0
         end do
         v(ip,ip) = 1.0d0
      end do
      do ip = 1, n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.0d0
      end do
c
c     perform the jacobi rotations
c
      do i = 1, maxrot
         sm = 0.0d0
         do ip = 1, n-1
            do iq = ip+1, n
               sm = sm + abs(a(ip,iq))
            end do
         end do
         if (sm .eq. 0.0d0)  goto 10
         if (i .lt. 4) then
            tresh = 0.2d0*sm / n**2
         else
            tresh = 0.0d0
         end if
         do ip = 1, n-1
            do iq = ip+1, n
               gg = 100.0d0 * abs(a(ip,iq))
               if (i.gt.4 .and. abs(d(ip))+gg.eq.abs(d(ip))
     :                    .and. abs(d(iq))+gg.eq.abs(d(iq))) then
                  a(ip,iq) = 0.0d0
               else if (abs(a(ip,iq)) .gt. tresh) then
                  h = d(iq) - d(ip)
                  if (abs(h)+gg .eq. abs(h)) then
                     t = a(ip,iq) / h
                  else
                     theta = 0.5d0*h / a(ip,iq)
                     t = 1.0d0 / (abs(theta)+sqrt(1.0d0+theta**2))
                     if (theta .lt. 0.0d0)  t = -t
                  end if
                  c  = 1.0d0 / sqrt(1.0d0+t**2)
                  ss = t*c
                  tau = ss/(1.0d0+c)
                  h = t * a(ip,iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip,iq) = 0.0d0
                  do j = 1, ip-1
                     gg = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = gg - ss*(h+gg*tau)
                     a(j,iq) = h + ss*(gg-h*tau)
                  end do
                  do j = ip+1, iq-1
                     gg = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = gg - ss*(h+gg*tau)
                     a(j,iq) = h + ss*(gg-h*tau)
                  end do
                  do j = iq+1, n
                     gg = a(ip,j)
                     h = a(iq,j)
                     a(ip,j) = gg - ss*(h+gg*tau)
                     a(iq,j) = h + ss*(gg-h*tau)
                  end do
                  do j = 1, n
                     gg = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = gg - ss*(h+gg*tau)
                     v(j,iq) = h + ss*(gg-h*tau)
                  end do
                  nrot = nrot + 1
               end if
            end do
         end do
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0d0
         end do
      end do
c
c     print warning if not converged
c
   10 continue
      if (nrot .eq. maxrot) then
         write (*,20)
   20    format (/,' JACOBI  --  Matrix Diagonalization not Converged')
      end if
c
c     sort the eigenvalues and vectors
c
      do i = 1, n-1
         k = i
         p = d(i)
         do j = i+1, n
            if (d(j) .lt. p) then
               k = j
               p = d(j)
            end if
         end do
         if (k .ne. i) then
            d(k) = d(i)
            d(i) = p
            do j = 1, n
               p = v(j,i)
               v(j,i) = v(j,k)
               v(j,k) = p
            end do
         end if
      end do

      return
      end
!=====================================================================
      subroutine movie(switch,frame)
!=====================================================================
      implicit none
      include 'var.inc'

      integer i,switch,pref,frame
      double precision rrx, rry, rrz
      character(len=10)  mname, fname

      pref = 4
!      mname(1:pref)='Maya'
!      mname(pref+5:pref+8)='.vtk'

      fname(1:pref)='cnf.'

      if(switch==0)then

!       open(22,file='vmd.xyz',status='unknown',action='write')
!       write(22,*)nat
!       write(22,*)
!       do i=1,nat
!        write(22,'(a,3f10.5)')'C ',rx(i),ry(i),rz(i)
!       enddo

      elseif(switch==1)then

!       write(22,*)nat
!       write(22,*)
!       do i=1,nat
!        write(22,'(a,3f10.5)')'C ',rx(i),ry(i),rz(i)
!       enddo

!       write(mname(pref+1:pref+4),'(i4.4)')frame
       write(fname(pref+1:pref+6),'(i6.6)')frame

!       open(3,file=mname,status='unknown')
!        write(3,*)nat,boxx,boxy,boxz
!        do i=1,nat
!         write(3,'(3f10.5))')rx(i),ry(i),rz(i)
!         !write(3,'(3f10.5))')ex(i),ey(i),ez(i)
!         write(3,'(3f10.5))')vx(i),vy(i),vz(i)
!        enddo
!        close(3)

       open(4,file=fname,status='unknown') !MGA
        write(4,*)nat
        write(4,*)boxx
        write(4,*)boxy
        write(4,*)boxz
        write(4,*)0,0
        if(le) then
         do i=1,nat
          rrx = rx(i) + hboxx
          rry = ry(i) + hboxy
          rrz = rz(i) + hboxz
          write(4,001)rrx,rry,rrz,0.0,0.0,0.0,ex(i),ey(i),ez(i),
     :                0.0,0.0,0.0,i
         enddo
        else
         do i=1,nat
          write(4,001)rx(i),ry(i),rz(i),0.0,0.0,0.0,
     :                ex(i),ey(i),ez(i),0.0,0.0,0.0,i
         enddo
        endif

        close(4)

      elseif(switch==2)then

!       close(22)

      endif
001   format(2(3f9.3,1X,3(f5.2,1X)),1X,i7)

      return
      end

!=====================================================================
!=====================================================================



program StressCalc

implicit none

   integer ox,jp
   real(8) nu
   integer, allocatable :: axyz(:,:,:)
   double precision, allocatable :: Vx(:,:,:) 
   double precision, allocatable :: temp(:,:,:)
   integer, allocatable :: a(:,:), ab(:,:)
   integer sum,len_ab,thr,n1,n2,n3,len_c
   integer i,j,k,x,y,z,rad,num,x0,y0,z0,j1,x1,y1,z1,z2,z3,dx,dy,dz,dirx1,dirx2,diry1,diry2,dirz1,dirz2
   double precision fraction,sr,sz,st,trz,ro,R,p1,p2,p3,fi,Umax,In1,In2,In3,a1,b1,c1,Q,S,i1,i2,i3,smax
   double precision k1,k2,theta,Utot,dU,old,pi,scf,ax,ay,az,p4
   logical flag  
   integer t1,t2,error
   character(20) name
   save temp,dx,dy,dz,x0,y0,z0,x,y,z,x1,y1,z1,dirx1,dirx2,diry1,diry2,dirz1,dirz2
   !$omp threadprivate(temp,dx,dy,dz)
   pi=2.*asin(1.)
   
   open(2, file='stress/parameters.txt')
   read(2,*) ox
   read(2,*) thr
   read(2,*) nu
   close(2)
   jp=ox*ox*ox-1
   
  allocate (a(0:jp,0:3), STAT=error)
!  if (error /= 0) STOP "Not enough memory for A!"
  allocate (axyz(0:ox-1,0:ox-1,0:ox-1), STAT=error)
!  if (error /= 0) STOP "Not enough memory for Axyz!"
  allocate (ab(0:jp,0:2), STAT=error)
!  if (error /= 0) STOP "Not enough memory for Ab!"
  
!  Loading a microstructure
  
   open(2, file='stress/source.txt')
!   read(2,400) axyz
!   close(2)
   x=0; sum=0; len_ab=0
   do k=0,ox-1
    do j=0,ox-1
      do i=0,ox-1
       read(2,400) axyz(i,j,k)
       a(x,0)=i; a(x,1)=j; a(x,2)=k
       a(x,3)=axyz(i,j,k)
       if (a(x,3)==0) then
          sum=sum+1
          ab(len_ab,0)=i; ab(len_ab,1)=j; ab(len_ab,2)=k; len_ab=len_ab+1 
       end if
       x=x+1
      enddo
    enddo
   enddo
  close(2)
  
    allocate (Vx(0:ox-1,0:ox-1,0:ox-1), STAT=error)
!  if (error /= 0) STOP "Not enough memory for Vx!"

   Vx=0. 
   n1=0;n2=0;n3=0
   fraction=dble(sum*100)/dble(ox*ox*ox)
   
!stress calculation
   rad=ox/2
   call omp_set_num_threads(thr)

!$omp parallel shared(rad,len_ab,a,ab,Vx,axyz,thr,n1,n2,n3) 
  allocate (temp(0:ox-1,0:ox-1,0:ox-1), STAT=error)
  if (error /= 0) STOP "Not enough memory for temp!"
!$omp do private(i,x,y,z,x0,y0,z0,x1,y1,j,k1,k2,z1,z2,z3,p1,p2,p3,i1,i2,i3,num,R,ro,theta,sr,st,sz,trz,In1,In2,In3,a1,b1,c1,Q,S,fi,scf,dirx1,dirx2,diry1,diry2,dirz1,dirz2)   
    do i=0, len_ab-1
      x0=ab(i,0); y0=ab(i,1); z0=ab(i,2)
      do j=0,jp                         !periodic boundaries
        if (a(j,3)/=0) then
          x=a(j,0); y=a(j,1); z=a(j,2)
          if (((x0-x)**2+(y0-y)**2+(z0-z)**2)>rad)then
             if (x<0) then 
                 x=ox+x 
             elseif (x>(ox-1)) then 
                 x=x-ox 
             end if
             if (y<0) then 
                 y=ox+y 
             elseif (y>(ox-1)) then 
                 y=y-ox 
             end if
             if (z<0) then 
                 z=ox+z 
             elseif (z>(ox-1)) then 
                 z=z-ox 
             end if
          endif             
                num=(x0-x)**2+(y0-y)**2+(z0-z)**2
                R=dsqrt(dble(num))
                ro=0.5/R
                theta=dacos(dble(z-z0)/R)
                sr=(dcos(theta))**2+ro*ro*ro/(7.0-5.0*nu)*(6.0*(1.0-ro*ro)-(5.0*(5.0-nu)-18.0*ro*ro)*(dcos(theta))**2)
                st=(dsin(theta))**2+ro*ro*ro/(7.0-5.0*nu)*(4.0-5.0*nu+9.0*ro*ro+(5.0*(1.0-2.0*nu)-21.0*ro*ro)*(dcos(theta))**2)
                sz=3.0*ro*ro*ro/(7.0-5.0*nu)*(-2.0+5.0*nu+ro*ro+5.0*(1.0-2.0*nu-ro*ro)*(dcos(theta))**2)
                trz=-dsin(2.*theta)*(0.5+(ro*ro*ro)/(2.*(7.-5.*nu))*(5.*(1.+nu)-12.*ro*ro))
                In1=sr+st+sz                                                        
                In2=sr*st+st*sz+sr*sz-trz*trz
                In3=sr*st*sz-trz*trz*sz
                a1=-In1; b1=In2; c1=-In3
                Q=(a1*a1-3.0*b1)/9.0
                R=(2.0*a1*a1*a1-9.0*a1*b1+27.0*c1)/54.0
                S=Q*Q*Q-R*R
                if (S>0) then
                    fi=dacos(R/dsqrt(Q*Q*Q))/3.0
                    i1=-2.0*dsqrt(Q)*dcos(fi)-a1/3.0
                    i2=-2.0*dsqrt(Q)*dcos(fi+2.0*pi/3.0)-a1/3.0
                    i3=-2.0*dsqrt(Q)*dcos(fi-2.0*pi/3.0)-a1/3.0
                 elseif (S<0) then
                    fi=dacosh(abs(R)/dsqrt(abs(Q*Q*Q)))/3.0
                    k1=1.0
                    i1=-2.0*sign(k1,R)*dsqrt(abs(Q))*dcosh(fi)-a1/3.0
                    i2=sign(k1,R)*dsqrt(abs(Q))*dcosh(fi)-a1/3.0
                    i3=sign(k1,R)*dsqrt(abs(Q))*dcosh(fi)-a1/3.0
                 else
                    if (R>=0) then
                        i1=-2.0*(R**(1/3))-a1/3.0
                        i2=(R**(1/3))-a1/3.0
                        i3=i2
                    else 
                        i1=2.0*(abs(R)**(1/3))-a1/3.0
                        i2=-(abs(R)**(1/3))-a1/3.0
                        i3=i2
                    endif
                 end if
                 temp(x,y,z)=max(i1,i2,i3)  
           endif
         enddo

        do j=0,jp                           !correction
            x=a(j,0); y=a(j,1); z=a(j,2)
            if (temp(x,y,z)>1.2)then
                dirz1=z-1
                if (dirz1<0) dirz1=ox-1
                do while (temp(x,y,dirz1)>1.2)
                    z1=dirz1; dirz1=z1-1
                    if (dirz1<0) dirz1=ox-1
                enddo
                k1=temp(x,y,dirz1)
                dirz2=z+1
                if (dirz2>(ox-1)) dirz2=0
                do while (temp(x,y,dirz2)>1.2)
                    z2=dirz2; dirz2=z2+1
                    if (dirz2>(ox-1)) dirz2=0
                enddo
                k2=temp(x,y,dirz2)
                temp(x,y,z)=(k1+k2)/2            
            endif
        enddo 
      temp(x0,y0,z0)=0. 

!$omp critical         
      Vx=Vx+temp; 
!$omp end critical 

    end do
!$omp end do nowait
    deallocate(temp) 
!$omp end parallel 

  where (axyz==0) Vx=0.
  open(2,file="stress/Vx_unc.txt")    
  do x=0,jp   
      i=a(x,0); j=a(x,1); k=a(x,2);
      write(2,100) a(x,3),Vx(i,j,k) 
  end do
  close(2)
    
! Calibration 

  k1=-44.63841-1.23756*dble(sum)        !new calibration up to 2
  k2=45.63837+1.23756*dble(sum)
  Vx=k1+k2*Vx/dble(sum)

  smax=3.*(9.-5.*nu)/(2.*(7.-5.*nu))
  do n1=0, len_ab-1                     !ellipsoid correction
     x0=ab(n1,0); y0=ab(n1,1); z0=ab(n1,2)           
                dy=1;  n2=1
                diry1=y0-1; if (diry1<0) diry1=ox-1
                do while ((axyz(x0,diry1,z0)==0).and.(dy<ox))
                    dz=1
                    dirz1=z0-1; if (dirz1<0) dirz1=ox-1
                    do while ((axyz(x0,diry1,dirz1)==0).and.(dz<ox))
                        dz=dz+1; z=dirz1 
                        dirz1=z-1; if (dirz1<0) dirz1=ox-1
                    end do
                    dirz2=z0+1; if (dirz2>(ox-1)) dirz2=0
                    do while ((axyz(x0,diry1,dirz2)==0).and.(dz<ox))
                        dz=dz+1; z=dirz2 
                        dirz2=z+1; if (dirz2>(ox-1)) dirz2=0
                    end do   
                    if (dz>n2) n2=dz
                    
                    dy=dy+1; y=diry1 
                    diry1=y-1; if (diry1<0) diry1=ox-1
                end do
                
                diry2=y0+1; if (diry2>(ox-1)) diry2=0
                do while ((axyz(x0,diry2,z0)==0).and.(dy<ox))
                    dz=1
                    dirz1=z0-1; if (dirz1<0) dirz1=ox-1
                    do while ((axyz(x0,diry2,dirz1)==0).and.(dz<ox))
                        dz=dz+1; z=dirz1 
                        dirz1=z-1; if (dirz1<0) dirz1=ox-1
                    end do
                    dirz2=z0+1; if (dirz2>(ox-1)) dirz2=0
                    do while ((axyz(x0,diry2,dirz2)==0).and.(dz<ox))
                        dz=dz+1; z=dirz2 
                        dirz2=z+1; if (dirz2>(ox-1)) dirz2=0
                    end do   
                    if (dz>n2) n2=dz

                    dy=dy+1; y=diry2 
                    diry2=y+1; if (diry2>(ox-1)) diry2=0
                end do
                
        dx=1; dirx1=x0-1 
        if (dirx1<0) dirx1=ox-1
        do while ((axyz(dirx1,y0,z0)==0).and.(dx<ox))
                    dz=1
                    dirz1=z0-1; if (dirz1<0) dirz1=ox-1
                   do while ((axyz(dirx1,y0,dirz1)==0).and.(dz<ox))
                        dz=dz+1; z=dirz1 
                        dirz1=z-1; if (dirz1<0) dirz1=ox-1
                    end do
                    dirz2=z0+1; if (dirz2>(ox-1)) dirz2=0
                    do while ((axyz(dirx1,y0,dirz2)==0).and.(dz<ox))
                        dz=dz+1; z=dirz2 
                        dirz2=z+1; if (dirz2>(ox-1)) dirz2=0
                    end do   
                    if (dz>n2) n2=dz            
            dx=dx+1; x=dirx1 
            dirx1=x-1; if (dirx1<0) dirx1=ox-1
        end do
        
        dirx2=x0+1; if (dirx2>(ox-1)) dirx2=0
        do while ((axyz(dirx2,y0,z0)==0).and.(dx<ox))
                    dz=1
                    dirz1=z0-1; if (dirz1<0) dirz1=ox-1
                    do while ((axyz(dirx2,y0,dirz1)==0).and.(dz<ox))
                        dz=dz+1; z=dirz1 
                        dirz1=z-1; if (dirz1<0) dirz1=ox-1
                    end do
                    dirz2=z0+1; if (dirz2>(ox-1)) dirz2=0
                    do while ((axyz(dirx2,y0,dirz2)==0).and.(dz<ox))
                        dz=dz+1; z=dirz2 
                        dirz2=z+1; if (dirz2>(ox-1)) dirz2=0
                    end do   
                    if (dz>n2) n2=dz            
            dx=dx+1; x=dirx2 
            dirx2=x+1
            if (dirx2>(ox-1)) dirx2=0
        end do

        dirx1=x0-1; dirx2=x0+1 
        diry1=y0-1; diry2=y0+1
        dirz1=z0-1; dirz2=z0+1
        if (dirx1<0) dirx1=ox-1; if (dirx2>(ox-1)) dirx2=0
        if (diry1<0) diry1=ox-1; if (diry2>(ox-1)) diry2=0
        if (dirz1<0) dirz1=ox-1; if (dirz2>(ox-1)) dirz2=0
               
        if ((axyz(dirx1,y0,z0)/=0).and.(axyz(dirx2,y0,z0)/=0).and.(axyz(x0,diry1,z0)/=0).and.(axyz(x0,diry2,z0)/=0).and.(axyz(x0,y0,dirz1)/=0).and.(axyz(x0,y0,dirz2)/=0)) then
         do i=dirx1,dirx2                                               !correction for 1 dot
             do j=diry1,diry2
                 do k=dirz1,dirz2
                     if (Vx(i,j,k)>smax) Vx(i,j,k)=smax
                 enddo
             enddo
         enddo
        else
        k2=-0.87939
        k1=2.72956-0.96978*(dble(dy)/dble(dx))
        scf=k1*(dble(n2)/dble(dy))**k2
        p1=Vx(dirx1,y0,z0)/smax; p2=Vx(dirx2,y0,z0)/smax; p3=Vx(x0,diry1,z0)/smax; p4=Vx(x0,diry2,z0)/smax
        if ((axyz(dirx1,y0,z0)/=0).and.(Vx(dirx1,y0,z0)>1.).and.(scf*p1>Vx(dirx1,y0,z0))) Vx(dirx1,y0,z0)=scf*p1
        if ((axyz(dirx2,y0,z0)/=0).and.(Vx(dirx2,y0,z0)>1.).and.(scf*p2>Vx(dirx2,y0,z0))) Vx(dirx2,y0,z0)=scf*p2
        if ((axyz(x0,diry1,z0)/=0).and.(Vx(x0,diry1,z0)>1.).and.(scf*p3>Vx(x0,diry1,z0))) Vx(x0,diry1,z0)=scf*p3      
        if ((axyz(x0,diry2,z0)/=0).and.(Vx(x0,diry2,z0)>1.).and.(scf*p4>Vx(x0,diry2,z0))) Vx(x0,diry2,z0)=scf*p4  
        
        k1=2.72956-0.96978*(dble(dx)/dble(dy))
        if ((axyz(dirx1,y0,z0)/=0).and.(Vx(dirx1,y0,z0)>1.).and.(scf*p1>Vx(dirx1,y0,z0))) Vx(dirx1,y0,z0)=scf*p1
        if ((axyz(dirx2,y0,z0)/=0).and.(Vx(dirx2,y0,z0)>1.).and.(scf*p2>Vx(dirx2,y0,z0))) Vx(dirx2,y0,z0)=scf*p2
        if ((axyz(x0,diry1,z0)/=0).and.(Vx(x0,diry1,z0)>1.).and.(scf*p3>Vx(x0,diry1,z0))) Vx(x0,diry1,z0)=scf*p3      
        if ((axyz(x0,diry2,z0)/=0).and.(Vx(x0,diry2,z0)>1.).and.(scf*p4>Vx(x0,diry2,z0))) Vx(x0,diry2,z0)=scf*p4   
        endif
  enddo
  where (axyz==0) Vx=0.
  where (Vx<0.) Vx=0.
  open(2,file="stress/scf.txt")    
!  do x=0,jp   
!      i=a(x,0); j=a(x,1); k=a(x,2);
!      write(2,100) a(x,3),Vx(i,j,k) 
!  end do
  do k=0,ox-1
      do j=0,ox-1
          do i=0,ox-1
              write(2,100) axyz(i,j,k),Vx(i,j,k)
          enddo
      enddo
  enddo
  close(2)
  
  deallocate(a); deallocate (axyz); deallocate(ab) 
  deallocate(Vx)
  
100 format(1 i4, 1 e30.20)
200 format(100e30.20)
300 format(100i)
400 format(1 i3)
500 format(1 e30.20)
  
end program StressCalc




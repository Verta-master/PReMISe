program Fracture

implicit none

   integer ox
   double precision nu,E,res
   integer, allocatable :: axyz(:,:,:)
   double precision, allocatable :: Vx1(:,:,:),Vx(:,:,:), Vxz(:,:,:),Vx2(:,:,:)
   double precision, allocatable :: temp(:,:,:), Uip(:,:,:), Ui(:,:,:), dUi(:,:,:) 
   integer, allocatable :: a(:,:), ab(:,:), acr(:,:)
   integer sum0,len_ab,len_acr,ck,n1,n2,def,delta,Rc,numU,num2,max_thr,thr,surf,cur,xmax,ymax,zmax,number,sum1,jp
   integer i,j,k,x,y,z,rad,num,x0,y0,z0,j1,x1,y1,z1,z2,z3,dx,dy,dz,dirx1,dirx2,diry1,diry2,dirz1,dirz2,x2,len
   double precision fraction,sr,sz,st,trz,ro,R,p1,p2,p3,p4,fi,Umax,In1,In2,In3,a1,b1,c1,Q,S,i1,i2,i3,fraction0,exp0,E0,xcM,dM,AcM,kcM,p,y01M,gminM,gmaxM,a1M,t1M,totP,pt,stm,xcF,dF,AcF,kcF,y01F,gminF,gmaxF,a1F,t1F
   double precision k1,k2,theta,Utot,dU,old,pi,Uav,scf,dU2,displ,Up,dUp,dUtot,Ur,Wf,Ucr,Ep0,pinf,step,km,V,dUav,Uapp,smax,Ua,Eloc,E1,Gcr,strain,volume,t1,maxp,cr_den,v_den,FF,MF,FP,MP
   logical flag, flag2,flag3,penny,fl_Vy,elliptic,notch,distrM,distrF  
   character (50) fname,lab
   integer :: clock_rate, clock_max,t2,error,weight
   integer, dimension (3) :: max_pos
   save temp,dx,dy,dz,x0,y0,z0,x,y,z,x1,y1,z1,dirz1,dirz2,dirx1,dirx2
   !$omp threadprivate(temp)
   pi=2.*asin(1.)
   
   open(2, file='fracture/parameters.txt')
   read(2,*) ox
   read(2,*) res
   read(2,*) max_thr
   read(2,*) i
   read(2,*) nu
   read(2,*) E0
   read(2,*) i
   read(2,*) Gcr 
!   read(2,*) FP                 !filler porosity
!   read(2,*) MP                 !matrix porosity
    read(2,*) cr_den
    read(2,*) v_den
   close(2)
  elliptic=.false.
   
   open(2, file='fracture/large_scale.txt')
   read(2,*) weight
   read(2,*) pt
   read(2,*) stm
   read(2,*) gminM
   read(2,*) gmaxM
   read(2,*) gminF
   read(2,*) gmaxF
   read(2,*) i
   if (i==1) then
       distrM=.true.
       read(2,*) xcM
       read(2,*) dM
       read(2,*) AcM
       read(2,*) kcM
   else
       distrM=.false.
       read(2,*) y01M
       read(2,*) a1M
       read(2,*) t1M
   endif
   read(2,*) i
   if (i==1) then
       distrF=.true.
       read(2,*) xcF
       read(2,*) dF
       read(2,*) AcF
       read(2,*) kcF
   else
       distrF=.false.
       read(2,*) y01F
       read(2,*) a1F
       read(2,*) t1F
   endif   
   close(2)
       
   jp=ox*ox*ox-1
  
  allocate (a(0:jp,0:3), STAT=error)
  allocate (axyz(0:ox-1,0:ox-1,0:ox-1), STAT=error)
  allocate (ab(0:jp,0:2), STAT=error)
  allocate (acr(0:jp,0:2), STAT=error)
  allocate (Ui(0:ox-1,0:ox-1,0:ox-1), STAT=error)
  allocate (dUi(0:ox-1,0:ox-1,0:ox-1), STAT=error)
  allocate (Uip(0:ox-1,0:ox-1,0:ox-1), STAT=error)
  allocate (Vxz(0:ox-1,0:ox-1,0:ox-1), STAT=error)
  allocate (Vx1(0:ox-1,0:ox-1,0:ox-1), STAT=error)
  allocate (Vx(0:ox-1,0:ox-1,0:ox-1), STAT=error)
  allocate (Vx2(0:ox-1,0:ox-1,0:ox-1), STAT=error)
  
! load initial data from file
   sum0=0; len_ab=0; x=0; num=0 
   axyz=1
   open(1, file="fracture/source.txt")
   do k=0,ox-1
    do j=0,ox-1
     do i=0,ox-1
       read(1,*) a(x,3),Vxz(i,j,k) 
       if (a(x,3)==3) num=num+1
       axyz(i,j,k)=a(x,3)
       a(x,0)=i; a(x,1)=j; a(x,2)=k 
       if (a(x,3)==0) then
           sum0=sum0+1
           ab(len_ab,0)=i; ab(len_ab,1)=j; ab(len_ab,2)=k; len_ab=len_ab+1 
       end if
       x=x+1
     end do
    enddo
   enddo
   close(1)
   TotP=1.-v_den/cr_den*(1.-weight/100.)            
   fraction=dble(sum0*100)/dble(ox*ox*ox)           !large pore fraction
!   FF=dble(num*100)/dble(ox*ox*ox)                  !filler fraction
!   MF=(100.-FF-fraction)/100.                                !matrix fraction
!   TotP=FP*FF/100.+MP*MF+fraction/100.                        !total porosity trough the fractions
!   open(1, file='fracture/fractions.txt')
!   write(1,*) 'Total porosity, %: ', TotP*100
!   write(1,*) 'Large pore fraction, %: ', fraction
!   write(1,*) 'Filler fraction, %: ', FF
!   write(1,*) 'Matrix fraction, %: ', MF*100
!   close(1)

   penny=.true. 
   
    Vx1=0.; Vx=0.
    k1=-44.63841-1.23756*dble(sum0)        
    k2=45.63837+1.23756*dble(sum0)
    Vx1=k1+k2*Vxz/dble(sum0)                                               !calibration
    where (axyz==0) Vx1=0.
    where (Vx1<0.) Vx1=0.
    Vx=Vx1; Vx2=Vx1
    if (elliptic) then                                                      !ellipsoid correction for z-tension
      smax=3.*(9.-5.*nu)/(2.*(7.-5.*nu))
      do n1=0, len_ab-1                                                           
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
                     if (Vx2(i,j,k)>smax) Vx2(i,j,k)=smax
                     if (Vx(i,j,k)>smax) Vx(i,j,k)=smax
                 enddo
             enddo
         enddo
        else
        k2=-0.87939
        k1=2.72956-0.96978*(dble(dy)/dble(dx))
        scf=k1*(dble(n2)/dble(dy))**k2
        if (penny) then
            p1=1.; p2=1.; p3=1.; p4=1.              !for regular ellipsoids
        else
            p1=Vx(dirx1,y0,z0)/smax; p2=Vx(dirx2,y0,z0)/smax; p3=Vx(x0,diry1,z0)/smax; p4=Vx(x0,diry2,z0)/smax
        endif
        if ((axyz(dirx1,y0,z0)/=0).and.(Vx(dirx1,y0,z0)>1.).and.(scf*p1>Vx(dirx1,y0,z0))) Vx(dirx1,y0,z0)=scf*p1
        if ((axyz(dirx2,y0,z0)/=0).and.(Vx(dirx2,y0,z0)>1.).and.(scf*p2>Vx(dirx2,y0,z0))) Vx(dirx2,y0,z0)=scf*p2
        if ((axyz(x0,diry1,z0)/=0).and.(Vx(x0,diry1,z0)>1.).and.(scf*p3>Vx(x0,diry1,z0))) Vx(x0,diry1,z0)=scf*p3
        if ((axyz(x0,diry2,z0)/=0).and.(Vx(x0,diry2,z0)>1.).and.(scf*p4>Vx(x0,diry2,z0))) Vx(x0,diry2,z0)=scf*p4
        
        k1=2.72956-0.96978*(dble(dx)/dble(dy))
        scf=k1*(dble(n2)/dble(dx))**k2
        if ((axyz(dirx1,y0,z0)/=0).and.(Vx(dirx1,y0,z0)>1.).and.(scf*p1>Vx(dirx1,y0,z0))) Vx(dirx1,y0,z0)=scf*p1
        if ((axyz(dirx2,y0,z0)/=0).and.(Vx(dirx2,y0,z0)>1.).and.(scf*p2>Vx(dirx2,y0,z0))) Vx(dirx2,y0,z0)=scf*p2
        if ((axyz(x0,diry1,z0)/=0).and.(Vx(x0,diry1,z0)>1.).and.(scf*p3>Vx(x0,diry1,z0))) Vx(x0,diry1,z0)=scf*p3
        if ((axyz(x0,diry2,z0)/=0).and.(Vx(x0,diry2,z0)>1.).and.(scf*p4>Vx(x0,diry2,z0))) Vx(x0,diry2,z0)=scf*p4
      endif
    enddo    
    where (axyz==0) Vx=0.
    where (Vx<0.) Vx=0.
    endif

    if (weight==0) pt=1.6                                       !virgin graphite calibration
    Ep0=E0*(1.-TotP)**4.12*pt*stm
    
    strain=1E-5;                                                !strain initial value 
    step=1E-5; pinf=strain*Ep0    
    volume=(dble(ox)*res)**3
    Eloc=pinf*pinf*sum(Vx*Vx)/(strain*strain*Ep0)*volume     !homogenisation    
    Uip=pinf*pinf*Vx*Vx/(2.*Eloc)
          
   rad=ox/2; 
   flag=.false.
   dU=0.; delta=0; dUp=0.; def=0
   open(2, file='fracture/fracture.txt')
   write(2,*) 'Strain Stress Ebulk Porosity'
   write(2,300) strain,pinf,Ep0,TotP*100.
   strain=strain+step; pinf=strain*Ep0; maxp=pinf

    call RANDOM_SEED ()
    old=TotP; fraction0=TotP; fraction=TotP
        
    ck=1
  do while ((1.-pinf/maxp)*100.<0.2)            ! % of damage
        len_acr=0
        Utot=0.; dU2=0.; Ui=0.; dUi=0.
        Vx1=0.
        num=sum0; Umax=0 
        Ui=Vx*Vx*pinf*pinf/(2.*Eloc); where (axyz==0) Ui=0.
        Utot=sum(Ui)
        dUi=(Ui-Uip)*res; Uip=Ui
        number=0
        do x2 = 0,jp
          if ((a(x2,3)==1).or.(a(x2,3)==3)) then
              x0=a(x2,0); y0=a(x2,1); z0=a(x2,2) 
              call RANDOM_NUMBER(p); 
              if (a(x2,3)==1) then
                  
                !Strain distribution in matrix                
                if (distrM) then
                   if (p>=AcM) p=AcM-0.001; Gcr=1./kcM*(-log(1.-p/AcM))**(1./dM)+xcM                                !Weibull distribution  
                else
                   Gcr=t1M*log((p-y01M)/a1M)                                                                        !ExpGro1 distribution
                endif
                if (Gcr<gminM) Gcr=gminM; if (Gcr>gmaxM) Gcr=gmaxM   
                                    
              else                  
                !Strain distribution in filler                
                if (distrF) then
                     if (p>=AcF) p=AcF-0.001; Gcr=1./kcF*(-log(1.-p/AcF))**(1./dF)+xcF                              !Weibull distribution
                else                    
                     Gcr=t1F*log((p-y01F)/a1F)                                                                      !ExpGro1 distribution                  
                endif
                if (Gcr<gminF) Gcr=gminF; if (Gcr>gmaxF) Gcr=gmaxF  
                  
              endif
                
                if (pinf*Vx(x0,y0,z0)/Ep0>=Gcr) then                        !critical strain criterion
                    max_pos=maxloc(Vx2)
                    x0=max_pos(1)-1; y0=max_pos(2)-1; z0=max_pos(3)-1
                    Vx2(x0,y0,z0)=0.
                    axyz(x0,y0,z0)=2
                    ab(len_ab,0)=x0; ab(len_ab,1)=y0; ab(len_ab,2)=z0
                    acr(len_acr,0)=x0; acr(len_acr,1)=y0; acr(len_acr,2)=z0
                    len_ab=len_ab+1; sum0=sum0+1; len_acr=len_acr+1; def=def+1
                    dU=dU+Ui(x0,y0,z0); dU2=dU2+Ui(x0,y0,z0); Ui(x0,y0,z0)=0
                    flag=.true.
            endif                   
          endif
        enddo 
        
        do x2=0,jp
            x0=a(x2,0); y0=a(x2,1); z0=a(x2,2)
            if (a(x2,3)/=2) a(x2,3)=axyz(x0,y0,z0)
            if (a(x2,3)==2) axyz(x0,y0,z0)=2
        enddo
        
       where (axyz==2) axyz=0  
        
       num=len_acr

       fraction=fraction+dble(num)/dble(ox*ox*ox)
        
       if (num>0) then
    
    Ep0=E0*(1.-fraction)**4.12*pt*stm

    write(2,300) strain,pinf,Ep0,fraction*100.          !write data
    
    Eloc=pinf*pinf*sum(Vx*Vx)/(strain*strain*Ep0)*volume      !homogenisation
            
    pinf=strain*Ep0; if (pinf>maxp) maxp=pinf
    strain=strain+step            
    
            if (num<max_thr) then
                thr=num
            else
                thr=max_thr
            endif
            call omp_set_num_threads(thr)
            
            !$omp parallel shared(rad,len_ab,a,ab,num,axyz,thr,acr,Vx1)
            allocate (temp(0:ox-1,0:ox-1,0:ox-1), STAT=error)
            if (error /= 0) STOP "Not enough memory for temp!"
           !$omp do private(x2,i,j,k,x0,y0,z0,x1,y1,z1,x,y,z,num2,R,ro,theta,sr,st,sz,trz,In1,In2,In3,a1,b1,c1,Q,S,fi,i1,i2,i3,k1,dirz1,dirz2,k2)
            
            do x2 = 0, num-1 
                i=acr(x2,0); j=acr(x2,1); k=acr(x2,2)
!                call Stress(i,j,k) z-tension
               
                    x0=i; y0=j; z0=k
                    do x1 = x0-rad, x0+rad-1               !periodic boundaries
                        do y1 = y0-rad, y0+rad-1 
                            do z1 = z0-rad, z0+rad-1 
                                if (x1<0) then 
                                    x=ox+x1 
                                elseif (x1>(ox-1)) then 
                                    x=x1-ox 
                                else 
                                    x=x1
                                end if
                                if (y1<0) then 
                                    y=ox+y1 
                                elseif (y1>(ox-1)) then 
                                    y=y1-ox 
                                else 
                                    y=y1
                                end if
                                if (z1<0) then 
                                    z=ox+z1 
                                elseif (z1>(ox-1)) then 
                                    z=z1-ox 
                                else 
                                    z=z1
                                end if
            if (axyz(x,y,z)/=0) then 
                num2=(x0-x)**2+(y0-y)**2+(z0-z)**2
                R=dsqrt(dble(num2))
                ro=0.5/R
                theta=dacos(dble(z-z0)/R)                                                                                  !uniaxial                
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
               p1=max(i1,i2,i3)
                temp(x,y,z)=p1 
              endif
            end do
            enddo
            enddo

            do j=0,jp                           !field correction
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
                temp(x,y,z)=(k1+k2)/2.     
              endif
            enddo 
             temp(x0,y0,z0)=0. 
           
        !$omp critical
            Vx1=Vx1+temp                      
        !$omp end critical
            enddo                   !end Stress z-tension
        !$omp end do nowait
        deallocate(temp)
        !$omp end parallel 
       where (axyz==0) Vx1=0.

       Vxz=Vxz+Vx1                  !save uncalibrated data with previous field

        k1=-44.63841-1.23756*dble(sum0)        
        k2=45.63837+1.23756*dble(sum0)
        Vx1=k1+k2*Vxz/dble(sum0)                                                !calibration
        where (axyz==0) Vx1=0.
        where (Vx1<0.) Vx1=0.
        Vx=Vx1; Vx2=Vx1
        
    if (elliptic) then
        do n1=0, len_ab-1                                                           !ellipsoid correction for z-tension
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
                     if (Vx2(i,j,k)>smax) Vx2(i,j,k)=smax
                     if (Vx(i,j,k)>smax) Vx(i,j,k)=smax
                 enddo
             enddo
         enddo
        else
        k2=-0.87939
        k1=2.72956-0.96978*(dble(dy)/dble(dx))
        scf=k1*(dble(n2)/dble(dy))**k2
        if (penny) then
            p1=1.; p2=1.; p3=1.; p4=1.              !for regular ellipsoids
        else
            p1=Vx(dirx1,y0,z0)/smax; p2=Vx(dirx2,y0,z0)/smax; p3=Vx(x0,diry1,z0)/smax; p4=Vx(x0,diry2,z0)/smax
        endif
        if ((axyz(dirx1,y0,z0)/=0).and.(Vx(dirx1,y0,z0)>1.).and.(scf*p1>Vx(dirx1,y0,z0))) Vx(dirx1,y0,z0)=scf*p1
        if ((axyz(dirx2,y0,z0)/=0).and.(Vx(dirx2,y0,z0)>1.).and.(scf*p2>Vx(dirx2,y0,z0))) Vx(dirx2,y0,z0)=scf*p2
        if ((axyz(x0,diry1,z0)/=0).and.(Vx(x0,diry1,z0)>1.).and.(scf*p3>Vx(x0,diry1,z0))) Vx(x0,diry1,z0)=scf*p3
        if ((axyz(x0,diry2,z0)/=0).and.(Vx(x0,diry2,z0)>1.).and.(scf*p4>Vx(x0,diry2,z0))) Vx(x0,diry2,z0)=scf*p4
        
        k1=2.72956-0.96978*(dble(dx)/dble(dy))
        scf=k1*(dble(n2)/dble(dx))**k2
        if ((axyz(dirx1,y0,z0)/=0).and.(Vx(dirx1,y0,z0)>1.).and.(scf*p1>Vx(dirx1,y0,z0))) Vx(dirx1,y0,z0)=scf*p1
        if ((axyz(dirx2,y0,z0)/=0).and.(Vx(dirx2,y0,z0)>1.).and.(scf*p2>Vx(dirx2,y0,z0))) Vx(dirx2,y0,z0)=scf*p2
        if ((axyz(x0,diry1,z0)/=0).and.(Vx(x0,diry1,z0)>1.).and.(scf*p3>Vx(x0,diry1,z0))) Vx(x0,diry1,z0)=scf*p3
        if ((axyz(x0,diry2,z0)/=0).and.(Vx(x0,diry2,z0)>1.).and.(scf*p4>Vx(x0,diry2,z0))) Vx(x0,diry2,z0)=scf*p4
        endif
    enddo  
    endif

        where (axyz==0) Vx=0.
        where (Vx<0.) Vx=0.
            
       endif                 !end if num>0; end stress
                      
        if (.not.(flag)) then
        strain=strain+step; pinf=strain*Ep0
     endif
      ck=ck+1
      flag=.false.
    enddo
    close(2)
    
    open(8,file='fracture/final_stress.txt')
    do x=0,jp   
        i=a(x,0); j=a(x,1); k=a(x,2);
        write(8,100) a(x,3),Vx(i,j,k)
    end do
    close(8)
    
    deallocate(a); deallocate(axyz); deallocate(ab); deallocate(acr)
    deallocate(Ui); deallocate(dUi); deallocate(Uip)
    deallocate(Vxz); deallocate(Vx1); deallocate(Vx); deallocate(Vx2)
100 format(1 i4, 2 e30.20)
200 format(100 e30.20)
300 format(4 e30.20)
400 format(2 i12, 3 e15.3)   
500 format(1 i1)
600 format(1 i1, 1e30.20)
    
end program Fracture   
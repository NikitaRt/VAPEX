!----------------------------------------------------------------!
!  Graphics output (TECPLOT, MESHTV, PARAVIEW, VisIt etc)        !
!----------------------------------------------------------------!
!  $Id: MOD_Graphics2D.f90 5 2014-01-13 11:44:50Z Sergey $
!----------------------------------------------------------------!
MODULE GRAPHICS_TEC2D
CONTAINS
  SUBROUTINE WriteTECPLOT(NPict, Time, &
       Alpha, Al3, Al3dr, Al3j, U1, V1, U2, V2, &
       P, T1, T2, T3, gam, Y_H2, &
       Ro_a, gamH2, &
       X, Y, N, M, N9, M9,&
       MsgLevel, FileMsgLevel, MUNIT)
    IMPLICIT NONE
    INTEGER:: NPict ! Number of picture
    INTEGER:: N,M,N9,M9
    INTEGER::i,j
    REAL :: UC1, VC1, UC2, VC2
    CHARACTER(len=60) filename
    CHARACTER(len=4) suff
    REAL(kind=8):: U1(N,M9), V1(N9,M), Alpha(N9,M9), X(*), Y(*)
    REAL(kind=8):: U2(N,M9), V2(N9,M), P(N9,M9), Y_H2(N9,M9)
    REAL(kind=8):: T1(N9,M9),T2(N9,M9),T3(N9,M9),gam(N9,M9), Ro_a(N9,M9)
    REAL(kind=8):: Al3(N9,M9), Al3dr(N9,M9), Al3j(N9,M9)
    REAL(kind=8):: gamH2(N9,M9)
    REAL(kind=8):: Time

    REAL(kind=8):: XC,YC
    INTEGER,INTENT(in):: MsgLevel, FileMsgLevel, MUNIT

999 FORMAT('VARIABLES= "X","Y","U1","V1","U2","V2","a_l","Alpha","Phi","Al3",'//&
         '"Al3j","Al3dr","P","T1","T2","T3","Gamma","Y_H2","Ro_H2","Gamma_H2"')
998 FORMAT('ZONE')
997 FORMAT(2X,'T= "',G12.5,'"',2X,'I=',I3,2X,'J=',I3,2X)

    WRITE(suff,'(i4)') NPict

    filename='#tec/0000.dat'
    filename(9:9) = suff(4:4)
    IF(suff(3:3).NE.' ') filename(8:8)=suff(3:3)
    IF(suff(2:2).NE.' ') filename(7:7)=suff(2:2)
    IF(suff(1:1).NE.' ') filename(6:6)=suff(1:1)

    OPEN(50,file=TRIM(filename))

    WRITE(50,999)
    WRITE(50,998)
    WRITE(50,997) Time,N9,M9
    DO  J=1,M9
       DO I=1,N9
          XC=X(i); YC=Y(j)
          IF(i.EQ.1) THEN
             UC1=U1(1,j); UC2=U2(1,j);XC=0.5*(X(1)+X(2))
          ELSE IF(i.EQ.N9) THEN
             UC1=U1(N,j); UC2=U2(N,j);XC=0.5*(X(i)+X(i-1))
          ELSE
             UC1=(U1(i,j)+U1(i-1,j))/2.
             UC2=(U2(i,j)+U2(i-1,j))/2.
          END IF

          IF(j.EQ.1) THEN
             VC1=V1(i,1); VC2=V2(i,1);YC=0.5*(Y(1)+Y(2))
          ELSE IF(j.EQ.M9) THEN
             VC1=V1(i,M); VC2=V2(i,M);YC=0.5*(Y(j)+Y(j-1))
          ELSE
             VC1=(V1(i,j)+V1(i,j-1))/2
             VC2=(V2(i,j)+V2(i,j-1))/2
          END IF

          WRITE (50,'(20(g15.6,1x))') &
               XC,YC, &
               UC1,VC1, &
               UC2,VC2, &
               1.d0-(Alpha(i,j)+Al3(i,j)),&
               Alpha(i,j), Alpha(i,j)/(1.d0-Al3(i,j)),&
               Al3(i,j),Al3j(i,j),Al3dr(i,j),&
               P(i,j), T1(i,j),T2(i,j),T3(i,j),gam(i,j),&
               Y_H2(i,j),Y_H2(i,j)*Alpha(i,j)*Ro_a(i,j),gamH2(i,j)
       END DO
    END DO
    CLOSE(50)

    IF(MsgLevel > 0) THEN
       WRITE(*,100) TRIM(filename),time
    ENDIF
    IF(FileMsgLevel > 0) THEN
       WRITE(MUNIT,100) TRIM(filename),time
    ENDIF
100 FORMAT(80('=')/' Writing graphics file ',a,'   at t = ',g10.3,' [s]'/80('='))

  END SUBROUTINE WriteTECPLOT


  SUBROUTINE WriteTECPLOTDisp(NPict, Time, x, y, u, v, temp, ll)

    IMPLICIT NONE
    INTEGER:: NPict ! Number of picture
    INTEGER:: ll ! number of particles
    CHARACTER(len=60) filename
    CHARACTER(len=4) suff
    REAL(kind=8):: X(*), Y(*),U(*),V(*)
    REAL(kind=8):: temp(*)
    REAL(kind=8):: Time

    INTEGER::i

999 FORMAT('VARIABLES= "X","Y","U","V","Temp"')
998 FORMAT('ZONE')
997 FORMAT(2X,'T= "',G12.5,'"',2X,'I=',I3,2X)

    WRITE(suff,'(i4)') NPict

    filename='#tec/P0000.dat'
    filename(10:10) = suff(4:4)
    IF(suff(3:3).NE.' ') filename(9:9)=suff(3:3)
    IF(suff(2:2).NE.' ') filename(8:8)=suff(2:2)
    IF(suff(1:1).NE.' ') filename(7:7)=suff(1:1)

    OPEN(50,file=filename)

    WRITE(50,999)
    WRITE(50,998)
    WRITE(50,997) Time,ll

    DO i = 1,ll
       WRITE (50,'(18(g15.6,1x))') &
            X(i),Y(i),U(i),V(i),Temp(i)
    END DO
    CLOSE(50)

  END SUBROUTINE WriteTECPLOTDisp


  SUBROUTINE WriteField(filename,f,n,m,x,y)
    IMPLICIT NONE

    CHARACTER(*) filename
    INTEGER:: n,m,i,j
    REAL(kind=8):: f(n,m),x(n),y(m)

    OPEN(50,file=filename)
    WRITE(50,'(i4,1x,i4)') n,m
    DO j = 1,m
       DO i = 1,n
          WRITE(50,'(3(d19.12,5x))') &
               x(i),y(j),f(i,j)
       ENDDO
    ENDDO
  END SUBROUTINE WriteField


  SUBROUTINE ReadField(filename,f,n,m,x,y)
    IMPLICIT NONE

    CHARACTER(*) filename
    INTEGER:: n,m,i,j
    REAL(kind=8),INTENT(out):: f(n,m)
    REAL(kind=8),INTENT(in) :: x(n),y(m)

    INTEGER:: nf, mf
    REAL(kind=8):: xf, yf, ff

    PRINT *,'Reading field from ',filename
    OPEN(50,file=filename)
    READ(50,*) nf, mf
    IF(n/=nf .OR. m/=mf) THEN
       PRINT *,' Different matrix sizes: ',n,nf,m,mf
       STOP
    ENDIF

    DO j = 1,m
       DO i = 1,n
          READ(50,*) &
               xf,yf,ff
          IF(ABS(xf-x(i)) > 1.d-8 .OR. ABS(yf-y(j)) > 1.d-8) THEN
             PRINT *,' Mismatching coordinates: ',i,j,x(i),xf,y(j),yf
             STOP
          ENDIF
          f(i,j) = ff
       ENDDO
    ENDDO
  END SUBROUTINE ReadField


  SUBROUTINE  WriteTECPLOTMesh(x,y,n,m)
    IMPLICIT NONE
    INTEGER:: n,m,i,j
    REAL(kind=8):: x(n),y(m)

999 FORMAT('VARIABLES= "X","Y"')
998 FORMAT('ZONE')
997 FORMAT(2X,'T= Mesh',2X,'I=',I3,2X,'J=',I3,2X)

    OPEN(50,file='#tec/mesh.dat')

    WRITE(50,999)
    WRITE(50,998)
    WRITE(50,997) n,m

    DO j = 1,m
       DO i = 1,n
          WRITE (50,'(18(g15.6,1x))') &
               X(i),Y(j)
       END DO
    ENDDO
    CLOSE(50)
  END SUBROUTINE WriteTECPLOTMesh
END MODULE GRAPHICS_TEC2D


MODULE GRAPHICS_PDB2D

CONTAINS
  SUBROUTINE WritePDB(NPict, Time, Al1, Al2, Al3,al3j,al3dr,al3de,&
       U1, V1, U2, V2, &
       ae1,ae2,aYa,&
       P, PH2O, T1, T2, T3, Ts, gam, Y_H2, Ro_a,&
       gamH2, & !r31c, r32c, &
       X, Y, N, M, N9, M9, &
       MsgLevel, FileMsgLevel,MUNIT)

    IMPLICIT NONE
    INTEGER,INTENT(in):: NPict ! Number of picture
    INTEGER,INTENT(in):: N,M,N9,M9
    REAL(kind=8),INTENT(in):: U1(N,M9), V1(N9,M), X(*), Y(*)
    REAL(kind=8),INTENT(in):: U2(N,M9), V2(N9,M), P(N9,M9), Y_H2(N9,M9)
    REAL(kind=8),INTENT(in):: T1(N9,M9),T2(N9,M9),T3(N9,M9),gam(N9,M9),Ro_a(N9,M9)
    REAL(kind=8),INTENT(in):: Al1(N9,M9),Al2(N9,M9), Ts(N9,M9)
    REAL(kind=8),INTENT(in):: Al3(N9,M9),al3j(N9,M9),al3dr(N9,M9),al3de(N9,M9)
    REAL(kind=8),INTENT(in):: gamH2(N9,M9), PH2O(N9,M9) !, r31c(N9,M9), r32c(N9,M9)
    REAL(kind=8),INTENT(in):: Time
    REAL(kind=8),INTENT(in):: ae1(N9,M9),ae2(N9,M9),aYa(N9,M9)
    INTEGER,INTENT(in):: MsgLevel, FileMsgLevel,MUNIT

    !    include 'silo.inc'
    INTEGER:: dbid,namelen,ierr
    INTEGER:: dims(2), ndims, id
    REAL(kind=8):: U1C(n9,m9), V1C(n9,m9)
    REAL(kind=8):: U2C(n9,m9), V2C(n9,m9)
    REAL(kind=8):: DeltaP(n9,m9), Pwall(m9)
    REAL(kind=8):: Superheat(n9,m9), Phi(n9,m9)
    INTEGER::i,j
    CHARACTER(len=60) filename
    CHARACTER(len=4) suff

    CHARACTER(len=80):: vectorvars
    INTEGER:: ndimvv

!!$    write(suff,'(i4)') NPict
!!$
!!$    filename='#pdb/fn0000.pdb'; namelen=LEN(trim(filename))
!!$    filename(11:11) = suff(4:4)
!!$    if(suff(3:3).ne.' ') filename(10:10)=suff(3:3)
!!$    if(suff(2:2).ne.' ') filename(9:9)=suff(2:2)
!!$    if(suff(1:1).ne.' ') filename(8:8)=suff(1:1)
!!$
!!$    if(MsgLevel > 0) then
!!$       write(*,100) trim(filename),time
!!$    endif
!!$    if(FileMsgLevel > 0) then
!!$       write(MUNIT,100) trim(filename),time
!!$    endif
!!$100 format(80('=')/' Writing graphics file ',a,' at t = ',g8.3,' [s]'/80('='))
!!$
!!$    ndims = 2
!!$    dims(1) = N9; dims(2) = M9
!!$
!!$    ierr = dbcreate(filename,namelen,DB_CLOBBER,DB_LOCAL, &
!!$         DB_F77NULL, 0, DB_PDB, dbid)
!!$    ierr = dbputqm(dbid,'Mesh',4, &
!!$         'R',1,'Z',1,'Z',1,&
!!$         X, Y, 0, dims, ndims, &
!!$         DB_DOUBLE, DB_COLLINEAR,&
!!$         DB_F77NULL, id)
!!$
!!$    U1C = 0;V1C = 0
!!$    U1C(2:n,:) = 0.5*(U1(2:n,:)+U1(1:n-1,:))
!!$    V1C(:,2:m) = 0.5*(V1(:,2:m)+V1(:,1:m-1))
!!$    U2C = 0;V2C = 0
!!$    U2C(2:n,:) = 0.5*(U2(2:n,:)+U2(1:n-1,:))
!!$    V2C(:,2:m) = 0.5*(V2(:,2:m)+V2(:,1:m-1))
!!$    Pwall(:) = P(n,:)
!!$    do i = 1,n9
!!$       DeltaP(i,:) = P(i,:)-Pwall(:)
!!$    enddo
!!$
!!$    Superheat=T1-Ts; 
!!$    Superheat(1,:)=Superheat(2,:);Superheat(n9,:)=Superheat(n,:)
!!$    Superheat(:,1)=Superheat(:,2);Superheat(:,m9)=Superheat(:,m)
!!$
!!$    Phi=Al2/max(Al1+Al2,1.d-10); 
!!$    Phi(1,:)=Phi(2,:);Phi(n9,:)=Phi(n,:)
!!$    Phi(:,1)=Phi(:,2);Phi(:,m9)=Phi(:,m)
!!$
!!$
!!$    !     Write variables
!!$    ierr = dbputqv1(dbid,'U_1',3,&
!!$         'Mesh',4,U1C,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'V_1',3,&
!!$         'Mesh',4,V1C,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'U_2',3,&
!!$         'Mesh',4,U2C,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'V_2',3,&
!!$         'Mesh',4,V2C,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$
!!$    vectorvars = 'Velocity_1 vector {U_1,V_1};Velocity_2 vector {U_2,V_2}'
!!$    ndimvv = LEN(TRIM(vectorvars))
!!$
!!$    ierr = dbwrite(dbid,'_meshtv_defvars',15,&
!!$         vectorvars,ndimvv,1,DB_CHAR)
!!$
!!$    ierr = dbputqv1(dbid,'T1',2,&
!!$         'Mesh',4,T1,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'T2',2,&
!!$         'Mesh',4,T2,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'T3',2,&
!!$         'Mesh',4,T3,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'T1_Ts',5,&
!!$         'Mesh',4,Superheat,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'Al1',3,&
!!$         'Mesh',4,Al1,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'Al2',3,&
!!$         'Mesh',4,Al2,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'Al3',3,&
!!$         'Mesh',4,Al3,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'Al3j',4,&
!!$         'Mesh',4,Al3j,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'Al3dr',5,&
!!$         'Mesh',4,Al3dr,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'Al3de',5,&
!!$         'Mesh',4,Al3de,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'Phi',3,&
!!$         'Mesh',4,Phi,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'P',1,&
!!$         'Mesh',4,P,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'PH2O',4,&
!!$         'Mesh',4,PH2O,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'dP',2,&
!!$         'Mesh',4,DeltaP,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'gamH2',5,&
!!$         'Mesh',4,gamH2,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'gam',3,&
!!$         'Mesh',4,gam,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    DeltaP = Al2*Ro_a*Y_H2
!!$    ierr = dbputqv1(dbid,'Ro_H2',5,&
!!$         'Mesh',4,DeltaP,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'Ro_a',4,&
!!$         'Mesh',4,Ro_a,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'ae1',3,&
!!$         'Mesh',4,ae1,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'ae2',3,&
!!$         'Mesh',4,ae2,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'aYa',3,&
!!$         'Mesh',4,aYa,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'r31c',4,&
!!$         'Mesh',4,r31c,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'r32c',4,&
!!$         'Mesh',4,r32c,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr = dbputqv1(dbid,'Tsv',3,&
!!$         'Mesh',4,Ts,dims,ndims,&
!!$         DB_F77NULL,0,&
!!$         DB_DOUBLE, DB_NODECENT,&
!!$         DB_F77NULL, id)
!!$    ierr=dbclose(dbid)
  END SUBROUTINE WritePDB

  SUBROUTINE WritePDBDisp(NPict, Time, xmesh, ymesh, n, m, x, y, u, v, temp, ind, ll,&
       MsgLevel, FileMsgLevel)

    IMPLICIT NONE
    INTEGER:: NPict ! Number of picture
    INTEGER:: ll ! number of particles
    INTEGER:: n,m
    REAL(kind=8):: X(*), Y(*),U(*),V(*)
    REAL(kind=8):: xmesh(n), ymesh(m)
    REAL(kind=8):: temp(*)
    REAL(kind=8):: Time
    INTEGER:: ind(*)
    INTEGER, INTENT(in):: MsgLevel, FileMsgLevel

    REAL(kind=8):: xbuf(ll), ybuf(ll), buf(ll)

    INTEGER::i, TYPE, numtype, ib
    INTEGER:: dbid,namelen,ierr,meshid,varid
    CHARACTER(len=60) filename
    CHARACTER(len=7) meshname
    CHARACTER(len=6) varname
    CHARACTER(len=4) suff
    INTEGER:: dims(2), ndims

!!$    include 'silo.inc'
!!$
!!$    write(suff,'(i4)') NPict
!!$    filename='#pdb/P0000.pdb'; namelen = LEN(TRIM(filename))
!!$    filename(10:10) = suff(4:4)
!!$    if(suff(3:3).ne.' ') filename(9:9)=suff(3:3)
!!$    if(suff(2:2).ne.' ') filename(8:8)=suff(2:2)
!!$    if(suff(1:1).ne.' ') filename(7:7)=suff(1:1)
!!$
!!$    ndims = 2
!!$    dims(1) = N; dims(2) = M
!!$
!!$    ierr = dbcreate(filename,namelen,DB_CLOBBER,DB_LOCAL, &
!!$         DB_F77NULL, 0, DB_PDB, dbid)
!!$    ierr = dbputqm(dbid,'Mesh',4, &
!!$         'R',1,'Z',1,'Z',1,&
!!$         Xmesh, Ymesh, 0, dims, ndims, &
!!$         DB_DOUBLE, DB_COLLINEAR,&
!!$         DB_F77NULL, meshid)
!!$
!!$    !...Write out the point mesh.
!!$
!!$
!!$    if(ll > 0) then
!!$       ierr = dbputpm (dbid, 'pmesh', 5, 2, x, y, 0, ll, DB_DOUBLE,&
!!$            DB_F77NULL, meshid)
!!$       !...Write out the point variables
!!$
!!$       ierr = dbputpv1 (dbid, 'Temp', 4, 'pmesh', 5, temp, ll, DB_DOUBLE,&
!!$            DB_F77NULL, varid)
!!$
!!$       do type=1,3
!!$          numtype = count(ind(1:ll)==type)
!!$          if(MsgLevel > 0) then
!!$             write(*,*)  'Type : ',type,' NumType = ',NumType
!!$          endif
!!$          if(FileMsgLevel > 0) then
!!$             write(16,*) 'Type : ',type,' NumType = ',NumType
!!$          endif
!!$          if(numtype /= 0) then
!!$             ib = 0
!!$             do i = 1,ll
!!$                if(ind(i) == type) then
!!$                   ib = ib+1
!!$                   xbuf(ib) = x(i); ybuf(ib) = y(i); buf(ib) = temp(i)
!!$                endif
!!$             enddo
!!$             meshname = 'pmesh_0'
!!$             varname =  'Temp_0'
!!$             write(meshname(7:7),'(i1)') type
!!$             write(varname(6:6),'(i1)') type
!!$             ierr = dbputpm (dbid, meshname, 7, 2, xbuf, ybuf, 0, ib, DB_DOUBLE,&
!!$                  DB_F77NULL, meshid)
!!$             ierr = dbputpv1 (dbid, varname, 6, meshname, 7, buf, ib, DB_DOUBLE,&
!!$                  DB_F77NULL, varid)
!!$          endif
!!$       enddo
!!$    endif
!!$    ierr=dbclose(dbid)
  END SUBROUTINE WritePDBDisp
END MODULE GRAPHICS_PDB2D

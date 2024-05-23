module m_Concrete
    use mMaterial
    USE mMechGauss
    use m_utilities

implicit none
! integer, parameter, private :: kreal = kind(0.0d0)

!Added by Shanthanu 2017/3/24

contains

subroutine calBondDmat(istatus,fstatus,matvar,D)
    integer :: istatus(8)
    real (kreal) :: fstatus(44)
    real(kind=kreal) ::  matvar(:)
    real (kreal) :: D(3,3)
    
    D=0.0D0
    if  (istatus(1)==0) call calBondInit (istatus,fstatus,matvar)
    D(1,1)= fstatus(7)  !EMB
    D(2,2)= matvar(5)
    return

end subroutine 

subroutine calBondInit(istatus,fstatus,matvar)
    integer :: istatus(8)
    real(kind=kreal) :: fstatus(44)
    real(kind=kreal) ::  matvar(:)
    
    real(kind=kreal) :: CBAI=2.0
    real(kind=kreal) ::ep_s,ep_t,fp_s,fp_t,hp_s,hp_t,hep_e
    real(kind=kreal) ::en_s,en_t,fn_s,fn_t,hn_s,hn_t,hen_e
    
    
    istatus(1)  =   1   !initflag
    istatus(7)  =   nint(matvar(14))     !icic7
    istatus(8)  =   nint(matvar(15))     !iconcyc
    
    fstatus(1)  =   matvar(6)    !circumference
    fstatus(8)  =   matvar(1)    !E1B
    fstatus(9)  =   matvar(7)    !E1BN
    fstatus(10) =   matvar(8)    !E2P
    fstatus(11) =   matvar(9)    !E2N
    fstatus(13) =   matvar(10)    !TCP
    fstatus(15) =   matvar(11)    !TCN
    fstatus(17) =   matvar(12)    !TYP
    fstatus(18) =   matvar(13)    !TYN
    
    
    istatus(2)  =   10          !IKB
    
    fstatus(7)  =   matvar(1)   !EMB(M,J)=E1B(K,J)
    fstatus(12) =   min(matvar(1),matvar(7))/100.0  !MIN(E1B(K,J),E1BN(K,J))/100.
    
    IF(istatus(7)==1) THEN
        ! TCPS(K,J)= TCP(K,J)/E1B(K,J)
        fstatus(14)=matvar(10)/matvar(1)
        ! TCNS(K,J)= TCN(K,J)/E1B(K,J)
        fstatus(16)=matvar(11)/matvar(1)
        CBAI=2.0
        CALL EFH( fstatus(14),fstatus(13),CBAI,fstatus(14),fstatus(13),  &
                 fstatus(8), EP_S , EP_T , FP_S , FP_T , HP_S , HP_T ,   &
                   HEP_E)
        CALL EFH(fstatus(16),fstatus(15),CBAI,fstatus(16),fstatus(15),   &
                fstatus(8), EN_S , EN_T , FN_S , FN_T , HN_S , HN_T ,    &
                   HEN_E)
! *** [ E点 ] ***
        fstatus(19)= fstatus(14)
        fstatus(20)= fstatus(13)
        fstatus(21)= fstatus(16)
        fstatus(22)= fstatus(15)
! *** [ X点 ] ***
        fstatus(23)= fstatus(14)
        fstatus(24)= fstatus(13)
        fstatus(25)= fstatus(16)
        fstatus(26)= fstatus(15)
! *** [ F点 ] ***
        fstatus(27)= FP_S
        fstatus(28)= FP_T
        fstatus(29)= FN_S
        fstatus(30)= FN_T
! *** [ H点 ] ***
        fstatus(31)= HP_S
        fstatus(32)= HP_T
        fstatus(33)= HN_S
        fstatus(34)= HN_T
! *** [ HE剛性 ] ***
        fstatus(39)= HEP_E
        fstatus(40)= HEN_E
    ENDIF
    

end subroutine 

subroutine calBondUpdate(dsl,istatus,fstatus)
    integer:: istatus(8)
    real(kind=kreal):: fstatus(44)
    
    integer :: initflag,ikb,inbf,itsp,itsn,ndik,icic7,iconcyc
    real (kreal) :: surb,dsl,sdb,tsb,sdbl,tsbl,emb,e1b, &
            e1bn,e2p,e2n,e1b_100,tcp,tcps,tcn,tcns,typ,tyn,esnt,etnt,   &
            fspt,ftpt,xspt,xtpt,ysnt,ytnt,aspt,atpt,bsnt,btnt,cspt,ctpt,&
            dsnt,dtnt,qsnt,qtnt,rspt,rtpt,gspt,gtpt,hsnt,htnt,back1,back2
            
            
    integer inh, inc,inc2,nst,numbf
    
    
        
    call read_fistr_bond (istatus,fstatus,surb,dsl,sdb,tsb,sdbl,tsbl,emb, &
                e1b,e1bn,e2p,e2n,e1b_100,tcp,tcps,tcn,tcns,typ,     &
                tyn,esnt,etnt,fspt,ftpt,xspt,xtpt,ysnt,ytnt,        &
                aspt,atpt,bsnt,btnt,cspt,ctpt,dsnt,dtnt,qsnt,       &
                qtnt,rspt,rtpt,gspt,gtpt,hsnt,htnt,back1,back2,     &
                initflag,ikb,inbf,itsp,itsn,ndik,icic7)
    
    sdb=sdb+dsl
    tsb=tsb+emb*dsl
! C *** CHECK BOND FAILURE BY YIELDING OR CRACKING

     INH=1
	 INC=1
	 INC2=1
    
		IF (ICONCYC.EQ.0) then
            
			CALL TSBOND (SDB,TSB,   &
     			EMB,E1B,E2P,E2N,    &   
     			TCP,TCN,TYP,TYN,    &  
     			IKB,INC,INH,INBF,NUMBF)
		ELSE
			! IF(ICIC7==1) GOTO 7521
			CALL TSBONDCYC(IKB,SDB,TSB,     &
                  EMB,E1B,E1B,E2P,E2N,      &
                  TCP,TCN,TYP,TYN,          &
                  INC,INC2,INH,INBF,NUMBF,  &
                  SDBL,TSBL,                &
                  ASPT,ATPT,BSNT,BTNT,      &
                  CSPT,CTPT,DSNT,DTNT,      &
                  ESNT,ETNT,FSPT,FTPT,      &
                  GSPT,GTPT,HSNT,HTNT,      &
                  QSNT,QTNT,RSPT,RTPT,      &
                  XSPT,XTPT,YSNT,YTNT,      &
                  ITSP,ITSN,                &
                  NDIK,BACK1,BACK2)
			GOTO 7520
7521			CALL TSBONDCYC_hj( IKB,             &
                    DSL,SDB,TSB,SDBL,TSBL,       &
                    EMB,E1B,E1BN,E2P,E2N,E1B_100,   &
                    TCP , TCPS , TCN , TCNS ,       &
                    TYP , TYN,                      &
                  ESNT , ETNT , FSPT , FTPT ,       & ! E
                  XSPT , XTPT , YSNT , YTNT ,   & ! X
                  ASPT , ATPT , BSNT , BTNT ,   & ! F
                  CSPT , CTPT , DSNT , DTNT ,   & ! H
                  QSNT , QTNT , RSPT , RTPT ,   & ! R
                  GSPT, GTPT,     & !HSNT,HTNT,
     			 NST )
		END IF

 7520  CONTINUE
    
    call write_fistr_bond (istatus, fstatus,surb,dsl,sdb,tsb,sdbl,tsbl,emb,     &
                e1b,e1bn,e2p,e2n,e1b_100,tcp,tcps,tcn,tcns,typ,     &
                tyn,esnt,etnt,fspt,ftpt,xspt,xtpt,ysnt,ytnt,        &
                aspt,atpt,bsnt,btnt,cspt,ctpt,dsnt,dtnt,qsnt,       &
                qtnt,rspt,rtpt,gspt,gtpt,hsnt,htnt,back1,back2,     &
                initflag,ikb,inbf,itsp,itsn,ndik,icic7)
    
end subroutine

subroutine read_fistr_bond (istatus, fstatus,surb,dsl,sdb,tsb,sdbl,tsbl,emb,     &
                e1b,e1bn,e2p,e2n,e1b_100,tcp,tcps,tcn,tcns,typ,     &
                tyn,esnt,etnt,fspt,ftpt,xspt,xtpt,ysnt,ytnt,        &
                aspt,atpt,bsnt,btnt,cspt,ctpt,dsnt,dtnt,qsnt,       &
                qtnt,rspt,rtpt,gspt,gtpt,hsnt,htnt,back1,back2,     &
                initflag,ikb,inbf,itsp,itsn,ndik,icic7)
        
real (kreal) :: fstatus(44)
integer :: istatus(8)

integer :: initflag,ikb,inbf,itsp,itsn,ndik,icic7,iconcyc
real (kreal) :: surb,dsl,sdb,tsb,sdbl,tsbl,emb,e1b, &
        e1bn,e2p,e2n,e1b_100,tcp,tcps,tcn,tcns,typ,tyn,esnt,etnt,   &
        fspt,ftpt,xspt,xtpt,ysnt,ytnt,aspt,atpt,bsnt,btnt,cspt,ctpt,&
        dsnt,dtnt,qsnt,qtnt,rspt,rtpt,gspt,gtpt,hsnt,htnt,back1,back2

surb        =   fstatus (   1   )
dsl         =   fstatus (   2   )
sdb         =   fstatus (   3   )
tsb         =   fstatus (   4   )
sdbl        =   fstatus (   5   )
tsbl        =   fstatus (   6   )
emb         =   fstatus (   7   )
e1b         =   fstatus (   8   )
e1bn        =   fstatus (   9   )
e2p         =   fstatus (   10  )
e2n         =   fstatus (   11  )
e1b_100     =   fstatus (   12  )
tcp         =   fstatus (   13  )
tcps        =   fstatus (   14  )
tcn         =   fstatus (   15  )
tcns        =   fstatus (   16  )
typ         =   fstatus (   17  )
tyn         =   fstatus (   18  )
esnt        =   fstatus (   19  )
etnt        =   fstatus (   20  )
fspt        =   fstatus (   21  )
ftpt        =   fstatus (   22  )
xspt        =   fstatus (   23  )
xtpt        =   fstatus (   24  )
ysnt        =   fstatus (   25  )
ytnt        =   fstatus (   26  )
aspt        =   fstatus (   27  )
atpt        =   fstatus (   28  )
bsnt        =   fstatus (   29  )
btnt        =   fstatus (   30  )
cspt        =   fstatus (   31  )
ctpt        =   fstatus (   32  )
dsnt        =   fstatus (   33  )
dtnt        =   fstatus (   34  )
qsnt        =   fstatus (   35  )
qtnt        =   fstatus (   36  )
rspt        =   fstatus (   37  )
rtpt        =   fstatus (   38  )
gspt        =   fstatus (   39  )
gtpt        =   fstatus (   40  )
hsnt        =   fstatus (   41  )
htnt        =   fstatus (   42  )
back1       =   fstatus (   43  )
back2       =   fstatus (   44  )


ikb         =   istatus (   2   )
inbf        =   istatus (   3   )
itsp        =   istatus (   4   )
itsn        =   istatus (   5   )
ndik        =   istatus (   6   )
icic7       =   istatus (   7   )
iconcyc     =   istatus (   8   )


end subroutine

subroutine write_fistr_bond (istatus, fstatus,surb,dsl,sdb,tsb,sdbl,tsbl,emb,     &
                e1b,e1bn,e2p,e2n,e1b_100,tcp,tcps,tcn,tcns,typ,     &
                tyn,esnt,etnt,fspt,ftpt,xspt,xtpt,ysnt,ytnt,        &
                aspt,atpt,bsnt,btnt,cspt,ctpt,dsnt,dtnt,qsnt,       &
                qtnt,rspt,rtpt,gspt,gtpt,hsnt,htnt,back1,back2,     &
                initflag,ikb,inbf,itsp,itsn,ndik,icic7)
        
real (kreal) :: fstatus(44)
integer :: istatus(8)


integer :: initflag,ikb,inbf,itsp,itsn,ndik,icic7,iconcyc
real (kreal) :: surb,dsl,sdb,tsb,sdbl,tsbl,emb,e1b, &
        e1bn,e2p,e2n,e1b_100,tcp,tcps,tcn,tcns,typ,tyn,esnt,etnt,   &
        fspt,ftpt,xspt,xtpt,ysnt,ytnt,aspt,atpt,bsnt,btnt,cspt,ctpt,&
        dsnt,dtnt,qsnt,qtnt,rspt,rtpt,gspt,gtpt,hsnt,htnt,back1,back2


fstatus (   1   )   =   surb     
fstatus (   2   )   =   dsl      
fstatus (   3   )   =   sdb      
fstatus (   4   )   =   tsb      
fstatus (   5   )   =   sdbl     
fstatus (   6   )   =   tsbl     
fstatus (   7   )   =   emb      
fstatus (   8   )   =   e1b      
fstatus (   9   )   =   e1bn     
fstatus (   10  )   =   e2p      
fstatus (   11  )   =   e2n      
fstatus (   12  )   =   e1b_100  
fstatus (   13  )   =   tcp      
fstatus (   14  )   =   tcps     
fstatus (   15  )   =   tcn      
fstatus (   16  )   =   tcns     
fstatus (   17  )   =   typ      
fstatus (   18  )   =   tyn      
fstatus (   19  )   =   esnt     
fstatus (   20  )   =   etnt     
fstatus (   21  )   =   fspt     
fstatus (   22  )   =   ftpt     
fstatus (   23  )   =   xspt     
fstatus (   24  )   =   xtpt     
fstatus (   25  )   =   ysnt     
fstatus (   26  )   =   ytnt     
fstatus (   27  )   =   aspt     
fstatus (   28  )   =   atpt     
fstatus (   29  )   =   bsnt     
fstatus (   30  )   =   btnt     
fstatus (   31  )   =   cspt     
fstatus (   32  )   =   ctpt     
fstatus (   33  )   =   dsnt     
fstatus (   34  )   =   dtnt     
fstatus (   35  )   =   qsnt     
fstatus (   36  )   =   qtnt     
fstatus (   37  )   =   rspt     
fstatus (   38  )   =   rtpt     
fstatus (   39  )   =   gspt     
fstatus (   40  )   =   gtpt     
fstatus (   41  )   =   hsnt     
fstatus (   42  )   =   htnt     
fstatus (   43  )   =   back1    
fstatus (   44  )   =   back2    


istatus (   2   )   =   ikb  
istatus (   3   )   =   inbf 
istatus (   4   )   =   itsp 
istatus (   5   )   =   itsn 
istatus (   6   )   =   ndik 
istatus (   7   )   =   icic7
istatus (   8   )   =   iconcyc

end subroutine

subroutine calEmbedBarStiff(stiff,gausses,ecoord)
USE elementInfo
REAL(kind=kreal), INTENT(IN)    :: ecoord(3, 8)
TYPE(tGaussStatus)              :: gausses(:) 
REAL(kind=kreal), INTENT(OUT)   :: stiff(24,24)
REAL(kind=kreal)                :: gp(2,3)! local gauss point
REAL(kind=kreal)                :: p1(3),p2(3)          ! global gauss point
REAL(kind=kreal)                :: rp1(3),rp2(3)        ! global reinforcement point
REAL(kind=kreal)                :: g1,g2                ! global point location
REAL(kind=kreal)                :: weight=1.0D0             ! global gauss weight
REAL(kind=kreal)                :: dcos(3)
REAL(kind=kreal)                :: length,det,dlength
REAL(kind=kreal)                :: gderiv(8, 3)
REAL(kind=kreal)                :: B(6, 24),DB(6,24)
REAL(kind=kreal)                :: D(6,6)           !Constitutive matrix
REAL(kind=kreal)                :: jacobian(3,3),inverse(3,3)
integer :: i,j,k,ii,jj
integer :: nn=8,ndof=3

g1=(0.5D0)-(1.0D0/(2.0D0*dsqrt(3.0D0)))
g2=(0.5D0)+(1.0D0/(2.0D0*dsqrt(3.0D0)))

do i =1,gausses(1)%istatus(60) !Number of bars
    if (gausses(i)%istatus(61).eq.1) then
        rp1=gausses(i)%fstatus(594:596)
        rp2=gausses(i)%fstatus(597:599)
        length=dsqrt(((rp2(1)-rp1(1))**2)+((rp2(2)-rp1(2))**2)+((rp2(3)-rp1(3))**2))
        dcos=(1/length)*(rp2-rp1)       !direction cosine
        p1=rp1+(dcos*g1*length)
        p2=rp1+(dcos*g2*length)
        call getIsoCoord(ecoord,p1,gp(1,:))
        call getIsoCoord(ecoord,p2,gp(2,:))
        gausses(i)%fstatus(600:602)=dcos
        gausses(i)%fstatus(603:605)=dcos
        gausses(i)%fstatus(606:608)=gp(1,:)
        gausses(i)%fstatus(609:611)=gp(2,:)
        gausses(i)%istatus(61)=2
    else
        dcos=gausses(i)%fstatus(600:602)
        gp(1,:)=gausses(i)%fstatus(606:608)
        gp(2,:)=gausses(i)%fstatus(609:611)
    end if
    call rotate_D_mat(gausses(i)%fstatus(561),gausses(i)%fstatus(593),dcos,D)
    do k = 1,2
        B=0.0D0
        DB=0.0D0
        gderiv=0.0D0
        call getglobalderiv(361, nn, gp(k,:), ecoord, det, gderiv)
        call getJacobian( 361, nn, gp(k,:), ecoord, det, jacobian, inverse )
        dlength=dsqrt(jacobian(1,1)*jacobian(1,1)+jacobian(1,2)*jacobian(1,2)+   &
                        jacobian(1,3)*jacobian(1,3))
        do j = 1,nn
            B(1, 3*j-2) = gderiv(j, 1)
            B(2, 3*j-1) = gderiv(j, 2)
            B(3, 3*j  ) = gderiv(j, 3)
            B(4, 3*j-2) = gderiv(j, 2)
            B(4, 3*j-1) = gderiv(j, 1)
            B(5, 3*j-1) = gderiv(j, 3)
            B(5, 3*j  ) = gderiv(j, 2)
            B(6, 3*j-2) = gderiv(j, 3)
            B(6, 3*j  ) = gderiv(j, 1)
        end do
        DB(1:6, 1:nn*ndof) = MATMUL( D, B(1:6, 1:nn*ndof) )
        FORALL( i=1:nn*ndof, j=1:nn*ndof )
          stiff(i, j) = stiff(i, j)+DOT_PRODUCT( B(:, i), DB(:, j) )*dlength*weight
        END FORALL 
    end do
end do


end subroutine 

subroutine calEmbedBarUpdate(qf,initdisp,totaldisp,gausses,ecoord,iter)
USE elementInfo
REAL(kind=kreal), INTENT(IN)    :: ecoord(3, 8)
TYPE(tGaussStatus)              :: gausses(:) 
REAL(kind=kreal)                :: gp(2,3)    ! local gauss point
REAL(kind=kreal)                :: rp1(3),rp2(3)  ! global reinforcement point
REAL(kind=kreal)                :: p1(3),p2(3)          ! global gauss point
REAL(kind=kreal)                :: weight=1.0D0             ! global gauss weight
REAL(kind=kreal)                :: length,det,strain,istrain,stress1d
REAL(kind=kreal)                :: dcos(3)
REAL(kind=kreal)                :: gderiv(2,8, 3)
REAL(kind=kreal)                :: qf(24)
REAL(kind=kreal)                :: gdispderiv(3, 3),igdispderiv(3, 3)
REAL(kind=kreal)                :: totaldisp(3, 8),initdisp(3, 8)
REAL(kind=kreal)                :: endstrain(2),iendstrain(2),dlength(2)
REAL(kind=kreal)                :: g1,g2                ! global point location
REAL(kind=kreal)                :: B(6, 24),BT(24,6)
REAL(kind=kreal)                :: stress(6)
REAL(kind=kreal)                :: disp(24),stiff(24,24)
REAL(kind=kreal)                :: jacobian(3,3),inverse(3,3)
INTEGER, INTENT(IN) :: iter

integer :: i,j,k,ii,jj
integer :: nn=8,ndof=3

g1=(0.5D0)-(1.0D0/(2.0D0*dsqrt(3.0D0)))
g2=(0.5D0)+(1.0D0/(2.0D0*dsqrt(3.0D0)))
gderiv=0.0D0

do i =1,gausses(1)%istatus(60) !Number of bars
    dcos=gausses(i)%fstatus(600:602)
    gp(1,:)=gausses(i)%fstatus(606:608)
    gp(2,:)=gausses(i)%fstatus(609:611)
    do k = 1,2
        call getglobalderiv(361, nn,gp(k,:), ecoord, det, gderiv(k,:,:))
        gdispderiv(1:ndof, 1:ndof) = MATMUL( totaldisp(1:ndof, 1:nn),   &
                                    gderiv(k,1:nn, 1:ndof) )
        igdispderiv(1:ndof, 1:ndof) = MATMUL( initdisp(1:ndof, 1:nn),   &
                                    gderiv(k,1:nn, 1:ndof) )
        endstrain(k)=dcos(1)*dcos(1)*gdispderiv(1, 1)+     &
                    dcos(2)*dcos(2)*gdispderiv(2, 2)+   &
                    dcos(3)*dcos(3)*gdispderiv(3, 3)+   &
                    dcos(1)*dcos(2)*( gdispderiv(1, 2)+gdispderiv(2, 1) )+   &
                    dcos(2)*dcos(1)*( gdispderiv(2, 3)+gdispderiv(3, 2) )+   &
                    dcos(1)*dcos(2)*( gdispderiv(3, 1)+gdispderiv(1, 3) )
        iendstrain(k)=dcos(1)*dcos(1)*igdispderiv(1, 1)+     &
                    dcos(2)*dcos(2)*igdispderiv(2, 2)+   &
                    dcos(3)*dcos(3)*igdispderiv(3, 3)+   &
                    dcos(1)*dcos(2)*( igdispderiv(1, 2)+igdispderiv(2, 1) )+   &
                    dcos(2)*dcos(1)*( igdispderiv(2, 3)+igdispderiv(3, 2) )+   &
                    dcos(1)*dcos(2)*( igdispderiv(3, 1)+igdispderiv(1, 3) )
        call getJacobian( 361, nn, gp(k,:), ecoord, det, jacobian, inverse )
        dlength(k)=dsqrt(jacobian(1,1)*jacobian(1,1)+jacobian(1,2)*jacobian(1,2)+   &
                        jacobian(1,3)*jacobian(1,3))
    end do
    strain=0.5D0*(endstrain(1)+endstrain(2))
    istrain=0.5D0*(iendstrain(1)+iendstrain(2))
    ! if (iter.eq.1) then
    gausses(i)%fstatus(566)=strain-istrain  !dstrain
    gausses(i)%fstatus(564)=strain
    call calRCSteelUpdate(gausses(i)%istatus(61:65),gausses(i)%fstatus(561:593))
    ! end if
    stress1d=gausses(i)%fstatus(565)
    stress(1)=stress1d*dcos(1)*dcos(1)*gausses(i)%fstatus(593)
    stress(2)=stress1d*dcos(2)*dcos(2)*gausses(i)%fstatus(593)
    stress(3)=stress1d*dcos(3)*dcos(3)*gausses(i)%fstatus(593)
    stress(4)=stress1d*dcos(1)*dcos(2)*gausses(i)%fstatus(593)
    stress(5)=stress1d*dcos(2)*dcos(3)*gausses(i)%fstatus(593)
    stress(6)=stress1d*dcos(3)*dcos(1)*gausses(i)%fstatus(593)
    disp=0.0D0
    do k=1,2
        B=0.0D0
        do j = 1,nn
            B(1, 3*j-2) = gderiv(k,j, 1)
            B(2, 3*j-1) = gderiv(k,j, 2)
            B(3, 3*j  ) = gderiv(k,j, 3)
            B(4, 3*j-2) = gderiv(k,j, 2)
            B(4, 3*j-1) = gderiv(k,j, 1)
            B(5, 3*j-1) = gderiv(k,j, 3)
            B(5, 3*j  ) = gderiv(k,j, 2)
            B(6, 3*j-2) = gderiv(k,j, 3)
            B(6, 3*j  ) = gderiv(k,j, 1)
        end do
      qf(1:nn*ndof) = qf(1:nn*ndof)                                           &
                     +MATMUL(stress(1:6), B(1:6, 1:nn*ndof) )*dlength(k)*weight
    end do
end do

end subroutine

subroutine rotate_D_mat(E,A,dcos,D)
real(kind=kreal)                :: D(6,6)!constitutive matrix
real(kind=kreal)                :: E,A,l,m,n
real(kind=kreal)                :: dcos(3)
integer :: i,j

D=0.0D0
l=dcos(1)
m=dcos(2)
n=dcos(3)
D(1,1)=E*A*l*l*l*l
D(1,2)=E*A*l*l*m*m
D(1,3)=E*A*l*l*n*n
D(1,4)=E*A*l*l*l*m
D(1,5)=E*A*m*n*l*l
D(1,6)=E*A*n*l*l*l
D(2,2)=E*A*m*m*m*m
D(2,3)=E*A*m*m*n*n
D(2,4)=E*A*l*m*m*m
D(2,5)=E*A*m*m*m*n
D(2,6)=E*A*m*m*l*n
D(3,3)=E*A*n*n*n*n
D(3,4)=E*A*l*m*n*n
D(3,5)=E*A*n*n*n*m
D(3,6)=E*A*n*n*n*l
D(4,4)=E*A*l*l*m*m
D(4,5)=E*A*l*m*m*n
D(4,6)=E*A*l*l*m*n
D(5,5)=E*A*m*m*n*n
D(5,6)=E*A*m*n*n*l
D(6,6)=E*A*l*l*n*n

do j =6,1,-1
    do i=1,j-1
            D(j,i)=D(i,j)
    end do
end do

end subroutine

subroutine getIsoCoord(ecoord,p,gp)
    real(kind=kreal), intent(in)    :: ecoord(3, 8)
    real(kind=kreal)                :: p(3),gp(3),dcos(3)
    real(kind=kreal)                :: p1(3),p7(3),factor(3)
    real(kind=kreal)                :: length
    integer::i
    
    p1=ecoord(:, 1)
    p7=ecoord(:, 7)
    factor=2.0D0/(p7-p1)
    do i=1,3
        gp(i)=-1.0D0+(p(i)-p1(i))*factor(i)
    end do
    
end subroutine

subroutine calRebarDmat(istatus,fstatus,matvar,D)
    real (kreal) :: D(3,3),Dtemp(3,3)
    real (kreal) :: fstatus(33)
    real (kreal) :: E,ang
    integer :: istatus(5)
    real(kind=kreal) ::  matvar(:)
    real (kreal) :: rotmat(3,3)
    integer :: i
    Dtemp=0.0D0
    D=0.0D0
    rotmat=0.0D0
    
    if  (istatus(1)==0) then
        call calRCSteelInit (E,istatus,fstatus, &
                matvar,0.0D0)
    else
        E=fstatus(1)
    end if
    Dtemp(1,1)=E
    ang=matvar(7)*4*atan(1.0D0)/180.0D0
    rotmat(1,1)=dcos(ang)
    rotmat(1,2)=dsin(ang)
    rotmat(2,1)=-dsin(ang)
    rotmat(2,2)=dcos(ang)
    rotmat(3,3)=1.0D0
    
    D=matmul(transpose(rotmat),matmul(Dtemp,rotmat))
    
end subroutine

    subroutine calRCSteelInit(E,istatus,fstatus,matvar,length)
    real (kreal) :: fstatus(33)
    integer :: istatus(5)
    real(kind=kreal) ::  matvar(:)
    real(kind=kreal) ::  E,sigy
    real(kind=kreal) ::  length
    E = matvar(1)
    sigy = matvar(5)
    
    istatus (1) =   1       !initialise flag
    istatus (2) =   nint(matvar(6))     !LOOP
    istatus (3) =   1                   !INHOOP
    
    
    
    fstatus (1) =   E               !EHOOP
    fstatus (2) =   sigy            !SIGY
    fstatus (3) =   sigy/E          !SMAXH
    fstatus (17)=   fstatus (3)     !HTBE
    fstatus (18)=   fstatus (2)     !HTBS
    fstatus (21)=   -fstatus(3)     !HCBE
    fstatus (22)=   -fstatus(2)     !HCBS
    fstatus (32)=   E               !EH
    fstatus (33)=   length
    
    end subroutine

    subroutine calRCSteelUpdate(istatus, fstatus)
    real (kreal) :: fstatus(33)
    integer :: istatus(5)
    
    real (kreal) :: DEHOOPREAL
    real (kreal) :: TEHOOP
    integer :: NUMHY,NELMH
    integer :: NIT
    real (kreal) :: EH,SIGY,SMAXH,TSHOOP,EHOOP
    integer :: INHOOP,LOOP
    real (kreal) :: OEHOOP,OSHOOP
    real (kreal) :: HCE,HCS,HRE,HRS,HXE,HXS
    real (kreal) :: HTAE,HTAS,HTBE,HTBS,HCAE,HCAS
    real (kreal) :: HCBE,HCBS,EH1,EH2,HXNE,HXPE
    real (kreal) :: HAE,HAS,XMAX,XMIN,XX2
    integer :: INTC

    call read_fistr_rc (istatus, fstatus,ehoop,sigy,smaxh,tshoop,       &
            tehoop,dehoopreal,oehoop,oshoop,hce,hcs,hre,hrs,hxe,hxs,    &
            htae,htas,htbe,htbs,hcae,hcas,hcbe,hcbs,eh1,eh2,hxne,hxpe,  &
            hae,has,xmax,xmin,xx2,eh,loop,inhoop,intc,numhy)
            
    
    numhy=0
    if (loop.eq.0) then
        call sehoop (eh,sigy,smaxh,tehoop,tshoop,   &
                 ehoop,inhoop,numhy)
    elseif (loop.eq.2.or.loop.eq.3) then
        call hoopz(dehoopreal,tehoop,tshoop,eh,             &
                oehoop,oshoop,inhoop,ehoop,                 &
                hce,hcs,hre,hrs,hxe,hxs,                    &
                htae,htas,htbe,htbs,hcae,hcas,hcbe,hcbs,    &
                 eh1,eh2,hxne,hxpe,loop)
    
    elseif (loop.eq.4) then
        call mmpts( dehoopreal,tehoop,tshoop,eh,oehoop,oshoop, &
                    inhoop,ehoop,hre,hrs,hxe,hxs,hae,has,  &
                    htbe,htbs,hcbe,hcbs,xmax,xmin,xx2,intc)     
                 
    else
    write(*,*)"ENTER PROPER VALUE of LOOP",loop
    stop
    end if

    call write_fistr_rc (istatus, fstatus,ehoop,sigy,smaxh,tshoop,      &
            tehoop,dehoopreal,oehoop,oshoop,hce,hcs,hre,hrs,hxe,hxs,    &
            htae,htas,htbe,htbs,hcae,hcas,hcbe,hcbs,eh1,eh2,hxne,hxpe,  &
            hae,has,xmax,xmin,xx2,eh,loop,inhoop,intc,numhy)
    
    end subroutine
    
    subroutine read_fistr_rc (istatus, fstatus,ehoop,sigy,smaxh,tshoop, &
            tehoop,dehoopreal,oehoop,oshoop,hce,hcs,hre,hrs,hxe,hxs,    &
            htae,htas,htbe,htbs,hcae,hcas,hcbe,hcbs,eh1,eh2,hxne,hxpe,  &
            hae,has,xmax,xmin,xx2,eh,loop,inhoop,intc,numhy)
    
    real (kreal) :: fstatus(33)
    integer :: istatus(5)
    
    real (kreal) :: dehoopreal
    real (kreal) :: tehoop
    integer :: numhy,nelmh
    integer :: nit
    real (kreal) :: eh,sigy,smaxh,tshoop,ehoop
    integer :: inhoop,loop
    real (kreal) :: oehoop,oshoop
    real (kreal) :: hce,hcs,hre,hrs,hxe,hxs
    real (kreal) :: htae,htas,htbe,htbs,hcae,hcas
    real (kreal) :: hcbe,hcbs,eh1,eh2,hxne,hxpe
    real (kreal) :: hae,has,xmax,xmin,xx2
    integer :: intc
    
    
    
    ehoop       =   fstatus (   1   )
    sigy        =   fstatus (   2   )
    smaxh       =   fstatus (   3   )
    tehoop      =   fstatus (   4   )
    tshoop      =   fstatus (   5   )
    dehoopreal  =   fstatus (   6   )
    oehoop      =   fstatus (   7   )
    oshoop      =   fstatus (   8   )
    hce         =   fstatus (   9   )
    hcs         =   fstatus (   10  )
    hre         =   fstatus (   11  )
    hrs         =   fstatus (   12  )
    hxe         =   fstatus (   13  )
    hxs         =   fstatus (   14  )
    htae        =   fstatus (   15  )
    htas        =   fstatus (   16  )
    htbe        =   fstatus (   17  )
    htbs        =   fstatus (   18  )
    hcae        =   fstatus (   19  )
    hcas        =   fstatus (   20  )
    hcbe        =   fstatus (   21  )
    hcbs        =   fstatus (   22  )
    eh1         =   fstatus (   23  )
    eh2         =   fstatus (   24  )
    hxne        =   fstatus (   25  )
    hxpe        =   fstatus (   26  )
    hae         =   fstatus (   27  )
    has         =   fstatus (   28  )
    xmax        =   fstatus (   29  )
    xmin        =   fstatus (   30  )
    xx2         =   fstatus (   31  )
    eh          =   fstatus (   32  )


    loop        =   istatus (   2   )
    inhoop      =   istatus (   3   )
    intc        =   istatus (   4   )
    numhy       =   istatus (   5   )
    
    end subroutine

    subroutine write_fistr_rc (istatus, fstatus,ehoop,sigy,smaxh,tshoop,    &
            tehoop,dehoopreal,oehoop,oshoop,hce,hcs,hre,hrs,hxe,hxs,    &
            htae,htas,htbe,htbs,hcae,hcas,hcbe,hcbs,eh1,eh2,hxne,hxpe,  &
            hae,has,xmax,xmin,xx2,eh,loop,inhoop,intc,numhy)
    
    real (kreal) :: fstatus(33)
    integer :: istatus(5)
    
    real (kreal) :: DEHOOPREAL
    real (kreal) :: TEHOOP
    integer :: NUMHY,NELMH
    integer :: NIT
    real (kreal) :: EH,SIGY,SMAXH,TSHOOP,EHOOP
    integer :: INHOOP,LOOP
    real (kreal) :: OEHOOP,OSHOOP
    real (kreal) :: HCE,HCS,HRE,HRS,HXE,HXS
    real (kreal) :: HTAE,HTAS,HTBE,HTBS,HCAE,HCAS
    real (kreal) :: HCBE,HCBS,EH1,EH2,HXNE,HXPE
    real (kreal) :: HAE,HAS,XMAX,XMIN,XX2
    integer :: INTC
    
    fstatus (   1   )   =   EHOOP       
    fstatus (   2   )   =   SIGY        
    fstatus (   3   )   =   SMAXH       
    fstatus (   4   )   =   TEHOOP      
    fstatus (   5   )   =   TSHOOP      
    fstatus (   6   )   =   DEHOOPREAL  
    fstatus (   7   )   =   OEHOOP      
    fstatus (   8   )   =   OSHOOP      
    fstatus (   9   )   =   HCE         
    fstatus (   10  )   =   HCS         
    fstatus (   11  )   =   HRE         
    fstatus (   12  )   =   HRS         
    fstatus (   13  )   =   HXE         
    fstatus (   14  )   =   HXS         
    fstatus (   15  )   =   HTAE        
    fstatus (   16  )   =   HTAS        
    fstatus (   17  )   =   HTBE        
    fstatus (   18  )   =   HTBS        
    fstatus (   19  )   =   HCAE        
    fstatus (   20  )   =   HCAS        
    fstatus (   21  )   =   HCBE        
    fstatus (   22  )   =   HCBS        
    fstatus (   23  )   =   EH1         
    fstatus (   24  )   =   EH2         
    fstatus (   25  )   =   HXNE        
    fstatus (   26  )   =   HXPE        
    fstatus (   27  )   =   HAE         
    fstatus (   28  )   =   HAS         
    fstatus (   29  )   =   XMAX        
    fstatus (   30  )   =   XMIN        
    fstatus (   31  )   =   XX2         
    fstatus (   32  )   =   EH
    
    
    istatus (   2   )   =   LOOP        
    istatus (   3   )   =   INHOOP     
    istatus (   4   )   =   INTC       
    istatus (   5   )   =   NUMHY     
        
        
    end subroutine

    subroutine calConcreteDmat(istatus,fstatus,matvar,D)
        real (kreal) :: D(6,6)
        real (kreal) :: fstatus(560)
        integer :: istatus(58),i,j
        real(kind=kreal) ::  matvar(:)
        real (kreal) :: E,P
        
        D=0.0d0
        

        if (istatus(58)==0) then 
            call calConcreteInit(istatus,fstatus,matvar,D)
            istatus(58)=1       
        else if (istatus(58)==1) then
            call calConcreteInit(istatus,fstatus,matvar,D)
            istatus(58)=2
        else
        D=0.0d0
        
        !Newton Rhapson

        D(1,:)  =   fstatus (   525 :   530 )
        D(2,:)  =   fstatus (   531 :   536 )
        D(3,:)  =   fstatus (   537 :   542 )
        D(4,:)  =   fstatus (   543 :   548 )
        D(5,:)  =   fstatus (   549 :   554 )
        D(6,:)  =   fstatus (   555 :   560 )
        
        end if

    end

    subroutine calConcreteInit(istatus,fstatus,matvar,D)
    
    
    real (kreal) :: fstatus(560)
    integer :: istatus(58)
    real (kreal) :: D(6,6)
    real(kind=kreal) :: E,P,CLC,PAI,ECR,VOLUME
    real(kind=kreal), intent(in)  :: matvar(:)
    
    D=0.0D0
    E = matvar(1) 
    P = matvar(2)
    
    D(1,1)=E*(1.d0-P)/(1.d0-2.d0*P)/(1.d0+P)
    D(1,2)=E*P/(1.d0-2.d0*P)/(1.d0+P)
    D(1,3)=D(1,2)
    D(2,1)=D(1,2)
    D(2,2)=D(1,1)
    D(2,3)=D(1,2)
    D(3,1)=D(1,3)
    D(3,2)=D(2,3)
    D(3,3)=D(1,1)
    D(4,4)=E/(1.d0+P)*0.5d0
    D(5,5)=E/(1.d0+P)*0.5d0
    D(6,6)=E/(1.d0+P)*0.5d0
    
    istatus(14)     =   nint(matvar(5) )        !IROT
    istatus(16)     =   nint(matvar(6) )        !ICIC(1)
    istatus(17)     =   nint(matvar(7) )        !ICIC(2)
    istatus(18)     =   nint(matvar(8) )        !ICIC(3)
    istatus(19)     =   nint(matvar(9) )        !ICIC(4)
    istatus(20)     =   nint(matvar(10))        !ICIC(5)
    istatus(21)     =   nint(matvar(11))        !ICIC(6)
    istatus(22)     =   nint(matvar(12))        !ICIC(7)
    istatus(23)     =   nint(matvar(13))        !ICIC(8)
    istatus(15)     =   nint(matvar(14))        !ICONCYC
            
    fstatus(176)    =   matvar(15)              !FC
    fstatus(175)    =   matvar(16)              !ECU
    fstatus(177)    =   matvar(17)              !FT
    fstatus(191)    =   matvar(18)              !EPCU
    fstatus(192)    =   matvar(19)              !EPCUS
    fstatus(211)    =   matvar(20)              !C_IC4 !0.6
    
   !*******************************************************

    fstatus(43:45)=E                        !EC                 
    fstatus(61:63)=P                        !PR 
    !fstatus(64:66)=E/(2*(1+P))             !RG     !d
    fstatus(178)=fstatus(177)               !SCT    !d  
    fstatus(179)=E                          !EO     !d
    
    !Mesh dependance parameter
    if (istatus(58)==0) then 
        PAI=ASIN(1.0D0)*2.0D0
        VOLUME=fstatus(180)/8
        CLC=2.0*(3.0*VOLUME/4/PAI)**(1/3.0D0)*10.
        ECR=matvar(17)/matvar(1) 
        ! EBUENG(M)=(0.23*FC*0.0981+136.)/(284.*SCT*0.0981*CLC(M)*1.0)+ECR
        fstatus(180)=(0.23*matvar(15)*0.0981+136.)/(284.*fstatus(178)*0.0981*CLC*1.0)+ECR
    end if
    fstatus(183)=(E+E-2*0.02*E- &
        (0.02*SQRT(E)+0.02*SQRT(E))**2)/4.  !G0     !d  
    
    fstatus(184)=fstatus(179)/100           !EC_100 !d
    fstatus(193:195)=fstatus(191)           !EPCU_  !d
    fstatus(196:198)=fstatus(192)           !EPCUS_ !d
    fstatus(199:201)=1.0D0                  !SC_3D
    fstatus(208:210)=fstatus(61)            !VUE    !d  
    fstatus(400) = 1.0D0                    !CYCN(1,27)
    fstatus(449) = 1.0D0                    !CYCN(1,27)
    fstatus(498) = 1.0D0                    !CYCN(1,27)
    fstatus(374) = 1                        !CYCN(1,27) !NST
    fstatus(423) = 1                        !CYCN(1,27) !NST
    fstatus(472) = 1                        !CYCN(1,27) !NST

    !CHECK
    fstatus (157:165)=-2.00D0               !DC_2
    fstatus(181)=0.6D0                      !VS 
    istatus(10:12)=1                        !IUL    
    istatus(13)=1                           !INCF   
    istatus(25)=1                           !IMM
    istatus(26)=1                           !NIT
    istatus(28)=1                           !NST
    istatus(52)=2                           !M
    istatus(53)=2                           !LX
    istatus(54)=2                           !LY
    istatus(55)=2                           !LZ 
    
    end subroutine

    subroutine calUpdateConcrete(dmatrix,istatus,fstatus)
    
    !fistr variables
    real (kreal) :: fstatus(560)
    integer :: istatus(58)
    real (kreal) :: stress(6), strain(6), dmatrix(6,6)
    
    
    !Concrete variables
    real(kreal) TSI(6),TEP(6),TEPLL(6),TEPC(6),TEPL(6)
    real(kreal) SHTEP(3),SHTEPL(3)
    real(kreal) EU(6),SN(6),SNCR(6)
    real(kreal) ET(3),SCOLD(3)
    real(kreal) EUOUT(6),SNCROUT(6)
    integer :: IULOUT(3),KULOUT(3)
    real(kreal) PR12,PR23,PR31,RG12,RG23,RG31
    real(kreal) DCC(3,3),DCCL(3,3),DCCR(3,3),DCCRL(3,3)
    real(kreal) DCCRCHG(3,3,3),DCCRCHGL(3,3,3)
    real(kreal) DC_1(3,3),DC_2(3,3)
    integer ICOORDCHGED,IACT
    real(kreal) ECU,FC,FT,SCT,EO,EBU,VS,SENS,G0,EC_100
    real(kreal) SC(3),EC(3)
    integer IIKEEPSC
    real(kreal) EPCU,EPCUS,EPCU_(3),EPCUS_(3)
    real(kreal) SC_3D(3),SC_RAM(3)
    real(kreal) EESC(3),VUE(3)
    integer IUL(3)
    integer INCF
    integer IROT,ICONCYC  
    integer ICIC(8)
    real(kreal) C_IC4
    integer ICR,IMM,NIT,IR,NST
    real(kreal) EEN(3),SEN(3),ERC(3),SRC(3),EXC(3),SXC(3)
    real(kreal) EBC(3),SBC(3),EJ(3),SJ(3),EJJ(3)
    real(kreal) ERT(3),SRT(3),ECP(3),SCP(3)
    real(kreal) ETP(3),STP(3),EPC(3),SPC(3),EPC1ST(3)
    real(kreal) EPT(3),ETMAX(3),STMAX(3)
    real(kreal) ECCR(3),SCCR(3)
    real(kreal) EUTMAX(3),EUOVER(3)
    integer ICOUNTLOOPT(3)
    real(kreal) E41_51(3),S41_51(3)
    integer ID_CRACK(3)
    real(kreal) EUOLD(6),EU_OLD(6),SNOLD(6),SNCROLD(6)
    real(kreal) EPPC(3),EPPT(3),EPEC(3),CCRACK(3)
    integer IVIRGIN(3)
    real(kreal) ELIMIT(3)
    integer NTCRACK
    integer LLCRACK(3)
    integer NTFAIL
    integer :: KUL(3),KULIN1(3),KULIN5(3)
    real(kreal) DS1(3),DE1(3),DS2(3),DE2(3)
    real(kreal) FS1(3),FE1(3),FS2(3),FE2(3)
    real(kreal) RS1(3),RE1(3),RS2(3),RE2(3)
    integer M,LX,LY,LZ
    real(kreal) CYCN(3,49)
    real(kreal) EUDL(3)
    integer Mhj,NG3C    
    real (kreal) :: GBETA
    
    !added by shan 19/01/2018
    real (kreal) :: SC_RAM1(3)
    
    
    real (kreal),dimension (6,6) :: TE,TET,UC,UTI,DMCP
    real (kreal) :: DDCC(3,3)
    real (kreal) :: E1,E2,E3,G12,G23,G31,EUC1,EUC2,EUC3,F,G0_100
    integer :: NDDCC,IJ,I,J

    DMCP=0.0d0
    dmatrix=0.0d0
    
    ! read the values of fstatus
    call read_fistr(istatus,fstatus,                                                        &
            TSI,TEP,TEPLL,SHTEP,SHTEPL,EU,SN,SNCR,ET,SCOLD,EUOUT,SNCROUT,IULOUT,KULOUT,     &
            PR23,PR31,PR12,RG23,RG31,RG12,DCC,DCCL,DCCR,DCCRL,DCCRCHG,DCCRCHGL,DC_1,DC_2,   &
            ICOORDCHGED,IACT,ECU,FC,FT,SCT,EO,EBU,VS,SENS,G0,EC_100,SC,EC,IIKEEPSC,         &
            EPCU,EPCUS,EPCU_,EPCUS_,SC_3D,SC_RAM,EESC,VUE,IUL,INCF,IROT,ICONCYC,            &
            ICIC,C_IC4,ICR,IMM,NIT,IR,NST,EEN,SEN,ERC,SRC,EXC,SXC,EBC,SBC,EJ,SJ,EJJ,        &
            ERT,SRT,ECP,SCP,ETP,STP,EPC,SPC,EPC1ST,EPT,ETMAX,STMAX,ECCR,SCCR,EUTMAX,EUOVER, &                                       
            ICOUNTLOOPT,E41_51,S41_51,ID_CRACK,EUOLD,EU_OLD,SNOLD,SNCROLD,EPPC,EPPT,EPEC,   &                                       
            CCRACK,IVIRGIN,ELIMIT,NTCRACK,LLCRACK,NTFAIL,KUL,KULIN1,KULIN5,DS1,DE1,DS2,DE2, &                                      
            FS1,FE1,FS2,FE2,RS1,RE1,RS2,RE2,M,LX,LY,LZ,CYCN,EUDL,Mhj,NG3C,GBETA) 
    
    
    if (nit.eq.1) SC_RAM=EU(1:3)
    
    
    ! IACTOLD=IACTIVEC !Check
    
    SHTEP   =   0.0D0
    SHTEPL  =   0.0D0
    SCOLD   =   0.0D0
    EUOUT   =   0.0D0
    SNCROUT =   0.0D0
    EPCU_   =   0.0D0
    EPCUS_  =   0.0D0
    SC_RAM1 =   0.0D0
    CYCN(1:3,30:33)=2 !element and gauss number
    
    INCF    =   1
    IULOUT  =   0
    KULOUT  =   0
    IR      =   0
    NTCRACK =   0
    LLCRACK =   0
    NTFAIL  =   0


        
    call CONCRT (                                                                           &   
            TSI,TEP,TEPLL,SHTEP,SHTEPL,EU,SN,SNCR,ET,SCOLD,EUOUT,SNCROUT,IULOUT,KULOUT      &
            ,PR23,PR31,PR12,RG23,RG31,RG12,DCC,DCCL,DCCR,DCCRL,DCCRCHG,DCCRCHGL,DC_1,DC_2,  &
            ICOORDCHGED,IACT,ECU,FC,FT,SCT,EO,EBU,VS,SENS,G0,EC_100,SC,EESC,IIKEEPSC,           &
            EPCU,EPCUS,EPCU_,EPCUS_,SC_3D,SC_RAM1,EESC,VUE,IUL,INCF,IROT,ICONCYC,           &
            ICIC,C_IC4,ICR,IMM,NIT,IR,NST,EEN,SEN,ERC,SRC,EXC,SXC,EBC,SBC,EJ,SJ,EJJ,        &
            ERT,SRT,ECP,SCP,ETP,STP,EPC,SPC,EPC1ST,EPT,ETMAX,STMAX,ECCR,SCCR,EUTMAX,EUOVER, &                                       
            ICOUNTLOOPT,E41_51,S41_51,ID_CRACK,EUOLD,EU_OLD,SNOLD,SNCROLD,EPPC,EPPT,EPEC,   &                                       
            CCRACK,IVIRGIN,ELIMIT,NTCRACK,LLCRACK,NTFAIL,KUL,KULIN1,KULIN5,DS1,DE1,DS2,DE2, &                                      
            FS1,FE1,FS2,FE2,RS1,RE1,RS2,RE2,M,LX,LY,LZ,CYCN,EUDL,Mhj,NG3C)    

    ! open (30,file='integerlist.txt')
    ! write(30,333)((nst-1)*2)+nit,nst,nit,iulout,kulout,icoordchged,iact,iikeepsc,iul,incf,irot, &
                    ! iconcyc,icr,imm,ir,icountloopt,        &
                    ! id_crack,ivirgin,ntcrack,llcrack,ntfail,kul,kulin1, &
                    ! kulin5,icic
    
    
    ! 333 format (i3,"nst",i3,"nit",i3,"iulout",3i3,"kulout",3i3, &
            ! "icoordchged",i3,"iact",i3,"iikeepsc",i3,"iul",3i3, &
            ! "incf",i3,"irot",i3,"iconcyc",i3,"icr",i3,"imm",i3,"ir",i3, &
            ! "icountloopt",3i3,"id_crack",3i3,"ivirgin",3i3,    &
            ! "ntcrack",i3,"llcrack",3i3,"ntfail",i3,"kul",3i3,"kulin1",  &
            ! 3i3,"kulin5",3i3,"icic",8i3)

            
            
! C *** 材料軸での応力σiを，�E体座標系でのσijに変換する ***

    IF(ICIC(6).EQ.0)        GOTO 3100            !     ICIC(6)= 0 : *回転ひび割れモチE��, 1 : *1方向�Eび割れモチE�� 1 : *多方向�Eび割れモチE��
    IF(NTCRACK.EQ.0)        GOTO 3100
        
    ! !ICIC(6)= 0 回転ひび割れモチE��の時�E、ここを通らなぁE��by sun
      ! CALL CHG3 (TSI(1),TSI(2),TSI(3),TSI(4),TSI(5),TSI(6),&
                ! DCCR(1,1),DCCR(1,2),DCCR(1,3),    &
                ! DCCR(2,1),DCCR(2,2),DCCR(2,3),    &
                ! DCCR(3,1),DCCR(3,2),DCCR(3,3),    &
                ! SNCR(1),  SNCR(2),  SNCR(3),  &
                ! SNCR(4),  SNCR(5),  SNCR(6),0)    !ICIC(6)= 0 回転ひび割れモチE��の時�E、ここを通らなぁE��by sun
      CALL CHG3 (TSI(1),TSI(2),TSI(3),TSI(4),TSI(5),TSI(6),&
                DCCR(1,1),DCCR(1,2),DCCR(1,3),  &
                DCCR(2,1),DCCR(2,2),DCCR(2,3),  &
                DCCR(3,1),DCCR(3,2),DCCR(3,3),  &
                SNCR(1),  SNCR(2),  SNCR(3),    &
                SNCR(4),  SNCR(5),  SNCR(6),0)             ! mark by sun 2017.1122
    
    
! C *** CHECK CONCRETE FAILURE AND STRAIN SOFTENING 
        GOTO 3110

3100  CALL CHG3 (TSI(1),TSI(2),TSI(3),TSI(4),TSI(5),TSI(6),&
                DCC(1,1),DCC(1,2),DCC(1,3),  &
                DCC(2,1),DCC(2,2),DCC(2,3),  &
                DCC(3,1),DCC(3,2),DCC(3,3),  &
                SN(1),SN(2),SN(3),  &
                SN(4),SN(5),SN(6),0)
    
! C *** FORM D-MATRIX OF INTEGRAL POINTS OF CONCRETE

 3110 CONTINUE

! C ***初期値めE.0に設定すめE**
       TE=0.0D0
       TET=0.0D0
       UC=0.0D0
       UTI=0.0D0
       dmatrix=0.0D0
       


    if (ET(1)<EC_100)ET(1)=EC_100
    if (ET(2)<EC_100)ET(2)=EC_100
    if (ET(3)<EC_100)ET(3)=EC_100

! C ***せん断剛性を算�Eする***
      F = 1.0D0-PR23*PR23-PR31*PR31-PR12*PR12-2.0D0*PR23*PR31*PR12
      E1 = ET(1)
      E2 = ET(2)
      E3 = ET(3)
      G12 = (E1+E2-2*PR12*SQRT(E1*E2)-                  &
           (PR23*SQRT(E1)+PR31*SQRT(E2))**2)/4./F
      G23 = (E2+E3-2*PR23*SQRT(E2*E3)-                  &
           (PR31*SQRT(E2)+PR12*SQRT(E3))**2)/4./F
      G31 = (E3+E1-2*PR31*SQRT(E3*E1)-                  &
           (PR12*SQRT(E3)+PR23*SQRT(E1))**2)/4./F

           
    EUC1 = EU(1)
    EUC2 = EU(2)
    EUC3 = EU(3)

    IF(ICIC(6).EQ.0)        GOTO 3510
    IF(NTCRACK.EQ.0)        GOTO 3510
      GO TO 3520

! C *** EQUIVALENT UNIAXIAL STRAIN MODEL BY DARWIN & PECKNOLD            !等価一軸ひずみモチE��

 3510   DDCC=DCC
        NDDCC=0

      UTI(1,1) = E1*(1.0D0-PR23*PR23)/F
      UTI(1,2) = SQRT(E1*E2)*(PR31*PR23+PR12)/F
      UTI(1,3) = SQRT(E1*E3)*(PR12*PR23+PR31)/F

      UTI(2,1) = UTI(1,2)
      UTI(2,2) = E2*(1.0D0-PR31*PR31)/F
      UTI(2,3) = SQRT(E2*E3)*(PR12*PR31+PR23)/F

      UTI(3,1) = UTI(1,3)
      UTI(3,2) = UTI(2,3)
      UTI(3,3) = E3*(1.0D0-PR12*PR12)/F

      UTI(4,4) = G12
      UTI(5,5) = G23
      UTI(6,6) = G31

      GO TO 3540


! C *** CRACK DIRECTION CORDINATE MODEL
! C ***二方向�Eび割れ�E場吁E**
 3520 CONTINUE

    DDCC=DCCR
    NDDCC=1

    CALL RG12RG23RG31(NTCRACK,LLCRACK,EUC1,EUC2,EUC3,G0,ICIC,RG12,RG23,RG31,GBETA)
    G0_100 = G0 / 100.
    IF( RG12 < G0_100 ) RG12 = G0_100
    IF( RG23 < G0_100 ) RG23 = G0_100
    IF( RG31 < G0_100 ) RG31 = G0_100

    
    IF(NTCRACK.EQ.1) GOTO 3523

 3522 UTI(1,1) = E1*(1.0D0-PR23*PR23)/F
      UTI(2,2) = E2*(1.0D0-PR31*PR31)/F
      UTI(3,3) = E3*(1.0D0-PR12*PR12)/F
      UTI(4,4) = RG12
      UTI(5,5) = RG23
      UTI(6,6) = RG31

      GOTO 3540

 3523 IF(LLCRACK(1).EQ.1) THEN
        UTI(1,1) = E1*(1.0D0-PR23*PR23)/F
        UTI(2,2) = E2/F
        UTI(2,3) = SQRT(E2*E3)*PR23/F
        UTI(3,2) = UTI(2,3)
        UTI(3,3) = E3/F
        UTI(4,4) = RG12
        UTI(5,5) =  G23
        UTI(6,6) = RG31
      ENDIF

      IF(LLCRACK(2).EQ.1) THEN
        UTI(1,1) = E1/F
        UTI(1,3) = SQRT(E1*E3)*PR31/F
        UTI(2,2) = E2*(1.0D0-PR31*PR31)/F
        UTI(3,1) = UTI(1,3)
        UTI(3,3) = E3/F
        UTI(4,4) = RG12
        UTI(5,5) = RG23
        UTI(6,6) =  G31
      ENDIF

      IF(LLCRACK(3).EQ.1) THEN
        UTI(1,1) = E1/F
        UTI(1,2) = SQRT(E1*E2)*PR12/F
        UTI(2,1) = UTI(1,2)
        UTI(2,2) = E2/F
        UTI(3,3) = E3*(1.0D0-PR12*PR12)/F
        UTI(4,4) = G12
        UTI(5,5) = RG23
        UTI(6,6) = RG31
      ENDIF

 3540 CONTINUE
 
    TE(1,1)=DDCC(1,1)*DDCC(1,1)
    TE(1,2)=DDCC(1,2)*DDCC(1,2)
    TE(1,3)=DDCC(1,3)*DDCC(1,3)
    TE(1,4)=DDCC(1,1)*DDCC(1,2)
    TE(1,5)=DDCC(1,2)*DDCC(1,3)
    TE(1,6)=DDCC(1,3)*DDCC(1,1)
    TE(2,1)=DDCC(2,1)*DDCC(2,1)
    TE(2,2)=DDCC(2,2)*DDCC(2,2)
    TE(2,3)=DDCC(2,3)*DDCC(2,3)
    TE(2,4)=DDCC(2,1)*DDCC(2,2)
    TE(2,5)=DDCC(2,2)*DDCC(2,3)
    TE(2,6)=DDCC(2,3)*DDCC(2,1)
    TE(3,1)=DDCC(3,1)*DDCC(3,1)
    TE(3,2)=DDCC(3,2)*DDCC(3,2)
    TE(3,3)=DDCC(3,3)*DDCC(3,3)
    TE(3,4)=DDCC(3,1)*DDCC(3,2)
    TE(3,5)=DDCC(3,2)*DDCC(3,3)
    TE(3,6)=DDCC(3,3)*DDCC(3,1)
    TE(4,1)=DDCC(2,1)*DDCC(3,1)*2
    TE(4,2)=DDCC(2,2)*DDCC(3,2)*2
    TE(4,3)=DDCC(2,3)*DDCC(3,3)*2
    TE(4,4)=DDCC(2,1)*DDCC(3,2)+DDCC(3,1)*DDCC(2,2)
    TE(4,5)=DDCC(2,2)*DDCC(3,3)+DDCC(3,2)*DDCC(2,3)
    TE(4,6)=DDCC(2,3)*DDCC(3,1)+DDCC(3,3)*DDCC(2,1)
    TE(5,1)=DDCC(3,1)*DDCC(1,1)*2
    TE(5,2)=DDCC(3,2)*DDCC(1,2)*2
    TE(5,3)=DDCC(3,3)*DDCC(1,3)*2
    TE(5,4)=DDCC(3,1)*DDCC(1,2)+DDCC(1,1)*DDCC(3,2)
    TE(5,5)=DDCC(3,2)*DDCC(1,3)+DDCC(1,2)*DDCC(3,3)
    TE(5,6)=DDCC(3,3)*DDCC(1,1)+DDCC(1,3)*DDCC(3,1)
    TE(6,1)=DDCC(1,1)*DDCC(2,1)*2
    TE(6,2)=DDCC(1,2)*DDCC(2,2)*2
    TE(6,3)=DDCC(1,3)*DDCC(2,3)*2
    TE(6,4)=DDCC(1,1)*DDCC(2,2)+DDCC(2,1)*DDCC(1,2)
    TE(6,5)=DDCC(1,2)*DDCC(2,3)+DDCC(2,2)*DDCC(1,3)
    TE(6,6)=DDCC(1,3)*DDCC(2,1)+DDCC(2,3)*DDCC(1,1)
    

      DO I=1,6 ; DO J=1,6
        TET(I,J)=TE(J,I)
    ENDDO ; ENDDO

        DO I=1,6 ; DO J=1,6 ; DO IJ=1,6
           UC(I,J)=UC(I,J)+TET(I,IJ)*UTI(IJ,J)
    ENDDO ; ENDDO ; ENDDO


        DO I=1,6 ;DO J=1,6 ; DO IJ=1,6
           DMCP(I,J)=DMCP(I,J)+UC(I,IJ)*TE(IJ,J)
    ENDDO ; ENDDO ; ENDDO
    

      ! open (31,file='floatlist.txt')
      ! ! write(31,'(36f10.0)')UTI
      ! write(31,'(2i3,3f10.2,6e15.6)')NST,NIT,E1,E2,E3,TEP
      
      
      
    dmatrix=DMCP
    fstatus(1:6)=TSI
    !UPDATE LAST STEP VALUES
    DCCRL=DCCR
    DCCL=DCC
    TEPLL=TEP
    EUOLD=EU
    SNOLD=SN
    SNCROLD=SNCR
    
    call write_fistr(istatus,fstatus,                                                       &
            TSI,TEP,TEPLL,SHTEP,SHTEPL,EU,SN,SNCR,ET,SCOLD,EUOUT,SNCROUT,IULOUT,KULOUT,     &
            PR23,PR31,PR12,RG23,RG31,RG12,DCC,DCCL,DCCR,DCCRL,DCCRCHG,DCCRCHGL,DC_1,DC_2,   &
            ICOORDCHGED,IACT,ECU,FC,FT,SCT,EO,EBU,VS,SENS,G0,EC_100,SC,EC,IIKEEPSC,         &
            EPCU,EPCUS,EPCU_,EPCUS_,SC_3D,SC_RAM,EESC,VUE,IUL,INCF,IROT,ICONCYC,            &
            ICIC,C_IC4,ICR,IMM,NIT,IR,NST,EEN,SEN,ERC,SRC,EXC,SXC,EBC,SBC,EJ,SJ,EJJ,        &
            ERT,SRT,ECP,SCP,ETP,STP,EPC,SPC,EPC1ST,EPT,ETMAX,STMAX,ECCR,SCCR,EUTMAX,EUOVER, &                                       
            ICOUNTLOOPT,E41_51,S41_51,ID_CRACK,EUOLD,EU_OLD,SNOLD,SNCROLD,EPPC,EPPT,EPEC,   &                                       
            CCRACK,IVIRGIN,ELIMIT,NTCRACK,LLCRACK,NTFAIL,KUL,KULIN1,KULIN5,DS1,DE1,DS2,DE2, &                                      
            FS1,FE1,FS2,FE2,RS1,RE1,RS2,RE2,M,LX,LY,LZ,CYCN,EUDL,Mhj,NG3C,GBETA)
    
    !change here for method
    
    ! if (NST==1) then  !Initial Stiffness Method
    fstatus (   525 :   530 ) = DMCP(1,:)
    fstatus (   531 :   536 ) = DMCP(2,:)
    fstatus (   537 :   542 ) = DMCP(3,:)
    fstatus (   543 :   548 ) = DMCP(4,:)
    fstatus (   549 :   554 ) = DMCP(5,:)
    fstatus (   555 :   560 ) = DMCP(6,:)
    ! end if

    end subroutine

!****************************
!read the values from fistr variables
!****************************   
    
    subroutine read_fistr(istatus,fstatus,                                                  &
            TSI,TEP,TEPLL,SHTEP,SHTEPL,EU,SN,SNCR,ET,SCOLD,EUOUT,SNCROUT,IULOUT,KULOUT,     &
            PR23,PR31,PR12,RG23,RG31,RG12,DCC,DCCL,DCCR,DCCRL,DCCRCHG,DCCRCHGL,DC_1,DC_2,   &
            ICOORDCHGED,IACT,ECU,FC,FT,SCT,EO,EBU,VS,SENS,G0,EC_100,SC,EC,IIKEEPSC,         &
            EPCU,EPCUS,EPCU_,EPCUS_,SC_3D,SC_RAM,EESC,VUE,IUL,INCF,IROT,ICONCYC,            &
            ICIC,C_IC4,ICR,IMM,NIT,IR,NST,EEN,SEN,ERC,SRC,EXC,SXC,EBC,SBC,EJ,SJ,EJJ,        &
            ERT,SRT,ECP,SCP,ETP,STP,EPC,SPC,EPC1ST,EPT,ETMAX,STMAX,ECCR,SCCR,EUTMAX,EUOVER, &                                       
            ICOUNTLOOPT,E41_51,S41_51,ID_CRACK,EUOLD,EU_OLD,SNOLD,SNCROLD,EPPC,EPPT,EPEC,   &                                       
            CCRACK,IVIRGIN,ELIMIT,NTCRACK,LLCRACK,NTFAIL,KUL,KULIN1,KULIN5,DS1,DE1,DS2,DE2, &                                      
            FS1,FE1,FS2,FE2,RS1,RE1,RS2,RE2,M,LX,LY,LZ,CYCN,EUDL,Mhj,NG3C,GBETA) 
    
        real (kreal) :: fstatus(560)
        integer :: istatus(58)
        
        real(kreal) TSI(6),TEP(6),TEPLL(6),TEPC(6),TEPL(6)
        real(kreal) SHTEP(3),SHTEPL(3)
        real(kreal) EU(6),SN(6),SNCR(6)
        real(kreal) ET(3),SCOLD(3)
        real(kreal) EUOUT(6),SNCROUT(6)
        integer :: IULOUT(3),KULOUT(3)
        real(kreal) PR12,PR23,PR31,RG12,RG23,RG31
        real(kreal) DCC(3,3),DCCL(3,3),DCCR(3,3),DCCRL(3,3)
        real(kreal) DCCRCHG(3,3,3),DCCRCHGL(3,3,3)
        real(kreal) DC_1(3,3),DC_2(3,3)
        integer ICOORDCHGED,IACT
        real(kreal) ECU,FC,FT,SCT,EO,EBU,VS,SENS,G0,EC_100
        real(kreal) SC(3),EC(3)
        integer IIKEEPSC
        real(kreal) EPCU,EPCUS,EPCU_(3),EPCUS_(3)
        real(kreal) SC_3D(3),SC_RAM(3)
        real(kreal) EESC(3),VUE(3)
        integer IUL(3)
        integer INCF
        integer IROT,ICONCYC  
        integer ICIC(8)
        real(kreal) C_IC4
        integer ICR,IMM,NIT,IR,NST
        real(kreal) EEN(3),SEN(3),ERC(3),SRC(3),EXC(3),SXC(3)
        real(kreal) EBC(3),SBC(3),EJ(3),SJ(3),EJJ(3)
        real(kreal) ERT(3),SRT(3),ECP(3),SCP(3)
        real(kreal) ETP(3),STP(3),EPC(3),SPC(3),EPC1ST(3)
        real(kreal) EPT(3),ETMAX(3),STMAX(3)
        real(kreal) ECCR(3),SCCR(3)
        real(kreal) EUTMAX(3),EUOVER(3)
        integer ICOUNTLOOPT(3)
        real(kreal) E41_51(3),S41_51(3)
        integer ID_CRACK(3)
        real(kreal) EUOLD(6),EU_OLD(6),SNOLD(6),SNCROLD(6)
        real(kreal) EPPC(3),EPPT(3),EPEC(3),CCRACK(3)
        integer IVIRGIN(3)
        real(kreal) ELIMIT(3)
        integer NTCRACK
        integer LLCRACK(3)
        integer NTFAIL
        integer :: KUL(3),KULIN1(3),KULIN5(3)
        real(kreal) DS1(3),DE1(3),DS2(3),DE2(3)
        real(kreal) FS1(3),FE1(3),FS2(3),FE2(3)
        real(kreal) RS1(3),RE1(3),RS2(3),RE2(3)
        integer M,LX,LY,LZ
        real(kreal) CYCN(3,49)
        real(kreal) EUDL(3)
        integer Mhj,NG3C    
        real (kreal) :: GBETA
    
    
        ! !FLOAT VALUES
    
    TSI             =   fstatus (   1   :   6   )       
    TEP             =   fstatus (   7   :   12  )       
    TEPLL           =   fstatus (   13  :   18  )       
    SHTEP           =   fstatus (   19  :   21  )      
    SHTEPL          =   fstatus (   22  :   24  )      
    EU              =   fstatus (   25  :   30  )      
    SN              =   fstatus (   31  :   36  )      
    SNCR            =   fstatus (   37  :   42  )      
    ET              =   fstatus (   43  :   45  )      
    SCOLD           =   fstatus (   46  :   48  )      
    EUOUT           =   fstatus (   49  :   54  )      
    SNCROUT         =   fstatus (   55  :   60  )      
    PR23            =   fstatus (   61  )      
    PR31            =   fstatus (   62  )      
    PR12            =   fstatus (   63  )      
    RG23            =   fstatus (   64  )      
    RG31            =   fstatus (   65  )      
    RG12            =   fstatus (   66  )      
    DCC(1,:)        =   fstatus (   67  :   69  )      
    DCC(2,:)        =   fstatus (   70  :   72  )      
    DCC(3,:)        =   fstatus (   73  :   75  )      
    DCCL(1,:)       =   fstatus (   76  :   78  )      
    DCCL(2,:)       =   fstatus (   79  :   81  )      
    DCCL(3,:)       =   fstatus (   82  :   84  )      
    DCCR(1,:)       =   fstatus (   85  :   87  )      
    DCCR(2,:)       =   fstatus (   88  :   90  )      
    DCCR(3,:)       =   fstatus (   91  :   93  )      
    DCCRL(1,:)      =   fstatus (   94  :   96  )      
    DCCRL(2,:)      =   fstatus (   97  :   99  )      
    DCCRL(3,:)      =   fstatus (   100 :   102 )      
    DCCRCHG(1,1,:)  =   fstatus (   103 :   105 )      
    DCCRCHG(1,2,:)  =   fstatus (   106 :   108 )      
    DCCRCHG(1,3,:)  =   fstatus (   109 :   111 )      
    DCCRCHG(2,1,:)  =   fstatus (   112 :   114 )      
    DCCRCHG(2,2,:)  =   fstatus (   115 :   117 )      
    DCCRCHG(2,3,:)  =   fstatus (   118 :   120 )      
    DCCRCHG(3,1,:)  =   fstatus (   121 :   123 )      
    DCCRCHG(3,2,:)  =   fstatus (   124 :   126 )      
    DCCRCHG(3,3,:)  =   fstatus (   127 :   129 )      
    DCCRCHGL(1,1,:) =   fstatus (   130 :   132 )      
    DCCRCHGL(1,2,:) =   fstatus (   133 :   135 )      
    DCCRCHGL(1,3,:) =   fstatus (   136 :   138 )      
    DCCRCHGL(2,1,:) =   fstatus (   139 :   141 )      
    DCCRCHGL(2,2,:) =   fstatus (   142 :   144 )      
    DCCRCHGL(2,3,:) =   fstatus (   145 :   147 )      
    DCCRCHGL(3,1,:) =   fstatus (   148 :   150 )      
    DCCRCHGL(3,2,:) =   fstatus (   151 :   153 )      
    DCCRCHGL(3,3,:) =   fstatus (   154 :   156 )      
    DC_1(1,:)       =   fstatus (   157 :   159 )      
    DC_1(2,:)       =   fstatus (   160 :   162 )      
    DC_1(3,:)       =   fstatus (   163 :   165 )      
    DC_2(1,:)       =   fstatus (   166 :   168 )      
    DC_2(2,:)       =   fstatus (   169 :   171 )      
    DC_2(3,:)       =   fstatus (   172 :   174 )      
    ECU             =   fstatus (   175 )      
    FC              =   fstatus (   176 )      
    FT              =   fstatus (   177 )      
    SCT             =   fstatus (   178 )      
    EO              =   fstatus (   179 )      
    EBU             =   fstatus (   180 )      
    VS              =   fstatus (   181 )      
    SENS            =   fstatus (   182 )      
    G0              =   fstatus (   183 )      
    EC_100          =   fstatus (   184 )      
    SC              =   fstatus (   185 :   187 )      
    EC              =   fstatus (   188 :   190 )      
    EPCU            =   fstatus (   191 )      
    EPCUS           =   fstatus (   192 )      
    EPCU_           =   fstatus (   193 :   195 )      
    EPCUS_          =   fstatus (   196 :   198 )      
    SC_3D           =   fstatus (   199 :   201 )      
    SC_RAM          =   fstatus (   202 :   204 )      
    EESC            =   fstatus (   205 :   207 )      
    VUE             =   fstatus (   208 :   210 )      
    C_IC4           =   fstatus (   211 )      
    EEN             =   fstatus (   212 :   214 )      
    SEN             =   fstatus (   215 :   217 )      
    ERC             =   fstatus (   218 :   220 )      
    SRC             =   fstatus (   221 :   223 )      
    EXC             =   fstatus (   224 :   226 )      
    SXC             =   fstatus (   227 :   229 )      
    EBC             =   fstatus (   230 :   232 )      
    SBC             =   fstatus (   233 :   235 )      
    EJ              =   fstatus (   236 :   238 )      
    SJ              =   fstatus (   239 :   241 )      
    EJJ             =   fstatus (   242 :   244 )      
    ERT             =   fstatus (   245 :   247 )      
    SRT             =   fstatus (   248 :   250 )      
    ECP             =   fstatus (   251 :   253 )      
    SCP             =   fstatus (   254 :   256 )      
    ETP             =   fstatus (   257 :   259 )      
    STP             =   fstatus (   260 :   262 )      
    EPC             =   fstatus (   263 :   265 )      
    SPC             =   fstatus (   266 :   268 )      
    EPC1ST          =   fstatus (   269 :   271 )      
    EPT             =   fstatus (   272 :   274 )      
    ETMAX           =   fstatus (   275 :   277 )      
    STMAX           =   fstatus (   278 :   280 )      
    ECCR            =   fstatus (   281 :   283 )      
    SCCR            =   fstatus (   284 :   286 )      
    EUTMAX          =   fstatus (   287 :   289 )      
    EUOVER          =   fstatus (   290 :   292 )      
    E41_51          =   fstatus (   293 :   295 )      
    S41_51          =   fstatus (   296 :   298 )      
    EUOLD           =   fstatus (   299 :   304 )      
    EU_OLD          =   fstatus (   305 :   310 )      
    SNOLD           =   fstatus (   311 :   316 )      
    SNCROLD         =   fstatus (   317 :   322 )      
    EPPC            =   fstatus (   323 :   325 )      
    EPPT            =   fstatus (   326 :   328 )      
    EPEC            =   fstatus (   329 :   331 )      
    CCRACK          =   fstatus (   332 :   334 )      
    ELIMIT          =   fstatus (   335 :   337 )      
    DS1             =   fstatus (   338 :   340 )      
    DE1             =   fstatus (   341 :   343 )      
    DS2             =   fstatus (   344 :   346 )      
    DE2             =   fstatus (   347 :   349 )      
    FS1             =   fstatus (   350 :   352 )      
    FE1             =   fstatus (   353 :   355 )      
    FS2             =   fstatus (   356 :   358 )      
    FE2             =   fstatus (   359 :   361 )      
    RS1             =   fstatus (   362 :   364 )      
    RE1             =   fstatus (   365 :   367 )      
    RS2             =   fstatus (   368 :   370 )      
    RE2             =   fstatus (   371 :   373 )      
    CYCN(1,:)       =   fstatus (   374 :   422 )       
    CYCN(2,:)       =   fstatus (   423 :   471 )       
    CYCN(3,:)       =   fstatus (   472 :   520 )       
    EUDL            =   fstatus (   521 :   523 )      
    GBETA           =   fstatus (   524 )               
                    
                    
                    
    IULOUT          =   istatus (   1   :   3   )       
    KULOUT          =   istatus (   4   :   6   )      
    ICOORDCHGED     =   istatus (   7   )      
    IACT            =   istatus (   8   )      
    IIKEEPSC        =   istatus (   9   )      
    IUL             =   istatus (   10  :   12  )       
    INCF            =   istatus (   13  )       
    IROT            =   istatus (   14  )       
    ICONCYC         =   istatus (   15  )       
    ICIC            =   istatus (   16  :   23  )       
    ICR             =   istatus (   24  )       
    IMM             =   istatus (   25  )       
    NIT             =   istatus (   26  )       
    IR              =   istatus (   27  )       
    NST             =   istatus (   28  )       
    ICOUNTLOOPT     =   istatus (   29  :   31  )       
    ID_CRACK        =   istatus (   32  :   34  )       
    IVIRGIN         =   istatus (   35  :   37  )       
    NTCRACK         =   istatus (   38  )       
    LLCRACK         =   istatus (   39  :   41  )       
    NTFAIL          =   istatus (   42  )       
    KUL             =   istatus (   43  :   45  )       
    KULIN1          =   istatus (   46  :   48  )       
    KULIN5          =   istatus (   49  :   51  )       
    M               =   istatus (   52  )       
    LX              =   istatus (   53  )       
    LY              =   istatus (   54  )       
    LZ              =   istatus (   55  )       
    Mhj             =   istatus (   56  )       
    NG3C            =   istatus (   57  )       
            
    end subroutine 

!****************************
!write the values the values from fistr variables
!****************************   
    
    subroutine write_fistr(istatus,fstatus,                                                 &
            TSI,TEP,TEPLL,SHTEP,SHTEPL,EU,SN,SNCR,ET,SCOLD,EUOUT,SNCROUT,IULOUT,KULOUT,     &
            PR23,PR31,PR12,RG23,RG31,RG12,DCC,DCCL,DCCR,DCCRL,DCCRCHG,DCCRCHGL,DC_1,DC_2,   &
            ICOORDCHGED,IACT,ECU,FC,FT,SCT,EO,EBU,VS,SENS,G0,EC_100,SC,EC,IIKEEPSC,         &
            EPCU,EPCUS,EPCU_,EPCUS_,SC_3D,SC_RAM,EESC,VUE,IUL,INCF,IROT,ICONCYC,            &
            ICIC,C_IC4,ICR,IMM,NIT,IR,NST,EEN,SEN,ERC,SRC,EXC,SXC,EBC,SBC,EJ,SJ,EJJ,        &
            ERT,SRT,ECP,SCP,ETP,STP,EPC,SPC,EPC1ST,EPT,ETMAX,STMAX,ECCR,SCCR,EUTMAX,EUOVER, &                                       
            ICOUNTLOOPT,E41_51,S41_51,ID_CRACK,EUOLD,EU_OLD,SNOLD,SNCROLD,EPPC,EPPT,EPEC,   &                                       
            CCRACK,IVIRGIN,ELIMIT,NTCRACK,LLCRACK,NTFAIL,KUL,KULIN1,KULIN5,DS1,DE1,DS2,DE2, &                                      
            FS1,FE1,FS2,FE2,RS1,RE1,RS2,RE2,M,LX,LY,LZ,CYCN,EUDL,Mhj,NG3C,GBETA)
    
        real (kreal) :: fstatus(:)
        integer :: istatus(:)
        
        real(kreal) TSI(6),TEP(6),TEPLL(6),TEPC(6),TEPL(6)
        real(kreal) SHTEP(3),SHTEPL(3)
        real(kreal) EU(6),SN(6),SNCR(6)
        real(kreal) ET(3),SCOLD(3)
        real(kreal) EUOUT(6),SNCROUT(6)
        integer :: IULOUT(3),KULOUT(3)
        real(kreal) PR12,PR23,PR31,RG12,RG23,RG31
        real(kreal) DCC(3,3),DCCL(3,3),DCCR(3,3),DCCRL(3,3)
        real(kreal) DCCRCHG(3,3,3),DCCRCHGL(3,3,3)
        real(kreal) DC_1(3,3),DC_2(3,3)
        integer ICOORDCHGED,IACT
        real(kreal) ECU,FC,FT,SCT,EO,EBU,VS,SENS,G0,EC_100
        real(kreal) SC(3),EC(3)
        integer IIKEEPSC
        real(kreal) EPCU,EPCUS,EPCU_(3),EPCUS_(3)
        real(kreal) SC_3D(3),SC_RAM(3)
        real(kreal) EESC(3),VUE(3)
        integer IUL(3)
        integer INCF
        integer IROT,ICONCYC  
        integer ICIC(8)
        real(kreal) C_IC4
        integer ICR,IMM,NIT,IR,NST
        real(kreal) EEN(3),SEN(3),ERC(3),SRC(3),EXC(3),SXC(3)
        real(kreal) EBC(3),SBC(3),EJ(3),SJ(3),EJJ(3)
        real(kreal) ERT(3),SRT(3),ECP(3),SCP(3)
        real(kreal) ETP(3),STP(3),EPC(3),SPC(3),EPC1ST(3)
        real(kreal) EPT(3),ETMAX(3),STMAX(3)
        real(kreal) ECCR(3),SCCR(3)
        real(kreal) EUTMAX(3),EUOVER(3)
        integer ICOUNTLOOPT(3)
        real(kreal) E41_51(3),S41_51(3)
        integer ID_CRACK(3)
        real(kreal) EUOLD(6),EU_OLD(6),SNOLD(6),SNCROLD(6)
        real(kreal) EPPC(3),EPPT(3),EPEC(3),CCRACK(3)
        integer IVIRGIN(3)
        real(kreal) ELIMIT(3)
        integer NTCRACK
        integer LLCRACK(3)
        integer NTFAIL
        integer :: KUL(3),KULIN1(3),KULIN5(3)
        real(kreal) DS1(3),DE1(3),DS2(3),DE2(3)
        real(kreal) FS1(3),FE1(3),FS2(3),FE2(3)
        real(kreal) RS1(3),RE1(3),RS2(3),RE2(3)
        integer M,LX,LY,LZ
        real(kreal) CYCN(3,49)
        real(kreal) EUDL(3)
        integer Mhj,NG3C
        real (kreal) :: GBETA       
        
        
        !FLOAT VALUES
    
        fstatus (   1   :   6   )   =   TSI 
        fstatus (   7   :   12  )   =   TEP 
        fstatus (   13  :   18  )   =   TEPLL
        fstatus (   19  :   21  )   =   SHTEP
        fstatus (   22  :   24  )   =   SHTEPL
        fstatus (   25  :   30  )   =   EU
        fstatus (   31  :   36  )   =   SN
        fstatus (   37  :   42  )   =   SNCR
        fstatus (   43  :   45  )   =   ET
        fstatus (   46  :   48  )   =   SCOLD
        fstatus (   49  :   54  )   =   EUOUT
        fstatus (   55  :   60  )   =   SNCROUT
        fstatus (   61  :   61  )   =   PR23
        fstatus (   62  :   62  )   =   PR31
        fstatus (   63  :   63  )   =   PR12
        fstatus (   64  :   64  )   =   RG23
        fstatus (   65  :   65  )   =   RG31
        fstatus (   66  :   66  )   =   RG12
        fstatus (   67  :   69  )   =   DCC(1,:)
        fstatus (   70  :   72  )   =   DCC(2,:)
        fstatus (   73  :   75  )   =   DCC(3,:)
        fstatus (   76  :   78  )   =   DCCL(1,:)
        fstatus (   79  :   81  )   =   DCCL(2,:)
        fstatus (   82  :   84  )   =   DCCL(3,:)
        fstatus (   85  :   87  )   =   DCCR(1,:)
        fstatus (   88  :   90  )   =   DCCR(2,:)
        fstatus (   91  :   93  )   =   DCCR(3,:)
        fstatus (   94  :   96  )   =   DCCRL(1,:)
        fstatus (   97  :   99  )   =   DCCRL(2,:)
        fstatus (   100 :   102 )   =   DCCRL(3,:)
        fstatus (   103 :   105 )   =   DCCRCHG(1,1,:)
        fstatus (   106 :   108 )   =   DCCRCHG(1,2,:)
        fstatus (   109 :   111 )   =   DCCRCHG(1,3,:)
        fstatus (   112 :   114 )   =   DCCRCHG(2,1,:)
        fstatus (   115 :   117 )   =   DCCRCHG(2,2,:)
        fstatus (   118 :   120 )   =   DCCRCHG(2,3,:)
        fstatus (   121 :   123 )   =   DCCRCHG(3,1,:)
        fstatus (   124 :   126 )   =   DCCRCHG(3,2,:)
        fstatus (   127 :   129 )   =   DCCRCHG(3,3,:)
        fstatus (   130 :   132 )   =   DCCRCHGL(1,1,:)
        fstatus (   133 :   135 )   =   DCCRCHGL(1,2,:)
        fstatus (   136 :   138 )   =   DCCRCHGL(1,3,:)
        fstatus (   139 :   141 )   =   DCCRCHGL(2,1,:)
        fstatus (   142 :   144 )   =   DCCRCHGL(2,2,:)
        fstatus (   145 :   147 )   =   DCCRCHGL(2,3,:)
        fstatus (   148 :   150 )   =   DCCRCHGL(3,1,:)
        fstatus (   151 :   153 )   =   DCCRCHGL(3,2,:)
        fstatus (   154 :   156 )   =   DCCRCHGL(3,3,:)
        fstatus (   157 :   159 )   =   DC_1(1,:)
        fstatus (   160 :   162 )   =   DC_1(2,:)
        fstatus (   163 :   165 )   =   DC_1(3,:)
        fstatus (   166 :   168 )   =   DC_2(1,:)
        fstatus (   169 :   171 )   =   DC_2(2,:)
        fstatus (   172 :   174 )   =   DC_2(3,:)
        fstatus (   175 :   175 )   =   ECU
        fstatus (   176 :   176 )   =   FC
        fstatus (   177 :   177 )   =   FT
        fstatus (   178 :   178 )   =   SCT
        fstatus (   179 :   179 )   =   EO
        fstatus (   180 :   180 )   =   EBU
        fstatus (   181 :   181 )   =   VS
        fstatus (   182 :   182 )   =   SENS
        fstatus (   183 :   183 )   =   G0
        fstatus (   184 :   184 )   =   EC_100
        fstatus (   185 :   187 )   =   SC
        fstatus (   188 :   190 )   =   EC
        fstatus (   191 :   191 )   =   EPCU
        fstatus (   192 :   192 )   =   EPCUS
        fstatus (   193 :   195 )   =   EPCU_
        fstatus (   196 :   198 )   =   EPCUS_
        fstatus (   199 :   201 )   =   SC_3D
        fstatus (   202 :   204 )   =   SC_RAM
        fstatus (   205 :   207 )   =   EESC
        fstatus (   208 :   210 )   =   VUE
        fstatus (   211 :   211 )   =   C_IC4
        fstatus (   212 :   214 )   =   EEN
        fstatus (   215 :   217 )   =   SEN
        fstatus (   218 :   220 )   =   ERC
        fstatus (   221 :   223 )   =   SRC
        fstatus (   224 :   226 )   =   EXC
        fstatus (   227 :   229 )   =   SXC
        fstatus (   230 :   232 )   =   EBC
        fstatus (   233 :   235 )   =   SBC
        fstatus (   236 :   238 )   =   EJ
        fstatus (   239 :   241 )   =   SJ
        fstatus (   242 :   244 )   =   EJJ
        fstatus (   245 :   247 )   =   ERT
        fstatus (   248 :   250 )   =   SRT
        fstatus (   251 :   253 )   =   ECP
        fstatus (   254 :   256 )   =   SCP
        fstatus (   257 :   259 )   =   ETP
        fstatus (   260 :   262 )   =   STP
        fstatus (   263 :   265 )   =   EPC
        fstatus (   266 :   268 )   =   SPC
        fstatus (   269 :   271 )   =   EPC1ST
        fstatus (   272 :   274 )   =   EPT
        fstatus (   275 :   277 )   =   ETMAX
        fstatus (   278 :   280 )   =   STMAX
        fstatus (   281 :   283 )   =   ECCR
        fstatus (   284 :   286 )   =   SCCR
        fstatus (   287 :   289 )   =   EUTMAX
        fstatus (   290 :   292 )   =   EUOVER
        fstatus (   293 :   295 )   =   E41_51
        fstatus (   296 :   298 )   =   S41_51
        fstatus (   299 :   304 )   =   EUOLD
        fstatus (   305 :   310 )   =   EU_OLD
        fstatus (   311 :   316 )   =   SNOLD
        fstatus (   317 :   322 )   =   SNCROLD
        fstatus (   323 :   325 )   =   EPPC
        fstatus (   326 :   328 )   =   EPPT
        fstatus (   329 :   331 )   =   EPEC
        fstatus (   332 :   334 )   =   CCRACK
        fstatus (   335 :   337 )   =   ELIMIT
        fstatus (   338 :   340 )   =   DS1
        fstatus (   341 :   343 )   =   DE1
        fstatus (   344 :   346 )   =   DS2
        fstatus (   347 :   349 )   =   DE2
        fstatus (   350 :   352 )   =   FS1
        fstatus (   353 :   355 )   =   FE1
        fstatus (   356 :   358 )   =   FS2
        fstatus (   359 :   361 )   =   FE2
        fstatus (   362 :   364 )   =   RS1
        fstatus (   365 :   367 )   =   RE1
        fstatus (   368 :   370 )   =   RS2
        fstatus (   371 :   373 )   =   RE2
        fstatus (   374 :   422 )   =   CYCN(1,:)
        fstatus (   423 :   471 )   =   CYCN(2,:)
        fstatus (   472 :   520 )   =   CYCN(3,:)
        fstatus (   521 :   523 )   =   EUDL
        fstatus (   524 :   524 )   =   GBETA
        
        
        
        
        istatus (   1   :   3   )   =   IULOUT
        istatus (   4   :   6   )   =   KULOUT
        istatus (   7   :   7   )   =   ICOORDCHGED
        istatus (   8   :   8   )   =   IACT
        istatus (   9   :   9   )   =   IIKEEPSC
        istatus (   10  :   12  )   =   IUL
        istatus (   13  :   13  )   =   INCF
        istatus (   14  :   14  )   =   IROT
        istatus (   15  :   15  )   =   ICONCYC
        istatus (   16  :   23  )   =   ICIC
        istatus (   24  :   24  )   =   ICR
        istatus (   25  :   25  )   =   IMM
        istatus (   26  :   26  )   =   NIT
        istatus (   27  :   27  )   =   IR
        istatus (   28  :   28  )   =   NST
        istatus (   29  :   31  )   =   ICOUNTLOOPT
        istatus (   32  :   34  )   =   ID_CRACK
        istatus (   35  :   37  )   =   IVIRGIN
        istatus (   38  :   38  )   =   NTCRACK
        istatus (   39  :   41  )   =   LLCRACK
        istatus (   42  :   42  )   =   NTFAIL
        istatus (   43  :   45  )   =   KUL
        istatus (   46  :   48  )   =   KULIN1
        istatus (   49  :   51  )   =   KULIN5
        istatus (   52  :   52  )   =   M
        istatus (   53  :   53  )   =   LX
        istatus (   54  :   54  )   =   LY
        istatus (   55  :   55  )   =   LZ
        istatus (   56  :   56  )   =   Mhj
        istatus (   57  :   57  )   =   NG3C
            
    
    
    end subroutine  


    !************************************************************************************** 
        SUBROUTINE CONCRT (                                                     &
            TSI,TEP,TEPLL,SHTEP,SHTEPL,                                         &
            EU,SN,SNCR,ET,SCOLD,                                                &
            EUOUT,SNCROUT,IULOUT,KULOUT,                                        &
            PR23,PR31,PR12,                                                     &
            RG23,RG31,RG12,                                                     &
            DCC,DCCL,DCCR,DCCRL,DCCRCHG,DCCRCHGL,                               &
            DC_1,DC_2,                                                          &
            ICOORDCHGED,IACT,                                                   &
            ECU,FC,FT,SCT,EO,EBU,VS,SENS,G0,EC_100,                             &
            SC,EC,IIKEEPSC,EPCU,EPCUS,EPCU_,EPCUS_,                             &
            SC_3D,SC_RAM,                                                       &
            EESC,VUE,IUL,INCF,IROT,ICONCYC,                                     &
            ICIC,C_IC4,ICR,IMM,NIT,IR,NST,                                      &
            EEN,SEN,ERC,SRC,EXC,SXC,                                            &
            EBC,SBC,EJ,SJ,EJJ,                                                  &
            ERT,SRT,ECP,SCP,                                                    &
            ETP,STP,EPC,SPC,EPC1ST,                                             &
            EPT,ETMAX,STMAX,ECCR,SCCR,                                          &
            EUTMAX,EUOVER,ICOUNTLOOPT,                                          &
            E41_51,S41_51,ID_CRACK,EUOLD,EU_OLD,SNOLD,SNCROLD,                  &
            EPPC,EPPT,EPEC,CCRACK,                                              &
            IVIRGIN,ELIMIT,NTCRACK,LLCRACK,NTFAIL,                              &
            KUL,KULIN1,KULIN5,                                                  &
            DS1,DE1,DS2,DE2,                                                    &
            FS1,FE1,FS2,FE2,                                                    &
            RS1,RE1,RS2,RE2,                                                    &
            M,LX,LY,LZ,CYCN,EUDL,Mhj,NG3C) 
        !**************************************************************************************
        ! *** STRESS - STRAIN RELATIONSHIP FOR CONCRETE ***
        !**************************************************************************************
    
        !add
            real(kreal) TSI(6),TEP(6),TEPLL(6)
            real(kreal) SHTEP(3),SHTEPL(3)
            real(kreal) EU(6),SN(6),SNCR(6)
            real(kreal) ET(3),SCOLD(3)
            real(kreal) EUOUT(6),SNCROUT(6)
            integer :: IULOUT(3),KULOUT(3)
            real(kreal) PR12,PR23,PR31,RG12,RG23,RG31
            real(kreal) DCC(3,3),DCCL(3,3),DCCR(3,3),DCCRL(3,3)
            real(kreal) DCCRCHG(3,3,3),DCCRCHGL(3,3,3)
            real(kreal) DC_1(3,3),DC_2(3,3)
            integer ICOORDCHGED,IACT
            real(kreal) ECU,FC,FT,SCT,EO,EBU,VS,SENS,G0,EC_100
            real(kreal) SC(3),EC(3)
            integer IIKEEPSC
            real(kreal) EPCU,EPCUS,EPCU_(3),EPCUS_(3)
            real(kreal) SC_3D(3),SC_RAM(3)
            real(kreal) EESC(3),VUE(3)
            integer IUL(3)
            integer INCF
            integer IROT,ICONCYC  
            integer ICIC(8)
            real(kreal) C_IC4
            integer ICR,IMM,NIT,IR,NST
            real(kreal) EEN(3),SEN(3),ERC(3),SRC(3),EXC(3),SXC(3)
            real(kreal) EBC(3),SBC(3),EJ(3),SJ(3),EJJ(3)
            real(kreal) ERT(3),SRT(3),ECP(3),SCP(3)
            real(kreal) ETP(3),STP(3),EPC(3),SPC(3),EPC1ST(3)
            real(kreal) EPT(3),ETMAX(3),STMAX(3)
            real(kreal) ECCR(3),SCCR(3)
            real(kreal) EUTMAX(3),EUOVER(3)
            integer ICOUNTLOOPT(3)
            real(kreal) E41_51(3),S41_51(3)
            integer ID_CRACK(3)
            real(kreal) EUOLD(6),EU_OLD(6),SNOLD(6),SNCROLD(6)
            real(kreal) EPPC(3),EPPT(3),EPEC(3),CCRACK(3)
            integer IVIRGIN(3)
            real(kreal) ELIMIT(3)
            integer NTCRACK
            integer LLCRACK(3)
            integer NTFAIL
            integer :: KUL(3),KULIN1(3),KULIN5(3)
            real(kreal) DS1(3),DE1(3),DS2(3),DE2(3)
            real(kreal) FS1(3),FE1(3),FS2(3),FE2(3)
            real(kreal) RS1(3),RE1(3),RS2(3),RE2(3)
            integer M,LX,LY,LZ
            real(kreal) CYCN(3,49)
            real(kreal) EUDL(3)
            integer Mhj,NG3C            
            
            ! additional variables
            real(kreal)SCHJ(3),TEPC(6),TEPL(6)

            !*****************NOT KNOWN
            real(kreal) ECR,A               
            real(kreal) BETA(3)
            real(kreal) ECOLD(3)
            real(kreal) EUDM(6),SNOLDDM(6),SNCRL(6),SNL(6)
            real(kreal) EURR9(3)

            real(kreal) S_NaN(99),N_NaN(99)
            CHARACTER*100 M_NaN
            
            integer I,J,K,IACTOLD,intI
            integer IREK,NTELS
            integer LLIUL,NUMCRACK,INTJ
            integer NNO,NNNOOO1,NNNOOO2,NNNOOO3,NNNOOO4,NNNOOO5,MM
            integer NNNOOO6,NNNOOO7,NNNOOO8,NNNEU,N1122,VUE0
            real(kreal)   SCSCEC,SCMIN,ECMIN
            real(kreal)   G12,G23,G31,G0_100,CN,VUE1L,VUE2L,VUE3L
            real(kreal) DEU(6),EUL(6)
        

            DO I=1,3;SCHJ(I)=SN(I);ENDDO

            !edit to be checked
            LLIUL       =   0  
            IREK            =   0
            IACTOLD         =   0
            NTELS           =   0
            NUMCRACK        =   0
            INTJ            =   0
            DEU         =   0.0D0
            EUL         =   0.0D0
            A           =   0.0D0
            BETA        =   0.0D0
            ECOLD       =   0.0D0
            EUDM        =   0.0D0
            SNOLDDM     =   0.0D0
            SNCRL       =   0.0D0
            SNL         =   0.0D0
            EURR9       =   0.0D0
            SCSCEC      =   0.0D0
            SCMIN       =   0.0D0
            ECMIN       =   0.0D0
            G12         =   0.0D0
            G23         =   0.0D0
            G31         =   0.0D0
            G0_100      =   0.0D0
            CN          =   0.0D0
            VUE1L       =   0.0D0
            VUE2L       =   0.0D0
            VUE3L       =   0.0D0
            
            

        !*** PARAMETERS OF FAILURE SURFACE ******************************************
        
            ECR=FT/EO
        
            DO I=1,6;EUDM(I)=EU(I);ENDDO    
            
            IACTOLD=IACT

 
            CALL sTCracks( DCC,DCCL,DCCR,DCCRL,DCCRCHG,DCCRCHGL,                        &
                           PR23,PR31,PR12,RG23,RG31,RG12,                               &
                           TSI,SN,SNCR,TEP,TEPLL,EUDM,ET,ECR,ECCR,ECU,EEN,              &
                           IUL,CCRACK,ICOORDCHGED,IACT,                                 &
                           M,I,J,K,NST,NIT,IROT,ICIC(6)                                 )
            

            CALL sGetIREK(ICIC(6),DCC,DCCL,DCCR,DCCRL,NST,IACT,IACTOLD,IREK)
            
                    
            CALL sChangeAxisTurn(DCC,DCCR,IACT,SNCR,IREK)
            
                       
            CALL CRACK_1(ICIC(6),FT,EO,EU,ETP,DC_1,DC_2,DCC,DCCR)          
            
                DO I=1,3
                    SCOLD(I)=SC(I)
                    ECOLD(I)=EC(I)
                ENDDO

        ! *** CONSTITUTIVE MODEL
        !
        !     <IROT=0>          : EQUIVALENT UNIAXIAL STRAIN MODEL
        !     <IROT=1 & ICR=0>  : EQUIVALENT UNIAXIAL STRAIN MODEL
        !     <IROT=1 & ICR=1>  : CRACK DIRECTION CORDINATE MODEL
        !
        !   CALL sCheckCMODEL(IROT,CCRACK,EU,ECR,ECCR,
        !     * ECU,LLCRACK,NTCRACK,NTFAIL,NTELS,LLIUL,NUMCRACK)
        !


            IF(ICIC(6).EQ.0)    GO TO 1200        !新固定一方向�Eび割れモチE��
            IF(IACT.NE.0)           GO TO 200
        !
        ! *** 1 : EQUIVALENT UNIAXIAL STRAIN MODEL ***********************
        !
        ! *** CALCULATE EQIVALENT UNIAXIAL STRAINS4
        !
    1200    continue
        !   ***主軸でのひずみの計箁E**
            CALL sEURote(EU,EUL,DEU,TEP,TEPLL,ET,DCC,DCCL, &
                    PR23,PR31,PR12,RG23,RG31,RG12)

                    
      GO TO 300
    
    200     continue
        

        CALL sEUFixed(EU,EUL,DEU,TEP,TEPLL,ET,DCCR,DCCRL, &
                            PR23,PR31,PR12,RG23,RG31,RG12)
                                        

                    
  300 CONTINUE
                    
        CALL sCheckCMODEL(IROT,CCRACK,EU,ECR,ECCR,ECU,EEN, &
            LLCRACK,NTCRACK,NTFAIL,NTELS,LLIUL,NUMCRACK)
            
    
        A=ABS(FT/20.)
        IF(ABS(SN(1))<A.AND.ABS(SN(2))<A.AND.ABS(SN(3))<A)THEN
            DO I=1,3
                SC(I)=FC 
                EC(I)=ECU
            ENDDO
        ELSE
! *** --- *********************************** ICIC(8):3D効极E:なぁE1:あり,2:新ルール

        IF( ICIC(8) == 0 )THEN
                DO I=1,3
                    SC(I)=FC
                    EC(I) = ECU
                ENDDO
            ELSE IF( ICIC(8) == 1 )THEN
                CALL sSCEC(SN,SNCR,EU,IUL,IROT,ICR,FC,ECU,EO,FT,SC,EC,      &
                        EEN,EBU,ICIC(1),ECCR,NST,CCRACK,EPCU,ICONCYC,SCHJ)
            ELSE IF( ICIC(8) == 2 .OR. ICIC(8) == 3 )THEN
                CALL sSCEC_hj(SN,SNCR,EU,IUL,IROT,ICR,FC,ECU,EO,FT,SC,EC,   &
                EEN,EBU,ICIC(1),ECCR,NST,CCRACK,EPCU,ICONCYC,SCHJ,ICIC(8))
            ENDIF
        END IF
                

        DO intI=1,3
            SCscec=1.0
            IF(SC(intI) > SCscec*FC ) THEN
                SC(intI)=SCscec*FC ; EC(intI)=ECU
            ENDIF
            SCscec=1.0
            IF(EC(intI) > SCscec*ECU ) THEN
                EC(intI)=ECU
            ENDIF
        ENDDO
        

        !CH---
    DO I=1,3 ; CYCN(I,35)=SCHJ(I) ; END DO ! BY HJ
    DO I=1,3 ; CYCN(I,36)=EC(I) ; END DO ! BY HJ

!c  IF (NTCRACK.NE.0) THEN
        SCMIN=MIN(SC(1),SC(2),SC(3))
        ECMIN=MIN(EC(1),EC(2),EC(3))
        DO intI=1,3
            IF(SC(intI).LT.FC)THEN
! c             SC(intI)=FC
! c             EC(intI)=ECU
            ENDIF
! c         SC(intI)=SCMIN
! c         EC(intI)=ECMIN
        ENDDO
! c ENDIF
    

    
    DO intI=1,3
        SC_3D(intI)=SC(intI)/FC
    ENDDO

! C     ***COMPRESSIVE REDUCTION FACTOR BY K.YONEZAWA (modified by tide)***
! C ***[圧縮低減係数βの計算]***
      BETA(1)=1.0D0
      BETA(2)=1.0D0
      BETA(3)=1.0D0
    
    
    ! C--- *** --- ***********************************
    IF(     ICIC(3) == 0 ) THEN             ![0:野口濱田式]
        CALL sRAMDA (EU,BETA,SC,ECU,FC,IUL,                 &
                 EURR9,EEN,EBU,SN,ECR,ECCR,NST,CCRACK)
    
      ELSE IF(ICIC(3) == 2 ) THEN             ![2:低減なしβ:1.0]
        GO TO 3110
      ELSE IF(ICIC(3) == 3 ) THEN             ![3:新野口濱田式]
         CALL sRAMDA_hj (EU,SC,EC,EBU,ETP,SC_RAM)
    ENDIF
    

 3110 CONTINUE
    
    ! CH---エラーチェチE��---
! C DO I=1,3
! C     IF (EEN(I)>=0.)THEN
! C         WRITE(1850,'(A50,I5,ES11.3)') 
! C     +           '●sRAMDA: EEN≧0.0 NST NIT LX LY LZ ',I,EEN(I)
! C     ENDIF
! C ENDDO
! C M_NaN='☁ERAMDA: [ EEN 1 2 3 ]  NST  NIT  LX  LY  LZ '
! C N_NaN(1) = NST ; N_NaN(2) = NIT
! C N_NaN(3) = LX ; N_NaN(4) = LY ; N_NaN(5) = LZ
! C S_NaN(1) = EEN(1) ; S_NaN(2) = EEN(2) ; S_NaN(3) = EEN(3)
! C CALL CCHJ(N_NaN, 5, S_NaN, 3 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
    DO I=1,3 ; CYCN(I,37)=EURR9(I) ; END DO ! BY HJ

    CALL CHG3(SNOLDDM(1),SNOLDDM(2),SNOLDDM(3),     &
        SNOLDDM(4),SNOLDDM(5),SNOLDDM(6),           &
        DCCRL(1,1),DCCRL(1,2),DCCRL(1,3),           &
        DCCRL(2,1),DCCRL(2,2),DCCRL(2,3),           &
        DCCRL(3,1),DCCRL(3,2),DCCRL(3,3),           &
        SNOLD(1),SNOLD(2),SNOLD(3),                 &
        SNOLD(4),SNOLD(5),SNOLD(6),0)

      CALL CHG3(SNL(1),SNL(2),SNL(3),               &
        SNL(4),SNL(5),SNL(6),                       &
        DCCR(1,1),DCCR(1,2),DCCR(1,3),              &
        DCCR(2,1),DCCR(2,2),DCCR(2,3),              &
        DCCR(3,1),DCCR(3,2),DCCR(3,3),              &
        SNOLDDM(1),SNOLDDM(2),SNOLDDM(3),           &
        SNOLDDM(4),SNOLDDM(5),SNOLDDM(6),1)
    
    CALL CHG3(SNOLDDM(1),SNOLDDM(2),SNOLDDM(3),     &
        SNOLDDM(4),SNOLDDM(5),SNOLDDM(6),           &
        DCCRL(1,1),DCCRL(1,2),DCCRL(1,3),           &
        DCCRL(2,1),DCCRL(2,2),DCCRL(2,3),           &
        DCCRL(3,1),DCCRL(3,2),DCCRL(3,3),           &
        SNCROLD(1),SNCROLD(2),SNCROLD(3),           &
        SNCROLD(4),SNCROLD(5),SNCROLD(6),0)

      CALL CHG3(SNCRL(1),SNCRL(2),SNCRL(3),         &
        SNCRL(4),SNCRL(5),SNCRL(6),                 &
        DCCR(1,1),DCCR(1,2),DCCR(1,3),              &
        DCCR(2,1),DCCR(2,2),DCCR(2,3),              &
        DCCR(3,1),DCCR(3,2),DCCR(3,3),              &
        SNOLDDM(1),SNOLDDM(2),SNOLDDM(3),           &
        SNOLDDM(4),SNOLDDM(5),SNOLDDM(6),1)

        
        
        
    IF(ICIC(6).EQ.0)    GOTO 3100
    IF(NTCRACK.EQ.0)    GOTO 3100
      GO TO 3200
    
    
    
    
    ! C     <IROT=0>                  : EQUIVALENT UNIAXIAL STRAIN MODEL
! C     <IROT=1 & ICR =0>         : EQUIVALENT UNIAXIAL STRAIN MODEL
! C
! C     *** NO CRACK ***
! C     ***[ひび割れなぁE�E��E輩の言ぁE��がどぁE��な�E�E]***
! C ***[固定�Eび割れモチE��を用ぁE��めEひび割れてぁE��ぁE��合では,回転ひび割れモチE��と見なして処琁E��る]***
 3100 CONTINUE
    
    ! C ***[回転ひび割れモチE��で応力と剛性を求める]***
        do intI=1,3
          IF (EC(intI).gt.0.0) THEN 
            CONTINUE
        ENDIF
        enddo


    CALL sSressStiff(NST,NIT,BETA,ICONCYC,              &
        DEU,EU,SN,                                      &
        ERC,SRC,EXC,SXC,                                &
        EBC,SBC,EJ,SJ,EJJ,                              &
        ERT,SRT,ECP,SCP,                                &
        ETP, STP,EPC,SPC,                               &
        EPC1ST,EPT,EEN,SEN,                             &
        ETMAX,STMAX,ECCR,SCCR,                          &
        EUTMAX,EUOVER,ICOUNTLOOPT,                      &
        E41_51,S41_51,ID_CRACK,                         &
        EUOLD,EU_OLD,SNOLD,EPCU,EPCUS,EPCU_,EPCUS_,     &
!C      EUL,SNL,EPCU,EPCUS,                             &
        EO,ET, EC,SC,FC,FT,EBU,ECU,EC_100,              &
        SCOLD,ECOLD,                                    &
        IUL,EPPC,EPPT,EPEC,CCRACK,                      &
        INCF,ICIC,C_IC4,SENS,VS,                        &
        IVIRGIN,ELIMIT,MM,IROT,ICR,0,CYCN,EUDL,Mhj)

        
    

        
        
        
    DO I=1,3
        SNCR(I)=SN(I)
    ENDDO
        
        ! CH---エラーチェチE��---
! C DO I=1,3
! C     IF (EEN(I)>=0.)THEN
! C         WRITE(1850,'(A50,I5,ES11.3)') 
! C     +           '●sSressStiff_1: EEN≧0.0 NST NIT LX LY LZ ',I,EEN(I)
! C     ENDIF
! C ENDDO
! C M_NaN='☁ESressStiff_1: [ EEN 1 2 3 ]  NST  NIT  LX  LY  LZ '
! C N_NaN(1) = NST ; N_NaN(2) = NIT
! C N_NaN(3) = LX ; N_NaN(4) = LY ; N_NaN(5) = LZ
! C S_NaN(1) = EEN(1) ; S_NaN(2) = EEN(2) ; S_NaN(3) = EEN(3)
! C CALL CCHJ(N_NaN, 5, S_NaN, 3 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
! C
! C ***[こ�E部刁E�E固定�Eび割れモチE��でひび割れてなぁE��合でもPSICRで出力できるためである]***
      IF(ICR.EQ.1) THEN
        DO intI=1,3
            KUL(intI)=0
        DO intJ=1,3
            DCCR(intI,intJ)=DCC(intI,intJ)
        ENDDO
        ENDDO
      ENDIF
    CALL sTransArrayReal(SN,SNCR,6)
    GOTO 678
        
        
! C ***[固定�Eび割れモチE��(ひび割れた場吁Eで応力と剛性を求める]***
 3200 CONTINUE
    CALL sSressStiff(NST,NIT,BETA,ICONCYC,                  &
        DEU,EU,SNCR,                                        &
        ERC,SRC,EXC,SXC,                                    &
        EBC,SBC,EJ,SJ,EJJ,                                  &
        ERT,SRT,ECP,SCP,                                    &
        ETP, STP,EPC,SPC,                                   &
        EPC1ST,EPT,EEN,SEN,                                 &
        ETMAX,STMAX,ECCR,SCCR,                              &
        EUTMAX,EUOVER,ICOUNTLOOPT,                          &
        E41_51,S41_51,ID_CRACK,                             &
        EUOLD,EU_OLD,SNCROLD,EPCU,EPCUS,EPCU_,EPCUS_,       &
! C     EUL,SNL,EPCU,EPCUS,                                 &
        EO,ET, EC,SC,FC,FT,EBU,ECU,EC_100,                  &
        SCOLD,ECOLD,                                        &
        IUL,EPPC,EPPT,  EPEC,CCRACK,                        &
        INCF,ICIC,C_IC4,SENS,VS,                            &
        IVIRGIN,ELIMIT,MM,IROT,ICR,1,CYCN,EUDL,Mhj)
        DO I=1,3;SN(I)=SNCR(I);ENDDO

! CH---エラーチェチE��---
! C DO I=1,3
! C     IF (EEN(I)>=0.)THEN
! C         WRITE(1850,'(A50,I5,ES11.3)') 
! C     +           '●sSressStiff_2: EEN≧0.0 NST NIT LX LY LZ ',I,EEN(I)
! C     ENDIF
! C ENDDO
! C M_NaN='☁ESressStiff_2: [ EEN 1 2 3 ]  NST  NIT  LX  LY  LZ '
! C N_NaN(1) = NST ; N_NaN(2) = NIT
! C N_NaN(3) = LX ; N_NaN(4) = LY ; N_NaN(5) = LZ
! C S_NaN(1) = EEN(1) ; S_NaN(2) = EEN(2) ; S_NaN(3) = EEN(3)
! C CALL CCHJ(N_NaN, 5, S_NaN, 3 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
  678 CONTINUE

        
        

        

  
    IF(NTCRACK.EQ.0)    GOTO 4000
    IF(ICIC(6).EQ.0)    GOTO 3300
! C
! C     *** SHEAR MODULUS : THE CONTACT DENSITY MODEL ***
! C     *** REFERENCE     : 鉁E��コンクリート�E非線形解析と構�E剁E��技報堂�E版！E***
! C
! C     G  : ELASTIC SHEAR MODULUS
! C     RG12 : REDUCED SHEAR MODULUS
! C     RG22 : REDUCED SHEAR MODULUS
! C     RG31 : REDUCED SHEAR MODULUS

    IF(ICIC(5).EQ.1)THEN
    
    G23 =0.5D0*EO/(1.0D0+0.2D0)
    G31 =0.5D0*EO/(1.0D0+0.2D0)
    G12 =0.5D0*EO/(1.0D0+0.2D0)

    CALL sTransRToA(EU(4),EU(5),EU(6),SHTEP,3)
    CALL sTransRToA(EUOLD(4),EUOLD(5),EUOLD(6),SHTEPL,3)

    
    CALL SHEARREDUCE(NST,G0,NG3C,KUL,KULIN1,KULIN5,         &
    SNCR,SNOLD,EU,EUOLD,ECCR,CCRACK,LLCRACK,NTCRACK,        &
! C  SNCR,SNOLD,EU,EUL,ECCR,CCRACK,LLCRACK,                 &
    FC,EO,FT,                                               &
    DS1,DE1,DS2,DE2,                                        &
    FS1,FE1,FS2,FE2,                                        &
    RS1,RE1,RS2,RE2,                                        &
    RG12,RG31,RG23)                                         

    G0_100 = G0 / 100.
    IF( RG12 < G0_100 ) RG12 = G0_100
    IF( RG23 < G0_100 ) RG23 = G0_100
    IF( RG31 < G0_100 ) RG31 = G0_100

! C
! CH---Gij 出力用 ------------------------------------
    !Mhj = XHJ(1)
    IF(M==Mhj) THEN
        NNO=6000+LX*100+LY*10+LZ
        CN=NST+NIT/10.
        IF( N1122<9.AND.NNNOOO1/=NNO.AND.NNNOOO2/=NNO.AND.NNNOOO3/=NNO      &
                    .AND.NNNOOO4/=NNO.AND.NNNOOO5/=NNO.AND.NNNOOO6/=NNO &
                    .AND.NNNOOO7/=NNO.AND.NNNOOO8/=NNO) THEN                
            ! WRITE(NNO,*) M80
            ! WRITE(NNO,'((A5,I5),1A1,3I1)')
     ! +            'MEL=' , M ,O,LX,LY,LZ
            WRITE(*,*) 'STEP,RG12,RG31,RG23,γ12,γ23,γ31,    &
     NCR,LCR1,LCR2,LCR3'
            N1122=N1122+1
            IF     (N1122==1)THEN
                NNNOOO1=NNO
            ELSE IF(N1122==2)THEN
                NNNOOO2=NNO
            ELSE IF(N1122==3)THEN
                NNNOOO3=NNO
            ELSE IF(N1122==4)THEN
                NNNOOO4=NNO
            ELSE IF(N1122==5)THEN
                NNNOOO5=NNO
            ELSE IF(N1122==6)THEN
                NNNOOO6=NNO
            ELSE IF(N1122==7)THEN
                NNNOOO7=NNO
            ELSE IF(N1122==8)THEN
                NNNOOO8=NNO
            ENDIF
        ENDIF
        ! WRITE(NNO,'(F10.1,A1,3(F10.3,A1),3(F10.0,A1),4(I2,A1))') &
                    ! CN,O,RG12/G0,O,RG31/G0,O,RG23/G0,O,            &
                    ! ((EU(J)*1000000.,O),J=4,6),                    &
                    ! NTCRACK,O,((LLCRACK(J),O),J=1,3)
    END IF
        GO TO 4000
    END IF

! C     <IROT=1 & ICR=1> : CRACK DIRECTION CORDINATE MODEL
! C      *** CRACK ***
    CONTINUE

! C     *** SHEAR MODULUS : AL-MAHAIDI MODEL ***

 3300 G23 =0.5D0*EO/(1.0D0+0.2D0)
      G31 =0.5D0*EO/(1.0D0+0.2D0)
      G12 =0.5D0*EO/(1.0D0+0.2D0)
      ECR=FT/EO
    RG12=G12
    RG23=G23
    RG31=G31
! C     *** SHEAR MODULUS : AL-MAHAIDI MODEL ***
    IF(ICONCYC.EQ.1)THEN
! C GOTO 988

! CH---修正1-----------------------------------------------------
      DO I=1,6
       NNNEU=EU(I)*10000000000.
       IF(EU(I)==0.)THEN
! c    IF(NNNEU==0)THEN
! c       EU(I)=0.00000001
       END IF
    END DO
! CH---修正1-----------------------------------------------------
        IF (LLCRACK(1).EQ.1.AND.CCRACK(1).EQ.4) THEN
            IF(EU(1)==0.)write(1861,*)'☁E��正1: EU(1)=',EU(1)
            RG12=G12*0.4D0/(EU(1)/ECR)
            RG31=G31*0.4D0/(EU(1)/ECR)
        ENDIF
        IF (LLCRACK(2).EQ.1.AND.CCRACK(1).EQ.4) THEN
            IF(EU(2)==0.)write(1861,*)'☁E��正1: EU(2)=',EU(2)
            RG12=G12*0.4D0/(EU(2)/ECR)
            RG23=G23*0.4D0/(EU(2)/ECR)
        ENDIF
        IF (LLCRACK(3).EQ.1.AND.CCRACK(1).EQ.4) THEN
            IF(EU(3)==0.)write(1861,*)'☁E��正1: EU(3)=',EU(3)
            RG31=G31*0.4D0/(EU(3)/ECR)
            RG23=G23*0.4D0/(EU(3)/ECR)
        ENDIF
        IF (LLCRACK(1).EQ.1.AND.CCRACK(1).EQ.5) THEN
            RG12=G12*0.4D0/((EU(1)-ECCR(1)+ECR)/ECR)
            RG31=G31*0.4D0/((EU(1)-ECCR(1)+ECR)/ECR)
        ENDIF
        IF (LLCRACK(2).EQ.1.AND.CCRACK(1).EQ.5) THEN
            RG12=G12*0.4D0/((EU(2)-ECCR(2)+ECR)/ECR)
            RG23=G23*0.4D0/((EU(2)-ECCR(2)+ECR)/ECR)
        ENDIF
        IF (LLCRACK(3).EQ.1.AND.CCRACK(1).EQ.5) THEN
            RG31=G31*0.4D0/((EU(3)-ECCR(3)+ECR)/ECR)
            RG23=G23*0.4D0/((EU(3)-ECCR(3)+ECR)/ECR)
        ENDIF
    GOTO 988
! C  988    CONTINUE
        IF (CCRACK(1).EQ.4.and.EU(1).GT.ECR) RG23=G23*0.4D0/(EU(1)/ECR)
        IF (CCRACK(2).EQ.4.and.EU(2).GT.ECR) RG31=G31*0.4D0/(EU(2)/ECR)
        IF (CCRACK(3).EQ.4.and.EU(3).GT.ECR) RG12=G12*0.4D0/(EU(3)/ECR)
        IF (CCRACK(1).EQ.5.and.EU(1).GT.ECCR(1))THEN
            RG23=G23*0.4D0/((EU(1)-ECCR(1)+ECR)/ECR)
        ENDIF
        IF (CCRACK(2).EQ.5.and.EU(2).GT.ECCR(2))THEN
            RG31=G31*0.4D0/((EU(2)-ECCR(2)+ECR)/ECR)
        ENDIF
        IF (CCRACK(3).EQ.5.and.EU(3).GT.ECCR(3))THEN
            RG12=G12*0.4D0/((EU(3)-ECCR(3)+ECR)/ECR)
        ENDIF
  988   CONTINUE
    ENDIF

 4000 CONTINUE

! C *** DECIDE POISSON'S RATIO

      VUE1L=VUE(1)
      VUE2L=VUE(2)
      VUE3L=VUE(3)

      VUE0=0.2D0

! C ***ポアソン比�E初期値を設定すめE**
    PR23=0.02D0
    PR31=0.02D0
    PR12=0.02D0

! C ***[一方向目が�Eび割れた場吁E***
      IF(LLCRACK(1).EQ.1) THEN
       PR12=0.00005D0
       PR31=0.00005D0
      END IF
! C ***[二方向目が�Eび割れた場吁E***
      IF(LLCRACK(2).EQ.1) THEN
       PR23=0.00005D0
       PR12=0.00005D0
      END IF
! C ***[二方向目が�Eび割れた場吁E***
      IF(LLCRACK(3).EQ.1) THEN
       PR23=0.00005D0
       PR31=0.00005D0
      END IF

! C ***非線形ポアソン比を求めめE**
! C CALL SPOISSON(  EU1,EU2,EU3,
! C     *               EC1,EC2,EC3,
! C     *               PR012,PR023,PR031,
! C     *               PR12,PR23,PR31)
! C 
! c WRITE(1820,'(2F5.0,TR10,3F3.0)')
! c     +CYCN(1,1),CYCN(1,49),CYCN(1,31),CYCN(1,32),CYCN(1,33)
! c WRITE(1820,*)'        SN      SNCR        EU'
! c DO intI=1,6
! c write(1820,'(I2,3ES10.2)')inti,SN(intI),SNCR(intI),EU(intI)
! c ENDDO
    IF(ICIC(6).EQ.0)THEN
        DO intI=4,6
            SN(intI)=0.0
            SNCR(intI)=0.0
            EU(intI)=0.0
        ENDDO
    ENDIF
    CALL sPSITurn(EU,SNCR,DCCR,EUOUT,SNCROUT,   &
        IUL,KUL,IULOUT,KULOUT)


        
        
      RETURN
      
                
      END
        

        !**************************************************************************************    
            SUBROUTINE sTCracks( DCC,DCCL,DCCR,DCCRL,DCCRCHG,DCCRCHGL,                  &
                                 PRDM23,PRDM31,PRDM12,RGDM23,RGDM31,RGDM12,                 &
                                 TSI,PSI,PSICR,TEP,TEPLL,EU,ET,ECR,ECCR,ECU,EEN,            &
                                 IUL,CCRACK,ICOORDCHGED,IACTIVEC,                           &
                                 M,I,J,K,NST,NIT,IROT,ICIC6)
                !**************************************************************************************
            
                real(kreal) DAIX(3,3),DAIXL(3,3)
                real(kreal) DCC(3,3),DCCL(3,3),DCCR(3,3),DCCRL(3,3),DCCDM(3,3)
                real(kreal) DCCRCHG(3,3,3),DCCRCHGL(3,3,3)
                real(kreal) PSICR(6),PSI(6),TSI(6),PSIL(6)
                real(kreal) TEPLL(6),TEP(6),EU(6),EULD(6),EUL(6),EEN(3) 
                real(kreal) DEU(6),ET(3),EUPSI(6),EUCR(6)      !sun:DEU(3)
                real(kreal) CCRACK(3)
                integer IUL(3),LCRACK(3)   
                real(kreal) euuu(3),tuuu(3)
                !add 
                real(kreal) ECCR(3),ECR,ECU
                integer NST,M,NIT,I,J,K,NTCRACK,NTFAIL,NTELS,NUMCRACK
                integer IACT,IACTOLD,ICOORDCHGED,IACTIVEC,IACTIVECOLD,IREK
                integer IJNC,IROT
                integer ICIC6
                real(kreal) PR12,PR23,PR31,RG12,RG23,RG31,RGDM12,RGDM23,RGDM31
                real(kreal) PRDM12,PRDM23,PRDM31
                !added by shan CHECK
                IREK=0
                EUPSI=0.0D0
                EUL=0.0D0
                DEU=0.0D0
                EUCR=0.0D0
                LCRACK=0
                NTCRACK=0
                NTFAIL=0
                NTELS=0
                NUMCRACK=0
                IJNC=0
                IACTOLD=IACT
                !********************************************************************************      
                
                PR23=PRDM23
                PR31=PRDM31
                PR12=PRDM12
            
                RG23=RGDM23
                RG31=RGDM31
                RG12=RGDM12
            
                  CALL PRIN (   TEP(1),TEP(2),TEP(3),                           &
                                TEP(4),TEP(5),TEP(6),                           &
                                EUPSI(1),EUPSI(2),EUPSI(3),                     &
                                EUPSI(4),EUPSI(5),EUPSI(6),                     &
                                DCC(1,1),DCC(1,2),DCC(1,3),                     &
                                DCC(2,1),DCC(2,2),DCC(2,3),                     &
                                DCC(3,1),DCC(3,2),DCC(3,3)      )

                IF(ICIC6.EQ.0)THEN
                    IACTIVEC=0
                    RETURN
                ENDIF

                IACTIVECOLD=IACTIVEC
                !********************************************************************************
                !!  ***[主軸でのひずみの計算]***
                !!  ***[主軸でのひずみの計算]***
                CALL sTransArrayreal(EU,EUPSI,6)
                CALL sGetIREK(ICIC6,DCC,DCCL,DCCR,DCCRL,NST,                     &
                                IACT,IACTOLD,IREK)
                CALL sChangeAxisPR(PR23,PR31,PR12,RG23,RG31,RG12,IREK)
            
                CALL sEURote(EUPSI,EUL,DEU,TEP,TEPLL,ET,DCC,DCCL,               & 
                             PR23,PR31,PR12,RG23,RG31,RG12) 
                !********************************************************************************
                !!  ***[ひび割れ�E判断]***
                !********************************************************************************
                CALL sCrackChk(CCRACK,EUPSI,ECR,ECCR,ECU,EEN,                   &
                                 LCRACK,NTCRACK,NTFAIL,NTELS,NUMCRACK)
                
                !********************************************************************************
                !!  ***[ひび割れと判断しなぁE��吁E主軸を用ぁE��]***
                !********************************************************************************
                IF(NTCRACK.EQ.0) THEN
                CALL sTransArrayreal33(DCC,DCCR)
                    IACTIVEC=0
                !!      CALL sTransArrayreal(EUPSI,EU,6)
                !!
                !!********************************************************************************
                !!! ***[主応力状態で新たなひび割れと判断する場吁E***
                !!********************************************************************************
                ELSE
                    SELECT CASE(ICOORDCHGED)
                !!********************************************************************************
                !!!     ***[ひび割れ暦がなぁE��吁E***
                !!********************************************************************************
                    CASE(0)
                !!********************************************************************************
                !!!         ***[一軸固定�Eび割れモチE��の適用]***
                !!********************************************************************************
                        ICOORDCHGED=1
                        IACTIVEC=1
                        CALL sTransArrayreal33(DCC,DCCR)
                        CALL sTransArrayreal33(DCCL,DCCRL)
                !!          CALL sCrack1ST(TSI,TEP,DCCR,PSICR,NTCRACK,LCRACK)
                !!********************************************************************************
                !!!         ***[ひび割れ暦の記録]***
                !!********************************************************************************
                        CALL sRembCrack(DCC,DCCL,DCCRCHG,DCCRCHGL,1)
                !!!
                !!********************************************************************************
                !!!     ***[ひび割れ暦ぁE回ある場吁E***
                !!********************************************************************************
                    CASE(1)
                !!********************************************************************************
                !!!         ***[新たなひび割れ発生�E条件に満たすか�E判断]***
                !!********************************************************************************
                        CALL sTransArrayreal(EU,EUCR,6)
                        CALL sJudgeNewCrack(DCC,DCCRCHG,ICOORDCHGED,        &
                            IJNC,LCRACK,EUPSI,EUCR,IACTIVEC,IACTOLD,        &
                            DEU,TEP,TEPLL,ET,DCCRL,DCCL,                    &
                            IROT,CCRACK,ECR,ECCR,ECU,EEN,                   &
                            PRDM23,PRDM31,PRDM12,RGDM23,RGDM31,RGDM12)
                !!!
                !!********************************************************************************
                !!!         ***[1方向�Eび割れモチE��]***
                !!********************************************************************************
                        IF(ICIC6.LT.2) THEN
                            IJNC=0
                        ENDIF
                !!********************************************************************************
                !!!         ***[ひび割れ暦がある場合�E処琁E��ECOORDCHGEDは�E�）]***
                !!********************************************************************************
                        CALL sCoordCHG1(DCCRCHG,DAIX,DAIXL,DCCR,DCCRL,DCC,DCCL,         &
                                        TSI,TEP,PSICR,M,I,J,K,NTCRACK,LCRACK,                       &
                                        IJNC,IACTIVEC,IACTIVECOLD,ICOORDCHGED)
                !********************************************************************************
                !!!     ***[ひび割れ暦ぁE回ある場吁E***
                !!********************************************************************************
                    CASE(2)
                !!********************************************************************************
                !!!         ***[新たなひび割れ発生�E条件に満たすか�E判断]***
                !!********************************************************************************
            
                        CALL sTransArrayreal(EU,EUCR,6)
                        CALL sJudgeNewCrack(DCC,DCCRCHG,ICOORDCHGED,                    &
                            IJNC,LCRACK,EUPSI,EUCR,IACTIVEC,IACTOLD,                    &
                            DEU,TEP,TEPLL,ET,DCCRL,DCCL,                                &
                            IROT,CCRACK,ECR,ECCR,ECU,EEN,                               &
                            PRDM23,PRDM31,PRDM12,RGDM23,RGDM31,RGDM12)
                !!********************************************************************************
                !!!             ***[ひび割れ暦がある場合�E処琁E��ECOORDCHGEDは2�E�]***
                !!********************************************************************************
                        CALL sCoordCHG2(DCCRCHG,DAIX,DAIXL,DCC,DCCL,DCCR,DCCRL,         &
                            TSI,TEP,PSICR,M,I,J,K,NTCRACK,LCRACK,                       &
                            IJNC,IACTIVEC,IACTIVECOLD,ICOORDCHGED)
                
                !********************************************************************************
                !!      ***[ひび割れ暦ぁE回ある場吁E***
                !********************************************************************************
                    CASE(3)
                !********************************************************************************
                !!          ***[新たなひび割れ発生�E条件に満たすか�E判断(ここではただどれに近いの判断)]***
                !********************************************************************************
                        CALL sTransArrayreal(EU,EUCR,6)
                        CALL sJudgeNewCrack(DCC,DCCRCHG,ICOORDCHGED,                    &
                            IJNC,LCRACK,EUPSI,EUCR,IACTIVEC,IACTOLD,                    &
                            DEU,TEP,TEPLL,ET,DCCRL,DCCL,                                &
                            IROT,CCRACK,ECR,ECCR,ECU,EEN,                               &
                            PRDM23,PRDM31,PRDM12,RGDM23,RGDM31,RGDM12)
                !********************************************************************************
                !!              ***[ひび割れ暦がある場合�E処琁E��ECOORDCHGEDは2�E�]***
                !********************************************************************************
                        CALL sCoordCHG3(DCCRCHG,DAIX,DAIXL,DCC,DCCL,DCCR,DCCRL,         &
                            TSI,TEP,PSICR,M,I,J,K,NTCRACK,LCRACK,                       &
                            IJNC,IACTIVEC,IACTIVECOLD)
                    END SELECT  
                    ENDIF
            END SUBROUTINE
            
        !***********************************************************************
            SUBROUTINE PRIN (SR1,SR2,SR3,                                               &
                                  SR4,SR5,SR6,                     &
                                  PS1,PS2,PS3,                     &
                                  PS4,PS5,PS6,                     &
                                  DC11,DC12,DC13,                  &
                                  DC21,DC22,DC23,                  &
                                  DC31,DC32,DC33)
                !!***********************************************************************
                !! *** CALCULATE PRINCIPAL STRESS, ITS DIRECTION AND MAXMUM SHEAR STRESS
                !!
                !!  !implicit none
                !!
                real(kreal)  A(3,3),V(3,3)
                integer L,M,NR,NR2,IERR
                real(kreal) EPS,AAA,AAA1,AAA2,AAA3,BBB,BBB1,BBB2,BBB3
                real(kreal) CCC,CCC1,CCC2,CCC3
                real(kreal) DC11,DC12,DC13,DC21,DC22,DC23,DC31,DC32,DC33
                real(kreal) PS1,PS2,PS3,PS4,PS5,PS6,SR1,SR2,SR3,SR4,SR5,SR6
            
                DATA L,M,NR,EPS / 3,3,100,1.D-6    /
                !add    

                  A(1,1)=SR1
                  A(1,2)=SR4
                  A(1,3)=SR6
                  A(2,1)=A(1,2)
                  A(2,2)=SR2
                  A(2,3)=SR5
                  A(3,1)=A(1,3)
                  A(3,2)=A(2,3)
                  A(3,3)=SR3
                
                CALL JACOBI (A,V,EPS,L,M,NR,NR2,IERR)
            
              630 FORMAT(1PE15.7)
              640 FORMAT(3(1PE15.7))
            
                  PS1=A(1,1)
                  PS2=A(2,2)
                  PS3=A(3,3)
            
                  DC11=V(1,1)
                  DC12=V(2,1)
                  DC13=V(3,1)
                  DC21=V(1,2)
                  DC22=V(2,2)
                  DC23=V(3,2)
                  DC31=V(1,3)
                  DC32=V(2,3)
                  DC33=V(3,3)
            
            
                !!     KUU=1
                !!     IF(KUU.EQ.1) GO TO 1000
                !!
                  IF(PS2.GT.PS1) THEN
                      AAA=PS2
                      PS2=PS1
                      PS1=AAA
                     AAA1=DC11
                     AAA2=DC12
                     AAA3=DC13
                     DC11=DC21
                     DC12=DC22
                     DC13=DC23
                     DC21=AAA1
                     DC22=AAA2
                     DC23=AAA3
                  ENDIF
                !!       PS1>PS2
                  IF(PS3.GT.PS1) THEN
                      BBB=PS3
                      PS3=PS1
                      PS1=BBB
                     BBB1=DC11
                     BBB2=DC12
                     BBB3=DC13
                     DC11=DC31
                     DC12=DC32
                     DC13=DC33
                     DC31=BBB1
                     DC32=BBB2
                     DC33=BBB3
                  ENDIF
                !!       PS1>PS3
                  IF(PS3.GT.PS2) THEN
                      CCC=PS3
                      PS3=PS2
                      PS2=CCC
                     CCC1=DC21
                     CCC2=DC22
                     CCC3=DC23
                     DC21=DC31
                     DC22=DC32
                     DC23=DC33
                     DC31=CCC1
                     DC32=CCC2
                     DC33=CCC3
                  ENDIF
                !!       PS2>PS3
            
             1000 PS4=0.0D0
                  PS5=0.0D0
                  PS6=0.0D0
                !!     PS4=ABS(PS1-PS2)*0.5D0
                !!     PS5=ABS(PS2-PS3)*0.5D0
                !!     PS6=ABS(PS3-PS1)*0.5D0
                  RETURN
            END 

        !**************************************************************************************
            SUBROUTINE sCrackChk( CCRACK,EU,ECR,ECCR,ECU,EEN,                           &
                                  LCRACK,NTCRACK,NTFAIL,NTELS,NUMCRACK )
                !**************************************************************************************
                !*引用名：sCheckCrack
                !*機�E�E��Eび割れをチェチE��する
                !*入力！EEU    :コンクリート�E三軸応力状態を老E�Eした等価一軸ひずみ
                !*      SN  :コンクリート�E主応力
                !*      IUL :
                !*出力！ENTCRACK   :ひび割れ総数
                !*      NUMCRACK:ひび割れ�E方向番号
                !*      LIUL    :ひび割れ�E持E��E
                !*目皁E��E
                !**************************************************************************************
                    
                real(kreal)  CCRACK(3),EU(:),ECCR(3),EEN(3)
                integer  NTCRACK,NTFAIL,intI,NUMCRACK
                integer NTELS,LCRACK(3)
                real(kreal) ECR,ECU

                NTCRACK=0
                NTFAIL=0

                DO intI=1,3
                    LCRACK(intI)=0
                    IF(NTCRACK.EQ.0) NUMCRACK=intI
            
                    IF(CCRACK(intI).EQ.0.0) THEN
                        IF(EEN(intI).EQ.0.0) THEN               
                            IF(EU(intI).GT.ECR) THEN
                                LCRACK(intI)=1
                                NTCRACK=NTCRACK+1
                            ENDIF
                            IF(EU(intI).LT.ECU) THEN
                            IF(LCRACK(intI).NE.1)LCRACK(intI)=2
                                NTFAIL=NTFAIL+1
                            ENDIF
                        ELSE
                            IF(EU(intI).GT.ECCR(intI)) THEN
                                LCRACK(intI)=1
                                NTCRACK=NTCRACK+1
                            ENDIF
                            IF(EU(intI).LT.ECU) THEN
                            IF(LCRACK(intI).NE.1)LCRACK(intI)=2
                                NTFAIL=NTFAIL+1
                            ENDIF
                        ENDIF
                    ENDIF
            
                    IF(CCRACK(intI).EQ.4.0) THEN
                        IF(EU(intI).GT.ECR) THEN
                            LCRACK(intI)=1
                            NTCRACK=NTCRACK+1
                        ENDIF
                        IF(EU(intI).LT.ECU) THEN
                            IF(LCRACK(intI).NE.1)LCRACK(intI)=2
                            NTFAIL=NTFAIL+1
                        ENDIF
                    ENDIF
                        
                    IF(CCRACK(intI).EQ.5.0) THEN
                        IF(EU(intI).GT.ECCR(intI)) THEN
                            LCRACK(intI)=1
                            NTCRACK=NTCRACK+1
                        ENDIF
                     !  IF(IUL(intI).EQ.8.OR.IUL(intI).EQ.9) THEN
                        IF(EU(intI).LT.ECU) THEN
                            IF(LCRACK(intI).NE.1)LCRACK(intI)=2
                            NTFAIL=NTFAIL+1
                        ENDIF
                    ENDIF
                ENDDO

                NTELS=3-NTCRACK-NTFAIL
            
            END SUBROUTINE
        !**************************************************************************************
            SUBROUTINE sTransIToA(I1,I2,I3,intA,NDIM)
                !IMPLICIT NONE
                !**************************************************************************************
                !*引用名：sTransAToR
                !*機�E�E�整数配�Eから引数に戻ぁE
                !*入力！ENARRAY    :
                !*出力！EN1,N2,N3  :
                !*目皁E��CODING量を減らし、簡単なMISSを行わなぁE��めE
                !********************************************************************
                integer intA(NDIM),I1,I2,I3
                integer NDIM
        
                intA(1)=I1
                intA(2)=I2
                intA(3)=I3

            END SUBROUTINE            
        !**************************************************************************************

            SUBROUTINE sTRANSARRAYreal33(DCC,DAIX)
                !********************************************************************
                !*引用名：STRANSARRAYreal(kreal)33            added by Kengo(2003.4.28)
                !*機�E�E��Eび割れをチェチE��する
                !*入力！EEU    :コンクリート�E三軸応力状態を老E�Eした等価一軸ひずみ
                !*      SN  :コンクリート�E主応力
                !*      IUL :
                !*出力！ENTCRACK   :ひび割れ総数
                !*      NUMCRACK:ひび割れ�E方向番号
                !*      LIUL    :ひび割れ�E持E��E
                !*目皁E��E
                !********************************************************************
                    !!IMPLICIT DOUBLE PRECISION(A-H,O-Z)
                 !     IMPLICIT NONE
      
                real(kreal) DCC(3,3),DAIX(3,3)
                integer intI,intJ
                
                DO intI=1,3
                    DO intJ=1,3
                      DAIX(intI,intJ)=DCC(intI,intJ)
                    ENDDO
                ENDDO
                
            END SUBROUTINE
            
        !**************************************************************************************      
            SUBROUTINE sTransArrayreal(RArraySource,RArrayTarget,IDimension)
      
                !********************************************************************
                !*引用名：sTransArrayreal
                !*機�E�E�実数配�Eの代入
                !*入力！ERArraySource  :ソース配�E
                !*      RArrayTarget    :目標�E刁E
                !*出力！ERArrayTarget  :目標�E刁E
                !*目皁E��E
                !********************************************************************
                  !IMPLICIT NONE
      
                integer intI,IDimension
                real(kreal) RArraySource(IDimension),RArrayTarget(IDimension)
      
    
                DO intI=1,IDimension
                    RArrayTarget(intI) = RArraySource(intI)
                ENDDO

            END SUBROUTINE            
        !**************************************************************************************
            SUBROUTINE sGetIREK(ICIC6,DCC,DCCL,DCCR,DCCRL,NST,                          &
                                IACT,IACTOLD,IREK)
                !********************************************************************
                !*引用名：sGetIREK
                !*機�E�E�SCの変化より特徴点の計算し直ぁE
                !*入力！E  FC  :
                !*出力！EPR12  :
                !*参老E��E
                !********************************************************************
                !IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        
              integer ICIC6,IACT,I,NST,IREK,IACTOLD,intI,intJ
        
                real(kreal) DCC(3,3),DCCL(3,3),DCCR(3,3),DCCRL(3,3) 
                real(kreal) TEP(6),TEPLL(6),TEPC(6),TEPL(6)
        
                !!  IF(NTCRACK.EQ.0) GOTO 3501
                IF(ICIC6.EQ.0)  GOTO 3501
                IF(IACT.EQ.0)       GOTO 3501
                !!************************************************************************
                !!********* CRACK **********************
                !!************************************************************************
               DO 3505 I=1,6
                 TEPL(I)=TEPC(I)
         3505  CONTINUE
                
        
              CALL CHG3 (TEPC(1),TEPC(2),TEPC(3),                       &
                TEPC(4),TEPC(5),TEPC(6),                                &
                DCCR(1,1),DCCR(1,2),DCCR(1,3),                          &
                DCCR(2,1),DCCR(2,2),DCCR(2,3),                          &
                DCCR(3,1),DCCR(3,2),DCCR(3,3),                          &
                TEP(1),TEP(2),TEP(3),                                   &
                TEP(4),TEP(5),TEP(6),3)
        
              IF(NST.GT.1) THEN
              CALL DIRM (                                               &
                DCCRL(1,1),DCCRL(1,2),DCCRL(1,3),                       &
                DCCRL(2,1),DCCRL(2,2),DCCRL(2,3),                       &
                DCCRL(3,1),DCCRL(3,2),DCCRL(3,3),                       &
                DCCR(1,1),DCCR(1,2),DCCR(1,3),                          &
                DCCR(2,1),DCCR(2,2),DCCR(2,3),                          &
                DCCR(3,1),DCCR(3,2),DCCR(3,3),                          &
                                        IREK)
              ENDIF
               GOTO 3508
                !!************************************************************************
                !!********* CRACK **********************
                !!************************************************************************
         3501  CONTINUE
        
    
                
               CALL CHG3 (TEPL(1),TEPL(2),TEPL(3),                      &
                TEPL(4),TEPL(5),TEPL(6),                                &
                DCC(1,1),DCC(1,2),DCC(1,3),                             &
                DCC(2,1),DCC(2,2),DCC(2,3),                             &
                DCC(3,1),DCC(3,2),DCC(3,3),                             &
                TEPLL(1),TEPLL(2),TEPLL(3),                             &
                TEPLL(4),TEPLL(5),TEPLL(6),3)
        
               CALL CHG3 (TEPC(1),TEPC(2),TEPC(3),                      &
                TEPC(4),TEPC(5),TEPC(6),                                &
                DCC(1,1),DCC(1,2),DCC(1,3),                             &
                DCC(2,1),DCC(2,2),DCC(2,3),                             &
                DCC(3,1),DCC(3,2),DCC(3,3),                             &
                TEP(1),TEP(2),TEP(3),                                   &
                TEP(4),TEP(5),TEP(6),3)
        
        
                IF(NST.GT.1) THEN
                    !!  ***[固定�Eび割れ座標系を経歴した場吁E***
                    IF(IACTOLD.NE.0) THEN 
                    DO intI = 1,3
                        DO intJ = 1,3
                            DCCL(inti,intJ) = DCCRL(inti,intJ)
                        ENDDO
                    ENDDO
               ENDIF
        
              CALL DIRM (                                               &
                DCCL(1,1),DCCL(1,2),DCCL(1,3),                          &
                DCCL(2,1),DCCL(2,2),DCCL(2,3),                          &
                DCCL(3,1),DCCL(3,2),DCCL(3,3),                          &
                DCC(1,1),DCC(1,2),DCC(1,3),                             &
                DCC(2,1),DCC(2,2),DCC(2,3),                             &
                DCC(3,1),DCC(3,2),DCC(3,3),                             &
                IREK)
                ENDIF
        
                !!  IREK=0
         3508 CONTINUE
    

            END SUBROUTINE
        !**************************************************************************************
            SUBROUTINE DIRM (                                                           &
                                  DC11,DC12,DC13,                       &
                                  DC21,DC22,DC23,                       &
                                  DC31,DC32,DC33,                       &
                                  DD11,DD12,DD13,                       &
                                  DD21,DD22,DD23,                       &
                                  DD31,DD32,DD33,IR                     )
                ! !***********************************************************************
                ! !
                ! ! *** CHECK OF ROTATION OF DIRECTION
                ! !
                    !IMPLICIT DOUBLE PRECISION(A-H,O-Z)
                  ! IMPLICIT NONE
            
                real(kreal) CK1A,CK2A,CK3A      
                real(kreal) CK1B,CK2B,CK3B      
                real(kreal) DC11,DC12,DC13      
                real(kreal) DC21,DC22,DC23      
                real(kreal) DC31,DC32,DC33      
                real(kreal) DD11,DD12,DD13      
                real(kreal) DD21,DD22,DD23      
                real(kreal) DD31,DD32,DD33
                !!
                  integer IR
                !!
                  CK1A=ABS(DC11*DD11+DC12*DD12+DC13*DD13)
                  CK2A=ABS(DC11*DD21+DC12*DD22+DC13*DD23)
                  CK3A=ABS(DC11*DD31+DC12*DD32+DC13*DD33)
                !!
                  CK1B=ABS(DC21*DD11+DC22*DD12+DC23*DD13)
                  CK2B=ABS(DC21*DD21+DC22*DD22+DC23*DD23)
                  CK3B=ABS(DC21*DD31+DC22*DD32+DC23*DD33)
                !!
                  IF(CK1A.GE.CK2A) THEN
                    IF(CK1A.GE.CK3A) THEN
                      GO TO 100
                    ELSE
                      GO TO 300
                    END IF
                  ELSE
                    IF(CK2A.GE.CK3A) THEN
                      GO TO 200
                    ELSE
                      GO TO 300
                    END IF
                  END IF
                !!
              100 IF(CK2B.GE.CK3B) THEN
                    GO TO 123
                  ELSE
                    GO TO 132
                  END IF
                  GO TO 500
                !!
              200 IF(CK1B.GE.CK3B) THEN
                    GO TO 213
                  ELSE
                    GO TO 231
                  END IF
                  GO TO 500
                !!
              300 IF(CK1B.GE.CK2B) THEN
                    GO TO 312
                  ELSE
                    GO TO 321
                  END IF
                !!
              500 CONTINUE
                ! !
                ! ! *** 1 ==> 1 ***
                ! !
                ! !   * 2 - 2 & 3 - 3
                ! !
              123 CONTINUE
                  IR=0
                  RETURN
                ! !
                ! !   * 2 - 3 & 3 - 2
                ! !
              132 CONTINUE
                  IR=1
                  RETURN
                ! !
                ! ! *** 1 ==> 2 ***
                ! !
                ! !   * 2 - 1 & 3 - 3
                ! !
              213 CONTINUE
                  IR=2
                  RETURN
                ! !
                ! !   * 2 - 3 & 3 - 1
                ! !
              231 CONTINUE
                  IR=3
                  RETURN
                ! !
                ! ! *** 1 ==> 3 ***
                ! !
                ! !   * 2 - 2 & 3 - 1
                ! !
              321 CONTINUE
                  IR=4
                  RETURN
                ! !
                ! !   * 2 - 1 & 3 - 2
                ! !
              312 CONTINUE
                  IR=5
                  RETURN
            END SUBROUTINE
            
        !**************************************************************************************
            SUBROUTINE CHG3 (RS1,RS2,RS3,                                               &
                             RS4,RS5,RS6,                                 &
                             D11,D12,D13,                                 &
                             D21,D22,D23,                                 &
                             D31,D32,D33,                                 &
                             SR1,SR2,SR3,                                 &
                             SR4,SR5,SR6,N)
                !!***********************************************************************
                !!
                !! *** CALCULATE CORDINTE TRANSLATION MATRIX
                !!     N = 0 : STRESS ,1-2-3 ==> X-Y-Z
                !!     N = 1 : STRESS ,X-Y-Z ==> 1-2-3
                !!     N = 2 : STRAIN ,1-2-3 ==> X-Y-Z
                !!     N = 3 : STRAIN ,X-Y-Z ==> 1-2-3
                !!
              !IMPLICIT DOUBLE PRECISION(A-H,O-Z)
              !IMPLICIT NONE
                ! add for implicit none
              real(kreal) SRE1(3,3),DC1(3,3),DC1T(3,3),SRED(3,3)
              real(kreal) SRE2(3,3)
              integer I,J,K,N
              real(kreal) D11,D12,D13,D21,D22,D23,D31,D32,D33
              real(kreal) RS1,RS2,RS3,RS4,RS5,SR1,SR2,SR3,SR4,SR5,SR6,RS6
                !!     real(kreal) RS(6),SRE2(3,3)
                !!
              IF(N.EQ.0.OR.N.EQ.1) GO TO 100
              IF(N.EQ.2.OR.N.EQ.3) GO TO 200
              GO TO 300
                !!
                !!  *** STRESS *************************************************
         100  SRE1(1,1)=SR1
              SRE1(1,2)=SR4
              SRE1(1,3)=SR6
              SRE1(2,1)=SR4
              SRE1(2,2)=SR2
              SRE1(2,3)=SR5
              SRE1(3,1)=SR6
              SRE1(3,2)=SR5
              SRE1(3,3)=SR3
            
              DC1(1,1)=D11
              DC1(1,2)=D12
              DC1(1,3)=D13
              DC1(2,1)=D21
              DC1(2,2)=D22
              DC1(2,3)=D23
              DC1(3,1)=D31
              DC1(3,2)=D32
              DC1(3,3)=D33
            
              DC1T(1,1)=D11
              DC1T(1,2)=D21
              DC1T(1,3)=D31
              DC1T(2,1)=D12
              DC1T(2,2)=D22
              DC1T(2,3)=D32
              DC1T(3,1)=D13
              DC1T(3,2)=D23
              DC1T(3,3)=D33
            
              IF(N.EQ.1) GO TO 150
            
                !!  ** 1-2-3 ==> X-Y-Z (STRESS)**
                !!
              DO 20 I=1,3
              DO 20 J=1,3
              SRED(I,J)=0.0D0
              DO 20 K=1,3
              SRED(I,J)=SRED(I,J)+DC1T(I,K)*SRE1(K,J)
           20 CONTINUE
            
              DO 30 I=1,3
              DO 30 J=1,3
              SRE2(I,J)=0.0D0
              DO 30 K=1,3
              SRE2(I,J)=SRE2(I,J)+SRED(I,K)*DC1(K,J)
           30 CONTINUE
            
              RS1=SRE2(1,1)
              RS2=SRE2(2,2)
              RS3=SRE2(3,3)
              RS4=SRE2(1,2)
              RS5=SRE2(2,3)
              RS6=SRE2(1,3)
                !!     WRITE(6,*) 'RS1',RS1,'RS2',RS2,'RS3',RS3
                !!     WRITE(6,*) 'RS4',RS4,'RS5',RS5,'RS6',RS6
              RETURN
                !!
                !!  ** X-Y-Z ==> 1-2-3 (STRESS)**
                !!
         150  DO 25 I=1,3
              DO 25 J=1,3
              SRED(I,J)=0.0D0
              DO 25 K=1,3
              SRED(I,J)=SRED(I,J)+DC1(I,K)*SRE1(K,J)
           25 CONTINUE
            
              DO 35 I=1,3
              DO 35 J=1,3
              SRE2(I,J)=0.0D0
              DO 35 K=1,3
              SRE2(I,J)=SRE2(I,J)+SRED(I,K)*DC1T(K,J)
           35 CONTINUE
            
              RS1=SRE2(1,1)
              RS2=SRE2(2,2)
              RS3=SRE2(3,3)
              RS4=SRE2(1,2)
              RS5=SRE2(2,3)
              RS6=SRE2(1,3)
              RETURN
            
                !! *** STRAIN *************************************************
         200  SRE1(1,1)=SR1
              SRE1(1,2)=SR4*0.5D0
              SRE1(1,3)=SR6*0.5D0
              SRE1(2,1)=SR4*0.5D0
              SRE1(2,2)=SR2
              SRE1(2,3)=SR5*0.5D0
              SRE1(3,1)=SR6*0.5D0
              SRE1(3,2)=SR5*0.5D0
              SRE1(3,3)=SR3
            
              DC1(1,1)=D11
              DC1(1,2)=D12
              DC1(1,3)=D13
              DC1(2,1)=D21
              DC1(2,2)=D22
              DC1(2,3)=D23
              DC1(3,1)=D31
              DC1(3,2)=D32
              DC1(3,3)=D33
            
              DC1T(1,1)=D11
              DC1T(1,2)=D21
              DC1T(1,3)=D31
              DC1T(2,1)=D12
              DC1T(2,2)=D22
              DC1T(2,3)=D32
              DC1T(3,1)=D13
              DC1T(3,2)=D23
              DC1T(3,3)=D33
            
              IF(N.EQ.3) GOTO 250
                !!
                !!  ** 1-2-3 ==> X-Y-Z (STRAIN)**
                !!
              DO 40 I=1,3
              DO 40 J=1,3
              SRED(I,J)=0.0D0
              DO 40 K=1,3
              SRED(I,J)=SRED(I,J)+DC1T(I,K)*SRE1(K,J)
           40 CONTINUE
            
              DO 50 I=1,3
              DO 50 J=1,3
              SRE2(I,J)=0.0D0
              DO 50 K=1,3
              SRE2(I,J)=SRE2(I,J)+SRED(I,K)*DC1(K,J)
           50 CONTINUE
            
              RS1=SRE2(1,1)
              RS2=SRE2(2,2)
              RS3=SRE2(3,3)
              RS4=SRE2(1,2)*2.0D0
              RS5=SRE2(2,3)*2.0D0
              RS6=SRE2(1,3)*2.0D0
              RETURN
                !!
                !!  ** X-Y-Z ==> 1-2-3 (STRAIN)**
                !!
         250  DO 60 I=1,3
              DO 60 J=1,3
              SRED(I,J)=0.0D0
              DO 60 K=1,3
              SRED(I,J)=SRED(I,J)+DC1(I,K)*SRE1(K,J)
           60 CONTINUE
            
              DO 70 I=1,3
              DO 70 J=1,3
              SRE2(I,J)=0.0D0
              DO 70 K=1,3
              SRE2(I,J)=SRE2(I,J)+SRED(I,K)*DC1T(K,J)
           70 CONTINUE
            
              RS1=SRE2(1,1)
              RS2=SRE2(2,2)
              RS3=SRE2(3,3)
              RS4=SRE2(1,2)*2.0D0
              RS5=SRE2(2,3)*2.0D0
              RS6=SRE2(1,3)*2.0D0
              RETURN
            
          300 WRITE(6,601)
          601 FORMAT('  *** ERROR IN CHG3 ***')
              RETURN
        END SUBROUTINE                          
        !**************************************************************************************
            SUBROUTINE sChangeAxisPR(PR23,PR31,PR12,RG23,RG31,RG12,IR)
                ! ********************************************************************
                !*引用名：sChangeAxis           added by TIDE(2003.4.28)
                !*機�E�E�回転ひび割れモチE���E�主軸�E�でのひずみの計箁E
                !*入力！EEU    :コンクリート�E三軸応力状態を老E�Eした等価一軸ひずみ
                !*      SN  :コンクリート�E主応力
                !*      IUL :
                !*出力！ENTCRACK   :ひび割れ総数
                !*      NUMCRACK:ひび割れ�E方向番号
                !*      LIUL    :ひび割れ�E持E��E
                !*目皁E��E
                !********************************************************************
                !add  
              integer IR
              real(kreal) PR12,PR23,PR31,RG12,RG23,RG31
            
              IF(IR.EQ.0) THEN
                !!   *** 1 ==> 1  &  2 ==> 2  &  3 ==> 3 ***
                ELSE
                CONTINUE
              ENDIF
            
               IF(IR.EQ.1) THEN
                !!   *** 2 <==> 3 ***
                CALL XRVS1(PR31,PR12)
                CALL XRVS1(RG31,RG12)
            
               ENDIF
            
               IF(IR.EQ.2) THEN
                !!   *** 1 <==> 2 ***
                CALL XRVS1(PR31,PR23)
                CALL XRVS1(RG31,RG23)
                
               ENDIF
            
               IF(IR.EQ.3) THEN
                !!   *** 1 ==> 2  &  2 ==> 3  &  3 ==> 1 ***
                CALL XRVS2(PR23,PR31,PR12)
                CALL XRVS2(RG23,RG31,RG12)
                
               ENDIF
            
               IF(IR.EQ.4) THEN
                !!   *** 1 <==> 3 ***
                CALL XRVS1(PR23,PR12)
                CALL XRVS1(RG23,RG12)
                
               ENDIF
            
               IF(IR.EQ.5) THEN
                !!   *** 1 ==> 3  &  2 ==> 1  &  3 ==> 2 ***
                CALL XRVS3(PR23,PR31,PR12)
                CALL XRVS3(RG23,RG31,RG12)
                
               ENDIF
            
        END SUBROUTINE   
        !**************************************************************************************
            SUBROUTINE XRVSI1 (IDATA1,IDATA2)
                !!***********************************************************************
                !!
                !IMPLICIT DOUBLE PRECISION(A-H,O-Z)
                  !IMPLICIT NONE
              real(kreal) IDUMMY,IDATA1,IDATA2,IDATA3
            
              IDUMMY=IDATA1
              IDATA1=IDATA2
              IDATA2=IDUMMY
              RETURN
            END SUBROUTINE
        !**************************************************************************************
            SUBROUTINE XRVSI2 (IDATA1,IDATA2,IDATA3)
        
                !!
                !!     1 ==> 2   2 ==> 3   3 ==> 1
                !!
                !IMPLICIT DOUBLE PRECISION(A-H,O-Z)
                  !IMPLICIT NONE
        
              real(kreal) IDUMMY,IDATA1,IDATA2,IDATA3
            
              IDUMMY=IDATA1
              IDATA1=IDATA2
              IDATA2=IDATA3
              IDATA3=IDUMMY
              RETURN
            END SUBROUTINE
        !**************************************************************************************
            SUBROUTINE XRVSI3 (IDATA1,IDATA2,IDATA3)
        
                !!     1 ==> 3   2 ==> 1   3 ==> 2
                !!
                !IMPLICIT DOUBLE PRECISION(A-H,O-Z)
              !IMPLICIT NONE
              real(kreal) IDUMMY,IDATA1,IDATA2,IDATA3
        
              IDUMMY=IDATA1
              IDATA1=IDATA3
              IDATA3=IDATA2
              IDATA2=IDUMMY
              RETURN
            END SUBROUTINE
        !**************************************************************************************
            SUBROUTINE XRVS1 (DATA1,DATA2)
        
                !IMPLICIT DOUBLE PRECISION(A-H,O-Z)
                  !IMPLICIT NONE
              real(kreal) DUMMY,DATA1,DATA2,DATA3
    
              DUMMY=DATA1
              DATA1=DATA2
              DATA2=DUMMY
              RETURN
            END SUBROUTINE
        !**************************************************************************************
            SUBROUTINE XRVS2 (DATA1,DATA2,DATA3)
        
                !!     1 ==> 2   2 ==> 3   3 ==> 1
                !!
                !IMPLICIT DOUBLE PRECISION(A-H,O-Z)
                  !IMPLICIT NONE
              real(kreal) DUMMY,DATA1,DATA2,DATA3
        
              DUMMY=DATA1
              DATA1=DATA2
              DATA2=DATA3
              DATA3=DUMMY
              RETURN
            END SUBROUTINE
        !**************************************************************************************
            SUBROUTINE XRVS3 (DATA1,DATA2,DATA3)
        
                !!     1 ==> 3   2 ==> 1   3 ==> 2
                !!
                !IMPLICIT DOUBLE PRECISION(A-H,O-Z)
                  !IMPLICIT NONE
              real(kreal) DUMMY,DATA1,DATA2,DATA3
        
              DUMMY=DATA1
              DATA1=DATA3
              DATA3=DATA2
              DATA2=DUMMY
              RETURN
            END SUBROUTINE                             
                                     
        !**************************************************************************************
            SUBROUTINE sEURote(EU,EUL,DEU,TEP,TEPLL,ET,DCC,DCCL,                        &
                               PR23,PR31,PR12,RG23,RG31,RG12)
                !!********************************************************************
                !*引用名：sEURote           added by Kengo(2003.4.28)
                !*機�E�E�回転ひび割れモチE���E�主軸�E�でのひずみの計箁E
                !*入力！EEU        :前stepの等価一軸ひずみ�E��E力として�E�E
                !*      TEP     :コンクリート�E応力�E��Eび割れ軸�E�E
                !*      DEP     :コンクリート�E応力増�E�E��Eび割れ軸�E�E
                !*      TEPL    :前stepコンクリート�E応力�E��Eび割れ軸�E�E
                !*      ET      :コンクリート�Eヤング係数
                !*      DCC     :主軸座標系
                !*      DCCL    :前主軸座標系
                !*      PR12    :ポアソン毁E
                !*出力！EEU        :前stepの等価一軸ひずみ�E��E力として�E�E
                !*　　　   DEU     :ひび割れ総数
                !*目皁E��E
                !********************************************************************
                !     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
                  !IMPLICIT NONE
            
                real(kreal) TEPLL(6),TEP(6),EU(6),EULD(6),EUL(6)
                real(kreal) TEPC(6),TEPL(6)
                real(kreal) DEU(6),DEP(6),ET(3)
                real(kreal) DCC(3,3),DCCL(3,3)
                !add 
                real(kreal) FAI,PR12,PR23,PR31,RG12,RG23,RG31
                
                DEU=0.0D0
                DEP=0.0D0
            
               CALL CHG3 (TEPL(1),TEPL(2),TEPL(3),                      &
                TEPL(4),TEPL(5),TEPL(6),                                &
                DCC(1,1),DCC(1,2),DCC(1,3),                             &
                DCC(2,1),DCC(2,2),DCC(2,3),                             &
                DCC(3,1),DCC(3,2),DCC(3,3),                             &
                TEPLL(1),TEPLL(2),TEPLL(3),                             &
                TEPLL(4),TEPLL(5),TEPLL(6),3)
            
               CALL CHG3 (TEPC(1),TEPC(2),TEPC(3),                      &
                TEPC(4),TEPC(5),TEPC(6),                                &
                DCC(1,1),DCC(1,2),DCC(1,3),                             &
                DCC(2,1),DCC(2,2),DCC(2,3),                             &
                DCC(3,1),DCC(3,2),DCC(3,3),                             &
                TEP(1),TEP(2),TEP(3),                                   &
                TEP(4),TEP(5),TEP(6),3)

                CALL CHG3(EULD(1),EULD(2),EULD(3),                          &
                      EULD(4),EULD(5),EULD(6),                          &
                    DCCL(1,1),DCCL(1,2),DCCL(1,3),                  &
                    DCCL(2,1),DCCL(2,2),DCCL(2,3),                  &
                    DCCL(3,1),DCCL(3,2),DCCL(3,3),                  &
                        EU(1),EU(2),EU(3),                              &
                        EU(4),EU(5),EU(6),2)
            
                CALL CHG3(EUL(1),EUL(2),EUL(3),                           &
                          EUL(4),EUL(5),EUL(6),                           &
                                DCC(1,1),DCC(1,2),DCC(1,3),             &
                                DCC(2,1),DCC(2,2),DCC(2,3),             &
                                DCC(3,1),DCC(3,2),DCC(3,3),             &
                         EULD(1),EULD(2),EULD(3),                       &
                         EULD(4),EULD(5),EULD(6),3)
            
              DEP(1)=TEPC(1)-TEPL(1)
              DEP(2)=TEPC(2)-TEPL(2)
              DEP(3)=TEPC(3)-TEPL(3)
              
            
              FAI=1.D0-PR12*PR12-PR23*PR23-PR31*PR31-2.D0*PR12*PR23*PR31
              
              DEU(1)=((1.0D0-PR23*PR23)*DEP(1)                          &
                     +(PR31*PR23+PR12)*SQRT(ET(2)/ET(1))*DEP(2)         &
                     +(PR12*PR23+PR31)*SQRT(ET(3)/ET(1))*DEP(3))/FAI
              
              DEU(2)=((PR31*PR23+PR12)*SQRT(ET(1)/ET(2))*DEP(1)         &
                    +(1.0D0-PR31*PR31)*DEP(2)                           &
                    +(PR12*PR31+PR23)*SQRT(ET(3)/ET(2))*DEP(3))/FAI
              
              DEU(3)=((PR12*PR23+PR31)*SQRT(ET(1)/ET(3))*DEP(1)         &
                    +(PR12*PR31+PR23)*SQRT(ET(2)/ET(3))*DEP(2)          &
                    +(1.0D0-PR12*PR12)*DEP(3))/FAI

                !c      DEU(1)=(1.0D0-PR23)*DEP(1)/(1.0D0-2.0D0*PR23)/(1.0D0+PR23)
                !c     *            +PR31*DEP(2)/(1.0D0-2.0D0*PR31)/(1.0D0+PR31)
                !c     *            +PR12*DEP(3)/(1.0D0-2.0D0*PR12)/(1.0D0+PR12)
                !c      DEU(2)=        PR23*DEP(1)/(1.0D0-2.0D0*PR23)/(1.0D0+PR23)
                !c     *    +(1.0D0-PR31)*DEP(2)/(1.0D0-2.0D0*PR31)/(1.0D0+PR31)
                !c     *            +PR12*DEP(3)/(1.0D0-2.0D0*PR12)/(1.0D0+PR12)
                !c      DEU(3)=        PR23*DEP(1)/(1.0D0-2.0D0*PR23)/(1.0D0+PR23)
                !c     *            +PR31*DEP(2)/(1.0D0-2.0D0*PR31)/(1.0D0+PR31)
                !c     *    +(1.0D0-PR12)*DEP(3)/(1.0D0-2.0D0*PR12)/(1.0D0+PR12)
                !!
              
              !delete by shan
              if (deu(1).gt.-1D-10.and.deu(1).lt.+1d-10)deu(1)=0.0D0
              if (deu(2).gt.-1D-10.and.deu(2).lt.+1d-10)deu(2)=0.0D0
              if (deu(3).gt.-1D-10.and.deu(3).lt.+1d-10)deu(3)=0.0D0

              EU(1)=EUL(1)+DEU(1)
              EU(2)=EUL(2)+DEU(2)
              EU(3)=EUL(3)+DEU(3)
              
              !delete by shan
              if (eu(1).gt.-1D-10.and.eu(1).lt.+1d-10)eu(1)=0.0D0
              if (eu(2).gt.-1D-10.and.eu(2).lt.+1d-10)eu(2)=0.0D0
              if (eu(3).gt.-1D-10.and.eu(3).lt.+1d-10)eu(3)=0.0D0
              EU(4)=0.0D0
              EU(5)=0.0D0
              EU(6)=0.0D0
            END SUBROUTINE
        
        !**************************************************************************************
            SUBROUTINE sRembCrack(DAIX,DAIXL,DCCRCHG,DCCRCHGL,N)
                !********************************************************************
                    !*引用名：sRembCrack
                    !*機�E�E��Eび割れ暦の記録
                    !*入力！EDAIX,DAIXL:
                    !*      N                   :アクチE��ブ�Eび割れ座標系の番号
                    !*出力！EDCCRCHG,DCCRCHGL:
                    !********************************************************************
                            
                real(kreal) DAIX(3,3),DAIXL(3,3)
                real(kreal) DCCRCHG(3,3,3),DCCRCHGL(3,3,3)
                !add
                  integer N,INTI,INTJ
                !!
                DO intI=1,3
                  DO intJ=1,3
                    DCCRCHG(N,intI,intJ)=DAIX(intI,intJ)
                    DCCRCHGL(N,intI,intJ)=DAIXL(intI,intJ)
                  ENDDO
                ENDDO
                !!
            END SUBROUTINE                          
                                      
        !**************************************************************************************
            SUBROUTINE sJudgeNewCrack(DCC,DCCRCHG,ICOORDCHGED,                          &
                        IJNC,LCRACK,EUPSI,EUCR,IACTIVEC,IACTOLD,                            &
                        DEU,TEP,TEPLL,ET,DCCRL,DCCL,                                        &
                        IROT,CCRACK,ECR,ECCR,ECU,EEN,                                       &
                        PRDM23,PRDM31,PRDM12,RGDM23,RGDM31,RGDM12)
                !c     *    PR23,PR31,PR12,RG23,RG31,RG12)
                    !********************************************************************
                    !*引用名：sJudgeNewCrack
                    !*機�E�E�新たなひび割れ発生�E条件に満たすか�E判断
                    !*入力！EDCC           :主軸の座標�EトリチE��ス
                    !*      DCCRCHG :ひび割れ暦を記録する座標�EトリチE��ス
                    !*      ICCHG       :総�Eび割れ暦数
                    !*      EUPSI       :主軸座標系でのひずみ
                    !*      EUCR        :ひび割れ座標系でのひずみ
                    !*      LCARCK  :座標系の中吁E��ひび割れてぁE��か�E持E��E
                    !*出力！EIJNC      :ひび割れ座標系が転換するかの持E��E
                    !*      IACTIVEC:アクチE��ブ�Eび割れ座標系の番号
                    !********************************************************************
                
                    real(kreal) DCCRCHG(3,3,3),DCC(3,3),DCCRDMY(3,3),DCCRL(3,3)
                    real(kreal) EUPSI(6),EUCR(6),EUL(6)
                    real(kreal) TEP(6),TEPLL(6)
                    real(kreal) DD(3,3),DDMAX(3),EUCHG(3)
                    real(kreal) CCRACK(3),ECCR(3)
                    real(kreal) ET(3),EEN(3)
                    real(kreal) EUCRDM(6),DEUDM(6)
                    real(kreal) AAAA(3,3),BBBB(3,3),CCCC(3),DCCRRRR(3,3)
                !add 
                  !c real(kreal) ECCR(3),ECR(3),ECU(3)
                    integer NST,NIT,NTFAIL,NTELS,NUMCRACK,LCRACK(3)
                  !integer IACT,IACTOLD,ICOORDCHGED,IACTIVEC,IACTIVECOLD,ICIC6,IREK
                    integer IACTOLD,IARGCLS,IREK,IROT,LIUL,NNN
                    integer  IMAINC(3),IJNCDMY(3)
                  !integer IJNC,IROT
                    real(kreal) PR12,PR23,PR31,RG12,RG23,RG31,RGDM12,RGDM23,RGDM31
                    real(kreal) PRDM12,PRDM23,PRDM31
                    integer INTI,INTJ,INTM,INTN,INTK,IACTIVEC,IACTIVECOLD,IJNC
                    integer I,J,K,M,NTCRACK,ICOORDCHGED
                    real(kreal) A,B,C,DCCL(3,3),DEU(6),ECR,ECU,DCCR(3,3),EUCRMAXTEMP,TEMPA
                !!
                IF(IACTIVEC.EQ.0) THEN
                    CALL sTRANSARRAYreal33(DCCL,DCCRL)
                ENDIF
                !!
                IJNC=0
                !!  ***[主軸1とひび割れ座標系履歴中の主なひび割れ方向�E角度]***
                DO intI=1,ICOORDCHGED   
                !!
                    IJNCDMY(intI)=0
                    IARGCLS=0
                    EUCHG(intI)=-1000.0
                !!
                    PR23=PRDM23
                    PR31=PRDM31
                    PR12=PRDM12
                !!
                    RG23=RGDM23
                    RG31=RGDM31
                    RG12=RGDM12
                !c
                    DO intM=1,3
                    DO intN=1,3
                        DCCRDMY(intM,intN)=DCCRCHG(intI,intM,intN)
                    ENDDO
                    ENDDO
                    CALL sTransArrayreal(EUCR,EUCRDM,6)
                    CALL sTransArrayreal(DEU,DEUDM,6)
                    
                !!          ***[ひび割れ座標系でのひずみの計算]***
                    CALL sGetIREK(2,DCC,DCCL,DCCRDMY,DCCRL,NST,                     &
                                    intI,IACTOLD,IREK)
                    CALL sChangeAxisPR(PR23,PR31,PR12,RG23,RG31,RG12,IREK)
                !!
                    CALL sEUFixed(EUCRDM,EUL,DEUDM,TEP,TEPLL,ET,DCCRDMY,DCCRL,      &
                                                    PR23,PR31,PR12,RG23,RG31,RG12)
                    CALL sCheckCMODEL(IROT,CCRACK,EUCRDM,ECR,ECCR,ECU,EEN,          &
                                      LCRACK,NTCRACK,NTFAIL,NTELS,LIUL,NUMCRACK)
                !!
                !!      CALL sCrack1ST(TEP,DCCRDMY,NTCRACK,LCRACK)
                !!
                    DO intJ=1,3
                !c          IF(LCRACK(intJ).EQ.1.AND.EUCR(intJ).GE.EUCHG(intI)) THEN
                        IF(EUCRDM(intJ).GE.EUCHG(intI)) THEN
                            EUCHG(intI)=EUCRDM(intJ)
                            IMAINC(intI)=intJ
                        ENDIF
                    ENDDO
                !!
                    DO intM=1,3
                    DO intN=1,3
                        DCCRRRR(intM,intN)=DCCRCHG(intI,intM,intN)
                !c          AAAA(intM,intN)=ACOSD(DCCRRRR(intM,intN))
                !c          BBBB(intM,intN)=ACOSD(DCC(intM,intN))
                    ENDDO
                    ENDDO
                    DO intJ=1,3
                        DD(intI,intJ)=DCC(1,1)*DCCRCHG(intI,intJ,1)                 &
                                     +DCC(1,2)*DCCRCHG(intI,intJ,2)                 &      
                                     +DCC(1,3)*DCCRCHG(intI,intJ,3)
                !c          CCCC(intJ)=ACOSD(DD(intI,intJ))
                    ENDDO
                !!  B=ACOS(DD(intI))*180.0/3.1415926
                !**********
                    A=1.2*EUCHG(intI)
                    B=COS(25.0*3.1415926/180)
                    !=COS(65.0*3.1415926/180)
                !!
                    IF(EUPSI(1).GT.A) THEN
                        IF(ABS(DD(intI,IMAINC(intI))).LT.B) THEN
                !c          IF(ABS(DD(intI,IMAINC(intI))).GT.!) THEN
                            IJNCDMY(intI)=1
                            DO intK=1,3
                                IF(intK.NE.IMAINC(intI)) THEN
                                    TempA=ABS(DD(intI,intK))
                !c                      IF(TempA.GT.B.OR.TempA.LT.!) IJNCDMY(intI)=0
                                    IF(TempA.GT.B) then 
                                        IJNCDMY(intI)=0
                                    endif
                                ENDIF                   
                            ENDDO
                        ENDIF 
                    ENDIF
                !!
                ENDDO
        
                NNN=IJNCDMY(1)+IJNCDMY(2)+IJNCDMY(3)
                !!  
                !!  ***[ICOORDCHGEDは1あるぁE�E2の場合、�E部の座標系履歴に対して条件満たさなぁE��IJNCは1にならない]***
                IF(NNN.EQ.ICOORDCHGED.AND.ICOORDCHGED.NE.3) THEN
                    IJNC=1
                    CALL sTRANSARRAYreal33(DCC,DCCR)
                    RETURN
                ENDIF
                !!
                    EUCRMAXTEMP=-10000.0
                DO intI=1,ICOORDCHGED   
                !!      IF (EUCRMAX(intI).GT.EUCRMAXTEMP.AND.IJNCDMY(intI).EQ.0) THEN
                    IF (EUCHG(intI).GT.EUCRMAXTEMP) THEN
                        EUCRMAXTEMP=EUCHG(intI)
                        IARGCLS=intI
                    ENDIF
                ENDDO
                !!
              888 CONTINUE
                IJNC=0
                IACTIVEC=IARGCLS
                RETURN
            END SUBROUTINE
        
        !**************************************************************************************
            SUBROUTINE sEUFixed(EU,EUL,DEU,TEP,TEPLL,ET,DCCR,DCCRL,                     &
                                PR23,PR31,PR12,RG23,RG31,RG12)
                !********************************************************************
                !*引用名：sEUFixed          added by Kengo(2003.4.28)
                !*機�E�E�固定�Eび割れモチE���E��Eび割れ軸�E�でのひずみの計箁E
                !*入力！EEU        :前stepの等価一軸ひずみ�E��E力として�E�E
                !*      TEP     :コンクリート�E応力�E��Eび割れ軸�E�E
                !*      DEP     :コンクリート�E応力増�E�E��Eび割れ軸�E�E
                !*      TEPL    :前stepコンクリート�E応力�E��Eび割れ軸�E�E
                !*      ET      :コンクリート�Eヤング係数
                !*      DCCR    :ひび割れ座標系
                !*      DCCRL   :前�Eび割れ座標系
                !*      PR12    :ポアソン毁E
                !*出力！EEU        :前stepの等価一軸ひずみ�E��E力として�E�E
                !*　　　   DEU     :ひび割れ総数
                !*目皁E��E
                !********************************************************************
            
              !IMPLICIT NONE
            
                real(kreal) TEP(6),TEPLL(6),EU(6),EULD(6),EUL(6)
                real(kreal) TEPC(6),TEPL(6)
                real(kreal) DEU(6),DEP(6),ET(3)
                real(kreal) DCCR(3,3),DCCRL(3,3)
                !add 
                  real(kreal) FAI,PR12,PR23,PR31,RG12,RG23,RG31
                    DEU=0.0D0
                    
                CALL CHG3(EULD(1),EULD(2),EULD(3),                              &
                          EULD(4),EULD(5),EULD(6),                            &
                                    DCCRL(1,1),DCCRL(1,2),DCCRL(1,3),           &
                                    DCCRL(2,1),DCCRL(2,2),DCCRL(2,3),           &
                                    DCCRL(3,1),DCCRL(3,2),DCCRL(3,3),           &
                           EU(1),EU(2),EU(3),                                  &
                           EU(4),EU(5),EU(6),2)
                   
                   CALL CHG3 (TEPL(1),TEPL(2),TEPL(3),                          &
                              TEPL(4),TEPL(5),TEPL(6),                                    &
                                DCCR(1,1),DCCR(1,2),DCCR(1,3),                              &
                                DCCR(2,1),DCCR(2,2),DCCR(2,3),                              &
                                DCCR(3,1),DCCR(3,2),DCCR(3,3),                              &
                             TEPLL(1),TEPLL(2),TEPLL(3),                                 &
                             TEPLL(4),TEPLL(5),TEPLL(6),3)
            
                  CALL CHG3 (TEPC(1),TEPC(2),TEPC(3),                               &
                             TEPC(4),TEPC(5),TEPC(6),                               &
                            DCCR(1,1),DCCR(1,2),DCCR(1,3),                          &
                            DCCR(2,1),DCCR(2,2),DCCR(2,3),                          &
                            DCCR(3,1),DCCR(3,2),DCCR(3,3),                          &
                             TEP(1),TEP(2),TEP(3),                                  &
                             TEP(4),TEP(5),TEP(6),3)
            
                  CALL CHG3(EUL(1),EUL(2),EUL(3),                                   &
                            EUL(4),EUL(5),EUL(6),                                   &
                                    DCCR(1,1),DCCR(1,2),DCCR(1,3),                  &
                                    DCCR(2,1),DCCR(2,2),DCCR(2,3),                  &
                                    DCCR(3,1),DCCR(3,2),DCCR(3,3),                  &
                           EULD(1),EULD(2),EULD(3),                                 &
                           EULD(4),EULD(5),EULD(6),3)
            
                      DEP(1)=TEPC(1)-TEPL(1)
                      DEP(2)=TEPC(2)-TEPL(2)
                      DEP(3)=TEPC(3)-TEPL(3)
                      DEP(4)=TEPC(4)-TEPL(4)
                      DEP(5)=TEPC(5)-TEPL(5)
                      DEP(6)=TEPC(6)-TEPL(6)
                      
                        
                    
                    
                    !c      DEU(1)=(1.0D0-PR23)*DEP(1)/(1.0D0-2.0D0*PR23)/(1.0D0+PR23)
                    !c     *            +PR31*DEP(2)/(1.0D0-2.0D0*PR31)/(1.0D0+PR31)
                    !c     *            +PR12*DEP(3)/(1.0D0-2.0D0*PR12)/(1.0D0+PR12)
                    !c      DEU(2)=        PR23*DEP(1)/(1.0D0-2.0D0*PR23)/(1.0D0+PR23)
                    !c     *    +(1.0D0-PR31)*DEP(2)/(1.0D0-2.0D0*PR31)/(1.0D0+PR31)
                    !c     *            +PR12*DEP(3)/(1.0D0-2.0D0*PR12)/(1.0D0+PR12)
                    !c      DEU(3)=        PR23*DEP(1)/(1.0D0-2.0D0*PR23)/(1.0D0+PR23)
                    !c     *            +PR31*DEP(2)/(1.0D0-2.0D0*PR31)/(1.0D0+PR31)
                    !c     *    +(1.0D0-PR12)*DEP(3)/(1.0D0-2.0D0*PR12)/(1.0D0+PR12)
                    !c
                  FAI=1.D0-PR12*PR12-PR23*PR23-PR31*PR31-2.D0*PR12*PR23*PR31        
                  
                 
                  
                  DEU(1)=((1.0D0-PR23*PR23)*DEP(1)                                  &
                        +(PR31*PR23+PR12)*SQRT(ET(2)/ET(1))*DEP(2)                  &
                        +(PR12*PR23+PR31)*SQRT(ET(3)/ET(1))*DEP(3))/FAI
                        
                  DEU(2)=((PR31*PR23+PR12)*SQRT(ET(1)/ET(2))*DEP(1)                 &
                        +(1.0D0-PR31*PR31)*DEP(2)                                   &
                        +(PR12*PR31+PR23)*SQRT(ET(3)/ET(2))*DEP(3))/FAI
                 
                  DEU(3)=((PR12*PR23+PR31)*SQRT(ET(1)/ET(3))*DEP(1)             &
                            +(PR12*PR31+PR23)*SQRT(ET(2)/ET(3))*DEP(2)          &
                            +(1.0D0-PR12*PR12)*DEP(3))/FAI
            
                    ! if (deu(2)<1.0D-15 .and. deu(2)>-1.0D-15)deu(2)=0.0D0
                    ! if (deu(3)<1.0D-15 .and. deu(3)>-1.0D-15)deu(3)=0.0D0
                    
                     
                  EU(1)=EUL(1)+DEU(1)
                  EU(2)=EUL(2)+DEU(2)
                  EU(3)=EUL(3)+DEU(3)
                  EU(4)=TEPC(4)
                  EU(5)=TEPC(5)
                  EU(6)=TEPC(6)
                !!      EU(4)=EUL(4)+DEP(4)
                !!      EU(5)=EUL(5)+DEP(5)
                !!      EU(6)=EUL(6)+DEP(6)
                !!
                
            END SUBROUTINE                              
                                      
        !**************************************************************************************
            SUBROUTINE sCrack1ST1(TSI,DCCCR,NTCRACK,LCRACK)
                ! ********************************************************************
                !*引用名：sCrack1ST
                !*機�E�E��Eび割れ暦ぁE回ある場合�E座標系の計箁E
                !*入力！ETSI           :ひび割れ暦を記録する座標�EトリチE��ス
                !*      DCCCR       :ひびわれ座標系
                !*      NTCRACK :ひび割れ総本数
                !*      LCRACK  :吁E��のひび割れ状況E
                !*出力！EDCCR      :前stepの等価一軸ひずみ�E��E力として�E�E
                !*　　　   PSICR       :ひび割れ座標系での応力
                !********************************************************************
                !IMPLICIT NONE
        
                real(kreal) DAIX(3,3),DAIXL(3,3),DCCCR(3,3)
                real(kreal) DCCCR2(3,3),DCCCR22(3,3),DCCCR3(3,3)
                real(kreal) TSI(6),TEP(3),PSICR(6),PSICR2(6)
            !add 
                  integer II,JJ,KK,INTI,INTJ,INTLCRACK,INTSAME,NTCRACK,LCRACK(3)
                  real(kreal) DXY,DXYZ,DYZ,DZX,PSICRDL
            
                  CALL CHG3 (PSICR(1),PSICR(2),PSICR(3),                            &
                             PSICR(4),PSICR(5),PSICR(6),                            &
                           DCCCR(1,1),DCCCR(1,2),DCCCR(1,3),                        &
                           DCCCR(2,1),DCCCR(2,2),DCCCR(2,3),                        &
                           DCCCR(3,1),DCCCR(3,2),DCCCR(3,3),                        &
                               TSI(1),TSI(2),TSI(3),                                &
                               TSI(4),TSI(5),TSI(6),1)
            
                IF(NTCRACK.GE.2) RETURN
                !!
                !!   *** CRACK IN 1 DIRECTION (MODIFYING DCCCR ON 2 real(kreal)) ***
                !!
                  IF(NTCRACK.EQ.1.AND.LCRACK(1).EQ.1) THEN
                  CALL PRIN (PSICR(1),PSICR(2),PSICR(3),                         &     !PSICR(1)
                             PSICR(4),PSICR(5),PSICR(6),                      &  !PSICR(4),PSICR(6)
                            PSICR2(1),PSICR2(2),PSICR2(3),                          &
                            PSICR2(4),PSICR2(5),PSICR2(6),                          &
                            DCCCR2(1,1),DCCCR2(1,2),DCCCR2(1,3),                    &
                            DCCCR2(2,1),DCCCR2(2,2),DCCCR2(2,3),                    &
                            DCCCR2(3,1),DCCCR2(3,2),DCCCR2(3,3))
            
                PSICRDL=PSICR(1)
                  ENDIF
                !!
                !!   *** CRACK IN 2 DIRECTION(MODIFYING DCCCR ON 2 real(kreal)) ***
                !!
                  IF(NTCRACK.EQ.1.AND.LCRACK(2).EQ.1) THEN
                  CALL PRIN (PSICR(1),PSICR(2),PSICR(3),                         &   !PSICR(2)
                           PSICR(4),PSICR(5),PSICR(6),                        &   !PSICR(4),PSICR(5)
                           PSICR2(1),PSICR2(2),PSICR2(3),                           &
                           PSICR2(4),PSICR2(5),PSICR2(6),                           &
                           DCCCR2(1,1),DCCCR2(1,2),DCCCR2(1,3),                     &
                           DCCCR2(2,1),DCCCR2(2,2),DCCCR2(2,3),                     &
                           DCCCR2(3,1),DCCCR2(3,2),DCCCR2(3,3))
            
                PSICRDL=PSICR(2)
                  ENDIF
                !!
                !!   *** CRACK IN 3 DIRECTION(MODIFYING DCCCR ON 2 real(kreal)) ***
                !!
                  IF(NTCRACK.EQ.1.AND.LCRACK(3).EQ.1) THEN
                  CALL PRIN (PSICR(1),PSICR(2),PSICR(3),                         &   !PSICR(3)
                             PSICR(4),PSICR(5),PSICR(6),                      &   !PSICR(5),PSICR(6)
                            PSICR2(1),PSICR2(2),PSICR2(3),                          &
                            PSICR2(4),PSICR2(5),PSICR2(6),                          &
                            DCCCR2(1,1),DCCCR2(1,2),DCCCR2(1,3),                    &
                            DCCCR2(2,1),DCCCR2(2,2),DCCCR2(2,3),                    &
                            DCCCR2(3,1),DCCCR2(3,2),DCCCR2(3,3))
        
                    PSICRDL=PSICR(3)
                  ENDIF
            
                INTSAME=0
                DO INTI=1,3
                    PSICR(INTI)=PSICR2(INTI)

                        DXY=    DCCCR2(INTI,1)*DCCCR2(INTI,2)
                        DYZ=    DCCCR2(INTI,2)*DCCCR2(INTI,3)
                        DZX=    DCCCR2(INTI,3)*DCCCR2(INTI,1)
                        DXYZ=   DCCCR2(INTI,1)+DCCCR2(INTI,2)+DCCCR2(INTI,3)
                    IF(DXY.EQ.0.0.AND.DYZ.EQ.0.0.AND.DZX.EQ.0.0.AND.DXYZ.EQ.1.0)THEN
                        INTLCRACK=INTI
                        INTSAME=INTSAME+1
                    ENDIF
                ENDDO
            
                IF(INTSAME.NE.3) THEN
                    PSICR(INTLCRACK)=PSICRDL
                END IF
            
                  DO 2305 II=1,3
                  DO 2305 JJ=1,3
                    DCCCR3(II,JJ)= 0.0
                  DO 2305 KK=1,3
                    DCCCR3(II,JJ)=DCCCR3(II,JJ)+                            &
                                  DCCCR22(II,KK)*DCCCR(KK,JJ)
             2305 CONTINUE
     
                DO intI=1,3
                DO intJ=1,3
                    DCCCR(intI,intJ) = DCCCR3(intI,intJ)
                ENDDO
                ENDDO
            
            END SUBROUTINE                               
                                      
        !**************************************************************************************
            SUBROUTINE sCrack1ST(TSI,TEP,DCCCR,PSICR,NTCRACK,LCRACK)
        
                
                  !IMPLICIT NONE
                !add 
                  integer II,JJ,KK,intI,intJ,INTLCRACK,INTSAME,NTCRACK,LCRACK(3)
                  real(kreal) DXY,DXYZ,DZX,EUDL,DYZ
            
                real(kreal) DAIX(3,3),DAIXL(3,3),DCCCR(3,3),DCCRTEMP(3,3)
                real(kreal) DCCCR2(3,3),DCCCR22(3,3),DCCCR3(3,3)
                real(kreal) TSI(6),TEP(6),PSICR(6)
                real(kreal) EU(6),EU2(6)
            
                !c  RETURN
                GOTO 8
                  CALL CHG3 (EU(1),EU(2),EU(3),                             &
                             EU(4),EU(5),EU(6),                             &
                           DCCCR(1,1),DCCCR(1,2),DCCCR(1,3),                &
                           DCCCR(2,1),DCCCR(2,2),DCCCR(2,3),                &
                           DCCCR(3,1),DCCCR(3,2),DCCCR(3,3),                &
                            TEP(1),TEP(2),TEP(3),                           &
                            TEP(4),TEP(5),TEP(6),3)

                DO intI=1,3
                DO intJ=1,3
                    DCCCR2(intI,intJ)=0.0
                ENDDO
                ENDDO
                DCCCR2(1,1)=1.0
                DCCCR2(2,2)=1.0
                DCCCR2(3,3)=1.0
            
                IF(NTCRACK.GE.2) GOTO 8
                !!
                !!   *** CRACK IN 1 DIRECTION (MODIFYING DCCCR ON 2 DIMENSION) ***
                !!
                  IF(NTCRACK.EQ.1.AND.LCRACK(1).EQ.1) THEN
                  CALL PRIN (EU(1),EU(2),EU(3),                              &   !EU(1)=0.D0
                            EU(4),EU(5),EU(6),                                        &   !EU(4)=0.D0,EU(6)=0.D0
                            EU2(1),EU2(2),EU2(3),                                   &
                            EU2(4),EU2(5),EU2(6),                                   &
                            DCCCR2(1,1),DCCCR2(1,2),DCCCR2(1,3),                    &
                            DCCCR2(2,1),DCCCR2(2,2),DCCCR2(2,3),                    &
                            DCCCR2(3,1),DCCCR2(3,2),DCCCR2(3,3))
            
                EUDL=EU(1)
                  ENDIF
                !!
                !!   *** CRACK IN 2 DIRECTION(MODIFYING DCCCR ON 2 DIMENSION) ***
                !!
                  IF(NTCRACK.EQ.1.AND.LCRACK(2).EQ.1) THEN
                  CALL PRIN (EU(1),EU(2),EU(3),                                      &   !EU(2)
                              EU(4),EU(5),EU(6),                                      &   !EU(4),EU(5)
                             EU2(1),EU2(2),EU2(3),                               &
                             EU2(4),EU2(5),EU2(6),                               &
                             DCCCR2(1,1),DCCCR2(1,2),DCCCR2(1,3),                &
                             DCCCR2(2,1),DCCCR2(2,2),DCCCR2(2,3),                &
                             DCCCR2(3,1),DCCCR2(3,2),DCCCR2(3,3))
            
                EUDL=EU(2)
                  ENDIF
                !!
                !!   *** CRACK IN 3 DIRECTION(MODIFYING DCCCR ON 2 DIMENSION) ***
                !!
                  IF(NTCRACK.EQ.1.AND.LCRACK(3).EQ.1) THEN
                  CALL PRIN (EU(1),EU(2),EU(3),                  &   !EU(3)
                               EU(4),EU(5),EU(6),                            &   !EU(5),EU(6)
                               EU2(1),EU2(2),EU2(3),                        &
                               EU2(4),EU2(5),EU2(6),                        &
                               DCCCR2(1,1),DCCCR2(1,2),DCCCR2(1,3),         &
                               DCCCR2(2,1),DCCCR2(2,2),DCCCR2(2,3),         &
                               DCCCR2(3,1),DCCCR2(3,2),DCCCR2(3,3))
        
                EUDL=EU(3)
                  ENDIF
            
                INTSAME=0
                DO INTI=1,3
                    EU(INTI)=EU2(INTI)

                        DXY=    DCCCR2(INTI,1)*DCCCR2(INTI,2)
                        DYZ=    DCCCR2(INTI,2)*DCCCR2(INTI,3)
                        DZX=    DCCCR2(INTI,3)*DCCCR2(INTI,1)
                        DXYZ=   DCCCR2(INTI,1)+DCCCR2(INTI,2)+DCCCR2(INTI,3)
                    IF(DXY.EQ.0.0.AND.DYZ.EQ.0.0.AND.DZX.EQ.0.0.AND.DXYZ.EQ.1.0)THEN
                        INTLCRACK=INTI
                        INTSAME=INTSAME+1
                    ENDIF
                ENDDO
            
                IF(INTSAME.NE.3) THEN
                    EU(INTLCRACK)=EUDL
                END IF
            
                IF(INTSAME.GE.3) THEN
                    CONTINUE
                END IF
            
                  DO 2305 II=1,3
                  DO 2305 JJ=1,3
                    DCCCR3(II,JJ)= 0.0
                  DO 2305 KK=1,3
                    DCCCR3(II,JJ)=DCCCR3(II,JJ)+                    &
                                  DCCCR2(II,KK)*DCCCR(KK,JJ)
             2305 CONTINUE
            
                DO intI=1,3
                DO intJ=1,3
                    IF(DCCCR(intI,intJ).NE.DCCCR3(intI,intJ))THEN
                        CONTINUE
                    ENDIF
                    DCCCR(intI,intJ)=DCCCR3(intI,intJ)
                ENDDO
                ENDDO
            
            
                8 CONTINUE
                  CALL CHG3 (PSICR(1),PSICR(2),PSICR(3),            &
                             PSICR(4),PSICR(5),PSICR(6),                     &
                             DCCCR(1,1),DCCCR(1,2),DCCCR(1,3),               &
                             DCCCR(2,1),DCCCR(2,2),DCCCR(2,3),               &
                             DCCCR(3,1),DCCCR(3,2),DCCCR(3,3),               &
                              TSI(1),TSI(2),TSI(3),                           &
                              TSI(4),TSI(5),TSI(6),1)
                !!
                !c    8 CONTINUE
                !!
                RETURN
            
        
               
            END SUBROUTINE
    
        !************************************************************************************** 
            SUBROUTINE sCheckCMODEL(IROT,CCRACK,EU,ECR,ECCR,ECU,EEN,                    &
                                    LCRACK,NTCRACK,NTFAIL,NTELS,LIUL,NUMCRACK)
                        !********************************************************************
                    !*引用名：sCheckCrack            added by Kengo(2003.4.28)
                    !*機�E�E��Eび割れをチェチE��する
                    !*入力！EEU    :コンクリート�E三軸応力状態を老E�Eした等価一軸ひずみ
                    !*      SN  :コンクリート�E主応力
                    !*      IUL :
                    !*出力！ENTCRACK   :ひび割れ総数
                    !*      NUMCRACK:ひび割れ�E方向番号
                    !*      LIUL    :ひび割れ�E持E��E
                    !*目皁E��E
                    !********************************************************************
                
                  !IMPLICIT NONE
            
                real(kreal) CCRACK(3),EU(6),ECCR(3),EEN(3)
                  real(kreal) ECR,ECU
                  integer NTCRACK,NTFAIL,NTELS,NUMCRACK,IROT,LIUL,LCRACK(3)
            
                !!  IF(IROT.EQ.0) RETURN
                CALL sCrackChk(CCRACK,EU,ECR,ECCR,ECU,EEN,                              &
                               LCRACK,NTCRACK,NTFAIL,NTELS,NUMCRACK)
            
            END SUBROUTINE

        !**************************************************************************************
            SUBROUTINE sCoordCHG1(DCCRCHG,DAIX,DAIXL,DCCR,DCCRL,DCC,DCCL,               &
                                  TSI,TEP,PSICR,M,I,J,K,NTCRACK,LCRACK,         &
                                  IJNC,IACTIVEC,IACTIVECOLD,ICOORDCHGED)
                !********************************************************************
                !*引用名：sCoordCHG1
                !*機�E�E��Eび割れ暦ぁE回ある場合�E座標系の計箁E
                !*入力！EDCCRCHG   :ひび割れ暦を記録する座標�EトリチE��ス
                !*      DCCR        :ひびわれ座標系
                !*      DCCRL       :前�Eびわれ座標系
                !*      DCC         :主軸座標系
                !*      DCCL        :前主軸座標系
                !*      NTCRACK :ひび割れ総本数
                !*      LCRACK  :吁E��のひび割れ状況E
                !*      IJNC        :座標系を変換されたかどぁE��の持E��E
                !*      IACTIVEC:アクチE��ブする座標系の番号
                !*出力！EDCCR      :前stepの等価一軸ひずみ�E��E力として�E�E
                !*      DAIX        :コンクリート�E応力増�E�E��Eび割れ軸�E�E
                !*      DAIXL       :前stepコンクリート�E応力�E��Eび割れ軸�E�E
                !*　　　   PSICR       :ひび割れ座標系での応力
                !********************************************************************
        
              !IMPLICIT NONE
        
                real(kreal) DAIX(3,3),DAIXL(3,3)
                real(kreal) DCC(3,3),DCCL(3,3),DCCR(3,3),DCCRL(3,3)
                real(kreal) DCCRCHG(3,3,3),DCCRCHGL(3,3,3)
                real(kreal) PSICR(6),TSI(6)
                ! add 
                  integer INTM,INTN,IACTIVEC,IACTIVECOLD,IJNC
                  integer I,J,K,M,NTCRACK,ICOORDCHGED,LCRACK(3)
                  real(kreal) TEP(6)
            
                !!  ***[古ぁE�Eび割れを用ぁE��場吁E***
                IF(IJNC.EQ.0) THEN
                
                    DO intM=1,3
                    DO intN=1,3
                        DCCR(intM,intN)=DCCRCHG(1,intM,intN)
                    ENDDO
                    ENDDO
            
                !!      ***[一軸固定�Eび割れモチE��の適用]***
                    IF(NTCRACK.EQ.1) THEN
                        CALL sCrack1ST(TSI,TEP,DCCR,PSICR,NTCRACK,LCRACK)
                !!          CALL sCrack1STTemp(TSI,TEP,PSICR,DCCR,NTCRACK,LCRACK)
                    ENDIF
            
                !!      ***[ひび割れ暦の記録]***
                    IACTIVEC=1
                    IF(IACTIVECOLD.EQ.0)THEN
                        CALL sTransArrayreal33(DCCL,DCCRL)
                    ENDIF
                    CALL sTransArrayreal33(DCCR,DAIX)
                    CALL sTransArrayreal33(DCCRL,DAIXL)
                !!      CALL sRembCrack(DCCR,DCCRL,DCCRCHG,DCCRCHGL,1)
                    RETURN
                ENDIF

               !!   ***[新たなひび割れが入る場吁E***
                CALL sTransArrayreal33(DCC,DCCR)
                IF(IACTIVECOLD.EQ.0)THEN
                    CALL sTransArrayreal33(DCCL,DCCRL)  
                ENDIF
                CALL sTransArrayreal33(DCC,DAIX)
                CALL sTransArrayreal33(DCCL,DAIXL)
                !!  ***[ひび割れ暦の記録]***
                ICOORDCHGED=2
                IACTIVEC=2
                CALL sRembCrack(DCCR,DCCRL,DCCRCHG,DCCRCHGL,2)
        
                RETURN
            END SUBROUTINE   
        !**************************************************************************************
            SUBROUTINE sCoordCHG2(DCCRCHG,DAIX,DAIXL,DCC,DCCL,DCCR,DCCRL,               &
                                      TSI,TEP,PSICR,M,I,J,K,NTCRACK,LCRACK,         &
                                      IJNC,IACTIVEC,IACTIVECOLD,ICOORDCHGED)
                    ! ********************************************************************
                    !*引用名：sCoordCHG2
                    !*機�E�E��Eび割れ暦ぁE回ある場合�E座標系の計箁E
                    !*入力！EDCCRCHG   :ひび割れ暦を記録する座標�EトリチE��ス
                    !*      DCCR        :ひびわれ座標系
                    !*      DCCRL       :前�Eびわれ座標系
                    !*      NTCRACK :ひび割れ総本数
                    !*      LCRACK  :吁E��のひび割れ状況E
                    !*      IJNC        :座標系を変換されたかどぁE��の持E��E
                    !*      IACTIVEC:アクチE��ブする座標系の番号
                    !*出力！EDCCR      :前stepの等価一軸ひずみ�E��E力として�E�E
                    !*      DAIX        :コンクリート�E応力増�E�E��Eび割れ軸�E�E
                    !*      DAIXL       :前stepコンクリート�E応力�E��Eび割れ軸�E�E
                    !*　　　   PSICR       :ひび割れ座標系での応力
                    !********************************************************************
                      !IMPLICIT NONE
            
                real(kreal) DAIX(3,3),DAIXL(3,3),DCCRACTIVE(3,3)
                real(kreal) DCC(3,3),DCCL(3,3),DCCR(3,3),DCCRL(3,3)
                real(kreal) DCCRCHG(3,3,3),DCCRCHGL(3,3,3)
                real(kreal) PSICR(6),TSI(6)
                ! add
                  integer INTM,INTN,IACTIVEC,IACTIVECOLD,IJNC
                  integer I,J,K,M,NTCRACK,ICOORDCHGED,LCRACK(3)
                  real(kreal) TEP(6)

                IF(IJNC.EQ.0) THEN
                    DO intM=1,3
                    DO intN=1,3
                        DCCR(intM,intN)=DCCRCHG(IACTIVEC,intM,intN)
                    ENDDO
                    ENDDO
            
                    IF(IACTIVECOLD.EQ.0)THEN
                        CALL sTransArrayreal33(DCCL,DCCRL)
                    ENDIF
                    IF(NTCRACK.EQ.1) THEN
                        CALL sCrack1ST(TSI,TEP,DCCR,PSICR,NTCRACK,LCRACK)
                    !!  CALL sCrack1STTemp(TSI,TEP,PSICR,DCCR,NTCRACK,LCRACK)
                    ENDIF
            
                    !!      ***[ひび割れ暦の記録]***
                    CALL sTransArrayreal33(DCCR,DAIX)
                    CALL sTransArrayreal33(DCCRL,DAIXL)
                    !!      CALL sRembCrack(DCCR,DCCRL,DCCRCHG,DCCRCHGL,IACTIVEC)
                    RETURN
                ENDIF
            
                    !!  ***[新たなひび割れが入る場吁E***
                CALL sTransArrayreal33(DCC,DCCR)
                    IF(IACTIVECOLD.EQ.0)THEN
                        CALL sTransArrayreal33(DCCL,DCCRL)
                    ENDIF
                CALL sTransArrayreal33(DCC,DAIX)
                CALL sTransArrayreal33(DCCL,DAIXL)
                    !!  ***[ひび割れ暦の記録]***
                ICOORDCHGED=3
                IACTIVEC=3
                CALL sRembCrack(DCCR,DCCRL,DCCRCHG,DCCRCHGL,3)
            
            END SUBROUTINE  
        
        !************************************************************************************** 
            SUBROUTINE sCoordCHG3(DCCRCHG,DAIX,DAIXL,DCC,DCCL,DCCR,DCCRL,               &
                                  TSI,TEP,PSICR,M,I,J,K,NTCRACK,LCRACK,                 &
                                  IJNC,IACTIVEC,IACTIVECOLD)
                    !********************************************************************
                    !*引用名：sCoordCHG3
                    !*機�E�E��Eび割れ暦ぁE回ある場合�E座標系の計箁E
                    !*入力！EDCCRCHG   :ひび割れ暦を記録する座標�EトリチE��ス
                    !*      DCCR        :ひびわれ座標系
                    !*      DCCRL       :前�Eびわれ座標系
                    !*      NTCRACK :ひび割れ総本数
                    !*      LCRACK  :吁E��のひび割れ状況E
                    !*      IJNC        :座標系を変換されたかどぁE��の持E��E
                    !*      IACTIVEC:アクチE��ブする座標系の番号
                    !*出力！EDCCR      :前stepの等価一軸ひずみ�E��E力として�E�E
                    !*      DAIX        :コンクリート�E応力増�E�E��Eび割れ軸�E�E
                    !*      DAIXL       :前stepコンクリート�E応力�E��Eび割れ軸�E�E
                    !*　　　   PSICR       :ひび割れ座標系での応力
                    !********************************************************************
                  !IMPLICIT NONE
            
                real(kreal) DAIX(3,3),DAIXL(3,3),DCCRACTIVE(3,3)
                real(kreal) DCC(3,3),DCCL(3,3),DCCR(3,3),DCCRL(3,3)
                real(kreal) DCCRCHG(3,3,3),DCCRCHGL(3,3,3)
                real(kreal) PSICR(6),TSI(6)
                ! add 
                  integer INTM,INTN,IACTIVEC,IACTIVECOLD,IJNC
                  integer I,J,K,M,NTCRACK,ICOORDCHGED,LCRACK(3)
                  real(kreal) TEP(6)
            
                    !!  ***[角度が小さぁE��座標系の選択]***
                IF(IACTIVECOLD.EQ.0)THEN
                    CALL sTransArrayreal33(DCCL,DCCRL)
                ENDIF
                CALL sCrack1ST(TSI,TEP,DCCR,PSICR,NTCRACK,LCRACK)
                !!      CALL sCrack1STTemp(TSI,TEP,PSICR,DCCR,NTCRACK,LCRACK)
          
                    DO intM=1,3
                    DO intN=1,3
                        DCCR(intM,intN)=DCCRCHG(IACTIVEC,intM,intN)
                    ENDDO
                    ENDDO
            
                !!      ***[ひび割れ暦の記録]***
                CALL sTransArrayreal33(DCCR,DAIX)
                CALL sTransArrayreal33(DCCRL,DAIXL)
                !!  CALL sRembCrack(DCCR,DCCRL,DCCRCHG,DCCRCHGL,IACTIVEC)
                RETURN
                !c  ENDIF
            
            END SUBROUTINE                 
        !************************************************************************************** 
            SUBROUTINE sChangeAxisTurn(DCC,DCCR,IACT,SNDM,IR)

                real (kreal)DCC(3,3),DCCR(3,3),SN(3),SNDM(3)
                real (kreal)  DCCTemp1(3,3),DCCTemp2(3,3)
                integer iact,ir,inti

                IF(IACT.EQ.0)THEN
                    CALL sTransArrayReal33(DCC,DCCTemp1)
                ELSE
                    CALL sTransArrayReal33(DCCR,DCCTemp1)
                ENDIF

                CALL sTransArrayReal33(DCCTemp1,DCCTemp2)

                IF(IR.EQ.0) THEN
            !  *** 1 ==> 1  &  2 ==> 2  &  3 ==> 3 ***
                ELSE
                CONTINUE
                ENDIF

              IF(IR.EQ.1) THEN
            !   *** 2 <==> 3 ***
                DO intI=1,3
                    DCCTemp2(2,intI)=DCCTemp1(3,intI)
                    DCCTemp2(3,intI)=DCCTemp1(2,intI)
                ENDDO

                CALL XRVS1(SN(2),SN(3))
              ENDIF

              IF(IR.EQ.2) THEN
        !   *** 1 <==> 2 ***
                DO intI=1,3
                    DCCTemp2(1,intI)=DCCTemp1(2,intI)
                    DCCTemp2(2,intI)=DCCTemp1(1,intI)
                ENDDO
            
                CALL XRVS1(SN(1),SN(2))
         
            ENDIF

               IF(IR.EQ.3) THEN
        !  *** 1 ==> 2  &  2 ==> 3  &  3 ==> 1 ***
                DO intI=1,3
                    DCCTemp2(1,intI)=DCCTemp1(2,intI)
                    DCCTemp2(2,intI)=DCCTemp1(3,intI)
                    DCCTemp2(3,intI)=DCCTemp1(1,intI)
                ENDDO
                CALL XRVS2(SN(1),SN(2),SN(3))
               ENDIF

               IF(IR.EQ.4) THEN
        !   *** 1 <==> 3 ***
                DO intI=1,3
                    DCCTemp2(1,intI)=DCCTemp1(3,intI)
                    DCCTemp2(3,intI)=DCCTemp1(1,intI)
                ENDDO
                CALL XRVS1(SN(1),SN(3))
               ENDIF

               IF(IR.EQ.5) THEN
        !   *** 1 ==> 3  &  2 ==> 1  &  3 ==> 2 ***
                DO intI=1,3
                    DCCTemp2(1,intI)=DCCTemp1(3,intI)
                    DCCTemp2(2,intI)=DCCTemp1(1,intI)
                    DCCTemp2(3,intI)=DCCTemp1(2,intI)
                ENDDO
                CALL XRVS3(SN(1),SN(2),SN(3))

               ENDIF
            
            IF(IACT.EQ.0)THEN
                CALL sTransArrayReal33(DCCTemp2,DCC)
            ELSE
                CALL sTransArrayReal33(DCCTemp2,DCCR)
              ENDIF

            END
        !************************************************************************************** 
            SUBROUTINE CRACK_1(ICIC6,FT,EO,EU,ETP,DC_1,DC_2,DCC,DCCR)

                !**** ICIC6のひび割れモチE��で�E�E新固定１方向�Eび割れモチE��となる　sun   
              
        
            real (kreal) DCC(3,3),DCCR(3,3),DC_1(3,3),DC_2(3,3),ETP(3),EU(6)
            integer ICIC6,i,j
            real (kreal) ft,eo
            real EUCR_,EUCR1,EUCR2,EUCR3
            
            IF(ICIC6 /= 3) GOTO 10

            IF(DC_1(1,1)==-2.)THEN
                EUCR_= FT/EO * (-0.0)
                EUCR1=ETP(1)-EUCR_ 
                EUCR2=ETP(2)-EUCR_ 
                EUCR3=ETP(3)-EUCR_
                IF(EU(1) > EUCR1 .OR. EU(2) > EUCR2 .OR. EU(3) > EUCR3)THEN
                    DO I=1,3 ; DO J=1,3 ; DC_1(I,J)=DCC(I,J) ; ENDDO;ENDDO
                    CALL InputDirection(DC_1,DC_2)
                ENDIF
            ENDIF

                IF(DC_1(1,1)/=-2.)THEN
                    DO I=1,3;DO J=1,3
                        DCC(I,J)=DC_1(I,J) ; DCCR(I,J)=DC_1(I,J)
                    ENDDO;ENDDO
                ENDIF

        10  CONTINUE
            END
        !************************************************************************************** 
                            
            SUBROUTINE InputDirection(DC_1,DC_2)
        
                real (kreal) DC_1(3,3),DC_2(3,3)
                real (kreal) c
                integer i,j
                
                    C=0.
                    DO I=1,3
                        DO J=1,3
                            !=!+DC_2(I,J)**2
                        ENDDO
                    ENDDO
                IF( C > 0.0 )THEN
                    DO I=1,3
                        DO J=1,3
                        DC_1(I,J)=DC_2(I,J)
                        ENDDO
                    ENDDO
                ENDIF
            10  END
        !**************************************************************************************             

            SUBROUTINE sSCEC(SND,SNCRD,EU,IUL,IROT,ICR,FC,ECU,EO,FT,SC,EC,  &
                    EEN,EBU,ICIC1,ECCR,NST,CCRACK,EPCU,ICONCYC,SCHJ)
     
    
                ! ********************************************************************
        ! *引用名：sSCEC                             added by KENGO(2003.5.15)
        ! *機�E�E�SC,ECを求めめE
        ! *入力！E     :
        ! *出力！E     :
        ! *参老E��E
        ! ********************************************************************
        ! *引数�E�SC     :コンクリート�E三軸応力状態を老E�Eした最大主応力�E�破壊曲面との接点の値�E�E
        ! *　　　EC     :コンクリート要素の3軸応力状態を老E�Eした最大圧縮応力時主ひずみ
        ! *　　　FC     :コンクリート要素の一軸圧縮強度
        ! *　　　ECU    :コンクリート要素の圧縮強度時�Eずみ
        ! *　　　EU     :コンクリート�E三軸応力状態を老E�Eした等価一軸ひずみ
        ! *　　　SN     :コンクリート�E主応力�E�ESI�E�E
        ! *　　　IUL    :コンクリート�E履歴持E��E
        ! *　　　LCRACK :弾性=0�E��Eび割めE1�E�圧壁E2
        ! *　　　NCR0   :弾性�E�　　LCRACK=0の数
        ! *　　　NCR1   :ひび割れ，LCRACK=1の数
        ! *　　　NCR2   :圧壊，　　LCRACK=2の数
        ! *********************************************************************

        real (kreal) :: SC(3),EC(3),SCC(3),SN(3),SNN(3),SCL(3),SNCR(3),EPT(3),  &
                EUU(3),EU(3),ECL(3),ECC(3),EEN(3),  &
                ECCR(3),CCRACK(3),STMAX(3),ETMAX(3)
        integer     ::  LPOSNEG(3),LCOMTEN(3),IUL(3),IROT,LCRACK(3)
        
        real(kreal)  ::  FT,EO,ECU,FC,EBU,EPCU
        real(kreal)  :: ALT,ALC,BET,BEC,FF30,FF65_1,FCC,SSS,SS
        real(kreal)  :: EUK1,EUK2,EUK3,ECR
        
        real (kreal) :: SND(3),SNCRD(3),SCHJ(3)
        real (kreal) :: FCCC
        real (kreal) :: SCDMY(3),ECDMY(3),eudm(6),SNK1,SNK2,SNK3
        real (kreal) :: SNKCR1,SNKCR2,SNKCR3,ECUUU,D
        real (kreal) :: AECU,AET,EET,ECUU,EE,EEE
        integer ::  I,ICR,ICIC1,NST,ICONCYC
        integer :: intI,NUMCRACK,NCOM,NTEN,NCR0,NCR1,NCR2,NNEG,NPOS
!       KUPHER
        ALT=ABS(FT/FC)
        ALC=1.15D0
        BET=ABS(FT/EO/ECU)
        BEC=1.7512D0
        FF30=1.0
        FF65_1=0.85
! C
! C      SCHICKERT
! C      ALT=ABS(FT/FC)
! C      ALC=1.21D0
! C      BET=ABS(FT/EO/ECU)
! C      BEC=1.3D0
! Cc
        DO intI=1,3
            SCDMY(intI)=SC(intI)
            ECDMY(intI)=EC(intI)
        ENDDO
! C
! C IF(SND(1).EQ.0.0.AND.SND(2).EQ.0.0.AND.SND(3).EQ.0.0) THEN
! C     DO intI=1,3
! C         SN(intI)=0.01
! C     ENDDO
! C ELSE
        DO intI=1,3
            SN(intI)=SND(intI)
        ENDDO
! C ENDIF
! C
! C IF(SNCRD(1).EQ.0.0.AND.SNCRD(2).EQ.0.0.AND.SNCRD(3).EQ.0.0) THEN
! C     DO intI=1,3
! C         SNCR(intI)=0.01
! C     ENDDO
! C ELSE
        DO intI=1,3
            SNCR(intI)=SNCRD(intI)
        ENDDO
! C ENDIF

        EUK1=EU(1)
        EUK2=EU(2)
        EUK3=EU(3)
        ECR=FT/EO

        CALL sCrackChk(CCRACK,EU,ECR,ECCR,ECU,EEN,  &
            LCRACK,NCR1,NCR2,NCR0,NUMCRACK)

        DO intI=1,3
            IF(IUL(intI).NE.4.AND.IUL(intI).NE.5.AND.IUL(intI).LT.40)THEN
                IF(LCRACK(intI).EQ.1)THEN
                    LCRACK(intI)=0
                    NCR1=NCR1-1
                    NCR0=NCR0+1
                ENDIF
            ENDIF
        ENDDO

        IF(IROT.EQ.0)     GOTO 305
        IF(NCR1.EQ.0)     GOTO 305
    ! c IF(ICR.EQ.0)      GOTO 305
        GOTO 310

    305 CONTINUE
        SNK1=SN(1)
        SNK2=SN(2)
        SNK3=SN(3)

! ! C***DECIDE MAXIMUM COMPRESSIVE STRESSES***

        CALL sComTen  (SN,LCRACK,LCOMTEN,NCOM,NTEN)

    FCCC=FF30*FC

 ! CCCCCCC/////�E�軸弾性/////
    IF (NCR0.EQ.3.and.NCR1.EQ.0.and.NCR2.EQ.0)THEN
        DO intI= 1,3
        SNN(intI)=SN(intI)
        SCL(intI)=SC(intI)
        ENDDO

    CALL XAGRIS(SNN(1),SNN(2),SNN(3),SC(1),SC(2),SC(3),FC,FT,   &
                        ALT,ALC,BET,BEC,1) 
    DO I=1,3 ; SCHJ(I)=SC(I) ; ENDDO
! C�E�－－－－３軸圧縮
        IF (NCOM.EQ.3.and.NTEN.EQ.0)THEN
            DO intI= 1,3
                IF(SC(intI).GE.FCCC) SC(intI)=FCCC
            ENDDO
            SSS=MIN(SC(1),SC(2),SC(3))
            FCC=FF65_1*FC

        ENDIF
! C�E�－－－－２軸圧縮�E�１軸引張
        IF (NCOM.EQ.2.and.NTEN.EQ.1)THEN
            IF(LCOMTEN(1).EQ.1.and.LCOMTEN(2).EQ.-1.and.LCOMTEN(3).EQ.-1)THEN
                SNN(1)=0.0D0
                SNN(2)=SN(2)
                SNN(3)=SN(3)
            CALL XAGRIS (SNN(1),SNN(2),SNN(3),SCC(1),SCC(2),SCC(3),FC,  &
            FT,ALT,ALC,BET,BEC,1) 

                SSS=MIN(SCC(2),SCC(3))
                FCC=FF65_1*SSS
                SS=MIN(SC(2),SC(3))
                    IF(SS.LT.FCC)THEN
                        SC(1)=FT
                    ELSE
                        SC(2)=SCC(2)*FF65_1
                        SC(3)=SCC(3)*FF65_1
                    ENDIF
                    IF(SC(2).GT.FCCC)SC(2)=FCCC
                    IF(SC(3).GT.FCCC)SC(3)=FCCC
        ENDIF

    IF(LCOMTEN(1).EQ.-1.and.LCOMTEN(2).EQ.1.and.LCOMTEN(3).EQ.-1)THEN
                SNN(1)=SN(1)
                SNN(2)=0.0D0
                SNN(3)=SN(3)
            CALL XAGRIS (SNN(1),SNN(2),SNN(3),SCC(1),SCC(2),SCC(3),     &
                        FC,FT,ALT,ALC,BET,BEC,1) 

                SSS=MIN(SCC(1),SCC(3))
                FCC=FF65_1*SSS
                SS=MIN(SC(1),SC(3))
                    IF(SS.LT.FCC)THEN
                        SC(2)=FT
                    ELSE
                        SC(1)=SCC(1)*FF65_1
                        SC(3)=SCC(3)*FF65_1
                    ENDIF
                    IF(SC(1).GT.FCCC)SC(1)=FCCC
                    IF(SC(3).GT.FCCC)SC(3)=FCCC
    ENDIF

    IF(LCOMTEN(1).EQ.-1.and.LCOMTEN(2).EQ.-1.and.LCOMTEN(3).EQ.1)THEN
                SNN(1)=SN(1)
                SNN(2)=SN(2)
                SNN(3)=0.0D0
            CALL XAGRIS (SNN(1),SNN(2),SNN(3),SCC(1),SCC(2),SCC(3), &
            FC,FT,ALT,ALC,BET,BEC,1) 
                SSS=MIN(SCC(1),SCC(2))
                FCC=FF65_1*SSS
                SS=MIN(SC(1),SC(2))
                    IF(SS.LT.FCC)THEN
                        SC(3)=FT
                    ELSE
                        SC(1)=SCC(1)*FF65_1
                        SC(2)=SCC(2)*FF65_1
                    ENDIF
                    IF(SC(1).GT.FCCC)SC(1)=FCCC
                    IF(SC(2).GT.FCCC)SC(2)=FCCC
    ENDIF

        ENDIF
! C�E�－－－－１軸圧縮�E�２軸引張
        IF (NCOM.EQ.1.and.NTEN.EQ.2)THEN
            SSS=FC
            FCC=FF65_1*FC

    IF(LCOMTEN(1).EQ.-1.and.LCOMTEN(2).EQ.1.and.LCOMTEN(3).EQ.1)THEN
            SS=SC(1)
        IF(SS.LT.FCC)THEN
            SC(2)=FT
            SC(3)=FT
        ELSE
            SC(1)=FCC
            SC(2)=MAX(SC(2),SC(3))
            SC(3)=MAX(SC(2),SC(3))
        ENDIF
    ENDIF

    IF(LCOMTEN(1).EQ.1.and.LCOMTEN(2).EQ.-1.and.LCOMTEN(3).EQ.1)THEN
            SS=SC(2)
        IF(SS.LT.FCC)THEN
            SC(1)=FT
            SC(3)=FT
        ELSE
            SC(2)=FCC
            SC(1)=MAX(SC(1),SC(3))
            SC(3)=MAX(SC(1),SC(3))
        ENDIF
    ENDIF

    IF(LCOMTEN(1).EQ.1.and.LCOMTEN(2).EQ.1.and.LCOMTEN(3).EQ.-1)THEN
            SS=SC(3)
        IF(SS.LT.FCC)THEN
            SC(1)=FT
            SC(2)=FT
        ELSE
            SC(3)=FCC
            SC(1)=MAX(SC(1),SC(2))
            SC(2)=MAX(SC(1),SC(2))
        ENDIF
    ENDIF

        ENDIF
!C�E�－－－－３軸引張
        IF(NCOM.EQ.0.and.NTEN.EQ.3)THEN
            SSS=FC
            FCC=FF65_1*FC
                DO intI= 1,3
                    SC(intI)=FT
                ENDDO
        ENDIF
        GOTO 350
    ENDIF

!CCCCCCC/////�E�軸弾性�E�１軸圧壁E////
    IF (NCR0.EQ.2.and.NCR1.EQ.0.and.NCR2.EQ.1)THEN
        DO intI= 1,3
        SNN(intI)=SN(intI)
        ENDDO

        DO intI= 1,3
            IF(LCRACK(intI).EQ.2)THEN
                SNN(intI)=0.0D0
                SCL(intI)=SC(intI)
            ENDIF
        ENDDO

    CALL XAGRIS(SNN(1),SNN(2),SNN(3),SC(1),SC(2),SC(3),FC,FT,&
                        ALT,ALC,BET,BEC,1) 
        DO intI= 1,3
            IF(LCRACK(intI).EQ.2)THEN
                SC(intI)=SCL(intI)
            ENDIF
        ENDDO

        SSS=FC
        FCC=FF65_1*FC
!C�E�－－－－１軸圧縮�E�１軸引張
        IF(NCOM.EQ.1.and.NTEN.EQ.1)THEN
            IF(LCRACK(1).EQ.2)THEN
                IF(LCOMTEN(2).EQ.1.and.LCOMTEN(3).EQ.-1)THEN
                    SS=SC(3)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(2)=FT
                    ELSE
                        SC(3)=FCC
                    ENDIF
                ENDIF
                IF(LCOMTEN(2).EQ.-1.and.LCOMTEN(3).EQ.1)THEN
                    SS=SC(2)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(3)=FT
                    ELSE
                        SC(2)=FCC
                    ENDIF
                ENDIF
            ENDIF

            IF(LCRACK(2).EQ.2)THEN
                IF(LCOMTEN(1).EQ.1.and.LCOMTEN(3).EQ.-1)THEN
                    SS=SC(3)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(1)=FT
                    ELSE
                        SC(3)=FCC
                    ENDIF
                ENDIF
                IF(LCOMTEN(1).EQ.-1.and.LCOMTEN(3).EQ.1)THEN
                    SS=SC(1)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(3)=FT
                    ELSE
                        SC(1)=FCC
                    ENDIF
                ENDIF
            ENDIF

            IF(LCRACK(3).EQ.2)THEN
                IF(LCOMTEN(1).EQ.1.and.LCOMTEN(2).EQ.-1)THEN
                    SS=SC(2)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(1)=FT
                    ELSE
                        SC(2)=FCC
                    ENDIF
                ENDIF
                IF(LCOMTEN(1).EQ.-1.and.LCOMTEN(2).EQ.1)THEN
                    SS=SC(1)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(2)=FT
                    ELSE
                        SC(1)=FCC
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
!C�E�－－－－２軸圧縮
        IF(NCOM.EQ.2.and.NTEN.EQ.0)THEN

            IF(LCRACK(1).EQ.2)THEN
                SS=MIN(SC(2),SC(3))
                IF(SS.GT.FCC)THEN
                SC(2)=FCC
                SC(3)=FCC
                ENDIF
                IF(SC(2).GT.FCCC) SC(2)=FCCC
                IF(SC(3).GT.FCCC) SC(3)=FCCC
            ENDIF

                IF(LCRACK(2).EQ.2)THEN
                SS=MIN(SC(1),SC(3))
                IF(SS.GT.FCC)THEN
                SC(1)=FCC
                SC(3)=FCC
                ENDIF
                IF(SC(1).GT.FCCC) SC(1)=FCCC
                IF(SC(3).GT.FCCC) SC(3)=FCCC
            ENDIF

                IF(LCRACK(3).EQ.2)THEN
                SS=MIN(SC(1),SC(2))
                IF(SS.GT.FCC)THEN
                SC(1)=FCC
                SC(2)=FCC
                ENDIF
                IF(SC(1).GT.FCCC) SC(1)=FCCC
                IF(SC(2).GT.FCCC) SC(2)=FCCC
            ENDIF
        ENDIF
!C�E�－－－－２軸引張
        IF(NCOM.EQ.0.and.NTEN.EQ.2)THEN
            DO intI=1,3
                IF(LCRACK(intI).EQ.0) THEN
                    SC(intI)=FT
                ENDIF
            ENDDO
        ENDIF
    GOTO 350
    ENDIF

!CCCCCCC/////�E�軸弾性�E�２軸圧壁E////
    IF (NCR0.EQ.1.and.NCR1.EQ.0.and.NCR2.EQ.2)THEN
        DO intI=1,3
            IF(LCRACK(intI).EQ.0)THEN
                IF(SN(intI).LE.0.0D0) SC(intI)=FC
                IF(SN(intI).GT.0.0D0) SC(intI)=FT
            ENDIF
        ENDDO
        SSS=FC
        FCC=FF65_1*FC
    GOTO 350
    ENDIF

!CCCCCCC/////3軸圧壁E////
    IF (NCR0.EQ.0.and.NCR1.EQ.0.and.NCR2.EQ.3)THEN
        SSS=FC
        FCC=FF65_1*FC
    GOTO 350
    ENDIF

!CCCCCCC/////�E�軸弾性�E�１軸ひび割めE////
    IF (NCR0.EQ.2.and.NCR1.EQ.1.and.NCR2.EQ.0)THEN
        DO intI= 1,3
        SNN(intI)=SN(intI)
        ENDDO

        DO intI= 1,3
            IF(LCRACK(intI).EQ.1)THEN
                SNN(intI)=0.0D0
                SCL(intI)=SC(intI)
            ENDIF
        ENDDO

        CALL XAGRIS(SNN(1),SNN(2),SNN(3),SC(1),SC(2),SC(3),FC,FT,   &
                        ALT,ALC,BET,BEC,1) 
        DO intI= 1,3
            IF(LCRACK(intI).EQ.1)THEN
                SC(intI)=SCL(intI)
            ENDIF
        ENDDO
        SSS=FC
            FCC=FF65_1*FC
            FCCC=FF30*FC
!C�E�－－－－１軸圧縮�E�１軸引張
        IF(NCOM.EQ.1.and.NTEN.EQ.1)THEN
            IF(LCRACK(1).EQ.1)THEN
                IF(LCOMTEN(2).EQ.1.and.LCOMTEN(3).EQ.-1)THEN
                    SS=SC(3)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(2)=FT
                    ELSE
                        SC(3)=FCC
                    ENDIF
                ENDIF
                IF(LCOMTEN(2).EQ.-1.and.LCOMTEN(3).EQ.1)THEN
                    SS=SC(2)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(3)=FT
                    ELSE
                        SC(2)=FCC
                    ENDIF
                ENDIF
            ENDIF

            IF(LCRACK(2).EQ.1)THEN
                IF(LCOMTEN(1).EQ.1.and.LCOMTEN(3).EQ.-1)THEN
                    SS=SC(3)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(1)=FT
                    ELSE
                        SC(3)=FCC
                    ENDIF
                ENDIF
                IF(LCOMTEN(1).EQ.-1.and.LCOMTEN(3).EQ.1)THEN
                    SS=SC(1)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(3)=FT
                    ELSE
                        SC(1)=FCC
                    ENDIF
                ENDIF
            ENDIF

            IF(LCRACK(3).EQ.1)THEN
                IF(LCOMTEN(1).EQ.1.and.LCOMTEN(2).EQ.-1)THEN
                    SS=SC(2)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(1)=FT
                    ELSE
                        SC(2)=FCC
                    ENDIF
                ENDIF
                IF(LCOMTEN(1).EQ.-1.and.LCOMTEN(2).EQ.1)THEN
                    SS=SC(1)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(2)=FT
                    ELSE
                        SC(1)=FCC
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
!C�E�－－－－２軸圧縮
        IF(NCOM.EQ.2.and.NTEN.EQ.0)THEN

            IF(LCRACK(1).EQ.1)THEN
                SS=MIN(SC(2),SC(3))
                IF(SS.GT.FCC)THEN
                SC(2)=FCC
                SC(3)=FCC
                ENDIF
                IF(SC(2).GT.FCCC) SC(2)=FCCC
                IF(SC(3).GT.FCCC) SC(3)=FCCC
            ENDIF

            IF(LCRACK(2).EQ.1)THEN
                SS=MIN(SC(1),SC(3))
                IF(SS.GT.FCC)THEN
                SC(1)=FCC
                SC(3)=FCC
                ENDIF
                IF(SC(1).GT.FCCC) SC(1)=FCCC
                IF(SC(3).GT.FCCC) SC(3)=FCCC
            ENDIF

            IF(LCRACK(3).EQ.1)THEN
                SS=MIN(SC(1),SC(2))
                IF(SS.GT.FCC)THEN
                SC(1)=FCC
                SC(2)=FCC
                ENDIF
                IF(SC(1).GT.FCCC) SC(1)=FCCC
                IF(SC(2).GT.FCCC) SC(2)=FCCC
            ENDIF
        ENDIF
!C�E�－－－－２軸引張
        IF(NCOM.EQ.0.and.NTEN.EQ.2)THEN
            DO intI=1,3
                IF(LCRACK(intI).EQ.0)THEN
                    SC(intI)=FT
                ENDIF
            ENDDO
        ENDIF
    GOTO 350
    ENDIF

!CCCCCCC/////�E�軸弾性�E�１軸圧壊，１軸ひび割めE////
    IF (NCR0.EQ.1.and.NCR1.EQ.1.and.NCR2.EQ.1)THEN
            DO intI=1,3
            IF(LCRACK(intI).EQ.0)THEN
                IF(SN(intI).LE.0.0D0) SC(intI)=FC
                IF(SN(intI).GT.0.0D0) SC(intI)=FT
            ENDIF
        ENDDO
        SSS=FC
        FCC=FF65_1*FC
    GOTO 350
    ENDIF
! c         DO intI=1,3
! c             SNN(intI)=SN(intI)
! c         ENDDO
! c
! c         DO intI=1,3
! c             IF(LCRACK(intI).EQ.1)THEN
! c                 SNN(intI)=0.0D0
! c             ENDIF
! c         ENDDO
! c         
! c         DO intI=1,3
! c             IF(LCRACK(intI).NE.0)THEN
! c                 SCL(intI)=SC(intI)
! c             ENDIF
! c         ENDDO
! c 
! c CALL XAGRIS(SNN(1),SNN(2),SNN(3),SC(1),SC(2),SC(3),FC,FT,
! c     *                   ALT,ALC,BET,BEC,1) 
! c         DO intI=1,3
! c             IF(LCRACK(intI).NE.0)THEN
! c                 SC(intI)=SCL(intI)
! c             ELSE
! c                 IF(SN(intI).LE.0.0D0) SC(intI)=FC
! c                 IF(SN(intI).GT.0.0D0) SC(intI)=FT
! c             ENDIF
! c         ENDDO
! c     SSS=FC
! c     FCC=0.65D0*FC
! c GOTO 350
! c ENDIF
! C
! CCCCCCC/////�E�軸圧壊，１軸ひび割めE////
    IF (NCR0.EQ.0.and.NCR1.EQ.1.and.NCR2.EQ.2)THEN
        SSS=FC
        FCC=FF65_1*FC
    GOTO 350
    ENDIF
! C
! CCCCCCC/////�E�軸弾性�E�２軸ひび割めE////
    IF (NCR0.EQ.1.and.NCR1.EQ.2.and.NCR2.EQ.0)THEN
        DO intI=1,3
            IF(LCRACK(intI).EQ.0)THEN
                IF(SN(intI).LE.0.0D0) SC(intI)=FC
                IF(SN(intI).GT.0.0D0) SC(intI)=FT
            ENDIF
        ENDDO
        SSS=FC
        FCC=FF65_1*FC
    GOTO 350
    ENDIF

! CCCCCCC/////�E�軸圧壊，２軸ひび割めE////
    IF (NCR0.EQ.0.and.NCR1.EQ.2.and.NCR2.EQ.1)THEN
        SSS=FC
        FCC=FF65_1*FC
    GOTO 350
    ENDIF

!CCCCCCC/////�E�軸ひび割めE////
    IF (NCR0.EQ.0.and.NCR1.EQ.3.and.NCR2.EQ.0)THEN
        SSS=FC
        FCC=FF65_1*FC
    GOTO 350
    ENDIF

 310    CONTINUE
! C***CRACK DIRECTION MODEL***
    SNKCR1=SNCR(1)
    SNKCR2=SNCR(2)
    SNKCR3=SNCR(3)

    CALL sCrackChk(CCRACK,EU,ECR,ECCR,ECU,EEN,  &
        LCRACK,NCR1,NCR2,NCR0,NUMCRACK)

    CALL sComTen  (SNCR,LCRACK,LCOMTEN,NCOM,NTEN)
! CCCCCCC
    FCCC=FF30*FC
! CCCCCCC
! C
! CCCCCCC/////�E�軸弾性�E�１軸ひび割めE////
    IF (NCR0.EQ.2.and.NCR1.EQ.1.and.NCR2.EQ.0)THEN
        DO intI= 1,3
        SNN(intI)=SNCR(intI)
        ENDDO

        DO intI= 1,3
            IF(LCRACK(intI).EQ.1)THEN
                SNN(intI)=0.0D0
                SCL(intI)=SC(intI)
            ENDIF
        ENDDO

        CALL XAGRIS(SNN(1),SNN(2),SNN(3),SC(1),SC(2),SC(3),FC,FT,   &
                       ALT,ALC,BET,BEC,1) 
        DO intI= 1,3
            IF(LCRACK(intI).EQ.1)THEN
                SC(intI)=SCL(intI)
            ENDIF
        ENDDO
        SSS=FC
            FCC=FF65_1*FC
            FCCC=FF30*FC
! C�E�－－－－１軸圧縮�E�１軸引張
        IF(NCOM.EQ.1.and.NTEN.EQ.1)THEN
            IF(LCOMTEN(2).EQ.1.and.LCOMTEN(3).EQ.-1)THEN
                SS=SC(3)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(2)=FT
                    ELSE
                        SC(3)=FCC
                    ENDIF
            ENDIF
            IF(LCOMTEN(2).EQ.-1.and.LCOMTEN(3).EQ.1)THEN
                SS=SC(2)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(3)=FT
                    ELSE
                        SC(2)=FCC
                    ENDIF
            ENDIF
            IF(LCOMTEN(1).EQ.1.and.LCOMTEN(3).EQ.-1)THEN
                SS=SC(3)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(1)=FT
                    ELSE
                        SC(3)=FCC
                    ENDIF
            ENDIF
            IF(LCOMTEN(1).EQ.-1.and.LCOMTEN(3).EQ.1)THEN
                SS=SC(1)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(3)=FT
                    ELSE
                        SC(1)=FCC
                    ENDIF
            ENDIF
            IF(LCOMTEN(1).EQ.-1.and.LCOMTEN(2).EQ.1)THEN
                SS=SC(1)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(2)=FT
                    ELSE
                        SC(1)=FCC
                    ENDIF
            ENDIF
            IF(LCOMTEN(1).EQ.1.and.LCOMTEN(2).EQ.-1)THEN
                SS=SC(2)
                    IF(SS.LT.FCC.and.SS.GE.FC)THEN
                        SC(1)=FT
                    ELSE
                        SC(2)=FCC
                    ENDIF
            ENDIF
        ENDIF
! C�E�－－－－２軸圧縮
        IF(NCOM.EQ.2.and.NTEN.EQ.0)THEN
            IF(LCRACK(1).EQ.1)THEN
                SS=MIN(SC(2),SC(3))
                IF(SS.GT.FCC)THEN
                SC(2)=FCC
                SC(3)=FCC
                ENDIF
                IF(SC(2).GT.FCCC) SC(2)=FCCC
                IF(SC(3).GT.FCCC) SC(3)=FCCC
            ENDIF

            IF(LCRACK(2).EQ.1)THEN
                SS=MIN(SC(1),SC(3))
                IF(SS.GT.FCC)THEN
                SC(1)=FCC
                SC(3)=FCC
                ENDIF
                IF(SC(1).GT.FCCC) SC(1)=FCCC
                IF(SC(3).GT.FCCC) SC(3)=FCCC
            ENDIF

            IF(LCRACK(3).EQ.1)THEN
                SS=MIN(SC(1),SC(2))
                IF(SS.GT.FCC)THEN
                SC(1)=FCC
                SC(2)=FCC
                ENDIF
                IF(SC(1).GT.FCCC) SC(1)=FCCC
                IF(SC(2).GT.FCCC) SC(2)=FCCC
            ENDIF
        ENDIF
! C�E�－－－－２軸引張
        IF(NCOM.EQ.0.and.NTEN.EQ.2)THEN
            DO intI=1,3
                IF(LCRACK(intI).EQ.0)THEN
                    SC(intI)=FT
                ENDIF
            ENDDO
        ENDIF
    GOTO 350
    ENDIF
! C
! CCCCCCC/////�E�軸弾性�E�１軸圧壊，１軸ひび割めE////
    IF (NCR0.EQ.1.and.NCR1.EQ.1.and.NCR2.EQ.1)THEN
            DO intI=1,3
            IF(LCRACK(intI).EQ.0)THEN
                IF(SNCR(intI).LE.0.0D0) SC(intI)=FC
                IF(SNCR(intI).GT.0.0D0) SC(intI)=FT
            ENDIF
        ENDDO
        SSS=FC
        FCC=FF65_1*FC
    GOTO 350
    ENDIF
! c         DO intI=1,3
! c             SNN(intI)=SNCR(intI)
! c         ENDDO
! c
! c         DO intI=1,3
! c             IF(LCRACK(intI).EQ.1)THEN
! c                 SNN(intI)=0.0D0
! c             ENDIF
! c         ENDDO
! c         
! c         DO intI=1,3
! c             IF(LCRACK(intI).NE.0)THEN
! c                 SCL(intI)=SC(intI)
! c             ENDIF
! c         ENDDO
! c         
! c CALL XAGRIS(SNN(1),SNN(2),SNN(3),SC(1),SC(2),SC(3),FC,FT,
! c     *                   ALT,ALC,BET,BEC,1) 
! c         DO intI=1,3
! c             IF(LCRACK(intI).NE.0)THEN
! c                 SC(intI)=SCL(intI)
! c             ENDIF
! c         ENDDO
! c
! c     DO intI=1,3
! c         IF(LCRACK(intI).EQ.0)THEN
! c             IF(SNCR(intI).LE.0.0D0) SC(intI)=FC
! c             IF(SNCR(intI).GT.0.0D0) SC(intI)=FT
! c         ENDIF
! c     ENDDO
! c     SSS=FC
! c     FCC=0.65D0*FC
! c GOTO 350
! c ENDIF
! C
! CCCCCCC/////�E�軸圧壊，１軸ひび割めE////
    IF (NCR0.EQ.0.and.NCR1.EQ.1.and.NCR2.EQ.2)THEN
        SSS=FC
        FCC=FF65_1*FC
    GOTO 350
    ENDIF
! C
! CCCCCCC/////�E�軸弾性�E�２軸ひび割めE////
    IF (NCR0.EQ.1.and.NCR1.EQ.2.and.NCR2.EQ.0)THEN
        DO intI=1,3
            IF(LCRACK(intI).EQ.0)THEN
                IF(SNCR(intI).LE.0.0D0) SC(intI)=FC
                IF(SNCR(intI).GT.0.0D0) SC(intI)=FCCC
            ENDIF
        ENDDO
        SSS=FC
        FCC=FF65_1*FC
    GOTO 350
    ENDIF
! C
! CCCCCCC/////�E�軸圧壊，２軸ひび割めE////
    IF (NCR0.EQ.0.and.NCR1.EQ.2.and.NCR2.EQ.1)THEN
        SSS=FC
        FCC=FF65_1*FC
    GOTO 350
    ENDIF
! C
! CCCCCCC/////�E�軸ひび割めE////
    IF (NCR0.EQ.0.and.NCR1.EQ.3.and.NCR2.EQ.0)THEN
        SSS=FC
        FCC=FF65_1*FC
    GOTO 350
    ENDIF
! C
 350    CONTINUE
    SN(1)=SNK1
    SN(2)=SNK2
    SN(3)=SNK3
    SNCR(1)=SNKCR1
    SNCR(2)=SNKCR2
    SNCR(3)=SNKCR3


! C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    DO INTI=1,3
        
    IF(SC(intI).GT.0.0) SC(intI)=FC

    IF(ABS(SC(INTI)).LT.-FC)THEN 
      EC(INTI)=ECU*(-1.6*(SC(INTI)/FC)**3+2.25*(SC(INTI)/FC)**2 &
                         +0.35*SC(INTI)/FC)
    ELSE
        EC(INTI)=ECU*(3.15*SC(INTI)/FC-2.15)
    ENDIF
    ENDDO
    GOTO 360
! C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    IF(ICONCYC.EQ.1)THEN
    DO intI=1,3
    IF(CCRACK(intI).EQ.4.and.LCRACK(intI).EQ.0)THEN
        EU(intI)=EU(intI)-ECR
    ENDIF
    IF(CCRACK(intI).EQ.5.and.LCRACK(intI).EQ.0)THEN
        EU(intI)=EU(intI)-ECCR(intI)
    ENDIF
    ENDDO
    ENDIF


! C *** DECIDE STRAIN FOR MAXIMUM COMPRESSIVE STRESSES ***
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CC  SEANZ MODEL OR FAFITIS-SHAH MODEL の選抁E
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    CALL sPOSNEG  (EU,LCRACK,LPOSNEG,NPOS,NNEG)
       IF(ICIC1.EQ.0) THEN
        ECUUU=((EO*ECU*ECU+2.0D0*FCCC*ECU-FCCC*EO*ECU*ECU/FC)-          &
             SQRT((EO*ECU*ECU+2.0D0*FCCC*ECU-FCCC*EO*ECU*ECU/FC)**2-    &
             4.0D0*FCCC*FCCC*ECU*ECU))*0.5D0/FCCC
 1116 FORMAT(1X,F16.9,1X)
       ELSE
        ECUUU=ECU*(1.0D0-(1.0D0-FCCC/FC)**(1.0D0/(EO*ECU/FC)))
 1114 FORMAT(1X,F16.9,1X)
       ENDIF

       EET=FT/EO

      AET =EET*10.0D0**5
      AECU=ECU*10.0D0**5

! CCCCCCC/////�E�軸弾性/////
    IF (NCR0.EQ.3.and.NCR1.EQ.0.and.NCR2.EQ.0)THEN
        DO intI=1,3
            EUU(intI)=EU(intI)*10.0D0**5
            ECL(intI)=EC(intI)
        ENDDO

        CALL XAGRIS(EUU(1),EUU(2),EUU(3),EC(1),EC(2),EC(3),AECU,AET,    &
                        ALT,ALC,BET,BEC,2) 

        DO intI=1,3
            EC(intI)=EC(intI)*10.0D0**(-5)
        ENDDO
! C�E�－－－－３軸圧縮
        IF(NPOS.EQ.0.and.NNEG.EQ.3)THEN
            DO intI=1,3
                IF(EC(intI).GE.ECUUU) EC(intI)=ECUUU
            ENDDO
            EEE=MIN(EC(1),EC(2),EC(3))
        ENDIF
!C�E�－－－－１軸引張�E�２軸圧縮
        IF(NPOS.EQ.1.and.NNEG.EQ.2)THEN
            IF(LPOSNEG(1).EQ.1)THEN
                EUU(1)=0.0D0
                EUU(2)=EU(2)*10.0D0**5
                EUU(3)=EU(3)*10.0D0**5
        CALL XAGRIS(EUU(1),EUU(2),EUU(3),ECC(1),ECC(2),ECC(3),AECU,AET, &
                        ALT,ALC,BET,BEC,2) 
                ECC(2)=ECC(2)*10.0D0**(-5)
                ECC(3)=ECC(3)*10.0D0**(-5)
                EEE=MIN(ECC(2),ECC(3))
            IF(ICIC1.EQ.0) THEN
                ECUU=((EO*EEE*EEE+2.0D0*FCC*EEE-FCC*EO*EEE*EEE/SSS)-    &
              SQRT((EO*EEE*EEE+2.0D0*FCC*EEE-FCC*EO*EEE*EEE/SSS)**2-    &
              4.0D0*FCC*FCC*EEE*EEE))*0.5D0/FCC
            ELSE
                ECUU=EEE*(1.0D0-(1.0D0-FCC/SSS)**(1.0D0/(EO*EEE/SSS)))
            ENDIF
                D=ECUU/EEE
                EE=MIN(EC(2),EC(3))
                IF(EE.LT.ECUU)THEN
                    EC(1)=EET
                ELSE
                    EC(2)=ECC(2)*D
                    EC(3)=ECC(3)*D
                ENDIF
                IF(EC(2).GT.ECUUU)EC(2)=ECUUU
                IF(EC(3).GT.ECUUU)EC(3)=ECUUU
            ENDIF

            IF(LPOSNEG(2).EQ.1)THEN
                EUU(2)=0.0D0
                EUU(1)=EU(2)*10.0D0**5
                EUU(3)=EU(3)*10.0D0**5
        CALL XAGRIS(EUU(1),EUU(2),EUU(3),ECC(1),ECC(2),ECC(3),AECU,AET, &
                        ALT,ALC,BET,BEC,2) 
                ECC(1)=ECC(1)*10.0D0**(-5)
                ECC(3)=ECC(3)*10.0D0**(-5)
                EEE=MIN(ECC(1),ECC(3))
            IF(ICIC1.EQ.0) THEN
                ECUU=((EO*EEE*EEE+2.0D0*FCC*EEE-FCC*EO*EEE*EEE/SSS)-    &
              SQRT((EO*EEE*EEE+2.0D0*FCC*EEE-FCC*EO*EEE*EEE/SSS)**2-    &
              4.0D0*FCC*FCC*EEE*EEE))*0.5D0/FCC
            ELSE
                ECUU=EEE*(1.0D0-(1.0D0-FCC/SSS)**(1.0D0/(EO*EEE/SSS)))
            ENDIF
                D=ECUU/EEE
                EE=MIN(EC(1),EC(3))
                IF(EE.LT.ECUU)THEN
                    EC(2)=EET
                ELSE
                    EC(1)=ECC(1)*D
                    EC(3)=ECC(3)*D
                ENDIF
                IF(EC(1).GT.ECUUU)EC(1)=ECUUU
                IF(EC(3).GT.ECUUU)EC(3)=ECUUU
            ENDIF

            IF(LPOSNEG(3).EQ.1)THEN
                EUU(3)=0.0D0
                EUU(2)=EU(2)*10.0D0**5
                EUU(1)=EU(3)*10.0D0**5
        CALL XAGRIS(EUU(1),EUU(2),EUU(3),ECC(1),ECC(2),ECC(3),AECU,AET, &
                        ALT,ALC,BET,BEC,2) 
                ECC(2)=ECC(2)*10.0D0**(-5)
                ECC(1)=ECC(1)*10.0D0**(-5)
                EEE=MIN(ECC(2),ECC(1))
            IF(ICIC1.EQ.0) THEN
                ECUU=((EO*EEE*EEE+2.0D0*FCC*EEE-FCC*EO*EEE*EEE/SSS)-    &
              SQRT((EO*EEE*EEE+2.0D0*FCC*EEE-FCC*EO*EEE*EEE/SSS)**2-    &
              4.0D0*FCC*FCC*EEE*EEE))*0.5D0/FCC
            ELSE
                ECUU=EEE*(1.0D0-(1.0D0-FCC/SSS)**(1.0D0/(EO*EEE/SSS)))
            ENDIF
                D=ECUU/EEE
                EE=MIN(EC(1),EC(2))
                IF(EE.LT.ECUU)THEN
                    EC(3)=EET
                ELSE
                    EC(1)=ECC(1)*D
                    EC(2)=ECC(2)*D
                ENDIF
                IF(EC(1).GT.ECUUU)EC(1)=ECUUU
                IF(EC(2).GT.ECUUU)EC(2)=ECUUU
            ENDIF
        ENDIF
! C�E�－－－－２軸引張�E�１軸圧縮
        IF(NPOS.EQ.2.and.NNEG.EQ.1)THEN
            DO intI=1,3
                IF(LPOSNEG(intI).EQ.-1) EE=EC(intI)
            ENDDO

            IF(ICIC1.EQ.0) THEN
                ECUU=((EO*ECU*ECU+2.0D0*FCC*ECU-FCC*EO*ECU*ECU/SSS)-    &
              SQRT((EO*ECU*ECU+2.0D0*FCC*ECU-FCC*EO*ECU*ECU/SSS)**2-    &
              4.0D0*FCC*FCC*ECU*ECU))*0.5D0/FCC
            ELSE
                ECUU=ECU*(1.0D0-(1.0D0-FCC/SSS)**(1.0D0/(EO*ECU/SSS)))
            ENDIF

            IF(LPOSNEG(1).EQ.-1)THEN
                IF(EE.LT.ECUU)THEN
                    EC(2)=EET
                    EC(3)=EET
                ELSE
                    EC(1)=ECUU
                    EC(2)=MAX(EC(2),EC(3))
                    EC(3)=MAX(EC(2),EC(3))
                ENDIF
            ENDIF

            IF(LPOSNEG(2).EQ.-1)THEN
                IF(EE.LT.ECUU)THEN
                    EC(1)=EET
                    EC(3)=EET
                ELSE
                    EC(2)=ECUU
                    EC(1)=MAX(EC(1),EC(3))
                    EC(3)=MAX(EC(1),EC(3))
                ENDIF
            ENDIF

            IF(LPOSNEG(3).EQ.-1)THEN
                IF(EE.LT.ECUU)THEN
                    EC(1)=EET
                    EC(2)=EET
                ELSE
                    EC(3)=ECUU
                    EC(1)=MAX(EC(1),EC(2))
                    EC(2)=MAX(EC(1),EC(2))
                ENDIF
            ENDIF
        ENDIF
! C�E�－－－－３軸引張
        IF(NPOS.EQ.3.and.NNEG.EQ.0)THEN
            DO intI=1,3
                EC(intI)=EET
            ENDDO
        ENDIF
    GOTO 360
    ENDIF

! CCCCCCC/////�E�軸弾性�E�１軸圧壁E////
    IF (NCR0.EQ.2.and.NCR1.EQ.0.and.NCR2.EQ.1)THEN
        DO intI=1,3
            EUU(intI)=EU(intI)*10.0D0**5
        ENDDO
        DO intI=1,3
            IF(LCRACK(intI).EQ.2)THEN
                EUU(intI)=0.0D0
                ECL(intI)=EC(intI)
            ENDIF
        ENDDO
     IF(ICIC1.EQ.0) THEN
        ECUU=((EO*ECU*ECU+2.0D0*FCC*ECU-FCC*EO*ECU*ECU/SSS)-        &
            SQRT((EO*ECU*ECU+2.0D0*FCC*ECU-FCC*EO*ECU*ECU/SSS)**2-  &
            4.0D0*FCC*FCC*ECU*ECU))*0.5D0/FCC
       ELSE
        ECUU=ECU*(1.0D0-(1.0D0-FCC/SSS)**(1.0D0/(EO*ECU/SSS)))
       ENDIF

      CALL XAGRIS(EUU(1),EUU(2),EUU(3),EC(1),EC(2),EC(3),AECU,AET,  &
                    ALT,ALC,BET,BEC,2)

        DO intI=1,3
            IF(LCRACK(intI).EQ.2)THEN
                EC(intI)=ECL(intI)
            ELSE
                EC(intI)=EC(intI)*10.0D0**(-5)
            ENDIF
        ENDDO
! C�E�－－－－１軸圧縮�E�１軸引張
        IF(NPOS.EQ.1.and.NNEG.EQ.1)THEN
            IF(LPOSNEG(2).EQ.1.and.LPOSNEG(3).EQ.-1)THEN
                EE=EC(3)
                IF(EE.LT.ECUU)THEN
                    EC(2)=EET
                ELSE
                    EC(3)=ECUU
                ENDIF
            ENDIF
            IF(LPOSNEG(2).EQ.-1.and.LPOSNEG(3).EQ.1)THEN
                EE=EC(2)
                IF(EE.LT.ECUU)THEN
                    EC(3)=EET
                ELSE
                    EC(2)=ECUU
                ENDIF
            ENDIF
            IF(LPOSNEG(1).EQ.1.and.LPOSNEG(3).EQ.-1)THEN
                EE=EC(3)
                IF(EE.LT.ECUU)THEN
                    EC(1)=EET
                ELSE
                    EC(3)=ECUU
                ENDIF
            ENDIF
            IF(LPOSNEG(1).EQ.-1.and.LPOSNEG(3).EQ.1)THEN
                EE=EC(1)
                IF(EE.LT.ECUU)THEN
                    EC(3)=EET
                ELSE
                    EC(1)=ECUU
                ENDIF
            ENDIF
                IF(LPOSNEG(2).EQ.1.and.LPOSNEG(1).EQ.-1)THEN
                EE=EC(1)
                IF(EE.LT.ECUU)THEN
                    EC(2)=EET
                ELSE
                    EC(1)=ECUU
                ENDIF
            ENDIF
            IF(LPOSNEG(2).EQ.-1.and.LPOSNEG(1).EQ.1)THEN
                EE=EC(2)
                IF(EE.LT.ECUU)THEN
                    EC(1)=EET
                ELSE
                    EC(2)=ECUU
                ENDIF
            ENDIF
        ENDIF
! C�E�－－－－２軸引張
        IF(NPOS.EQ.2.and.NNEG.EQ.0)THEN
            DO intI=1,3
                IF(LCRACK(intI).NE.2)THEN
                    EC(intI)=EET
                ENDIF
            ENDDO
        ENDIF
! C�E�－－－－２軸圧縮
        IF(NPOS.EQ.0.and.NNEG.EQ.2)THEN
            IF(LCRACK(1).EQ.2)THEN
                EE=MIN(EC(2),EC(3))
                IF(EE.GT.ECUU)THEN
                    EC(2)=ECUU
                    EC(3)=ECUU
                ENDIF
                    IF(EC(2).GT.ECUUU) EC(2)=ECUUU
                    IF(EC(3).GT.ECUUU) EC(3)=ECUUU
            ENDIF

            IF(LCRACK(2).EQ.2)THEN
                EE=MIN(EC(1),EC(3))
                IF(EE.GT.ECUU)THEN
                    EC(1)=ECUU
                    EC(3)=ECUU
                ENDIF
                    IF(EC(1).GT.ECUUU) EC(1)=ECUUU
                    IF(EC(3).GT.ECUUU) EC(3)=ECUUU
            ENDIF

                IF(LCRACK(3).EQ.2)THEN
                EE=MIN(EC(1),EC(2))
                IF(EE.GT.ECUU)THEN
                    EC(1)=ECUU
                    EC(2)=ECUU
                ENDIF
                    IF(EC(1).GT.ECUUU) EC(1)=ECUUU
                    IF(EC(2).GT.ECUUU) EC(2)=ECUUU
            ENDIF
        ENDIF
    GOTO 360
    ENDIF
! C
! CCCCCCC/////�E�軸弾性�E�２軸圧壁E////
    IF (NCR0.EQ.1.and.NCR1.EQ.0.and.NCR2.EQ.2)THEN
        DO intI=1,3
            IF(LCRACK(intI).EQ.0)THEN
                IF(EU(intI).LE.0.0D0) EC(intI)=ECU
                IF(EU(intI).GT.0.0D0) EC(intI)=EET
            ENDIF
        ENDDO
    GOTO 360
    ENDIF
! C
! CCCCCCC/////�E�軸圧壁E////
    IF (NCR0.EQ.0.and.NCR1.EQ.0.and.NCR2.EQ.3)GOTO 360
! C
! CCCCCCC/////�E�軸弾性�E�１軸ひび割めE////
    IF (NCR0.EQ.2.and.NCR1.EQ.1.and.NCR2.EQ.0)THEN
        DO intI=1,3
            EUU(intI)=EU(intI)*10.0D0**5
        ENDDO
        DO intI=1,3
            IF(LCRACK(intI).EQ.1)THEN
                EUU(intI)=0.0D0
                ECL(intI)=EC(intI)
            ENDIF
        ENDDO
      CALL XAGRIS(EUU(1),EUU(2),EUU(3),EC(1),EC(2),EC(3),AECU,AET,  &
                   ALT,ALC,BET,BEC,2)
        DO intI=1,3
            IF(LCRACK(intI).EQ.1)THEN
                EC(intI)=ECL(intI)
            ELSE
                EC(intI)=EC(intI)*10.0D0**(-5)
            ENDIF
        ENDDO
     IF(ICIC1.EQ.0) THEN
        ECUU=((EO*ECU*ECU+2.0D0*FCC*ECU-FCC*EO*ECU*ECU/SSS)-        &
            SQRT((EO*ECU*ECU+2.0D0*FCC*ECU-FCC*EO*ECU*ECU/SSS)**2-  &
            4.0D0*FCC*FCC*ECU*ECU))*0.5D0/FCC
       ELSE
        ECUU=ECU*(1.0D0-(1.0D0-FCC/SSS)**(1.0D0/(EO*ECU/SSS)))
       ENDIF
! C�E�－－－－１軸引張�E�１軸圧縮
        IF(NPOS.EQ.1.and.NNEG.EQ.1)THEN
            IF(LPOSNEG(2).EQ.1.and.LPOSNEG(3).EQ.-1)THEN
                EE=EC(3)
                IF(EE.LT.ECUU)THEN
                    EC(2)=EET
                ELSE
                    EC(3)=ECUU
                ENDIF
            ENDIF
            IF(LPOSNEG(2).EQ.-1.and.LPOSNEG(3).EQ.1)THEN
                EE=EC(2)
                IF(EE.LT.ECUU)THEN
                    EC(3)=EET
                ELSE
                    EC(2)=ECUU
                ENDIF
            ENDIF
            IF(LPOSNEG(1).EQ.1.and.LPOSNEG(3).EQ.-1)THEN
                EE=EC(3)
                IF(EE.LT.ECUU)THEN
                    EC(1)=EET
                ELSE
                    EC(3)=ECUU
                ENDIF
            ENDIF
            IF(LPOSNEG(1).EQ.-1.and.LPOSNEG(3).EQ.1)THEN
                EE=EC(1)
                IF(EE.LT.ECUU)THEN
                    EC(3)=EET
                ELSE
                    EC(1)=ECUU
                ENDIF
            ENDIF
            IF(LPOSNEG(2).EQ.1.and.LPOSNEG(1).EQ.-1)THEN
                EE=EC(1)
                IF(EE.LT.ECUU)THEN
                    EC(2)=EET
                ELSE
                    EC(1)=ECUU
                ENDIF
            ENDIF
            IF(LPOSNEG(2).EQ.-1.and.LPOSNEG(1).EQ.1)THEN
                EE=EC(2)
                IF(EE.LT.ECUU)THEN
                    EC(1)=EET
                ELSE
                    EC(2)=ECUU
                ENDIF
            ENDIF
        ENDIF
! C�E�－－－－２軸引張
        IF(NPOS.EQ.2.and.NNEG.EQ.0)THEN
            DO intI=1,3
                IF(LCRACK(intI).NE.1)THEN
                    EC(intI)=EET
                ENDIF
            ENDDO
        ENDIF
! C�E�－－－－２軸圧縮
        IF(NPOS.EQ.0.and.NNEG.EQ.2)THEN
            IF(LCRACK(1).EQ.1)THEN
                EE=MIN(EC(2),EC(3))
                IF(EE.GT.ECUU)THEN
                    EC(2)=ECUU
                    EC(3)=ECUU
                ENDIF
                    IF(EC(2).GT.ECUUU) EC(2)=ECUUU
                    IF(EC(3).GT.ECUUU) EC(3)=ECUUU
            ENDIF

            IF(LCRACK(2).EQ.1)THEN
                EE=MIN(EC(1),EC(3))
                IF(EE.GT.ECUU)THEN
                    EC(1)=ECUU
                    EC(3)=ECUU
                ENDIF
                    IF(EC(1).GT.ECUUU) EC(1)=ECUUU
                    IF(EC(3).GT.ECUUU) EC(3)=ECUUU
            ENDIF

                IF(LCRACK(3).EQ.1)THEN
                EE=MIN(EC(1),EC(2))
                IF(EE.GT.ECUU)THEN
                    EC(1)=ECUU
                    EC(2)=ECUU
                ENDIF
                    IF(EC(1).GT.ECUUU) EC(1)=ECUUU
                    IF(EC(2).GT.ECUUU) EC(2)=ECUUU
            ENDIF
        ENDIF
    GOTO 360
    ENDIF

! CCCCCCC/////�E�軸弾性�E�１軸圧壊，１軸ひび割めE////
    IF (NCR0.EQ.1.and.NCR1.EQ.1.and.NCR2.EQ.1)THEN
            DO intI=1,3
            IF(LCRACK(intI).EQ.0)THEN
                IF(EU(intI).LE.0.0D0) EC(intI)=ECU
                IF(EU(intI).GT.0.0D0) EC(intI)=EET
            ENDIF
        ENDDO
    GOTO 360
    ENDIF
! c     DO intI=1,3
! c         IF(LCRACK(intI).EQ.1)THEN
! c             EUU(intI)=0.0D0
! c         ELSE
! c             EUU(intI)=EU(intI)*10.0D0**5
! c             ECL(intI)=EC(intI)
! c         ENDIF
! c     ENDDO
! c
! c      CALL XAGRIS(EUU(1),EUU(2),EUU(3),EC(1),EC(2),EC(3),AECU,AET,
! c     *               ALT,ALC,BET,BEC,2)
! c     DO intI=1,3
! c         IF(LCRACK(intI).NE.0)THEN
! c             EC(intI)=ECL(intI)
! c         ENDIF
! c     ENDDO
! c
! c     DO intI=1,3
! c         IF(LCRACK(intI).EQ.0)THEN
! c             EC(intI)=EC(intI)*10.0D0**(-5)
! c                 IF(EU(intI).LE.0.0D0) EC(intI)=ECU
! c                 IF(EU(intI).GT.0.0D0) EC(intI)=EET
! c         ENDIF
! c     ENDDO
! c GOTO 360
! c ENDIF

! CCCCCCC/////�E�軸ひび割れ，２軸圧壁E////
    IF (NCR0.EQ.0.and.NCR1.EQ.1.and.NCR2.EQ.2)GOTO 360

! CCCCCCC/////�E�軸弾性�E�２軸ひび割めE////
    IF (NCR0.EQ.1.and.NCR1.EQ.2.and.NCR2.EQ.0)THEN
        DO intI=1,3
            IF(LCRACK(intI).EQ.0)THEN
                IF(EU(intI).LE.0.0D0) EC(intI)=ECU
                IF(EU(intI).GT.0.0D0) EC(intI)=EET
            ENDIF
        ENDDO
    GOTO 360
    ENDIF
! C
! CCCCCCC/////�E�軸圧壊，２軸ひび割めE////
    IF (NCR0.EQ.0.and.NCR1.EQ.2.and.NCR2.EQ.1)GOTO 360
! C
! CCCCCCC/////�E�軸ひび割めE////
    IF (NCR0.EQ.0.and.NCR1.EQ.3.and.NCR2.EQ.0)GOTO 360
! C
360 CONTINUE

    EU(1)=EUK1
    EU(2)=EUK2
    EU(3)=EUK3
! CCCCCCCCCCCCCCCCC
! C IF(ICONCYC.EQ.1)THEN
    DO INTI=1,3
    IF(iul(intI).ne.1)THEN
        SC(intI)=SCDMY(intI)
        EC(intI)=ECDMY(intI)
! c     SC(intI)=FC
! c   EC(intI)=ECU
    ENDIF
    ENDDO
! c     IF(LCRACK(intI).EQ.1)THEN
! c     SC(INTI)= FC
! c     IF(ABS(SC(INTI)).LT.-FC)THEN
! c      EC(INTI)=ECU*(-1.6*(SC(INTI)/FC)**3+2.25*(SC(INTI)/FC)**2
! c     *                    +0.35*SC(INTI)/FC)
! c      ELSE
! c         EC(INTI)=ECU*(3.15*SC(INTI)/FC-2.15)
! c     ENDIF
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.42.or.IUL(intI).EQ.46.or.
! c     *     IUL(intI).EQ.52.or.IUL(intI).EQ.56)THEN
! c         IF(LCRACK(intI).EQ.0.and.SN(intI).GT.0.0D0.and.
! c     *       EU(intI).GT.0.0D0)THEN
! c             SC(INTI)= FC
! c             IF(ABS(SC(INTI)).LT.-FC)THEN
! c              EC(INTI)=ECU*(-1.6*(SC(INTI)/FC)**3+2.25*(SC(INTI)/FC)**2
! c     *                    +0.35*SC(INTI)/FC)
! c              ELSE
! c                 EC(INTI)=ECU*(3.15*SC(INTI)/FC-2.15)
! c             ENDIF
! c         ENDIF
! c     ENDIF
    
! C IF(SC(intI).GT.0.0) THEN
! C     SC(intI)=FC
! C     EC(INTI)=ECU
! C ENDIF
! C IF(EC(intI).GT.0.0) THEN
! C     SC(intI)=FC
! C     EC(INTI)=ECU
! C ENDIF
! C IF(LCRACK(intI).EQ.1) THEN
! C     SC(intI)=FC
! C     EC(INTI)=ECU
! C ENDIF
! C ENDDO
    
! C ENDIF
! CCCCCCCCCCCCCCCCC
! CCCCCCCC     ADDED BY KENGO (2003.10.02)start
! c DO intI=1,3
! c IF(ICR.EQ.0)THEN
! c
! c     IF(IUL(intI).EQ.41)THEN
! c         EC(intI)=ETMAX(intI)
! c         SC(intI)=STMAX(intI)
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.42.and.LCRACK(intI).EQ.1.
! c     *and.SN(intI).LE.0.0D0)THEN
! c             EC(intI)=ECU
! c             SC(intI)=FC
! c             EC(intI)=EC(intI)-ECR
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.42.and.LCRACK(intI).EQ.1.
! c     *and.SN(intI).GT.0.0D0)THEN
! c             EC(intI)=EPT(intI)
! c             SC(intI)=0.0D0
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.46.and.LCRACK(intI).EQ.1.
! c     *and.SN(intI).LE.0.0D0)THEN
! c             EC(intI)=EPT(intI)
! c             SC(intI)=0.0D0
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.46.and.LCRACK(intI).EQ.1.
! c     *and.SN(intI).GT.0.0D0)THEN
! c     EC(intI)=ETMAX(intI)
! c     SC(intI)=STMAX(intI)
! C     EC(intI)=ECU
! C SC(intI)=FC
! C     EC(intI)=EC(intI)-ECR
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.51)THEN
! c         EC(intI)=ETMAX(intI)
! c         SC(intI)=STMAX(intI)
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.52.and.LCRACK(intI).EQ.1.
! c     *and.SN(intI).LE.0.0D0)THEN
! c             EC(intI)=ECU
! c             SC(intI)=FC
! c             EC(intI)=EC(intI)-ECCR(intI)
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.52.and.LCRACK(intI).EQ.1.
! c     *and.SN(intI).GT.0.0D0)THEN
! c             EC(intI)=EPT(intI)
! c             SC(intI)=0.0D0
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.56.and.LCRACK(intI).EQ.1.
! c     *and.SN(intI).LE.0.0D0)THEN
! c             EC(intI)=EPT(intI)
! c             SC(intI)=0.0D0
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.56.and.LCRACK(intI).EQ.1.
! c     *and.SN(intI).GT.0.0D0)THEN
! c     EC(intI)=ETMAX(intI)
! c     SC(intI)=STMAX(intI)
! C     EC(intI)=ECU
! C     SC(intI)=FC
! C     EC(intI)=EC(intI)-ECCR(intI)
! c     ENDIF
! c
! c ELSE
! c
! c     IF(IUL(intI).EQ.41)THEN
! c         EC(intI)=ETMAX(intI)
! c         SC(intI)=STMAX(intI)
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.42.and.LCRACK(intI).EQ.1.
! c     *and.SNCR(intI).LE.0.0D0)THEN
! c             EC(intI)=ECU
! c             SC(intI)=FC
! c             EC(intI)=EC(intI)-ECR
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.42.and.LCRACK(intI).EQ.1.
! c     *and.SNCR(intI).GT.0.0D0)THEN
! c             EC(intI)=EPT(intI)
! c             SC(intI)=0.0D0
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.46.and.LCRACK(intI).EQ.1.
! c     *and.SNCR(intI).LE.0.0D0)THEN
! c             EC(intI)=EPT(intI)
! c             SC(intI)=0.0D0
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.46.and.LCRACK(intI).EQ.1.
! c     *and.SNCR(intI).GT.0.0D0)THEN
! c     EC(intI)=ETMAX(intI)
! c     SC(intI)=STMAX(intI)
! C     EC(intI)=ECU
! C     SC(intI)=FC
! C     EC(intI)=EC(intI)-ECR
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.51)THEN
! c         EC(intI)=ETMAX(intI)
! c         SC(intI)=STMAX(intI)
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.52.and.LCRACK(intI).EQ.1.
! c     *and.SNCR(intI).LE.0.0D0)THEN
! c             EC(intI)=ECU
! c             SC(intI)=FC
! c             EC(intI)=EC(intI)-ECCR(intI)
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.52.and.LCRACK(intI).EQ.1.
! c     *and.SNCR(intI).GT.0.0D0)THEN
! c             EC(intI)=EPT(intI)
! c             SC(intI)=0.0D0
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.56.and.LCRACK(intI).EQ.1.
! c     *and.SNCR(intI).LE.0.0D0)THEN
! c             EC(intI)=EPT(intI)
! c             SC(intI)=0.0D0
! c     ENDIF
! c
! c     IF(IUL(intI).EQ.56.and.LCRACK(intI).EQ.1.
! c     *and.SNCR(intI).GT.0.0D0)THEN       EC(intI)=ETMAX(intI)
! c     SC(intI)=STMAX(intI)
! C     EC(intI)=ECU
! C     SC(intI)=FC
! C     EC(intI)=EC(intI)-ECCR(intI)
! c     ENDIF
! c ENDIF
! c ENDDO

! CCCCCCCC     ADDED BY KENGO (2003.10.02)end
! C Do intI=1,3
! C IF(IUL(intI).EQ.42)THEN
! C CALL CHECK_SC_EC(IUL(1),IUL(2),IUL(3),SN(1),SN(2),SN(3),
! C     *                   SC(1),SC(2),SC(3),EC(1),EC(2),EC(3),NST,LCOMTEN,LPOSNEG)
! C ENDIF
! C ENDDO
    END

        !************************************************************************************** 
            SUBROUTINE sSCEC_hj(SND,SNCRD,EU,IUL,IROT,ICR,FC,ECU,EO,FT,SC,EC,       &
        EEN,EBU,ICIC1,ECCR,NST,CCRACK,EPCU,ICONCYC,SCHJ,ICIC8)
    
    real (kreal):: SC(3),EC(3),SCC(3),SN(3),SNN(3),SCL(3),SNCR(3),EPT(3),   &
            EUU(3),EU(3),ECL(3),ECC(3),EEN(3),                              &
            ECCR(3),CCRACK(3),STMAX(3),ETMAX(3)
    real(kreal) :: SND(3),SNCRD(3),SCHJ(3)
    integer :: LCRACK(3),IUL(3),LPOSNEG(3),LCOMTEN(3)
    real(kreal) :: SCDMY(3),ECDMY(3),eudm(6)
    
    integer :: irot, icr, icic1,nst,iconcyc,icic8
    real(kreal) :: fc,ecu,eo,ft,ebu,epcu
    
    real (kreal) :: ALT,ALC,BET,BEC,FF30,FF65_1,XBAI
    integer :: I
    !C       KUPHER
     ALT=ABS(FT/FC)
     ALC=1.15D0
     BET=ABS(FT/EO/ECU)
     BEC=1.7512D0

    FF30=1.0
    FF65_1=0.85

    CALL XAGRIS(SND(1),SND(2),SND(3),SC(1),SC(2),SC(3),FC,FT,       &
                        ALT,ALC,BET,BEC,1)

    XBAI=4.
    DO I=1,3
        IF(SC(I)<XBAI*FC) THEN
            SC(I)=XBAI*FC
        ENDIF
    ENDDO
    DO I=1,3
        IF(     ICIC8==2) THEN
            EC(I)=ECU*SC(I)/FC
        ELSE IF(ICIC8==3) THEN
            EC(I)=ECU
        ENDIF
    ENDDO

    END

        !********************************************************************
            SUBROUTINE sComTen (SN,LCRACK,LCOMTEN,NCOM,NTEN)
    ! C ***********************************************************************
    ! *引数名：sComTen            added by Kengo(2003.4.16)
    ! *機�E�E�圧縮か引張か判別
    ! *入力：EU�E�E
    ! *　　　SN�E�E
    ! *　　　IUL�E�E
    ! *出力：LCOMTEN�E�引張のとぁE�E�圧縮のとぁE1
    ! *目皁E��E
    ! ***********************************************************************

    real (kreal) SN(3)
    integer NCOM,NTEN,intI,LCOMTEN(3),LCRACK(3)
    
    NCOM=0
    NTEN=0

    DO intI=1,3
        IF(LCRACK(intI).EQ.0)THEN
            IF(SN(intI).GT.0)THEN
                    LCOMTEN(intI) = 1
                    NTEN=NTEN+1
            ELSE
                    LCOMTEN(intI) = -1
                    NCOM=NCOM+1
            ENDIF
        ENDIF
    ENDDO
    END
    
        !********************************************************************
            SUBROUTINE XAGRIS(SNN1,SNN2,SNN3,SC1,SC2,SC3,FC,FT, &
                        ALT,ALC,BET,BEC,JJ)


! C
! C *** DECIDE MAXIMUM COMPRESSIVE STRESSE  S ***
! C                                       WITH ARGYRIS FAILURE SURFACE
! C                               <1991.8.27> <K.UCHIDA AND A.AMEMIYA>
! C
! c  ***  Original Kupher's concrete data ***
! C      CONCRETE(FC=-326.89,ECU=-0.00215,FT=29.74699,ET=0.0000991)
! C
! C  *** DATA For unknown tension strength or max strain ***
! c          DATA ALT,ALC/0.091,1.15/
! c          DATA BET,BEC/0.0461,1.7512/
! C
! C  *** using input data of tension strength or max strain ***
! C      ALT=ABS(FT/FC)
! C      ALC=1.15
! C      BET=ABS(FT/330457./.00215)
! C      BEC=1.7512
! C*******************************************************************
! C
! C       write(6,*) 'ALT=',ALT,'   ALC=',ALC,'  BET=',BET,'  BEC=',BEC
! C       WRITE(6,*) 'FC=',FC,'   FT=',FT
! C
! C     <JJ=1>  :STRESS
! C     <JJ=2>  :STRAIN
! C
! C     IF(ABS(SNN1)<1.E-25.AND.ABS(SNN2)<1.E-25.AND.
! C     *                      ABS(SNN3)<1.E-25)THEN !
! CH--- *** ---
    real(kreal) ::  SNN1,SNN2,SNN3,SC1,SC2,SC3,FC,FT,ALT,ALC,BET,BEC
    integer     ::  JJ

    real(kreal) :: A,AAA,BBB,CCC,SMAX,SSNN,SSNN1,SSNN2,SSNN3
    real(kreal) :: ZK1,ZK2,ZK3,PSIM,PSIMF,TM,CC,TMF
    real(kreal) :: GZA,RO1,RO2,A0,A1,A2,a1a1a0,GZA0,B0,B1,B2,R1,R2,R21,R    
    real(kreal) :: S,DD,EE,FF,DSN1,DSN2,DSN3
    integer     :: III,JJJ,KKK,ICA,I,K
    


    A=ABS(FT/20.)
    IF(ABS(SNN1)<A.AND.ABS(SNN2)<A.AND.ABS(SNN3)<A)THEN !
        DO I=1,3
            SC1=FC ; SC2=FC ; SC3=FC
        END DO
        GOTO 999
    ENDIF
    IF(SNN1==SNN2.AND.SNN2==SNN3.AND.SNN1==SNN3) SNN1 = SNN1 * 0.999

      III=0
      JJJ=0
      KKK=0

      IF(SNN2.GT.SNN1) THEN
          AAA=SNN2
          SNN2=SNN1
          SNN1=AAA
          III=1
      ENDIF
!C       SNN1>SNN2
      IF(SNN3.GT.SNN1) THEN
          BBB=SNN3
          SNN3=SNN1
          SNN1=BBB
          JJJ=1
      ENDIF
!C       SNN1>SNN3
      IF(SNN3.GT.SNN2) THEN
          CCC=SNN3
          SNN3=SNN2
          SNN2=CCC
          KKK=1
      ENDIF
! C       SNN2>SNN3
! C *** HIRITU
        SMAX=MAX(ABS(SNN1),ABS(SNN2),ABS(SNN3))
!CH---修正2--------------------------------------------------------
      SSNN1=ABS(SNN1);SSNN2=ABS(SNN2);SSNN3=ABS(SNN3)
    SSNN=SSNN1+SSNN2+SSNN3
!   IF(SSNN == 0.0)write(1861,*)'●●●修正2: SSNN=',SSNN
! C      IF(SSNN == 0.0) SMAX=1.0
! CH---修正2--------------------------------------------------------
! C
        ZK1=SNN1/SMAX
        ZK2=SNN2/SMAX
        ZK3=SNN3/SMAX

! C *** INITIAL VALUE

        SNN1=ZK1*10.D0
        SNN2=ZK2*10.D0
        SNN3=ZK3*10.D0

        ICA=0

      DO 8000 K=1,5

        DSN1=ZK1*10.D0**(1-K)
        DSN2=ZK2*10.D0**(1-K)
        DSN3=ZK3*10.D0**(1-K)

      DO 1000 I=1,100000
        SNN1=SNN1+DSN1
        SNN2=SNN2+DSN2
        SNN3=SNN3+DSN3

        ICA=ICA+1

  100   PSIM=(SNN1+SNN2+SNN3)/3.D0
        PSIMF=(-1.D0)*PSIM/FC

    !A = ISNAN(SNN1)+ISNAN(SNN2)+ISNAN(SNN2)
!   IF(A < 0.) WRITE(1861,*)'☁E��正2 SNN1,SNN2,SNN3',SNN1,SNN2,SNN3
! C IF(A < 0.) THEN
! C     SNN1=1.;SNN2=1.;SNN3=1.
! C END IF
! CH-------------
        TM=SQRT((SNN1-SNN2)**2+(SNN2-SNN3)**2+(SNN3-SNN1)**2)   &
         /SQRT(15.D0)
        TMF=(-1.D0)*TM/FC

        CC=(2.D0*SNN1-SNN2-SNN3)/                               &
            SQRT(2.D0*((SNN1-SNN2)**2+(SNN2-SNN3)**2+(SNN3-SNN1)**2))
! CH---修正3--------------------------------------------------------
      SSNN1=ABS(SNN1-SNN2);SSNN2=ABS(SNN2-SNN3);SSNN3=ABS(SNN3-SNN1)
    SSNN=SSNN1+SSNN2+SSNN3
! C A = ISNAN(SNN1) ; IF(A<0.) SNN1=0.0
! C A = ISNAN(SNN2) ; IF(A<0.) SNN2=0.0
! C A = ISNAN(SNN3) ; IF(A<0.) SNN3=0.0
! C A = ISNAN(CC) ; IF(A<0.) CC=0.0
      IF(SSNN == 0.0) THEN
!   write(1861,*)'●●●修正3: SSNN=',SSNN,CC
    ENDIF
!C      IF(SSNN == 0.0) CC=0.0
!CH---修正3--------------------------------------------------------
        IF(CC.GT.1.0D0) CC=1.0D0
        IF(CC.LT.0.5D0) CC=0.5D0
! C
! C*****************************************************************
! C       MODULUS FOR STRESS FOR KUPHER DATA
! C        IF(JJ.EQ.2) GOTO 157
! C        A2=-0.038292
! C        A1=-0.513326
! C        A0=0.048835
! C
! C        B2=-0.232785
! C        B1=-0.909151
! C        B0=0.087963
! C
! C        GOTO 155
! C*****************************************************************
! C       MODULUS FOR STRAIN FOR KUPHER DATA
! C  157 CONTINUE
! C        A2=-0.010786
! C        A1=-0.538802
! C        A0=0.025115
! C
! C        B2=-0.230689
! C        B1=-1.027311
! C        B0=0.048340
! C
! C        GOTO 155
! C*****************************************************************
! C  155 GZA0=((-1)*A1-SQRT(A1*A1-4.*A0*A2))/2./A2
! C
! C      GZA=((SQRT(1.2)+3.*B1*GZA0+3.*B2*GZA0-B2-B1)-
! C     *     SQRT((B2-B1-SQRT(1.2)-3.*B1*GZA0-3.*B2*GZA0)**2
! C     *     -4.*(3.*B2*GZA0+4.*B2)*(B1*GZA0+SQRT(1.2)/3.
! C     *     +B1/3.-B2/9.+2/3*B2*GZA0)))/2./(3.*B2*GZA0+4.*B2)
! C
! C      GZAP=((SQRT(1.2)+3.*B1*GZA0+3.*B2*GZA0-B2-B1)+
! C     *     SQRT((B2-B1-SQRT(1.2)-3.*B1*GZA0-3.*B2*GZA0)**2
! C     *     -4.*(3.*B2*GZA0+B2+4.*B2)*(B1*GZA0+SQRT(1.2)/3.
! C     *     +B1/3.-B2/9.+2/3*B2*GZA0)))/2./(3.*B2*GZA0+4.*B2)
! C
! C      GZAM=((SQRT(1.2)+3.*B1*GZA0+3.*B2*GZA0-B2-B1)-
! C     *     SQRT((B2-B1-SQRT(1.2)-3.*B1*GZA0-3.*B2*GZA0)**2
! C     *     -4.*(3.*B2*GZA0+B2+4.*B2)*(B1*GZA0+SQRT(1.2)/3.
! C     *     +B1/3.-B2/9.+2/3*B2*GZA0)))/2./(3.*B2*GZA0+4.*B2)
! C
! C      IF(JJ.EQ.1) THEN
! C
! C      RO1=(A2*(2.*ALC+ALT)*(3.*(-GZA)-2.*ALC)*(3.*(-GZA)+ALT)/9.-
! C     *    SQRT(1.2)*(-GZA)*(ALT-ALC)+SQRT(1.2)*ALT*ALC)/
! C     *    (2.*ALC+ALT)
! C      ELSE
! C      RO1=(A2*(2.*BEC+BET)*(3.*(-GZA)-2.*BEC)*(3.*(-GZA)+BET)/9.-
! C     *    SQRT(1.2)*(-GZA)*(BET-BEC)+SQRT(1.2)*BET*BEC)/
! C     *    (2.*BEC+BET)
! C      ENDIF
! C
! C     RO2=B2*((-GZA)+GZA0)*(3.*(-GZA)-1.)*(3.*GZA0+1.)/(9.*GZA0+3.)+
! C     *    SQRT(2./15.)*(GZA0+(-GZA))/(GZA0+1/3.)
! C*******************************************************************
! C STRESS FAILURE SURFACE  GZA=  -8.935134E-01
! C STRESS FAILURE SURFACE  RO1=   4.769774E-01
! C STRESS FAILURE SURFACE  RO2=   7.144530E-01
! C STRAIN FAILURE SURFACE  GZA=      -1.071775
! C STRAIN FAILURE SURFACE  RO1=   5.902000E-01
! C STRAIN FAILURE SURFACE  RO2=   8.844008E-01
! C*******************************************************************
! C
      IF(JJ.EQ.1) THEN
! C
! C **< kupher >*************************
! C
        GZA=-8.935134D-01
        RO1= 4.769774D-01
        RO2= 7.144530D-01
! C
! C **< richart >************************
! C
! c      GZA=-2.5D0
! c      RO1= 1.1D0
! c      RO2= 1.5D0
! C
! C **< schickert >**********************
! C
! C       GZA=-8.232883D-01
! C       RO1= 4.493055D-01
! C       RO2= 6.632722D-01
! C
      A2=9.D0*(SQRT(1.2D0)*(-1.D0)*GZA*(ALT-ALC)        &
             -SQRT(1.2D0)*ALT*ALC+RO1*(2.D0*ALC+ALT))   &
             /(2.D0*ALC+ALT)                            &
             /(3.D0*(-1.D0)*GZA-2.D0*ALC)               &
             /(3.D0*(-1.D0)*GZA+ALT)                    
      A1=(2.D0*ALC-ALT)*A2/3.D0+SQRT(1.2D0)*(ALT-ALC)/(2.D0*ALC+ALT)
      A0=2.D0*ALC*A1/3.D0-4.D0*ALC*ALC*A2/9.D0+SQRT(2.D0/15.D0)*ALC

      ELSE
! C
! C **< kupher >*************************
! C
        GZA=-1.071775D0
        RO1= 5.902000D-01
        RO2= 8.844008D-01
! c
! C **< richart >************************
! C
! c      GZA=-6.8925D0
! c      RO1= 4.5D0
! c      RO2= 7.0D0
! C
! C:1     GZA=-6.8925D0
! C       RO1= 5.0D0
! C       RO2= 7.5D0
! C
! C:2     GZA=-6.8925D0
! C       RO1= 6.0D0
! C       RO2= 9.0D0
! C
! C:3     GZA=-6.8925D0
! C       RO1= 4.0D0
! C       RO2= 6.0D0
! C
! C:4     GZA=-6.8925D0
! C       RO1= 3.5D0
! C       RO2= 5.5D0
! C
! C:5     GZA=-6.8925D0
! C       RO1= 4.0D0
! C       RO2= 6.0D0
! C
! C:6     GZA=-3.0D0
! C       RO1= 1.364D0
! C       RO2= 1.86D0
! C
! C:7     GZA=-6.8925D0
! C       RO1= 4.531D0
! C       RO2= 6.1793D0
! C
! C **< schickert >**********************
! C
! C       GZA=-5.420829D0
! C       RO1= 2.24377D0
! C       RO2= 4.050877D0
! C
      A2=9.D0*(SQRT(1.2D0)*(-1.D0)*GZA*(BET-BEC)        &
             -SQRT(1.2D0)*BET*BEC+RO1*(2.D0*BEC+BET))   &
             /(2.D0*BEC+BET)                            &
             /(3.D0*(-1.D0)*GZA-2.D0*BEC)               &
             /(3.D0*(-1.D0)*GZA+BET)                    
      A1=(2.D0*BEC-BET)*A2/3.D0+SQRT(1.2D0)*(BET-BEC)/(2.D0*BEC+BET)
      A0=2.D0*BEC*A1/3.D0-4.D0*BEC*BEC*A2/9.D0+SQRT(2.D0/15.D0)*BEC

      ENDIF

!C      GZA0=((-1.D0)*A1-SQRT(A1*A1-4.D0*A0*A2))/2.D0/A2
    a1a1a0=A1*A1-4.D0*A0*A2
    if(a1a1a0<0.)then
        write(1861,*)'☁E��正4,SQRT=',a1a1a0
        a1a1a0=0.
    endif
       GZA0=((-1.D0)*A1-SQRT(a1a1a0))/2.D0/A2
!CH-----------------------------------



      B2=9.*(RO2*(GZA0+1.D0/3.D0)-SQRT(2.D0/15.D0)*(GZA0+(-1)*GZA)) &
          /((-1.D0)*GZA+GZA0)                                       &   
          /(3.D0*(-1.D0)*GZA-1.D0)                                  &
          /(3.D0*GZA0+1.D0)                                     
      B1=((-1.D0)*GZA+1.D0/3.D0)*B2                                 &
          +(SQRT(1.2D0)-3.D0*RO2)/(3.D0*(-1.D0)*GZA-1.D0)
      B0=-GZA0*B1-GZA0*GZA0*B2

  150   R1=A2*(PSIMF)**2+A1*(PSIMF)+A0
        R2=B2*(PSIMF)**2+B1*(PSIMF)+B0

        R21=R2**2-R1**2

!CH--------------------------------
    !A = ISNAN(CC)
    IF(A < 0.) WRITE(1861,*)'☁E��正3 CC=',CC

    a1a1a0=4.D0*R21*CC**2+5.D0*R1**2-4.D0*R1*R2
    if(a1a1a0<0.)then
        write(1861,*)'☁E��正5,SQRT=',a1a1a0
        a1a1a0=0.
    endif
       GZA0=((-1.D0)*A1-SQRT(a1a1a0))/2.D0/A
! C        R=(2.D0*R2*R21*CC
! C     *    +R2*(2.D0*R1-R2)*SQRT(4.D0*R21*CC**2+5.D0*R1**2-4.D0*R1*R2))
! C     *    /(4.D0*R21*CC**2+(R2-2.D0*R1)**2)\
        R=(2.D0*R2*R21*CC                       &
         +R2*(2.D0*R1-R2)*SQRT(a1a1a0))         &
         /(4.D0*R21*CC**2+(R2-2.D0*R1)**2)
!CH-----------------------------------------

        S=R-TMF

        IF(S.LT.0.D0) THEN
        SNN1=SNN1-DSN1
        SNN2=SNN2-DSN2
        SNN3=SNN3-DSN3
        IF(K.EQ.5) GOTO 2000
        GO TO 8000
        ENDIF
 1000   CONTINUE

 8000 CONTINUE

 2000   SC1=SNN1
        SC2=SNN2
        SC3=SNN3

        IF(KKK.EQ.1) THEN
        DD=SC2
        SC2=SC3
        SC3=DD
        ENDIF
        IF(JJJ.EQ.1) THEN
        EE=SC1
        SC1=SC3
        SC3=EE
        ENDIF
       IF(III.EQ.1) THEN
        FF=SC1
        SC1=SC2
        SC2=FF
        ENDIF

! C       WRITE(6,*) 'ICA=',ICA,' 1 ',SC1,' 2 ',SC2,' 3 ',SC3
999     RETURN

      END

        !********************************************************************
            SUBROUTINE sPosNeg (EU,LCRACK,LPOSNEG,NPOS,NNEG)
    
        ! ***********************************************************************
! *引数名：sPosNeg            added by Kengo(2003.5.12)
! *機�E�E��Eずみの正負を判宁E
! *入力：EU�E�E
! *出力：LPOSNEG�E�＋�EとぁE(引張側)�E�－�EとぁE1(圧縮側)
! *目皁E��E
! ***********************************************************************

    real(kreal) :: EU(3)
    integer ::  LCRACK(3),LPOSNEG(3)
    integer :: NPOS,NNEG,intI

    NPOS=0
    NNEG=0

    DO intI=1,3
        IF(LCRACK(intI).EQ.0)THEN
            IF(EU(intI).GT.0)THEN
                    LPOSNEG(intI) = 1
                    NPOS=NPOS+1
            ELSE
                    LPOSNEG(intI) = -1
                    NNEG=NNEG+1
            ENDIF
        ENDIF
    ENDDO
    END

        ! ********************************************************************
            SUBROUTINE sRAMDA (EU,BETA,SC,ECU,FC,IUL,           &
                EURR9,EEN,EBU,SN,ECR,ECCR,NST,CCRACK)

        ! **************************************************************
        ! *引用名：sRamda                               Added By KENGO(2003.6.3)
        ! *機�E�E�圧縮強度低減係数λの計箁E
        ! *入力！E
        ! *出力！E
        ! *参老E��E
        ! *********************************************************************


! C     COMPRESSIVE REDUCTION FACTOR (NOGUCHI EQUATION)


    real (kreal) ::  SC(3),EU(3),BETA(3),EUR(3),EURR9(3),ECCR(3),   &
                CCRACK(3),EEN(3),SN(6)
    
    integer :: LCRACK(3),IUL(3),LCRA(3),NST
    
    real (kreal) :: ECU,FC,EBU,ECR,BETABETA,EUREUR
    
    integer :: intI,NTCRACK,NTELS,NUMCRACK,NTFAIL
    
    
    DO intI=1,3
        LCRACK(intI)=0
        EUR(intI)=EU(intI)
        EURR9(intI)=EU(intI)
    ENDDO

    CALL sCrackChk(CCRACK,EU,ECR,ECCR,ECU,EEN,      &
        LCRACK,NTCRACK,NTFAIL,NTELS,NUMCRACK)
    DO intI=1,3
        IF(LCRACK(intI).EQ.2) LCRACK(intI)=0
    ENDDO
    DO intI=1,3
        IF(CCRACK(intI).EQ.5)THEN
            EUR(intI)=EUR(intI)-ECCR(intI)+ECR
        ENDIF
    ENDDO
! c DO intI=1,3
! c     IF(IUL(intI).EQ.4)THEN
! c         LCRACK(intI)=1
! c     ENDIF
! c     IF(IUL(intI).EQ.5)THEN
! c         EUR(intI)=EUR(intI)-ECCR(intI)+ECR
! c         LCRACK(intI)=1
! c     ENDIF
! c     IF(IUL(intI).EQ.41.and.EUR(intI).GT.0.0D0)THEN
! c         LCRACK(intI)=1
! c     ENDIF
! c     IF(IUL(intI).EQ.42.and.EUR(intI).GT.ECR)THEN
! c         LCRACK(intI)=1
! c     ENDIF
! c     IF(IUL(intI).EQ.46.and.EUR(intI).GT.ECR)THEN
! c         LCRACK(intI)=1
! c     ENDIF
! c     IF(IUL(intI).EQ.51.and.EUR(intI).GT.0.0D0)THEN
! c         LCRACK(intI)=1
! c     ENDIF
! c     IF(IUL(intI).EQ.52.and.EUR(intI).GT.ECCR(intI))THEN
! c         EUR(intI)=EUR(intI)-ECCR(intI)+ECR
! c         LCRACK(intI)=1
! c     ENDIF
! c     IF(IUL(intI).EQ.56.and.EUR(intI).GT.ECCR(intI))THEN
! c         EUR(intI)=EUR(intI)-ECCR(intI)+ECR
! c         LCRACK(intI)=1
! c     ENDIF
! c ENDDO
! CCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(LCRACK(1).EQ.0.and.LCRACK(2).EQ.0.and.LCRACK(3).EQ.0)THEN
        GO TO 500
    ENDIF
 
    IF(LCRACK(1).EQ.1.and.LCRACK(2).EQ.0.and.LCRACK(3).EQ.0)THEN
        IF(EUR(1).LT.0.0D0) GO TO 500
        IF(EUR(2).LT.0.0D0) THEN
          BETA(2)=1.0D0/EXP(-0.2*(-EUR(1)/ECU)**0.5*(-FC/250))
          IF(BETA(2).LT.1.0D0) BETA(2)=1.0D0
          SC(2)=SC(2)/BETA(2)

    BETABETA=BETA(2)
    EUREUR=-EUR(1)/ECU

        END IF
        IF(EUR(3).LT.0.0D0) THEN
          BETA(3)=1.0D0/EXP(-0.2*(-EUR(1)/ECU)**0.5*(-FC/250))
          IF(BETA(3).LT.1.0D0) BETA(3)=1.0D0
          SC(3)=SC(3)/BETA(3)
    
    BETABETA=BETA(3)
    EUREUR=-EUR(1)/ECU

        END IF
        GO TO 500
      END IF
    IF(LCRACK(1).EQ.0.and.LCRACK(2).EQ.1.and.LCRACK(3).EQ.0)THEN
        IF(EUR(2).LT.0.0D0) GO TO 500
        IF(EUR(1).LT.0.0D0) THEN
          BETA(1)=1.0D0/EXP(-0.2*(-EUR(2)/ECU)**0.5*(-FC/250))
          IF(BETA(1).LT.1.0D0) BETA(1)=1.0D0
          SC(1)=SC(1)/BETA(1)

    BETABETA=BETA(1)
    EUREUR=-EUR(2)/ECU

        END IF
        IF(EUR(3).LT.0.0D0) THEN
          BETA(3)=1.0D0/EXP(-0.2*(-EUR(2)/ECU)**0.5*(-FC/250))
          IF(BETA(3).LT.1.0D0) BETA(3)=1.0D0
          SC(3)=SC(3)/BETA(3)

    BETABETA=BETA(3)
    EUREUR=-EUR(2)/ECU

        END IF
        GO TO 500
      END IF

    IF(LCRACK(1).EQ.0.and.LCRACK(2).EQ.0.and.LCRACK(3).EQ.1)THEN
        IF(EUR(3).LT.0.0D0) GO TO 500
        IF(EUR(1).LT.0.0D0) THEN
          BETA(1)=1.0D0/EXP(-0.2*(-EUR(3)/ECU)**0.5*(-FC/250))
          IF(BETA(1).LT.1.0D0) BETA(1)=1.0D0
          SC(1)=SC(1)/BETA(1)

    BETABETA=BETA(1)
    EUREUR=-EUR(3)/ECU
 
        END IF
        IF(EUR(2).LT.0.0D0) THEN
          BETA(2)=1.0D0/EXP(-0.2*(-EUR(3)/ECU)**0.5*(-FC/250))
          IF(BETA(2).LT.1.0D0) BETA(2)=1.0D0
          SC(2)=SC(2)/BETA(2)    

    BETABETA=BETA(2)
    EUREUR=-EUR(3)/ECU
 
        END IF
        GO TO 500
      END IF

    IF(LCRACK(1).EQ.1.and.LCRACK(2).EQ.1.and.LCRACK(3).EQ.0)THEN
        IF(EUR(1).LT.0.0D0.AND.EUR(2).LT.0.0D0) GO TO 500
        IF(EUR(3).LT.0.0D0) THEN
          IF(EUR(1).GT.EUR(2)) THEN
                BETA(3)=1.0D0/EXP(-0.2*(-EUR(1)/ECU)**0.5*(-FC/250))

    BETABETA=BETA(3)
    EUREUR=-EUR(1)/ECU
 
          ELSE
                BETA(3)=1.0D0/EXP(-0.2*(-EUR(2)/ECU)**0.5*(-FC/250))

    BETABETA=BETA(3)
    EUREUR=-EUR(2)/ECU
 
          END IF
          IF(BETA(3).LT.1.0D0) BETA(3)=1.0D0
          SC(3)=SC(3)/BETA(3)
        END IF
        GO TO 500
      END IF

    IF(LCRACK(1).EQ.0.and.LCRACK(2).EQ.1.and.LCRACK(3).EQ.1)THEN
        IF(EUR(2).LT.0.0D0.AND.EUR(3).LT.0.0D0) GO TO 500
        IF(EUR(1).LT.0.0D0) THEN
          IF(EUR(2).GT.EUR(3)) THEN
         BETA(1)=1.0D0/EXP(-0.2*(-EUR(2)/ECU)**0.5*(-FC/250))

    BETABETA=BETA(1)
    EUREUR=-EUR(2)/ECU

          ELSE
         BETA(1)=1.0D0/EXP(-0.2*(-EUR(3)/ECU)**0.5*(-FC/250))

    BETABETA=BETA(1)
    EUREUR=-EUR(3)/ECU
 
          END IF
          IF(BETA(1).LT.1.0D0) BETA(1)=1.0D0
          SC(1)=SC(1)/BETA(1)
        END IF
        GO TO 500
      END IF

    IF(LCRACK(1).EQ.1.and.LCRACK(2).EQ.0.and.LCRACK(3).EQ.1)THEN
      IF(EUR(1).LT.0.0D0.AND.EUR(3).LT.0.0D0) GO TO 500
        IF(EUR(2).LT.0.0D0) THEN
          IF(EUR(1).GT.EUR(3)) THEN
          BETA(2)=1.0D0/EXP(-0.2*(-EUR(1)/ECU)**0.5*(-FC/250))

    BETABETA=BETA(2)
    EUREUR=-EUR(1)/ECU

         ELSE
         BETA(2)=1.0D0/EXP(-0.2*(-EUR(3)/ECU)**0.5*(-FC/250))

    BETABETA=BETA(2)
    EUREUR=-EUR(3)/ECU
 
         END IF
         IF(BETA(2).LT.1.0D0) BETA(2)=1.0D0
         SC(2)=SC(2)/BETA(2)
       END IF
       GO TO 500
     END IF

500 CONTINUE
    !CALL CHECK_RAMDA(BETABETA,EUREUR)
    RETURN
    END 

        ! ********************************************************************

            SUBROUTINE sRAMDA_hj (EU,SC,EC,EBU,ETP,SC_RAM)
    
 
        real(kreal) :: SC(3),EU(3),C_RAMUDA(3),C_RAMU(3),EC(3),     &
                        C_K(3),ETP(3),E_CR(3),SC_RAM(3)
     
        integer :: I,ID
        
        real (kreal) :: EBU,EUmax

    DO I=1,3
        E_CR(I)=EU(I)-ETP(I)
    ENDDO
    EUmax=MAX(E_CR(1),E_CR(2),E_CR(3))
    DO I=1,3
        IF(EUmax==E_CR(I)) ID=I
    ENDDO

    IF(EUmax<=0.) THEN
        DO I=1,3
            C_RAMUDA(I)=1.0
        ENDDO
        GOTO 999
    ENDIF

    DO I=1,3
        C_K(I) = ABS(EUmax/EC(I))
    END DO

    DO I=1,3
        C_RAMU(I) = 1/(0.27+0.96*C_K(I)**0.167)
    ENDDO
    IF(     ID==1) THEN
        C_RAMUDA(1)=1.0
        C_RAMUDA(2)=C_RAMU(2)
        C_RAMUDA(3)=C_RAMU(3)
    ELSE IF(ID==2) THEN
        C_RAMUDA(1)=C_RAMU(1)
        C_RAMUDA(2)=1.0
        C_RAMUDA(3)=C_RAMU(3)
    ELSE IF(ID==3) THEN
        C_RAMUDA(1)=C_RAMU(1)
        C_RAMUDA(2)=C_RAMU(2)
        C_RAMUDA(3)=1.0
    ENDIF

    DO I=1,3
        IF(C_RAMUDA(I) > 1.0) C_RAMUDA(I) =1.0
        IF(C_RAMUDA(I) < 0.4) C_RAMUDA(I) =0.4
    ENDDO

999 DO I=1,3
        SC(I)=SC(I)*C_RAMUDA(I)
        SC_RAM(I)=C_RAMUDA(I)
    ENDDO

    END
        ! ********************************************************************
    SUBROUTINE sSressStiff(NST,NIT,BETA,ICONCYC,                        &
                            DEU,EU,SN,                                  &
                            ERC,SRC,EXC,SXC,                            &
                            EBC,SBC,EJ,SJ,EJJ,                          &
                            ERT,SRT,ECP,SCP,                            &
                            ETP, STP,EPC,SPC,                           &
                            EPC1ST,EPT,EEN,SEN,                         &
                            ETMAX,STMAX,ECCR,SCCR,                      &
                            EUTMAX,EUOVER,ICOUNTLOOPT,                  &
                            E41_51,S41_51,ID_CRACK,                     &
                            EUOLD,S_OLD,SNOLD,EPCU,EPCUS,EPCU_,EPCUS_,  &
                            EO,ET, EC,SC,FC,FT,EBU,ECU,EC_100,          &
                            SCOLD,ECOLD,                                &
                            IUL,EPPC,EPPT,EPEC,CCRACK,                  &
                            INCF,ICIC,C_IC4,SENS,VS,                    &
                            IVIRGIN,ELIMIT,MM,IROT,ICR,LCRACK,          &
                             CYCN,EUDL,Mhj)

        ! ********************************************************************
        ! *引用名：sSressStiff
        ! *機�E�E��Eずみから応力と剛性を求めめE
        ! *入力！E  FC :
        ! *出力！EPR12 :
        ! *参老E��E
        ! ********************************************************************
        ! ********************************************************************
        ! *引用名：sSressStiff
        ! *機�E�E��Eずみから応力と剛性を求めめE
        ! *入力！E  FC :
        ! *出力！EPR12 :
        ! *参老E��E
        ! ********************************************************************
! c      INCLUDE 'com_hj.f'  
        real(kreal) ::  XHJ(1),C_Infinity(1)
        
        real(kreal) ::  BETA(3),DEU(6),DEUTEMP(6),EUDL(3),EU(6),SN(6)
        real(kreal) ::  SNCR(6),ET(3),SCOLD(3),ECOLD(3),SC(3)
        real(kreal) ::  EC(3),EPCU_(3),EPCUS_(3),EESC(3),VUE(3)
        real(kreal) ::  PRO(3),EEN(3),SEN(3),ERC(3),SRC(3)
        real(kreal) ::  EXC(3),SXC(3),EBC(3),SBC(3),EJ(3),SJ(3),EJJ(3)
        real(kreal) ::  ERT(3),SRT(3),ECP(3),SCP(3),ETP(3),STP(3),EPC(3),SPC(3)
        real(kreal) ::  EPC1ST(3),EPT(3),ETMAX(3),STMAX(3),ECCR(3),SCCR(3)
        real(kreal) ::  E41_51(3),S41_51(3), EUOLD(6),S_OLD(6),SNOLD(6)
        real(kreal) ::  EPPC(3),EPPT(3),EPEC(3),CCRACK(3),ELIMIT(3)
        real(kreal) ::  EUTMAX(3),EUOVER(3),CYCN(3,49),CYCNDM(49)
        real(kreal) :: RRHL(99)
        real(kreal) :: EPCU,EPCUS,EO,FC,FT,EBU,ECU,EC_100
        real(kreal) :: C_IC4,SENS,VS
        
        real(kreal) :: copx,ebudm,ecr,eudmdm,euolddm,euover1ratio
        real(kreal) :: euovertemp1,euovertemp2
        
        integer :: i,i1,i2,i3,iflgtemp1,inti,ninte,nreal
        
        integer :: IUL(3), ICIC(8),ID_CRACK(3),IVIRGIN(3)
        integer :: IIHL(99),ICOUNTLOOPT(3)
        integer :: NST,NIT,ICONCYC,INCF,MM,IROT,ICR,LCRACK,MHJ
        integer :: LX,LY,LZ,M,NITT
        CHARACTER MRRHL(99)*10,MIIHL(99)*10
        real(kreal) ::S_NaN(99),N_NaN(99)
        CHARACTER*100 M_NaN
        
        !added by shan 12/12/2017
        DEUTEMP=0.0D0



    LX = CYCNDM(31) ; LY = CYCNDM(32) ; LZ = CYCNDM(33)
    NITT = CYCNDM(49) ; M = CYCNDM(30)
! CH---***---
! C ***[単調載荷の場吁E***
! c SC(intI) = FC
    IF (ICONCYC.EQ.0) THEN
        DO intI=1,3
            CALL XSECON (EU(intI),SN(intI),ET(intI),EC(intI),SC(intI),  &
                            FC,FT,EO,ECU,EPCU,EPCUS,BETA(intI),         &
! C     *            DCBEF11,DCBEF12,DCBEF13,
! C     *            DCBEF21,DCBEF22,DCBEF23,
! C     *            DCBEF31,DCBEF32,DCBEF33,
! C     *            DCAFT11,DCAFT12,DCAFT13,
! C     *            DCAFT21,DCAFT22,DCAFT23,
! C     *            DCAFT31,DCAFT32,DCAFT33,
                 CCRACK(intI),IUL(intI),LCRACK,INCF,IROT,               &
                            ICIC(1),ICIC(2),ICIC(3),ICIC(4),ICR,intI)
        IF (EC(intI).gt.0.0) THEN 
            CONTINUE
        ENDIF
        ENDDO

        CALL sTransArrayReal(SN,SNCR,3)
        
! C ***[繰返し載荷の場吁E***
    
    ELSE IF (ICONCYC.GE.1) THEN
! C
! C     **** DECIDE SC/SCOLD ****
! C     ***[もしSCが変化しなぁE��吁E特徴点を計算し直さない]***
! C     ***[まあまあ、今だってSCの計算してぁE��ぁE��ら変化しなぁE��ろう]***
        ECR=FT/EO
        
        DO intI=1,3
            DEUTEMP(intI)=EU(intI)-EUOLD(intI)
        ENDDO
        
        
        DO intI = 1,3

! C         EUOVER=0.0
            IFLGTemp1=0
            EUOver1Ratio=0.94/6.3
            EUOverTemp1=EUOver(intI)

            IF(EU(intI).GT.EUTMAX(intI).AND.EU(intI).GT.ECR)THEN
                EUTMAX(intI)=EU(intI)
! C             IFLGTemp1=1
            ENDIF

            IF(IFLGTemp1.EQ.0)THEN
                IF(CCRACK(intI).EQ.4.0)THEN
                    EUOverTemp2=(EUTMAX(intI)-ECR)*EUOver1Ratio*1.
                ELSEIF(CCRACK(intI).EQ.5.0)THEN
                    EUOverTemp2=(EUTMAX(intI)-ECCR(intI))*EUOver1Ratio*0.0
                ENDIF

                IF(EUOverTemp2.LT.0.0)THEN
                    EUOverTemp2=0.0
                ENDIF
            ENDIF
! CC
            IF(EUOverTemp2.GT.EUOverTemp1)THEN
                EUOVER(intI)=EUOverTemp2
            ENDIF

            EBUDM=0.002D0+EUOVER(intI)
            IF(EU(intI).GT.EBUDM)THEN
                EUOverTemp2=0.0
            ELSEIF(EUOLD(intI).GT.EBUDM)THEN
! C             EUOverTemp2=EUOVER(intI)+EUOverTemp2*0.5**ICOUNTLOOPT(intI)
! C             EUTMAX(intI)=0.0    
                ICOUNTLOOPT(intI)=ICOUNTLOOPT(intI)+1           
            ENDIF

            EUOVER(intI)=0.0

            IF(EUOverTemp2.GT.0.0.AND.intI.EQ.1)THEN
                CONTINUE
            ENDIF
            EUOverTemp2=0.0

            EUDMDM=EU(intI)-EUOVER(intI)
            EUOLDDM=EUOLD(intI)-EUOVER(intI)
! C         EUDMDM=EU(intI)
! C         EUOLDDM=EUOLD(intI)

            IF(IUL(intI).EQ.31)THEN
                CONTINUE
            ENDIF

            DO I=1,49 ; CYCNDM(I)=CYCN(intI,I) ; ENDDO
            CYCNDM(34)=intI

! C         ***CHANGE FOR REDUCTION***
! C WRITE(*,*)EBUDM
! CH----------
    EBUDM=EBU ! NO コメンチE:BY HJ
! C ICONCYC=1
! CH----------  

             IF( ICONCYC == 2 ) GOTO 30              ! by sun need check
! C         ***[SCの変化より特徴点の計算し直し]***
            CALL sChangeForReduc(FC,FT,EO,ECU,ECR,ICIC(1),NST,              &
                            EPCU,EPCUS,EPCU_(intI),EPCUS_(intI),            &
                            SC(intI),EC(intI),SCOLD(intI),                  &
                            EEN(intI),SEN(intI),ECP(intI),SCP(intI),        &
                            ETP(intI),STP(intI),EPC(intI),SPC(intI),        &
                            EJ(intI),SJ(intI),EXC(intI),SXC(intI),          &
                        ETMAX(intI),EPEC(intI),CCRACK(intI),IUL(intI),      &
                            NITT,M,LX,LY,LZ,intI)
            I=1
            IF(I==1) GOTO 10
            CALL XSECONCYC (DEUTEMP(intI),EUDMDM,SN(intI),                  &
                        ERC(intI),SRC(intI),EXC(intI),SXC(intI),            &
                        EBC(intI),SBC(intI),EJ(intI),SJ(intI),EJJ(intI),    &
                        ERT(intI),SRT(intI),ECP(intI),SCP(intI),            &
                        ETP(intI), STP(intI),EPC(intI),SPC(intI),           &
                        EPC1ST(intI),EPT(intI),EEN(intI),SEN(intI),         &
                        ETMAX(intI),STMAX(intI),ECCR(intI),SCCR(intI),      &
                        E41_51(intI),S41_51(intI),                          &
! C                     EUOLD(intI),SNOLD(intI),EPCU,EPCUS,                 &
                        EUOLDDM,SNOLD(intI),EPCU,EPCUS,                     &
                        EO,ET(intI), EC(intI),SC(intI),FC,FT,EBUDM,ECU,     &
                        IUL(intI),EPPC(intI),EPPT(intI),  EPEC(intI),       &
                        CCRACK(intI),INCF,ICIC,SENS,VS,                     &
                        IVIRGIN(intI),ELIMIT(intI),MM,IROT,ICR,             &
                        LCRACK)
            GOTO 20
10          CALL XSECONCYC_hj (DEUTEMP(intI),EUDMDM,SN(intI),               &
                    ERC(intI),SRC(intI),EXC(intI),SXC(intI),                &
                    EBC(intI),SBC(intI),EJ(intI),SJ(intI),EJJ(intI),        &
                    ERT(intI),SRT(intI),ECP(intI),SCP(intI),                &
                    ETP(intI), STP(intI),EPC(intI),SPC(intI),               &
                    EPC1ST(intI),EPT(intI),EEN(intI),SEN(intI),             &
                    ETMAX(intI),STMAX(intI),ECCR(intI),SCCR(intI),          &
                    E41_51(intI),S41_51(intI),                              &
                    EUOLDDM,SNOLD(intI),EPCU,EPCUS,                         &
                    EO,ET(intI), EC(intI),SC(intI),FC,FT,EBUDM,ECU,         &
                    EPCU_(intI),EPCUS_(intI),                               &
                    IUL(intI),EPPC(intI),EPPT(intI),  EPEC(intI),           &
                    CCRACK(intI),INCF,ICIC,C_IC4,SENS,VS,                   &
                    IVIRGIN(intI),ELIMIT(intI),MM,IROT,ICR,                 &
                    LCRACK,                                                 &
                     CYCNDM, EUDL(intI),COPX,NIT)
                
                            
            
                     
            DO I=1,49 ; CYCN(intI,I)=CYCNDM(I) ; ENDDO
            CYCN(intI,13) = EXC(intI)
            CYCN(intI,14) = SXC(intI)
            CYCN(intI,22) = EPC(intI)
            CYCN(intI,23) = SPC(intI)
! C         
        
20          GOTO 40
! CH---エラーチェチE��---
    M_NaN='☁EXSECONCYC_hj: [ SN EEN SEN ]  NST  NIT  LX  LY  LZ  intI '
    N_NaN(1) = NST ; N_NaN(2) = NITT
    N_NaN(3) = LX ; N_NaN(4) = LY ; N_NaN(5) = LZ ; N_NaN(6) = intI
    S_NaN(1) = SN(intI) ; S_NaN(2) = EEN(intI) ; S_NaN(3) = SEN(intI)
    ! CALL CCHJ(N_NaN, 6, S_NaN, 3 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
    LX = CYCNDM(31) ; LY = CYCNDM(32) ; LZ = CYCNDM(33)
    NITT = CYCNDM(49) ; M = CYCNDM(30)
! CH-------------PARA-------------------------------------
            I=0
            I=I+1 ; IIHL(I) = CYCNDM( 1) ; MIIHL(I) = 'NST' 
            I=I+1 ; IIHL(I) = CYCNDM(49) ; MIIHL(I) = 'NIT'
            I=I+1 ; IIHL(I) = CYCNDM(30) ; MIIHL(I) = 'M_element'
            I=I+1 ; IIHL(I) = LX*100+LY*10+LZ ; MIIHL(I) = 'In'
            I=I+1 ; IIHL(I) = CYCNDM(34) ; MIIHL(I) = 'intI'
            I=I+1 ; IIHL(I) = CYCNDM(18)*2. ; MIIHL(I) = 'IVCYCN'
            I=I+1 ; IIHL(I) = CYCNDM(19)    ; MIIHL(I) = '転換回数'
            Ninte = I

            I=0
            I=I+1 ; RRHL(I) = CYCNDM(1)+CYCNDM(49)/10. ;MRRHL(I)='Ns'
            I=I+1 ; RRHL(I) = IUL(intI)
                    IF( RRHL(I) > 9. ) RRHL(I) = RRHL(I)/10. + 0.05
                                              MRRHL(I) = 'IUL'
            I=I+1 ; RRHL(I) = EU (intI)     ; MRRHL(I) = 'EU_εu'
            I=I+1 ; RRHL(I) = EUOLDDM       ; MRRHL(I) = 'EUOLD'
            I=I+1 ; RRHL(I) = SN (intI)     ; MRRHL(I) = 'SN_ ρE'
            I=I+1 ; RRHL(I) = SNOLD (intI)      ; MRRHL(I) = 'SNOLD'
            I=I+1 ; RRHL(I) = DEUTEMP(intI) ; MRRHL(I) = 'DEU_dεu'
            I=I+1 ; RRHL(I) = EUDL(intI)    ; MRRHL(I) = 'DEUL' 
            I=I+1 ; RRHL(I) = EXC (intI)    ; MRRHL(I) = 'EXC_εx'
            I=I+1 ; RRHL(I) = SXC (intI)    ; MRRHL(I) = 'SXC_σx'
            I=I+1 ; RRHL(I) = EPC (intI)    ; MRRHL(I) = 'EPC_εp'
            I=I+1 ; RRHL(I) = SPC (intI)    ; MRRHL(I) = 'SPC_σp'
            I=I+1 ; RRHL(I) = COPX          ; MRRHL(I) = 'COPX_Epx'
            I=I+1 ; RRHL(I) = CYCNDM(24)    ; MRRHL(I) = 'COPXL'
            I=I+1 ; RRHL(I) = CYCNDM(25)    ; MRRHL(I) = 'EXC0'
            I=I+1 ; RRHL(I) = CYCNDM(26)    ; MRRHL(I) = 'SXC0'
            I=I+1 ; RRHL(I) = CYCNDM(27)-1  ; MRRHL(I) = 'NCYCN'
            I=I+1 ; RRHL(I) = CYCN(intI,35)/FC  ; MRRHL(I) = 'SC_L'
            I=I+1 ; RRHL(I) = CYCN(intI,36) ; MRRHL(I) = 'EC_L'
            I=I+1 ; RRHL(I) = CYCN(intI,37) ; MRRHL(I) = 'LCRA'
            I=I+1 ; RRHL(I) = SC (intI)/FC      ; MRRHL(I) = 'SC'
            I=I+1 ; RRHL(I) = EC (intI)     ; MRRHL(I) = 'EC'
            I=I+1 ; RRHL(I) = EEN (intI)    ; MRRHL(I) = 'EEN'
            I=I+1 ; RRHL(I) = SEN (intI)    ; MRRHL(I) = 'SEN'
            I=I+1 ; RRHL(I) = EPC (intI)    ; MRRHL(I) = 'EPC'
            I=I+1 ; RRHL(I) = ECP (intI)    ; MRRHL(I) = 'ECP'
            Nreal = I

    Mhj=XHJ(1)
    DO I1=1,2
    DO I2=1,2
    DO I3=1,2
    I=I1*1000+I2*100+I3*10
! c        CALL HJ_CHECKpara_1(NITT,M,LX,LY,LZ,intI,I,Mhj,I1,I2,I3,
! c     +                       RRHL,Nreal,IIHL,Ninte,MRRHL,MIIHL)
    ENDDO
    ENDDO
    ENDDO
! C
! CH-----------------------------------------------------------------------
! C
! C10           CALL XSECONCYC_hj (EUDMDM,                           
! C     *                                                     
! C     *                                                                          
! C     *               EPC1ST(intI),                     
! C     *               ECCR(intI),SCCR(intI),
! C     *               EUOLDDM,SNOLD(intI),EPCU,EPCUS,
! C     *                ECU,
! C     *                         EPPC(intI),EPPT(intI),  EPEC(intI),
! C     *               CCRACK(intI),                SENS,VS,
! C     *               IVIRGIN(intI),ELIMIT(intI),MM,IROT,ICR,
! C     *               LCRACK,
! C     *                CYCNDM, EUDL(intI),COPX)
 30 CALL CONCYC_hj( FC,FT,EO,EC_100,IUL(intI),EU(intI),DEUTEMP(intI),           &
                S_OLD(intI) , SNOLD(intI) ,                                     &
                EC(intI) , SC(intI) , ETP(intI) , STP(intI) ,                   &
                EEN(intI) , SEN(intI), EXC(intI) , SXC(intI),                   &
                E41_51(intI),S41_51(intI),ETMAX(intI),STMAX(intI),              &
                EPCU_(intI),EPCUS_(intI),ECP(intI) , SCP(intI) ,                &
                EJ(intI) , SJ(intI), ERT(intI) , SRT(intI) ,                    &
                ERC(intI) , SRC(intI) ,EJJ(intI),EPT(intI),                     &
                EPC(intI) , EBC(intI) , SBC(intI), SPC(intI) ,                  &
                ID_CRACK(intI),ICIC,C_IC4,                                      &
                SN(intI),ET(intI),INCF,                                         &
                NST , NITT , M ,LX,LY,LZ ,intI)                                 
                    EUOLDDM=S_OLD(intI)                                         
                                                                                
            I=0
            I=I+1 ; IIHL(I) = CYCNDM( 1) ; MIIHL(I) = 'NST' 
            I=I+1 ; IIHL(I) = IUL(intI) ; MIIHL(I) = 'IUL' 
            I=I+1 ; IIHL(I) = ID_CRACK(intI) ; MIIHL(I) = 'ID_CRACK'
            Ninte = I

            I=0
            I=I+1 ; RRHL(I)= CYCNDM(1)+CYCNDM(49)/10. ;MRRHL(I)='Ns'
            I=I+1 ; RRHL(I)= EU(intI)     ;MRRHL(I)= 'S'
            I=I+1 ; RRHL(I)= SN(intI)     ;MRRHL(I)= 'T'
            I=I+1 ; RRHL(I)= ET(intI)     ;MRRHL(I)= 'E'
            I=I+1 ; RRHL(I)= DEUTEMP(intI);MRRHL(I)= 'DS'
            I=I+1 ; RRHL(I)= EC(intI)     ;MRRHL(I)= 'O_S'
            I=I+1 ; RRHL(I)= SC(intI)     ;MRRHL(I)= 'O_T'
            I=I+1 ; RRHL(I)= ETP(intI)    ;MRRHL(I)= 'T_S'
            I=I+1 ; RRHL(I)= STP(intI)    ;MRRHL(I)= 'T_T'
            I=I+1 ; RRHL(I)= EEN(intI)    ;MRRHL(I)= 'E_S'
            I=I+1 ; RRHL(I)= SEN(intI)    ;MRRHL(I)= 'E_T'
            I=I+1 ; RRHL(I)= E41_51(intI) ;MRRHL(I)= 'F_S'
            I=I+1 ; RRHL(I)= S41_51(intI) ;MRRHL(I)= 'F_T'
! C         I=I+1 ; RRHL(I)= ETMAX(intI)  ;MRRHL(I)= 'Y_S'
! C         I=I+1 ; RRHL(I)= STMAX(intI)  ;MRRHL(I)= 'Y_T'
            I=I+1 ; RRHL(I)= EPCU_(intI)  ;MRRHL(I)= 'U_S'
            I=I+1 ; RRHL(I)= EPCUS_(intI) ;MRRHL(I)= 'U_T'
            I=I+1 ; RRHL(I)= EJ(intI)     ;MRRHL(I)= 'V_S'
            I=I+1 ; RRHL(I)= SJ(intI)     ;MRRHL(I)= 'V_T'
            I=I+1 ; RRHL(I)= ERT(intI)     ;MRRHL(I)= 'W_S'
            I=I+1 ; RRHL(I)= SRT(intI)     ;MRRHL(I)= 'W_T'
! C         I=I+1 ; RRHL(I)= ECP(intI) ;MRRHL(I)= 'C_S'
! C         I=I+1 ; RRHL(I)= SCP(intI) ;MRRHL(I)= 'C_T'
! C         I=I+1 ; RRHL(I)= EXC(intI) ;MRRHL(I)= 'X_S'
! C         I=I+1 ; RRHL(I)= SXC(intI) ;MRRHL(I)= 'X_T'
            I=I+1 ; RRHL(I)= ERC(intI) ;MRRHL(I)= 'R1_S'
            I=I+1 ; RRHL(I)= SRC(intI) ;MRRHL(I)= 'R1_T'
            I=I+1 ; RRHL(I)= EJJ(intI) ;MRRHL(I)= 'R2_S'
            I=I+1 ; RRHL(I)= EPT(intI) ;MRRHL(I)= 'R2_T'
            I=I+1 ; RRHL(I)= EPC(intI) ;MRRHL(I)= 'P_S'
            I=I+1 ; RRHL(I)= EBC(intI) ;MRRHL(I)= 'H_S'
! C         I=I+1 ; RRHL(I)= SBC(intI) ;MRRHL(I)= 'G_S'
! C         I=I+1 ; RRHL(I)= SPC(intI) ;MRRHL(I)= 'CR_S'
! C         I=I+1 ; RRHL(I)=  ;MRRHL(I)= ''
            Nreal = I
! C
    Mhj=XHJ(1)
    DO I1=1,2
    DO I2=1,2
    DO I3=1,2
    I=I1*1000+I2*100+I3*10
! c        CALL HJ_CHECKpara_1(NITT,M,LX,LY,LZ,intI,I,Mhj,I1,I2,I3,
! c     +                       RRHL,Nreal,IIHL,Ninte,MRRHL,MIIHL)
    ENDDO
    ENDDO
    ENDDO
! C
40  CONTINUE
! C NST=CYCNDM(1)

! CH---------------------------------------------------
        ENDDO

    END IF
! C IF(NST==215) THEN
! C      WRITE(98,'(2I5,3I1,3es12.3)') NST,NIT,LX,LY,LZ,(SC(I),I=1,3)
! C ENDIF
! C
    END

    !***********************************************************************
      SUBROUTINE XSECON (EU,SN,ET,EC,SC,FC,FT,EO,ECU,EPCU,EPCUS,BETA,       &
! C     *             DCC11,DCC12,DCC13,
! C     *             DCC21,DCC22,DCC23,
! C     *             DCC31,DCC32,DCC33,
! C     *             DCR11,DCR12,DCR13,
! C     *             DCR21,DCR22,DCR23,
! C     *             DCR31,DCR32,DCR33,
                  CCRACK,IUL,NIUL,INCF,IROT,                                &
                             ICIC1,ICIC2,ICIC3,ICIC4,ICR,N)
! C***********************************************************************
! C
! C *** STRESS - STRAIN RELATIONSHIP FOR CONCRETE
! C     ORIGINAL PROGRAM = SUBROUTINE SECONC(CYCLIC RULE)
! C     THIS PROGRAM IS MODIFIED FOR MONOTONIC LOADING IN ORDER TO
! C     APPLY TO CALCULATE UNBLANCED STRESS < 1991.8.9 K.UCHIDA & A.AMEMIY
! C
 
        real (kreal) :: EU,SN,ET,EC,SC,FC,FT,EO,ECU,EPCU,EPCUS,BETA,CCRACK

        integer  :: IUL,NIUL,INCF,IROT,ICIC1,ICIC2,ICIC3,ICIC4,ICR,N
 
        real (kreal) :: dum,dum1,ebu,ecr,eex1,emax,epcuec,epcusc,es
        real (kreal) :: smax,ssx1
 ! C
! C *** BOND LIMIT STRAIN
! C
      EBU=0.002D0
! c IF(ICIC(4)==3.OR.ICIC(4)==4)  EBU=2000.0
! C
! C *** CRACK STRAIN
! C
      ECR=FT/EO
! C
! C     ELASTIC (TENSION & COMPRESSION)  IUL = 1 ------------ >  1
! C     TENSION STIFENNIG                IUL = 4 ------------ >  4
! C     COMPRESSIVE SOFTENNING           IUL = 8 ------------ >  8
! C     COMPRESSIVE SOFTENNING (EU>EPCU) IUL = 9 ------------ >  9
! C
! C     GO TO (1,2,3,4,5,6,7,8,9),IUL
! C
! CCCCCCCC
! C       WRITE(6,*) 'EU=',EU
! C       WRITE(6,*) 'EC=',EC
! C      WRITE(6,*) 'ECR=',ECR
! C      WRITE(6,*) 'EPCU=',EPCU
! C       WRITE(6,*) 'SC=',SC
! CCCCCCC
! C      IF(EC.GE.0.0D0) THEN
! C        IF(EU.LT.ECR)                GOTO 1
! C        IF(EU.GE.ECR)                GOTO 4
! C      ELSEIF
! CCCCCCC
      IF(EU.GE.0.0D0) THEN
        IF(EU.LT.ECR)                GOTO 1
        IF(EU.GE.ECR)                GOTO 4
      ELSE

        IF(SC.LT.FC) THEN
          SMAX=SC
        ELSE
          SMAX=FC/BETA
        ENDIF
        IF(EC.LT.ECU) THEN
          EMAX=EC
        ELSE
          EMAX=ECU
        ENDIF
! C
! CH---***---
! c SMAX=FC*1.0 ; EMAX=ECU
! CH---***---

! C
! C *** COMPRESSION LIMIT
! C 
        EPCUEC=EPCU*EMAX/ECU
        EPCUSC=EPCUS*SMAX/FC
! C
! C*******************************
        IF(EU.GE.EC)                  GOTO 1
        IF(EU.LT.EC.AND.EU.GE.EPCUEC) GOTO 8
        IF(EU.LT.EPCUEC)              GOTO 9
      ENDIF
! C
! C *** IUL=1
! C
    1 CONTINUE
! CCC
! C      WRITE(6,*)'****IUL=1******'
! CCC
      IUL=1

      IF(EU)  101,100,100
! C
! C  ** IN TENSION ZONE
! C
  100 CONTINUE
      ET=EO
      SN=ET*EU
! CCCC
! C      IF(SC.LE.0.0D0) SC=FT
! C      IF(EC.LE.0.0D0) EC=ECR
! CCCC
! C       WRITE(6,*) 'SN=',SN
! CCCC
! CC      IF(SN.GT.FT) WRITE(6,6001) IUL
      RETURN
! C
! C  ** IN COMPRESSION ZONE
! C
  101 CONTINUE
! CCCCCCCCCC
! C      IF(SC.GT.0.0D0) SC=FC
! C      IF(EC.GT.0.0D0) EC=ECU
! CCCCCCCCCC
      IF(ICIC1.EQ.0) THEN
! C
! C     *** SAENZ EQUATION ***
! C
! CCCCCCCCCCCCCCCCCC
! C     ET=EO
! C     SN=ET*EU
! CCCCCCCCCCCCCCCCCCC
      DUM=1.0D0+(EO*EMAX/SMAX-2.0D0)*EU/EMAX+(EU/EMAX)**2
      ET=EO*(1.0D0-(EU/EMAX)**2)/DUM**2
      SN=EO*EU/DUM
! CCCCCCCCCCCCCCCCC
      IF(SN.LT.SC) THEN
      EC=((EO*EMAX*EMAX+2.0D0*SC*EMAX-SC*EO*EMAX*EMAX/SMAX)-            &
         SQRT((EO*EMAX*EMAX+2.0D0*SC*EMAX-SC*EO*EMAX*EMAX/SMAX)**2-     &
         4.0D0*SC*SC*EMAX*EMAX))*0.5D0/SC
      GOTO 8
      ENDIF
! CCCCCCCCCCCCCCCCCC
! C       WRITE(6,*) 'SN=',SN
! CCC
! C
! C     *** FAFITIS-SHAH EQUATION ***
! C
      ELSE
      DUM1=EU/EC
      DUM=EO*EC/SC
      ET=EO*(1.0D0-DUM1)**(DUM-1.0D0)
      SN=SC*(1.0D0-(1.0D0-DUM1)**DUM)

      END IF
      IF(ET.LT.EO/100.0D0) ET=EO/100.0D0
! CC      IF(SN.LT.SC) WRITE(6,6001) IUL
      RETURN
! C
! C *** IUL=4
! C
    4 CONTINUE
! CCC
! C      WRITE(6,*)'****IUL=4****'
! CCC
    CCRACK=4.0
      IUL=4
      INCF=3
    ICR = 1
      ET=EO/100.0D0
! C
! C     <IROT=0>                  : EQUIVALENT UNIAXIAL STRAIN MODEL
! C     <IROT=1 & ICR =0>         : EQUIVALENT --->> CRACK
! C     <IROT=1 & ICR =1>         : CRACK DIRECTION CORDINATE MODEL
! C
! C
! C     < ICIC4=0 > : TENSION STIFENNING (SIRAI EQUATION)
! C     < ICIC4=1 > : NONE
! C     < ICIC4=2 > : TENSION CUT OFF
! C
  104 CONTINUE

      IF(ICIC4.EQ.0.OR.ICIC4.EQ.1) THEN
      IF(EU.LE.EBU) THEN
       EEX1=(EU-ECR)/(EBU-ECR)
       SSX1=1.0D0-2.748D0*EEX1+2.654D0*EEX1**2-0.906D0*EEX1**3
       SN=SSX1*FT
      ELSE
       SN=0.0D0
      END IF
      END IF

      IF(ICIC4.EQ.2) THEN
       SN=0.0D0
      END IF
! CCCCCC
! c      ET=EO
! c      SN=ET*EU
! CCCCCC
      RETURN
! C
! C *** IUL=8
! C
    8 CONTINUE
! CCC
! C      WRITE(6,*) ' *****IUL=8******'
! CCC
      IUL=8
      INCF=3
      ET=EO/100.0D0
! CCCCCCCCCCCCCCCCCCCC
! C      IF(EPCUS.GT.0.0D0) THEN
! C      SN=SC+ET*(EU-EC)
! C         IF(SC.LT.SN) THEN
! C           WRITE(6,*) 'SC=',SC
! C           WRITE(6,*) 'SN=',SN
! C           stop
! C         ENDIF
! C      ELSE
! C      ES=(EPCUS-FC)/(EPCU-ECU)
! C      SN=SC+ES*(EU-EC)
! C      ENDIF
! C
      IF(EPCUS.GT.0.0D0) THEN
      SN=SC+ET*(EU-EC)
         IF(SC.LT.SN) THEN
           WRITE(6,*) 'SC=',SC
           WRITE(6,*) 'SN=',SN
           stop
         ENDIF
      ELSE
! C        EPCUEC=EPCU*EC/ECU
! C        EPCUSC=EPCUS*SMAX/FC
        es=(epcusc-Fc)/(epcuec-ecU)
        sn=sc+es*(eu-ec)
! C        ES=(EPCUSC-SMAX)/(EPCUEC-EMAX)
! C        SN=SMAX+ES*(EU-EC)
      ENDIF
! C
! CCCCCCCCCCCCCCCCCCC
! C
! C     WRITE(6,*) '(IUL=8) SC=',SC
! C     WRITE(6,*) '(IUL=8) SN=',SN
! C     WRITE(6,*) '(IUL=8) EC=',EC
! C     WRITE(6,*) '(IUL=8) ES=',ES
      RETURN
! C
! C *** IUL=9
! C
    9 CONTINUE
! CCC
! C      WRITE(6,*) '*****IUL=9********'
! CCC
      IUL=9
      INCF=3
      ET=EO/100.0D0
      SN=EPCUS
      RETURN
! C
! C *** ERROR MESSAGE
! C
 6001 FORMAT(' *** ERROR IN SECON3 (IUL=',I2,') ***')

      END
    
    !***********************************************************************
    SUBROUTINE sChangeForReduc(FC,FT,EO,ECU,ECR,ICIC1,NST,      &
                            EPCU,EPCUS,EPCU_,EPCUS_,            &
                            SC,EC,SCOLD,                        &
                            EEN,SEN,ECP,SCP,                    &
                            ETP,STP,EPC,SPC,                    &
                            EJ,SJ,EXC,SXC,                      &
                            ETMAX,EPEC,CCRACK,IUL,              &
                            NIT,M,LX,LY,LZ,intI)                
                                                                
     ! ********************************************************************
! *引用名：sChangeForReduc                                        
! *機�E�E�SCの変化より特徴点の計算し直ぁE                              
! *入力！E  FC :                                               
! *出力！EPR12 :                                           
! *参老E��E                                                       
! ********************************************************************
     
     
    real (kreal) :: FC,FT,EO,ECU,ECR,EPCU,EPCUS,EPCU_,EPCUS_
    real (kreal) :: SC,EC,SCOLD,EEN,SEN,ECP,SCP,ETP,STP,EPC,SPC
    real (kreal) :: EJ,SJ,EXC,SXC,ETMAX,EPEC,CCRACK
    real (kreal) :: DUMMY,ES,ESDM,PRO,ECCR,EPPC,SCCR
    integer :: ICIC1,NST,NIT,M,LX,LY,LZ,intI,IUL
    integer :: ncall
! C
! CH---圧壁E収斂)点応力ひずみの計箁E----
! C こ�E部刁E�E本来は�E�PRO /= 1.の時，実行すべぁE
    ES = (FC-EPCUS) / (ECU-EPCU)
    EPCUS_= EPCUS
    EPCU_ = EC-(SC-EPCUS)/ES
! C EPCUS = EPCUS_
! C EPCU  = EPCU_
! CH------------------------------
! C
    IF(NST.EQ.1)THEN
        PRO=1.0
    ELSE
        PRO=SC/SCOLD
    ENDIF
    IF(PRO.NE.1.0) THEN
! CH---圧壁E収斂)点応力ひずみの計箁E----
! CH------------------------------
! C         ***圧縮域で除荷したことある場吁E**
        IF(EEN.NE.0.0) THEN
            IF(EEN.GT.ECU)THEN
                CALL CompressionZone(ICIC1,EO,DUMMY,ECU,SC,EEN,SEN)
            ELSE    
                ESDM=(EPCUS_-SC)/(EPCU_-ECU)
                SEN=SC+ESDM*(EEN-ECU)   
            ENDIF
! C
! c         CALL ECTPJ_hj(EEN,SEN,ECP,SCP,
! c     *                   ETP,STP,EPC,SPC,
! c     *                   EJ,SJ,ECCR, SCCR,EPCU_, EPPC,
! c     *                   EO,SC,FT,EC,CCRACK,ICIC1,
! c     *                   EC,ETMAX,ECR,EPEC,IVCYCN ,NCYCNOFF,
! c     +           NST,NIT,M,LX,LY,LZ,intI,NCALL )
            NCALL = 99999 
            CALL ECTPJ(EEN,SEN,ECP,SCP,                 &
                        ETP,STP,EPC,SPC,                &
                        EJ,SJ,ECCR, SCCR,EPCU_, EPPC,   &
                        EO,SC,FT,EC,CCRACK,ICIC1,       &
                        EC,ETMAX,ECR,EPEC,              &
                NST,NIT,M,LX,LY,LZ,intI,NCALL )         !FROM KOBAYASHI
                                     
            CALL CROSS(EEN,SEN,EPPC,ECP,SCP,            &
                        EXC,SXC,EO,SC,EPC,              &
                        ECU,IUL,EPCU_,EPCUS_,ICIC1) 
            ! CALL ROOT_CHECK(2611, 10)

        ELSE 
! C         ***圧縮域で除荷したことなぁE��吁E**
            CALL CompressionZone(ICIC1,EO,DUMMY,ECU,SC,EJ,SJ)
            !CALL ROOT_CHECK(2610, 10)
        ENDIF

    ENDIF
    END

    !*******************************************************************
      SUBROUTINE CompressionZone(ICIC1,EO,DUMMY,ECU,SC,EEN,SEN)
! *******************************************************************
! ****                          hj . 2008.8.30                   ****
! *******************************************************************
    real (kreal) :: EO,DUMMY,ECU,SC,EEN,SEN
    integer :: ICIC1
    
    IF(     ICIC1 == 1) THEN
        CALL FS(EO,DUMMY,ECU,SC,EEN,SEN)
    ELSE IF(ICIC1 == 2) THEN
        CALL Kimura(EO,DUMMY,ECU,SC,EEN,SEN)
    ELSE
        CALL SAENZ(EO,DUMMY,ECU,SC,EEN,SEN)
    ENDIF
    END

    !*******************************************************************
      SUBROUTINE Kimura( E0,E,O_S,O_T,S,T)
!*******************************************************************

! c S: ひずみ
! C T: 応力
! C
    real (kreal) :: E0,E,O_S,O_T,S,T,Esec,r,X,g,gd,h,hd

    Esec = O_T / O_S
    r = E0 / ( E0 - Esec )
    X = S / O_S

    g = r - 1. + X**r
    gd = X**r / S * r
    h = O_T * X * r
    hd = Esec * r

    E = ( hd * g - h * gd ) / g**2
    T = h / g

      END
    
    !*******************************************************************
      SUBROUTINE SAENZ(EO,ET,ECU,SC,EPS,SIG)
! *******************************************************************
! ****  ADDED BY JOE ZAKKY ON 1996.8.14 FOR BECOMING MORE EASILY ****
! **** TO ADJUST, CHECK OR CHANGE INTO OTHER EQUATIONS.          ****
! ! *******************************************************************
! *引数�E�　 EO  �E�弾性初期接線係数
! *     ET  �E�当応力点で接線剛性    
! *     ECU �E�最大圧縮応力時�E等価一軸ひずみ 
! *     SC  �E�最大圧縮応力   
! *     EPS �E�条件としてのひずみ    
! *     SIG �E�求める応力        
! ! *******************************************************************
! *$INCLUDE"D:\ANALYSIS\ORIGINAL\TS\H\FILE.H"
! *$INCLUDE'/USER2/ZAKKY/Z9/H/FILE.H'
      !CALL ROOT_CHECK(0000, 56)
        real (kreal) :: EO,ET,ECU,SC,EPS,SIG,Z,SIGLINE
      Z  =1.+(EO*ECU/SC-2.)*EPS/ECU+(EPS/ECU)**2
      SIG=EO*EPS/Z
      ET =EO*(1.-(EPS/ECU)**2)/Z**2
        SIGLINE=EPS*EO
        IF(SIG.LT.SIGLINE) THEN
            SIG=SIGLINE
            ET=EO
        ENDIF
      !CALL ROOT_CHECK(9999, 56)
      RETURN
      END

    !*******************************************************************
      SUBROUTINE ECTPJ(  EEN,  SEN,  ECP,  SCP,             &
                        ETP,  STP,  EPC,  SPC,              &
                         EJ,   SJ, ECCR, SCCR,              &
                       EPCU, EPPC,   EO,   FC,              &
                         FT,   ECU,CRACK, ICIC,             &
                         EC,ETMAX,  ECR, EPEC,              &
                NST,NIT,M,LX,LY,LZ,intI,NCALL )
! *******************************************************************
! ****  SUBROUTINE FOR SOLVING POINT LOCATION OF                 ***
! ****  E, C, T, P, J AND T'                                     ****
! *******************************************************************

! *$INCLUDE"D:\ANALYSIS\ORIGINAL\TS\H\FILE.H"
! *$INCLUDE'/USER2/ZAKKY/Z9/H/FILE.H'
        real(kreal) :: EEN,SEN,ECP,SCP,ETP,STP,EPC,SPC,EJ,SJ,ECCR,SCCR
        real(kreal) :: EPCU,EPPC,EO,FC,FT,ECU,CRACK,EC,ETMAX,ECR,EPEC
        integer ::  ICIC,NST,NIT,M,LX,LY,LZ,intI,NCALL
        
        real(kreal) ::aaaaaa,alfa,bbb,epecc,epep
        integer id_nan,n
        
        real(kreal) S_NaN(99),N_NaN(99)
        CHARACTER*100 M_NaN

! CH---老E��EECTPJ------------------------------------------------
    IF(EEN >= 0.0) THEN
        ! WRITE(1850,'(A72,8I6,ES11.3)')'●●●老E��_ECTPJ [ EEN ] :
     ! + ≧0.0,EEN NCALL NST NIT M LX LY LZ intI'
     ! +                ,NCALL,NST,NIT,M,LX,LY,LZ,intI,EEN
        GOTO 9999
    ENDIF
    M_NaN='☁EECTPJ [ EEN ] : NCALL NST NIT M LX LY LZ intI'
    N_NaN(1) = NCALL ;N_NaN(2) =  NST ;N_NaN(3) =  NIT ;N_NaN(4) =  M
    N_NaN(5) = LX  ;N_NaN(6) = LY  ;N_NaN(7) = LZ  ;N_NaN(8) = intI
    S_NaN(1) = EEN
    ID_NaN=0
    !CALL CCHJ(N_NaN, 8, S_NaN, 1 ,M_NaN,ID_NaN)
    IF(ID_NaN==1) GOTO 9999

! CH---老E��EECTPJ------------------------------------------------
      ! CALL ROOT_CHECK(0000, 23)
      ! CALL ROOT_CHECKR1('EEN', EEN)
      ! CALL ROOT_CHECKR1('SEN', SEN)

      IF(DABS(EEN).LT.DABS(4.*EC)) THEN
        EPC=(0.145*(EEN**2.)/EC)+0.127*EEN            !�E�点ひずみ
        SPC=0.0                                    !�E�点応力
    ELSE
        EPC=(EEN/EC-2.828)*EC/5.                      !�E�点ひずみ
        SPC=0.0                                    !�E�点応力
      ENDIF
      ! CALL ROOT_CHECKR1('EPC', EPC)
      EPEC=1.5*SEN/(EEN-EPC)                     !EC間�E剛性
! CH---修正5_hj------------------------------------------------
    IF(EEN >= EPC) THEN
        ! WRITE(1850,'(A60,8I6,2ES11.3)')
     ! +'●●●修正_51:EEN≧EPC NCALL NST NIT M LX LY LZ intI',
     ! +                NCALL,NST,NIT,M,LX,LY,LZ,intI,EEN,EPC
        EPEP=EO/1.5
    ENDIF
    M_NaN='☁EECTPJ_52:[SEN EEN EPC ] NCALL NST NIT M LX LY LZ intI'
    N_NaN(1) = NCALL ;N_NaN(2) =  NST ;N_NaN(3) =  NIT ;N_NaN(4) =  M
    N_NaN(5) = LX  ;N_NaN(6) = LY  ;N_NaN(7) = LZ  ;N_NaN(8) = intI
    S_NaN(1) = SEN ; S_NaN(2) = EEN ; S_NaN(3) = EPC
    ID_NaN=0
    !CALL CCHJ(N_NaN, 8, S_NaN, 3 ,M_NaN,ID_NaN)
    IF(ID_NaN==1) EPEP=EO/1.5
    IF(ID_NaN==1) EPEC=EO
! CH---修正5_hj------------------------------------------------

      ! CALL ROOT_CHECKR1('EPEC1', EPEC)

      IF(EPEC.GT.EO) EPEC=EO                     !EC間�E剛性は初期剛性以丁E

      ! CALL ROOT_CHECKR1('EPEC2', EPEC)

      EPECC=2*SEN/(EEN-EPC)

      IF(EPEC.GT.EPECC) EPEC=EPECC               ! 点�E�で剛性ぁE以下にならなぁE��ぁE��する処置、E

      IF(EEN.GE.EC) THEN
        SCP=0.833333*SEN                           !�E�点応力1
! c     SCP=1.0*SEN 
    ELSE
        SCP=MIN(0.66667*SEN,SEN-0.166667*FC)       !�E�点応力2E
      ENDIF

      BBB=SEN-EPEC*EEN
! C      ECP=(SCP-BBB)/EPEC                         !�E�点ひずみ
    ECP=EEN-(SEN-SCP)/EPEC
      ETP=EEN-SEN/2.0/EPEC                       !�E�点ひずみ
      STP=SEN-SEN/2.0                            !�E�点応力
      EPPC=SCP/(ECP-EPC)               !IUL=6の剛性
! CH---修正4_hj------------------------------------------------
    IF( ECP >= EPC ) THEN
        ! WRITE(1850,'(A60,8I6,2ES11.3)')
     ! +'●●●修正_41:ECP≧EPC NCALL NST NIT M LX LY LZ intI',
     ! +                NCALL,NST,NIT,M,LX,LY,LZ,intI,ECP,EPC
        EPPC=EO
    ELSE
        EPPC=SCP/(ECP-EPC)               !IUL=6の剛性
    ENDIF
    M_NaN='☁EECTPJ_42:[SCP ECP EPC] NCALL NST NIT M LX LY LZ intI'
    N_NaN(1) = NCALL ;N_NaN(2) =  NST ;N_NaN(3) =  NIT ;N_NaN(4) =  M
    N_NaN(5) = LX  ;N_NaN(6) = LY  ;N_NaN(7) = LZ  ;N_NaN(8) = intI
    S_NaN(1) = SCP ; S_NaN(2) = ECP ; S_NaN(3) = EPC
    ID_NaN=0
    ! CALL CCHJ(N_NaN, 8, S_NaN, 3 ,M_NaN,ID_NaN)
    IF(ID_NaN==1) EPPC=EO
! CH---修正4_hj------------------------------------------------
      ! CALL ROOT_CHECKR1('EPPC', EPPC)
      IF(EO.LT.EPPC)THEN
        EPPC=EO
        EPC=EEN-SEN/EO
      ENDIF

    IF (EPC.GT.0.0) THEN
        EPC=0.0
        EPPC=SEN/EEN
    ENDIF   

      ! CALL ROOT_CHECKR1('EPC1', EPC)
  ! 100 CALL ROOT_CHECK( 100, 23)
      IF(CRACK.EQ.4.0)THEN
        IF(EEN.EQ.0.0)THEN
  ! 200 CALL ROOT_CHECK( 200, 23)
!         SJ=-FT*1.5                             !IUL=2を未経騁E
! C          SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)*15      !FROM KOBAYASHI(MODIFIED BY SAKURAI)
          SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)*1     !FROM KOBAYASHI(MODIFIED BY TIDE)
! *         EJ=SJ/EPPC+EPC                         !FROM KOBAYASHI
          IF(ICIC.EQ.0)THEN
  ! 300 CALL ROOT_CHECK( 300, 23)
            CALL SAENZ_REVERS(EJ,SJ,EO,ECU,FC)    !SAENZ式送E��数
          ELSE IF(ICIC == 1) THEN
  ! 400       CALL ROOT_CHECK( 400, 23)
            CALL SHAH_REVERS (EJ,SJ,EO,ECU,FC)    !FAFITIS-SHAH式送E��数
          ELSE IF(ICIC == 2) THEN
            CALL Kimura_REVERS(EJ,SJ,EO,ECU,FC)
          ENDIF
        ELSE
  ! 500 CALL ROOT_CHECK( 500, 23)
! *         SJ=SEN*0.1                             !IUL=2を経騁E
! C          SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)*15      !FROM KOBAYASHI(MODIFIED BY SAKURAI)
          SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)*1       !FROM KOBAYASHI(MODIFIED BY TIDE)
          EJ=SJ/EPPC+EPC
        ENDIF
      ELSE
  ! 600   CALL ROOT_CHECK( 600, 23)
! C        SJ=SEN*0.1*15                               !IUL=2を経騁E
        SJ=SEN*0.1*1                             !IUL=2を経騁E
        EJ=SJ/EPPC+EPC
      ENDIF
    IF(EJ.GE.0.0)THEN
        AAAAAA=1.
    ENDIF
  ! 700 CALL ROOT_CHECK( 700, 23)
! c      EPCU=-.01
      IF(CRACK.NE.0.0)GOTO 9999                  !圧縮効果を老E�Eしたひび割れ発生点の算�E
      ALFA=0.0
      ALFA=FT/(FT/EO-EPCU)
! C      ECCR=(EPPC*EPC-ALFA*EPCU)/(EPPC-ALFA)
    ECCR=EPC+(1-EPC/EPCU)*(FT/EO)
      SCCR=EPPC*ALFA*(EPCU-EPC)/(ALFA-EPPC)
! CH---修正6------------------------------------------------
    IF( ALFA==EPPC ) THEN
        ! WRITE(1850,'(A60,8I6,2ES11.3)')
     ! +'●●●修正6_hj:ALFA=EPPC NCALL NST NIT M LX LY LZ intI',
     ! +                NCALL,NST,NIT,M,LX,LY,LZ,intI,ALFA,EPPC
        SCCR=FT
    ELSE
    ENDIF
    M_NaN='☁E ECTPJ_hj_62:[ EPPC ALFA EPCU_ EPC ]   &
      NCALL NST NIT M LX LY LZ intI'
    N_NaN(1) = NCALL ;N_NaN(2) =  NST ;N_NaN(3) =  NIT ;N_NaN(4) =  M
    N_NaN(5) = LX  ;N_NaN(6) = LY  ;N_NaN(7) = LZ  ;N_NaN(8) = intI
    S_NaN(1)=EPPC; S_NaN(2)=ALFA ; S_NaN(3) = EPCU;S_NaN(4) = EPC
    N=0
    !CALL CCHJ(N_NaN, 8, S_NaN, 4 ,M_NaN,N)
    IF(N==1) SCCR=FT
! CH---修正6------------------------------------------------
    IF(SJ.LT.-100.0) THEN 
    CONTINUE
    ENDIF
 9999 RETURN
      END

    !***********************************************************************
      SUBROUTINE CROSS (EEN,SEN,EPP,ECP,SCP,EX,SX,EO,SC,EPC,    &
                       EC,IUL,EPCU,EPCUS,ICIC)
! ***********************************************************************
! * 機�E         �E�EX.SXを決めるSUB
! * 引数�E�EEPP       :PC線�E剛性�E�EUL=6の剛性�E�！EPPC�E�E
! *         ECP     :C点の主ひずみ、主応力
! *         SCP     :C点の主ひずみ、主応力    
! *         EX      :X点の主ひずみ、主応力    
! *         SX      :X点の主ひずみ、主応力    
! *         EO      :コンクリート要素の初期剛性
! *         SC      :コンクリート�E三軸応力状態を老E�Eした最大主ひずみ、応力�E�破壊曲面との接点の値�E�E
! *         EPC     :P点の主ひずみ、主応力
! *         EC      :コンクリート�E三軸応力状態を老E�Eした最大主ひずみ、応力�E�破壊曲面との接点の値�E�E
! *         IUL     :コンクリート要素の吁E���E点の2つの等価一軸方向�E応力状態を表す指樁E
! *         EPCU    :コンクリート�E圧壊が生じた後�E収斂点ひずみ
! *         EPCUS   :コンクリート�E圧壊が生じた後�E収斂点応力
! *         ICIC    :圧縮応力上�E曲線オプション　0�E�SAENZ式　1�E�FAFITIS-SHAH弁E
! *         EEN     :E点の主ひずみ
! *         SEN     :E点の主応力
! ***********************************************************************

    real (kreal) :: EEN,SEN,EPP,ECP,SCP,EX,SX,EO,SC,EPC,EC,EPCU,EPCUS
    integer :: IUL,ICIC
    
    real (kreal) :: EY,R,ET
    
      !CALL ROOT_CHECK(0000, 14)

      EY=(SC-SCP)/EPP+ECP
      IF(EY-EC) 10,10,20
! ****X点がIUL=8丁E
   10 R=(SC-EPCUS)/(EC-EPCU)
      EX=(SC-R*EC+EPP*EPC)/(EPP-R)
      SX=R*(EX-EC)+SC
      !CALL ROOT_CHECK(9999, 14)
      RETURN
! ****X点が応力上�E曲線丁E
   20 R=-EC/100.
    CALL CompressionZone(ICIC,EO,ET,EC,SC,EEN,SEN)
    IF(ET.EQ.EPP) THEN
        EX=EEN
    ELSE    
        EX=(ET*EEN-SEN-EPP*EPC)/(ET-EPP)
    ENDIF

    CALL CompressionZone(ICIC,EO,ET,EC,SC,EX,SX)

      RETURN
      END
    
    
    !*******************************************************************
      SUBROUTINE XSECONCYC (    DEU, EU,    SN,                         &
                 ERC,   SRC,    EXC,   SXC,                             &
                 EBC,   SBC,    EJ,    SJ,    EJJ,                      &
                 ERT,   SRT,    ECP,   SCP,                             &
                 ETP,   STP,    EPC,   SPC,                             &
                            EPC1ST,EPT,    EEN,   SEN,                  &
                 ETMAX, STMAX,  ECCR,  SCCR,                            &
                            E41_51,S41_51,                              &
                            EUOLD, SNOLD,  EPCU,  EPCUS,                &
                            EO,ET, EC,SC,  FC,FT, EBU,   ECU,           &
                 IUL,   EPPC,   EPPT,  EPEC,                            &
                 CRACK, INCF,   ICIC , SENS,    VS,                     &
                            IVIRGIN,ELIMIT,MM,IROT,ICR,LLLCRACK)        
! C     *                 DCC11,DCC12,DCC13,                            
! C     *                       DCC21,DCC22,DCC23,
! C     *                       DCC31,DCC32,DCC33,
! C     *                       DCR11,DCR12,DCR13,
! C     *                       DCR21,DCR22,DCR23,
! C     *                       DCR31,DCR32,DCR33)
! **********************************************************************
! * 引数�E�　DEU           �E�EI)吁E��チE��プ�E主ひずみ増�Eの解
! * 　　EU�E�SN       �E�EI,O)現STEPで計算する主ひずみ、主応力
! *     EEN,SEN     �E�EI)E点の主ひずみ、主応力
! *     EXC,SXC     �E�EI)X点の主ひずみ、主応力
! *     ECP,SCP     �E�EI)C点の主ひずみ、主応力
! *     ETP,STP     �E�EI)T点の主ひずみ、主応力
! *     EPC,SPC     �E�EI)P点の主ひずみ、主応力
! *     EJ,SJ       �E�EI)J点の主ひずみ、主応力
! *     ERC,SRC     �E�EI)R点(除荷曲線上！ETの間）�E再裁荷点)の主ひずみ、主応力
! *     EBC,SBC     �E�EI)B点(除荷曲線上！ETの間）�E再裁荷した時にPC線との交点)の主ひずみ、主応力
! *     EJJ         �E�EI)J点�E�剛性復活点�E��E接線剛性
! *     ERT�E�SRT   �E�EI)R'点の主ひずみ、主応力
! *     EPC1ST      �E�EI)初めて圧縮老E�E履歴でひび割れが入ったとき�E�E�点のひずみを記録�E�EU, LUモチE��で使ぁE��E
! *     EPT,SPT     �E�EI)応力が０になるとき�E(P'点)ひずみ�E�応力
! *     ECR,SCR     �E�EI)コンクリート�E引張強度時！E点�E�主ひずみ、応力
! *     ETMAX,STMAX �E�EI)履歴中経験した最大�E�E'点�E�主ひずみ、応力
! *     ECCR,  SCCR �E�EI)圧縮側から除荷し、�Eび割れた時�Eひずみ、応力
! *     E41_51,S41_51�E�EI)IUL=41OR51の除荷点ひずみ、応力
! *     EUOLD, SNOLD�E�EI)直前スチE��プ�E等価一軸ひずみ,応力
! *     EPCU        �E�EI)コンクリート�E圧壊が生じた後�E収斂点ひずみ
! *     EPCUS       �E�EI)コンクリート�E圧壊が生じた後�E収斂点応力
! *     EO          �E�EI)コンクリート要素の初期剛性
! *     ET          �E�EI/O)当応力点で接線剛性
! *     EC,SC       �E�EI)コンクリート�E三軸応力状態を老E�Eした最大主ひずみ、応力�E�破壊曲面との接点の値�E�E
! *     FC,FT       �E�EI)コンクリート要素の一軸圧縮、引張強度
! *     EBU         �E�EI)チE��ションスチE��フニング効果有効限界ひずみ
! *     ECU         �E�EI)�E�ECU�E�最大圧縮応力時�E等価一軸ひずみ
! *     IUL         �E�EI/O)コンクリート要素の吁E���E点の2つの等価一軸方向�E応力状態を表す指樁E
! *     EPPC        �E�EI)PC線�E剛性�E�EUL=6の剛性�E�E再載荷剛性)
! *     EPPT        �E�EI)除荷剛性�E�引張方向かな�E�？！E
! *     EPEC        �E�EI)EC間�E剛性
! *     ANGL        �E�EI)1スチE��プ前の主応力方向�E角度�E�第�E�！E,3�E�主方向角度�E�E
! *     ANGCR       �E�EO)第�E�！E,3�E��Eび割れ角度
! *     CRACK       �E�EI)ひび割れ状態�E持E��！E
! *     INCF        �E�EO)コンクリート要素の残差力による収斂計算指樁E
! *     ICIC        �E�EI)圧縮応力上�E曲線オプション　0�E�SAENZ式　1�E�FAFITIS-SHAH弁E
! *     SENS        �E�EI)繰返し用の剛性復活点Jを決めるための感度係数�E�現在は使ってぁE��ぁE��E
! *     VS          �E�EI)繰返し履歴用でIUL=4,5へ戻るため�E感度係数
! *     IVIRGIN     �E�EI)除荷が経験したかの持E��E
! *     ELIMIT      �E�EI)
! *     MM          �E�EI)コンクリート要素番号
! **********************************************************************
! **** STRESS - STRAIN RELATIONSHIP FOR CONCRETE                    ****
! **** ORIGINAL PROGRAM = SUBROUTINE SECONC(CYCLIC RULE)            ****
! **** THIS PROGRAM IS MODIFIED    < 1995.12.12 S.TORIZO >          ****
! ****                                                              ****
! **** コンクリート�E吁E���E点でのひずみの増�E量によって除荷また�E   ****
! **** 載荷判定をするプログラム        < 1996.08    J.ZAKKY >       ****
! **********************************************************************
! C
! C *** STRESS - STRAIN RELATIONSHIP FOR CONCRETE
! C     ORIGINAL PROGRAM = SUBROUTINE SECONC(CYCLIC RULE)
! C     THIS PROGRAM IS MODIFIED FOR MONOTONIC LOADING IN ORDER TO
! C     APPLY TO CALCULATE UNBLANCED STRESS < 1991.8.9 K.UCHIDA & A.AMEMIY
! C
    real (kreal) :: DEU,EU,SN,ERC,SRC,EXC,SXC,EBC,SBC,EJ,SJ,EJJ
    real (kreal) :: ERT,SRT,ECP,SCP,ETP,STP,EPC,SPC,EPC1ST,EPT,EEN,SEN
    real (kreal) :: ETMAX,STMAX,ECCR,SCCR,E41_51,S41_51,EUOLD,SNOLD,EPCU,EPCUS 
    real (kreal) :: EO,ET,EC,SC,FC,FT,EBU,ECU,EPPC,EPPT,EPEC,CRACK,SENS,VS
    real (kreal) :: ELIMIT,MODEL42_52,ECR,SCR
    
    real (kreal) :: DUMMY_SHIRAI,C_IC4,DUMMY_41,DUM41,DUM5,DUM51,DUMMY_51
    real (kreal) :: ec4,dum,dum4,ec5,es,eu4,eu5,saa,sbb,sc4,sc5,sjdammy
    
    integer :: IUL,INCF,IVIRGIN,MM,IROT,ICR,LLLCRACK
    integer :: inti,lx,ly,lz,m,ncall,nst,nit
    integer :: ICIC(8)
    real (kreal) :: CYCN(49)
! c IF(EBU.EQ.0.0) EBU=0.005D0
    IF(EBU.EQ.0.0) EBU=0.002D0
    MODEL42_52=1
    VS=1.0

      ECR=FT/EO
      SCR=FT
      IF(IUL.EQ. 1)  GOTO  1
      IF(IUL.EQ. 2)  GOTO  2
      IF(IUL.EQ. 3)  GOTO  3
      IF(IUL.EQ. 31)  GOTO 3
      IF(IUL.EQ. 4)  GOTO  4
      IF(IUL.EQ.41)  GOTO 41
      IF(IUL.EQ.42)  GOTO 42
      IF(IUL.EQ.43)  GOTO 43
      IF(IUL.EQ.46)  GOTO 46
      IF(IUL.EQ. 5)  GOTO  5
      IF(IUL.EQ.51)  GOTO 51
      IF(IUL.EQ.52)  GOTO 52
      IF(IUL.EQ.53)  GOTO 53
      IF(IUL.EQ.56)  GOTO 56
      IF(IUL.EQ. 6)  GOTO  6
      IF(IUL.EQ. 7)  GOTO  7
      IF(IUL.EQ. 8)  GOTO  8
      IF(IUL.EQ. 9)  GOTO  9
! ***************
! **** IUL=1 ****
! ***************
    1 CALL ROOT_CHECK(   1, 57)
      IF(EU)       1001,1002,1002
 1001 CALL ROOT_CHECK(1001, 57)
      IF(DEU)      1010,1010,1003
 1003 CALL ROOT_CHECK(1003, 57)
      IF(IVIRGIN)   901,1010,1004
 1010 CALL ROOT_CHECK(1010, 57)
      IF(EU.LT.EPCU) GOTO 9000
      IF(EU.LE.EC)   GOTO 8000
      IF(EC.LT.EU)   GOTO 1000
 1002 CALL ROOT_CHECK(1002, 57)
    IF(CRACK.EQ.4.0) THEN
        IF(EU.GE.EBU) GOTO 4090
        CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
        CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
        IF(DUMMY_SHIRAI.LT.DUMMY_41) THEN
            GOTO 4000
        ELSE
            GOTO 4100
        ENDIF
    ENDIF
      IF(ECR.LT.EU)  GOTO 4002
      IF(EU.LE.ECR)  GOTO 1040
! ****IUL=1からの除荷の準備  IUL=2を経由
! C 1004 CALL ROOT_CHECK(1004, 57)
! C
1004    CONTINUE
      IVIRGIN=1
      EEN=EUOLD ! EUCOLD(*,*,*,*,*,*)
      SEN=SNOLD ! PSIOLD(*,*,*,*,*,*)
      CALL ECTPJ( EEN, SEN, ECP, SCP, ETP, STP,         &
                 EPC, SPC,  EJ,  SJ,ECCR,SCCR,          &
                EPCU,EPPC,                              &
                  EO,  SC,  FT,  ECU,CRACK,ICIC(1),     &
                  EC,ETMAX,ECR,EPEC,                    &
                NST,NIT,M,LX,LY,LZ,intI,NCALL )
      ET=EPEC
      IF(EPPC.EQ.EO)THEN
        EXC=EEN
        SXC=SEN
      ELSE
        CALL CROSS1(EPPC,ECP,SCP,EXC,SXC,EO,SC,EPC,     &
                  EC,IUL,EPCU,EPCUS,ICIC(1))
      ENDIF
      IF(EPPC.EQ.EO)THEN
        IF(EU.LT.EPC) GOTO 6000                        !もし、現ひずみが残留ひずみより小さければIUL=3決定、これより下�E引張域、E
        IF(CRACK.EQ.0.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          IF(ECCR.LE.EU)  GOTO 5000
          IF(EU.LT.ECCR)  GOTO 3010
        ENDIF
        IF(CRACK.EQ.4.0)THEN
          IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
        ENDIF
        IF(CRACK.EQ.5.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
            GOTO 5000
          ELSE
            GOTO 5100
          ENDIF
        ENDIF
      ELSE
        IF(EU.LE.ETP)  GOTO 2000
       IF(EU.LT.EPC)  GOTO 3000
        IF(CRACK.EQ.0.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          IF(ECCR.LE.EU)  GOTO 5000
          IF(EU.LT.ECCR)  GOTO 3010
        ENDIF
        IF(CRACK.EQ.4.0)THEN
          IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
        ENDIF
        IF(CRACK.EQ.5.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
            GOTO 5000
          ELSE
            GOTO 5100
          ENDIF
        ENDIF
      ENDIF
! ****応力上�E曲線丁E
 1000 CALL ROOT_CHECK(1000, 57)
      IUL=1
      INCF=1
    CALL CompressionZone(ICIC(1),EO,ET,ECU,SC,EU,SN)
      IF(ET.LT.EO/100.0) ET=EO/100.0
! C     IF(SN.LT.SC) GOTO 903
      IF(IVIRGIN.EQ.0)THEN
        IF(SN.LT.SC/3.)IVIRGIN=1
      ENDIF
      RETURN
! ****引張弾性埁E
 1040 CALL ROOT_CHECK(1040, 57)
      IUL=1
      INCF=1
      ET=EO
      SN=EO*EU
      RETURN
! ***************
! **** IUL=2 ****
! ***************
    2 CALL ROOT_CHECK(   2, 57)
      IF(DEU)2001,2001,2002
 2001 CALL ROOT_CHECK(2001, 57)
    IF(EEN.GT.EC)THEN
        IF(EU.LT.EPCU) GOTO 9000
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EEN)  GOTO 1000
        IF(EEN.LT.EU)  GOTO 2000
    ELSE
        IF(EU.LT.EPCU) GOTO 9000
        IF(EU.LE.EEN)  GOTO 8000
        IF(EEN.LT.EU)  GOTO 2000        
    ENDIF
 2002 CALL ROOT_CHECK(2002, 57)
      IF(EU.LE.ETP)  GOTO 2000
      IF(EU.LE.EPC)  GOTO 3000
      IF(CRACK.EQ.0.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        IF(ECCR.LT.EU)  GOTO 5000
        IF(EU.LE.ECCR)  GOTO 3010
      ENDIF
      IF(CRACK.EQ.4.0)THEN
        IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
      ENDIF
      IF(CRACK.EQ.5.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
        CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
        IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
          GOTO 5000
        ELSE
          GOTO 5100
        ENDIF
      ENDIF
 2000 CALL ROOT_CHECK(2000, 57)
      IUL=2
      INCF=1
      ET=EPEC
      SN=SEN+(EU-EEN)*EPEC
      RETURN
! ***************
! **** IUL=3 ****
! ***************
    3 CALL ROOT_CHECK(   3, 57)
    
    IF(EPEC.EQ.EPPC) GOTO 6

      IF(EPC.LT.EUOLD) GOTO 3003                                               !直線部にジャンチE
      IF(DEU) 3001,3002,3002
 3001 CALL ROOT_CHECK(3001, 57)
      ERC=EUOLD
      SRC=SNOLD

      EBC=(EPEC*ERC-SRC-EPPC*EPC)/(EPEC-EPPC)
      SBC=EPPC*(EBC-EPC)

      IF(EU.LT.EPCU) GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)  GOTO 8000
      ELSE
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EXC)  GOTO 1000
      ENDIF
      IF(EU.LT.EBC)  GOTO 6000
      IF(EBC.LT.EU)  GOTO 7000
 3002 CALL ROOT_CHECK(3002, 57)
      IF(EU.LE.EPC)  GOTO 3000
      IF(CRACK.EQ.0.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        IF(ECCR.LT.EU)  GOTO 5000
        IF(EU.LE.ECCR)  GOTO 3010
      ENDIF
      IF(CRACK.EQ.4.0)THEN
        IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
      ENDIF
      IF(CRACK.EQ.5.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
        CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
        IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
          GOTO 5000
        ELSE
          GOTO 5100
        ENDIF
      ENDIF
 3003 CALL ROOT_CHECK(3003, 57)
      IF(DEU) 3004,3005,3005
 3004 CALL ROOT_CHECK(3004, 57)
      IF(EU.LT.EPCU) GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)  GOTO 8000
      ELSE
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EXC)  GOTO 1000
      ENDIF
      IF(EU.LT.EPC)  GOTO 6000
      IF(EPC.LE.EU)  GOTO 3010
 3005 CALL ROOT_CHECK(3005, 57)
      IF(EBU.LT.EU)   GOTO 5090
      IF(ECCR.LT.EU)  GOTO 5000
      IF(EU.LE.ECCR)  GOTO 3010
! ****IUL=3 曲線丁E
 3000 CALL ROOT_CHECK(3000, 57)
      IUL=3
      INCF=1
      IF(EPPC.EQ.EO)THEN
        ET=EO
        SN=(EU-EPC)*EO
      ELSE
        CALL CURVE3(EU,SN,EPC,ETP,STP,EPEC,ET )
      ENDIF
      RETURN
! ****IUL=3 直線丁E
 3010 CALL ROOT_CHECK(3010, 57)
      IUL=31
      INCF=1
      ET=EO
      SN=EO*(EU-EPC)
      RETURN
! ***********************************************
! **** IUL=4 (圧縮域を経験してぁE��いT/S領域)****
! ***********************************************
    4 CALL ROOT_CHECK(   4, 57)
      IF(DEU) 4001,4004,4004
 4001 CALL ROOT_CHECK(4001, 57)                                                !除荷発甁E
      IF(EBU.LT.EU) GOTO 4090                                                  !EBU:T/S限界ひずみ
      IF(EBU.LE.EUOLD)THEN
        ETMAX=EBU
        STMAX=0.0
      ELSE
        ETMAX=EUOLD
        STMAX=SNOLD
      ENDIF
! ***************************************
      IF(MODEL42_52.EQ.0)THEN                                                  !鳥山垁E
        EPPT=EO*ECR/ETMAX                                                      !出所不�E
      ENDIF
      IF(MODEL42_52.EQ.1.OR.        &                                              !直線モチE��
       MODEL42_52.EQ.4.OR.          &                                           !鳥山�E�E��線　混合モチE��
       MODEL42_52.EQ.5.OR.          &                                         !改良鳥山モチE��
       MODEL42_52.EQ.7)THEN                                                  !櫻井モチE��用
        EPPT=(STMAX-SJ)/(ETMAX-EJ)                                             !圧縮域での除荷を経験してぁE��ぁE�Eで、E��荷剛性は引張強度をもとに決める、E
      ENDIF
      IF(MODEL42_52.EQ.2)THEN                                                  !大乁E��論文
        EPPT=1.5*EO*ECR/ETMAX
      ENDIF
      IF(MODEL42_52.EQ.3)THEN                                                  !初期剛性モチE��
        EPPT=EO
      ENDIF
! C      IF(MODEL42_52.EQ.7)THEN                                                  !長沼モチE��
! C        EPPT=EO*ECR/ETMAX
! C      ENDIF
! ***************************************
! C      IF(MODEL42_52.EQ.7.AND.EPT.LT.5.0E-5)THEN
! C        EPPT=EO
! C      ENDIF
        EPT=ETMAX-STMAX/EPPT
      IF(EU.LT.EPCU) GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)  GOTO 8000
        IF(EU.LE.EJ)   GOTO 6000
        IF(EJ.LT.EU)   GOTO 4210
      ELSE
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EJ)THEN
          IF(EEN.EQ.0.0)THEN
            GOTO 1000
          ELSE
            IF(EU.LE.EXC) GOTO 1000
            IF(EXC.LT.EU) GOTO 6000
          ENDIF
        ELSE
          GOTO 4210
        ENDIF
      ENDIF
 4004 CALL ROOT_CHECK(4004, 57)
      IF(EBU.LE.EU) GOTO 4090
      IF(EU.LE.EBU) GOTO 4000
! ****TENSION STIFFENING 進衁E
! ****IUL=1>>4  ひび発甁ETENSION STIFFENING
 4002 CALL ROOT_CHECK(4002, 57)
! C      IF(ANGCR.LT.900.)GOTO 4003
! C      ANGCR=ANGL
 4003 CALL ROOT_CHECK(4003, 57)
      IVIRGIN=1
      ET=EO/100.
      IF(EEN.EQ.0.0)THEN
! *       SJ=-FT*SENS                                                            !IUL=2を未経験（感度解析部刁E��E
! C        SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)*15                                      !FROM KOBAYASHI(MODIFIED BY SAKURAI)
        SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)                                      !FROM TIDE
      CALL CompressionZone_REVERS(ICIC(1),FC,SJ,EO,EJJ,ECU,SC,EJ,SJDAMMY)
      ELSE
! *       SJ=SEN*0.1                                                             !IUL=2を経騁E
! C        SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)*15                                      !FROM KOBAYASHI(MODIFIED BY SAKURAI)
        SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)                                     !FROM TIDE
        EJ=SJ/EPPC+EPC
      ENDIF
! ****収斂回数�E�回目以陁E
 4005 CALL ROOT_CHECK(4005, 57)
      IF(EU.LT.EBU) GOTO 4000
      IF(EBU.LE.EU) GOTO 4090
! ****白井式曲線丁E
 4000 CALL ROOT_CHECK(4000, 57)
      IUL=4
      INCF=3
      CRACK=4.0
      ET=EO/100.
      CALL SHIRAI(EU,SN,ECR,SCR,EBU,ICIC(4),C_IC4)
    ICR = 1
      RETURN
! ****さらに允E
 4090 CALL ROOT_CHECK(4090, 57)
      IUL=4
      INCF=3
      CRACK=4.0
      ET=EO/100.
      SN=0.0
    ICR = 1
      RETURN
! ****************
! **** IUL=41 ****
! ****************
   41 CALL ROOT_CHECK(  41, 57)
      IF(DEU)4101,4102,4102
 4101 CALL ROOT_CHECK(4101, 57)
      IF(EU.LT.EPCU)    GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)   GOTO 8000
      ELSE
        IF(EU.LE.EC)    GOTO 8000
        IF(EU.LE.EXC)   GOTO 1000
      ENDIF
      IF(EU.LE.EJ)      GOTO 6000
      IF(EJ.LT.EU)      GOTO 4303
 4102 CALL ROOT_CHECK(4102, 57)
      IF(EBU.LT.EU)    GOTO 4090
      CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
      CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
      IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
        GOTO 4000
      ELSE
        GOTO 4100
      ENDIF
 4100 CALL ROOT_CHECK(4100, 57)
      IUL=41
      INCF=1
      CALL STRAIGHT41_51(EU,SN,ET,ETMAX,STMAX,EPC,VS)
      RETURN
! **************************************************************************************
! **** IUL=42 (圧縮の影響を受けてぁE��ぁE��ンション・スチE��フニング領域からの除荷曲緁E****
! **************************************************************************************
   42 CALL ROOT_CHECK(  42, 57)
      IF(DEU) 4201,4201,4202
 4201 CALL ROOT_CHECK(4201, 57)                                                !そ�Eまま進行、E
      IF(EU.LT.EPCU)   GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)   GOTO 8000
        IF(EU.LE.EJ)    GOTO 6000
        IF(EJ.LT.EU)    GOTO 4200
      ELSE
        IF(EU.LE.EC)    GOTO 8000
        IF(EU.LE.EJ)THEN
          IF(EEN.EQ.0.0)THEN
          IVIRGIN=0
! C       IVIRGIN=1
            GOTO 1000
          ELSE
            IF(EU.LE.EXC) GOTO 1000
            IF(EXC.LT.EU) GOTO 6000
          ENDIF
        ELSE
          GOTO 4200
        ENDIF
      ENDIF
 4202 CALL ROOT_CHECK(4202, 57)                                                !停滞（増�Eひずみ�E�E�E�もしくは、�E載荷、E
      IF(EBU.LT.EU)   GOTO 4090
      IF(ETMAX.LT.EU) GOTO 4000
      IF(EU.LE.ETMAX) GOTO 4610
 4210 CALL ROOT_CHECK(4210, 57)
      ETMAX=EUOLD
      STMAX=SNOLD
 4200 CALL ROOT_CHECK(4200, 57)
      IUL=42
      INCF=1
! ***********************************************************
      IF(MODEL42_52.EQ.0.OR.    &                                                 !鳥山垁E
        MODEL42_52.EQ.5)THEN                                                  !改良鳥山垁E
        IF(EPT.LT.EU)THEN
          EU4=ETMAX-EU
          EC4=ETMAX-EPT
          DUM4=1.+(EO*EC4/STMAX-2.)*EU4/EC4+(EU4/EC4)**2
          ET=EO*(1.-(EU4/EC4)**2)/DUM4**2
          SN=STMAX-EO*EU4/DUM4
        ELSE
          EU4=EU-EJ
          EC4=EPT-EJ
          SC4=-SJ
          DUM41=EU4/EC4
          DUM4=EO*EC4/SC4
          ET=EO*(1.-DUM41)**(DUM4-1.)
          SN=SC4*(1.-(1.-DUM41)**DUM4)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.1)THEN                                                  !直線モチE��
        ET=EPPT
        SN=ET*(EU-EJ)+SJ
      ENDIF
      IF(MODEL42_52.EQ.2)THEN                                                  !大乁E��論文
        IF(EPT.LT.EU)THEN
          ET=EPPT
          SN=ET*(EU-EPT)
        ELSE
          ET=-SJ*(EPT-EJ)
          SN=ET*(EU-EJ)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.3)THEN                                                  !初期剛性モチE��
        IF(EPT.LT.EU)THEN
          ET=EPPT
          SN=ET*(EU-EPT)
        ELSE
          ET=-SJ*(EPT-EJ)
          SN=ET*(EU-EJ)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.4)THEN                                                  !曲線！E��線　混合型
        IF(EPT.LT.EU)THEN
          EU4=ETMAX-EU
          EC4=ETMAX-EPT
          DUM4=1.+(EO*EC4/STMAX-2.)*EU4/EC4+(EU4/EC4)**2
          ET=EO*(1.-(EU4/EC4)**2)/DUM4**2
          SN=STMAX-EO*EU4/DUM4
        ELSE
          ET=EPPT
          SN=ET*(EU-EJ)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.7)THEN                                                  !長沼モチE��

       IF(SJ.GE.0.)THEN
         WRITE(6,*)'SJの応力が決まってません'
         STOP
       ENDIF

        IF(EPT.LT.EU)THEN
         ET=EPPT
         SN=ET*(EU-EJ)+SJ
        ELSE
          IF(EPT.LT.1.5*ECR)THEN
         ET=EPPT
         SN=ET*(EU-EJ)+SJ
          ELSE
           IF(EEN.LT.0.0) THEN
            CALL NN42_52( SJ, EJ, EPPC, EPT, EU, SN, ET, MM)
           ELSE
            CALL NN42_52( SJ, EJ,  EJJ, EPT, EU, SN, ET, MM)
           ENDIF
          ENDIF
        ENDIF

! *      WRITE(1114,*)'来てますよ',SJ,EJ

      ENDIF
      RETURN
! ****************
! **** IUL=43 ****
! ****************
   43 CALL ROOT_CHECK(  43, 57)
      IF(DEU)4301,4302,4302
 4301 CALL ROOT_CHECK(4301, 57)
      IF(EU.LT.EPCU)    GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)   GOTO 8000
      ELSE
        IF(EU.LE.EC)    GOTO 8000
        IF(EU.LE.EXC)   GOTO 1000
      ENDIF
      IF(EU.LE.EJ)      GOTO 6000
      IF(EJ.LT.EU)      GOTO 4300
 4302 CALL ROOT_CHECK(4302, 57)
      IF(EBU.LT.EU)     GOTO 4090
      CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
      CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
      IF(DUMMY_SHIRAI.LT.DUMMY_41)GOTO 4000
      IF(E41_51.LT.EU)  GOTO 4100
      IF(EU.LE.E41_51)  GOTO 4300
 4303 CALL ROOT_CHECK(4303, 57)
      E41_51=EUOLD
      S41_51=SNOLD
 4300 CALL ROOT_CHECK(4300, 57)
      IUL=43
      INCF=1
! *      IF(MODEL42_52.EQ.7)THEN
! *       IF(EU.GT.0.001)THEN
! *         GOTO4001
! *         CALL NN43_53( EO, SJ, EJ, EPPC, S41_51, E41_51, EU, SN, ET, MM)
! *       ENDIF
! *      ELSE
       ET=(S41_51-SJ)/(E41_51-EJ)
       SN=ET*(EU-EJ)+SJ
! *      ENDIF
      RETURN
! ****************
! **** IUL=46 ****
! ****************
   46 CALL ROOT_CHECK(  46, 57)
      IF(DEU)4601,4602,4602
 4601 CALL ROOT_CHECK(4601, 57)
      IF(EU.LT.EPCU)  GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC) GOTO 8000
        IF(EU.LE.EJ)  GOTO 6000
      ELSE
        IF(EU.LE.EC)    GOTO 8000
        IF(EU.LE.EJ)THEN
          IF(EEN.EQ.0.0)THEN
            GOTO 1000
          ELSE
            IF(EU.LE.EXC) GOTO 1000
            IF(EXC.LT.EU) GOTO 6000
          ENDIF
        ENDIF
      ENDIF
      IF(EU.LE.ERT)   GOTO 4200
      IF(ERT.LT.EU)   GOTO 4600
 4602 CALL ROOT_CHECK(4602, 57)
      IF(EBU.LT.EU)   GOTO 4090
      IF(ETMAX.LT.EU) GOTO 4000
      IF(EU.LE.ETMAX) GOTO 4600
 4610 CALL ROOT_CHECK(4610, 57)
      ERT=EUOLD
      SRT=SNOLD
 4600 CALL ROOT_CHECK(4600, 57)
      IUL=46
      INCF=1
    IF(MODEL42_52.EQ.0.OR.      &                                                  !鳥山垁E
        MODEL42_52.EQ.5)THEN                                                  !改良鳥山垁E
        IF(EPT.LT.EU)THEN
          EU4=ETMAX-EU
          EC4=ETMAX-EPT
          DUM4=1.+(EO*EC4/STMAX-2.)*EU4/EC4+(EU4/EC4)**2
          ET=EO*(1.-(EU4/EC4)**2)/DUM4**2
          SN=STMAX-EO*EU4/DUM4
        ELSE
          EU4=EU-EJ
          EC4=EPT-EJ
          SC4=-SJ
          DUM41=EU4/EC4
          DUM4=EO*EC4/SC4
          ET=EO*(1.-DUM41)**(DUM4-1.)
          SN=SC4*(1.-(1.-DUM41)**DUM4)+SJ
        ENDIF
      ELSE IF(MODEL42_52.EQ.7)THEN
        IF(EU.GE.EPT)THEN
         ET=(STMAX-SRT)/(ETMAX-ERT)
         SN=SRT+(EU-ERT)*ET
        ELSE
          IF(EPT.LT.1.5*ECR)THEN
         ET=(STMAX-SRT)/(ETMAX-ERT)
         SN=SRT+(EU-ERT)*ET
          ELSE
           IF(EEN.LT.0.0) THEN
            CALL NN42_52( SJ, EJ, EPPC, EPT, EU, SN, ET, MM)
           ELSE
            CALL NN42_52( SJ, EJ,  EJJ, EPT, EU, SN, ET, MM)
           ENDIF
          ENDIF
        ENDIF
      ELSE
       ET=(STMAX-SRT)/(ETMAX-ERT)
       SN=SRT+(EU-ERT)*ET
      ENDIF
      RETURN
! ****************
! **** IUL=5  ****
! ****************
    5 CALL ROOT_CHECK(   5, 57)
      IF(DEU) 5001,5003,5003
 5001 CALL ROOT_CHECK(5001, 57)
      IF(EBU.LE.EU)  GOTO 5090
      IF(EBU.LT.EUOLD)THEN
        ETMAX=EBU
        STMAX=0.0
      ELSE
        ETMAX=EUOLD
        STMAX=SNOLD
      ENDIF
! ***************************************
      IF(MODEL42_52.EQ.0)THEN                                                  !鳥山垁E
 5099   CALL ROOT_CHECK(5099, 57)
        EPPT=EO*ECCR/ETMAX                                                     !出所不�E
      ENDIF
      IF(MODEL42_52.EQ.1.OR.                &                                  !直線モチE��
        MODEL42_52.EQ.4.OR.                &                                  !鳥山�E�E��線　混合モチE��
        MODEL42_52.EQ.5.OR.                &                                  !改良鳥山モチE��
        MODEL42_52.EQ.7)THEN                                                  !櫻井モチE��用
 5098   CALL ROOT_CHECK(5098, 57)
        EPPT=(STMAX-SJ)/(ETMAX-EJ)
      ENDIF
      IF(MODEL42_52.EQ.2)THEN                                                  !大乁E��論文
 5097   CALL ROOT_CHECK(5097, 57)
        EPPT=1.5*EO*ECCR/ETMAX
      ENDIF
      IF(MODEL42_52.EQ.3)THEN                                                  !初期剛性モチE��
 5096   CALL ROOT_CHECK(5096, 57)
        EPPT=EO
      ENDIF
      IF(MODEL42_52.EQ.4)THEN                                                  !鳥山�E�E��線　混合モチE��
 5095   CALL ROOT_CHECK(5095, 57)
        EPPT=(STMAX-SJ)/(ETMAX-EJ)
      ENDIF
! C      IF(MODEL42_52.EQ.7)THEN                                                  !長沼モチE��
! C 5094   CALL ROOT_CHECK(5094, 57)
! C        EPPT=EO*ECCR/ETMAX
! C      ENDIF
! ***************************************
! C      IF(MODEL42_52.EQ.7.AND.EPT.LT.(EPT+5.0E-5))THEN
! C        EPPT=EO
! C      ENDIF
        EPT=ETMAX-STMAX/EPPT
      IF(EU.LT.EPCU) GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)  GOTO 8000
      ELSE
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EXC)  GOTO 1000
      ENDIF
      IF(EU.LE.EJ)   GOTO 6000
      IF(EJ.LT.EU)   GOTO 5210
 5003 CALL ROOT_CHECK(5003, 57)
      IF(EBU.LT.EU)  GOTO 5090
      IF(EU.LE.EBU)  GOTO 5000
! *****収斂�E�回目以陁E
 5004 CALL ROOT_CHECK(5003, 57)
      IF(EBU.LT.EU)  GOTO 5090
      IF(EU.LE.EBU)  GOTO 5002
 5000 CALL ROOT_CHECK(5000, 57)
! C      IF(ANGCR.LT.900.)THEN                                                    !ひずみ軟化領域
! C        EPC1ST=EPC                                                             !初めて圧縮老E�E履歴でひび割れが入ったとき�E�E�点のひずみを記録�E�EU, LUモチE��で使ぁE��E
! C        GOTO 5002
! C      ENDIF
! C      ANGCR=ANGL
 5002 CALL ROOT_CHECK(5002, 57)
      IUL=5
      INCF=3
      CRACK=5.0
      ET=EO/100.
      CALL SHIRAI(EU,SN,ECCR,SCCR,EBU,ICIC(4),C_IC4)
    ICR = 1
      RETURN
 5090 CALL ROOT_CHECK(5090, 57)
      IUL=5                                                                    !チE��ションスチE��フニングを趁E��た領域
      INCF=3
      CRACK=5.0
      ET=EO/100.
      SN=0.0
    ICR = 1
      RETURN
! ****************
! **** IUL=51 ****
! ****************
   51 CALL ROOT_CHECK(  51, 57)
      IF(DEU)5101,5102,5102
 5101 CALL ROOT_CHECK(5101, 57)
      IF(EU.LT.EPCU)  GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC) GOTO 8000
      ELSE
        IF(EU.LE.EC)  GOTO 8000
        IF(EU.LE.EXC) GOTO 1000
      ENDIF
      IF(EU.LE.EJ)    GOTO 6000
      IF(EJ.LT.EU)    GOTO 5303
 5102 CALL ROOT_CHECK(5102, 57)
      IF(EBU.LT.EU)   GOTO 5090
      CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
      CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
      IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
        GOTO 5000
      ELSE
        GOTO 5100
      ENDIF
 5100 CALL ROOT_CHECK(5100, 57)
      IUL=51
      INCF=1
      CALL STRAIGHT41_51(EU,SN,ET,ETMAX,STMAX,EPC,VS)
      RETURN
! ****************
! **** IUL=52 ****
! ****************
   52 CALL ROOT_CHECK(  52, 57)
      IF(DEU)5201,5202,5202
 5201 CALL ROOT_CHECK(5201, 57)
      IF(EU.LT.EPCU)  GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)  GOTO 8000
      ELSE
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EXC)  GOTO 1000
      ENDIF
      IF(EU.LE.EJ)    GOTO 6000
      IF(EJ.LT.EU)    GOTO 5200
 5202 CALL ROOT_CHECK(5202, 57)
      IF(EBU.LT.EU)   GOTO 5090
      IF(ETMAX.LT.EU) GOTO 5000
      IF(EU.LE.ETMAX) GOTO 5610
 5210 CALL ROOT_CHECK(5210, 57)
      ETMAX=EUOLD
      STMAX=SNOLD
 5200 CALL ROOT_CHECK(5200, 57)
      IUL=52
      INCF=1
! ***********************************************************
      IF(MODEL42_52.EQ.0.OR.    &                                              !鳥山垁E
        MODEL42_52.EQ.5)THEN                                                  !改良鳥山垁E
        IF(EPT.LT.EU)THEN
          EU5=ETMAX-EU
          EC5=ETMAX-EPT
          DUM5=1.+(EO*EC5/STMAX-2.)*EU5/EC5+(EU5/EC5)**2
          ET=EO*(1.-(EU5/EC5)**2)/DUM5**2
          SN=STMAX-EO*EU5/DUM5
        ELSE
          EU5=EU-EJ
          EC5=EPT-EJ
          SC5=-SJ
          DUM51=EU5/EC5
          DUM5=EO*EC5/SC5
          ET=EO*(1.-DUM51)**(DUM5-1.)
          SN=SC5*(1.-(1.-DUM51)**DUM5)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.1)THEN                                                  !直線モチE��
! C        EPPT=(STMAX-SJ)/(ETMAX-EJ)
        ET=EPPT
        SN=ET*(EU-EJ)+SJ
      ENDIF
      IF(MODEL42_52.EQ.2)THEN                                                  !大乁E��論文
        IF(EPT.LT.EU)THEN
          ET=EPPT
          SN=ET*(EU-EPT)
        ELSE
          ET=-SJ*(EPT-EJ)
          SN=ET*(EU-EJ)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.3)THEN                                                  !初期剛性モチE��
        IF(EPT.LT.EU)THEN
          ET=EPPT
          SN=ET*(EU-EPT)
        ELSE
          ET=-SJ*(EPT-EJ)
          SN=ET*(EU-EJ)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.4)THEN                                                  !鳥山�E�E��線　混合型
        IF(EPT.LT.EU)THEN
          EU5=ETMAX-EU
          EC5=ETMAX-EPT
          DUM5=1.+(EO*EC5/STMAX-2.)*EU5/EC5+(EU5/EC5)**2
          ET=EO*(1.-(EU5/EC5)**2)/DUM5**2
          SN=STMAX-EO*EU5/DUM5
        ELSE
          ET=EPPT
          SN=ET*(EU-EJ)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.7)THEN                                                  !長沼モチE��
        IF(EPT.LT.EU)THEN
         ET=EPPT
         SN=ET*(EU-EJ)+SJ
        ELSE
          IF(EPT.LT.1.0E-3)THEN
         ET=EPPT
         SN=ET*(EU-EJ)+SJ
          ELSE
           CALL NN42_52( SJ, EJ, EPPC, EPT, EU, SN, ET, MM)
         ENDIF
        ENDIF
      ENDIF
      RETURN
! ****************
! **** IUL=53 ****
! ****************
   53 CALL ROOT_CHECK(  53, 57)
      IF(DEU)5301,5302,5302
 5301 CALL ROOT_CHECK(5301, 57)
      IF(EU.LT.EPCU)    GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)   GOTO 8000
      ELSE
        IF(EU.LE.EC)    GOTO 8000
        IF(EU.LE.EXC)   GOTO 1000
      ENDIF
      IF(EU.LE.EJ)      GOTO 6000
      IF(EJ.LT.EU)      GOTO 5300
 5302 CALL ROOT_CHECK(5302, 57)
      IF(EBU.LT.EU)     GOTO 5090
      CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
      CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
      IF(DUMMY_SHIRAI.LT.DUMMY_51)GOTO 5000
      IF(E41_51.LT.EU)  GOTO 5100
      IF(EU.LE.E41_51)  GOTO 5300
 5303 CALL ROOT_CHECK(5303, 57)
      E41_51=EUOLD
      S41_51=SNOLD
 5300 CALL ROOT_CHECK(5300, 57)
      IUL=53
      INCF=1
! C      IF(MODEL42_52.EQ.7)THEN
! C       IF(EU.GT.0.001)THEN
! C        CALL NN43_53( EO, SJ, EJ, EPPC, S41_51, E41_51, EU, SN, ET, MM)
! C       ENDIF
! C      ELSE
       ET=(S41_51-SJ)/(E41_51-EJ)
       SN=ET*(EU-EJ)+SJ
! C      ENDIF
      RETURN
! ****************
! **** IUL=56 ****
! ****************
   56 CALL ROOT_CHECK(  56, 57)
      IF(DEU)5601,5602,5602
 5601 CALL ROOT_CHECK(5601, 57)
      IF(EU.LT.EPCU)  GOTO 9000
      IF(EXC.LT.EU)THEN
        IF(EU.LE.EXC)   GOTO 8000
      ELSE
        IF(EU.LE.EC)    GOTO 8000
        IF(EU.LE.EXC)   GOTO 1000
      ENDIF
      IF(EU.LE.EJ)    GOTO 6000
      IF(EU.LE.ERT)   GOTO 5200
      IF(ERT.LT.EU)   GOTO 5600
 5602 CALL ROOT_CHECK(5602, 57)
      IF(EBU.LT.EU)   GOTO 5090
      IF(ETMAX.LT.EU) GOTO 5000
      IF(EU.LE.ETMAX) GOTO 5600
 5610 CALL ROOT_CHECK(5610, 57)
      ERT=EUOLD
      SRT=SNOLD
 5600 CALL ROOT_CHECK(5600, 57)
      IUL=56
      INCF=1
      IF(MODEL42_52.EQ.7)THEN
        IF(EU.GE.EPT)THEN
         ET=(STMAX-SRT)/(ETMAX-ERT)
         SN=SRT+(EU-ERT)*ET
        ELSE
          IF(EPT.LT.1.0E-3)THEN
         ET=(STMAX-SRT)/(ETMAX-ERT)
         SN=SRT+(EU-ERT)*ET
          ELSE
            CALL NN42_52( SJ, EJ, EPPC, EPT, EU, SN, ET, MM)
         ENDIF
        ENDIF
      ELSE
       ET=(STMAX-SRT)/(ETMAX-ERT)
       SN=SRT+(EU-ERT)*ET
      ENDIF
      RETURN
! ***************
! **** IUL=6 ****
! ***************
    6 CALL ROOT_CHECK(   6, 57)
      IF(DEU)6001,6002,6002
 6001 CALL ROOT_CHECK(6001, 57)
      IF(EU.LT.EPCU) GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)  GOTO 8000
      ELSE
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EXC)  GOTO 1000
      ENDIF
      IF(EXC.LT.EU) GOTO 6000
 6002 CALL ROOT_CHECK(6002, 57)
      IF(EUOLD.LT.ECP) GOTO 6004
      IF(EPPC.EQ.EO)THEN
        IF(EU.LT.EPC) GOTO 6000
        IF(EPC.LE.EU) GOTO 6003
      ELSE
        EBC=EUOLD
        SBC=SNOLD
        SBB=SBC+(EU-EBC)*EPEC
        CALL CURVE3(EU,SAA,EPC,ETP,STP,EPEC,ET)
        IF(SBB.LT.SAA) GOTO 7000
      ENDIF
      IF(EU.LE.EPC) GOTO 3000
 6003 CALL ROOT_CHECK(6003, 57)
      IF(CRACK.EQ.0.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        IF(ECCR.LT.EU)  GOTO 5000
        IF(EU.LE.ECCR)  GOTO 3010
      ENDIF
      IF(CRACK.EQ.4.0)THEN
        IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
      ENDIF
      IF(CRACK.EQ.5.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
        CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
        IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
          GOTO 5000
        ELSE
          GOTO 5100
        ENDIF
      ENDIF
 6004 CALL ROOT_CHECK(6004, 57)
      CALL CROSS1(EO,EUOLD,SNOLD,EEN,SEN,EO,SC,-SNOLD/EO+EUOLD,     &
                EC,IUL,EPCU,EPCUS,ICIC(1))
      CALL ECTPJ( EEN, SEN, ECP, SCP, ETP, STP,             &
                EPC, SPC,  EJ,  SJ,ECCR,SCCR,               &
               EPCU,EPPC,                                   &
                 EO,  SC,  FT,  ECU,CRACK,ICIC(1),          &
                 EC,ETMAX,ECR,EPEC,                         &
                NST,NIT,M,LX,LY,LZ,intI,NCALL )
      CALL CROSS1(EPPC,ECP,SCP,EXC,SXC,EO,SC,EPC,           &
                EC,IUL,EPCU,EPCUS,ICIC(1))
      IF(EPPC.EQ.EO)THEN
        IF(EU.LT.EPC) GOTO 6000
        IF(CRACK.EQ.0.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          IF(ECCR.LE.EU)  GOTO 5000
          IF(EU.LT.ECCR)  GOTO 3010
        ENDIF
        IF(CRACK.EQ.4.0)THEN
          IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
        ENDIF
        IF(CRACK.EQ.5.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
            GOTO 5000
          ELSE
            GOTO 5100
          ENDIF
        ENDIF
      ELSE
        IF(EU.LE.ETP)  GOTO 2000
        IF(EU.LE.EPC)  GOTO 3000
        IF(CRACK.EQ.0.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          IF(ECCR.LE.EU)  GOTO 5000
          IF(EU.LT.ECCR)  GOTO 3010
        ENDIF
        IF(CRACK.EQ.4.0)THEN
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
        ENDIF
        IF(CRACK.EQ.5.0)THEN
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
            GOTO 5000
          ELSE
            GOTO 5100
          ENDIF
        ENDIF
      ENDIF
 6000 CALL ROOT_CHECK(6000, 57)
      IUL=6
      INCF=1
      ET=EPPC
      SN=EPPC*(EU-EPC)
      RETURN
! ***************
! **** IUL=7 ****
! ***************
    7 CALL ROOT_CHECK(   7, 57)
      IF(DEU)7001,7002,7002
 7001 CALL ROOT_CHECK(7001, 57)
      IF(EU.LT.EPCU) GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)  GOTO 8000
      ELSE
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EXC)  GOTO 1000
      ENDIF
      IF(EU.LE.EBC)  GOTO 6000
      IF(EBC.LT.EU)  GOTO 7000
 7002 CALL ROOT_CHECK(7002, 57)
      SBB=SBC+(EU-EBC)*EPEC
      CALL CURVE3(EU,SAA,EPC,ETP,STP,EPEC,ET )
      IF(SBB.LT.SAA) GOTO 7000
      IF(EU.LE.EPC)  GOTO 3000
      IF(CRACK.EQ.0.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        IF(ECCR.LT.EU)  GOTO 5000
        IF(EU.LE.ECCR)  GOTO 3010
      ENDIF
      IF(CRACK.EQ.4.0)THEN
        IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
      ENDIF
      IF(CRACK.EQ.5.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
        CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
        IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
          GOTO 5000
        ELSE
          GOTO 5100
        ENDIF
      ENDIF
 7000 CALL ROOT_CHECK(7000, 57)
      IUL=7
      INCF=1
      ET=EPEC
      SN=SBC+(EU-EBC)*EPEC
      RETURN
! ***************
! **** IUL=8 ****
! ***************
    8 CALL ROOT_CHECK(   8, 57)
      IF(DEU)8001,8002,8002
 8001 CALL ROOT_CHECK(8001, 57)                                                !除荷が起きてぁE��ぁE��E
      IF(EU.LT.EPCU) GOTO 9000
      IF(EPCU.LE.EU) GOTO 8000
 8002 CALL ROOT_CHECK(8002, 57)                                                !除荷が発生した、E
      EEN=EUOLD
      SEN=SNOLD
      ET=EO
      CALL ECTPJ( EEN, SEN, ECP, SCP, ETP, STP,         &
                 EPC, SPC,  EJ,  SJ,ECCR,SCCR,          &
                EPCU,EPPC,                              &
                  EO,  SC,  FT,  ECU,CRACK,ICIC(1),     &
                  EC,ETMAX,ECR,EPEC,                    &
                NST,NIT,M,LX,LY,LZ,intI,NCALL )
      IF(EPPC.EQ.EO)THEN
        EXC=EEN
        SXC=SEN
      ELSE
        CALL CROSS1(EPPC,ECP,SCP,EXC,SXC,EO,SC,EPC,     &
                 EC,IUL,EPCU,EPCUS,ICIC(1))
      ENDIF
      IF(EPPC.EQ.EO)THEN
        IF(EU.LT.EPC) GOTO 6000                                                !圧縮域�Eまま�E�EUL=6�E�E
        IF(CRACK.EQ.0.0)THEN                                                   !引張域へ突�E
          IF(EBU.LT.EU)   GOTO 5090
          IF(ECCR.LE.EU)  GOTO 5000
          IF(EU.LT.ECCR)  GOTO 3010
        ENDIF
        IF(CRACK.EQ.4.0)THEN
          IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
        ENDIF
        IF(CRACK.EQ.5.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
            GOTO 5000
          ELSE
            GOTO 5100
          ENDIF
        ENDIF
      ELSE
        IF(EU.LE.ETP)  GOTO 2000
        IF(EU.LT.EPC)  GOTO 3000
        IF(CRACK.EQ.0.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          IF(ECCR.LE.EU)  GOTO 5000
          IF(EU.LT.ECCR)  GOTO 3010
        ENDIF
        IF(CRACK.EQ.4.0)THEN
          IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
        ENDIF
        IF(CRACK.EQ.5.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
            GOTO 5000
          ELSE
            GOTO 5100
          ENDIF
        ENDIF
      ENDIF
 8000 CALL ROOT_CHECK(8000, 57)
      IUL=8
      INCF=3
      ES=(EPCUS-FC)/(EPCU-EC)
      SN=SC+ES*(EU-EC)
      ET=EO/100.
      RETURN
! ***************
! **** IUL=9 ****
! ***************
    9 CALL ROOT_CHECK(   9, 57)
      IF(DEU)9001,9002,9002
 9001 CALL ROOT_CHECK(9001, 57)
      IF(ELIMIT.EQ.0.0) GOTO 9000
      IF(EU.LE.ELIMIT)  GOTO 9000
      IF(ELIMIT.LT.EU)  GOTO 9002
 9000 CALL ROOT_CHECK(9000, 57)
      IUL=9
      INCF=3
      ELIMIT=EU
      SN=EPCUS
      ET=EO/100.
      RETURN
 9002 CALL ROOT_CHECK(9002, 57)
      SN=0.0001
      ET=EO/100.
      RETURN
! *******************************************************************
! * ERROR MESSAGES                                                  *
! *******************************************************************
  901 WRITE(6,902)
  902 FORMAT(' IVIRGIN ERROR')
      STOP
  903 WRITE(6,904)
  904 FORMAT(' ***ERROR IN SECON2 (IUL=',I2,')***')

      STOP
      END
    
    
    subroutine ROOT_CHECK(A,B)
    integer :: A,B

    end
    
    !*******************************************************************
      SUBROUTINE SHIRAI(EU,SN,ECR,SCR,EBU,ICIC4,C_IC4)
! *******************************************************************
! ****  ADDED BY JOE ZAKKY ON 1996.8.20 FOR BECOMING MORE EASILY ****
! ****  TO ADJUST, CHECK OR CHANGE INTO OTHER EQUATIONS.         ****
! *******************************************************************
    real (kreal) :: EU,SN,ECR,SCR,EBU,C_IC4,EEX1,SSX1,EU_CR
    integer :: ICIC4
! *$real (kreal) ::INCLUDE"D:\ANALYSIS\ORIGINAL\TS\H\FILE.H"
! *$INCLUDE'/USER2/ZAKKY/Z9/H/FILE.H'

 !     CALL ROOT_CHECK(0000, 76)

    IF( ICIC4 == 0 ) GOTO 1000
    IF( ICIC4 == 3 ) GOTO 3000
    IF( ICIC4 == 4 ) GOTO 4000

1000    EEX1= ( EU - ECR ) / ( EBU - ECR )
    SSX1=1.0-2.748*EEX1+2.654*EEX1**2-0.906*EEX1**3
    SN= SSX1 * SCR
    GOTO 9999

3000    EU_CR=EU-ECR
    IF( EU_CR <= 0. ) THEN
        SN=0.0
        GOTO 9999
    ENDIF
    SN=(0.0002/EU_CR)**C_IC4*SCR
    !SN=(0.0002/EU)**C_IC4*SCR          !Edited by shan

    IF(SN < 0.) SN = 0.
    IF(SN > SCR ) SN = SCR
    GOTO 9999

4000    EU_CR=EU-ECR
    SN = C_IC4 * EU_CR + SCR
    GOTO 9999

9999    RETURN
      END
    
    !*******************************************************************
      SUBROUTINE STRAIGHT41_51(EU,SN,ET,ETMAX,STMAX,EPC,VS)
! *******************************************************************
! *入力頁E��  EU   : 等価一軸ひずみ                                  *
! *入力頁E��  EU   : 等価一軸ひずみ                                  *
! *          ETMAX: チE��ションスチE��フニング域での除荷開始点ひずみ  *
! *          STMAX: チE��ションスチE��フニング域での除荷開始点応力    *
! *          EPC  : 圧縮域残留ひずみ                                *
! *          VS   : 感度係数                                        *
! *                                                                 *
! *出力頁E��  SN   :履歴に対応する応力                               *
! ! *          ET   :最新接線剛性                                     *
! ! *******************************************************************
! *$INCLUDE"D:\ANALYSIS\ORIGINAL\TS\H\FILE.H"
! *$INCLUDE'/USER2/ZAKKY/Z9/H/FILE.H'
        real (kreal) :: EU,SN,ET,ETMAX,STMAX,EPC,VS
      ! CALL ROOT_CHECK(0000, 78)

      ET=VS*STMAX/(ETMAX-EPC)
      SN=ET*(EU-EPC)

      RETURN
      END

    !***********************************************************************
      SUBROUTINE CROSS1 (EPP,ECP,SCP,EX,SX,EO,SC,EPC,           &
                       EC,IUL,EPCU,EPCUS,ICIC)
! ***********************************************************************
! * 機�E         �E�EX.SXを決めるSUB
! * 引数�E�EEPP       :PC線�E剛性�E�EUL=6の剛性�E�！EPPC�E�E
! *         ECP     :C点の主ひずみ、主応力
! *         SCP     :C点の主ひずみ、主応力    
! *         EX      :X点の主ひずみ、主応力    
! *         SX      :X点の主ひずみ、主応力    
! *         EO      :コンクリート要素の初期剛性
! *         SC      :コンクリート�E三軸応力状態を老E�Eした最大主ひずみ、応力�E�破壊曲面との接点の値�E�E
! *         EPC     :P点の主ひずみ、主応力
! *         EC      :コンクリート�E三軸応力状態を老E�Eした最大主ひずみ、応力�E�破壊曲面との接点の値�E�E
! *         IUL     :コンクリート要素の吁E���E点の2つの等価一軸方向�E応力状態を表す指樁E
! *         EPCU    :コンクリート�E圧壊が生じた後�E収斂点ひずみ
! *         EPCUS   :コンクリート�E圧壊が生じた後�E収斂点応力
! *         ICIC    :圧縮応力上�E曲線オプション　0�E�SAENZ式　1�E�FAFITIS-SHAH弁E
! ***********************************************************************

        real (kreal) :: EPP,ECP,SCP,EX,SX,EO,SC,EPC,EC,EPCU,EPCUS
        integer :: ICIC,IUL,I
        
        real (kreal) :: EENOLD,SENOLD,EY,R,O,B,P,ET



! *$INCLUDE"D:\ANALYSIS\ORIGINAL\TS\H\FILE.H"
! *$INCLUDE"D:\ANALYSIS\ORIGINAL\TS\H\URAWAZA.H"

      ! CALL ROOT_CHECK(0000, 14)

    EENOLD=EX
    SENOLD=SX

      EY=(SC-SCP)/EPP+ECP
      IF(EY-EC) 10,10,20
! ****X点がIUL=8丁E
   10 R=(SC-EPCUS)/(EC-EPCU)
      EX=(SC-R*EC+EPP*EPC)/(EPP-R)
    IF(EENOLD.LE.EX) EX=EENOLD
      SX=R*(EX-EC)+SC
     ! CALL ROOT_CHECK(9999, 14)
      RETURN
! ****X点が応力上�E曲線丁E
   20 R=-EC/100.
   30 O=0.9*EC
      DO 40 I=1,50
    CALL CompressionZone(ICIC,EO,ET,EC,SC,O,B)
         P=-(B-O*ET+EPP*EPC)/(ET-EPP)
         IF(ABS(P-O).LT.R)  GOTO 50
         O=P
   40 CONTINUE
      IF(EPP.EQ.EO) GOTO 60
      EPP=EO
      GOTO 30
   50 EX=P
    CALL CompressionZone(ICIC,EO,ET,EC,SC,EX,SX)
      RETURN
! *** CROSS POINT NOT DEFINED
   60 P=EX-SX/EO
      GOTO 50
! C
  ! 101 FORMAT(' SUBROUTINE CROSS でトラブル発甁E SAENZ '/
     ! *       ' STAGE=',I6,' ITR=',I2,' ELEMENT#',I5,
     ! *       ' INT#',I1,' DIR#',I1,' IUL=',I3,
     ! *       ' EPS=',1PE15.6, ' SIG=',1PE15.6)
  ! 102 FORMAT(' SUBROUTINE CROSS でトラブル発甁E FAFITIS-SHAH '/
     ! *       ' STAGE=',I6,' ITR=',I2,' ELEMENT#',I5,
     ! *       ' INT#',I1,' DIR#',I1,' IUL=',I3,
     ! *       ' EPS=',1PE15.6, ' SIG=',1PE15.6)

      END
    
    !*******************************************************************
      SUBROUTINE CURVE3( EU,SN,EPC,ETP,STP,EPEC,ET )
! *******************************************************************
! ****  FOR CALCULATING THE VALUE OF POINT ON IUL=3 ON EU        ****
! *******************************************************************
    real (kreal) ::  EU,SN,EPC,ETP,STP,EPEC,ET
    real (kreal) :: NAA,NAB,NAC
! *$INCLUDE"D:\ANALYSIS\ORIGINAL\TS\H\FILE.H"
! *$INCLUDE'/USER2/ZAKKY/Z9/H/FILE.H'
    
      ! CALL ROOT_CHECK(0000, 15)

! *長沼さんの論文を参老E��以下�E式を使用することにしました。！E櫻亁E2001 �E�E

      NAA=(EPEC*(ETP-EPC)-STP)/(ETP-EPC)**2
      NAB=EPEC-2*NAA*ETP
      NAC=STP-NAA*(ETP**2)-NAB*ETP

      ET=2*NAA*EU+NAB
      SN=NAA*(EU**2)+NAB*EU+NAC

! *野崎さん�E方程式　こ�E式だとひずみが大きくなるにつれて曲線が正側に食い込んでしまぁE��E
! *EXCLEで確認済み、E
! *ちなみに、E��崎さん（小林さん�E��E論文にある方程式と下�E方程式�E違うので注意を要する、E
! *     IF(CRACK.EQ.0.0) ETW=EPPC/2.0
! *     IF(CRACK.NE.0.0) ETW=STMAX/(ETMAX-EPC)
! *     EW=EPC-ETP
! *     SSC=-STP
! *     EEW=SSC*EW/(SSC-EW*(ETW*EO)**0.5)
! *     SSW=EO*EEW/((EEW*EO/SSC-EEW/EW-EW/EEW)+2.0)
! *     EE=EU-ETP
! *     DUMW=1.0+(EO*EEW/SSW-2.0)*EE/EEW+(EE/EEW)**2
! *     ET=EO*(1.0-(EE/EEW)**2)/DUMW**2
! *     SN=EO*EE/DUMW+STP

      RETURN
      END
    
    !*******************************************************************
      SUBROUTINE NN42_52( SJ, EJ, EJJ, EPT, EU, SN, ET, MM)
! *******************************************************************
! *入力頁E��  EU   : 等価一軸ひずみ                                  *
! *          EJ   : J点�E�剛性復活点�E��Eひずみ                       *
! *          SJ   : J点�E�剛性復活点�E��E応力                         *
! *          EJJ  : J点�E�剛性復活点�E��E接線剛性                     *
! *          EPT  : 応力ぁEのとき�Eひずみ                           *
! *                                                                 *
! *出力頁E��  SN   :履歴に対応する応力                               *
! *          ET   :最新接線剛性                                     *
! *******************************************************************

! *$INCLUDE"D:\ANALYSIS\ORIGINAL\TS\H\FILE.H"
    real (kreal) :: SJ,EJ,EJJ,EPT,EU,SN,ET
    integer :: MM,I,J
    real (kreal) ::SE,SS,SS1,SS2,SS3,SEMAX,SEMIN,EE
    real (kreal) :: ANAGA,BNAGA,CNAGA
    
    
      SE=0.0
      SS=0.0
       IF(EJJ.EQ.0.0)THEN
! C         WRITE(NFOUT,401)
         STOP
       ENDIF
      SE=SJ/EJJ
      SS1=-EJ
      SS2=1.0                                                                   !収斂計算�Eための初期値の設宁E
      SS3=1.0                                                                   !収斂計算�Eための初期値の設宁E

      SEMAX=SE+1.0E-8
      SEMIN=SE-1.0E-8

      DO 200 I=1,10
       DO 100 J=1,10
           EE=1.0/10.**I
           SS2=SS2-EE
            IF(SS2.LE.SS1) THEN
              SS2=SS3
              GOTO 200
            ENDIF

            IF((EPT+SS2).EQ.0.0)THEN
! C             WRITE(NFOUT,402)
             STOP
            ENDIF

           SS=(EJ+SS2)*LOG((EJ+SS2)/(EPT+SS2))
            IF(SS.LE.SEMAX.AND.SS.GE.SEMIN) GOTO 300                            !SSが篁E��冁E��入っぁE
            IF(SS.GT.SEMAX) THEN
             SS2=SS3
             GOTO 200
             ENDIF

           SS3=SS2
  100  CONTINUE
  200 CONTINUE
  300 CONTINUE
         IF(SS.GT.SEMAX.OR.SS.LT.SEMIN) THEN
! C           WRITE(NFOUT,400)
! C           WRITE(NFOUT,*)MM,SS,SEMAX,SEMIN
           STOP
         ENDIF

       ANAGA=SS2                                                                !係数Aの決宁E
       BNAGA=-1.0*LOG(EPT+ANAGA)                                                !係数Bの決宁E
       CNAGA=EJJ*(EJ+ANAGA)                                                     !係数Cの決宁E

       SN=CNAGA*(LOG(EU+ANAGA)+BNAGA)
       ET=CNAGA/(EU+ANAGA)

  ! 400 FORMAT( 'SS42_52が収束しませんでした'、E )
  ! 401 FORMAT( 'NN42_52でEJJが０でぁE、E )
  ! 402 FORMAT( 'NN42_52で数学皁E��ラーでぁE、E )
      RETURN
      END
    
        SUBROUTINE CONCYC_hj( FC , FT , E0 , E0_100 , IUL , S , DS ,    &
                        S_OLD , T_OLD ,                                 &
                        O_S , O_T , T_S , T_T ,                         &
                        E_S , E_T , X_S , X_T ,                         &
                        F_S , F_T , Y_S , Y_T ,                         &
                        U_S , U_T , C_S , C_T ,                         &
                        V_S , V_T , W_S , W_T ,                         &
                        R1_S , R1_T , R2_S , R2_T ,                     &
                        P_S , H_S , G_S , CR_S ,                        &
                        ID_CRACK , ICIC , C_IC4 ,                       &
                        T , E , INCF,                                   &
                        NST , NIT , M ,LX,LY,LZ,intI )

    
    
    real (kreal) :: FC,FT,E0,E0_100,S,DS,S_OLD,T_OLD,O_S,O_T,T_S,T_T
    real (kreal) :: E_S,E_T,X_S,X_T,F_S,F_T,Y_S,Y_T,U_S,U_T,C_S,C_T
    real (kreal) :: V_S,V_T,W_S,W_T,R1_S,R1_T,R2_S,R2_T,P_S,H_S,G_S,CR_S 
    real (kreal) :: C_IC4,T,E ,Scr,Xdd,ZERO
    
    integer :: ID_CRACK,INCF,NST,NIT,M,LX,LY,LZ,intI,I,IUL
    
    integer :: ICIC(8), NN(999)

    DO I=1,342 ; NN(I)=0 ; ENDDO
    ZERO= 0.0
    IF(IUL /= 1 ) GOTO 1
    IF(S <= 0.) THEN ;  IUL=11 ; ELSE ; IUL=31 ; ENDIF

1   IF( IUL == 11 ) GOTO 11
    IF( IUL == 31 ) GOTO 31

    IF( IUL == 12 ) GOTO 12
    IF( IUL == 32 ) GOTO 32

    IF( IUL == 13 ) GOTO 13

    IF( IUL == 141 ) GOTO 141
    IF( IUL == 142 ) GOTO 142
    IF( IUL == 20 ) GOTO 20

    IF( IUL == 341 ) GOTO 341
    IF( IUL == 342 ) GOTO 342
    IF( IUL == 40 ) GOTO 40

    IF( IUL == 24 ) GOTO 24
    IF( IUL == 42 ) GOTO 42

! *****************
! **** IUL=11  ****
! *****************
11  IUL= 11
    CALL NN_C( NST,NIT,M,LX,LY,LZ,intI,IUL,NN )

    IF( S < U_S ) GOTO 13
    IF( S < O_S ) GOTO 12

    IF( ID_CRACK == 0 ) THEN
        IF( S > T_S ) GOTO 32
        IF( S > P_S ) GOTO 31
    ELSE
        IF( S > F_S ) GOTO 32
        IF( S > P_S ) GOTO 20
    ENDIF


    IF( DS > 0. .AND. T_OLD < O_T / 3. )    &
     CALL EXPT(S_OLD , T_OLD , E0 , FT , O_S , ID_CRACK, & ! 収斂後�E計算追加
                  E_S , E_T , X_S , X_T , V_S , V_T , P_S , T_S , T_T)
    IF( DS > 0. .AND. T_OLD < O_T / 3. ) GOTO 141

     INCF= 1
    CALL  CompressZone( ICIC(1) , E0 , O_S , O_T , S , T , E )
    GOTO 9999
! *****************
! **** IUL=31  ****
! *****************
31  IUL= 31
    CALL NN_C( NST,NIT,M,LX,LY,LZ,intI,IUL,NN )

    IF( S < E_S ) THEN
        IF( S < U_S ) GOTO 13
        IF( S < O_S ) GOTO 12
                      GOTO 11
    ENDIF
    IF( S < P_S ) THEN
        IF( P_S ==0. ) THEN
            GOTO 11
        ELSE
            GOTO 342
        ENDIF
    ENDIF

    IF( S > T_S ) GOTO 32
! C IF( ID_CRACK == 1 ) GOTO 20 !老E��E

    INCF= 1
    E= E0
    Scr= S - P_S
    T= E * Scr
    GOTO 9999
! *****************
! **** IUL=12  ****
! *****************
12  IUL= 12
    CALL NN_C( NST,NIT,M,LX,LY,LZ,intI,IUL,NN )

    IF( S < U_S ) GOTO 13

    IF( ID_CRACK == 0 ) THEN
        IF( S > T_S ) GOTO 32
        IF( S > P_S ) GOTO 31
    ELSE
        IF( S > F_S ) GOTO 32
        IF( S > P_S ) GOTO 20
    ENDIF
    IF( S > O_S ) GOTO 11

    IF( DS > 0. ) &
     CALL EXPT(  S_OLD , T_OLD , E0 , FT , O_S , ID_CRACK , & ! 収斂後�E計算追加
              E_S , E_T , X_S , X_T , V_S , V_T , P_S , T_S , T_T)
    IF( DS > 0. ) GOTO 141

    INCF= 3
    E= 0.0
    !if( O_S == U_S ) then ; write(*,*) '☁EY_Point:',IUL ; endif
    CALL XY_Point( O_S , O_T , U_S , U_T , S , T , Xdd , 2)
    GOTO 9999
! *****************
! **** IUL=32  ****
! *****************
32  IUL= 32
    CALL NN_C( NST,NIT,M,LX,LY,LZ,intI,IUL,NN )

    IF( S < E_S ) THEN
        IF( S < U_S ) GOTO 13
        IF( S < O_S ) GOTO 12
                      GOTO 11       
    ENDIF

    IF( S < P_S ) THEN
        IF( P_S ==0. ) THEN
            GOTO 11
        ELSE
            GOTO 342
        ENDIF
    ENDIF

    IF( ID_CRACK == 0 ) THEN
        IF( S < T_S ) GOTO 31
    ELSE
        IF( S < F_S ) GOTO 40
    ENDIF

    IF( DS < 0. )   &
     CALL FHW( S_OLD , T_OLD , E0 , FT , P_S , T_S , T_T ,  & ! 収斂後�E計算追加
            F_S , F_T , Y_S , Y_T , W_S , W_T , H_S )
! C IF(NST==10)WRITE(1850,*)S,S_OLD,DS,F_S
    IF( DS < 0. ) GOTO 341

    INCF= 3
    E= 0.0
    CALL STIF( T_S , T_T , G_S , ICIC(4) , C_IC4 , S , T )

    IF( ID_CRACK==0 )THEN
        ID_CRACK= 1
        T_T= FT
        T_S= T_T / E0 + P_S
        CALL FHW( S , T , E0 , FT , P_S , T_S , T_T ,   &  ! 収斂後�E計算追加
                F_S , F_T , Y_S , Y_T , W_S , W_T , H_S )
    ENDIF

    GOTO 9999
! *****************
! **** IUL=13  ****
! *****************
13  IUL= 13
    CALL NN_C( NST,NIT,M,LX,LY,LZ,intI,IUL,NN )

    IF( ID_CRACK == 0 ) THEN
        IF( S > T_S ) GOTO 32
        IF( S > P_S ) GOTO 31
    ELSE
        IF( S > F_S ) GOTO 32
        IF( S > P_S ) GOTO 20
    ENDIF
    IF( S > O_S ) GOTO 11
    IF( S > U_S ) GOTO 12

    IF( DS > 0. )   &
     CALL EXPT( S_OLD , T_OLD , E0 , FT , O_S , ID_CRACK , &  ! 収斂後�E計算追加
              E_S , E_T , X_S , X_T , V_S , V_T , P_S , T_S , T_T)
    IF( DS > 0. ) GOTO 141

    INCF= 1
    E= 0.0
    T= U_T
    GOTO 9999
! *****************
! **** IUL=141 ****
! *****************
141 IUL= 141
    CALL NN_C( NST,NIT,M,LX,LY,LZ,intI,IUL,NN )

    IF( S < E_S ) GOTO 11

    IF( S > P_S ) GOTO 20

    IF( DS < 0. ) THEN
        GOTO 142
    ENDIF

    INCF= 1
    if( E_S == P_S )  then ; write(1850,*) '☁EY_Point:',IUL ; endif
    CALL XY_Point( E_S , E_T , P_S , ZERO , S , T , E , 2)
    GOTO 9999
! *****************
! **** IUL=142 ****
! *****************
142 IUL= 142
    CALL NN_C( NST,NIT,M,LX,LY,LZ,intI,IUL,NN )

    IF( S < X_S ) GOTO 11
    IF( S > P_S ) GOTO 20
    IF( E_S == 0. ) GOTO 11

    IF( DS > 0. ) THEN
        GOTO 141
    ENDIF

    INCF= 1
    if( P_S == X_S )  then ; write(1850,*) '☁EY_Point:',IUL ; endif
    CALL XY_Point( P_S , ZERO , X_S , X_T , S , T , E , 2)
    GOTO 9999
! *****************
! **** IUL=20 ****
! *****************
20  IUL= 20
    CALL NN_C( NST,NIT,M,LX,LY,LZ,intI,IUL,NN )

    IF( S < P_S ) GOTO 142

    IF( S > W_S ) GOTO 342

    IF( ID_CRACK == 0 ) GOTO 31
    IF( DS < 0. ) THEN
        CALL R1R2( S_OLD , T_OLD ,  &
                     F_S , F_T , P_S , ZERO , H_S , ZERO , V_S , V_T ,  &
                    R1_S , R1_T , R2_S , R2_T )
        GOTO 24
    ENDIF

    INCF= 1
    if( P_S == W_S )  then ; write(1850,*) '☁EY_Point:',IUL ; endif
    CALL XY_Point( P_S , ZERO , W_S , W_T , S , T , E , 2)
    GOTO 9999
! *****************
! **** IUL=341 ****
! *****************
341 IUL= 341
    CALL NN_C( NST,NIT,M,LX,LY,LZ,intI,IUL,NN )

    IF( S > F_S ) GOTO 32

    IF( S < H_S ) GOTO 40

    IF( DS > 0. ) THEN
        GOTO 342
    ENDIF

    INCF= 1
    if( F_S == H_S )  then ; write(1850,*) '☁EY_Point:',IUL ; endif
    CALL XY_Point( F_S , F_T , H_S , ZERO , S , T , E , 2)
    GOTO 9999
! *****************
! **** IUL=342 ****
! *****************
342 IUL= 342
    CALL NN_C( NST,NIT,M,LX,LY,LZ,intI,IUL,NN )

    IF( S > F_S ) GOTO 32

    IF( S < H_S ) GOTO 40
    IF( ID_CRACK == 0 ) GOTO 31

    IF( DS < 0. ) THEN
        GOTO 341
    ENDIF

    INCF= 1
    if( Y_S == H_S )  then ; write(1850,*) '☁EY_Point:',IUL ; endif
    CALL XY_Point( Y_S , Y_T , H_S , ZERO , S , T , E , 2)
    GOTO 9999
! *****************
! **** IUL=40 ****
! *****************
40  IUL= 40
    CALL NN_C( NST,NIT,M,LX,LY,LZ,intI,IUL,NN )

    IF( S < V_S ) GOTO 142
    IF( S > H_S ) GOTO 342
    IF( DS > 0. ) THEN
        CALL R1R2( S_OLD , T_OLD ,  &
                H_S , ZERO , V_S , V_T , F_S , F_T , P_S , ZERO ,   &
                R1_S , R1_T , R2_S , R2_T )
        GOTO 42
    ENDIF

    INCF= 1
    if( H_S == V_S )  then ; write(1850,*) '☁EY_Point:',IUL ; endif
    CALL XY_Point( H_S , ZERO , V_S , V_T , S , T , E , 2)
    GOTO 9999
! *****************
! **** IUL=24 ****
! *****************
24  IUL= 24
    CALL NN_C( NST,NIT,M,LX,LY,LZ,intI,IUL,NN )

    IF( S > R1_S ) GOTO 20
    IF( S < R2_S ) GOTO 40

    INCF= 1
    if( R1_S == R2_S )  then ; write(1850,*) '☁EY_Point:',IUL ; endif
2400    CALL XY_Point( R1_S , R1_T , R2_S , R2_T , S , T , E , 2)
    GOTO 9999
! *****************
! **** IUL=42 ****
! *****************
42  IUL= 42
    CALL NN_C( NST,NIT,M,LX,LY,LZ,intI,IUL,NN )

    INCF= 1
    IF( S < R1_S ) GOTO 40
    IF( S > R2_S ) GOTO 20
    GOTO 2400
! *********************************************
! *********************************************
! C
9999    CONTINUE
    S_OLD= S
    T_OLD= T
    IF(E > E0 ) E= E0
    IF(E < E0_100 ) E= E0_100
    END
    ! *********************************************

    SUBROUTINE NN_C( NST,NIT,M,LX,LY,LZ,intI,IUL,NN )
   integer ::  NN(999)
   integer :: NST,NIT,M,LX,LY,LZ,intI,IUL
    integer :: NID
    NN(IUL)=NN(IUL)+1
    IF(NN(IUL) == 1000000 ) THEN
        ! WRITE(1850,'(A40,I5,I2,I5,TR1,3I1,I2,I6)')            &
     ! "☁EUL:NST,NIT,M,LX,LY,LZ,intI,IUL",NST,NIT,M,LX,LY,LZ,intI,IUL
        NID=NID+1
    ENDIF
    IF(NID==4) STOP

    END
    
    !*******************************************************************
    SUBROUTINE EXPT(  S_OLD , T_OLD , E0 , FT , O_S , ID_CRACK ,    &
                E_S , E_T , X_S , X_T , V_S , V_T , P_S , T_S ,T_T )
    
    real (kreal) ::  S_OLD,T_OLD,E0,FT,O_S,E_S,E_T,X_S
    real (kreal) :: X_T,V_S,V_T,P_S,T_S,T_T,ZERO,Edd,xbai
    integer :: ID_CRACK
    
    ZERO= 0.0

! *** [ E点 ] ***
    E_S= S_OLD
    E_T= T_OLD
! *** [ X点 ] ***
    X_S= E_S
    X_T= E_T
! *** [ P点 ] ***
! C xbai= 0. !原点持E��
! C CALL X0_point( E_S , E_T , E0 , xbai ) !初期剛性E0載除荷
    CALL Karsan( E_S , E_T , O_S , E0 , xbai )
    P_S= xbai * E_S
    IF( P_S > 0. ) P_S = 0.
! *** [ V点 ] ***
    IF(E_S==0.) THEN ; V_T=0. ; V_S=0. ; GOTO 10 ; ENDIF
    V_T= E_T / 20. !刁E��修正
    if( P_S == E_S ) then; write(1850,*) '☁EY_Point,EXPT:'; endif
    CALL XY_Point(P_S , ZERO , E_S , E_T , V_T , V_S , Edd , 1 )
10  CONTINUE
! *** [ T点 ] ***
    IF( ID_CRACK==0 )THEN
        T_T= FT
        T_S= T_T / E0 + P_S
    ENDIF

    END
    
    !*******************************************************************
    SUBROUTINE Karsan( E_S , E_T , O_S , E0 , xbai )

    real (kreal) :: E_S,E_T,O_S,E0,xbai,E

    xbai= ( 0.145 * E_S / O_S ) + 0.127
    E= E_T /  (E_S * ( 1. - xbai ))

    IF( E > E0 .OR. E < 0. ) THEN
        !IF( E > E0 ) WRITE(*,*) '●Kaesan�E�P_E > E0'
        !IF( E < 0. ) WRITE(*,*) '●Kaesan�E�E < 0.'
        CALL X0_point( E_S , E_T , E0 , xbai )
    ENDIF
        
    END
    
    SUBROUTINE X0_point( x1 , y1 , E , xbai )
! C
! C 点�E�E1,y1�E�を通り係数Eを持つ直線�EX軸との交点、X座樁E
! C
    real(kreal) :: x1,y1,E,xbai


    xbai= 1. - y1 / ( E * x1 )

    END
    
!*******************************************************************
    SUBROUTINE  XY_Point( x1 , y1 , x2 , y2 , xy , YX , E , N )
! *******************************************************************
! * 点(x1,y1)、Ex2,y2)を通る直線上で、E
! *       剛性 E および
! * 1:  xy にあためEX を求める、E
! * 2:  xy にあためEY を求める、E
    real(kreal) :: x1,y1,x2,y2,xy,YX,E,dx,dy,c
    integer :: i,id,N

    dx = x1 - x2
    dy = y1 - y2

    do i=1,2
        if(i==1 ) then ; id= 6 ; else ;id= 1850 ; endif
        IF( dx ==0. ) write(id,*) '☁EXY_Point�E�Ex1 = x2'
        IF( N==1 ) THEN
            IF( dy ==0. ) write(id,*) '☁EXY_Point�E�Ey1 = y2'
        ENDIF
    end do

    E= dy / dx
    c  = x1 * y2 - x2 * y1
    IF( N==1 ) THEN
        YX = dx / dy * xy - c / dy
    ELSEIF( N==2 ) THEN
        YX = dy / dx * xy + c / dx
    ENDIF

    END
    
    !*******************************************************************
    SUBROUTINE  FHW( S_OLD , T_OLD , E0 , FT , P_S , T_S , T_T ,    &
                    F_S , F_T , Y_S , Y_T , W_S , W_T , H_S )
    
        real(kreal) ::  S_OLD,T_OLD,E0,FT,P_S,T_S,T_T
        real (kreal) :: F_S,F_T,Y_S,Y_T,W_S,W_T,H_S,ZERO,EDD
! *** [ F点 ] ***
    F_S= S_OLD
    F_T= T_OLD
! *** [ Y点 ] ***
    Y_S= F_S
    Y_T= F_T
! *** [ H点 ] ***
! C H_S= T_S - T_T / E0 ! 原点持E��
! C H_S= F_S - F_T / E0 ! 初期剛性載除荷
    CALL Naganuma( F_S , F_T , T_S , E0 , H_S )
! *** [ W点 ] ***
    IF( F_T == 0. ) THEN ; W_S= S_OLD ; W_T= T_OLD ; GOTO 10 ; ENDIF
    W_T=  F_T / 2.0 !刁E��修正
    if( H_S == F_S ) then; write(1850,*) '☁EY_Point,FHW:'; endif
    CALL XY_Point(H_S , ZERO , F_S , F_T , W_T , W_S , Edd , 1 )
10  CONTINUE

    END
    
    
    SUBROUTINE Naganuma( F_S , F_T , T_S , E0 , H_S )
    real (kreal) :: F_S,F_T,T_S,E0,H_S,Scr,E

    Scr= 100. / 10.**6.
    E= Scr / (Scr + F_S - T_S ) * E0
    H_S= F_S - F_T / E

    END
    
    !*******************************************************************
    SUBROUTINE STIF( T_S , T_T ,G_S , IC4 , C_IC4 , S , T )

    real (kreal) ::T_S,T_T,G_S,C_IC4,S,T,SX,TX,EU_CR
    integer :: IC4

    IF( IC4 == 0 ) GOTO 1000 ! 白井弁E
    IF( IC4 == 3 ) GOTO 3000 ! 持E��弁E
! C IF( IC4 == 4 ) GOTO 4000 ! 一定剛性

! *******************
! **** IC4=1000  ****
! *******************
1000    SX= ( S - T_S ) / ( G_S - T_S )
    TX=1.0 - 2.748 * SX +2.654 * SX**2 - 0.906 * SX**3
    T= TX * T_T
    GOTO 9999

! *******************
! **** IC4=3000  ****
! *******************
3000    EU_CR= S - T_S
    IF( EU_CR <= 0. ) THEN
        T= 0.0
        GOTO 9999
    ENDIF
    T= ( 0.0002 / EU_CR )**C_IC4 * T_T
    IF(T < 0.) T= 0.
    IF(T > T_T ) T= T_T
    GOTO 9999
! C
! C4000 EU_CR= S - T_S
! C T = C_IC4 * EU_CR +  T_T
! C GOTO 9999
! c
9999    CONTINUE
    END
    
    !*******************************************************************
    SUBROUTINE  R1R2( x , y , &
                x1_1 , y1_1 , x2_1 , y2_1 , x1_2 , y1_2 , x2_2 , y2_2, &
                X_1 , Y_1 , X_2 , Y_2 )
    
    real (kreal) ::  x,y,x1_1,y1_1,x2_1,y2_1,x1_2,y1_2,x2_2,y2_2
    real (kreal) :: X_1 , Y_1 , X_2 , Y_2 ,R, EDD

    X_1= x
    Y_1= y

    IF(x2_1 == x1_1) THEN
        !WRITE(*   ,*)'☁E1R2�E�x2_1 == x1_1'
        !WRITE(*,*)'☁E1R2�E�x2_1 == x1_1'
    ENDIF
    R= (X_1 - x1_1 ) / ( x2_1 - x1_1 )
    X_2= (x2_2 - x1_2 ) * R + x1_2
    !if(x1_2 == x2_2) then; write(*,*) '☁EY_Point,R1R2:'; endif
    CALL XY_Point( x1_2 , y1_2 , x2_2 , y2_2 , X_2 , Y_2 , Edd , 2 )

    END
    
    !********************************************************************

    SUBROUTINE sTransRToA(R1,R2,R3,RA,NDIM)
    ! *引用名：sTransAToR
! *機�E�E�整数配�Eから引数に戻ぁE
! *入力！ENARRAY   :
! *出力！EN1,N2,N3 :
! *目皁E��CODING量を減らし、簡単なMISSを行わなぁE��めE
! ********************************************************************
    real(kreal) R1,R2,R3
    integer NDIM
    real (kreal) RA(NDIM)
    
    RA(1)=R1
    RA(2)=R2
    RA(3)=R3

    END
    
    ! ********************************************************************
    SUBROUTINE sTransAToR(R1,R2,R3,RA,NDIM)
    ! *引用名：SBACKFROMAINT
! *機�E�E�整数配�Eから引数に戻ぁE
! *入力！ENARRAY   :
! *出力！EN1,N2,N3 :
! *目皁E��CODING量を減らし、簡単なMISSを行わなぁE��めE
! ********************************************************************
    real(kreal) R1,R2,R3
    integer NDIM
    real (kreal) RA(NDIM)
    
    R1=RA(1)
    R2=RA(2)
    R3=RA(3)

    END
    
    !*******************************************************************
    SUBROUTINE CompressZone( IC1 , E0 , O_S , O_T , S , T , E )
    real(kreal) :: E0,O_S,O_T,S,T,E 
    integer :: IC1
    
    real(kreal) :: Z,ZZ

    IF( IC1 > 2 .AND. IC1 < 1 ) GOTO  100 ! Saenz 弁E
    IF( IC1 == 1 ) GOTO 1000 ! Fafitis-Shag 弁E
    IF( IC1 == 2 ) GOTO 2000
! *******************
! **** IC4=100  ****
! *******************
100 Z  =1.+( E0 * O_S / O_T -2. )* S / O_S +( S / O_S )**2.
    T= E0 * S / Z
    E =E0*(1.-( S / O_S )**2. )/ Z**2.
    GOTO 9999
! *******************
! **** IC4=1000  ****
! *******************
1000  ZZ = S / O_S
      Z  = E0 * O_S / O_T
      T  = O_T * ( 1.-( 1.-ZZ )**Z )
      E  = E0 * ( 1.-ZZ )**( Z-1. )
    GOTO 9999
! *******************
! **** IC4=2000  ****
! *******************
2000    CONTINUE
    Z  =1.+( E0 * O_S / O_T -2. )* S / O_S +( S / O_S )**2.
    T= E0 * S / Z
    E =E0*(1.-( S / O_S )**2. )/ Z**2.
! C CALL Kimura(EO,DUMMY,ECU,SC,EEN,SEN)
    GOTO 9999

9999    CONTINUE
    END

    
    !*******************************************************************
      SUBROUTINE CompressionZone_REVERS(ICIC1,FC,SJ,            &
                                    EO,EJJ,ECU,SC,EJ,SJDAMMY)
! *******************************************************************
! ****                          hj . 2008.8.30                   ****
! *******************************************************************
    real (kreal) :: FC,SJ,EO,EJJ,ECU,SC,EJ,SJDAMMY
    integer :: ICIC1

    IF(     ICIC1==1) THEN
        CALL SHAH_REVERS(EJ,SJ,EO,ECU,FC)     !FAFITIS-SHAH式送E��数
        CALL FS(EO,EJJ,ECU,SC,EJ,SJDAMMY)     !J点の接線剛性
    ELSE IF(ICIC1==2) THEN
        CALL Kimura_REVERS(EJ,SJ,EO,ECU,FC)
        CALL Kimura(EO,EJJ,ECU,SC,EJ,SJDAMMY)
    ELSE
        CALL SAENZ_REVERS(EJ,SJ,EO,ECU,FC)    !SAENZ式送E��数
        CALL SAENZ(EO,EJJ,ECU,SC,EJ,SJDAMMY)  !J点の接線剛性
      ENDIF

    END
    
    
    !*********************************************************************
      SUBROUTINE SAENZ_REVERS(E,S,EO,ECU,FC)
! *********************************************************************
! ***SAENZ式�E送E��数                                                ***
! *********************************************************************
! *$INCLUDE"D:\ANALYSIS\ORIGINAL\TS\H\FILE.H"
! *$INCLUDE'/USER2/ZAKKY/Z9/H/FILE.H'
    real (kreal) :: E,S,EO,ECU,FC
    real (kreal) :: ALFA,ELINE
      ! CALL ROOT_CHECK(0000, 94)

      ALFA=0.0

      ALFA=ECU/S*EO-EO*ECU/FC+2.0
      E   =ECU*0.5*(ALFA-(ALFA**2-4.0)**0.5)
    
    ELINE=S/EO
    IF(E.GT.ELINE) THEN
        E=ELINE
    ENDIF

      RETURN
      END
    
    !*********************************************************************
      SUBROUTINE SHAH_REVERS(E,S,EO,ECU,FC)
! *********************************************************************
! ***SHAH式�E送E��数                                                 ***
! *********************************************************************
! *$INCLUDE"D:\ANALYSIS\ORIGINAL\TS\H\FILE.H"
! *$INCLUDE'/USER2/ZAKKY/Z9/H/FILE.H'
    real (kreal) :: E,S,EO,ECU,FC
    real (kreal) :: ALFA,ZZZ
      ! CALL ROOT_CHECK(0000, 95)

      ALFA=0.0
      ZZZ =0.0

      ZZZ =EO*ECU/FC
      ALFA=(LOG(1.-S/FC))/ZZZ
      E   =ECU*(1.0-EXP(ALFA))

      RETURN
      END
        
    !*********************************************************************
      SUBROUTINE Kimura_REVERS( S , T , E0 , O_S , O_T )
! *********************************************************************
    real (kreal) :: S,T,E0,O_S,O_T
    real (kreal) :: R0,Smax,Smin,R,T0,dle
    
    R0 = ABS ( O_T / 1000. )
    Smax = 0.
    Smin = O_S
    S = ( Smax + Smin ) / 2.

10  CALL Kimura( E0 , dle , O_S , O_T , S , T0 )
    R = ABS( T - T0 )
    IF( R > R0 ) THEN
        IF( T > T0 ) THEN
            Smin = S
            S = ( Smax + Smin ) / 2.
        ELSE
            Smax = S
            S = ( Smax + Smin ) / 2.
        ENDIF
        GOTO 10
    ENDIF

      END
    
    !*******************************************************************
      SUBROUTINE SHEARREDUCE(NST,G0,NG3C,KUL,KULIN1,KULIN5,         &
        PSICR,PSICROLD,EU,EUOLD,ECCR,CCRACK,LCR,NTCRACK,            &                               
        FC,EO,FT,                                                   &
        DS1,DE1,DS2,DE2,                                            &
        FS1,FE1,FS2,FE2,                                            &
        RS1,RE1,RS2,RE2,                                            &
        RG12,RG13,RG23)
! **********************************************************************
! C KUL,KULIN1,KULIN5,          !履歴番号、履歴1を通ったとぁE��持E��、履歴2を通ったとぁE��持E��E
! C     *PSICR,EU,ECCR,CCRACK,      !等価�E�軸ひずみ、�Eび割れ発生時ひずみ、引張先行OR圧縮先行�E判宁E
! C     *LCR,                       !ひび割れが発生してるかどぁE��の持E��E
! C     *TE,TELL,                   !現スチE��プせん断ひずみ、前スチE��プせん断ひずみ
! C     *FC,EO,FT,                      !コンクリート圧縮強度、ヤング係数、引張強度
! C     *DS1,DE1,DS2,DE2,               !負�E�正�E�域D点せん断応力、せん断ひずみ
! C     *FS1,FE1,FS2,FE2,               !正�E�負�E�域F点せん断応力、せん断ひずみ
! C     *RS1,RE1,RS2,RE2,               !正�E�負�E�域Rせん断応力、せん断ひずみ
! C     *RG12,RG13,RG23)                !せん断剛性�E�アウト�EチE���E�E
! **********************************************************************
! ****              G  : ELASTIC SHEAR MODULUS                      ****
! ****              RG : REDUCED SHEAR MODULUS                      ****
! ****                                                              ****
! ****コンクリート骨材�E噛合ぁE��果を老E�Eしたせん断剛性繰返しモチE��  ****
! ****�E�鉄筋�Eダボ作用効果を老E�Eしたせん断伝達モチE��は、現時点で、E ****
! ****   コンクリート�E積層要素がなぁE��めモチE��化が困難である。！E  ****
! ****                                                              ****
! ****                     REPRODUCED BY RYOTARO HIGASHI 2003.12.18 ****
! **********************************************************************
    real (kreal) :: G(3),RG(3),DR(3),CCRACK(3),DUMMY(3)
    real (kreal) :: TE(3),TELL(3),EU(6),EUOLD(6),ECR(3),ECCR(3),CCRACKDM(3)
    real (kreal) :: PSIT(3),PSITL(3),PSICR(6),PSICROLD(6),RS1(3),RE1(3),RS2(3),RE2(3)
    real (kreal) :: DS1(3),DE1(3),DS2(3),DE2(3),FS1(3),FE1(3),FS2(3),FE2(3)
    real (kreal) :: G0,FC,EO,FT,RG12,RG13,RG23,ATemp,eecr
    integer :: KUL(3),KULIN1(3),KULIN5(3)
    
    integer :: LCR(3),LCRDM(3),NST,NG3C,NTCRACK
    
    real (kreal) :: G23,G13,G12
    integer :: inti,ii
    
    
    NG3C=0
! C ***[DEBUG用]***
    IF(LCR(1).EQ.1.OR.LCR(2).EQ.1.OR.LCR(3).EQ.1) THEN
        CONTINUE
    ENDIF

! C ***[ひび割れ�EなぁE��のせん断剛性�E�せん断弾性剛性�E�]***
    G23=0.5D0*EO/(1.0D0+0.2D0)
    G13=0.5D0*EO/(1.0D0+0.2D0)
    G12=0.5D0*EO/(1.0D0+0.2D0)

    RG23=G23
    RG13=G13
    RG12=G12

! C ***[ここの計算�E実際はEUおよびPSICRの4-6につぁE��の計算]***
! C ***[LCRの頁E��を基準にTE,TELL,PSITの頁E��を調整する、要するに12ↁE、E3ↁE、E1ↁE]***
! C ***[ひび割れてぁE��ぁE��向�E計算�EしなぁETSIから変換した数値をそのまま用ぁE��]***
    PSIT(1)=PSICR(5)
    PSIT(2)=PSICR(6)
    PSIT(3)=PSICR(4)
    PSITL(1)=PSICROLD(5)
    PSITL(2)=PSICROLD(6)
    PSITL(3)=PSICROLD(4)

    TE(1)=EU(5)
    TE(2)=EU(6)
    TE(3)=EU(4)

    TELL(1)=EUOLD(5)
    TELL(2)=EUOLD(6)
    TELL(3)=EUOLD(4)

    LCRDM(1)=LCR(1)
    LCRDM(2)=LCR(2)
    LCRDM(3)=LCR(3)

    CCRACKDM(1)=CCRACK(1)
    CCRACKDM(2)=CCRACK(2)
    CCRACKDM(3)=CCRACK(3)

! C ***[!もしひび割れが閉じてぁE��ら、剛性を復活させる]***
    DO intI=1,3
        IF(NTCRACK.EQ.1.AND.LCRDM(intI).EQ.1) THEN        
! C     IF(LCRDM(intI).NE.1) THEN             
            KUL(intI)=0
        ENDIF
    ENDDO
    IF (LCRDM(1).NE.1.AND.LCRDM(2).NE.1.AND.LCRDM(3).NE.1) GOTO 9000

! C     ***[!ひび割れ発生時ひずみ�E�単調載荷の場合�EこれでぁE��が繰り返し載荷の場合には適用不可�E�]***
    EECR=FT/EO     
    DO II=1,3
        IF(CCRACKDM(II).EQ.4.0) ECR(II)=EECR
        IF(CCRACKDM(II).EQ.5.0) ECR(II)=ECCR(II)
    END DO

! C ***[!増�Eせん断ひずみの計算]***
    DO II=1,3
        IF(KUL(II).EQ.0)THEN
            TELL(II)=0.0
        ENDIF
        DR(II)=TE(II)-TELL(II)        
    END DO

    CALL sTransRToA(G23,G13,G12,G,3)
    CALL sTransRToA(RG23,RG13,RG12,RG,3)

     DO 40 II=1,3

! C     ***[ひび割れてぁE��ぁE��向�E計算�EしなぁETSIから変換した数値をそのまま用ぁE��]***
        IF(NTCRACK.GE.2.OR.NTCRACK.EQ.0)THEN
            CONTINUE
        ELSEIF(LCRDM(II).EQ.1)THEN
! C     ELSEIF(LCRDM(II).NE.1)THEN
            KUL(II)=0.0
            PSIT(II)=0.0
! C         PSIT(II)=PSITL(II)+DR(II)*G(II)
            GOTO 30         
        ENDIF

        IF(DR(II).EQ.0.0) GOTO 30

        IF(KUL(II).EQ.1) GOTO 1
        IF(KUL(II).EQ.2) GOTO 2
        IF(KUL(II).EQ.3) GOTO 3
        IF(KUL(II).EQ.4) GOTO 4
        IF(KUL(II).EQ.5) GOTO 5
        IF(KUL(II).EQ.6) GOTO 6
        IF(KUL(II).EQ.7) GOTO 7
        IF(KUL(II).EQ.8) GOTO 8

        IF(TE(II).GE.0.0.AND.KULIN1(II).NE.0) THEN
            RE1(II)=0.0
            GOTO 4
        ENDIF
        IF(TE(II).LT.0.0.AND.KULIN5(II).NE.0) THEN
            RE2(II)=0.0
            GOTO 8
        ENDIF

! C     ***[KULは0の場合増�Eせん断ひずみによる1ぁEかにぁE��]***
        IF(DR(II).GT.0.0) GOTO 1020
        IF(DR(II).LT.0.0) GOTO 5010


! C     ***[!KUL=1を通過したとぁE��持E��]***
    1       KULIN1(II)=1                     
        IF(DR(II)) 1010,1020,1020
! C     ***[!正埁EKUL=1)で増�Eが負]***
 1010       CONTINUE                         
        CALL G3_DF(FC,ECR(II),EU(II),TELL(II),G(II),    &
            DS1(II),DE1(II),FS1(II),FE1(II))

        IF(TE(II).GT.FE1(II)) GOTO 2010
        IF(TE(II).GT.0.0)    GOTO 3010
        IF(KULIN5(II).EQ.0)THEN
             GOTO 5010
        ELSE
        RE2(II)=0.0
         IF(TE(II).GT.DE2(II)) GOTO 8010
         IF(TE(II).LE.DE2(II)) GOTO 5010
        ENDIF

 1020       CONTINUE                                                                !正埁EKUL=1)で増�Eが正
        KUL(II)=1
        CALL G3_15(NST,G0,NG3C,KUL(II),FC,EU(II),ECR(II),TE(II),    &
            G(II) ,RG(II),PSIT(II))
        GOTO 40

    2       IF(DR(II)) 2010,2020,2020
 2010       CONTINUE                                            !正埁EKUL=2)で増�Eが負
        IF(TE(II).GT.FE1(II))THEN
         KUL(II)=2
        CALL G3_2468(KUL(II),FC, EU(II),   ECR(II),  TE(II),    &
            DS1(II),DE1(II),FE1(II),G(II),RG(II),PSIT(II))
          GOTO 40
        ENDIF
        IF(TE(II).GT.0.0) GOTO 3010
        IF(KULIN5(II).EQ.0)THEN
           GOTO 5010
        ELSE
          RE2(II)=0.0
         IF(TE(II).GT.DE2(II)) GOTO 8010
         IF(TE(II).LE.DE2(II)) GOTO 5010
        ENDIF

 2020  CONTINUE                                          !正埁EKUL=2)で増�Eが正
        IF(TE(II).GT.DE1(II))THEN
       GOTO 1020
        ELSE
         KUL(II)=2
            CALL G3_2468(KUL(II),FC, EU(II),   ECR(II),  TE(II),    &
                DS1(II),DE1(II),FE1(II),G(II),RG(II),PSIT(II))
          GOTO 40
        ENDIF

    3   IF(DR(II)) 3010,3020,3020
 3010  CONTINUE                                                                !正埁EKUL=3)で増�Eが負
        IF(TE(II).GT.0.0)THEN
            KUL(II)=3
            PSIT(II)=0.0
            RG(II)=0.1
            GOTO 40
        ENDIF
        IF(KULIN5(II).EQ.0)THEN
           GOTO 5010
        ELSE
          RE2(II)=0.0
         IF(TE(II).GT.DE2(II)) GOTO 8010
         IF(TE(II).LE.DE2(II)) GOTO 5010
        ENDIF

 3020  CONTINUE                                                                !正埁EKUL=3)で増�Eが正
        RE1(II)=TELL(II)
        RS1(II)=0.0
         IF(TE(II).GT.DE1(II))THEN
           GOTO 1020
         ELSE
           GOTO 4020
         ENDIF

    4   IF(DR(II)) 4010,4020,4020
 4010  CONTINUE                                                                !正埁EKUL=4)で増�Eが負
         IF(TE(II).GT.RE1(II))THEN
         KUL(II)=4
            CALL G3_2468(KUL(II),FC, EU(II),   ECR(II),  TE(II),    &
                DS1(II),DE1(II),RE1(II),G(II),RG(II),PSIT(II))
          GOTO 40
        ENDIF
         IF(TE(II).GT.0.0) GOTO 3010
         IF(KULIN5(II).EQ.0)THEN
            GOTO 5010
         ELSE
          RE2(II)=0.0
          IF(TE(II).GT.DE2(II)) GOTO 8010
          IF(TE(II).LE.DE2(II)) GOTO 5010
         ENDIF

 4020  CONTINUE                                                                !正埁EKUL=4)で増�Eが正
         IF(TE(II).GT.DE1(II))THEN
            GOTO 1020
         ELSE
         KUL(II)=4
            CALL G3_2468(KUL(II),FC, EU(II),   ECR(II),  TE(II),    &
                DS1(II),DE1(II),RE1(II),G(II),RG(II),PSIT(II))
          GOTO 40
        ENDIF

    5  KULIN5(II)=1                                                            !KUL=5を通過したとぁE��持E��E
        IF(DR(II)) 5010,5020,5020
 5010  CONTINUE                                                                !負埁EKUL=5)で増�Eが負
        KUL(II)=5
        CALL G3_15(NST,G0,NG3C,KUL(II),FC,EU(II),ECR(II),TE(II),    &
            G(II),RG(II),PSIT(II))
      GOTO 40

 5020  CONTINUE                                                                !負埁EKUL=5)で増�Eが正
        CALL G3_DF(   FC,  EU(II),  ECR(II), TELL(II),    G(II),    &
                    DS2(II),  DE2(II),  FS2(II),  FE2(II))
         IF(TE(II).LT.FE2(II)) GOTO 6020
         IF(TE(II).LT.0.0) GOTO 7020
         IF(KULIN1(II).EQ.0)THEN
            GOTO 1020
         ELSE
          RE1(II)=0.0
          IF(TE(II).LT.DE1(II)) GOTO 4020
          IF(TE(II).GE.DE1(II)) GOTO 1020
         ENDIF

    6  IF(DR(II)) 6010,6020,6020
 6010   CONTINUE                                                               !負埁EKUL=6)で増�Eが負
        IF(TE(II).LT.DE2(II))THEN
            GOTO 5010
        ELSE
         KUL(II)=6
            CALL G3_2468(KUL(II),FC, EU(II),ECR(II),TE(II),     &
                DS2(II),DE2(II),FE2(II),G(II),RG(II),PSIT(II))
          GOTO 40
        ENDIF

 6020   CONTINUE                                                               !負埁EKUL=6)で増�Eが正
        IF(TE(II).LT.FE2(II))THEN
         KUL(II)=6
            CALL G3_2468(KUL(II),FC, EU(II),ECR(II),TE(II),         &
                DS2(II),DE2(II),FE2(II),G(II),RG(II),PSIT(II))
          GOTO 40
        ENDIF
        IF(TE(II).LT.0.0) GOTO 7020
        IF(KULIN1(II).EQ.0)THEN
           GOTO 1020
        ELSE
         RE1(II)=0.0
         IF(TE(II).LT.DE1(II)) GOTO 4020
         IF(TE(II).GE.DE1(II)) GOTO 1020
        ENDIF

    7   IF(DR(II)) 7010,7020,7020
 7010   CONTINUE                                                               !負埁EKUL=7)で増�Eが負
         RE2(II)=TELL(II)
         RS2(II)=0.0
        IF(TE(II).LT.DE2(II))THEN
           GOTO 5010
        ELSE
           GOTO 8010
        ENDIF

 7020   CONTINUE                                                               !負埁EKUL=7)で増�Eが正
        IF(TE(II).LT.0.0)THEN
            KUL(II)=7
            PSIT(II)=0.0
            RG(II)=0.1  
            GOTO 40

        ENDIF
        IF(KULIN1(II).EQ.0)THEN
           GOTO 1020
        ELSE
         RE1(II)=0.0
         IF(TE(II).LT.DE1(II)) GOTO 4020
         IF(TE(II).GE.DE1(II)) GOTO 1020
        ENDIF

    8   IF(DR(II)) 8010,8020,8020
 8010   CONTINUE                                                               !負埁EKUL=8)で増�Eが負
        IF(TE(II).LT.DE2(II))THEN
           GOTO 5010
        ELSE
         KUL(II)=8
            CALL G3_2468(KUL(II),FC, EU(II),ECR(II),TE(II),     &
                DS2(II),DE2(II),RE2(II),G(II),RG(II),PSIT(II))
          GOTO 40
        ENDIF

 8020   CONTINUE                                                               !負埁EKUL=8)で増�Eが正
        IF(TE(II).LT.RE2(II))THEN
         KUL(II)=8
            CALL G3_2468(KUL(II),FC, EU(II),ECR(II),TE(II),         &
                DS2(II),DE2(II),RE2(II),G(II),RG(II),PSIT(II))
          GOTO 40
        ENDIF
        
        GOTO 7020

        IF(KULIN1(II).EQ.0) THEN
           GOTO 1020
        ELSE
         RE1(II)=0.0
         IF(TE(II).LT.DE1(II)) GOTO 4020
         IF(TE(II).GE.DE1(II)) GOTO 1020
        ENDIF

        TELL(II)=TE(II)

   30 CONTINUE

   40 CONTINUE

    PSICR(4)=PSIT(3)
    PSICR(5)=PSIT(1)
    PSICR(6)=PSIT(2)

    CALL sTransAToR(G23,G13,G12,G,3)
    CALL sTransAToR(RG23,RG13,RG12,RG,3)
!C  WRITE(6111,*)RG23,RG13,RG12

    DO intI=1,3
        ATemp=PSIT(intI)*TE(intI)
        IF(ATemp.LT.0.0.AND.LCRDM(intI).NE.0)THEN
            CONTINUE
        ENDIF
    ENDDO
 9000   CONTINUE
        RETURN
        END 
    
    !*******************************************************************
      SUBROUTINE G3_DF(   FC,   ECR,     EU,  TELL,   G,    &
                     DS1_5, DE1_5,   FS1_5,  FE1_5)
! *******************************************************************
! **** こ�EサブルーチンはKUL=1 OR 5で除荷されたときに呼び出され,****
! **** D点とF点のひずみを計算する、E                            ****
! *******************************************************************
! *     FC        :コンクリーチE軸圧縮強度
! *    EU        :ひび割れ方向等価�E�軸ひずみ
! *    ECR        :ひび割れ発生時ひずみ
! *   TELL　      :せん断ひずみ
! *    G12        :弾性せん断剛性
! *******************************************************************
    real (kreal) ::  FC,ECR,EU,TELL,G,DS1_5,DE1_5,FS1_5,FE1_5
    
    real (kreal) :: A,B,C,STS,BETA
     B=1./3.
     C=(-1*FC)**B
       STS=18.01024*C     
       A=EU-ECR !ひび割れ幁E��現スチE��プ�E等価�E�軸ひずみーひび割れ発生時ひずみ�E�E
    A=0.005  
    BETA=2.0*TELL/A

! C TS=18.01024*C*(BETA**2)/(1.0+BETA**2)
       IF(TELL.GT.0.0)THEN
         DS1_5=(STS*BETA**2)/(1.0+BETA**2)                !D点応力
       ELSE
         DS1_5=-1.*(STS*BETA**2)/(1.0+BETA**2)                !D点応力
       ENDIF
         DE1_5=TELL                                       !D点ひずみ

      FE1_5=0.5*DE1_5                                 !F点ひずみ
      FS1_5=0.0                                       !F点応力

! C        CHECK=DE1_5-4.0*DS1_5/G
! C         IF(ABS(FE1_5).GT.ABS(CHECK)) FE1_5=CHECK

       RETURN
       END

    !*******************************************************************
      SUBROUTINE G3_15(NST,G0,NG3C,KUL,FC,EU,ECR,RNOW,G,G3,TS)
! *******************************************************************
! ****                                                           ****
! ****         KUL=1 OR 5を進むとき�E、履歴上�EG3を求めめE         ****
! ****                                                           ****
! *******************************************************************
     real (kreal) :: G0,FC,EU,ECR,RNOW,G,G3,TS
     integer :: NST,NG3C,KUL
     
     real (kreal) ::A,DRNOW,B,C,STS,BETA,G3C

     
      A=EU-ECR
! C A=0.005
    DRNOW=ABS(RNOW)                                    !せん断ひずみ
    B=1./3.
    C=(-1*FC)**B
    STS=18.01024*C                                         !KG/CM2単佁E
    BETA=2.0*DRNOW/A                                !ひび割れ方向変佁E

! C ***[ひび割れ後�Eせん断剛性]***
    IF (DRNOW.NE.0.0)THEN
        G3C=STS*(BETA**2)/((1.0+BETA**2)*DRNOW)    
! C     IF(NST > 121.AND.NST<281)  THEN
! C         BETA_=1000.
! C         DRNOW_ = 0.0000001
! C         BETA_=BETA
! C         DRNOW_ = DRNOW / 1000.
! C         G3C=STS*(BETA_**2)/((1.0+BETA_**2)*DRNOW_)  
! C         NG3C=1
! C     ELSE IF(NST > 281.AND.NST<381)  THEN
! C         NG3C=2
! C     ELSE IF(NST > 381.AND.NST<481)  THEN
! C         NG3C=3
! C     ENDIF
    ELSE
        TS=0.0
        RETURN
    ENDIF

! C ***[ひび割れてぁE��ぁE��ンクリートとの連続性を確保]***
      G3C  =G3C*G/(G3C+G)                              
      IF(G3C.GT.G) G3C=G

      G3=G3C
      IF(G3.GT.G) G3=G
      IF(G3.LT.0.0) G3=0.00001

    TS=18.01024*C*(BETA**2)/(1.0+BETA**2)
    IF(KUL.GE.5) THEN
        TS=-TS
    ENDIF

      RETURN
      END
    
!*******************************************************************
      SUBROUTINE G3_2468(KUL,FC,EU,ECR,R,DS,DE,FE,G,G3,TS)      
! C                       (KUL(II),FC, EU(II),ECR(II),TE(II),       
! C     *           DS2(II),DE2(II),RE2(II),G(II),RG(II),PSIT(II))
! *******************************************************************
! ****                                                           ****
! ****        KUL=2 OR 4 OR 6 OR 8を進むとき�E、履歴上�EG3を求めめE    ****
! ****                                                           ****
! *******************************************************************
     real (kreal) :: FC,EU,ECR,R,DS,DE,FE,G,G3,TS
    real (kreal) :: FFA,FFB,G3C
    integer :: KUL
    
         IF(DE.EQ.FE) GOTO 9999
         FFA=DS/((DE-FE)**4)                        !係数A
         FFB=FE                                     !係数B

         G3C=4.0D0*FFA*(R-FFB)**3.   !コンクリート�E噛合ぁE��って生じるせん断剛性

         IF(G3C.GT.G) G3C=G

         G3=G3C
        IF(G3.GT.G) G3=G
        IF(G3.LT.0.0) G3=0.0
    
        TS=FFA*(R-FFB)**4
        IF(KUL.GT.5) THEN
            IF(FE.GE.0.0)THEN
                CONTINUE                
            ENDIF
            IF(TS.GE.0.0)THEN
                CONTINUE
            ENDIF
        ENDIF
    
        IF(KUL.LT.5) THEN
            IF(FE.LT.0.0)THEN
                CONTINUE                
            ENDIF
        ENDIF

 9999    RETURN
         END
    
    
    ! ********************************************************************

    SUBROUTINE sPSITurn(EU,PSICR,DCCR,EUOUT,PSICROUT,   &
        IUL,KUL,IULOUT,KULOUT)
     
     ! *引用名：sCoordCHG2
! *機�E�E��Eび割れ暦ぁE回ある場合�E座標系の計箁E
! *入力！EDCCRCHG  :ひび割れ暦を記録する座標�EトリチE��ス
! *     DCCR        :ひびわれ座標系
! *     DCCRL       :前�Eびわれ座標系
! *     NTCRACK :ひび割れ総本数
! *     LCRACK  :吁E��のひび割れ状況E
! *     IJNC        :座標系を変換されたかどぁE��の持E��E
! *     IACTIVEC:アクチE��ブする座標系の番号
! *出力！EDCCR     :前stepの等価一軸ひずみ�E��E力として�E�E
! *     DAIX        :コンクリート�E応力増�E�E��Eび割れ軸�E�E
! *     DAIXL       :前stepコンクリート�E応力�E��Eび割れ軸�E�E
! *　　　  PSICR       :ひび割れ座標系での応力
! ********************************************************************

    real (kreal) :: DCCR(3,3),DCCRTEMP(3,3),EU(6),EUOUT(6),EUTEMP(6)
    real (kreal) :: PSICR(6),PSICROUT(6),PSICRTEMP(6)
    integer :: IUL(3),KUL(3),IULOUT(3),KULOUT(3)
    
    
    real (kreal) :: PSIMAX,PSIMIN
    integer :: III,inti, INTPSIMAX, INTPSIMID,INTPSIMIN
    
    

! C RETURN
    III=0
    IF(III.EQ.0)THEN
        DO intI=1,3
            IULOUT(intI)=IUL(intI)
            KULOUT(intI)=KUL(intI)
        ENDDO
        
        DO intI=1,6
            EUOUT(intI)=EU(intI)
            PSICROUT(intI)=PSICR(intI)
        ENDDO
        
        RETURN
    ENDIF
    PSIMAX=MAX(PSICR(1),PSICR(2),PSICR(3))
    PSIMIN=MIN(PSICR(3),PSICR(2),PSICR(1))
    DO INTI=1,3
        IF (PSICR(INTI).EQ.PSIMAX) THEN
            INTPSIMAX=INTI
            EXIT    
        ENDIF
    ENDDO
    DO INTI=3,1,-1
        IF (PSICR(INTI).EQ.PSIMIN) THEN
            INTPSIMIN=INTI
            EXIT    
        ENDIF
    ENDDO
    INTPSIMID = 6-INTPSIMAX-INTPSIMIN

    DO INTI=1,3
        DCCRTEMP(1,INTI)=DCCR(INTPSIMAX,INTI)
        DCCRTEMP(2,INTI)=DCCR(INTPSIMID,INTI)
        DCCRTEMP(3,INTI)=DCCR(INTPSIMIN,INTI)
    ENDDO

    IULOUT(1)=IUL(INTPSIMAX)
    IULOUT(2)=IUL(INTPSIMID)
    IULOUT(3)=IUL(INTPSIMIN)

    KULOUT(1)=KUL(INTPSIMAX)
    KULOUT(2)=KUL(INTPSIMID)
    KULOUT(3)=KUL(INTPSIMIN)

      CALL CHG3(EUTEMP(1),EUTEMP(2),EUTEMP(3),  &
        EUTEMP(4),EUTEMP(5),EUTEMP(6),          &
        DCCR(1,1),DCCR(1,2),DCCR(1,3),          &
        DCCR(2,1),DCCR(2,2),DCCR(2,3),          &
        DCCR(3,1),DCCR(3,2),DCCR(3,3),          &
        EU(1),EU(2),EU(3),                      &
        EU(4),EU(5),EU(6),2)

      CALL CHG3(EUOUT(1),EUOUT(2),EUOUT(3),             &
        EUOUT(4),EUOUT(5),EUOUT(6),                     &
        DCCRTEMP(1,1),DCCRTEMP(1,2),DCCRTEMP(1,3),      &
        DCCRTEMP(2,1),DCCRTEMP(2,2),DCCRTEMP(2,3),      &
        DCCRTEMP(3,1),DCCRTEMP(3,2),DCCRTEMP(3,3),      &
        EUTEMP(1),EUTEMP(2),EUTEMP(3),                  &
        EUTEMP(4),EUTEMP(5),EUTEMP(6),3)

      CALL CHG3(PSICRTEMP(1),PSICRTEMP(2),PSICRTEMP(3), &
        PSICRTEMP(4),PSICRTEMP(5),PSICRTEMP(6),         &
        DCCR(1,1),DCCR(1,2),DCCR(1,3),                  &
        DCCR(2,1),DCCR(2,2),DCCR(2,3),                  &
        DCCR(3,1),DCCR(3,2),DCCR(3,3),                  &
        PSICR(1),PSICR(2),PSICR(3),                     &
        PSICR(4),PSICR(5),PSICR(6),2)

      CALL CHG3(PSICROUT(1),PSICROUT(2),PSICROUT(3),    &
        PSICROUT(4),PSICROUT(5),PSICROUT(6),            &
        DCCRTEMP(1,1),DCCRTEMP(1,2),DCCRTEMP(1,3),      &
        DCCRTEMP(2,1),DCCRTEMP(2,2),DCCRTEMP(2,3),      &
        DCCRTEMP(3,1),DCCRTEMP(3,2),DCCRTEMP(3,3),      &
        PSICRTEMP(1),PSICRTEMP(2),PSICRTEMP(3),         &
        PSICRTEMP(4),PSICRTEMP(5),PSICRTEMP(6),3)

    
    END
    
    
    !*******************************************************************
      SUBROUTINE XSECONCYC_hj ( DEU, EU,    SN,                         &
                 ERC,   SRC,    EXC,   SXC,                             &
                 EBC,   SBC,    EJ,    SJ,    EJJ,                      &
                 ERT,   SRT,    ECP,   SCP,                             &
                 ETP,   STP,    EPC,   SPC,                             &
                            EPC1ST,EPT,    EEN,   SEN,                  &
                 ETMAX, STMAX,  ECCR,  SCCR,                            &
                            E41_51,S41_51,                              &
                            EUOLD, SNOLD,  EPCU,  EPCUS,                &
                            EO,ET, EC,SC,  FC,FT, EBU,   ECU,           &
                            EPCU_,EPCUS_,                               &
                 IUL,   EPPC,   EPPT,  EPEC,                            &
                 CRACK, INCF,   ICIC , C_IC4,  SENS,    VS,             &
                            IVIRGIN,ELIMIT,MM,IROT,ICR,LLLCRACK,        &
                 CYCN, EUDL,COPX,NIT )

! C     *                 DCC11,DCC12,DCC13,
! C     *                       DCC21,DCC22,DCC23,
! C     *                       DCC31,DCC32,DCC33,
! C     *                       DCR11,DCR12,DCR13,
! C     *                       DCR21,DCR22,DCR23,
! C     *                       DCR31,DCR32,DCR33)
! **********************************************************************
! * 引数�E�　DEU           �E�EI)吁E��チE��プ�E主ひずみ増�Eの解
! * 　　EU�E�SN       �E�EI,O)現STEPで計算する主ひずみ、主応力
! *     EEN,SEN     �E�EI)E点の主ひずみ、主応力
! *     EXC,SXC     �E�EI)X点の主ひずみ、主応力
! *     ECP,SCP     �E�EI)C点の主ひずみ、主応力
! *     ETP,STP     �E�EI)T点の主ひずみ、主応力
! *     EPC,SPC     �E�EI)P点の主ひずみ、主応力
! *     EJ,SJ       �E�EI)J点の主ひずみ、主応力
! *     ERC,SRC     �E�EI)R点(除荷曲線上！ETの間）�E再裁荷点)の主ひずみ、主応力
! *     EBC,SBC     �E�EI)B点(除荷曲線上！ETの間）�E再裁荷した時にPC線との交点)の主ひずみ、主応力
! *     EJJ         �E�EI)J点�E�剛性復活点�E��E接線剛性
! *     ERT�E�SRT   �E�EI)R'点の主ひずみ、主応力
! *     EPC1ST      �E�EI)初めて圧縮老E�E履歴でひび割れが入ったとき�E�E�点のひずみを記録�E�EU, LUモチE��で使ぁE��E
! *     EPT,SPT     �E�EI)応力が０になるとき�E(P'点)ひずみ�E�応力
! *     ECR,SCR     �E�EI)コンクリート�E引張強度時！E点�E�主ひずみ、応力
! *     ETMAX,STMAX �E�EI)履歴中経験した最大�E�E'点�E�主ひずみ、応力
! *     ECCR,  SCCR �E�EI)圧縮側から除荷し、�Eび割れた時�Eひずみ、応力
! *     E41_51,S41_51�E�EI)IUL=41OR51の除荷点ひずみ、応力
! *     EUOLD, SNOLD�E�EI)直前スチE��プ�E等価一軸ひずみ,応力
! *     EPCU        �E�EI)コンクリート�E圧壊が生じた後�E収斂点ひずみ
! *     EPCUS       �E�EI)コンクリート�E圧壊が生じた後�E収斂点応力
! *     EO          �E�EI)コンクリート要素の初期剛性
! *     ET          �E�EI/O)当応力点で接線剛性
! *     EC,SC       �E�EI)コンクリート�E三軸応力状態を老E�Eした最大主ひずみ、応力�E�破壊曲面との接点の値�E�E
! *     FC,FT       �E�EI)コンクリート要素の一軸圧縮、引張強度
! *     EBU         �E�EI)チE��ションスチE��フニング効果有効限界ひずみ
! *     ECU         �E�EI)�E�ECU�E�最大圧縮応力時�E等価一軸ひずみ
! *     IUL         �E�EI/O)コンクリート要素の吁E���E点の2つの等価一軸方向�E応力状態を表す指樁E
! *     EPPC        �E�EI)PC線�E剛性�E�EUL=6の剛性�E�E再載荷剛性)
! *     EPPT        �E�EI)除荷剛性�E�引張方向かな�E�？！E
! *     EPEC        �E�EI)EC間�E剛性
! *     ANGL        �E�EI)1スチE��プ前の主応力方向�E角度�E�第�E�！E,3�E�主方向角度�E�E
! *     ANGCR       �E�EO)第�E�！E,3�E��Eび割れ角度
! *     CRACK       �E�EI)ひび割れ状態�E持E��！E
! *     INCF        �E�EO)コンクリート要素の残差力による収斂計算指樁E
! *     ICIC        �E�EI)圧縮応力上�E曲線オプション　0�E�SAENZ式　1�E�FAFITIS-SHAH弁E
! *     SENS        �E�EI)繰返し用の剛性復活点Jを決めるための感度係数�E�現在は使ってぁE��ぁE��E
! *     VS          �E�EI)繰返し履歴用でIUL=4,5へ戻るため�E感度係数
! *     IVIRGIN     �E�EI)除荷が経験したかの持E��E
! *     ELIMIT      �E�EI)
! *     MM          �E�EI)コンクリート要素番号
! **********************************************************************
! **** STRESS - STRAIN RELATIONSHIP FOR CONCRETE                    ****
! **** ORIGINAL PROGRAM = SUBROUTINE SECONC(CYCLIC RULE)            ****
! **** THIS PROGRAM IS MODIFIED    < 1995.12.12 S.TORIZO >          ****
! ****                                                              ****
! **** コンクリート�E吁E���E点でのひずみの増�E量によって除荷また�E   ****
! **** 載荷判定をするプログラム        < 1996.08    J.ZAKKY >       ****
! **********************************************************************
! C
! C *** STRESS - STRAIN RELATIONSHIP FOR CONCRETE
! C     ORIGINAL PROGRAM = SUBROUTINE SECONC(CYCLIC RULE)
! C     THIS PROGRAM IS MODIFIED FOR MONOTONIC LOADING IN ORDER TO
! C     APPLY TO CALCULATE UNBLANCED STRESS < 1991.8.9 K.UCHIDA & A.AMEMIY
! C

    real (kreal) :: DEU,EU,SN,ERC,SRC,EXC,SXC,EBC,SBC,EJ,SJ,EJJ
    real (kreal) :: ERT,SRT,ECP,SCP,ETP,STP,EPC,SPC,EPC1ST,EPT,EEN,SEN
    real (kreal) :: ETMAX,STMAX,ECCR,SCCR,E41_51,S41_51,EUOLD,SNOLD,EPCU,EPCUS 
    real (kreal) :: EO,ET,EC,SC,FC,FT,EBU,ECU,EPPC,EPPT,EPEC,CRACK,SENS,VS
    real (kreal) :: ELIMIT,MODEL42_52,ECR,SCR,EPCU_,EPCUS_
    
    real (kreal) :: DUMMY_SHIRAI,C_IC4,DUMMY_41,DUM41,DUM5,DUM51,DUMMY_51
    real (kreal) :: ec4,dum,dum4,ec5,es,eu4,eu5,saa,sbb,sc4,sc5,sjdammy
    real (kreal) :: eudl, copx, COPXL, EEBA, EEBAL,EECH,EECHL,EETO 
    real (kreal) :: EEUUOOLLDD,EPCL,ESHL,ESHLN,EUDL50,EXC0,EXCL,EXCLL
    real (kreal) :: spcl,ssba,ssbal,ssch,sschangera,sschl,sschll,ssnnoolldd
    real (kreal) :: sxc0,sxcl,sxcll,y,eeuu,sncycnex0,ssnn,ssnnexc0,xdd,xkei
    
    integer :: IUL,INCF,MM,IROT,ICR,LLLCRACK,IVIRGIN, NIT
    integer :: inti,lx,ly,lz,m
    integer :: ICIC(8),IIHL(99)
    integer :: iva,ivb,ivcycn,ivl,melem,ncycn,ncycnoff,nroot
    integer :: nroot1,nroot6,nroot8,nst,ncall,ID_NAN
    real (kreal) :: CYCN(49)
    real (kreal) ::  RRHL(99),CKSXC0(99)
    CHARACTER M10*10
! c      INCLUDE 'com_hj.f'
    real (kreal) :: XHJ(999),C_Infinity(1),S_NaN(99),N_NaN(99)

    CHARACTER*100 M_NaN


    IF(EBU.EQ.0.0) EBU=0.002D0
    IF(ICIC(4)==3.OR.ICIC(4)==4)  EBU=2000.0
    MODEL42_52=1
    VS=1.0
    
    !CHECK
    ivcycn=0
! CH------------------------------
! C *** 問題点 ***
! C *** 1) EEUU > 0.
! C *** 2) EXC,SXCの安宁EX点)
! C *** 3) EPC,SPCの安宁EP点)

    NROOT    = 1 ; NROOT1 = 1 ; NROOT6 = 1 ; NROOT8 = 1
    NCYCNOFF = 0

    SSCH  = CYCN( 4)
    SSCHL = CYCN(20)
    EECH  = CYCN( 5)
    EECHL = CYCN(21)
    SSBA  = CYCN( 6)
    SSBAL = CYCN(16)
    EEBA  = CYCN( 7)
    EEBAL = CYCN(17)
    EETO  = CYCN( 8)
    ESHL  = CYCN( 9)
    ESHLN = CYCN(19)
    IVA   = CYCN(10)
    IVB   = CYCN(11)
    IVL   = CYCN(12)
    SXCL  = CYCN(13)
    SXCLL = CYCN(14)
    EPCL  = CYCN(22)
    SPCL  = CYCN(23)
    COPXL = CYCN(24)
    EXC0  = CYCN(25)
    SXC0  = CYCN(26)
    NCYCN = CYCN(27)
    EXCL  = CYCN(28)
    EXCLL = CYCN(29)

    NST   = CYCN( 1)
    MELEM = CYCN(30) 
    M=MELEM
    LX    = CYCN(31) 
    LY    = CYCN(32) 
    LZ    = CYCN(33) 
    intI  = CYCN(34) 
    NIT   = CYCN(49)
! C
! C SC=CYCN(35) EC=CYCN(36) BETA =CYCN(37)
! C IF(NST>=58 )THEN
! C     WRITE(1820,'(1A20,7I5)')'XSECONCYC_hj',NST,NIT,MELEM,LX,LY,LZ,intI
! C ENDIF
! C---------TAKAHASHI----------------
    IF( EXC0 < ECU ) THEN
        CKSXC0(1)=1.0;CKSXC0(2)=0.85 ;CKSXC0(3)=0.80;CKSXC0(4)=0.77
        CKSXC0(5)=0.75;CKSXC0(6)=0.74;CKSXC0(7)=0.735;CKSXC0(8)=0.73
! C DO I=1,8 ; CKSXC0(I)=CKSXC0(I)-0.1 ; END DO
    ELSE
        CKSXC0(1)=1.0;CKSXC0(2)=0.95 ;CKSXC0(3)=0.90;CKSXC0(4)=0.85
        CKSXC0(5)=0.81;CKSXC0(6)=0.78;CKSXC0(7)=0.76;CKSXC0(8)=0.75
        CKSXC0(9)=0.74;CKSXC0(10)=0.72;CKSXC0(11)=0.69;CKSXC0(12)=0.64
       CKSXC0(13)=0.58;CKSXC0(14)=0.51;CKSXC0(15)=0.44;CKSXC0(16)=0.34
! C DO I=1,8 ; CKSXC0(I)=CKSXC0(I)-0.1 ; END DO
    ENDIF

    IF( NCYCNOFF == 0 ) GOTO 10
    EUDL50 = DABS(EUDL/50.)
    Y=EUDL * DEU
    IF( Y < 0. .AND. DABS(DEU) > EUDL50 )THEN
        IVCYCN=1 ; ELSE ; IVCYCN=0
    ENDIF
    IF(IVCYCN==1) THEN
        ESHL=ESHL+1.
        IF( SC==0. ) THEN;M10='SC=';CALL HJ_C_1(SC,M10);ENDIF
        IF( SNOLD > 0. ) THEN
            SsNnOoLlDd = 0.
            EeUuOoLlDd = EPC
        ELSE
            SsNnOoLlDd = SNOLD
            EeUuOoLlDd = EUOLD
        ENDIF
        SSchangeRA = DABS( (SsNnOoLlDd-SSCH)/SC )
        IF( SSchangeRA > 9./40. ) THEN ! 9./40. ) THEN
! C         IF( DEU > 0. .AND. SsNnOoLlDd > 1./3. * SC ) GOTO 20
            SSCHLL = SSCHL
            SSCHL  = SSCH
            SSCH   = SsNnOoLlDd
            EECHL  = EECH
            EECH   = EeUuOoLlDd
            EXCLL  = EXCL
            EXCL   = EXC
            SXCLL  = SXCL
            SXCL   = SXC
            ESHLN  = ESHLN + 1.
            IF( EUOLD <= EXC * 0.99 ) THEN
                EXC0  = EUOLD
! c             SXC0  = SNOLD
                CALL EPSSIG(EO,ECU,SC,FC,EC,EPCUS_,EPCU_,EXC0,SXC0)
                NCYCN = 1
            ENDIF
            IF( DEU > 0. ) THEN ; IVCYCN=-1 ; ELSE ; IVCYCN=1 ; ENDIF
! C20           CONTINUE
        ELSE
            IVCYCN = 0
        ENDIF
        IF( IVCYCN == -1 ) THEN
            IF( SSCH > 1./3. * SC ) THEN
                IVCYCN = 0
                SSCH   = SSCHL
                SSCHL  = SSCHLL
                ESHLN  = ESHLN - 1.
            ENDIF 
        ENDIF
    ENDIF

    IF( IVCYCN == 0 ) GOTO 10
! ***      EXC,SXCの計算し直ぁE  ***
    IF( IVCYCN == -1 ) THEN
        IF(ESHLN == 1. ) THEN
            EXC = EUOLD
! c         SXC = SNOLD
            CALL EPSSIG(EO,ECU,SC,FC,EC,EPCUS_,EPCU_,EXC,SXC)
        ENDIF
        IF( 0.99*EXC < EUOLD .AND. EUOLD < EECHL ) THEN
            XKEI = ( EUOLD - EECHL ) / (EXC - EECHL )
            IF( 0.99*EXC < EXCLL ) THEN
                EXC = EXCLL + ( EXC - EXCLL ) * XKEI
                CALL EPSSIG(EO,ECU,SC,FC,EC,EPCUS_,EPCU_,EXC,SXC)
! c     IF(nst >= 223 .AND. MELEM == 290) 
! c     +write(1850,*) 'end EPSSIG'
            ENDIF
        END IF
    ELSEIF( IVCYCN == 1 ) THEN
        IF(EXC==0.OR.EXC-EPC==0.) THEN !臨時対策，微小区閁EIVCYCN 問顁E
            EXC = EUOLD
! c         SXC = SNOLD
            CALL EPSSIG(EO,ECU,SC,FC,EC,EPCUS_,EPCU_,EXC,SXC)
        ELSE
            IF( SNOLD > 0. ) THEN
                SSNN = 0.
            ELSE
                SSNN = SNOLD
            ENDIF
            IF ( ESHLN > 0. )   &
                CALL XY_Point(EPCL,SPCL,EXCL,SXCL,SSNN,EEUU,Xdd,1)
! c     IF(nst >= 223 .AND. MELEM == 290) 
! c     +write(1850,*) 'end XY_Point'
           IF(EEUU>=0.)THEN;M10='EEUU=';CALL HJ_C_1(EUOLD,M10);ENDIF !臨時対筁E
            IF(EEUU > 0.0 )  EEUU = EUOLD  ! 臨時対筁E
            IF(EEUU > 0.0 )  EEUU = 0.0     !臨時対筁E
            IF(EXC0 < EPC) THEN
                COPX = ( SXC0 * CKSXC0(NCYCN))/( EXC0-EPC );COPXL=COPX
                IF(COPX > EO ) COPX = EO
            ELSE
                COPX = EO
            ENDIF
            SSNNEXC0  = SSNN + COPX * ( EXC0 - EEUU )
            SNCYCNEX0 = SXC0 * CKSXC0(NCYCN)
            IF( SSNNEXC0 > 1.02*SNCYCNEX0 ) THEN
                NCYCN = NCYCN + 1
                IF( EXC0 < ECU ) THEN
                    IF( NCYCN > 8 ) NCYCN = 8
                ELSE
                    IF( NCYCN > 16 ) NCYCN = 16
                ENDIF
            IF(EXC0 < EPC) THEN
                COPX = ( SXC0 * CKSXC0(NCYCN)) / (EXC0-EPC);COPXL=COPX
                IF(COPX > EO ) COPX = EO
            ELSE
                COPX = EO
            ENDIF
        ENDIF
! c IF(nst >= 223 .AND. MELEM == 290) write(1850,'(14ES10.3)')EEUU,
! c     +SSNN,EO,ECU,SC,FC,EC,EPCUS,EPCU,COPX,SXC0,EXC0,EPC,CKSXC0(NCYCN)
           IF(EXC-EPC==0.)THEN;M10='COPX=';CALL HJ_C_1(COPX,M10);ENDIF
            CALL XPOINT(EEUU,SSNN,EO,ECU,SC,FC,EC,EPCUS_,EPCU_,COPX,    &
                        EXC,SXC)
! c     IF(nst >= 223 .AND. MELEM == 290) 
! c     +write(1850,*) 'end XPOINT'
    CONTINUE
        ENDIF
    ENDIF

    IF( EXC > EXCLL .AND. SXCLL < 1./4.*SC ) THEN
        EXC = EXCLL
        SXC = SXCLL
    END IF
    SEN = SXC
    EEN = EXC
! ***   点の計算し直ぁE  ***
    NCALL = 100000
    CALL ECTPJ_hj( EXC, SXC, ECP, SCP, ETP, STP,        &
               EPC, SPC,  EJ,  SJ,ECCR,SCCR,            &
               EPCU_,EPPC,                              &
                EO,  SC,  FT,  ECU,CRACK,ICIC(1),       &
                EC,ETMAX,ECR,EPEC,IVCYCN ,NCYCNOFF,     &
                NST,NIT,M,LX,LY,LZ,intI,NCALL )
    EPCL=EPC
! ***********************
10  CONTINUE

    CYCN( 4) = SSCH   
    CYCN( 5) = EECH   
    CYCN( 6) = SSBA   
    CYCN(16) = SSBAL  
    CYCN( 7) = EEBA   
    CYCN(17) = EEBAL  
    CYCN( 8) = EETO   
    CYCN( 9) = ESHL  
    CYCN(18) =IVCYCN 
    CYCN(19) = ESHLN  
    CYCN(10) = IVA
    CYCN(11) = IVB 
    CYCN(12) = IVL
    CYCN(20) = SSCHL
    CYCN(21) = EECHL
    CYCN(24) = COPXL
    CYCN(25) = EXC0
    CYCN(26) = SXC0
    CYCN(27) = NCYCN
    CYCN(28) = EXCL
    CYCN(29) = EXCLL
    CYCN(13) = SXCL 
    CYCN(14) = SXCLL
! C
! CH-------------------------------
! C

      ECR=FT/EO
      SCR=FT
    IF(NIT==1) THEN
        IF(IUL.EQ. 1)  GOTO  1
        IF(IUL.EQ. 2)  GOTO  2
        IF(IUL.EQ. 3)  GOTO  3
        IF(IUL.EQ. 31)  GOTO 3
        IF(IUL.EQ. 4)  GOTO  4
        IF(IUL.EQ.41)  GOTO 41
        IF(IUL.EQ.42)  GOTO 42
        IF(IUL.EQ.43)  GOTO 43
        IF(IUL.EQ.46)  GOTO 46
        IF(IUL.EQ. 5)  GOTO  5
        IF(IUL.EQ.51)  GOTO 51
        IF(IUL.EQ.52)  GOTO 52
        IF(IUL.EQ.53)  GOTO 53
        IF(IUL.EQ.56)  GOTO 56
        IF(IUL.EQ. 6)  GOTO  6
        IF(IUL.EQ. 7)  GOTO  7
        IF(IUL.EQ. 8)  GOTO  8
        IF(IUL.EQ. 9)  GOTO  9
    ELSE
        IF( IUL .EQ.  1 )  GOTO  1
        IF( IUL .EQ.  2 )  GOTO  2
        IF( IUL .EQ.  3 )  GOTO  3
        IF( IUL .EQ. 31 )  GOTO  3
        IF( IUL .EQ.  4 )  GOTO  4
        IF( IUL .EQ. 41 )  GOTO  41
        IF( IUL .EQ. 42 )  GOTO  42
        IF( IUL .EQ. 43 )  GOTO  43
        IF( IUL .EQ. 46 )  GOTO  46
        IF( IUL .EQ.  5 )  GOTO  5
        IF( IUL .EQ. 51 )  GOTO  51
        IF( IUL .EQ. 52 )  GOTO  52
        IF( IUL .EQ. 53 )  GOTO  53
        IF( IUL .EQ. 56 )  GOTO  56
        IF( IUL .EQ.  6 )  GOTO  6
        IF( IUL .EQ.  7 )  GOTO  7
        IF( IUL .EQ.  8 )  GOTO  8
        IF( IUL .EQ.  9 )  GOTO  9
    ENDIF
! ***************
! **** IUL=1 ****
! ***************
    1 CALL ROOT_CHECK(   1, 57)
      IF(EU)       1001,1002,1002
 1001 CALL ROOT_CHECK(1001, 57)
      IF(DEU)      1010,1010,1003
 1003 CALL ROOT_CHECK(1003, 57)
      IF(IVIRGIN)   901,1010,1004
 1010 CALL ROOT_CHECK(1010, 57)
! C       IF(EU.LT.EPCU_) WRITE(*,*)EPCU
      IF(EU.LT.EPCU_) GOTO 9000
      IF(EU.LE.EC)   GOTO 8000
      IF(EC.LT.EU)   GOTO 1000
 1002 CALL ROOT_CHECK(1002, 57)
    IF(CRACK.EQ.4.0) THEN
        IF(EU.GE.EBU) GOTO 4090
        CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
        CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
        IF(DUMMY_SHIRAI.LT.DUMMY_41) THEN
            GOTO 4000
        ELSE
            GOTO 4100
        ENDIF
    ENDIF
    
      IF(ECR.LT.EU)  GOTO 4002
      IF(EU.LE.ECR)  GOTO 1040
! ****IUL=1からの除荷の準備  IUL=2を経由
! C 1004 CALL ROOT_CHECK(1004, 57)

1004    CONTINUE
      IVIRGIN=1
      EEN=EUOLD ! EUCOLD(*,*,*,*,*,*)
      SEN=SNOLD ! PSIOLD(*,*,*,*,*,*)
! CH--------------------------------------------------
    IF(NROOT==0) GOTO 10010
    IF(NROOT1==0) GOTO 10010
    NCALL = 11111
      CALL ECTPJ_hj( EEN, SEN, ECP, SCP, ETP, STP,      &
                 EPC, SPC,  EJ,  SJ,ECCR,SCCR,          &
                EPCU_,EPPC,                             &
                  EO,  SC,  FT,  ECU,CRACK,ICIC(1),     &
                  EC,ETMAX,ECR,EPEC,IVCYCN ,NCYCNOFF,   &
                NST,NIT,M,LX,LY,LZ,intI,NCALL )
      ET=EPEC
      IF(EPPC.EQ.EO)THEN
        EXC=EEN
        SXC=SEN
      ELSE
        CALL CROSS1(EPPC,ECP,SCP,EXC,SXC,EO,SC,EPC,     &
                  EC,IUL,EPCU_,EPCUS_,ICIC(1))
      ENDIF
10010   CONTINUE

        
! CH--------------------------------------------------
      IF(EPPC.EQ.EO)THEN
        IF(EU.LT.EPC) GOTO 6000                        !もし、現ひずみが残留ひずみより小さければIUL=3決定、これより下�E引張域、E
        IF(CRACK.EQ.0.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          IF(ECCR.LE.EU)  GOTO 5000
          IF(EU.LT.ECCR)  GOTO 3010
        ENDIF
        IF(CRACK.EQ.4.0)THEN
          IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
        ENDIF
        IF(CRACK.EQ.5.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
            GOTO 5000
          ELSE
            GOTO 5100
          ENDIF
        ENDIF
      ELSE
        IF(EU.LE.ETP)  GOTO 2000
       IF(EU.LT.EPC)  GOTO 3000
        IF(CRACK.EQ.0.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          IF(ECCR.LE.EU)  GOTO 5000
          IF(EU.LT.ECCR)  GOTO 3010
        ENDIF
        IF(CRACK.EQ.4.0)THEN
          IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
        ENDIF
        IF(CRACK.EQ.5.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
            GOTO 5000
          ELSE
            GOTO 5100
          ENDIF
        ENDIF
      ENDIF
! ****応力上�E曲線丁E
 1000 CALL ROOT_CHECK(1000, 57)
      IUL=1
      INCF=1
    CALL CompressionZone(ICIC(1),EO,ET,EC,SC,EU,SN)
      IF(ET.LT.EO/100.0) ET=EO/100.0
! C     IF(SN.LT.SC) GOTO 903
      IF(IVIRGIN.EQ.0)THEN
        IF(SN.LT.SC/3.)IVIRGIN=1
      ENDIF
! CH---エラーチェチE��---
    M_NaN='☁E1000: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ****引張弾性埁E
 1040 CALL ROOT_CHECK(1040, 57)
      IUL=1
      INCF=1
      ET=EO
      SN=EO*EU
! CH---エラーチェチE��---
    M_NaN='☁E1040: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ***************
! **** IUL=2 ****
! ***************
    2 CALL ROOT_CHECK(   2, 57)
      IF(DEU)2001,2001,2002
 2001 CALL ROOT_CHECK(2001, 57)
    IF(EEN.GT.EC)THEN
        IF(EU.LT.EPCU_) GOTO 9000
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EEN)  GOTO 1000
        IF(EEN.LT.EU)  GOTO 2000
    ELSE
        IF(EU.LT.EPCU_) GOTO 9000
        IF(EU.LE.EEN)  GOTO 8000
        IF(EEN.LT.EU)  GOTO 2000
    ENDIF
 2002 CALL ROOT_CHECK(2002, 57)
      IF(EU.LE.ETP)  GOTO 2000
      IF(EU.LE.EPC)  GOTO 3000
      IF(CRACK.EQ.0.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        IF(ECCR.LT.EU)  GOTO 5000
        IF(EU.LE.ECCR)  GOTO 3010
      ENDIF
      IF(CRACK.EQ.4.0)THEN
        IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
      ENDIF
      IF(CRACK.EQ.5.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
        CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
        IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
          GOTO 5000
        ELSE
          GOTO 5100
        ENDIF
      ENDIF
 2000 CALL ROOT_CHECK(2000, 57)
      IUL=2
      INCF=1
      ET=EPEC
      SN=SEN+(EU-EEN)*EPEC
! CH---エラーチェチE��---
    M_NaN='☁E2000: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
 2100 CONTINUE
      IUL=2
      INCF=1
      ET=(SN-SPC)/(EU-EPC)
      SN=SN+DEU*ET
! CH---エラーチェチE��---
    M_NaN='☁E2100: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ***************
! **** IUL=3 ****
! ***************
    3 CALL ROOT_CHECK(   3, 57)

    IF(EPEC.EQ.EPPC) GOTO 6

      IF(EPC.LT.EUOLD) GOTO 3003                                               !直線部にジャンチE
      IF(DEU) 3001,3002,3002
 3001 CALL ROOT_CHECK(3001, 57)
      ERC=EUOLD
      SRC=SNOLD
        
      EBC=(EPEC*ERC-SRC-EPPC*EPC)/(EPEC-EPPC)
      SBC=EPPC*(EBC-EPC)
      IF(EU.LT.EPCU_) GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)  GOTO 8000
      ELSE
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EXC)  GOTO 1000
      ENDIF
      IF(EU.LT.EBC)  GOTO 6000
      IF(EBC.LT.EU)  GOTO 7000
 3002 CALL ROOT_CHECK(3002, 57)
      IF(EU.LE.EPC)  GOTO 3000
      IF(CRACK.EQ.0.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        IF(ECCR.LT.EU)  GOTO 5000
        IF(EU.LE.ECCR)  GOTO 3010
      ENDIF
      IF(CRACK.EQ.4.0)THEN
        IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
      ENDIF
      IF(CRACK.EQ.5.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
        CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
        IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
          GOTO 5000
        ELSE
          GOTO 5100
        ENDIF
      ENDIF
 3003 CALL ROOT_CHECK(3003, 57)
      IF(DEU) 3004,3005,3005
 3004 CALL ROOT_CHECK(3004, 57)
      IF(EU.LT.EPCU_) GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)  GOTO 8000
      ELSE
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EXC)  GOTO 1000
      ENDIF
      IF(EU.LT.EPC)  GOTO 6000
      IF(EPC.LE.EU)  GOTO 3010
 3005 CALL ROOT_CHECK(3005, 57)
      IF(EBU.LT.EU)   GOTO 5090
      IF(ECCR.LT.EU)  GOTO 5000
      IF(EU.LE.ECCR)  GOTO 3010
! ****IUL=3 曲線丁E
 3000 CALL ROOT_CHECK(3000, 57)
      IUL=3
      INCF=1
      IF(EPPC.EQ.EO)THEN
        ET=EO
        SN=(EU-EPC)*EO
      ELSE
        CALL CURVE3(   EU,   SN,    &
                     EPC,           &
                    ETP,  STP,  &
                   EPEC,   ET )
      ENDIF
! CH---エラーチェチE��---
    M_NaN='☁E3000: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    ! S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ****IUL=3 直線丁E
 3010 CALL ROOT_CHECK(3010, 57)
      IUL=31
      INCF=1
      ET=EO
      SN=EO*(EU-EPC)
! CH---エラーチェチE��---
    M_NaN='☁E3010: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ***********************************************
! **** IUL=4 (圧縮域を経験してぁE��いT/S領域)****
! ***********************************************
    4 CALL ROOT_CHECK(   4, 57)
      IF(DEU) 4001,4004,4004
 4001 CALL ROOT_CHECK(4001, 57)                                                !除荷発甁E
      IF(EBU.LT.EU) GOTO 4090                                                  !EBU:T/S限界ひずみ
      IF(EBU.LE.EUOLD)THEN
        ETMAX=EBU
        STMAX=0.0
      ELSE
        ETMAX=EUOLD
        STMAX=SNOLD
      ENDIF
! ***************************************
      IF(MODEL42_52.EQ.0)THEN                                                  !鳥山垁E
        EPPT=EO*ECR/ETMAX                                                      !出所不�E
      ENDIF
      IF(MODEL42_52.EQ.1.OR.    &                                                  !直線モチE��
            MODEL42_52.EQ.4.OR.    &                                               !鳥山�E�E��線　混合モチE��
            MODEL42_52.EQ.5.OR.    &                                              !改良鳥山モチE��
            MODEL42_52.EQ.7)THEN                                                  !櫻井モチE��用
        EPPT=(STMAX-SJ)/(ETMAX-EJ)                                             !圧縮域での除荷を経験してぁE��ぁE�Eで、E��荷剛性は引張強度をもとに決める、E
      ENDIF
      IF(MODEL42_52.EQ.2)THEN                                                  !大乁E��論文
        EPPT=1.5*EO*ECR/ETMAX
      ENDIF
      IF(MODEL42_52.EQ.3)THEN                                                  !初期剛性モチE��
        EPPT=EO
      ENDIF
! C      IF(MODEL42_52.EQ.7)THEN                                                  !長沼モチE��
! C        EPPT=EO*ECR/ETMAX
! C      ENDIF
! ***************************************
! C      IF(MODEL42_52.EQ.7.AND.EPT.LT.5.0E-5)THEN
! C        EPPT=EO
! C      ENDIF
        EPT=ETMAX-STMAX/EPPT
      IF(EU.LT.EPCU_) GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)  GOTO 8000
        IF(EU.LE.EJ)   GOTO 6000
        IF(EJ.LT.EU)   GOTO 4210
      ELSE
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EJ)THEN
          IF(EEN.EQ.0.0)THEN
            GOTO 1000
          ELSE
            IF(EU.LE.EXC) GOTO 1000
            IF(EXC.LT.EU) GOTO 6000
          ENDIF
        ELSE
          GOTO 4210
        ENDIF
      ENDIF
 4004 CALL ROOT_CHECK(4004, 57)
      IF(EBU.LE.EU) GOTO 4090
      IF(EU.LE.EBU) GOTO 4000
! ****TENSION STIFFENING 進衁E
! ****IUL=1>>4  ひび発甁ETENSION STIFFENING
 4002 CALL ROOT_CHECK(4002, 57)
! C      IF(ANGCR.LT.900.)GOTO 4003
! C      ANGCR=ANGL
 4003 CALL ROOT_CHECK(4003, 57)
      IVIRGIN=1
      ET=EO/100.
      IF(EEN.EQ.0.0)THEN
! *       SJ=-FT*SENS                                                            !IUL=2を未経験（感度解析部刁E��E
! C        SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)*15                                      !FROM KOBAYASHI(MODIFIED BY SAKURAI)
        SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)                                      !FROM TIDE
      CALL CompressionZone_REVERS(ICIC(1),FC,SJ,    &
                    EO,EJJ,ECU,SC,EJ,SJDAMMY)
      ELSE
! *       SJ=SEN*0.1                                                             !IUL=2を経騁E
! C        SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)*15                                      !FROM KOBAYASHI(MODIFIED BY SAKURAI)
        SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)                                     !FROM TIDE
        EJ=SJ/EPPC+EPC
      ENDIF
! ****収斂回数�E�回目以陁E
 4005 CALL ROOT_CHECK(4005, 57)
      IF(EU.LT.EBU) GOTO 4000
      IF(EBU.LE.EU) GOTO 4090
! ****白井式曲線丁E
 4000 CALL ROOT_CHECK(4000, 57)
      IUL=4
      INCF=3
      CRACK=4.0
      ET=EO/100.
      CALL SHIRAI(EU,SN,ECR,SCR,EBU,ICIC(4),C_IC4)
    ICR = 1
! CH---エラーチェチE��---
    M_NaN='☁E4000: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ****さらに允E
 4090 CALL ROOT_CHECK(4090, 57)
      IUL=4
      INCF=3
      CRACK=4.0
      ET=EO/100.
      SN=0.0
    ICR = 1
      RETURN
! ****************
! **** IUL=41 ****
! ****************
   41 CALL ROOT_CHECK(  41, 57)
      IF(DEU)4101,4102,4102
 4101 CALL ROOT_CHECK(4101, 57)
      IF(EU.LT.EPCU_)    GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)   GOTO 8000
      ELSE
        IF(EU.LE.EC)    GOTO 8000
        IF(EU.LE.EXC)   GOTO 1000
      ENDIF
      IF(EU.LE.EJ)      GOTO 6000
      IF(EJ.LT.EU)      GOTO 4303
 4102 CALL ROOT_CHECK(4102, 57)
      IF(EBU.LT.EU)    GOTO 4090
      CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
      CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
      IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
        GOTO 4000
      ELSE
        GOTO 4100
      ENDIF
 4100 CALL ROOT_CHECK(4100, 57)
      IUL=41
      INCF=1
      CALL STRAIGHT41_51(EU,SN,ET,ETMAX,STMAX,EPC,VS)
! CH---エラーチェチE��---
    M_NaN='☁E4100: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! **************************************************************************************
! **** IUL=42 (圧縮の影響を受けてぁE��ぁE��ンション・スチE��フニング領域からの除荷曲緁E****
! **************************************************************************************
   42 CALL ROOT_CHECK(  42, 57)
      IF(DEU) 4201,4201,4202

 4201 CALL ROOT_CHECK(4201, 57)                                                !そ�Eまま進行、E
      IF( EU < EPCU_ )   GOTO 9000
      IF( EXC < EC ) THEN
! C WRITE(1850,'(3ES11.3,7I5)') EU,EXC,EC,NST,NIT,M,LX,LY,LZ,INTI
        IF( EU <= EXC )    GOTO 8000
        IF( EU <= EJ  )    GOTO 6000
        IF( EJ <  EU  )    GOTO 4200
      ELSE
! C WRITE(1850,'(3ES11.3,7I5)') EU,EXC,EC,NST,NIT,M,LX,LY,LZ,INTI
        IF( EU <= EC )    GOTO 8000
        IF( EU <= EJ )THEN
            IF(EEN.EQ.0.0)THEN
                IVIRGIN=0
! C             IVIRGIN=1
                GOTO 1000
            ELSE
                IF(EU.LE.EXC) GOTO 1000
                IF(EXC.LT.EU) GOTO 6000
            ENDIF
        ELSE
            GOTO 4200
        ENDIF
      ENDIF

 4202 CALL ROOT_CHECK(4202, 57)                                                !停滞（増�Eひずみ�E�E�E�もしくは、�E載荷、E
      IF(EBU.LT.EU)   GOTO 4090
      IF(ETMAX.LT.EU) GOTO 4000
      IF(EU.LE.ETMAX) GOTO 4610
 4210 CALL ROOT_CHECK(4210, 57)
      ETMAX=EUOLD
      STMAX=SNOLD
 4200 CALL ROOT_CHECK(4200, 57)
      IUL=42
      INCF=1
! ***********************************************************
      IF(MODEL42_52.EQ.0.OR.   &                                                !鳥山垁E
        MODEL42_52.EQ.5)THEN                                                  !改良鳥山垁E
        IF(EPT.LT.EU)THEN
          EU4=ETMAX-EU
          EC4=ETMAX-EPT
          DUM4=1.+(EO*EC4/STMAX-2.)*EU4/EC4+(EU4/EC4)**2
          ET=EO*(1.-(EU4/EC4)**2)/DUM4**2
          SN=STMAX-EO*EU4/DUM4
        ELSE
          EU4=EU-EJ
          EC4=EPT-EJ
          SC4=-SJ
          DUM41=EU4/EC4
          DUM4=EO*EC4/SC4
          ET=EO*(1.-DUM41)**(DUM4-1.)
          SN=SC4*(1.-(1.-DUM41)**DUM4)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.1)THEN                                                  !直線モチE��
        ET=EPPT
        SN=ET*(EU-EJ)+SJ
      ENDIF
      IF(MODEL42_52.EQ.2)THEN                                                  !大乁E��論文
        IF(EPT.LT.EU)THEN
          ET=EPPT
          SN=ET*(EU-EPT)
        ELSE
          ET=-SJ*(EPT-EJ)
          SN=ET*(EU-EJ)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.3)THEN                                                  !初期剛性モチE��
        IF(EPT.LT.EU)THEN
          ET=EPPT
          SN=ET*(EU-EPT)
        ELSE
          ET=-SJ*(EPT-EJ)
          SN=ET*(EU-EJ)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.4)THEN                                                  !曲線！E��線　混合型
        IF(EPT.LT.EU)THEN
          EU4=ETMAX-EU
          EC4=ETMAX-EPT
          DUM4=1.+(EO*EC4/STMAX-2.)*EU4/EC4+(EU4/EC4)**2
          ET=EO*(1.-(EU4/EC4)**2)/DUM4**2
          SN=STMAX-EO*EU4/DUM4
        ELSE
          ET=EPPT
          SN=ET*(EU-EJ)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.7)THEN                                                  !長沼モチE��

       IF(SJ.GE.0.)THEN
         WRITE(6,*)'SJの応力が決まってません'
         STOP
       ENDIF

        IF(EPT.LT.EU)THEN
         ET=EPPT
         SN=ET*(EU-EJ)+SJ
        ELSE
          IF(EPT.LT.1.5*ECR)THEN
         ET=EPPT
         SN=ET*(EU-EJ)+SJ
          ELSE
           IF(EEN.LT.0.0) THEN
            CALL NN42_52( SJ, EJ, EPPC, EPT, EU, SN, ET, MM)
           ELSE
            CALL NN42_52( SJ, EJ,  EJJ, EPT, EU, SN, ET, MM)
           ENDIF
          ENDIF
        ENDIF

! *      WRITE(1114,*)'来てますよ',SJ,EJ

      ENDIF
! CH---エラーチェチE��---
    M_NaN='☁E4200: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ****************
! **** IUL=43 ****
! ****************
   43 CALL ROOT_CHECK(  43, 57)
      IF(DEU)4301,4302,4302
 4301 CALL ROOT_CHECK(4301, 57)
      IF(EU.LT.EPCU_)    GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)   GOTO 8000
      ELSE
        IF(EU.LE.EC)    GOTO 8000
        IF(EU.LE.EXC)   GOTO 1000
      ENDIF
      IF(EU.LE.EJ)      GOTO 6000
      IF(EJ.LT.EU)      GOTO 4300
 4302 CALL ROOT_CHECK(4302, 57)
      IF(EBU.LT.EU)     GOTO 4090
      CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
      CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
      IF(DUMMY_SHIRAI.LT.DUMMY_41)GOTO 4000
      IF(E41_51.LT.EU)  GOTO 4100
      IF(EU.LE.E41_51)  GOTO 4300
 4303 CALL ROOT_CHECK(4303, 57)
      E41_51=EUOLD
      S41_51=SNOLD
 4300 CALL ROOT_CHECK(4300, 57)
      IUL=43
      INCF=1
! *      IF(MODEL42_52.EQ.7)THEN
! *       IF(EU.GT.0.001)THEN
! *         GOTO4001
! *         CALL NN43_53( EO, SJ, EJ, EPPC, S41_51, E41_51, EU, SN, ET, MM)
! *       ENDIF
! *      ELSE
       ET=(S41_51-SJ)/(E41_51-EJ)
       SN=ET*(EU-EJ)+SJ
! *      ENDIF
! CH---エラーチェチE��---
    M_NaN='☁E4300: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ****************
! **** IUL=46 ****
! ****************
   46 CALL ROOT_CHECK(  46, 57)
      IF(DEU)4601,4602,4602
 4601 CALL ROOT_CHECK(4601, 57)
      IF(EU.LT.EPCU_)  GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC) GOTO 8000
        IF(EU.LE.EJ)  GOTO 6000
      ELSE
        IF(EU.LE.EC)    GOTO 8000
        IF(EU.LE.EJ)THEN
          IF(EEN.EQ.0.0)THEN
            GOTO 1000
          ELSE
            IF(EU.LE.EXC) GOTO 1000
            IF(EXC.LT.EU) GOTO 6000
          ENDIF
        ENDIF
      ENDIF
      IF(EU.LE.ERT)   GOTO 4200
      IF(ERT.LT.EU)   GOTO 4600
 4602 CALL ROOT_CHECK(4602, 57)
      IF(EBU.LT.EU)   GOTO 4090
      IF(ETMAX.LT.EU) GOTO 4000
      IF(EU.LE.ETMAX) GOTO 4600
 4610 CALL ROOT_CHECK(4610, 57)
      ERT=EUOLD
      SRT=SNOLD
 4600 CALL ROOT_CHECK(4600, 57)
      IUL=46
      INCF=1
    IF(MODEL42_52.EQ.0.OR.  &                                                   !鳥山垁E
        MODEL42_52.EQ.5)THEN                                                  !改良鳥山垁E
        IF(EPT.LT.EU)THEN
          EU4=ETMAX-EU
          EC4=ETMAX-EPT
          DUM4=1.+(EO*EC4/STMAX-2.)*EU4/EC4+(EU4/EC4)**2
          ET=EO*(1.-(EU4/EC4)**2)/DUM4**2
          SN=STMAX-EO*EU4/DUM4
        ELSE
          EU4=EU-EJ
          EC4=EPT-EJ
          SC4=-SJ
          DUM41=EU4/EC4
          DUM4=EO*EC4/SC4
          ET=EO*(1.-DUM41)**(DUM4-1.)
          SN=SC4*(1.-(1.-DUM41)**DUM4)+SJ
        ENDIF
      ELSE IF(MODEL42_52.EQ.7)THEN
        IF(EU.GE.EPT)THEN
         ET=(STMAX-SRT)/(ETMAX-ERT)
         SN=SRT+(EU-ERT)*ET
        ELSE
          IF(EPT.LT.1.5*ECR)THEN
         ET=(STMAX-SRT)/(ETMAX-ERT)
         SN=SRT+(EU-ERT)*ET
          ELSE
           IF(EEN.LT.0.0) THEN
            CALL NN42_52( SJ, EJ, EPPC, EPT, EU, SN, ET, MM)
           ELSE
            CALL NN42_52( SJ, EJ,  EJJ, EPT, EU, SN, ET, MM)
           ENDIF
          ENDIF
        ENDIF
      ELSE
       ET=(STMAX-SRT)/(ETMAX-ERT)
       SN=SRT+(EU-ERT)*ET
      ENDIF
! CH---エラーチェチE��---
    M_NaN='☁E4600: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ****************
! **** IUL=5  ****
! ****************
    5 CALL ROOT_CHECK(   5, 57)
      IF(DEU) 5001,5003,5003
 5001 CALL ROOT_CHECK(5001, 57)
      IF(EBU.LE.EU)  GOTO 5090
      IF(EBU.LT.EUOLD)THEN
        ETMAX=EBU
        STMAX=0.0
      ELSE
        ETMAX=EUOLD
        STMAX=SNOLD
      ENDIF
! ***************************************
      IF(MODEL42_52.EQ.0)THEN                                                  !鳥山垁E
 5099   CALL ROOT_CHECK(5099, 57)
        EPPT=EO*ECCR/ETMAX                                                     !出所不�E
      ENDIF
      IF(MODEL42_52.EQ.1.OR.   &                                                !直線モチE��
        MODEL42_52.EQ.4.OR.   &                                                !鳥山�E�E��線　混合モチE��
        MODEL42_52.EQ.5.OR.    &                                               !改良鳥山モチE��
        MODEL42_52.EQ.7)THEN                                                  !櫻井モチE��用
 5098   CALL ROOT_CHECK(5098, 57)
        EPPT=(STMAX-SJ)/(ETMAX-EJ)
      ENDIF
      IF(MODEL42_52.EQ.2)THEN                                                  !大乁E��論文
 5097   CALL ROOT_CHECK(5097, 57)
        EPPT=1.5*EO*ECCR/ETMAX
      ENDIF
      IF(MODEL42_52.EQ.3)THEN                                                  !初期剛性モチE��
 5096   CALL ROOT_CHECK(5096, 57)
        EPPT=EO
      ENDIF
      IF(MODEL42_52.EQ.4)THEN                                                  !鳥山�E�E��線　混合モチE��
 5095   CALL ROOT_CHECK(5095, 57)
        EPPT=(STMAX-SJ)/(ETMAX-EJ)
      ENDIF
! C      IF(MODEL42_52.EQ.7)THEN                                                  !長沼モチE��
! C 5094   CALL ROOT_CHECK(5094, 57)
! C        EPPT=EO*ECCR/ETMAX
! C      ENDIF
! ***************************************
! C      IF(MODEL42_52.EQ.7.AND.EPT.LT.(EPT+5.0E-5))THEN
! C        EPPT=EO
! C      ENDIF
        EPT=ETMAX-STMAX/EPPT
      IF(EU.LT.EPCU_) GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)  GOTO 8000
      ELSE
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EXC)  GOTO 1000
      ENDIF
      IF(EU.LE.EJ)   GOTO 6000
      IF(EJ.LT.EU)   GOTO 5210
 5003 CALL ROOT_CHECK(5003, 57)
      IF(EBU.LT.EU)  GOTO 5090
      IF(EU.LE.EBU)  GOTO 5000
! *****収斂�E�回目以陁E
 5004 CALL ROOT_CHECK(5003, 57)
      IF(EBU.LT.EU)  GOTO 5090
      IF(EU.LE.EBU)  GOTO 5002
 5000 CALL ROOT_CHECK(5000, 57)
! C      IF(ANGCR.LT.900.)THEN                                                    !ひずみ軟化領域
! C        EPC1ST=EPC                                                             !初めて圧縮老E�E履歴でひび割れが入ったとき�E�E�点のひずみを記録�E�EU, LUモチE��で使ぁE��E
! C        GOTO 5002
! C      ENDIF
! C      ANGCR=ANGL
 5002 CALL ROOT_CHECK(5002, 57)
      IUL=5
      INCF=3
      CRACK=5.0
      ET=EO/100.
      CALL SHIRAI(EU,SN,ECCR,SCCR,EBU,ICIC(4),C_IC4)
    ICR = 1
! CH---エラーチェチE��---
    M_NaN='☁E5002: [ EU,SN,ECCR,SCCR,EBU ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = EU ; S_NaN(2) = SN; S_NaN(3) = ECCR; S_NaN(4) = SCCR
    S_NaN(5) =EBU 
    ! CALL CCHJ(N_NaN, 7, S_NaN, 5 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
 5090 CALL ROOT_CHECK(5090, 57)
      IUL=5                                                                    !チE��ションスチE��フニングを趁E��た領域
      INCF=3
      CRACK=5.0
      ET=EO/100.
      SN=0.0
    ICR = 1
! CH---エラーチェチE��---
    M_NaN='☁E5090: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ****************
! **** IUL=51 ****
! ****************
   51 CALL ROOT_CHECK(  51, 57)
      IF(DEU)5101,5102,5102
 5101 CALL ROOT_CHECK(5101, 57)
      IF(EU.LT.EPCU_)  GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC) GOTO 8000
      ELSE
        IF(EU.LE.EC)  GOTO 8000
        IF(EU.LE.EXC) GOTO 1000
      ENDIF
      IF(EU.LE.EJ)    GOTO 6000
      IF(EJ.LT.EU)    GOTO 5303
 5102 CALL ROOT_CHECK(5102, 57)
      IF(EBU.LT.EU)   GOTO 5090
      CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
      CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
      IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
        GOTO 5000
      ELSE
        GOTO 5100
      ENDIF
 5100 CALL ROOT_CHECK(5100, 57)
      IUL=51
      INCF=1
      CALL STRAIGHT41_51(EU,SN,ET,ETMAX,STMAX,EPC,VS)
! CH---エラーチェチE��---
    M_NaN='☁E5100: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ****************
! **** IUL=52 ****
! ****************
   52 CALL ROOT_CHECK(  52, 57)
      IF(DEU)5201,5202,5202
 5201 CALL ROOT_CHECK(5201, 57)
      IF(EU.LT.EPCU_)  GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)  GOTO 8000
      ELSE
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EXC)  GOTO 1000
      ENDIF
      IF(EU.LE.EJ)    GOTO 6000
      IF(EJ.LT.EU)    GOTO 5200
 5202 CALL ROOT_CHECK(5202, 57)
      IF(EBU.LT.EU)   GOTO 5090
      IF(ETMAX.LT.EU) GOTO 5000
      IF(EU.LE.ETMAX) GOTO 5610
 5210 CALL ROOT_CHECK(5210, 57)
      ETMAX=EUOLD
      STMAX=SNOLD
 5200 CALL ROOT_CHECK(5200, 57)
      IUL=52
      INCF=1
! ***********************************************************
      IF(MODEL42_52.EQ.0.OR.    &                                                !鳥山垁E
        MODEL42_52.EQ.5)THEN                                                  !改良鳥山垁E
        IF(EPT.LT.EU)THEN
          EU5=ETMAX-EU
          EC5=ETMAX-EPT
          DUM5=1.+(EO*EC5/STMAX-2.)*EU5/EC5+(EU5/EC5)**2
          ET=EO*(1.-(EU5/EC5)**2)/DUM5**2
          SN=STMAX-EO*EU5/DUM5
        ELSE
          EU5=EU-EJ
          EC5=EPT-EJ
          SC5=-SJ
          DUM51=EU5/EC5
          DUM5=EO*EC5/SC5
          ET=EO*(1.-DUM51)**(DUM5-1.)
          SN=SC5*(1.-(1.-DUM51)**DUM5)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.1)THEN                                                  !直線モチE��
! C        EPPT=(STMAX-SJ)/(ETMAX-EJ)
        ET=EPPT
        SN=ET*(EU-EJ)+SJ
      ENDIF
      IF(MODEL42_52.EQ.2)THEN                                                  !大乁E��論文
        IF(EPT.LT.EU)THEN
          ET=EPPT
          SN=ET*(EU-EPT)
        ELSE
          ET=-SJ*(EPT-EJ)
          SN=ET*(EU-EJ)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.3)THEN                                                  !初期剛性モチE��
        IF(EPT.LT.EU)THEN
          ET=EPPT
          SN=ET*(EU-EPT)
        ELSE
          ET=-SJ*(EPT-EJ)
          SN=ET*(EU-EJ)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.4)THEN                                                  !鳥山�E�E��線　混合型
        IF(EPT.LT.EU)THEN
          EU5=ETMAX-EU
          EC5=ETMAX-EPT
          DUM5=1.+(EO*EC5/STMAX-2.)*EU5/EC5+(EU5/EC5)**2
          ET=EO*(1.-(EU5/EC5)**2)/DUM5**2
          SN=STMAX-EO*EU5/DUM5
        ELSE
          ET=EPPT
          SN=ET*(EU-EJ)+SJ
        ENDIF
      ENDIF
      IF(MODEL42_52.EQ.7)THEN                                                  !長沼モチE��
        IF(EPT.LT.EU)THEN
         ET=EPPT
         SN=ET*(EU-EJ)+SJ
        ELSE
          IF(EPT.LT.1.0E-3)THEN
         ET=EPPT
         SN=ET*(EU-EJ)+SJ
          ELSE
           CALL NN42_52( SJ, EJ, EPPC, EPT, EU, SN, ET, MM)
         ENDIF
        ENDIF
      ENDIF
! CH---エラーチェチE��---
    M_NaN='☁E5200: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ****************
! **** IUL=53 ****
! ****************
   53 CALL ROOT_CHECK(  53, 57)
      IF(DEU)5301,5302,5302
 5301 CALL ROOT_CHECK(5301, 57)
      IF(EU.LT.EPCU_)    GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)   GOTO 8000
      ELSE
        IF(EU.LE.EC)    GOTO 8000
        IF(EU.LE.EXC)   GOTO 1000
      ENDIF
      IF(EU.LE.EJ)      GOTO 6000
      IF(EJ.LT.EU)      GOTO 5300
 5302 CALL ROOT_CHECK(5302, 57)
      IF(EBU.LT.EU)     GOTO 5090
      CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
      CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
      IF(DUMMY_SHIRAI.LT.DUMMY_51)GOTO 5000
      IF(E41_51.LT.EU)  GOTO 5100
      IF(EU.LE.E41_51)  GOTO 5300
 5303 CALL ROOT_CHECK(5303, 57)
      E41_51=EUOLD
      S41_51=SNOLD
 5300 CALL ROOT_CHECK(5300, 57)
      IUL=53
      INCF=1
! C      IF(MODEL42_52.EQ.7)THEN
! C       IF(EU.GT.0.001)THEN
! C        CALL NN43_53( EO, SJ, EJ, EPPC, S41_51, E41_51, EU, SN, ET, MM)
! C       ENDIF
! C      ELSE
       ET=(S41_51-SJ)/(E41_51-EJ)
       SN=ET*(EU-EJ)+SJ
! C      ENDIF
! CH---エラーチェチE��---
    M_NaN='☁E5300: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ****************
! **** IUL=56 ****
! ****************
   56 CALL ROOT_CHECK(  56, 57)
      IF(DEU)5601,5602,5602
 5601 CALL ROOT_CHECK(5601, 57)
      IF(EU.LT.EPCU_)  GOTO 9000
      IF(EXC.LT.EU)THEN
        IF(EU.LE.EXC)   GOTO 8000
      ELSE
        IF(EU.LE.EC)    GOTO 8000
        IF(EU.LE.EXC)   GOTO 1000
      ENDIF
      IF(EU.LE.EJ)    GOTO 6000
      IF(EU.LE.ERT)   GOTO 5200
      IF(ERT.LT.EU)   GOTO 5600
 5602 CALL ROOT_CHECK(5602, 57)
      IF(EBU.LT.EU)   GOTO 5090
      IF(ETMAX.LT.EU) GOTO 5000
      IF(EU.LE.ETMAX) GOTO 5600
 5610 CALL ROOT_CHECK(5610, 57)
      ERT=EUOLD
      SRT=SNOLD
 5600 CALL ROOT_CHECK(5600, 57)
      IUL=56
      INCF=1
      IF(MODEL42_52.EQ.7)THEN
        IF(EU.GE.EPT)THEN
         ET=(STMAX-SRT)/(ETMAX-ERT)
         SN=SRT+(EU-ERT)*ET
        ELSE
          IF(EPT.LT.1.0E-3)THEN
         ET=(STMAX-SRT)/(ETMAX-ERT)
         SN=SRT+(EU-ERT)*ET
          ELSE
            CALL NN42_52( SJ, EJ, EPPC, EPT, EU, SN, ET, MM)
         ENDIF
        ENDIF
      ELSE
       ET=(STMAX-SRT)/(ETMAX-ERT)
       SN=SRT+(EU-ERT)*ET
      ENDIF
! CH---エラーチェチE��---
    M_NaN='☁E5600: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ***************
! **** IUL=6 ****
! ***************
    6 CALL ROOT_CHECK(   6, 57)
      IF(DEU)6001,6002,6002
 6001 CALL ROOT_CHECK(6001, 57)
      IF(EU.LT.EPCU_) GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)  GOTO 8000
      ELSE
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EXC)  GOTO 1000
      ENDIF
      IF(EXC.LT.EU) GOTO 6000
 6002 CALL ROOT_CHECK(6002, 57)
      IF(EUOLD.LT.ECP) GOTO 6004
      IF(EPPC.EQ.EO)THEN
        IF(EU.LT.EPC) GOTO 6000
        IF(EPC.LE.EU) GOTO 6003
      ELSE
        EBC=EUOLD
        SBC=SNOLD
        SBB=SBC+(EU-EBC)*EPEC
        CALL CURVE3(   EU,  SAA,    &
                     EPC,           &
                     ETP,  STP, &
                    EPEC,   ET )    
        IF(SBB.LT.SAA) GOTO 7000
      ENDIF
      IF(EU.LE.EPC) GOTO 3000
 6003 CALL ROOT_CHECK(6003, 57)
      IF(CRACK.EQ.0.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        IF(ECCR.LT.EU)  GOTO 5000
        IF(EU.LE.ECCR)  GOTO 3010
      ENDIF
      IF(CRACK.EQ.4.0)THEN
        IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
      ENDIF
      IF(CRACK.EQ.5.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
        CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
        IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
          GOTO 5000
        ELSE
          GOTO 5100
        ENDIF
      ENDIF
 6004 CALL ROOT_CHECK(6004, 57)
        
! CH--------------------------------------------------
    IF(NROOT==0) GOTO 10020
    IF(NROOT6==0) GOTO 10020
    IF(EUOLD > EC .AND. SNOLD > 0.5*SC ) GOTO 10020 ! BY HJ 
    IF(EEN >= 0.0) GOTO 10020                       !  老E��E6666
      CALL CROSS1(EO,EUOLD,SNOLD,EEN,SEN,EO,SC,-SNOLD/EO+EUOLD, &
                EC,IUL,EPCU_,EPCUS_,ICIC(1))
    IF(EEN >= 0.0) GOTO 10020                       !  老E��E6666
    NCALL = 66666
    

      CALL ECTPJ_hj( EEN, SEN, ECP, SCP, ETP, STP,  &
                EPC, SPC,  EJ,  SJ,ECCR,SCCR,       &
               EPCU_,EPPC,                      &
                 EO,  SC,  FT,  ECU,CRACK,ICIC(1),&
                 EC,ETMAX,ECR,EPEC,IVCYCN ,NCYCNOFF,    &
                NST,NIT,M,LX,LY,LZ,intI,NCALL )

      CALL CROSS1(EPPC,ECP,SCP,EXC,SXC,EO,SC,EPC,   &
                EC,IUL,EPCU_,EPCUS_,ICIC(1))
10020   CONTINUE
! CH--------------------------------------------------

      IF(EPPC.EQ.EO)THEN
        IF(EU.LT.EPC) GOTO 6000
        IF(CRACK.EQ.0.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          IF(ECCR.LE.EU)  GOTO 5000
          IF(EU.LT.ECCR)  GOTO 3010
        ENDIF
        IF(CRACK.EQ.4.0)THEN
          IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
        ENDIF
        IF(CRACK.EQ.5.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
            GOTO 5000
          ELSE
            GOTO 5100
          ENDIF
        ENDIF
      ELSE
        IF(EU.LE.ETP)  GOTO 2000
        IF(EU.LE.EPC)  GOTO 3000
        IF(CRACK.EQ.0.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          IF(ECCR.LE.EU)  GOTO 5000
          IF(EU.LT.ECCR)  GOTO 3010
        ENDIF
        IF(CRACK.EQ.4.0)THEN
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
        ENDIF
        IF(CRACK.EQ.5.0)THEN
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
            GOTO 5000
          ELSE
            GOTO 5100
          ENDIF
        ENDIF
      ENDIF
 6000 CALL ROOT_CHECK(6000, 57)
      IUL=6
      INCF=1
      ET=EPPC
      SN=EPPC*(EU-EPC)
! CH---エラーチェチE��---
    M_NaN='☁E6000: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ***************
! **** IUL=7 ****
! ***************
    7 CALL ROOT_CHECK(   7, 57)
      IF(DEU)7001,7002,7002
 7001 CALL ROOT_CHECK(7001, 57)
      IF(EU.LT.EPCU_) GOTO 9000
      IF(EXC.LT.EC)THEN
        IF(EU.LE.EXC)  GOTO 8000
      ELSE
        IF(EU.LE.EC)   GOTO 8000
        IF(EU.LE.EXC)  GOTO 1000
      ENDIF
      IF(EU.LE.EBC)  GOTO 6000
      IF(EBC.LT.EU)  GOTO 7000
 7002 CALL ROOT_CHECK(7002, 57)
      SBB=SBC+(EU-EBC)*EPEC
      CALL CURVE3(   EU,  SAA,  &
                   EPC,     &
                   ETP,  STP,   &
                  EPEC,   ET )
      IF(SBB.LT.SAA) GOTO 7000
      IF(EU.LE.EPC)  GOTO 3000
      IF(CRACK.EQ.0.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        IF(ECCR.LT.EU)  GOTO 5000
        IF(EU.LE.ECCR)  GOTO 3010
      ENDIF
      IF(CRACK.EQ.4.0)THEN
        IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
      ENDIF
      IF(CRACK.EQ.5.0)THEN
        IF(EBU.LT.EU)   GOTO 5090
        CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
        CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
        IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
          GOTO 5000
        ELSE
          GOTO 5100
        ENDIF
      ENDIF
 7000 CALL ROOT_CHECK(7000, 57)
      IUL=7
      INCF=1
      ET=EPEC
      SN=SBC+(EU-EBC)*EPEC
! CH---エラーチェチE��---
    M_NaN='☁E7000: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ***************
! **** IUL=8 ****
! ***************
    8 CALL ROOT_CHECK(   8, 57)
      IF(DEU)8001,8002,8002
 8001 CALL ROOT_CHECK(8001, 57)                                  !除荷が起きてぁE��ぁE��E
      IF(EU.LT.EPCU_) GOTO 9000
      IF(EPCU_.LE.EU) GOTO 8000
 8002 CALL ROOT_CHECK(8002, 57)                                                !除荷が発生した、E
      EEN=EUOLD
      SEN=SNOLD
      ET=EO
! CH--------------------------------------------------
    IF(NROOT==0) GOTO 10030
    IF(NROOT8==0) GOTO 10030
    NCALL = 88888
      CALL ECTPJ_hj( EEN, SEN, ECP, SCP, ETP, STP,      &
                 EPC, SPC,  EJ,  SJ,ECCR,SCCR,          &
                EPCU_,EPPC,                         &
                  EO,  SC,  FT,  ECU,CRACK,ICIC(1), &
                  EC,ETMAX,ECR,EPEC,IVCYCN ,NCYCNOFF,   &
                NST,NIT,M,LX,LY,LZ,intI,NCALL )
      IF(EPPC.EQ.EO)THEN
        EXC=EEN
        SXC=SEN
      ELSE
        CALL CROSS1(EPPC,ECP,SCP,EXC,SXC,EO,SC,EPC, &
                  EC,IUL,EPCU_,EPCUS_,ICIC(1))
      ENDIF
10030   CONTINUE
! CH--------------------------------------------------
      IF(EPPC.EQ.EO)THEN
        IF(EU.LT.EPC) GOTO 6000                                                !圧縮域�Eまま�E�EUL=6�E�E
        IF(CRACK.EQ.0.0)THEN                                                   !引張域へ突�E
          IF(EBU.LT.EU)   GOTO 5090
          IF(ECCR.LE.EU)  GOTO 5000
          IF(EU.LT.ECCR)  GOTO 3010
        ENDIF
        IF(CRACK.EQ.4.0)THEN
          IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
        ENDIF
        IF(CRACK.EQ.5.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
            GOTO 5000
          ELSE
            GOTO 5100
          ENDIF
        ENDIF
      ELSE
        IF(EU.LE.ETP)  GOTO 2000
        IF(EU.LT.EPC)  GOTO 3000
        IF(CRACK.EQ.0.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          IF(ECCR.LE.EU)  GOTO 5000
          IF(EU.LT.ECCR)  GOTO 3010
        ENDIF
        IF(CRACK.EQ.4.0)THEN
          IF(EBU.LT.EU)   GOTO 4090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECR,SCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_41,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_41)THEN
            GOTO 4000
          ELSE
            GOTO 4100
          ENDIF
        ENDIF
        IF(CRACK.EQ.5.0)THEN
          IF(EBU.LT.EU)   GOTO 5090
          CALL SHIRAI(EU,DUMMY_SHIRAI,ECCR,SCCR,EBU,ICIC(4),C_IC4)
          CALL STRAIGHT41_51(EU,DUMMY_51,DUM,ETMAX,STMAX,EPC,VS)
          IF(DUMMY_SHIRAI.LT.DUMMY_51)THEN
            GOTO 5000
          ELSE
            GOTO 5100
          ENDIF
        ENDIF
      ENDIF
 8000 CALL ROOT_CHECK(8000, 57)
      IUL=8
      INCF=3
    CALL SOFTENING(FC,ECU,EC,SC,EPCUS_,EPCU_,EU,SN)
      ET=EO/100.
! CH---エラーチェチE��---
    M_NaN='☁E8000: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
      RETURN
! ***************
! **** IUL=9 ****
! ***************
    9 CALL ROOT_CHECK(   9, 57)
      IF(DEU)9001,9002,9002
 9001 CALL ROOT_CHECK(9001, 57)
      IF(ELIMIT.EQ.0.0) GOTO 9000
      IF(EU.LE.ELIMIT)  GOTO 9000
      IF(ELIMIT.LT.EU)  GOTO 9002
 9000 CALL ROOT_CHECK(9000, 57)
      IUL=9
      INCF=3
      ELIMIT=EU
      SN=EPCUS_
      ET=EO/100.
      RETURN
! CH---エラーチェチE��---
    M_NaN='☁E9000: [ SN ET ]  NST NIT M LX LY LZ intI '
    N_NaN(1) = NST ; N_NaN(2) = NIT ; N_NaN(3) = M
    N_NaN(4) = LX ; N_NaN(5) = LY ; N_NaN(6) = LZ ; N_NaN(7) = intI
    S_NaN(1) = SN ; S_NaN(2) = ET
    ! CALL CCHJ(N_NaN, 7, S_NaN, 2 ,M_NaN,ID_NaN)
! CH--- --- --- --- --- --- ---
 9002 CALL ROOT_CHECK(9002, 57)
      SN=0.0001
! C      SN=EPCUS_ ! BY HJ
      ET=EO/100.
      RETURN
! *******************************************************************
! * ERROR MESSAGES                                                  *
! *******************************************************************
  901 WRITE(6,902)
  902 FORMAT(' IVIRGIN ERROR')
      STOP
  903 WRITE(6,904)
  904 FORMAT(' ***ERROR IN SECON2 (IUL=',I2,')***')

      STOP
      END


!***********************************************************************
      SUBROUTINE HJ_C_1(C,M10)
      ! IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      ! INCLUDE 'par.f'
      ! INCLUDE 'common.f'
    CHARACTER M10*10
    real (kreal) C

    ! DO I=1,2
        ! IF(I==1) THEN ; N=1820 ; ELSE ; N=6 ; ENDIF
        ! WRITE(N,'(1A6,1I5,1A10,1ES20.5)') '●●◁E,NST,ADJUSTR(M10),C
    ! ENDDO

    END

!*******************************************************************
      SUBROUTINE EPSSIG(EO,ECU,SC,FC,EC,EPCUS_,EPCU_,EPS,SIG)
! *******************************************************************
! * 匁E��線{SAENZ(EO,ECU,SC) or SOFTENING(FC,EC,SC,EPCUS,EPCU)}上�E�E�EPSに対応するSIGを求める、E
    
    real (kreal) :: EO,ECU,SC,FC,EC,EPCUS_,EPCU_,EPS,SIG,ET
      

    IF( EPS <= EC .AND. EPS >= EPCU_ ) THEN
        CALL SOFTENING(FC,ECU,EC,SC,EPCUS_,EPCU_,EPS,SIG)
    ELSEIF( EPS > EC .AND. EPS < 0.0 ) THEN
        CALL SAENZ(EO,ET,EC,SC,EPS,SIG)
    ELSEIF( EPS < EPCU_ ) THEN
        SIG = EPCUS_
    ELSE
        SIG=0.0
    ENDIF

    RETURN
    END


    !*******************************************************************
      SUBROUTINE SOFTENING(FC,ECU,EC,SC,EPCUS_,EPCU_,EU,SN)

! *******************************************************************
      real (kreal) :: FC,ECU,EC,SC,EPCUS_,EPCU_,EU,SN,ES

      ES=(EPCUS_-SC)/(EPCU_-EC)
      SN=SC+ES*(EU-EC)

    END

    
    ! *******************************************************************
      SUBROUTINE XPOINT(EEUU,SSNN,EO,ECU,SC,FC,EC,EPCUS_,EPCU_,COPX,    &
                       EXC,SXC)
! *******************************************************************
! * 点(EU,SN)を通り�E�剛性Eを持つ直線と匁E��線{SAENZ(EO,ECU,SC) or SOFTENING(FC,EC,SC,EPCUS,EPCU)}との
! * 交点(EXC,SXC)を求める、E
! *引数�E�　 EO  �E�弾性初期接線係数
! *     ET  �E�当応力点で接線剛性    
! *     ECU �E�最大圧縮応力時�E等価一軸ひずみ 
! *     SC  �E�最大圧縮応力   
! *     EPS �E�条件としてのひずみ    
! *     SIG �E�求める応力
! C
     real(kreal) :: EEUU,SSNN,EO,ECU,SC,FC,EC,EPCUS_,EPCU_,COPX,EXC,SXC
     real(kreal) :: R,EU,SN,E,EPS,EK,EKI,EE,ET,H,RL,SIG,SS,ES
     integer :: N
    R=SC/100000000.
! C
5   N=0
    EU=EEUU
    SN=SSNN
    E=COPX
    EPS=EU
    EK = (SC-SN)/(EC-EU)
    EKI= (EPCUS_-SN)/(EPCU_-EU)

20  IF( EU > EC ) THEN
        IF( EK*0.999 > E ) THEN
            CALL SOFTENING(FC,ECU,EC,SC,EPCUS_,EPCU_,EPS,SIG)
            CALL vertical(E,EU,SN,EPS,SIG,EE,SS,H)
        ELSEIF( 1.001*EK < E ) THEN
            CALL SAENZ(EO,ET,EC,SC,EPS,SIG)
            CALL vertical(E,EU,SN,EPS,SIG,EE,SS,H)
        ELSE
            SIG = SC
            H=0.0
        ENDIF
    ELSEIF( EU <= EC .AND. EU >= EPCU_ ) THEN
        IF( EKI*0.999 > E ) THEN
            SIG = EPCUS_
            H=0.0
        ELSEIF( 1.001*EKI < E ) THEN
            CALL SOFTENING(FC,ECU,EC,SC,EPCUS_,EPCU_,EPS,SIG)
            CALL vertical(E,EU,SN,EPS,SIG,EE,SS,H)
        ELSE
            SIG = EPCUS_
            H=0.0
        ENDIF
    ELSEIF( EU < EPCU_ ) THEN
        SIG = EPCUS_
        H=0.0
    ENDIF

    IF( DABS( H ) < DABS( R ) )  THEN
        GOTO 10
    ELSE
        N=N+1
        IF( N > 200 ) GOTO 30
        IF( N==1 .OR. RL > DABS( H ) ) RL=DABS( H )
        EPS=EE
        GOTO 20
    END IF
30  CONTINUE
    R = RL * 1.001
    GOTO 5

10  CONTINUE

    EXC=EPS
    SXC=SIG

    RETURN
    END
    
    !*******************************************************************
      SUBROUTINE vertical(E,x1,y1,x2,y2,X,Y,H)
! *******************************************************************
! * 点(x2,y2)から�E�方向余弦(1,E)で点(x1,y1)を通る直線への投影座樁EX,Y)と垂直距離 H を求める、E
! *
      real (kreal) :: E,x1,y1,x2,y2,X,Y,H,a,b,a1,b1,c1,a2,b2,c2
       real (kreal) :: z1,z2

     a = 1.
     b = E
! C 平衁E
     a1 = b
     b1 = - a
     c1 = a * y1 - b * x1
!  C    垂直
     a2 = a
     b2 = b
     c2 = - a * x2 - b * y2

     z1 = c1 * b2 / b1 - c2
     z2 = a2 - a1 * b2 / b1

     X = z1 / z2
     Y = ( - a1 * X - c1 ) / b1
     H = sqrt( ( x2 - X )**2 + ( y2 - Y )**2 )

    END
    
    ! *******************************************************************
      SUBROUTINE FS(EO,ET,ECU,SC,EPS,SIG)
! *******************************************************************
! ****  ADDED BY JOE ZAKKY ON 1996.8.14 FOR BECOMING MORE EASILY ****
! **** TO ADJUST, CHECK OR CHANGE INTO OTHER EQUATIONS.          ****
! *******************************************************************

! *$INCLUDE"D:\ANALYSIS\ORIGINAL\TS\H\FILE.H"
! *$INCLUDE'/USER2/ZAKKY/Z9/H/FILE.H'
    real (kreal) :: EO,ET,ECU,SC,EPS,SIG
    real (kreal) :: ZZ, Z 
      ! CALL ROOT_CHECK(0000, 29)

      ZZ =EPS/ECU
      Z  =EO*ECU/SC
      SIG=SC*(1.-(1.-ZZ)**Z)
      ET =EO*(1.-ZZ)**(Z-1.)

      ! CALL ROOT_CHECK(9999, 29)
      RETURN
      END
    
    
    ! *******************************************************************
      SUBROUTINE ECTPJ_hj(  EEN,  SEN,  ECP,  SCP,                  &
                        ETP,  STP,  EPC,  SPC,                      &
                         EJ,   SJ, ECCR, SCCR,                      &
                      EPCU_, EPPC,   EO,   FC,                      &
                         FT,   ECU,CRACK, ICIC1,                    &
                         EC,ETMAX,  ECR, EPEC, IVCYCN ,NCYCNOFF,    &
                        NST,NIT,M,LX,LY,LZ,intI,NCALL)!NCYCNOFF
! *******************************************************************
! ****  SUBROUTINE FOR SOLVING POINT LOCATION OF                  ***
! ****  E, C, T, P, J AND T'                                     ****
! *******************************************************************

    real (kreal) :: EEN,SEN,ECP,SCP,ETP,STP,EPC,SPC,EJ,SJ,ECCR,SCCR
    real (kreal) :: EPCU_,EPPC,EO,FC,FT,ECU,CRACK,EC,ETMAX,ECR,EPEC
    integer :: ICIC1,IVCYCN,NCYCNOFF,NST,NIT,M,LX,LY,LZ,intI,NCALL
    
    real (kreal) :: XHJ(999),C_Infinity(1),S_NaN(99),N_NaN(99)
    real (kreal) :: aaaaaa,alfa,bbb,epecc,epep
    integer :: id_nan,n
! c      INCLUDE 'com_hj.f'

    CHARACTER*100 M_NaN
! CH---老E��EECTPJ------------------------------------------------
    IF(EEN >= 0.0) THEN
        WRITE(*,'(A72,8I6,ES11.3)')'●●●老E��_ECTPJ_hj [ EEN ] :  &
      ≧0.0,EEN NCALL NST NIT M LX LY LZ intI'                       &
                    ,NCALL,NST,NIT,M,LX,LY,LZ,intI,EEN
        GOTO 9999
    ENDIF
    M_NaN='☁EECTPJ_hj [ EEN ] : NCALL NST NIT M LX LY LZ intI'
    N_NaN(1) = NCALL ;N_NaN(2) =  NST ;N_NaN(3) =  NIT ;N_NaN(4) =  M
    N_NaN(5) = LX  ;N_NaN(6) = LY  ;N_NaN(7) = LZ  ;N_NaN(8) = intI
    S_NaN(1) = EEN
    ID_NaN=0
    ! CALL CCHJ(N_NaN, 8, S_NaN, 1 ,M_NaN,ID_NaN)
    IF(ID_NaN==1) GOTO 9999

! CH---老E��EECTPJ------------------------------------------------

    IF(IVCYCN == 1) GOTO 10

      IF(DABS(EEN).LT.DABS(4.*EC)) THEN
        EPC=(0.145*EEN**2./EC)+0.127*EEN ; EPC=EPC !�E�点ひずみ
        SPC=0.0                                    !�E�点応力
    ELSE
        EPC=(EEN/EC-2.828)*EC/5.                      !�E�点ひずみ
        SPC=0.0                                    !�E�点応力
    ENDIF
10  CONTINUE

    EPEP=SEN/(EEN-EPC)
! CH---修正5_hj------------------------------------------------
    IF(EEN >= EPC.OR.EEN >=0.0.OR.EPC >=0.0) THEN
        WRITE(1850,'(A80,8I6,2ES11.3)')                                 &
     '●●●修正51_hj:EEN≧EPC,EEN+,EPC+,NCALL NST NIT M LX LY LZ intI'   &
                    ,NCALL,NST,NIT,M,LX,LY,LZ,intI,EEN,EPC
        EPEP=EO/1.5
    ENDIF
    M_NaN='☁EECTPJ_hj_52:[SEN EEN EPC ] NCALL NST NIT M LX LY LZ intI'
    N_NaN(1) = NCALL ;N_NaN(2) =  NST ;N_NaN(3) =  NIT ;N_NaN(4) =  M
    N_NaN(5) = LX  ;N_NaN(6) = LY  ;N_NaN(7) = LZ  ;N_NaN(8) = intI
    S_NaN(1) = SEN ; S_NaN(2) = EEN ; S_NaN(3) = EPC
    ID_NaN=0
    ! CALL CCHJ(N_NaN, 8, S_NaN, 3 ,M_NaN,ID_NaN)
    IF(ID_NaN==1) EPEP=EO/1.5
! CH---修正5_hj------------------------------------------------
      EPEC=1.5*EPEP                     !EC間�E剛性

      IF(EPEC.GT.EO) EPEC=EO                     !EC間�E剛性は初期剛性以丁E
      EPECC=2*SEN/(EEN-EPC)
      IF(EPEC.GT.EPECC) EPEC=EPECC               ! 点�E�で剛性ぁE以下にならなぁE��ぁE��する処置、E
! C IF(EPEC < EO/2.) EPEC=EO/2.
      IF(EEN.GE.EC) THEN
        SCP=0.833333*SEN                           !�E�点応力1
    ELSE
        SCP=MIN(0.66667*SEN,SEN-0.166667*FC)       !�E�点応力2E
      ENDIF

      BBB=SEN-EPEC*EEN
    ECP=EEN-(SEN-SCP)/EPEC
      ETP=EEN-SEN/2.0/EPEC                       !�E�点ひずみ
      STP=SEN-SEN/2.0                            !�E�点応力
    IF(NCYCNOFF==1)THEN
        SCP=SEN !変更-C点をX点と一致させ�E�C点を廁E��E��EJ
        ECP=EEN !変更-C点をX点と一致させ�E�C点を廁E��E��EJ
    ENDIF
    EPPC=SCP/(ECP-EPC)               !IUL=6の剛性

! CH---修正4_hj------------------------------------------------
    IF( ECP >= EPC ) THEN
        WRITE(1850,'(A60,8I6,2ES11.3)')                     &
     '●●●修正41_hj:ECP≧EPC NCALL NST NIT M LX LY LZ intI',    &
                    NCALL,NST,NIT,M,LX,LY,LZ,intI,ECP,EPC
        EPPC=EO
    ELSE
        EPPC=SCP/(ECP-EPC)               !IUL=6の剛性
    ENDIF   
    
    M_NaN='☁EECTPJ_hj_42:                                   &
     [SCP ECP EPC SEN EEN EPEC] NCALL NST NIT M LX LY LZ intI'
    N_NaN(1) = NCALL ;N_NaN(2) =  NST ;N_NaN(3) =  NIT ;N_NaN(4) =  M
    N_NaN(5) = LX  ;N_NaN(6) = LY  ;N_NaN(7) = LZ  ;N_NaN(8) = intI
    S_NaN(1) = SCP ; S_NaN(2) = ECP ; S_NaN(3) = EPC ; S_NaN(4) = SEN
    S_NaN(5) = EEN ; S_NaN(6) = EPEC
    ID_NaN=0
    ! CALL CCHJ(N_NaN, 8, S_NaN, 6 ,M_NaN,ID_NaN)
    IF(ID_NaN==1) EPPC=EO
! CH---修正4_hj------------------------------------------------
      !CALL ROOT_CHECKR1('EPPC', EPPC)
      IF(EO.LT.EPPC)THEN
        EPPC=EO
        EPC=EEN-SEN/EO
      ENDIF
    IF (EPC.GT.0.0) THEN
        EPC=0.0
        EPPC=SEN/EEN
    ENDIF   

100 CONTINUE
      IF(CRACK.EQ.4.0)THEN
        IF(EEN.EQ.0.0)THEN
200 CONTINUE
! *         SJ=-FT*1.5                             !IUL=2を未経騁E
          SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)*15      !FROM KOBAYASHI(MODIFIED BY SAKURAI)
          SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)*1     !FROM KOBAYASHI(MODIFIED BY TIDE)
! *         EJ=SJ/EPPC+EPC                         !FROM KOBAYASHI
          IF(ICIC1.EQ.0)THEN
  300 CALL ROOT_CHECK( 300, 23)
            CALL SAENZ_REVERS(EJ,SJ,EO,ECU,FC)    !SAENZ式送E��数
          ELSE IF(ICIC1 == 1 )THEN
  400       CALL ROOT_CHECK( 400, 23)
            CALL SHAH_REVERS (EJ,SJ,EO,ECU,FC)    !FAFITIS-SHAH式送E��数
          ELSE IF(ICIC1 == 2 )THEN
            CALL Kimura_REVERS(EJ,SJ,EO,ECU,FC)
          ENDIF
        ELSE
500 CONTINUE
! *         SJ=SEN*0.1                             !IUL=2を経騁E
! C          SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)*15      !FROM KOBAYASHI(MODIFIED BY SAKURAI)
          SJ=-FT*(1.0+0.02*(ETMAX-ECR)/ECR)*1       !FROM KOBAYASHI(MODIFIED BY TIDE)
          EJ=SJ/EPPC+EPC
        ENDIF
      ELSE
  600   CALL ROOT_CHECK( 600, 23)

! C        SJ=SEN*0.1*15                               !IUL=2を経騁E
        SJ=SEN*0.1*1                             !IUL=2を経騁E
        EJ=SJ/EPPC+EPC
      ENDIF
    IF(EJ.GE.0.0)THEN
        AAAAAA=1.
    ENDIF
  700 CALL ROOT_CHECK( 700, 23)
! c      EPCU_=-.01
      IF(CRACK.NE.0.0)GOTO 9999                  !圧縮効果を老E�Eしたひび割れ発生点の算�E
      ALFA=0.0
      ALFA=FT/(FT/EO-EPCU_)
! C      ECCR=(EPPC*EPC-ALFA*EPCU_)/(EPPC-ALFA)
    ECCR=EPC+(1-EPC/EPCU_)*(FT/EO)
      SCCR=EPPC*ALFA*(EPCU_-EPC)/(ALFA-EPPC)
! CH---修正6------------------------------------------------
    IF( ALFA==EPPC ) THEN
        WRITE(1850,'(A60,8I6,2ES11.3)')                         &
     '●●●修正6_hj:ALFA=EPPC NCALL NST NIT M LX LY LZ intI',   &
                    NCALL,NST,NIT,M,LX,LY,LZ,intI,ALFA,EPPC
        SCCR=FT
    ELSE
    ENDIF
    M_NaN='☁E ECTPJ_hj_62:[ EPPC ALFA EPCU_ EPC ]       &
     NCALL NST NIT M LX LY LZ intI'
    N_NaN(1) = NCALL ;N_NaN(2) =  NST ;N_NaN(3) =  NIT ;N_NaN(4) =  M
    N_NaN(5) = LX  ;N_NaN(6) = LY  ;N_NaN(7) = LZ  ;N_NaN(8) = intI
    S_NaN(1)=EPPC; S_NaN(2)=ALFA ; S_NaN(3) = EPCU_;S_NaN(4) = EPC
    N=0
    ! CALL CCHJ(N_NaN, 8, S_NaN, 4 ,M_NaN,N)
    IF(N==1) SCCR=FT
! CH---修正6------------------------------------------------
    IF(SJ.LT.-100.0) THEN 
    CONTINUE
    ENDIF
 9999 RETURN
      END
    
        ! !************************************************************************
    ! SUBROUTINE CCHJ(N_NaN,NN,S_NaN,NS,M_NaN,ID_NaN)
    ! USE DFLIB 

! ! c      INCLUDE 'com_hj.f'
      ! DIMENSION C_Infinity(1)

    ! DIMENSION S_NaN(99),N_NaN(99)
    ! CHARACTER*100 M_NaN

    ! result = LEN_TRIM(M_NaN)
    ! Y=0. ; Y1=0. ; Y2=0.
    ! DO I=1,NS
        ! IF(ABS(S_NaN(I))==C_Infinity(1)) Y1=1.           !modify by sun 
        ! Y2=ISNAN(S_NaN(I))
        ! IF(Y1 == 1. .OR. Y2 < 0. ) Y=1.
    ! ENDDO

    ! IF( Y == 1. ) THEN
        ! ID_NaN=1
        ! WRITE(1850,'(A3,A<result>,TR2,<NN>I6,TR2,<NS>ES11.3)')    &
                ! '☁E',M_NaN,(N_NaN(I),I=1,NN),(S_NaN(I),I=1,NS)
            ! WRITE(*,'(A3,A<result>,TR2,<NN>I6,TR2,<NS>ES11.3)') 
                ! '☁E',M_NaN,(N_NaN(I),I=1,NN),(S_NaN(I),I=1,NS)    
    ! ENDIF

    ! END
    
    !C***********************************************************************
      SUBROUTINE RG12RG23RG31(NTCRACK,LCRACK,EUC1,EUC2,EUC3,G0,ICIC,RG12,RG23,RG31,GBETA)
! C***********************************************************************
    real (kreal) :: EUC1,EUC2,EUC3,G0,RG12,RG23,RG31,GBETA
    integer :: LCRACK(3),ICIC(8)
    integer :: NTCRACK

    real (kreal):: GKK,EU_1,EU_2,EU_3,GBAI,BETA
    real(kreal) :: GGGG_1,GGGG_2,GGGG_3,GGGG_bai,RRGG
    integer :: NG3C
    
    
    
      ! INCLUDE 'par.f'
      ! INCLUDE 'common.f'


    GKK=G0
    EU_1 = EUC1
    EU_2 = EUC2
    EU_3 = EUC3

    IF (ICIC(5)==0) GOTO 10
    IF (ICIC(5)==1) GOTO 100
    IF (ICIC(5)==2) GOTO 200
    IF (ICIC(5)==3) GOTO 300
! CH------------ AL --------------
   10 CONTINUE
      RETURN
! CH------------前巁E-------------
  100 CONTINUE
    NG3C=0
    IF(NG3C==1) THEN
        GBAI=0.99
    ELSE IF(NG3C==2) THEN
        GBAI=0.2
    ELSE IF(NG3C==3) THEN
        GBAI=0.02
    ENDIF

    IF(NTCRACK.EQ.1) GOTO 110

    IF     (LCRACK(1)==1.AND.LCRACK(2)==1.AND.LCRACK(3)/=1) THEN
        RG12 = 1./(1./GKK+1./RG12+1./RG12)
        RG23 = 1./(1./GKK+1./RG23)
        RG31 = 1./(1./GKK+1./RG31)
        IF(NG3C/=0) THEN
            RG12 = GKK*GBAI/2.
            RG23 = GKK*GBAI
            RG31 = GKK*GBAI
        ENDIF
    ELSE IF(LCRACK(1)/=1.AND.LCRACK(2)==1.AND.LCRACK(3)==1) THEN
        RG12 = 1./(1./GKK+1./RG12)
        RG23 = 1./(1./GKK+1./RG23+1./RG23)
        RG31 = 1./(1./GKK+1./RG31)
        IF(NG3C/=0) THEN
            RG12 = GKK*GBAI
            RG23 = GKK*GBAI/2.
            RG31 = GKK*GBAI
        ENDIF
    ELSE IF(LCRACK(1)==1.AND.LCRACK(2)/=1.AND.LCRACK(3)==1) THEN
        RG12 = 1./(1./GKK+1./RG12)
        RG23 = 1./(1./GKK+1./RG23)
        RG31 = 1./(1./GKK+1./RG31+1./RG31)
        IF(NG3C/=0) THEN
            RG12 = GKK*GBAI
            RG23 = GKK*GBAI
            RG31 = GKK*GBAI/2.
        ENDIF
    ELSE IF(LCRACK(1)==1.AND.LCRACK(2)==1.AND.LCRACK(3)==1) THEN
        RG12 = 1./(1./GKK+1./RG12+1./RG12)
        RG23 = 1./(1./GKK+1./RG23+1./RG23)
        RG31 = 1./(1./GKK+1./RG31+1./RG31)
        IF(NG3C/=0) THEN
            RG12 = GKK*GBAI/2.
            RG23 = GKK*GBAI/2.
            RG31 = GKK*GBAI/2.
        ENDIF
    ELSE
        RG12 = 1./(1./GKK+1./RG12+1./RG12)
        RG23 = 1./(1./GKK+1./RG23+1./RG23)
        RG31 = 1./(1./GKK+1./RG31+1./RG31)
        IF(NG3C/=0) THEN
            RG12 = GKK*GBAI/2.
            RG23 = GKK*GBAI/2.
            RG31 = GKK*GBAI/2.
        ENDIF
    ENDIF
      RETURN

  110 IF(LCRACK(1).EQ.1) THEN
        RG12 = 1./(1./GKK+1./RG12)
        RG31 = 1./(1./GKK+1./RG31)
        IF(NG3C/=0) THEN
            RG12 = GKK*GBAI
            RG31 = GKK*GBAI
        ENDIF
    RETURN
      ENDIF

      IF(LCRACK(2).EQ.1) THEN
        RG12 = 1./(1./GKK+1./RG12)
        RG23 = 1./(1./GKK+1./RG23)
        IF(NG3C/=0) THEN
            RG12 = GKK*GBAI
            RG23 = GKK*GBAI
        ENDIF
    RETURN
      ENDIF

      IF(LCRACK(3).EQ.1) THEN
        RG23 = 1./(1./GKK+1./RG23)
        RG31 = 1./(1./GKK+1./RG31)
        IF(NG3C/=0) THEN
            RG23 = GKK*GBAI
            RG31 = GKK*GBAI
        ENDIF
    RETURN
      ENDIF
! CH------------ 青柳 -----------
  200 CONTINUE
    IF( EU_1 < 0.0001 ) EU_1 = 0.0001
    IF( EU_2 < 0.0001 ) EU_2 = 0.0001
    IF( EU_3 < 0.0001 ) EU_3 = 0.0001
    GGGG_1 = 360./EU_1
    GGGG_2 = 360./EU_2
    GGGG_3 = 360./EU_3
    GGGG_bai = 0.0

    IF(NTCRACK.EQ.1) GOTO 210

    IF     (LCRACK(1)==1.AND.LCRACK(2)==1.AND.LCRACK(3)/=1) THEN
        RG12 = 1./ ( 1./GKK + 1./GGGG_1 + 1./GGGG_2 )
        RG23 = 1./ ( 1./GKK +             1./GGGG_2 )
        RG31 = 1./ ( 1./GKK + 1./GGGG_1             )
    ELSE IF(LCRACK(1)/=1.AND.LCRACK(2)==1.AND.LCRACK(3)==1) THEN
        RG12 = 1./ ( 1./GKK + 1./GGGG_2             )
        RG23 = 1./ ( 1./GKK + 1./GGGG_2 + 1./GGGG_3 )
        RG31 = 1./ ( 1./GKK +             1./GGGG_3 )
    ELSE IF(LCRACK(1)==1.AND.LCRACK(2)/=1.AND.LCRACK(3)==1) THEN
        RG12 = 1./ ( 1./GKK +             1./GGGG_1 )
        RG23 = 1./ ( 1./GKK + 1./GGGG_3             )
        RG31 = 1./ ( 1./GKK + 1./GGGG_3 + 1./GGGG_1 )
    ELSE IF(LCRACK(1)==1.AND.LCRACK(2)==1.AND.LCRACK(3)==1) THEN
        RG12 = 1./ ( 1./GKK + 1./GGGG_1 + 1./GGGG_2 )
        RG23 = 1./ ( 1./GKK + 1./GGGG_2 + 1./GGGG_3 )
        RG31 = 1./ ( 1./GKK + 1./GGGG_3 + 1./GGGG_1 )
    ELSE
        RG12 = 1./ ( 1./GKK + 1./GGGG_1 + 1./GGGG_2 )
        RG23 = 1./ ( 1./GKK + 1./GGGG_2 + 1./GGGG_3 )
        RG31 = 1./ ( 1./GKK + 1./GGGG_3 + 1./GGGG_1 )
    ENDIF
    RG12 = RG12+(GKK-RG12) * GGGG_bai
    RG23 = RG23+(GKK-RG23) * GGGG_bai
    RG31 = RG31+(GKK-RG31) * GGGG_bai
      RETURN

  210 IF(LCRACK(1).EQ.1) THEN
        RRGG = 1./(1./GKK+1./GGGG_1)
        RRGG = RRGG+(GKK-RRGG) * GGGG_bai
        RG12 = RRGG
        RG31 = RRGG
    RETURN
      ENDIF

      IF(LCRACK(2).EQ.1) THEN
        RRGG = 1./(1./GKK+1./GGGG_2)
        RRGG = RRGG+(GKK-RRGG) * GGGG_bai
        RG12 = RRGG
        RG23 = RRGG
    RETURN
      ENDIF

      IF(LCRACK(3).EQ.1) THEN
        RRGG = 1./(1./GKK+1./GGGG_3)
        RRGG = RRGG+(GKK-RRGG) * GGGG_bai
        RG23 = RRGG
        RG31 = RRGG
    RETURN
      ENDIF
! CH -----------一定保有せん断係数----------------
  300 CONTINUE
    BETA=GBETA
    IF( BETA < 0.1 ) BETA=0.1
    IF( BETA > 1.0 ) BETA=1.0
    RRGG = GKK * BETA

    IF(NTCRACK.EQ.1) GOTO 310

    IF     (LCRACK(1)==1.AND.LCRACK(2)==1.AND.LCRACK(3)/=1) THEN
        RG12 = RRGG / 2.
        RG23 = RRGG
        RG31 = RRGG
    ELSE IF(LCRACK(1)/=1.AND.LCRACK(2)==1.AND.LCRACK(3)==1) THEN
        RG12 = RRGG
        RG23 = RRGG / 2.
        RG31 = RRGG
    ELSE IF(LCRACK(1)==1.AND.LCRACK(2)/=1.AND.LCRACK(3)==1) THEN
        RG12 = RRGG
        RG23 = RRGG
        RG31 = RRGG / 2.
    ELSE IF(LCRACK(1)==1.AND.LCRACK(2)==1.AND.LCRACK(3)==1) THEN
        RG12 = RRGG / 2.
        RG23 = RRGG / 2.
        RG31 = RRGG / 2.
    ELSE
        RG12 = RRGG / 2.
        RG23 = RRGG / 2.
        RG31 = RRGG / 2.
    ENDIF
      RETURN

  310 IF(LCRACK(1).EQ.1) THEN
        RG12 = RRGG
        RG31 = RRGG
    RETURN
      ENDIF

      IF(LCRACK(2).EQ.1) THEN
        RG12 = RRGG
        RG23 = RRGG
    RETURN
      ENDIF

      IF(LCRACK(3).EQ.1) THEN
        RG23 = RRGG
        RG31 = RRGG
    RETURN
      ENDIF

    END
    
    ! C***********************************************************************
      SUBROUTINE JACOBI (A,V,EPS,L,M,NR,NR2,ERR)
! C***********************************************************************
! C 参�E�E�FORTRAN77による　数値計算ソフトウェア　P300
! C 渡部　力、名取　亮、小国　劁E
! C �E�丸喁E��式会社�E�E
        integer :: L,M,NR,NR2
        real (kreal) ::  A(L,M),V(L,M),EPS
        
        
        real (kreal) :: A1,A2,A3,CO,SI,T1,T2,T3,T4,TA,W,WMAX
        integer :: i,j,k,ll,m1,m2,n
      INTEGER :: ERR


! C4      INTEGER I,J,K,L,LL,M,M1,M2,N,NR,ERR
! C4      REAL A(L,M),V(L,M),EPS,W,WMAX,A1,A2,A3,
! C4     1     T1,T2,T3,T4,TA,SI,CO
! C
! C     WRITE(6,*)
! CCCCC
! C     WRITE(6,*) '***** IN JACOBI *****'
! C     WRITE(6,*) 'EIGENVALUE'
! C     DO 21 I=1,M
! C  21 WRITE(6,639) A(I,I)
! C     WRITE(6,*) 'EIGENVECTORS'
! C     DO 31 I=1,M
! C  31 WRITE(6,649) (V(I,J),J=1,M)
! C
! C     WRITE(6,*)'L=',L,'M=',M,'NR=',NR,'NR2=',NR2
! C     WRITE(6,*)'EPS=',EPS,'ERR=',ERR
! C
  639 FORMAT(1PE15.7)
  649 FORMAT(3(1PE15.7))

      IF(M.LT.2.OR.NR.LE.0.OR.EPS.LE.0.0D0) THEN
        ERR=999
        NR2=0
      ELSE
        ERR=0
        N=0                         !n:計算回数
        M1=M-1
        DO 20 I=1,M1
          M2=I+1
          DO 10 J=M2,M
            IF (A(I,J).NE.A(J,I)) THEN
              ERR=999
              NR2=0
              RETURN
            ENDIF
   10     CONTINUE
   20   CONTINUE
        DO 40 I=1,M
          DO 30 J=1,M
            V(I,J)=0.0D0
   30     CONTINUE
          V(I,I)=1.0D0
   40   CONTINUE
   50   WMAX=0.0D0
        DO 70 I=1,M1
          M2=I+1
          DO 60 J=M2,M
            W=ABS(A(I,J))
            IF (W.GT.WMAX) THEN
              WMAX=W                    !a(i,j)�E�角度の頁E���E最大値を探ぁE
              K=I                       !i,jを記録する
              LL=J
            ENDIF
   60     CONTINUE
   70   CONTINUE
        IF(WMAX.LE.EPS) THEN
          NR2=N
          RETURN
        ENDIF
        IF(N.GE.NR) THEN
! CCCCCCC
! C         WRITE(6,*) '++++++++++'
! C         WRITE(6,*) 'N=',N,'NR=',NR,'NR2=',NR2
! C         WRITE(6,*) '++++++++++'
! CCCCCCC
          ERR=1
          NR2=N
          RETURN
        ENDIF
        N=N+1
        A1=A(K,K)
        A2=A(LL,LL)
        A3=A(K,LL)
        T1=ABS(A1-A2)
        T2=T1*T1
        T3=4.0D0*A3*A3
        TA=2.0D0*A3/(T1+SQRT(T2+T3))
        IF(A1.LT.A2) TA=-TA
        T4=TA*TA+1.0D0
        CO=SQRT(1.0D0/T4)
        SI=TA*CO
        DO 80 I=1,M
          W=V(I,K)
          V(I,K)=W*CO+V(I,LL)*SI
          V(I,LL)=-W*SI+V(I,LL)*CO
   80   CONTINUE
        DO 90 I=1,M
          IF(I.NE.K.AND.I.NE.LL) THEN
            W=A(I,K)
            A(I,K)=W*CO+A(I,LL)*SI
            A(I,LL)=-W*SI+A(I,LL)*CO
            A(K,I)=A(I,K)
            A(LL,I)=A(I,LL)
          ENDIF
   90   CONTINUE
        A(K,K)=A1*CO*CO+A3*(CO+CO)*SI+A2*SI*SI
        A(LL,LL)=A1+A2-A(K,K)
        A(K,LL)=0.0D0
        A(LL,K)=0.0D0
        GOTO 50
      ENDIF
      END
    
    subroutine oldcalConcreteDmat(strain,matvar,fail,D)
    !********************************************************************************************
    ! Calculates the Elastic Element D matrix before cracking
    !           Shanthanu 2017/13/3
    !********************************************************************************************
        real(kind=kreal), intent(in)  :: strain(6)
        real(kind=kreal), intent(in)  :: matvar(:)
        integer, intent(in) :: fail
        real(kind=kreal),intent(out) :: D(6,6)
        real(kind=kreal) :: eigval(3),EE,PP     
        real(kind=kreal) :: princ(3,3)    
        integer :: i,j
                        
        D=0.0d0
        PP=matvar(2)
        call eigen3(strain,eigval,princ)            
        if(minval(eigval).ge.0) then 
            EE=matvar(1)
        else
            call oldcalConcreteCompMod(-minval(eigval),matvar(5),matvar(6),EE)
        end if
        
        if (fail==1) PP=0.0d0
        D(1,1)=EE*(1.d0-PP)/(1.d0-2.d0*PP)/(1.d0+PP)
        D(1,2)=EE*PP/(1.d0-2.d0*PP)/(1.d0+PP)
        D(1,3)=D(1,2)
        D(2,1)=D(1,2)
        D(2,2)=D(1,1)
        D(2,3)=D(1,2)
        D(3,1)=D(1,3)
        D(3,2)=D(2,3)
        D(3,3)=D(1,1)
        D(4,4)=EE/(1.d0+PP)*0.5d0
        D(5,5)=EE/(1.d0+PP)*0.5d0
        D(6,6)=EE/(1.d0+PP)*0.5d0
    end subroutine
    
    subroutine oldcalUpdateConcrete(gauss)
    !********************************************************************************************
    ! Calculates the Stress in Cracked concrete
    !           Shanthanu 2017/13/3
    !********************************************************************************************
    type( tGaussStatus ), intent(inout) :: gauss                    !> material properties
    real(kind=kreal) :: Dmat(6,6)
    real(kind=kreal) :: strainvec(6),stressvec(6),e0,fcr,Ec,fck
    real(kind=kreal) :: rotmat(3,3),strainmat(3,3),stressmat(3,3)
    integer :: i,j
    
    strainvec=gauss%strain(1:6)
    Dmat=0.d0
    stressvec=0.d0
    if (gauss%istatus(1)==0) then
        call oldcalConcreteDmat(strainvec,gauss%pMaterial%variables,&
                        gauss%istatus(1),Dmat)
        stressvec=matmul(Dmat,strainvec)    
    
    else 
    
        rotmat=0.d0
        fck=gauss%fstatus(1)
        e0=gauss%fstatus(2)
        fcr=gauss%fstatus(3)
        Ec=gauss%pMaterial%variables(1)
        rotmat(:,1)=gauss%fstatus(4:6)
        rotmat(:,2)=gauss%fstatus(7:9)
        rotmat(:,3)=gauss%fstatus(10:12)    
        stressmat=0.d0
        strainmat=0.d0

        call oldvectomat(strainvec,strainmat)
        strainmat=matmul(transpose(rotmat),matmul(strainmat,rotmat))        
        call oldmattovec(strainmat,strainvec)
        
        do i = 1,3 
            if (strainvec(i).le.0) then
                call oldcalConcreteCompMod(-strainvec(i),fck,e0,Dmat(i,i))      
            else
                call oldcalConcreteTensMod(strainvec(i),fcr,Ec,Dmat(i,i))
            end if
        end do
                    
        Dmat(6,6)=(Dmat(1,1)*Dmat(2,2))/(Dmat(1,1)+Dmat(2,2))
        Dmat(4,4)=(Dmat(2,2)*Dmat(3,3))/(Dmat(2,2)+Dmat(3,3))
        Dmat(5,5)=(Dmat(3,3)*Dmat(1,1))/(Dmat(3,3)+Dmat(1,1))
            
        stressvec=matmul(Dmat,strainvec)
        call oldvectomat(stressvec,stressmat)
        stressmat=matmul(rotmat,matmul(stressmat,transpose(rotmat)))
        call oldmattovec (stressmat,stressvec)      
    end if
    gauss%stress(1:6)=stressvec
    end subroutine 
    
    subroutine oldcalConcreteFailCheck(gauss)
    
    !********************************************************************************************
    !       This Subroutine checks whether the concrete has cracked along the principle 
    !           direction or not
    !           Shanthanu 2017/13/3
    !********************************************************************************************   
    
    
        type( tGaussStatus ), intent(inout) :: gauss      !> material properties
        real(kind=kreal) :: strain(6)
        real(kind=kreal) :: eigval(3)     !< vector containing the eigvalues
        real(kind=kreal) :: princ(3,3)    !< matrix containing the three principal column vectors
        real(kind=kreal) :: fcr,failstrain,fck,E
        integer :: e1,e2,e3
        e2=2;
        strain=gauss%strain(1:6)
        E=gauss%pMaterial%variables(1)
        fck=gauss%pMaterial%variables(5)
        fcr=0.33*((fck/1.0e6)**0.5)*1.0e6
        failstrain=fcr/E
        
        call eigen3(strain,eigval,princ)
        
        gauss%istatus(1)=0
        
        e1=maxloc(eigval,1)
        e3=minloc(eigval,1)
        if (e1*e3==6)e2=1
        if (e1*e3==2)e2=3
        
        gauss%istatus(1)=1
        gauss%fstatus(1)=gauss%pMaterial%variables(5)
        gauss%fstatus(2)=gauss%pMaterial%variables(6)
        gauss%fstatus(3)=fcr
        gauss%fstatus(4:6)= princ(:,e1)
        gauss%fstatus(7:9)= princ(:,e2)
        gauss%fstatus(10:12)= princ(:,e3)
        
    end subroutine 
    
    subroutine oldcalConcreteCompMod(e,fck,e0,mod)
        real(kind=kreal), intent(in) :: e,fck,e0
        real(kind=kreal), intent(out) :: mod
        real(kind=kreal) :: fc
        
        fc= fck * ((2*e/e0)-((e/e0)**2))
        if (e==0.0 ) then
        mod=2*fck/e0
        else
        mod = fc/e
        end if

    end subroutine
    
    subroutine oldcalConcreteTensMod(e,fcr,Ec,mod)
        real(kind=kreal), intent(in) :: e,fcr,Ec
        real(kind=kreal), intent(out) :: mod
        real(kind=kreal) :: ft
        
        if (e.le.(fcr/Ec))then
            mod=Ec
        else        
            ft= fcr/(1+sqrt(200*e))
            mod=ft/e
        end if
    
    end subroutine 
    
    subroutine oldvectomat(vec,mat)
        real(kind=kreal), intent(in) :: vec(6)
        real(kind=kreal), intent(out) :: mat(3,3)
        
        mat(1,1)=vec(1); mat(2,2)=vec(2); mat(3,3)=vec(3)
        mat(1,2)=vec(4); mat(2,1)=mat(1,2)
        mat(2,3)=vec(5); mat(3,2)=mat(2,3)
        mat(3,1)=vec(6); mat(1,3)=mat(3,1)
                
    end subroutine
    
    subroutine oldmattovec(mat,vec)
        real(kind=kreal), intent(out) :: vec(6)
        real(kind=kreal), intent(in) :: mat(3,3)
        
        vec(1)=mat(1,1);vec(2)=mat(2,2);vec(3)=mat(3,3)
        vec(4)=mat(1,2);vec(5)=mat(2,3);vec(6)=mat(3,1)         
            
        
    end subroutine
    
    
        
! C***********************************************************************
      SUBROUTINE SEHOOP (EO,SIGY,SMAX,E,S,EH,IN,NUMHY)
! ! C***********************************************************************
! C
! C *** STRESS - STRAIN RELATIONSHIP FOR CONCRETE
! C     ORIGINAL PROGRAM = SUBROUTINE SEHOOP(CYCLIC RULE)
! C     THIS PROGRAM IS MODIFIED FOR MONOTONIC LOADING IN ORDER TO
! C     APPLY TO CALCULATE UNBLANCED STRESS < 1991.08.01 K.U and A.A >
! C
    implicit none
    real (kreal) :: EO,SIGY,SMAX,E,S,EH
    integer :: IN,NUMHY
    real (kreal) :: SMIN,E3
    

! C
! C     *** BI-LINEAR MODEL ***
! C
      SMIN=SMAX-2.0D0*SIGY/EO
! CCCCCCCCCCCCCCCC
      E3=EO/100.0D0
! c       E3=EO/10.0D0
! C       E3=EO/5.0D0
! CCCCCCCCCCCCCCCC
      IF(E.GT.SMIN.AND.E.LT.SMAX) GO TO 1
      IF(E.GE.SMAX)               GO TO 102
      IF(E.LE.SMIN)               GO TO 202
! C     GO TO (1,2,3),IN
! C
! C *** IN=1: WITHIN ELASTIC RANGE
! C
    1 CONTINUE
      IN=1
      EH=EO
      S=EO*E
      RETURN
! C
! C *** IN=2 : FROM ELASTIC TO PLASTIC
! C
  102 CONTINUE
      IF(IN.GE.2) GO TO 103
      IN=2
      EH=E3
      S=SIGY+E3*(E-SMAX)
      NUMHY=NUMHY+1
      RETURN

  202 CONTINUE
      IF(IN.GE.2) GO TO 203
      IN=2
      EH=E3
      S=-SIGY+E3*(E-SMIN)
      NUMHY=NUMHY+1
      RETURN
! C
! C *** IN=3 : WITHIN PLASTIC RANGE
! C
  103 CONTINUE
      IN=3
      EH=E3
      S=SIGY+E3*(E-SMAX)
      NUMHY=NUMHY+1
      RETURN

  203 CONTINUE
      IN=3
      EH=E3
      S=-SIGY+E3*(E-SMIN)
      NUMHY=NUMHY+1
      RETURN

      END

      
      ! ! *********************************************************************
      SUBROUTINE HOOPZ(   DEHOOP, TEHOOP,TSHOOP,EH,     &
                         OEHOOP, OSHOOP,                &
                         INHOOP,   EHOOP,               &
                         HCE,   HCS,   HRE,             &
                         HRS,   HXE,   HXS,             &
                         HTAE,  HTAS,                   &
                         HTBE,  HTBS,                   &
                         HCAE,  HCAS,                   &
                         HCBE,  HCBS,                   &
                         CK,   TK,                      &
                         HXNE,  HXPE,                   &
                         LOOP)                  
! **********************************************************************
         !DEHOOP:増分ひずみ
       !TEHOOP:ひずみ
       !TSHOOP:応力
       !EH:初期剛性
       !OEHOOP:前ステップでのひずみ
       !OSHOOP:前ステップでの応力
       !INHOOP:鉄筋指標
       !EHOOP:接線剛性
       !HCE:Ｃ点ひずみ
       !HCS:Ｃ点応力
       !HRE:Ｒ点ひずみ
       !HRS:Ｒ点応力
       !HXE:Ｘ点ひずみ
       !HXS:Ｘ点応力
       !HTAE:フープ鉄筋のトリリニアモデルでの引張域仮降伏ひずみ（第１折点）バイリニアモデルの場合には0.0
       !HTAS:フープ鉄筋のトリリニアモデルでの引張域仮降伏応力   　　　　　 バイリニアモデルの場合には0.0
       !HTBE:フープ鉄筋のトリリニアモデルでの引張域　降伏ひずみ（第２折点）
       !HTBS:フープ鉄筋のトリリニアモデルでの引張域　降伏応力
       !HCAE:フープ鉄筋のトリリニアモデルでの圧縮域仮降伏ひずみ（第１折点）バイリニアモデルの場合には0.0
       !HCAS:フープ鉄筋のトリリニアモデルでの圧縮域仮降伏応力   　　　　　 バイリニアモデルの場合には0.0
       !HCBE:フープ鉄筋のトリリニアモデルでの圧縮域　降伏ひずみ（第２折点）
       !HCBS:フープ鉄筋のトリリニアモデルでの圧縮域　降伏応力
       !CK:入力時フープ引張側鉄筋第２剛性
       !TK:入力時フープ圧縮側鉄筋第２剛性
       !HXNE:フープ鉄筋のトリリニアモデルでの指標番号31特別区域におけるε軸切片（ひずみ）
       !HXPE:フープ鉄筋のトリリニアモデルでの指標番号31特別区域におけるε軸切片（応力）                                                                       Ｘ点：除荷直線とＸ軸との交点
! *********************************************************************
! **** バイリニアモデルの繰返対応鉄筋履歴モデル                    ****
! ****                            PRODUCED BY TATSU SEP. 2002      ****
! *********************************************************************
        real (kreal) ::  DEHOOP,TEHOOP,TSHOOP,EH,OEHOOP,OSHOOP,EHOOP
        real (kreal) ::  HCE,HCS,HRE,HRS,HXE,HXS,HTAE,HTAS,HTBE,HTBS,HCAE,HCAS 
        real (kreal) ::  HCBE,HCBS,CK,TK,HXNE,HXPE
        integer :: INHOOP,LOOP
        
        real (kreal) :: E30
        
        

      ! CALL ROOT_CHECK(0000, 32)
! C EH=200000.0
      E30=EH/100. 
        !降伏後（履歴30）での剛性
      if(LOOP.eq.2)goto 2                                                     !バイリニアモデル対応ルーチンへ
      if(LOOP.eq.3)goto 3                                                     !トリリニアも出る対応ルーチンへ
! ****************************
! ***** バイリニアモデル *****
! ****************************
    2 IF(INHOOP.EQ.10)GOTO 210
      IF(INHOOP.EQ.30)GOTO 230
      IF(INHOOP.EQ.31)GOTO 231
      IF(INHOOP.EQ.32)GOTO 232
! *********************
! ***** INHOOP=10 *****
! *********************
  210 CONTINUE
      IF(  TEHOOP.LT.HCBE)GOTO 2306 
        !圧縮降伏しました
      IF(HTBE.LT.TEHOOP  )GOTO 2307 
        !引張降伏しました

      INHOOP=10 
        !弾性域(10)と判定されました
      TSHOOP =TEHOOP*EH
      EHOOP=EH
      GOTO 9000
! *********************
! ***** INHOOP=30 *****
! *********************
  230 CONTINUE
      IF(OEHOOP)2301,2304,2304

 2301 IF(DEHOOP)2306,2303,2303 
          !直前ステップが圧縮ひずみの場合
 2303 CALL RXC(  EH, E30, HTBE, HTBS,      &
              OEHOOP,OSHOOP,               &
              HRE,  HRS,  HXE,             &
              HXS,  HCE,  HCS)
        !除荷
      IF(HTBE.LT.TEHOOP  )GOTO 2307
      IF( HCE.LT.TEHOOP  )GOTO 2328
      IF(  TEHOOP.LE.HCE )GOTO 2317

 2304 IF(DEHOOP)2305,2307,2307 
          !直前ステップが引張ひずみの場合
 2305 CALL RXC(  EH, E30, HCBE, HCBS,   &
              OEHOOP,OSHOOP,            &
              HRE,  HRS,  HXE,          &
              HXS,  HCE,  HCS)
          !除荷 
      IF(  TEHOOP.LT.HCBE)GOTO 2306
      IF(  TEHOOP.LT.HCE )GOTO 2327
      IF( HCE.LE.TEHOOP  )GOTO 2317

 2306 INHOOP=30      
          !圧縮側30と判定されました
      TSHOOP =HCBS+(TEHOOP-HCBE)*E30
      EHOOP=E30
      GOTO 9000

 2307 INHOOP=30                       
          !引張側30と判定されました
      TSHOOP =HTBS+(TEHOOP-HTBE)*E30
      EHOOP=E30
      GOTO 9000
! *********************
! ***** INHOOP=31 *****
! *********************
  231 CONTINUE
      IF(HXE)2311,2314,2314

 2311 IF(DEHOOP)2312,2313,2313
 2312 IF(HCE.LT.HRE)THEN
        IF(  TEHOOP.LT.HCE )GOTO 2306  
        !下向き31でＸ切片が負、増分ひずみが負の場合
        IF( HCE.LE.TEHOOP  )GOTO 2317
      ELSE
        IF(  TEHOOP.LT.HRE )GOTO 2306      
        !上向き31でＸ切片が負、増分ひずみが負の場合
        IF( HRE.LE.TEHOOP  )GOTO 2317
      ENDIF
 2313 IF(HCE.LT.HRE)THEN
        IF(HTBE.LT.TEHOOP  )GOTO 2307    
        !下向き31でＸ切片が負、増分ひずみが正の場合
        IF( HRE.LT.TEHOOP  )GOTO 2328
        IF(  TEHOOP.LE.HRE )GOTO 2317
      ELSE
        IF(HTBE.LT.TEHOOP  )GOTO 2307   
        !上向き31でＸ切片が負、増分ひずみが正の場合
        IF( HCE.LT.TEHOOP  )GOTO 2328
        IF(  TEHOOP.LE.HCE )GOTO 2317
      ENDIF

 2314 IF(DEHOOP)2315,2316,2316
 2315 IF(HCE.LT.HRE)THEN
        IF(  TEHOOP.LT.HCBE)GOTO 2306  
         !下向き31でＸ切片が正、増分ひずみが正の場合
        IF(  TEHOOP.LT.HCE )GOTO 2327
        IF( HCE.LE.TEHOOP  )GOTO 2317
      ELSE
        IF(  TEHOOP.LT.HCBE)GOTO 2306 
         !下向き31でＸ切片が正、増分ひずみが正の場合
        IF(  TEHOOP.LT.HRE )GOTO 2327
        IF( HRE.LE.TEHOOP  )GOTO 2317
      ENDIF
 2316 IF(HCE.LT.HRE)THEN
        IF( HRE.LT.TEHOOP  )GOTO 2307 
         !下向き31でＸ切片が正、増分ひずみが負の場合
        IF(  TEHOOP.LE.HRE )GOTO 2317
      ELSE
        IF( HCE.LT.TEHOOP  )GOTO 2307    
         !上向き31でＸ切片が正、増分ひずみが負の場合
        IF(  TEHOOP.LE.HCE )GOTO 2317
      ENDIF

 2317 INHOOP=31 
           !31と判定されました
      TSHOOP =EH*TEHOOP+HRS-EH*HRE
      EHOOP=EH
      GOTO 9000
! *********************
! ***** INHOOP=32 *****
! *********************
  232 CONTINUE
      IF(OSHOOP)2321,2324,2324

 2321 IF(DEHOOP)2322,2323,2323 
          !圧縮応力域での32です
 2322 IF(  TEHOOP.LT.HCBE)GOTO 2306
      IF(HCBE.LE.TEHOOP  )GOTO 2327
 2323 CALL RXC(  EH, E30, HTBE, HTBS,  &
             OEHOOP,OSHOOP,          &
             HRE,  HRS,  HXE,        &
             HXS,  HCE,  HCS)
      IF( HCE.LT.TEHOOP  )GOTO 2307
      IF(  TEHOOP.LE.HCE )GOTO 2317

 2324 IF(DEHOOP)2325,2326,2326  
          !引張応力域での32です
 2325 CALL RXC(  EH, E30, HCBE, HCBS,    &
              OEHOOP,OSHOOP,            &
              HRE,  HRS,  HXE,          &
              HXS,  HCE,  HCS)
      IF(  TEHOOP.LT.HCE )GOTO 2306
      IF( HCE.LE.TEHOOP  )GOTO 2317
 2326 IF(HTBE.LT.TEHOOP  )GOTO 2307
      IF(  TEHOOP.LE.HTBE)GOTO 2328

 2327 INHOOP=32                       
          !圧縮域での32と判定されました
      TSHOOP =HCBS+(TEHOOP-HCBE)*E30
      EHOOP=E30
      GOTO 9000

 2328 INHOOP=32                               
          !引張域での32と判定されました
      TSHOOP =HTBS+(TEHOOP-HTBE)*E30
      EHOOP=E30
      GOTO 9000

! ****************************
! ***** トリリニアモデル *****
! ****************************
    3 if(INHOOP.eq.10)goto 310
      if(INHOOP.eq.20)goto 320
      if(INHOOP.eq.21)goto 321
      if(INHOOP.eq.22)goto 322
      if(INHOOP.eq.30)goto 330
      if(INHOOP.eq.31)goto 331
      if(INHOOP.eq.32)goto 332
! C      write(NFOUT,5000)
      stop
! *********************
! ***** INHOOP=10 *****
! *********************
  310 continue                                     !前回のステップは弾性域でした
      if(  TEHOOP.lt.HCBE)goto 3306                !圧縮降伏しました
      if(  TEHOOP.lt.HCAE)goto 3207                !圧縮第２剛性に入りました
      if(HTBE.lt.TEHOOP  )goto 3307                !引張降伏しました
      if(HTAE.lt.TEHOOP  )goto 3208                !引張第２剛性にはいました

      INHOOP=10                                    !弾性域と判定されました
      TSHOOP =TEHOOP*EH
      EHOOP=EH
      goto 9000
! *********************
! ***** INHOOP=20 *****
! *********************
  320 continue
      if(OEHOOP)3201,3204,3204
 3201 if(DEHOOP)3202,3203,3203                       !圧縮域です（DEHOOPが負）
 3202 if(  TEHOOP.lt.HCBE)goto 3306                  !圧縮降伏しました
      if(HCBE.le.TEHOOP  )goto 3207
 3203 call RXC(  EH,  TK, HTAE, HTAS,   &
                OEHOOP,OSHOOP,  &
                HRE,  HRS,  HXE,  HXS,  HCE,  HCS)
      if(HTBE.lt.TEHOOP  )goto 3307
      if(HTAE.lt.TEHOOP  )goto 3208
      if( HCE.lt.TEHOOP  )goto 3228
      if(  TEHOOP.le.HCE )goto 3217
 3204 if(DEHOOP)3205,3206,3206                        !引張域です（DEHOOPが正）
 3205 call RXC(  EH,  CK, HCAE, HCAS,   &
                OEHOOP,OSHOOP,          &
                HRE,  HRS,  HXE,  HXS,  HCE,  HCS)
      if(  TEHOOP.lt.HCBE)goto 3306
      if(  TEHOOP.lt.HCAE)goto 3207
      if(  TEHOOP.lt.HCE )goto 3227
      if( HCE.le.TEHOOP  )goto 3217
 3206 if(HTBE.lt.TEHOOP  )goto 3307
      if(  TEHOOP.le.HTBE)goto 3208

 3207 INHOOP=20                                        !圧縮域20と判定されました
      TSHOOP =HCAS+(TEHOOP-HCAE)*CK
      EHOOP=CK
      goto 9000

 3208 INHOOP=20                                        !引張域20と判定されました
      TSHOOP =HTAS+(TEHOOP-HTAE)*TK
      EHOOP=TK
      goto 9000
! *********************
! ***** INHOOP=21 *****
! *********************
  321 continue
      if(HXE)3211,3214,3214

 3211 if(DEHOOP)3212,3213,3213
 3212 if(HCE.lt.HRE)then
        if(  TEHOOP.lt.HCBE)goto 3306
        if(  TEHOOP.lt.HCE )goto 3207
        if( HCE.le.TEHOOP  )goto 3217
      else
        if(  TEHOOP.lt.HCBE)goto 3306
        if(  TEHOOP.lt.HRE )goto 3207
        if( HRE.le.TEHOOP  )goto 3217
      endif
 3213 if(HCE.lt.HRE)then
        if(HTBE.lt.TEHOOP  )goto 3307
        if(HTAE.lt.TEHOOP  )goto 3208
        if( HRE.lt.TEHOOP  )goto 3228
        if(  TEHOOP.le.HRE )goto 3217
      else
        if(HTBE.lt.TEHOOP  )goto 3307
        if(HTAE.lt.TEHOOP  )goto 3208
        if( HCE.lt.TEHOOP  )goto 3228
        if(  TEHOOP.le.HCE )goto 3217
      endif
 3214 if(DEHOOP)3215,3216,3216
 3215 if(HCE.lt.HRE)then
        if(  TEHOOP.lt.HCBE)goto 3306
        if(  TEHOOP.lt.HCAE)goto 3207
        if(  TEHOOP.lt.HCE )goto 3227
        if( HCE.le.TEHOOP  )goto 3217
      else
        if(  TEHOOP.lt.HCBE)goto 3306
        if(  TEHOOP.lt.HTBS)goto 3207
        if(  TEHOOP.lt.HRE )goto 3227
        if( HRE.le.TEHOOP  )goto 3217
      endif
 3216 if(HCE.lt.HRE)then
        if(HTBE.lt.TEHOOP  )goto 3307
        if( HRE.lt.TEHOOP  )goto 3208
        if(  TEHOOP.le.HRE )goto 3217
      else
        if(HTBE.lt.TEHOOP  )goto 3307
        if( HCE.lt.TEHOOP  )goto 3208
        if(  TEHOOP.le.HCE )goto 3217
      endif

 3217 INHOOP=21                                                        !21と判定されました
      TSHOOP =EH*TEHOOP+HRS-EH*HRE
      EHOOP=EH
      goto 9000
! *********************
! ***** INHOOP=22 *****
! *********************
  322 continue
      if(OSHOOP)3221,3224,3224

 3221 if(DEHOOP)3222,3223,3223
 3222 if(  TEHOOP.lt.HCBE)goto 3306
      if(  TEHOOP.lt.HTBS)goto 3207
      if(HCAE.le.TEHOOP  )goto 3227
 3223 call RXC(  EH,  TK, HTAE, HTAS,      &
              OEHOOP,OSHOOP,               &
                HRE,  HRS,  HXE,  HXS,  HCE,  HCS)
      if(HTBE.lt.TEHOOP  )goto 3307
      if( HCE.lt.TEHOOP  )goto 3208
      if(  TEHOOP.le.HCE )goto 3217
 3224 if(DEHOOP)3225,3226,3226
 3225 call RXC(  EH,  CK, HCAE, HCAS,      &
              OEHOOP,OSHOOP,               &
                HRE,  HRS,  HXE,  HXS,  HCE,  HCS)
      if(  TEHOOP.lt.HCBE)goto 3306
      if(  TEHOOP.lt.HCE )goto 3207
      if( HCE.le.TEHOOP  )goto 3217
 3226 if(HTBE.lt.TEHOOP  )goto 3307
      if(HTAE.lt.TEHOOP  )goto 3208
      if(  TEHOOP.le.HTAE)goto 3228

 3227 INHOOP=22                                                                   !圧縮応力側22と判定されました
      TSHOOP =HCAS+(TEHOOP-HCAE)*CK
      EHOOP=CK
      goto 9000

 3228 INHOOP=22                                                                     !引張応力側22と判定されました
      TSHOOP =HTAS+(TEHOOP-HTAE)*TK
      EHOOP=TK
      goto 9000
! *********************
! ***** INHOOP=30 *****
! *********************
  330 continue
      if(OEHOOP)3301,3304,3304

 3301 if(DEHOOP)3306,3303,3303                                                 !直前ステップが圧縮ひずみの場合
 3303 call RXC(  EH, E30, HTBE, HTBS,        &                                
              OEHOOP,OSHOOP,                 &                                
                HRE,  HRS,  HXE,  HXS,  HCE,  HCS)                            
      if(HTBE.lt.TEHOOP  )goto 3307                                           
      if( HCE.lt.TEHOOP  )goto 3328                                           
      if(  TEHOOP.lt.HCE )goto 3317                                           
                                                                              
 3304 if(DEHOOP)3305,3307,3307                                                 !直前ステップが引張ひずみの場合
 3305 call RXC(  EH, E30, HCBE, HCBS,      &                                  
              OEHOOP,OSHOOP,               &                                  
                HRE,  HRS,  HXE,  HXS,  HCE,  HCS)                            
      if(  TEHOOP.lt.HCBE)goto 3306                                           
      if(  TEHOOP.lt.HCE )goto 3327                                           
      if( HCE.lt.TEHOOP  )goto 3317                                           
                                                                              
 3306 INHOOP=30                                                                !圧縮側30と判定されました
      TSHOOP =HCBS+(TEHOOP-HCBE)*E30                                          
      EHOOP=E30                                                               
      goto 9000                                                               
                                                                              
 3307 INHOOP=30                                                                !引張側30と判定されました
      TSHOOP =HTBS+(TEHOOP-HTBE)*E30
      EHOOP=E30
      goto 9000
! *********************
! ***** INHOOP=31 *****
! *********************
  331 continue
      if( HXE.lt.HXNE)goto 3311
      if(HXPE.lt.HXE )goto 3314
      goto 3318                                                                !31特別区域内

 3311 if(DEHOOP)3312,3313,3313
 3312 if(HCE.lt.HRE)then
        if(  TEHOOP.lt.HCE )goto 3306                                                !下向き31でＸ切片が負、増分ひずみが負の場合
        if( HCE.le.TEHOOP  )goto 3317
      else
        if(  TEHOOP.lt.HRE )goto 3306                                                !上向き31でＸ切片が負、増分ひずみが負の場合
        if( HRE.le.TEHOOP  )goto 3317
      endif
 3313 if(HCE.lt.HRE)then
        if(HTBE.lt.TEHOOP  )goto 3307                                                !下向き31でＸ切片が負、増分ひずみが正の場合
        if( HRE.lt.TEHOOP  )goto 3328
        if(  TEHOOP.le.HRE )goto 3317
      else
        if(HTBE.lt.TEHOOP  )goto 3307                                                !上向き31でＸ切片が負、増分ひずみが正の場合
        if( HCE.lt.TEHOOP  )goto 3328
        if(  TEHOOP.le.HCE )goto 3317
      endif

 3314 if(DEHOOP)3315,3316,3316
 3315 if(HCE.lt.HRE)then
        if(  TEHOOP.lt.HCBE)goto 3306                                                !下向き31でＸ切片が正、増分ひずみが正の場合
        if(  TEHOOP.lt.HCE )goto 3327
        if( HCE.le.TEHOOP  )goto 3317
      else
        if(  TEHOOP.lt.HCBE)goto 3306                                                !下向き31でＸ切片が正、増分ひずみが正の場合
        if(  TEHOOP.lt.HRE )goto 3327
        if( HRE.le.TEHOOP  )goto 3317
      endif
 3316 if(HCE.lt.HRE)then
        if( HRE.lt.TEHOOP  )goto 3307                                                !下向き31でＸ切片が正、増分ひずみが負の場合
        if(  TEHOOP.le.HRE )goto 3317
      else
        if( HCE.lt.TEHOOP  )goto 3307                                                !上向き31でＸ切片が正、増分ひずみが負の場合
        if(  TEHOOP.le.HCE )goto 3317
      endif

 3317 INHOOP=31                                                                    !31と判定されました
      TSHOOP =EH*TEHOOP+HRS-EH*HRE
      EHOOP=EH
      goto 9000

 3318 if(DEHOOP)3319,3310,3310                                                     !特別区域内処理
 3319 if(HCE.lt.HRE)then
        if(  TEHOOP.lt.HCBE)goto 3306
        if(  TEHOOP.lt.HCE )goto 3327
        if( HCE.le.TEHOOP  )goto 3317
      else
        if(  TEHOOP.lt.HCBE)goto 3306
        if(  TEHOOP.lt.HRE )goto 3327
        if( HRE.le.TEHOOP  )goto 3317
      endif
 3310 if(HCE.lt.HRE)then
        if(HTBE.lt.TEHOOP  )goto 3307
        if( HRE.lt.TEHOOP  )goto 3328
        if(  TEHOOP.le.HRE )goto 3317
      else
        if(HTBE.lt.TEHOOP  )goto 3307
        if( HCE.lt.TEHOOP  )goto 3328
        if(  TEHOOP.le.HCE )goto 3317
      endif
! *********************
! ***** INHOOP=32 *****
! *********************
  332 continue
      if(OSHOOP)3321,3324,3324
 3321 if(DEHOOP)3322,3323,3323
 3322 if(  TEHOOP.lt.HCBE)goto 3306
      if(HCBE.le.TEHOOP  )goto 3327
 3323 call RXC(  EH, E30, HTBE, HTBS,         &
              OEHOOP,OSHOOP,                  &
                HRE,  HRS,  HXE,  HXS,  HCE,  HCS)
      if( HCE.lt.HTBE)then
        if(HTBE.lt.TEHOOP  )goto 3307
        if( HCE.lt.TEHOOP  )goto 3328
        if(  TEHOOP.le.HCE )goto 3317
      else
        if( HCE.lt.TEHOOP  )goto 3307
        if(  TEHOOP.le.HCE )goto 3317
      endif
 3324 if(DEHOOP)3325,3326,3326
 3325 call RXC(  EH, E30, HCBE, HCBS,       &
              OEHOOP,OSHOOP,                &
                HRE,  HRS,  HXE,  HXS,  HCE,  HCS)
      if(HCBE.lt.HCE )then
        if(  TEHOOP.lt.HCBE)goto 3306
        if(  TEHOOP.lt.HCE )goto 3327
        if( HCE.le.TEHOOP  )goto 3317
      else
        if(  TEHOOP.lt.HCE )goto 3306
        if( HCE.le.TEHOOP  )goto 3317
      endif
 3326 if(HTBE.lt.TEHOOP  )goto 3307
      if(  TEHOOP.le.HTBE)goto 3328

 3327 INHOOP=32                                                                    !圧縮側32と判定されました
      TSHOOP =HCBS+(TEHOOP-HCBE)*E30
      EHOOP=E30
      goto 9000

 3328 INHOOP=32                                                                    !引張側32と判定されました
      TSHOOP =HTBS+(TEHOOP-HTBE)*E30
      EHOOP=E30
      goto 9000
! **********************
! ***** Ｒ点の準備 *****   !バイリニアモデル、トリリニアモデル共通
! **********************
 9000 OEHOOP=TEHOOP    
          !除荷が発生してもいいように記録する
      OSHOOP=TSHOOP
      RETURN
      END
      
      
    ! **********************************************************************
      SUBROUTINE RXC(  EH, EHOOP, TE, TS,       &
              OEHOOP,OSHOOP,                    &
              HRE,  HRS,  HXE,                  &
              HXS,  HCE,  HCS)
! **********************************************************************
       !EH:初期剛性
       !EHOOP:第2剛性
       !TE:
       !TS:
       !OEHOOP:前スチE��プでのひずみ
       !OSHOOP:前スチE��プでの応力
       !HCE:�E�点ひずみ
       !HCS:�E�点応力
       !HRE:�E�点ひずみ
       !HRS:�E�点応力
       !HXE:�E�点ひずみ
       !HXS:�E�点応力
! **********************************************************************
! ****  こ�Eサブルーチンでの E はひずみではなく最新の鉁E��剛性       ****
! **********************************************************************
        implicit none
        real (kreal) :: EH, EHOOP, TE, TS,OEHOOP,OSHOOP
        real (kreal) :: HRE,HRS,HXE,HXS,HCE,HCS

      !CALL ROOT_CHECK(0000, 77)

      HRE=OEHOOP    
        !�E�点の座標記録
      HRS=OSHOOP
      HXE=HRE-HRS/EH                                                              !�E�刁E��算�E�E�Ｘ点�E�E
      HXS=0.0
      HCE=(EH*HRE-HRS-EHOOP*TE+TS)/(EH-EHOOP)                                             !�E�点の座標算�E
      HCS=EHOOP*(HCE-TE)+TS

      RETURN
      END  

      ! *********************************************************************************
      SUBROUTINE MMPTS(   DEHOOP,    TEHOOP,   TSHOOP,     &
                        EH, OEHOOP, OSHOOP,                &
                        INHOOP,   EHOOP,                   &
                    HRE,HRS,                               &
                        HXE,   HXS,                        &
                        HAE,   HAS,                        &
                        HTBE,  HTBS,                       &
                        HCBE,  HCBS,                       &
                         XMAX, XMIN,  XX2,INTC) 
         !�E�点�E�除荷点
! **********************************************************************
       !DEHOOP:増�Eひずみ
       !TEHOOP:ひずみ
       !TSHOOP:応力
       !EH:鉁E���E期剛性
       !OEHOOP:前スチE��プでのひずみ
       !OSHOOP:前スチE��プでの応力
       !INHOOP:鉁E��指樁E
       !EHOOP:接線剛性
       !HRE:�E�点ひずみ
       !HRS:�E�点応力
       !HXE:塑性ひずみ
       !HXS:塑性ひずみ点応力
       !HAE:�E�点から第2折線へ初期剛性で下ろしたとき�E交点座標�Eひずみ
       !HAS:�E�点から第2折線へ初期剛性で下ろしたとき�E交点座標�E応力
       !HTBE:�E��E�点ひずみ
       !HTBS:�E��E�点応力
       !HCBE:�E��E�点ひずみ
       !HCBS:�E��E�点応力
       !XMAX:正側最大塑性ひずみ
       !XMIN:負側最大塑性ひずみ
         !XX2:最大塑性ひずみ
! *********************************************************************************
! ****   バイリニアモチE��にCIAMPIが提案しぁEMODIFIED MENEGOTTO-PINTO MODELめE  ****
! ****   対応させた繰返し鉁E��履歴モチE��                                        ****
! ****                                                                         ****
! ****                                PRODUCED BY TATSUYA WATANABE OCT. 2002   ****
! ****                                IMPROVED BY AKIHISA HORIBE   OCT. 2003   ****
! *********************************************************************************
        implicit none
        real (kreal) :: DEHOOP,TEHOOP,TSHOOP,EH,OEHOOP,OSHOOP,EHOOP
        real (kreal) :: HRE,HRS,HXE,HXS,HAE,HAS,HTBE,HTBS,HCBE,HCBS,XMAX,XMIN,XX2
        real (kreal) :: E30
        
        integer :: INTC,INHOOP
      

! C       open(1116,FILE='t-s.csv')
! c       write(1116,*)DEHOOP
     
       E30=EH/100.
! *       E30=50027.0
! ****************************************
! ***** バイリニアモチE��でのみ使用可 *****
! ****************************************
      IF(INHOOP.EQ.10)GOTO 210           
        !弾性埁E
      IF(INHOOP.EQ.30)GOTO 230   
        !塑性埁E
      IF(INHOOP.EQ.31)GOTO 231   
        !除荷埁E �E�曲線！E
      IF(INHOOP.EQ.32)GOTO 232  
        !再載荷域（曲線！E
! *********************
! ***** INHOOP=10 *****
! *********************
  210 CONTINUE
      IF(  TEHOOP.LT.HCBE)GOTO 2306        
        !圧縮降伏しました
      IF(HTBE.LT.TEHOOP  )GOTO 2307                  
        !引張降伏しました

      INHOOP=10                                            
          !弾性埁E10)と判定されました
      TSHOOP =TEHOOP*EH
      EHOOP=EH
      GOTO 9000
! *********************
! ***** INHOOP=30 *****
! *********************
  230 CONTINUE
      IF(OEHOOP)2301,2304,2304

 2301 IF(DEHOOP)2306,2306,2303                        
          !直前スチE��プが圧縮ひずみの場吁E
! ************************************************************************
! ****   MAMI:今までに経験した最大塑性ひずみを求めるサブルーチン      ****
! ************************************************************************
 2303 CALL MAMI(   EH, OEHOOP, OSHOOP,   HRE,   HRS,    &  !圧縮降伏してぁE��増�Eが正のとぁE
                  HXE,  HXS, XMAX, XMIN,  XX2,          &
                  INHOOP )
! ************************************************************************
! ****  　   　 　CROSS_A:HAEとHASを求めるサブルーチン　    　　　　  ****
! ************************************************************************
      CALL CROSS_A(   EH,  E30,  HTBE,  HTBS,   HRE,   HRS,     &
                     HAE,   HAS )

      INTC=-1

      GOTO 2316

 2304 IF(DEHOOP)2305,2307,2307
! ************************************************************************
! ****   MAMI:今までに経験した最大塑性ひずみを求めるサブルーチン      ****
! ************************************************************************
 2305 CALL MAMI(   EH, OEHOOP, OSHOOP,   HRE,   HRS,            &
                  HXE,   HXS, XMAX, XMIN,  XX2,                 &
                  INHOOP )
        !引張降伏してぁE��増�Eが負のとぁE
! ************************************************************************
! ****  　　　    CROSS_A:HAEとHASを求めるサブルーチン　　　　　      ****
! ************************************************************************
      CALL CROSS_A(   EH,  E30,  HCBE,  HCBS,   HRE,   HRS,     &
                    HAE,   HAS )

      INTC=1

      GOTO 2316

 2306 INHOOP=30  
          !圧縮降伏してぁE��増�Eが負のとぁE
      TSHOOP =HCBS+(TEHOOP-HCBE)*E30      
    EHOOP=E30
      GOTO 9000

 2307 INHOOP=30    
          !引張降伏してぁE��増�Eが正のとぁE
      TSHOOP =HTBS+(TEHOOP-HTBE)*E30
      EHOOP=E30
      GOTO 9000
! *********************
! ***** INHOOP=31 *****
! *********************
  231 CONTINUE
      
    IF(INTC.EQ.-1) GOTO 2311  !圧縮先衁E
      IF(INTC.EQ. 1) GOTO 2312  !引張先衁E   

 2311 IF(DEHOOP)2315,2316,2316
 2312 IF(DEHOOP)2316,2316,2317

 2315 CONTINUE
! ************************************************************************
! ****   MAMI:今までに経験した最大塑性ひずみを求めるサブルーチン      ****
! ************************************************************************
      CALL MAMI(   EH, OEHOOP, OSHOOP,   HRE,   HRS,            &
                  HXE,   HXS, XMAX, XMIN,  XX2,             &
                  INHOOP )
! ************************************************************************
! ****  　    　　CROSS_A:HAEとHASを求めるサブルーチン　    　　　　  ****
! ************************************************************************
      CALL CROSS_A(   EH,  E30,  HCBE,  HCBS,   HRE,   HRS,     &
                     HAE,   HAS)
      GOTO 2325

 2316 CONTINUE
      INHOOP=31
! ************************************************************************
! ****  PRIMA:除荷・再載荷履歴中の応力と接線剛性を求めるサブルーチン  ****
! ************************************************************************
      CALL PRIMA(    TEHOOP,    TSHOOP,   EHOOP,   EH,  E30,    &
                  HRE,   HRS,   HAE,   HAS,                     &
                XX2 )
      GOTO 9000

 2317 CONTINUE
! ************************************************************************
! ****   MAMI:今までに経験した最大塑性ひずみを求めるサブルーチン      ****
! ************************************************************************
      CALL MAMI(   EH, OEHOOP, OSHOOP,   HRE,   HRS,            &
                  HXE,   HXS, XMAX, XMIN,  XX2,                 &
                  INHOOP )
! ************************************************************************
! ****  　    　　CROSS_A:HAEとHASを求めるサブルーチン　    　　　　  ****
! ************************************************************************
      CALL CROSS_A(   EH,  E30,  HTBE,  HTBS,   HRE,   HRS,     &
                     HAE,   HAS)
      GOTO 2325

! *********************
! ***** INHOOP=32 *****
! *********************
  232 CONTINUE

      IF(INTC.EQ.-1) GOTO 2321  !圧縮先衁E
      IF(INTC.EQ. 1) GOTO 2322  !引張先衁E

 2321 IF(DEHOOP)2325,2325,2326
 2322 IF(DEHOOP)2327,2325,2325

 2325 CONTINUE
      INHOOP=32
! ************************************************************************
! ****  PRIMA:除荷・再載荷履歴中の応力と接線剛性を求めるサブルーチン  ****
! ************************************************************************
      CALL PRIMA(   TEHOOP,    TSHOOP,   EHOOP,   EH,  E30,     &
                  HRE,   HRS,   HAE,   HAS,                     &
                 XX2 )
      GOTO 9000

 2326 CONTINUE
! ************************************************************************
! ****   MAMI:今までに経験した最大塑性ひずみを求めるサブルーチン      ****
! ************************************************************************
      CALL MAMI(   EH, OEHOOP, OSHOOP,   HRE,   HRS,            &
                  HXE,   HXS, XMAX, XMIN,  XX2,             &
                  INHOOP )
! ************************************************************************
! ****  　　　CROSS_A:HAEとHASを求めるサブルーチン　　　　　  ****
! ************************************************************************
      CALL CROSS_A(   EH,  E30,  HTBE,  HTBS,   HRE,   HRS,     &
                     HAE,   HAS)
      GOTO 2316

 2327 CONTINUE
! ************************************************************************
! ****   MAMI:今までに経験した最大塑性ひずみを求めるサブルーチン      ****
! ************************************************************************
      CALL MAMI(   EH, OEHOOP, OSHOOP,   HRE,   HRS,            &
                  HXE,   HXS, XMAX, XMIN,  XX2,             &
                  INHOOP )
! ************************************************************************
! ****  　　　CROSS_A:HAEとHASを求めるサブルーチン　　　　　  ****
! ************************************************************************
      CALL CROSS_A(   EH,  E30,  HCBE,  HCBS,   HRE,   HRS,     &
                     HAE,   HAS)
      GOTO 2316

! *********************
! *****    終亁E  *****
! *********************
 9000 OEHOOP=TEHOOP
      OSHOOP=TSHOOP

      RETURN
      END

! ************************************************************************
      SUBROUTINE MAMI(   EH, OEHOOP, OSHOOP,   HRE,   HRS,      &
                        HXE,   HXS, XMAX, XMIN,  XX2,           &
                        INHOOP )
! ************************************************************************
       !EH:初期剛性
       !OEHOOP:前スチE��プでのひずみ
       !OSHOOP:前スチE��プでの応力
       !HRE:�E�点ひずみ
       !HRS:�E�点応力
       !HXE:X点ひずみ
       !HXS:X点応力
       !XMAX:正側最大塑性ひずみ
       !XMIN:負側最大塑性ひずみ
       !XX2:最大塑性ひずみ
       !INHOOP:鉁E��指樁E
! ************************************************************************
! ****   MAMI:今までに経験した最大塑性ひずみを求めるサブルーチン      ****
! ************************************************************************
        implicit none
        real (kreal) :: EH,OEHOOP,OSHOOP,HRE,HRS,HXE,HXS,XMAX,XMIN,XX2
        integer :: INHOOP

        HRE=OEHOOP
        HRS=OSHOOP
        HXE=HRE-HRS/EH
        HXS=0.0
        IF(INHOOP.EQ.30) THEN
         IF(HXE.LT.0.0) THEN
           XMIN=HXE
           XMAX=0.0
         ELSE
           XMAX=HXE
           XMIN=0.0
         ENDIF
        ENDIF

        IF(HXE.LT.0.0.AND.HXE.LT.XMIN) XMIN=HXE
        IF(HXE.GE.0.0.AND.HXE.GT.XMAX) XMAX=HXE
        XX2=XMAX+ABS(XMIN)

        RETURN
        END
! ************************************************************************
      SUBROUTINE CROSS_A(   EH,  E30,  TCBE,  TCBS,   HRE,  HRS,        &
                           HAE,   HAS )
! ************************************************************************
       !EH:鉁E���E期剛性
       !E30:鉁E��第2剛性
       !TCBE:TCB点ひずみ
       !TCBS:TCB点応力
       !HRE:�E�点ひずみ
       !HRS:�E�点応力
       !HAE:�E�点から第2折線へ初期剛性で下ろしたとき�E交点座標�Eひずみ
       !HAS:�E�点から第2折線へ初期剛性で下ろしたとき�E交点座標�E応力
! ************************************************************************
! ****  　　　CROSS_A:HAEとHASを求めるサブルーチン　　　　　  ****
! ************************************************************************
        implicit none
        real (kreal) :: EH,E30,TCBE,TCBS,HRE,HRS,HAE,HAS
      
       HAE=(E30*TCBE-EH*HRE+HRS-TCBS)/(E30-EH)
       HAS=EH*(HAE-HRE)+HRS

       RETURN
       END
! ************************************************************************
      SUBROUTINE PRIMA(    TEHOOP,    TSHOOP,   EHOOP,   EH,  E30,      &
                          HRE, HRS, HAE,   HAS,                         &
                          XX2 )
! ************************************************************************
       !TEHOOP:ひずみ
       !TSHOOP:応力
       !EHOOP:接線剛性
       !EH:鉁E���E期剛性
       !E30:鉁E��第2剛性
       !HRE:�E�点ひずみ
       !HRS:�E�点応力
       !HAE:�E�点から第2折線へ初期剛性で下ろしたとき�E交点座標�Eひずみ
       !HAS:�E�点から第2折線へ初期剛性で下ろしたとき�E交点座標�E応力
         !XX2:最大塑性ひずみ
! ************************************************************************
! ****  PRIMA:除荷・再載荷履歴中の応力と接線剛性を求めるサブルーチン  ****
! ************************************************************************
      implicit none
      real (kreal) ::   TEHOOP,TSHOOP,EHOOP,EH,E30,HRE,HRS,HAE,HAS,XX2
      real (kreal) ::   HH,RR,EPY,RON1,RON2,PRES,ZZ,YY
      


        HH=E30/EH
        RR=20.0-18.5*XX2/(0.00015+XX2)

         IF(HAE.EQ.HRE)THEN
! C         WRITE(NFOUT,5002)
         STOP
         ENDIF
! C
! C         IF(ABS(HAE-HRE).LT.0.0002)THEN
! C            IF((HAE-HRE).GE.0.0.AND.(HAE-HRE).LT.0.0002) 
! C     *         HAE=HRE+0.0002
! C            IF((HAE-HRE).LE.0.0.AND.(HAE-HRE).GT.-0.0002) 
! C     *         HAE=HRE-0.0002
! C         ENDIF
! C
        EPY=(TEHOOP-HRE)/(HAE-HRE)
        RON1=1.0+(ABS(EPY))**RR
        RON2=RON1**(1.0/RR)
         IF(RON2.EQ.0.0)THEN
! C          WRITE(NFOUT,5002)
           STOP
         ENDIF

        PRES=HH*EPY+(1-HH)*EPY/RON2
        TSHOOP=PRES*(HAS-HRS)+HRS                                                      !これが履歴上�E応力

        IF(TEHOOP.GT.HRE)THEN
         ZZ=RR*(ABS((TEHOOP-HRE)/(HAE-HRE))**(RR-1))
        ELSE
         ZZ=-RR*(ABS((TEHOOP-HRE)/(HAE-HRE))**(RR-1))
        ENDIF

         YY=RON1**(-1.0/RR)-(TEHOOP-HRE)*ZZ             &
               *(RON1**(-1.0/RR-1.0))/RR
       EHOOP=ABS((HH/(HAE-HRE)+(1-HH)*YY/(HAE-HRE))     &
               *(HAS-HRS))  
           !これが履歴上�E接線剛性

       IF(EHOOP.GT.EH)EHOOP=EH

       RETURN
! *******************
! ***** FORMATS *****
! *******************
 ! 5002 FORMAT('     エラー発生！強制終亁E��まぁE                   &
          ! /'  エラー発生場所  SUBROUTINE PRIMA'              &
          ! /'RON2ぁEになりました。�E毁EにしてどぁE��るつもり�E�E)
       END
 
! ************************************************************************
  SUBROUTINE EFH( S_ , T_ , CBAI , TC_S , TC_T ,  &
        E1 , E_S , E_T , F_S , F_T , H_S , H_T ,  &
        HE_E )
! ************************************************************************
      real (kreal) :: S_ , T_ , CBAI , TC_S , TC_T ,  &
            E1 , E_S , E_T , F_S , F_T , H_S , H_T ,  &
            HE_E 
        real (kreal) :: E,Htc_S,X_S,X_T
! *** [ E点 ] ***
	E= CBAI * E1
	E_S= S_
	E_T= T_
! *** [ X点 ] ***
	X_S= E_S
	X_T= E_T
! *** [ F点 ] ***
	F_T= -0.18 * TC_T
	F_S=  E_S - ( E_T - F_T ) / E
! *** [ H点 ] ***
	H_T=  0.18 * TC_T
	Htc_S= H_T / E1
	H_S=  E_S - ( E_T - H_T ) / E1
	H_S=  (H_S + Htc_S) / 2.0
! *** [ EH剛性 ] ***
	HE_E= ( E_T - H_T ) / ( E_S - H_S )
	HE_E= ABS( HE_E )
	END
    
! C***********************************************************************
      SUBROUTINE TSBOND (S,T,E,E1B,E2P,E2N,         &
                        TCP,TCN,TYP,TYN,           &
                        IK,INC,INH,INBF,NUMBF) 
! ! C***********************************************************************
! C
! C *** BOND STRESS - SLIP DISPLACEMENT RELATIONSHIP
! C     ORIGINAL PROGRAM = SUBROUTINE TSBOND(CYCLIC RULE)
! C     THIS PROGRAM IS MODIFIED FOR MONOTONIC LOADING IN ORDER TO
! C     APPLY TO CALCULATE UNBLANCED STRESS
! C     < 1991.08.17 K.UCHIDA & A.AMEMIYA>
! C
! ************************************************************************
! *引数：	S	:	   	   
! *		T	:	   	   
! *		E	:	   	   
! *    E1B       : ボンドリンク要素の初期剛性係数
! *    E2N       : ボンドリンク要素のすべり方向が負側の2次剛性係数
! *    E2P       : ボンドリンク要素のすべり方向が正側の2次剛性係数	   	   
! *		TCP	:	   	   
! *		TCN	:	   	   
! *		TYP	:	   	   
! *		TYN	:	   	   
! *		IK	:	   	   
! *		INC	:周りのコンクリートが降伏したかの判断指標
! *		INH	:周りの鉄筋が降伏したかの判断指標	   	   
! *		INBF:	   	   
! *		NUMBF:			
! ************************************************************************
    integer :: IK,INC,INH,INBF,NUMBF
    real (kreal) ::S,T,E,E1B,E2P,E2N,TCP,TCN,TYP,TYN
    
    real (kreal)::SCN,SCN2,SCP,SCP2
! C
! C *** CALCULATE INDEPENDENT VALUE
! C
      IF(E1B.EQ.0.0D0) RETURN
      SCP=TCP/E1B
      SCP2=SCP+(TYP-TCP)/E2P
      SCN=TCN/E1B
      SCN2=SCN+(TYN-TCN)/E2N
      IF(INH.EQ.2) GO TO 12
      IF(INC.EQ.2) GO TO 20
      IF(IK.EQ.12) GO TO 12
      IF(IK.EQ.13) GO TO 13
      IF(IK.EQ.20) GO TO 20
      IF(S.GT.SCN.AND.S.LT.SCP)  GO TO 1
      IF(S.GE.SCP.AND.S.LT.SCP2) GO TO 101
      IF(S.GT.SCN2.AND.S.LE.SCN) GO TO 201
      IF(S.GE.SCP2)              GO TO 111
      IF(S.LE.SCN2)              GO TO 211
! C
! C *** IK=1: PRIMARY CURVE ****************************
! C
! C   * WITHIN ELASTIC RANGE ( TCN<T<TCP, E=E1B )
! C
    1 CONTINUE
      IK=1
      E=E1B
      T=E1B*S
      RETURN
! C
! C   * WITHIN PLASTIC RANGE ( TCP<T<TYP, E=E2P )
! C
  101 CONTINUE
      IK=1
      E=E2P
      T=TCP+E2P*(S-SCP)
      RETURN
! C
! C   * WITHIN PLASTIC RANGE ( TYN<T<TCN, E=E2N )
! C
  201 CONTINUE
      IK=1
      E=E2N
      T=TCN+E2N*(S-SCN)
      RETURN
! C
! C *** IK=11: BOND IS IN FAILURE
! C
! C   * ( S>SCP2 )
! C
  111 CONTINUE
      IF(IK.GE.11) GO TO 14
      IK=11
      E=0.0D0
      T=TYP/2.0D0
      INBF=3
      NUMBF=NUMBF+1
      RETURN
! C
! C   * ( S<SCN2 )
! C
  211 CONTINUE
      IF(IK.GE.11) GO TO 14
      IK=11
      E=0.0D0
      T=TYN/2.0D0
      INBF=3
      NUMBF=NUMBF+1
      RETURN
! C
! C *** IK=12: REINFORCEMENT IS YIELDED
! C
   12 CONTINUE
      IF(IK.GE.11) GO TO 13
      IK=12
      E=0.0D0
      IF(T.GE.0.0D0) THEN
       IF(T.GE.TYP) T=TYP/2.0D0
       IF(T.LT.TYP) T=T/2.0D0
      ELSE
       IF(T.LE.TYN) T=TYN/2.0D0
       IF(T.GT.TYN) T=T/2.0D0
      END IF
      INBF=3
      NUMBF=NUMBF+1
      RETURN
! C
! C *** IK=13: REINFORCEMENT WAS YIELDED
! C
   13 CONTINUE
      IK=13
      E=0.0D0
      INBF=3
      NUMBF=NUMBF+1
      RETURN
! C
! C *** IK=14: BOND WAS IN FAILUER
! C
   14 CONTINUE
      IK=14
      E=0.0D0
      INBF=3
      NUMBF=NUMBF+1
      RETURN
! C
! C *** IK=20: BOND WAS IN FAILURE DUE TO CONCRETE CRACKING
! C
   20 CONTINUE
      IK=20
      E=0.0D0
      T=0.0D0
      INBF=3
      NUMBF=NUMBF+1
      RETURN
! C
! C *** ERROR FORMAT
! C
 6001 FORMAT(' *** ERROR IN SUB. TSBON3 *** < IK =',I2,'>')
  
  END SUBROUTINE

! ************************************************************************
      SUBROUTINE TSBONDCYC (   IK,    S,    T,               &
                            E,  E1P,  E1N,  E2P,  E2N,       &  
                          TCP,  TCN,  TYP,  TYN,             &
                          INC, INC2,  INH, INBF,NUMBF,       &
     										 SOLD, TOLD,     &      !１ステップ前のすべりと応力
                          ASP,  ATP,  BSN,  BTN,             &       !A点、B点
                          CSP,  CTP,  DSN,  DTN,             &       !C点、D点
                          ESN,  ETN,  FSP,  FTP,             &       !E点、F点
                          GSP,  GTP,  HSN,  HTN,             &       !G点、H点
                          QSN,  QTN,  RSP,  RTP,             &       !I点、J点
                          XSP,  XTP,  YSN,  YTN,             &       !K点、L点
                         ITSP, ITSN,                         &       !IK=12で2度応力半減しないための指標
                         NDIK,BACK1,BACK2)                          !IK=1 OR 11を通過したという指標
! ************************************************************************
! **** BOND STRESS - SLIP DISPLACEMENT RELATIONSHIP                   ****
! **** ORIGINAL PROGRAM = SUBROUTINE TSBOND(CYCLIC RULE)              ****
! **** THIS PROGRAM IS MODIFIED FOR MONOTONIC LOADING IN ORDER TO     ****
! **** APPLY TO CALCULATE UNBLANCED STRESS < 1991.01.29 K.UCHIDA >    ****
! ************************************************************************
! ****   張愛暉の博士論文を一部参考にして、森田らのΤ-S 繰返しモデル  ****
! ****   を作成                                                       ****
! ****                                          <2002.1.10 T.SAKURAI> ****
! ************************************************************************
    integer :: IK,ITSP,ITSN,INC,INC2,INH,INBF,NUMBF,NDIK
    real (kreal) ::S,T,E,E1P,E1N,E2P,E2N,TCP,TCN,TYP,TYN,SOLD,TOLD
    real (kreal) ::ASP,ATP,BSN,BTN,CSP,CTP,DSN,DTN,ESN,ETN,FSP,FTP,GSP,GTP,HSN,HTN
    real (kreal) ::QSN,QTN,RSP,RTP,XSP,XTP,YSN,YTN,BACK1,BACK2
    
    real (kreal) ::AAA,AASP,AATP,BB1,BB2,BBSN,BBTN,DS,E1P_100,SCN,SCN2,SCP,SCP2
    
      CALL ROOT_CHECK(0000, 69)
! **** CALCULATE INDEPENDENT VALUE
! *CCCC
	E1P_100=E1P/100.
      IF(E1P.EQ.0.0) RETURN
      IF(E1N.EQ.0.0) RETURN
      SCP=TCP/E1P                                                               !正側第1折れ点すべり
      SCP2=SCP+(TYP-TCP)/E2P                                                    !正側第2折れ点すべり
      SCN=TCN/E1N                                                               !負側第1折れ点すべり
      SCN2=SCN+(TYN-TCN)/E2N                                                    !負側第2折れ点すべり

      DS=S-SOLD                                                                 !増分すべり

! **** JUDGEMENT 1

      IF(IK.EQ. 2) GOTO  2
      IF(IK.EQ. 3) GOTO  3
      IF(IK.EQ. 4) GOTO  4
      IF(IK.EQ. 5) GOTO  5
      IF(IK.EQ. 6) GOTO  6
      IF(IK.EQ. 7) GOTO  7
      IF(IK.EQ. 8) GOTO  8
      IF(IK.EQ.14) GOTO 14

      IF(IK.EQ.20) GOTO 20
      IF(IK.EQ.30) GOTO 30

! ***** JUDGEMENT 2

      IF(INC.EQ.2)  GOTO 20
      IF(INC2.EQ.2) GOTO 30
      IF(INH.EQ.2)  GOTO 12

! ***** JUDGEMENT 3

      IF(IK.EQ. 1) GOTO  1
      IF(IK.EQ.11) GOTO  1
      IF(IK.EQ.31) GOTO 31
      IF(IK.EQ.32) GOTO 32

! ***** IK=1: PRIMARY CURVE ****************************

! *   * WITHIN ELASTIC RANGE ( SCN<S<SCP, E=E1B )

    1 CONTINUE
         IF(NDIK.EQ.0)GOTO15
         IF(NDIK.EQ.1)THEN                                                     !IK= 1上の除荷の確認
           IF(IK.EQ.1.AND.DS.LT.0.0)THEN
             ASP=SOLD
             ATP=TOLD
             BACK1=ASP
             GOTO25
           ENDIF
           IF(IK.EQ.11.AND.DS.GT.0.0)THEN                                       !IK=11上の除荷の確認
             BSN=SOLD
             BTN=TOLD 
             BACK2=BSN
             GOTO35
           ENDIF
         ENDIF

         IF(S.GE.SCP) GOTO309
         IF(S.LE.SCN) GOTO319

   15 NDIK=10                                                                  !IK=1 OR 11を通過したという指標
         IF(S.GE.0.0)THEN
      IK=1
           E=E1P
           T=E*S
         ELSE
      IK=11
           E=E1N
           T=E*S
         ENDIF
         GOTO 9999
         RETURN

! *   * SCP<S<SCP2,E=E2P の包絡線上の領域

   31 CONTINUE

         IF(DS.LT.0.0)THEN                                                      !増分が負のとき
             ASP=SOLD                                                           !A点,正側包絡線からの除荷点ひずみ
             ATP=TOLD                                                           !A点,正側包絡線からの除荷点応力
             BACK1=ASP
             GOTO25
         ENDIF
  309    IF(S.GT.SCP2) GOTO141
  310 IK=31
           E=E2P
           T=TCP+E*(S-SCP)
         GOTO 9999
         RETURN

! *   * SCN2<S<SCN,E=E2N の包絡線上の領域

   32 CONTINUE

         IF(DS.GT.0.0)THEN                                                      !増分が正のとき
             BSN=SOLD                                                           !B点,負側包絡線からの除荷点ひずみ
             BTN=TOLD                                                           !B点,負側包絡線からの除荷点応力
             BACK2=BSN
             GOTO35
         ENDIF
  319    IF(S.LT.SCN2) GOTO142
  320 IK=32
           E=E2N
           T=TCN+E*(S-SCN)
         GOTO 9999
         RETURN

! *   * IK=2; IK=31で除荷されたとき

    2 CONTINUE
         IF(DS.GT.0.0)THEN                                                      !増分が正のとき
            IF(S.GT.SCP2) GOTO141
            IF(S.GT.ASP)THEN
               IF(ASP.GT.SCP2) GOTO141
               IF(ASP.LE.SCP2) GOTO310
            ENDIF
         ENDIF

   25 CONTINUE
         CTP=-0.18*ATP                                                          !C点
         CSP=ASP-(ATP-CTP)/(2*E1P)                                              !C点
         DTN=0.0                                                                !D点のクリア
         DSN=0.0                                                                !D点のクリア
          IF(S.LT.CSP) GOTO45

      IK=2
           E=2.0*E1P
           T=ATP-E*(ASP-S)
         GOTO 9999
         RETURN

! *   * IK=3; IK=32で除荷されたとき

    3 CONTINUE
         IF(DS.LT.0.0)THEN
            IF(S.LT.SCN2) GOTO142
            IF(S.LT.BSN) THEN
               IF(BSN.LT.SCN2) GOTO142
               IF(BSN.GE.SCN2) GOTO320
            ENDIF
         ENDIF

   35 CONTINUE
           DTN=-0.18*BTN                                                        !D点
           DSN=BSN-(BTN-DTN)/(2*E1N)                                            !D点
           CTP=0.0                                                              !C点のクリア
           CSP=0.0                                                              !C点のクリア
            IF(S.GT.DSN) GOTO55

      IK=3
           E=2.0*E1N
           T=BTN-E*(BSN-S)
         GOTO 9999
         RETURN

! *   * IK=4; IK=2から行き着く剛性ゼロ領域

    4 CONTINUE
         IF(DS.GT.0.0)THEN
           ESN=SOLD                                                             !E点
           ETN=TOLD                                                             !E点
           FTP=-0.18*ETN                                                        !F点
           FSP=ESN+(FTP-ETN)/(2*E1P)                                            !F点
           BACK2=ESN
           GOTO66
         ENDIF
      GOTO46

   45 CONTINUE
         IF(BSN.EQ.0.0)THEN
           IF(CTP.NE.0.0)HTN=CTP                                                !H点
           IF(ETN.NE.0.0)HTN=ETN                                                !H点
           HSN=HTN/E1N                                                          !H点
         ELSE
           IF(CSP.GT.0.0)THEN
             HSN=(3*BSN+BACK1)/4.0                                              !H点
           ELSE
             HSN=BSN/2.0                                                        !H点
           ENDIF
           IF(HSN.GT.0.0)HSN=0.0
           IF(CTP.NE.0.0)HTN=CTP                                                !H点
           IF(ETN.NE.0.0)HTN=ETN                                                !H点
         ENDIF

   46 CONTINUE

            IF(S.LT.HSN)THEN
              IF(BSN.EQ.0.0)GOTO15

              IF(BSN.GT.SCN)THEN
                 YTN=0.9*BTN                                                    !Y点
                 YSN=BSN-0.1*BTN/(2.0*E1N)                                      !Y点
                 GOTO85
              ENDIF

              IF(BSN.LE.SCN.AND.BSN.GE.SCN2)THEN
                 BB1=TCN-E2N*SCN
                 AAA=(0.9*BTN-HTN)/(BSN-HSN)
                 BB2=HTN-AAA*HSN

                 YSN=(BB1-BB2)/(AAA-E2N)                                        !Y点
                 YTN=E2N*YSN+BB1                                                !Y点
                   IF(YSN.LT.SCN2)THEN
                      YSN=SCN2
                      YTN=AAA*YSN+BB2
                   ENDIF
                 GOTO85
              ENDIF

              IF(BSN.LT.SCN2)THEN
                 AAA=(0.9*BTN-HTN)/(BSN-HSN)
                 BB2=HTN-AAA*HSN

                 YSN=(BTN-BB2)/AAA                                              !Y点
                 YTN=BTN                                                        !Y点
                 GOTO85
              ENDIF
            ENDIF

      IK=4
          E=E1P_100
          IF(CTP.NE.0.0)T=CTP
          IF(ETN.NE.0.0)T=ETN
         GOTO 9999
         RETURN

! *   * IK=5; IK=3から行き着く剛性ゼロ領域 

    5 CONTINUE
         IF(DS.LT.0.0)THEN
             FSP=SOLD                                                           !F点
             FTP=TOLD                                                           !F点
             ETN=-0.18*FTP                                                      !E点
             ESN=FSP+(ETN-FTP)/(2*E1P)                                          !E点
             BACK1=FSP
             GOTO65
         ENDIF
      GOTO56

   55 CONTINUE
         IF(ASP.EQ.0.0)THEN
             IF(DTN.NE.0.0)GTP=DTN                                              !G点
             IF(FTP.NE.0.0)GTP=FTP                                              !G点
             GSP=GTP/E1P                                                        !G点
         ELSE
             IF(DSN.LT.0.0)THEN
               GSP=(3.0*ASP+BACK2)/4.0
             ELSE
               GSP=ASP/2.0
             ENDIF
             IF(GSP.LT.0.0) GSP=0.0
             IF(DTN.NE.0.0)GTP=DTN                                              !G点
             IF(FTP.NE.0.0)GTP=FTP                                              !G点
         ENDIF

   56 CONTINUE

         IF(S.GT.GSP)THEN
            IF(ASP.EQ.0.0)GOTO15

            IF(ASP.LT.SCP)THEN
               XTP=0.9*ATP                                                      !X点
               XSP=ASP-0.1*ATP/(2.0*E1P)                                        !X点
               GOTO75
            ENDIF

            IF(ASP.GE.SCP.AND.ASP.LE.SCP2)THEN
               BB1=TCP-E2P*SCP
               AAA=(0.9*ATP-GTP)/(ASP-GSP)
               BB2=GTP-AAA*GSP

               XSP=(BB1-BB2)/(AAA-E2P)                                          !X点
               XTP=E2P*XSP+BB1                                                  !X点

                  IF(XSP.GT.SCP2)THEN
                      XSP=SCP2
                      XTP=AAA*XSP+BB2
                  ENDIF
               GOTO75
            ENDIF

            IF(ASP.GT.SCP2)THEN
               AAA=(0.9*ATP-GTP)/(ASP-GSP)
               BB2=GTP-AAA*GSP

               XSP=(ATP-BB2)/AAA                                                !X点
               XTP=ATP                                                          !X点
               GOTO75
            ENDIF
         ENDIF

      IK=5
         E=E1P_100
         IF(DTN.NE.0.0)T=DTN
         IF(FTP.NE.0.0)T=FTP
         GOTO 9999
         RETURN

! *   * IK=6; IK=4と5を結ぶ直線領域

    6 CONTINUE
          IF(DS) 65,66,66

   65 CONTINUE
          IF(S.LT.ESN.AND.FSP.LE.BSN/2.0)THEN
             HSN=ESN
             HTN=ETN
             GOTO46
          ENDIF
          IF(S.LT.ESN) GOTO45
          GOTO67
   66 CONTINUE
          IF(S.GT.FSP.AND.ESN.GE.ASP/2.0)THEN
             GSP=FSP
             GTP=FTP
             GOTO56
          ENDIF
          IF(S.GT.FSP) GOTO55
          GOTO67
   67 CONTINUE
      IK=6
          E=2*E1P
          T=FTP-E*(FSP-S)
          GOTO 9999
          RETURN

! *   * IK=7; IK=5から正側の包絡線に戻る領域

    7 CONTINUE
          IF(DS.LT.0.0)THEN
            RSP=SOLD                                                            !R点
            RTP=TOLD                                                            !R点
            CTP=-0.18*RTP                                                       !C点を新しく算出
            CSP=RSP-(RTP-CTP)/(2*E1P)                                           !C点を新しく算出
            BACK1=RSP

             IF(ASP.LT.SCP)THEN                                                 !以下、A点の入れ換え
               ASP=2.0*RSP-RTP/E1P
               ATP=E1P*ASP
               GOTO25
             ENDIF

             IF(ASP.LE.SCP2)THEN
               BB1=TCP-E2P*SCP
               BB2=RTP-2.0*E1P*RSP

               AASP=(BB1-BB2)/(2.0*E1P-E2P)
               AATP=E2P*AASP+BB1

                IF(AASP.GT.SCP2)THEN
                   ASP=AASP
                   ATP=TYP/2.0
                   GOTO25
                ENDIF

                IF(AASP.GT.ASP)THEN
                   ASP=AASP
                   ATP=AATP
                   GOTO25
                ENDIF
                   GOTO25
             ELSE
               BB2=RTP-2.0*E1P*RSP
               AASP=(ATP-BB2)/(2.0*E1P)
               AATP=ATP

                IF(AASP.GT.ASP)THEN
                   ASP=AASP
                   ATP=AATP
                   GOTO25
                ENDIF
                   GOTO25
             ENDIF
          ENDIF

   75 CONTINUE                                                                  !ここからＤＳがゼロ以上
            IF(S.GT.XSP)THEN
                IF(XSP.GE.SCP2)GOTO141
                IF(XSP.GE.SCP)THEN
                    GOTO309
                ELSE
                    GOTO 15
                ENDIF
            ENDIF

   76 CONTINUE
      IK=7
          E=(XTP-GTP)/(XSP-GSP) ; IF(E<=E1P_100) E=E1P_100
          T=GTP+(S-GSP)*E
          GOTO 9999
          RETURN

! *   * IK=8; IK=4から負側の包絡線に戻る領域

    8 CONTINUE
            IF(DS.GT.0.0)THEN
               QSN=SOLD                                                         !Q点
               QTN=TOLD                                                         !Q点
               DTN=-0.18*QTN                                                    !D点を新しく算出
               DSN=QSN-(QTN-DTN)/(2*E1P)                                        !D点を新しく算出
               BACK2=QSN

                 IF(BSN.GT.SCN)THEN                                             !以下、B点の入れ換え
                     BSN=2.0*QSN-QTN/E1N
                     BTN=E1N*BSN
                     GOTO35
                 ENDIF

                 IF(BSN.GE.SCN2)THEN
                     BB1=TCN-E2N*SCN
                     BB2=QTN-2.0*E1P*QSN

                     BBSN=(BB1-BB2)/(2*E1P-E2N)
                     BBTN=E2N*BBSN+BB1

                       IF(BBSN.LT.SCN2)THEN
                          BSN=BBSN
                          BTN=TYN/2.0
                          GOTO35
                       ENDIF

                       IF(BBSN.LT.BSN)THEN
                          BSN=BBSN
                          BTN=BBTN
                          GOTO35
                       ENDIF
                          GOTO35
                 ELSE
                     BB2=QTN-2.0*E1P*QSN
                     BBSN=(BTN-BB2)/(2*E1P)
                     BBTN=BTN

                       IF(BBSN.LT.BSN)THEN
                          BSN=BBSN
                          BTN=BBTN
                          GOTO35
                       ENDIF
                          GOTO35
                 ENDIF
            ENDIF

   85 CONTINUE                                                                  !ここからＤＳがゼロ以下
           IF(S.LT.YSN)THEN
              IF(YSN.LE.SCN2)GOTO142
              IF(YSN.LE.SCN)THEN
                  GOTO319
              ELSE
                  GOTO15
              ENDIF
           ENDIF

   86 CONTINUE
      IK=8
          E=(YTN-HTN)/(YSN-HSN) ; IF(E<=E1P_100) E=E1P_100
          T=HTN+E*(S-HSN)
          GOTO 9999
          RETURN

! ***** IK=12: REINFORCEMENT IS YIELDED
! ***** 鉄筋が降伏した時のみ通る履歴で、IK=12は存在しない。
! ***** 降伏応力の変更を行う。

   12 CONTINUE
           IF(IK.EQ.1)THEN
              IF(S.GE.0.0)THEN
                IF(ITSP.EQ.0) GOTO121
                 TYP=TOLD
                 TCP=TOLD
                 SCP=TCP/E1P
                 SCP2=SCP

                 ITSP=1

  121            IF(DS.GE.0.0) GOTO141
                 IF(DS.LT.0.0) GOTO 15
              ELSE
                IF(ITSN.EQ.0) GOTO122
                 TYN=TOLD
                 TCN=TOLD
                 SCN=TCN/E1N
                 SCN2=SCN

                 ITSN=1

  122            IF(DS.LE.0.0) GOTO142
                 IF(DS.GT.0.0) GOTO 15
              ENDIF
           ENDIF

           IF(IK.EQ.31)THEN
                IF(ITSP.EQ.0) GOTO123
                 TYP=TOLD
                 SCP2=SCP+(TYP-TCP)/E2P

                 ITSP=1

  123            IF(DS.GE.0.0) GOTO141
                 IF(DS.LT.0.0) GOTO 31
           ENDIF

           IF(IK.EQ.32)THEN
                IF(ITSN.EQ.0) GOTO124
                 TYN=TOLD
                 SCN2=SCN+(TYN-TCN)/E2N

                 ITSN=1

  124            IF(DS.LE.0.0) GOTO142
                 IF(DS.GT.0.0) GOTO 32
           ENDIF
          RETURN

! ***** IK=14: BOND WAS IN FAILUER

   14 CONTINUE
           IF(S.GT.0.0)THEN
              IF(DS.GE.0.0) GOTO141
              IF(DS.LT.0.0)THEN
                 ASP=SOLD
                 ATP=TOLD
                 GOTO25
              ENDIF
           ELSE
              IF(DS.LE.0.0) GOTO142
              IF(DS.GT.0.0)THEN
                 BSN=SOLD
                 BTN=TOLD
                 GOTO35
              ENDIF
           ENDIF

  141 CONTINUE                                                                  !Sがゼロ以上で、DSもゼロ以上のとき
      IK=14
          E=E1P_100
          T=TYP/2.0
          INBF=3
          NUMBF=NUMBF+1
          GOTO 9999

  142 CONTINUE                                                                  !Sがゼロ以下で、DSもゼロ以下のとき
      IK=14
          E=E1P_100
          T=TYN/2.0
          INBF=3
          NUMBF=NUMBF+1
          GOTO 9999
          RETURN

! * *** IK=20: BOND WAS IN FAILURE DUE TO DISCRETE CRACK

   20 CONTINUE
      IK=20
      E=E1P_100
      T=0.0
      INBF=3
      NUMBF=NUMBF+1
      RETURN

! * *** IK=30: BOND WAS IN FAILURE DUE TO SMEARED CRACK

   30 CONTINUE
      IK=30
       T=0.0
       E=E1P_100
      INBF=3
      NUMBF=NUMBF+1
      RETURN

 9999 CONTINUE
       SOLD=S
       TOLD=T
	 IF(E<=E1P_100) E=E1P_100
      RETURN

! * *** ERROR FORMAT

 6001 FORMAT(' *** ERROR IN SUB. TSBON2 *** < IK =',I2,'>')

      END


! ************************************************************************
! *    SDB       : 各ボンドリング要素の二つのバネのズレ
! *    DSL
! *    EMB       : 各ボンドリング要素の二つのバネの接線剛性
! *    TSB       : 各ボンドリング要素の二つのバネの持つ力
! *    E1B       : ボンドリンク要素の初期剛性係数
! *    E2P       : ボンドリンク要素のすべり方向が正側の2次剛性係数
! *    E2N       : ボンドリンク要素のすべり方向が負側の2次剛性係数
! *    TCP       : ボンドリンク要素のすべり方向が正側の第1折れ点
! *    TCN       : ボンドリンク要素のすべり方向が負側の第1折れ点
! *    TYP       : ボンドリンク要素のすべり方向が正側の最大点
! *    TYN       : ボンドリンク要素のすべり方向が負側の最大点
! *    IKB       : ボンドリング要素の応力状態の番号

! ************************************************************************
      SUBROUTINE TSBONDCYC_hj(IK ,                              &
     						DS , S , T , S_OLD , T_OLD ,		&								!IKB(M,I),SDB(M,I),TSB(M,I),
     						E , E1P , E1N , E2P , E2N , E_100 ,	&					!EMB(M,I),E1B(K,I),E1B(K,I),E2P(K,I),E2N(K,I),
     						TCP_T ,	TCP_S , TCN_T , TCN_S ,		&			!TCP(K,I),TCN(K,I),TYP(K,I),TYN(K,I)
     						TYP , TYNK ,                        &
     						EP_S , EP_T , EN_S , EN_T ,         &
     						XP_S , XP_T , XN_S , XN_T ,         &
     						FP_S , FP_T , FN_S , FN_T ,         &
     						HP_S , HP_T , HN_S , HN_T ,         &
     						R0_S , R0_T , R1_S , R1_T ,         &
     						HEP_E , HEN_E ,                     &
     						NST )
    integer :: IK,NST
    real (kreal) ::DS,S,T,S_OLD,T_OLD,E,E1P,E1N,E2P,E2N,E_100
    real (kreal) ::TCP_T,TCP_S,TCN_T,TCN_S,TYP,TYNK,EP_S,EP_T,EN_S,EN_T
    real (kreal) ::XP_S,XP_T,XN_S,XN_T,FP_S,FP_T,FN_S,FN_T,HP_S,HP_T,HN_S,HN_T
    real (kreal) ::R0_S,R0_T,R1_S,R1_T,HEP_E,HEN_E
! C
    real (kreal) ::CBAI,ALFA,XH
    
	XH= S * DS
	CBAI= 2.0
	ALFA= 0.18
! C
	IF(IK==10) GOTO 10

	IF(IK==12) GOTO 12
	IF(IK==32) GOTO 32

	IF(IK==13) GOTO 13
	IF(IK==33) GOTO 33

	IF(IK==14) GOTO 14
	IF(IK==34) GOTO 34

	IF(IK==15) GOTO 15
	IF(IK==35) GOTO 35

	IF(IK==16) GOTO 16
	IF(IK==36) GOTO 36

	IF(IK==17) GOTO 17
	IF(IK==37) GOTO 37
! C
! ***************
! **** IK=10  ****
! ***************
10	CONTINUE
	IF( S > TCP_S ) THEN ; IK=12 ; GOTO 1200 ; ENDIF
	IF( S < TCN_S ) THEN ; IK=32 ; GOTO 3200 ; ENDIF
! C
	IF( S > 0. ) THEN
		E= E1P
	ELSE
		E= E1N
	ENDIF
	T= S * E
	GOTO 9999
! C
! ***************
! **** IK=12  ****
! ***************
12	CONTINUE
	IF( DS < 0. ) THEN
		CALL EFH( S_OLD , T_OLD , CBAI , TCP_S , TCP_T ,    &
     			  E1P , EP_S , EP_T , FP_S , FP_T , HP_S , HP_T ,   &
     			  HEP_E)
		GOTO 13
	ENDIF
1200	CALL TE_20( E2P , TCP_S , TCP_T , S , E , T )
	GOTO 9999
! ***************
! **** IK=32  ****
! ***************
32	CONTINUE
	IF( DS > 0. ) THEN
		CALL EFH( S_OLD , T_OLD , CBAI , TCN_S , TCN_T ,            &
     			  E1N , EN_S , EN_T , FN_S , FN_T , HN_S , HN_T ,   &
     			  HEN_E)
		GOTO 33
	ENDIF
3200	CALL TE_20( E2N , TCN_S , TCN_T , S , E , T )
	GOTO 9999
! C
! ***************
! **** IK=13  ****
! ***************
13	IK= 13
! C
	IF( S > EP_S ) THEN ; IK= 12 ; GOTO 1200 ; ENDIF

	IF( S < EN_S ) THEN ; IK= 32 ; GOTO 3200 ; ENDIF
	IF( S < HN_S ) GOTO 35
	IF( S < FP_S ) GOTO 14
! C
	CALL TE_30( E1P , CBAI , EP_S , EP_T , S , E , T )
	GOTO 9999
! ***************
! **** IK=33  ****
! ***************
33	IK= 33
! C
	IF( S < EN_S ) THEN ; IK= 32 ; GOTO 3200 ; ENDIF

	IF( S > EP_S ) THEN ; IK= 12 ; GOTO 1200 ; ENDIF
	IF( S > HP_S ) GOTO 15
	IF( S > FN_S ) GOTO 34
! C
	CALL TE_30( E1N , CBAI , EN_S , EN_T , S , E , T )
	GOTO 9999
! ***************
! **** IK=14  ****
! ***************
14	IK= 14
! C
	IF( S > EP_S ) THEN ; IK= 12 ; GOTO 1200 ; ENDIF
	IF( S > FP_S ) GOTO 13

	IF( S < EN_S ) THEN ; IK= 32 ; GOTO 3200 ; ENDIF
	IF( S < HN_S ) GOTO 35
! C
	IF( DS > 0. ) THEN
		CALL R01( S_OLD,T_OLD,E1P,CBAI,ALFA , HP_S , HP_T ,EP_S ,EP_T,  &
     			R0_S , R0_T , R1_S , R1_T )
		GOTO 16
	ENDIF
! C
	E= 0.0
	T= FP_T
	GOTO 9999
! ***************
! **** IK=34  ****
! ***************
34	IK= 34
! C
	IF( S < EN_S ) THEN ; IK= 32 ; GOTO 3200 ; ENDIF
	IF( S < FN_S ) GOTO 33

	IF( S > EP_S ) THEN ; IK= 12 ; GOTO 1200 ; ENDIF
	IF( S > HP_S ) GOTO 15
! C
	IF( DS < 0. ) THEN
		CALL R01(S_OLD,T_OLD,E1N,CBAI, ALFA , HN_S , HN_T ,EN_S ,EN_T,  &
     			R0_S , R0_T , R1_S , R1_T )
		GOTO 36
	ENDIF
! C
	E= 0.0
	T= FN_T
	GOTO 9999
! C
! ***************
! **** IK=15  ****
! ***************
15	IK= 15
	IF( S > EP_S ) THEN ; IK= 12 ; GOTO 1200 ; ENDIF
! C
	IF( S < EN_S ) THEN ; IK= 32 ; GOTO 3200 ; ENDIF
	IF( S < FN_S ) GOTO 33
	IF( S < HP_S ) GOTO 34
! c
	IF( DS < 0. ) THEN
		CALL R01_5( S_OLD , T_OLD , E1P , CBAI , HN_T , &
     			R0_S , R0_T , R1_S , R1_T )
		GOTO 17
	ENDIF
! C
	E= HEP_E
	CALL TE_50(E , HP_S , HP_T , S , T )
	GOTO 9999
! ***************
! **** IK=35  ****
! ***************
35	IK= 35
! C
	IF( S < EN_S ) THEN ; IK= 32 ; GOTO 3200 ; ENDIF
! C
	IF( S > EP_S ) THEN ; IK= 12 ; GOTO 1200 ; ENDIF
	IF( S > FP_S ) GOTO 13
	IF( S > HN_S ) GOTO 14
! c
	IF( DS > 0. ) THEN
		CALL R01_5( S_OLD , T_OLD , E1N , CBAI , HP_T , &
     			R0_S , R0_T , R1_S , R1_T )
		GOTO 37
	ENDIF
! C
	E= HEN_E
	CALL TE_50(E , HN_S , HN_T , S , T )
	GOTO 9999

! ***************
! **** IK=16  ****
! ***************
16	IK= 16
! C
	IF( S < R0_S ) GOTO 14
	IF( S > R1_S ) THEN
		IF(R1_S > HP_S ) THEN
			GOTO 15
		ELSE
			GOTO 34
		ENDIF
	ENDIF
! C
	E= CBAI * E1P
	CALL TE_50(E , R0_S , R0_T , S , T)
	GOTO 9999
! ***************
! **** IK=36  ****
! ***************
36	IK= 36
! C
	IF( S > R0_S ) GOTO 34
	IF( S < R1_S ) THEN
		IF(R1_S < HN_S ) THEN
			GOTO 35
		ELSE
			GOTO 14
		ENDIF
	ENDIF
! C
	E= CBAI * E1N
	CALL TE_50(E , R0_S , R0_T , S , T)
	GOTO 9999
! ***************
! **** IK=17  ****
! ***************
17	IK= 17
! C
	IF( S > R0_S ) GOTO 15
	IF( S < R1_S ) GOTO 14
	E= CBAI * E1P
	CALL TE_50(E , R0_S , R0_T , S , T)
	GOTO 9999
! ***************
! **** IK=37  ****
! ***************
37	IK= 37
	IF( S < R0_S ) GOTO 35
	IF( S > R1_S ) GOTO 34
	E= CBAI * E1N
	CALL TE_50(E , R0_S , R0_T , S , T)
	GOTO 9999
! ***************
! ***************
! ***************
9999	CONTINUE
	S_OLD= S 
	T_OLD= T
	IF( E < E_100 ) E = E_100
	END
! ************************************************************************
      SUBROUTINE TE_20( E2 , TC_S , TC_T , S , E , T )
! ************************************************************************
    real (kreal) ::E2 , TC_S , TC_T , S , E , T
	E= E2
	T= ( S - TC_S ) * E + TC_T
	END
! ************************************************************************
      SUBROUTINE TE_30( E1 , CBAI , E_S , E_T , S , E , T )
! ************************************************************************
    real (kreal) ::E1 , CBAI , E_S , E_T , S , E , T 
	E= CBAI * E1
	T= ( S - E_S ) * E + E_T
	END
! ************************************************************************
      SUBROUTINE TE_50( E , H_S , H_T , S , T )
! ************************************************************************
    real (kreal) ::E , H_S , H_T , S , T
	T= ( S - H_S ) * E + H_T
	END
! ************************************************************************
      SUBROUTINE R01_5( S , T , E1 , CBAI , H_T ,   &
     			R0_S , R0_T , R1_S , R1_T )
! ************************************************************************
    real (kreal) ::S,T,E1,CBAI,H_T,R0_S,R0_T,R1_S,R1_T
    real (kreal) ::E
	E= E1 * CBAI
	R0_S= S
	R0_T= T

	R1_T= H_T
	R1_S= R0_S + ( R1_T - R0_T ) / E
	END
! ************************************************************************
      SUBROUTINE R01( S , T , E1 , CBAI , ALFA , H_S , H_T , E_S , E_T ,    &
     			R0_S , R0_T , R1_S , R1_T )
! ***********************************************************************
    real (kreal) ::S,T,E1,CBAI,ALFA,H_S,H_T,E_S,E_T,R0_S,R0_T,R1_S,R1_T
    real (kreal) ::E,XA1,YA1,XA2,YA2,XB1,YB1,XB2,YB2,X,Y
	E= E1 * CBAI
	R0_S= S
	R0_T= T
	R1_T= H_T
	R1_S= R0_S + ( R1_T - R0_T ) / E
	XA1= R0_S ; YA1= R0_T ; XA2= R1_S ; YA2= R1_T
	XB1= H_S  ; YB1= H_T  ; XB2= E_S  ; YB2= E_T
	CALL intersectXY( XA1,YA1 , XA2,YA2 , XB1,YB1 , XB2,YB2 , X,Y )
	IF( H_T >=0. ) THEN
		IF( Y > H_T ) THEN
			R1_S= X ; R1_T= Y
		ENDIF
	ELSE
		IF( Y < H_T ) THEN
			R1_S= X ; R1_T= Y
		ENDIF
	ENDIF
	END
! ************************************************************************
      SUBROUTINE intersectXY( XA1,YA1,XA2,YA2,XB1,YB1,XB2,YB2,X,Y )
! C
! C	点（XA1,YA1）と点（XA2,YA2）を結ぶ直線Aと、
! C	点（XB1,YB1）と点（XB2,YB2）を結ぶ直線Bの
! C	交点（X,Y）を求める。
! C
    real (kreal) ::XA1,YA1,XA2,YA2,XB1,YB1,XB2,YB2,X,Y
    
    real (kreal) ::A,A12,B,B12
	A12= ( YA1 - YA2 ) / ( XA1 - XA2 )
	B12= ( YB1 - YB2 ) / ( XB1 - XB2 )
	A= YB2 - YA2 + XB2 * B12 - XA2 * A12
	B= B12 - A12
	X= A / B
	Y= ( X - XA2 ) * A12 + YA2

	END SUBROUTINE




end module m_Concrete

! *************************************************************************
! *************************************************************************
! *************************************************************************
! *    AREAC     : 各コンクリート要素の面積
! *    AREAH     : 各種フープ要素材料の断面積
! *    AREAH     : フープ要素断面積
! *    AREAS     : 各面材鉄筋要素の面積
! *    ASHOOP    : 各フープ要素の体積
! *    COSH      : 各フープ要素各ステップのおいてX軸となす角度の余
! *    CRDSI     : クラックリンク要素の二つの節点の各ステップにおいての応力増分
! *    CRDSI     : クラックリンク要素2つの節点の応力増分
! *    CRPSI     : クラックリンク要素の二つの節点の主応力
! *    CRTSI     : クラックリンク要素の二つの節点の応力
! *    CRTSI     : クラックリンク要素の全体座標系(X-Y-Z)の各積分点の応力
! *    DCCCR     : ひび割れが入った後のDCC
! *    DCCL      : 1ステップ前のDCC
! *    DCC       : 座標or応力変換マトリックス　座標か応力か、または局所系か全体系かはそのときどきで異なる
! *    DDIS      : 1ステップごとの増分変位
! *    DEP       : コンクリート要素の各積分点の各ステップのひずみ増分
! *    DEPCS     : 積層の要素各積分点の二方向の鉄筋のひずみの増分
! *    DEPF      : 接合要素の各積分点におけるひずみ増分
! *    DEPS      : 鉄筋面材要素の各積分点の各ステップのひずみ増分
! *    DEPSL     : シェル要素の各積分点におけるひずみ増分
! *    DEPS      : 鉄骨要素の各積分点におけるひずみ増分
! *    DEP       : コンクリート要素の各積分点における全体座標系ひずみ増分
! *    DEU       : コンクリートの各主軸方向の増分等価一軸ひずみ
! *    DF        : 各ステップの荷重増分配列
! *    DFK       : 各節点の荷重係数
! *    DFK       : 荷重係数。増分荷重に荷重係数を掛けたものが実際の荷重になる。
! *    DF        : SOLVEで求めた各ステップの荷重増分配列
! *    DIVANG    : 現在使っていない
! *    DLOAD     : 1ステップごとの増分荷重
! *    DMCP1     : 積層要素のXY座標系の「D]マトリックス
! *    DMCP      : コンクリート要素のXY座標系の［D]マトリックス
! *    DMSP      : 鉄筋面材要素の各積分点のXY座標系の「Ｄ」
! *    DSI       : コンクリート要素の各積分点の各ステップの応力増分
! *    DSIF      : 接合要素の各積分点における応力増分
! *    DSIS      : 鉄筋面材要素の各積分点のX各ステップ応力増分
! *    DSISL     : シェル要素の各積分点における応力増分
! *    DSIS      : 鉄骨要素の各積分点における応力増分
! *    DSI       : コンクリート要素の各積分点における応力増分
! *    DU        : 各ステップの変位増分の解
! *    DUK       : 変位係数。増分変位に変位係数を掛けたものが実際の変位になる。
! *    DU        : SOLVEで求めた各ステップの変位増分配列
! *    E1B       : ボンドリンク要素の初期剛性係数
! *    E2N       : ボンドリンク要素のすべり方向が負側の2次剛性係数
! *    E2P       : ボンドリンク要素のすべり方向が正側の2次剛性係数
! *    EBU       : テンションスティフニングがなくなるときのひずみ
! *    EC*       : コンクリート要素の3軸応力状態を考慮した最大圧縮応力時主ひずみ
! *    EC*       : コンクリート要素の初期剛性　DMAT内ではEOに変化するので注意
! *    EC        : 各種類コンクリートの初期剛性
! *    ECB       : 各種類コンクリートの一軸圧縮強度時のひずみ
! *    ECB       : コンクリート要素の圧縮強度時ひずみ
! *    ECL       : 1ステップ前のEC
! *    ECR       : コンクリートの引張強度時ひずみ
! *    ECS       : 各積分点の材料定数(入力データの作成を参照）
! *    ECSEYA    : 積層要素の各積分点の降伏後の除荷点のひずみ
! *    ECSSYA    : 積層要素の各積分点の降伏後の除荷点の応力ly
! *    ECUUU     : コンクリート要素の0.3FC時の等価一軸主ひずみ
! *    ECUU      : コンクリート要素の0.65FC時の等価一軸主ひずみ
! *    ECU       : ECBと同じ意味　プログラムの階層によって呼び名が変わるのだ
! *    EDS       : 各種類鉄筋面材材料の剛性マトリックス
! *    EECS      : 積層の要素各積分点の二方向の鉄筋の接線剛性
! *    EESC      : コンクリート要素の3軸応力状態を考慮した最大主ひずみ　DMAT内ではECに変化するので注意が必要
! *    EH(EHOOP) : フープ要素の初期剛性
! *    EH1       : 各種フープ要素の除荷点の応力
! *    EH2       : 各種フープ要素の除荷点のひずみ
! *    EH        : 各種フープ要素材料の初期剛性
! *    EHOOP     : 各フープ要素の各ステップの接線剛性
! *    EINT      : 接合要素の初期剛性
! *    EK1       : クラックと平行方向のバネの接線剛性
! *    EK1E      : クラックリンク要素材料の閉じている場合のひび割れに平行方向のバネ剛性
! *    EK1E      : クラックリンク要素の初期剛性
! *    EK1Y      : クラックリンク要素材料の開いている場合のひび割れに平行方向のバネ剛性
! *    EK1Y      : クラックリンク要素のクラック発生後のダミー剛性
! *    EK2       : クラックと鉛直方向のバネの接線剛性
! *    EK2E      : クラックリンク要素材料の閉じている場合のひび割れに鉛直方向のバネ剛性
! *    EK2Y      : クラックリンク要素材料の開いている場合のひび割れに鉛直方向のバネ剛性
! *    EMB       : 各ボンドリング要素の二つのバネの接線剛性
! *    EO        : コンクリート要素の初期剛性
! *    EPCU      : コンクリートの圧壊が生じた後の収斂点ひずみ
! *    EPCUS     : コンクリートの圧壊が生じた後の収斂点ひずみ
! *    EPCUS     : コンクリート要素ピーク時以降の下降域の収斂点応力
! *    EPCU      : コンクリート要素ピーク時以降の下降域の収斂点ひずみ
! *    EPSYA     : 降伏したフープ要素の除荷点のひずみ
! *    ES        : 各種類鉄筋面材材料の初期剛性
! *    ESL       : シェル要素の初期剛性
! *    ES        : 鉄骨要素の初期剛性
! *    ETC       : コンクリート要素の各積分点の等価一軸方向の接線剛性
! *    EUC       : コンクリート要素の各積分点の等価一軸ひずみ
! *    EUL       : 1ステップ前のEU
! *    EUDL      : 1ステップ前の等価一軸増分ひずみ　by HJ
! *    EUU       : EU×105
! *    EU        : コンクリートの三軸応力状態を考慮した等価一軸ひずみ
! *    FC        : 各種類コンクリートの一軸圧縮強度
! *    FCCC      : コンクリートSCの下限値　0.3×FCの値
! *    FCC       : 0.6×FCの値
! *    FC        : コンクリート要素の一軸圧縮強度
! *    FINT      : 接合要素の最大強度
! *    FOECE     : 各ステージの全体荷重
! *    FT        : コンクリートの引張強度
! *    ICHK      : データチェックフラッグ　0 : チェックなし　1 : チェック有り
! *    ICIC(1)   : 圧縮強度上昇曲線　0 : SAENZ式　1 : FAFITIS-SHAH式
! *    ICIC(2)   : 圧縮応力ひずみ曲線ピーク以降のモデル　0 : 線形下降モデル
! *    ICIC(3)   : 圧縮強度低減係数モデル
! *    ICIC(4)   : テンションスティフニングモデル
! *    ICODK     : クラックリンク要素の開閉状態指標
! *    ICR(IICR) : ひび割れが入っているかどうかの指標
! *    IESMB     : 各ボンドリング要素の剛性マトリックスを作成するための節点関係マトリックス
! *    IESMC     : 各コンクリート要素の剛性マトリックスを作成するための節点関係マトリックス
! *    IESMH     : 各フープ要素の剛性マトリクックスを作成するの節点関係マトリックス
! *    IESMK     : 各クラックリンク要素の剛性マトリックスを作成するための節点関係マトリックス
! *    IESMS     : 各面材鉄筋要素の剛性マトリックスを作成するための節点関係マトリックス
! *    IIREK     : 主軸の回転に関する指標　SUBROUTINE DIRMで求める
! *    IKB       : ボンドリング要素数の応力状態
! *    INBF      : ボンドリング要素の残差力による収斂計算指標
! *    INCF      : コンクリート要素の残差力による収斂計算指標
! *    INCR      : クラックリンク要素による収斂計算指標
! *    INECS(in) : 積層要素の各積分点の二方向の鉄筋の応力状態(1  : 弾性状態である､２ : 降伏し始めている、３、降伏している
! *    INHOOP    : 各フープ要素の各ステップにおいての応力状態
! *    INK       : 現在使っていない
! *    INS       : 鉄筋面材要素の各積分点の応力状態
! *    INSY      : 鉄筋面材要素による収斂計算の判定指標
! *    IO        : 出力情報フラッグ
! *    IREK      : IIREKと同じ
! *    IROT      : せん断伝達モデルの指標（１ : 固定ひび割れモデル､　０ : 回転ひび割れモデル）
! *    IROT      : ひび割れモデル　0 : ひびわれ回転モデル　1 : ひび割れ座標系モデル（未対応）
! *    IR        : IIREKと同じ
! *    ITEM      : ウェブ。フロント法による剛性マトリックスを作成するための各節点の関係を示すマトリックス
! *    IUL       : コンクリート要素の各積分点の二つの等価一軸方向に応力状態
! *    IUL       : コンクリートの履歴指標　単調の場合は1、4、8、9のどれか
! *    JBO       : 読み込む境界条件の節点数
! *    LIMIT     : 総ステップ数
! *    LIMNIT    : 収斂回数
! *    LK        : 荷重種類　1 : 荷重制御　2 : 変位制御
! *    LW        : NOUTで指定した出力番号
! *    MATB      : ボンドの材料種類数
! *    MATB      : ボンドリンク要素材料種類数
! *    MATC      : コンクリートの材料種類数
! *    MATC      : コンクリート要素材料種類数
! *    MATF      : 接合要素材料種類数
! *    MATH      : フープ要素の材料種類数
! *    MATH      : フープ要素材料種類数
! *    MATK      : クラックリンク要素の材料種類数
! *    MATK      : クラックリンク要素材料種類数
! *    MATS      : 各種類鉄筋面材材料の材料種類数
! *    MATSL     : シェル要素材料種類数
! *    MATS      : 鉄骨要素材料種類数
! *    MBO       : 境界条件を与えた節点番号
! *    MCONB     : 各ボンドリング要素を構成する節点の番号
! *    MCONB     : ボンドリンクを構成する節点番号
! *    MCONC     : 各コンクリート要素を構成する節点番号
! *    MCONC     : コンクリートを構成する節点番号
! *    MCONF     : 接合要素を構成する節点番号
! *    MCONH     : 各フープ要素を構成する節点番号
! *    MCONH     : フープを構成する節点番号
! *    MCONK     : クラックリンク要素を構成する節点の番号
! *    MCONK     : クラックリンクを構成する節点番号
! *    MCONS     : 各面材鉄筋要素を構成する節点の番号
! *    MCONSL    : シェルを構成する節点番号
! *    MCONS     : 鉄骨を構成する節点番号
! *    MCR       : クラックリンクの両サイドのコンクリート要素番号
! *    MCRACK    : クラックリンクの両サイドのコンクリート要素番号
! *    MCRACKP   : クラックリンクの節点がそのコンクリート要素の中の順番
! *    NBAND2    : NBANDの２倍
! *    NBAND     : PARAMETER文でマトリックスの大きさを宣言する最大節点数
! *    NBOND     : PARAMETER文で配列の大きさを宣言する最大ボンドリンク要素数
! *    NBOND     : パラメータ文で指定したボンドリンク要素数
! *    NCONC     : PARAMETER文で配列の大きさを宣言する最大コンクリート数
! *    NCONC     : パラメータ文で指定したコンクリート要素数
! *    NCR       : 各ステップの支持状況（収斂計算の時に、節点に残差力があるため、NCRとNCRLが違う可能性がある）
! *    NCRACK    : PARAMETER文で配列の大きさを宣言する最大クラックリンク要素数
! *    NCRACK    : パラメータ文で指定したクラックリンク要素数
! *    NCRCP     : クラックリンクの両サイドのコンクリート要素のクラック要素に一番近い積分点の番号
! *    NCRL      : 各荷重種類の解析対象の実際の支持状況
! *    NCRL      : 境界条件　0 : フリー　1 : 荷重制御　2 : 変位制御
! *    NELMB     : ボンドリング要素数
! *    NELMB     : ボンドリンク要素数（2節点ボンドリンク要素）
! *    NELMC     : 解析対象のコンクリート要素数
! *    NELMC     : コンクリート要素数（8節点ソリッド要素）
! *    NELMF     : 接合要素（8節点接合要素）
! *    NELMH     : 解析対象のフープ要素数
! *    NELMH     : フープ要素数（2節点線要素）
! *    NELMK     : クラックリンク要素数（2節点クラックリンク要素）
! *    NELMK     : クラックリンク要素の数
! *    NELMS     : 解析対象の面材鉄筋要素数
! *    NELMSL    : シェル要素数（4節点シェル要素）
! *    NELMS     : 鉄骨要素数（8節点ソリッド要素）
! *    NEND      : 最終荷重ステージ数
! *    NENDL     : 各荷重種類の荷重ステージ数
! *    NENDL     : 1つの荷重種類における増分ステップの回数
! *    NEND      : 解析終了までの総ステップ数
! *    NEXTJ     : 複数の解析を連続して行うときの解析回数（未対応）
! *    NFC       : PARAMETER文で配列の大きさを宣言する最大荷重ステージ数
! *    NFC       : パラメータ文で指定している荷重種類ごとのステップ数の最大値
! *    NFD       : 継続JOB用出力装置番号
! *    NFP       : 図形データ（節点変位、要素応力、ひずみ）の出力装置番号。ファイル名はPLU1
! *    NFRE      : 節点の自由度数。シェルが関与していなければ3、関与していれば6になる
! *    NFR       : 継続JOB用入力装置番号
! *    NGA       : 8節点立方体要素のガウス積分点の数　2×2×2＝8
! *    NGL       : シェル要素のガウス積分点の数
! *    NHBOND    : ボンドリング要素と連結している鉄筋要素とクラックリンク要素の番号
! *    NHOOP     : パラメータ文で指定したフープ要素数
! *    NINT  
! *    NIT       : 現在計算中の収斂計算の回数
! *    NLD       : 現在計算している荷重種類番号
! *    NLK       : PARAMETER文で配列の大きさを宣言する最大荷重種類数
! *    NLK       : パラメータ文で指定した荷重種類数
! *    NLOA      : 解析対象の荷重種類数
! *    NLOA      : 荷重種類数
! *    NMAT      : PARAMETER文で配列の大きさを宣言する最大材料種類数
! *    NMAT      : パラメータ文で指定した材料種類数
! *    NNOD      : NNODの２倍
! *    NNOD      : PARAMETER文で配列の大きさを宣言する最大節点数
! *    NNOD      : パラメータ文で指定した総節点数
! *    NODE      : 解析対象の総節点数
! *    NODE      : 総節点数
! *    NORF      : 拘束節点の番号
! *    NOUT      : 解析結果を出力するステージの数、実際の番号はLWに入る
! *    NRFC      : PARAMETER文で配列の大きさを宣言する最大反力を求める節点数
! *    NSHEL     : パラメータ文で指定したシェル要素数
! *    NST       : 現在計算している荷重ステージ番号
! *    NSTEEL    : PARAMETER文で配列の大きさを宣言する最大鉄筋面材要素数
! *    NSTEEL    : パラメータ文で指定した鉄骨要素数
! *    NSTL      : 現在計算している荷重種類下のステージ番号
! *    NSYM      : PARAMETER文で配列の大きさを宣言する最大点対象節点数
! *    NSYM      : パラメータ文で指定した点対称の組み合わせの数
! *    NSZF
! *    NUMBF     : 降伏したボンドリンク要素の数
! *    NUMCF     : コンクリート要素のひずみ軟化域に入った積分点数（現在使っていない）
! *    NUMCR     : 開いているクラックリンク要素数
! *    NUMHY     : 降伏したフープ要素の本数
! *    NUMSY     : 降伏した積分点の数
! *    NW        : 解析結果を出力WRITE文の番号
! *    NWF       : 基本図形データ出力装置番号、ファイル名はPLU2
! *    NWRST     : 継続JOB用出力ステップの指定
! *    NW        : 解析結果を出力するWRIRE文の番号
! *    OAN       : コンクリート要素の各積分点の等価一軸の最初の方向
! *    OSM       : 解析対象の全体剛性マトリックス
! *    PRC       : コンクリート要素の各積分点の等価一軸方向のパアソン比
! *    PRC       : コンクリートのポアソン比
! *    PRS       : 各種類鉄筋面材材料のパアソン比
! *    PRSL      : シェル要素のポアソン比
! *    PRS       : 鉄骨要素のポアソン比
! *    PR        : PRCと同じ
! *    PSI       : コンクリート要素の各積分点の主応力
! *    PSICR     : コンクリートのひび割れ後の主応力
! *    PSIS      : 鉄筋面材要素の各積分点の主応力
! *    PSI       : コンクリートの主応力
! *    QB        : 各ボンドリング要素を通すフープ要素がＸ軸となす角度
! *    QK        : クラックリンク要素がX軸となす角度
! *    RESIDE    : 各ステップの収斂計算指標
! *    SC1～3    : コンクリートの三軸応力状態を考慮した最大主応力（破壊曲面との接点の値）
! *    SC2       : コンクリート要素の各積分点のの二つの等価一軸方向の最大圧縮強度
! *    SCT       : 各種類コンクリートの一軸引張強度
! *    SCT       : FTと同じ。コンクリートの引張強度
! *    SDB       : 各ボンドリング要素の二つのバネのズレ
! *    SECS      : 積層の要素各積分点の各ステップの降伏応力
! *    SFT       : コンクリート要素の各積分点の二軸応力状態下の引張強度
! *    SIGE      : 鉄筋面材要素の各積分点の二軸応力下降伏応力
! *    SIGH      : 各フープ要素の各ステップにおいてのX軸となす角度の正弦値
! *    SIGY1     : 各種フープ要素の除荷点の引張降伏し始める応力
! *    SIGY2     : 各種フープ要素の除荷点の圧縮降伏し始める応力
! *    SIGY      : 各種フープ要素材料の降伏応力
! *    SIGYA     : 降伏したフープ要素の除荷点の応力
! *    SIGY      : フープ要素の降伏応力
! *    SMAXH     : 各フープ要素各ステップにおいての降伏応力
! *    SN1～3    : コンクリートの主応力（PSI）
! *    SNCR      : コンクリートのひび割れ方向の主応力（固定ひび割れ使用時）
! *    SURB      : 各ボンドリング要素を構成する節点の負担面積
! *    SURK      : クラックリンク要素の負担面積
! *    TCN       : ボンドリンク要素のすべり方向が負側の第1折れ点
! *    TCP       : ボンドリンク要素のすべり方向が正側の第1折れ点
! *    TCW       : クラックの二つのバネのズレ
! *    TEP       : コンクリート要素の各積分点のひすみ
! *    TEPC      : コンクリート要素の各積分点のひずみ
! *    TEPCS     : 積層の要素各積分点の二方向の鉄筋のひずみ
! *    TEPC      : コンクリート要素の局所座標系(1-2-3)の各積分点のひずみ　（ε1,ε2,ε3,γ12,γ23,γ31）　（現ステップ）
! *    TEPF      : 接合要素の全体座標系(X-Y-Z)の各積分点のひずみ
! *    TEPLL     : 1ステップ前のTEP
! *    TEPL      : 1ステップ前のTEPC
! *    TEPS      : 鉄筋面材要素の各積分点のひずみ
! *    TEPSH     : シェル要素の全体座標系(X-Y-Z)の各積分点のひずみ
! *    TEPS      : 鉄骨要素の全体座標系(X-Y-Z)の各積分点のひずみ
! *    TEP       : コンクリート要素の全体座標系(X-Y-Z)の各積分点のひずみ　（εx,εy,εz,γxy,γyz,γzx）　（現ステップ）
! *    TF        : 各節点の全体荷重
! *    TF        : SOLVEで求めた各節点の全体荷重
! *    THC       : 各種類コンクリートの厚さ
! *    THKSL     : シェル要素の厚さ
! *    THS       : 各種類鉄筋面材材料の厚さ
! *    TSB       : 各ボンドリング要素の二つのバネの持つ力
! *    TSHOOP    : 各フープ要素の応力
! *    TSI       : コンクリート要素の各積分点の応力
! *    TSIF      : 接合要素の全体座標系(X-Y-Z)の各積分点の応力
! *    TSIS      : 鉄筋面材要素の各積分点の応力
! *    TSISL     : シェル要素の全体座標系(X-Y-Z)の各積分点の応力
! *    TSIS      : 鉄骨要素の全体座標系(X-Y-Z)の各積分点の応力
! *    TSI       : コンクリート要素の全体座標系(X-Y-Z)の各積分点の応力　(σx,σy,σz,τxy,τyz,τzx)
! *    TSK       : クラックの二つのバネが持っている力
! *    TSPCS     : 積層の要素各積分点の二方向の鉄筋の応力
! *    TU        : 各節点の全体変位
! *    TU        : SOLVEで求めた各節点の全体変位
! *    TYN       : ボンドリンク要素のすべり方向が負側の最大点
! *    TYP       : ボンドリンク要素のすべり方向が正側の最大点
! *    X         : 節点座標マトリックス
! *    YST       : 各種類鉄筋面材材料の降伏強度
! *    YSTL      : シェル要素の降伏応力
! *    YST       : 鉄骨要素の降伏応力
! *    angll     : 1ステップ前の主応力方向の角度
! *    e1p,e1n   : ボンドリンクの第1勾配（正負）
! *    e2p,e2n   : ボンドリンクの第2勾配（正負）
! *    ebu       : テンションスティフニング効果有効限界ひずみ
! *    edc       : コンクリート1要素の初期XY座標系の〔D〕マトリックス
! *    ehoop     : 各フープ要素の各ステップでの接線剛性
! *    ft        : 各コンクリートの一軸引張強度
! *    icic      : 圧縮応力上昇曲線オプション　0  : saenz式　1  : fafitis-shah式
! *    io        : 各データを各ステップで出力するorしない
! *    issc      : せん断剛性低減モデルオプション
! *    istor     : 積層要素および履歴（単調or繰返し）オプション
! *    jbo       : 境界条件を与える節点数
! *    lamda     : 圧縮強度低減係数モデルオプション
! *    limnit    : 最大ステップ数（収斂計算を含む）
! *    lk        : 荷重種類　変位or荷重制御
! *    loop      : 繰返し履歴用の鉄筋モデルオプション
! *    model42_52: 繰返し履歴上のiul=42,52のモデルオプション
! *    ncons     : 面材鉄筋要素を構成する節点の番号
! *    nfc       : paremeter文で配列の大きさを宣言する最大荷重ステップ数
! *    ngp       : コンクリート要素の積分点の数
! *    nline     : クラックライン数
! *    nout      : 指定詳細出力ステージの個数
! *    npsym     : 点対称節点の組の数
! *    numcp     : 総隅節点数
! *    sens      : 繰返し用の剛性復活点jを決めるための感度係数（現在は使っていない）
! *    tcp,tcn   : ボンドリンクの第1折れ点応力（正負）
! *    tehoop    : 各フープ要素のひずみ
! *    typ,tyn   : ボンドリンクの第2折れ点応力（正負）
! *    vs        : 繰返し履歴用でiul=4,5へ戻るための感度係数
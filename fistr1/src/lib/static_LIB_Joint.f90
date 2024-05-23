      MODULE m_static_LIB_Joint
!####################################################################
!
!     Gaku Hashimoto, The University of Tokyo, 2014/02/13
!
!####################################################################
      
      USE hecmw, ONLY : kint, kREAL
      USE elementInfo
      
!--------------------------------------------------------------------
      
      IMPLICIT NONE
      
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      
      CONTAINS
      
      
!####################################################################
      SUBROUTINE STF_Joint                           &
                 (etype, nn, ecoord, gausses, stiff, fstrDYNAMIC) 
!####################################################################
      
      USE m_fstr
      USE mMechGauss
      USE m_MatMatrix
      USE m_common_struct
      USE mMaterial
      
!-------------------------------------------------------------------
      
      INTEGER(KIND = kint), INTENT(IN) :: etype                   ! element type
      INTEGER(KIND = kint), INTENT(IN) :: nn                      ! number of elemental nodes
      REAL(KIND = kREAL), INTENT(IN) :: ecoord(3, nn)             ! coordinates of elemental nodes
      TYPE(tGaussStatus), INTENT(IN) :: gausses(:)                ! status of Gaussian qudrature points
      REAL(KIND = kREAL), INTENT(OUT) :: stiff(:,:)               ! stiffness matrix
      type ( fstr_dynamic ), optional, INTENT(IN) :: fstrDYNAMIC
      
!--------------------------------------------------------------------
      
      TYPE( tMaterial ), POINTER :: matl
      
      INTEGER(KIND = kint), PARAMETER :: ndof = 3
      INTEGER(KIND = kint) :: i, j
      INTEGER(KIND = kint) :: isize, jsize, ksize, lsize
      INTEGER(KIND = kint) :: mtype
      !INTEGER(KIND = kint) :: table(12)
      
      REAL(KIND = kREAL) :: k1, k2, k3 !, kk1, kk2, kk3
      REAL(KIND = kREAL) :: e11_hat(3), e12_hat(3), e13_hat(3)
      REAL(KIND = kREAL) :: e1_hat0(3), e2_hat0(3), e11_hat0(3), e12_hat0(3)
      REAL(KIND = kREAL) :: t(3, 3)
      REAL(KIND = kREAL) :: emat(12, 12)
      REAL(KIND = kREAL) :: area
      REAL(KIND = kREAL) :: lmat(8, 12)
      REAL(KIND = kREAL) :: IIw(8, 8)
      REAL(KIND = kREAL) :: IIw_dp(12, 12)
      REAL(KIND = kREAL) :: Vs, Vp, rho
      REAL(KIND = kREAL) :: rayleigh_alpha, rayleigh_beta, delta_t
      REAL(KIND = kREAL) :: newmark_beta, newmark_gamma
      REAL(KIND = kREAL) :: coeff_damp
      REAL(KIND = kREAL) :: lmat_dp(12, 12)
      
      TYPE elem_shape
      
      endtype
     
!--------------------------------------------------------------------
      
      stiff(:, :) = 0.0D0
      
!--------------------------------------------------------------------
      
      matl=>gausses(1)%pMaterial
      
      mtype = matl%mtype
      
!--------------------------------------------------------------------
      coeff_damp = 1.0d0
      
      ! Goodman joint element
      IF( ( mtype .EQ. GOODMAN ) .OR.              &
          ( mtype .EQ. SLIP_WEAKENING     )) THEN 
              
        call Set_stiff_Joint(gausses, mtype, matl, k1, k2, k3)

        call Set_LocalCood_Joint(nn, ecoord, t, area, e1_hat0, e2_hat0, &
                                e11_hat, e12_hat, e11_hat0, e12_hat0)
       
        lmat=0.0d0
        call Set_Bmat_Joint(nn, t, ecoord, lmat)

        IIw = 0.0D0
        call Set_Dmat_Joint(nn, ecoord, t, area, k1, k2, k3, e1_hat0, e2_hat0, &
                           e11_hat, e12_hat, e11_hat0, e12_hat0, IIw)
       
        emat = MATMUL( TRANSPOSE( lmat ) , MATMUL( IIW, lmat))
      
       !--------------------------------------------------------
             
      else IF( ( mtype .EQ. DASHPOD_ELEMENT ) ) THEN 
                
        call Set_param_Joint_dp(matl, Vp, Vs, rho, rayleigh_alpha, rayleigh_beta, &
                                   delta_t, newmark_beta, newmark_gamma, coeff_damp, fstrDYNAMIC)
        
        call Set_LocalCood_Joint(nn, ecoord, t, area, e1_hat0, e2_hat0, &
                                e11_hat, e12_hat, e11_hat0, e12_hat0)
       
        lmat_dp=0.0d0
        call Set_Bmat_Joint_dp(nn, t, ecoord, lmat_dp)
        
        IIw_dp = 0.0D0
        Call Set_Dmat_Joint_dp(area, rho, Vs, Vp, IIw_dp)       
        
        emat = MATMUL( TRANSPOSE( lmat_dp ) , MATMUL( IIW_dp, lmat_dp))
       
      END IF
      
      DO jsize = 1, 12
        
        DO isize = 1, 12
         
         ! ここでnewmark beta法の部分でつじつまが合うような係数coeff_dampをかける
          stiff(isize   , jsize   ) =  emat(isize, jsize) * coeff_damp
          stiff(isize+12, jsize   ) = -emat(isize, jsize) * coeff_damp
          stiff(isize   , jsize+12) = -emat(isize, jsize) * coeff_damp
          stiff(isize+12, jsize+12) =  emat(isize, jsize) * coeff_damp
         
        END DO
        
      END DO
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE STF_Joint
!####################################################################
      
      
!####################################################################
      SUBROUTINE Update_Joint                                  &
                 (icel, etype, nn, ecoord, u, du, qf, gausses, &
                  !strain_hat, stress_hat, strain_pla, strain_pla0, &
                  !strain_increment, &
                  !strain_init, stress_init, &
                  !flag_initset, &
                  iter, fstrDYNAMIC)                      
!####################################################################
      
      USE m_fstr
      USE mMaterial
      USE mMechGauss
      USE m_MatMatrix
      USE m_ElastoPlastic
      USE mHyperElastic
      USE m_utilities
      
!--------------------------------------------------------------------
      
      INTEGER(KIND = kint), INTENT(IN) :: etype            ! element type
      INTEGER(KIND = kint), INTENT(IN) :: nn               ! number of elemental nodes
      REAL(KIND = kREAL), INTENT(IN) :: ecoord(3, nn)      ! coordinates of elemental nodes
      REAL(KIND = kREAL), INTENT(IN) :: u(3, nn)           ! nodal dislplacements 
      REAL(KIND = kREAL), INTENT(IN) :: du(3, nn)          ! nodal displacement increment
      REAL(KIND = kREAL), INTENT(OUT) :: qf(nn*3)          ! internal force    
      TYPE(tGaussStatus), INTENT(INOUT) :: gausses(:)      ! status of Gaussian qudrature points
      type ( fstr_dynamic        ), INTENT(IN), optional :: fstrDYNAMIC
      !REAL(KIND = kREAL), INTENT(INOUT) :: strain_hat(8)   ! local strain
      !REAL(KIND = kREAL), INTENT(INOUT) :: stress_hat(8)   ! local stress
      !
      !REAL(KIND = kREAL), INTENT(INOUT) :: strain_pla0(8)   ! plastic strain(previous load step)
      !REAL(KIND = kREAL), INTENT(INOUT) :: strain_pla(8)   ! plastic strain
      !REAL(KIND = kREAL), INTENT(INOUT) :: strain_increment(8)   ! increment strain of present load step
      !
      !REAL(KIND = kREAL), INTENT(INOUT) :: strain_init(8)   ! plastic strain(previous load step)
      !REAL(KIND = kREAL), INTENT(INOUT) :: stress_init(8)   ! plastic strain
      !
      !INTEGER(KIND = kint), INTENT(INOUT) :: flag_initset   ! increment strain of present load step
      
      INTEGER(KIND = kint), INTENT(IN) :: iter            ! iter
      
!--------------------------------------------------------------------
      
      TYPE( tMaterial ), POINTER :: matl
      
      INTEGER(KIND = kint), PARAMETER :: ndof = 3
      INTEGER(KIND = kint) :: i, j
      INTEGER(KIND = kint) :: na, nb, nc, nd
      INTEGER(KIND = kint) :: isize, jsize, ksize, lsize
      INTEGER(KIND = kint) :: mtype
      !INTEGER(KIND = kint) :: table(12)
      
      REAL(KIND = kREAL) :: k1, k2, k3
      REAL(KIND = kREAL) :: phi, c
      REAL(KIND = kREAL) :: epsilon
      REAL(KIND = kREAL) :: e11_hat(3), e12_hat(3), e13_hat(3)
      REAL(KIND = kREAL) :: e1_hat0(3), e2_hat0(3), e11_hat0(3), e12_hat0(3)
      REAL(KIND = kREAL) :: t(3, 3), t_transpose(3, 3)

      !!<YM, KKE, 2014/03/05>
      REAL(KIND = kREAL) :: strain_ela(8)
      REAL(KIND = kREAL) :: ks0
      REAL(KIND = kREAL) :: strain_tildash, sig_c_temp, sig_tilde
      REAL(KIND = kREAL) :: strain_pla_recon(8)
      !!<YM, KKE, 2014/03/05 END>
      !!<YM, KKE, 2014/03/24>
      REAL(KIND = kREAL) :: strain_ela0(8)
      REAL(KIND = kREAL) :: coeff_a, coeff_b, coeff_c, root_beta
      !!<YM, KKE, 2014/03/24 END>
      REAL(KIND = kREAL) :: d_hat(8), dsig_hat(8)
      REAL(KIND = kREAL) :: emat(12, 12)
      REAL(KIND = kREAL) :: area
      REAL(KIND = kREAL) :: lmat(8, 12)
      REAL(KIND = kREAL) :: IIw(8, 8)
      REAL(KIND = kREAL) :: bmat(ndof*nn, 8)
      REAL(KIND = kREAL) :: pi
      INTEGER(KIND = kint) :: icel
      
      !!<YM, KKE, 2014/03/24>
      REAL(KIND = kREAL) :: ave_cood_z
      REAL(KIND = kREAL) :: sigma_init(3, 3)
      !!<YM, KKE, 2014/03/24 END>
      !!<YM, KKE, 2014/06/05>
      REAL(KIND = kREAL), POINTER :: stress_tildash0
      REAL(KIND = kREAL), POINTER :: strain_tildash0
      REAL(KIND = kREAL), POINTER :: stiff_coeff
      !!<YM, KKE, 2014/06/05 END>
      !!<YM, KKE, 2014/06/05>
      !REAL(KIND = kREAL) :: du_hour_hat_0, dv_hour_hat_0
      !!<YM, KKE, 2014/06/05 END>
      !!<YM, KKE, 2015/08/10>
      REAL(KIND = kREAL) :: sig_0, dc, stress_init_tildash, strain_yield, stress_abs
      !!<YM, KKE, 2015/08/10 END>
      
      INTEGER(KIND = kint), POINTER :: flag_initset   ! increment strain of present load step
      INTEGER(KIND = kint), POINTER :: istatus
      INTEGER(KIND = kint), POINTER :: istatus0
      
      REAL(KIND = kREAL), POINTER :: strain_hat(:)   ! local strain
      REAL(KIND = kREAL), POINTER :: stress_hat(:)   ! local stress
      
      REAL(KIND = kREAL), POINTER :: strain_pla0(:)   ! plastic strain(previous load step)
      REAL(KIND = kREAL), POINTER :: strain_pla(:)   ! plastic strain
      REAL(KIND = kREAL), POINTER :: strain_increment(:)   ! increment strain of present load step
      
      REAL(KIND = kREAL), POINTER :: strain_init(:)   ! plastic strain(previous load step)
      REAL(KIND = kREAL), POINTER :: stress_init(:)   ! plastic strain

!--------------------------------------------------------------------

      !  pointer設定 - 不要なものはコメントアウトする
      istatus          => gausses(1)%istatus(1)
      istatus0         => gausses(1)%istatus(2)
      flag_initset     => gausses(1)%istatus(3)
      
      strain_tildash0  => gausses(1)%fstatus( 1   )
      stress_tildash0  => gausses(1)%fstatus( 2   )
      stiff_coeff      => gausses(1)%fstatus( 3   )
      strain_hat       => gausses(1)%fstatus( 4:11)
      stress_hat       => gausses(1)%fstatus(12:19)
      strain_pla0      => gausses(1)%fstatus(20:27)
      strain_pla       => gausses(1)%fstatus(28:35)
      strain_increment => gausses(1)%fstatus(36:43)
      strain_init      => gausses(1)%fstatus(44:51)
      stress_init      => gausses(1)%fstatus(52:59)
      
!--------------------------------------------------------------------
      
      pi = 4.0D0*DATAN( 1.0D0 )
      
!--------------------------------------------------------------------
      
      qf(:) = 0.0D0
      d_hat(:) = 0.0d0
      dsig_hat(:) = 0.0d0
      
!--------------------------------------------------------------------
      
      matl=>gausses(1)%pMaterial
      
      mtype = matl%mtype
      
!--------------------------------------------------------------------
      
      ! Goodman joint element
      IF( ( mtype .EQ. GOODMAN ) .OR.              &
          ( mtype .EQ. SLIP_WEAKENING ) ) THEN 
                
       !--------------------------------------------------------
       
       ! 鉛直座標の平均値を取得
       ave_cood_z = 0.125d0 * (  ecoord(3, 1) + ecoord(3, 2) + ecoord(3, 3) + ecoord(3, 4) &
                               + ecoord(3, 5) + ecoord(3, 6) + ecoord(3, 7) + ecoord(3, 8) )
       
       !--------------------------------------------------------
       
       ! 載荷ステップ毎に更新する変数はここで処理する
       if(iter .eq. 1) then
         strain_pla0(:)=strain_pla(:)
         !gausses(1)%istatus0(1)=gausses(1)%istatus(1)
         istatus0=istatus
         strain_increment(:)=0.0d0
        
         !IF( mtype .EQ. SLIP_WEAKENING ) then
         !  write(997, '(A30, 2(I5, 1X), 100(E12.5, 1X) )')  matl%name, iter, istatus, ave_cood_z, stress_init(1:3), stress_hat(1:3), strain_init(1:3), strain_hat(1:3)
         !endif
       endif
       IF( mtype .EQ. SLIP_WEAKENING ) then
         !if(mod(iter, 10) == 0) then
         if(iter == 1) then
           ! write(997, '(A30, 2(I5, 1X), 100(E12.5, 1X) )')  matl%name, iter, istatus, ave_cood_z, stress_init(1:3), stress_hat(1:3), strain_init(1:3), strain_hat(1:3)
         endif
       endif

       !--------------------------------------------------------
       
       call Set_stiff_Joint(gausses, mtype, matl, k1, k2, k3)
       phi = matl%variables(M_GOODMAN_PHI)
       c   = matl%variables(M_GOODMAN_C) 
       epsilon = matl%variables(M_GOODMAN_EPS)
        
       !--------------------------------------------------------
       
       call Set_LocalCood_Joint(nn, ecoord, t, area, e1_hat0, e2_hat0, &
                                e11_hat, e12_hat, e11_hat0, e12_hat0)
       
       lmat=0.0d0
       call Set_Bmat_Joint(nn, t, ecoord, lmat)

       IIw = 0.0D0
       call Set_Dmat_Joint(nn, ecoord, t, area, k1, k2, k3, e1_hat0, e2_hat0, &
                           e11_hat, e12_hat, e11_hat0, e12_hat0, IIw)
                           
       ! 以下，単位面積当たりの判定を実施するため，並進成分はareaで除する
       ! 対応してqfの計算時にareaを乗じる
       do i=1, 3
         do j = 1, 3
           IIw(i, j) = IIw(i, j) / area
         enddo
       enddo
       
       !do i=4, 8
       !  do j = 4, 8
       !    IIw(i, j) = 0.0d0
       !  enddo
       !enddo
                 
       !--------------------------------------------------------

       DO i = 1, 3
         do j=1,8
           d_hat(j) = d_hat(j) - lmat(j, i    ) * du(i, 1) &
                               - lmat(j, i + 3) * du(i, 2) &
                               - lmat(j, i + 6) * du(i, 3) &
                               - lmat(j, i + 9) * du(i, 4) &
                               + lmat(j, i    ) * du(i, 5) &
                               + lmat(j, i + 3) * du(i, 6) &
                               + lmat(j, i + 6) * du(i, 7) &
                               + lmat(j, i + 9) * du(i, 8)
         end do
       end do

       dsig_hat = matmul(IIw, d_hat)
             
       !--------------------------------------------------------
       !!<YM, KKE, 2014/04/17>
       ! 初期応力の設定
       if(flag_initset .EQ. 0) then
       
         call Set_initStress_Joint(nn, matl, ecoord, t, &
                                     sigma_init, stress_init, strain_init, flag_initset)
                                     
       endif
       ! 初期応力の設定終了
       !!<YM, KKE, 2014/04/17 END>
       
       ! !Added by Shan 14/07/2017
       ! if(flag_initset .EQ. 2) then
         ! stress_init(1:3) = gausses(1)%stress(1:3)
         ! call Set_initStress_Joint_restart(nn,k1,k2,k3,t,stress_init,         &
                                            ! strain_init,flag_initset)              
       ! endif
             
       
       
       if(iter .eq. 1) then
       ! 接線剛性を計算するために、前回ステップの応力値を計算しておく
        strain_tildash0=DSQRT((strain_hat(1)+strain_init(1)-strain_pla0(1))*(strain_hat(1)+strain_init(1)-strain_pla0(1)) &
                             +(strain_hat(2)+strain_init(2)-strain_pla0(2))*(strain_hat(2)+strain_init(2)-strain_pla0(2)) )
        stress_tildash0=DSQRT((stress_hat(1)+stress_init(1))*(stress_hat(1)+stress_init(1)) &
                             +(stress_hat(2)+stress_init(2))*(stress_hat(2)+stress_init(2)) )
       endif
       
       strain_hat(1:8) = strain_hat(1:8) + d_hat(1:8) + strain_init(1:8)     
       
       ! 仮の要素応力の算出
       ! "仮"というのは、この後の接触・滑動判定で更新するから
       stress_hat(1:8) = stress_hat(1:8) + dsig_hat(1:8) + stress_init(1:8)      
       
       ! 各載荷ステップでの増分ひずみの更新
       strain_increment(1:8) = strain_increment(1:8) + d_hat(1:8)
       
       ! 通常のジョイント要素の場合
       IF( mtype .EQ. GOODMAN ) THEN
           ! 増分後の要素応力からモールクーロンの判定を実施する
           If(istatus .EQ. 3) then         !破壊時の判定
            strain_pla(1:2)=strain_hat(1:2)
            istatus=3
            stress_hat(1)=0.0d0
            stress_hat(2)=0.0d0
            stress_hat(6)=0.0d0
           else if(epsilon < strain_hat(3)) then         !剥離の判定
       
            stress_hat(:)=0.0d0                      !剥離時には応力を負担しない
            ! normal: 0, separation: 1, sliding: 2
            istatus=1
       
           else
        
            ! 弾性ひずみの計算
            strain_ela(1:6)=strain_hat(1:6)-strain_pla0(1:6)
            strain_tildash=DSQRT(strain_ela(1)*strain_ela(1)   &
                                +strain_ela(2)*strain_ela(2) )
            sig_c_temp=c-stress_hat(3)*DTAN( phi/180.0*pi )
            ks0=0.5d0*(matl%variables(M_GOODMAN_K1)+matl%variables(M_GOODMAN_K2))
            !write(998, '( 2(I5, 1X), 100(E12.5, 1X) )') icel, iter, c, phi, stress_hat(3), sig_c_temp
            If(istatus0 .EQ. 1) then
             !前回載荷ステップが剥離だった場合（再接触時）は塑性ひずみを計算する
             !この後の滑動・非滑動の判定はここで定義する塑性ひずみを用いて行う
             strain_pla_recon(1:6)=strain_hat(1:6)-(strain_hat(3)/strain_increment(3))*strain_increment(1:6)
         
             !弾性ひずみは更新した塑性ひずみから再計算する
             strain_pla(1:6) = strain_pla_recon(1:6)
             strain_ela(1:6)=strain_hat(1:6)-strain_pla_recon(1:6)
             strain_tildash=DSQRT(strain_ela(1)*strain_ela(1)   &
                                 +strain_ela(2)*strain_ela(2) )
            endif

            !! この部分の判定は、前回載荷ステップからのstrain_tildashで行う
            If(sig_c_temp<ks0*strain_tildash) then
             !滑動状態
         
             If(.true.) then
         
              !今回ステップの法線方向から計算
              coeff_b=strain_tildash-(sig_c_temp/ks0)
              root_beta=coeff_b/strain_tildash
              strain_pla(1:2)=strain_pla0(1:2)+root_beta*strain_ela(1:2)
          
             else
         
              !前回ステップの法線方向から計算
              !増分塑性ひずみの計算
              !（幾何関係より計算）
              strain_ela0(1:6)=(strain_hat(1:6)-strain_increment(1:6))-strain_pla0(1:6)
              coeff_a=strain_ela0(1)*strain_ela0(1)+strain_ela0(2)*strain_ela0(2)
              coeff_b=strain_ela0(1)*strain_increment(1)+strain_ela0(2)*strain_increment(2)
              coeff_c=strain_increment(1)*strain_increment(1)+strain_increment(2)*strain_increment(2) &
                        -(sig_c_temp/ks0)*(sig_c_temp/ks0)
          
              ! 解の公式（sig_c_tempの変化が大きいときに解を持たない時がある）
              !root_beta=(-coeff_b+(coeff_b*coeff_b-coeff_a*coeff_c)**0.5d0)/coeff_a
          
              ! 解を持たない場合の処理
              If(coeff_a .EQ. 0.0d0) then
          
               ! 初期ステップでせん断破壊する場合など
               ! strain_ela0 が0.0d0のため解の公式は使えない
               ! ⇔ 流れ則の法線方向が定義できない
               coeff_a=(strain_hat(1)*strain_hat(1)+strain_hat(2)*strain_hat(2))**0.5
               strain_pla(1:2)=(sig_c_temp/ks0)*strain_hat(1:2)/coeff_a
          
              else if (coeff_b*coeff_b-coeff_a*coeff_c .GT. 0.0d0) then
          
               root_beta=(-coeff_b+(coeff_b*coeff_b-coeff_a*coeff_c)**0.5d0)/coeff_a
               ! 面内せん断ひずみのみの更新
               strain_pla(1:2)=strain_pla0(1:2)+(1-root_beta)*strain_ela0(1:2)
          
              else
          
               ! とりあえずひずみ増分＝塑性ひずみ増分としておく
               strain_pla(1:2)=strain_pla0(1:2)+strain_increment(1:2)
          
              endif
          
             endif

             ! 面内回転は、塑性ひずみ＝全ひずみ（滑ると力を負担しない）
             strain_pla(6)=strain_hat(6)
             !strain_pla(6)=strain_pla0(6)+strain_increment(6)
         
             sig_tilde=sig_c_temp
             istatus=2
         
             ! 破壊の判定
             ! 塑性ひずみの大きさが弾性ひずみのx倍（仮）を超えたら破壊とする
             coeff_a = strain_pla(1)*strain_pla(1)+strain_pla(2)*strain_pla(2)
             if(coeff_a .GT. 1.1d10*(sig_c_temp/ks0)*(sig_c_temp/ks0)) then
             ! if(coeff_a .GT. 0.0d10*(sig_c_temp/ks0)*(sig_c_temp/ks0)) then
              strain_pla(1:2)=strain_hat(1:2)
              istatus=3
             endif
         
            else
             !非滑動状態
             sig_tilde=ks0*strain_tildash
             istatus=0
            endif
            !sig_tildeが設定されたので応力を再計算
            !この部分の応力は、更新後の塑性ひずみから計算しなければならない
            !弾性ひずみを更新する
            strain_ela(1:6)=strain_hat(1:6)-strain_pla(1:6)
            strain_tildash=DSQRT(strain_ela(1)*strain_ela(1)   &
                                +strain_ela(2)*strain_ela(2) )
            !弾性ひずみの方向余弦により分力  
            If(strain_tildash .EQ. 0.0d0) then
             stress_hat(1)=0.0d0
             stress_hat(2)=0.0d0
             stress_hat(6)=0.0d0
            else
             stress_hat(1)=sig_tilde*strain_ela(1)/strain_tildash
             stress_hat(2)=sig_tilde*strain_ela(2)/strain_tildash
             stress_hat(6)=IIw(6,6)*strain_ela(6)
            endif
        
           endif
       
           ! 接線剛性を計算するための係数を更新する
           !gausses(1)%stiff_coeff=(sig_tilde-gausses(1)%stress_tildash0)/(strain_tildash-gausses(1)%strain_tildash0)/ks0
           stiff_coeff=DABS(sig_tilde-stress_tildash0)/coeff_b/ks0
           if(stiff_coeff .LT. 0.366e+0) then
            stiff_coeff=0.366e+0
           elseif(stiff_coeff .GT. 0.1D0) then
            stiff_coeff=1.0d0
           endif
       
       Else if(mtype .EQ. SLIP_WEAKENING) THEN
           !                                               / if C < 0
           !matval(M_GOODMAN_K1)  = fval(1)
           !matval(M_GOODMAN_K2)  = fval(2)
           !matval(M_GOODMAN_K3)  = fval(3)
           !matval(M_GOODMAN_C)   = fval(4)     粘着力C
           !matval(M_PLCONST2)    = fval(5)     σ0(残留応力)
           !matval(M_PLCONST3)    = fval(6)     dc
           !matval(M_DASHPOD_VP)  = fval(7)     σinit_(1) / a
           !matval(M_PLCONST1)    = fval(8)     σinit_(2) / b
           !matval(M_DASHPOD_VS)  = fval(9)     σinit_(3) / c
           !matval(M_GOODMAN_PHI) = fval(10)    μ
           !matval(M_GOODMAN_EPS) = fval(11)               / d
           ! σinit_(3) = a*z + d
           ! σinit_(1) = b * σinit_(3)
           ! σinit_(2) = c * σinit_(3)
       
           ! C + μ･z
           c = dabs(matl%variables(M_GOODMAN_C)) + matl%variables(M_GOODMAN_PHI) * ave_cood_z
           
           sig_0 = dabs(matl%variables(M_PLCONST2))
           dc = dabs(matl%variables(M_PLCONST3))
           
           ks0=0.5d0*(matl%variables(M_GOODMAN_K1)+matl%variables(M_GOODMAN_K2))
           sig_tilde = DSQRT(stress_hat(1)*stress_hat(1)   &
                                +stress_hat(2)*stress_hat(2) )
           strain_tildash=DSQRT(strain_hat(1)*strain_hat(1)   &
                                +strain_hat(2)*strain_hat(2) )
           stress_init_tildash = DSQRT(stress_init(1)*stress_init(1)   &
                                +stress_init(2)*stress_init(2) )
           !strain_yield = (c - stress_init_tildash) / ks0
           strain_yield = c / ks0
           
           ! 始めから降伏している場合の処理
           if(stress_init_tildash > c)then
             istatus = 1    
             !strain_yield = 0.0d0
           endif
           
           
           ! 線形領域
           IF((istatus .eq. 0) .and. (strain_tildash < strain_yield)) then
             ! 線形領域は何もしなくてOK
               
           else if((istatus .ne. 0) .and. (strain_tildash < strain_yield)) then
             stress_hat(1) = (c / sig_tilde) * stress_hat(1)
             stress_hat(2) = (c / sig_tilde) * stress_hat(2)
           ! 負勾配の領域
           else if(strain_tildash < dc ) then
             stress_abs = c + (strain_tildash - strain_yield) * (sig_0 - c) / (dc - strain_yield)
             !stress_hat(1) = (stress_abs / sig_tilde) * stress_hat(1)
             !stress_hat(2) = (stress_abs / sig_tilde) * stress_hat(2)
             if(strain_tildash .ne. 0.0d0)then
                stress_hat(1) = (stress_abs / strain_tildash) * strain_hat(1)
                stress_hat(2) = (stress_abs / strain_tildash) * strain_hat(2)
             else
                stress_hat(1) = (stress_abs / sig_tilde) * stress_hat(1)
                stress_hat(2) = (stress_abs / sig_tilde) * stress_hat(2) 
             endif
             istatus = 1
             !write(999, '( 100(E12.5, 1X) )') c, strain_tildash, strain_yield, sig_0, dc, stress_abs, stress_hat(1), stress_hat(2)
           
           ! 一定の領域 
           else
             !stress_hat(1) = (sig_0 / sig_tilde) * stress_hat(1)
             !stress_hat(2) = (sig_0 / sig_tilde) * stress_hat(2)
             stress_hat(1) = (sig_0 / strain_tildash) * strain_hat(1)
             stress_hat(2) = (sig_0 / strain_tildash) * strain_hat(2)
             istatus = 2
           endif

           !write(998, '( 2(I5, 1X), 100(E12.5, 1X) )')  iter, gausses(1)%istatus(1), ave_cood_z, c, stress_init(1:3), stress_hat(1:3), strain_hat(1:3), strain_yield, strain_tildash
       ENDIF
           
       ! 結果出力用
       do j=1,3
            gausses(1)%stress(j) = stress_hat(j)
            gausses(1)%strain(j) = strain_hat(j)
       enddo
       
       ! 判定が終わったので戻す
       do j=1,8
        strain_hat(j) = strain_hat(j) - strain_init(j)
        stress_hat(j) = stress_hat(j) - stress_init(j)
       end do
       
       do j=4,6
            gausses(1)%stress(j) = stress_hat(j-3)
            gausses(1)%strain(j) = strain_pla(j-3)
       enddo
       ! gausses(1)%stress(4) = sig_c_temp
       ! gausses(1)%stress(5) = stress_init(1)
       ! gausses(1)%stress(6) = stress_init(3)
       gausses(1)%strain_out(1:6) = gausses(1)%strain(1:6)
       gausses(1)%stress_out(1:6) = gausses(1)%stress(1:6)
              
       bmat = 0.0D0
       
       DO jsize = 1, 8
         
        DO isize = 1, 12
         ! bmat = TRANSPOSE( lmat )
         bmat(   isize, jsize) = -lmat(jsize, isize)
         bmat(12+isize, jsize) =  lmat(jsize, isize)
         
        END DO
        
       END DO
       
       DO isize = 1, nn*3
         
        ! 並進の成分のみ面積を乗じる（bマトに応力→力の変換を含めていない）
        qf(isize) = qf(isize)+bmat(isize, 1)*stress_hat(1)*area
        qf(isize) = qf(isize)+bmat(isize, 2)*stress_hat(2)*area
        qf(isize) = qf(isize)+bmat(isize, 3)*stress_hat(3)*area
        qf(isize) = qf(isize)+bmat(isize, 4)*stress_hat(4)
        qf(isize) = qf(isize)+bmat(isize, 5)*stress_hat(5)
        qf(isize) = qf(isize)+bmat(isize, 6)*stress_hat(6)
        !qf(isize) = qf(isize)+bmat(isize, 7)*stress_hat(7)*area
        !qf(isize) = qf(isize)+bmat(isize, 8)*stress_hat(8)*area
       END DO
       
      END IF
      ! dashpodの処理は記載しないが、qfが初期化されているので何もしないとゼロになる
      
!--------------------------------------------------------------------      
      
      RETURN
      
!####################################################################
    END SUBROUTINE Update_Joint
!####################################################################

                  
!####################################################################
      SUBROUTINE Update_Joint_damp                                  &
                 (icel, etype, nn, ecoord, u, du, qf, gausses, &
                  acc, vec3, ddux, &
                  iter, fstrDYNAMIC)                      
!####################################################################
      
      USE m_fstr
      USE mMaterial
      USE mMechGauss
      USE m_MatMatrix
      USE m_ElastoPlastic
      USE mHyperElastic
      USE m_utilities
      
!--------------------------------------------------------------------
      
      INTEGER(KIND = kint), INTENT(IN) :: etype            ! element type
      INTEGER(KIND = kint), INTENT(IN) :: nn               ! number of elemental nodes
      REAL(KIND = kREAL), INTENT(IN) :: ecoord(3, nn)      ! coordinates of elemental nodes
      REAL(KIND = kREAL), INTENT(IN) :: u(3, nn)           ! nodal dislplacements 
      REAL(KIND = kREAL), INTENT(IN) :: du(3, nn)          ! nodal displacement increment
      REAL(KIND = kREAL), INTENT(OUT) :: qf(nn*3)          ! internal force    
      TYPE(tGaussStatus), INTENT(INOUT) :: gausses(:)      ! status of Gaussian qudrature points
      
      REAL(KIND = kREAL), INTENT(IN) :: acc(3, nn)
      REAL(KIND = kREAL), INTENT(IN) :: vec3(3, nn)
      REAL(KIND = kREAL), INTENT(IN) :: ddux(3, nn)
      type ( fstr_dynamic ), optional, INTENT(IN) :: fstrDYNAMIC
      
      INTEGER(KIND = kint), INTENT(IN) :: iter            ! iter
      INTEGER(KIND = kint) :: icel
      
!--------------------------------------------------------------------
      
      TYPE( tMaterial ), POINTER :: matl
      
      INTEGER(KIND = kint), PARAMETER :: ndof = 3
      INTEGER(KIND = kint) :: i, j
      INTEGER(KIND = kint) :: na, nb, nc, nd
      INTEGER(KIND = kint) :: isize, jsize, ksize, lsize
      INTEGER(KIND = kint) :: mtype
      !INTEGER(KIND = kint) :: table(12)
      
      REAL(KIND = kREAL) :: k1, k2, k3 !, kk1, kk2, kk3
      !REAL(KIND = kREAL) :: l, m
      REAL(KIND = kREAL) :: e11_hat(3), e12_hat(3), e13_hat(3)
      REAL(KIND = kREAL) :: e1_hat0(3), e2_hat0(3), e11_hat0(3), e12_hat0(3)
      REAL(KIND = kREAL) :: t(3, 3)
      REAL(KIND = kREAL) :: emat(12, 12)
      REAL(KIND = kREAL) :: area
      REAL(KIND = kREAL) :: IIw_dp(12, 12)
      REAL(KIND = kREAL) :: Vs, Vp, rho
      REAL(KIND = kREAL) :: rayleigh_alpha, rayleigh_beta, delta_t
      REAL(KIND = kREAL) :: newmark_beta, newmark_gamma
      REAL(KIND = kREAL) :: param_a3, param_b3, coeff_damp
      REAL(KIND = kREAL) :: lmat_dp(12, 12)
      REAL(KIND = kREAL) :: stiff(ndof*nn, ndof*nn)     
      
!--------------------------------------------------------------------
      
      ! 載荷ステップ毎に更新する変数はここで処理する
      if(iter .eq. 1) then

      endif
      
!--------------------------------------------------------------------
     
      qf(1:nn*3) = 0.0D0
      
!--------------------------------------------------------------------
      
      matl=>gausses(1)%pMaterial
      
      mtype = matl%mtype
      
!--------------------------------------------------------------------
      
      ! dashpod element
      IF( ( mtype .EQ. DASHPOD_ELEMENT ) ) THEN 
         
         call Set_param_Joint_dp(matl, Vp, Vs, rho, rayleigh_alpha, rayleigh_beta, &
                                 delta_t, newmark_beta, newmark_gamma, coeff_damp, fstrDYNAMIC) 
       
       call Set_LocalCood_Joint(nn, ecoord, t, area, e1_hat0, e2_hat0, &
                                e11_hat, e12_hat, e11_hat0, e12_hat0)
                                
       lmat_dp=0.0d0
       call Set_Bmat_Joint_dp(nn, t, ecoord, lmat_dp)
       
       IIw_dp = 0.0D0
       Call Set_Dmat_Joint_dp(area, rho, Vs, Vp, IIw_dp)    
             
       emat = MATMUL( TRANSPOSE( lmat_dp ) , MATMUL( IIW_dp, lmat_dp))
             
       !--------------------------------------------------------
      
       DO jsize = 1, 12
        
        DO isize = 1, 12
         
         stiff(isize   , jsize   ) =  emat(isize, jsize)
         stiff(isize+12, jsize   ) = -emat(isize, jsize)
         stiff(isize   , jsize+12) = -emat(isize, jsize)
         stiff(isize+12, jsize+12) =  emat(isize, jsize)
         
        END DO
        
       END DO
       
       DO ksize = 1, nn
       
        DO jsize = 1, 3
         
         DO isize = 1, nn*3
          
          !qf(isize) = qf(isize) + stiff(isize, (ksize-1)*3 + jsize ) * &
          !                ! ダミー剛性の寄与打消し分
          !              ( - u(jsize, ksize) - du(jsize, ksize) + &
          !                ! 正しい有効荷重
          !                 (param_b3 * du(jsize, ksize) - param_b1 * acc(jsize, ksize) - param_b2 * vel(jsize, ksize)) / coeff_damp)
          !                 !(-ddux(jsize, ksize)) / coeff_damp)
          qf(isize) = qf(isize) + stiff(isize, (ksize-1)*3 + jsize ) * &
                        !((param_b3 * du(jsize, ksize) - param_b1 * acc(jsize, ksize) - param_b2 * vel(jsize, ksize)) / coeff_damp)
                        !(-1.0d0 * ddux(jsize, ksize))
                        (-1.0d0 + rayleigh_beta * coeff_damp) * ddux(jsize, ksize)

         END DO
         
        END DO
       
       END DO
       
      END IF
      
      
!--------------------------------------------------------------------      
      
      RETURN
      
!####################################################################
      END SUBROUTINE Update_Joint_damp
!####################################################################
      
!####################################################################                  
       SUBROUTINE Set_stiff_Joint                           &            
                  (gausses, mtype, matl, k1, k2, k3)
!####################################################################                  
                  
       ! USE
       USE mMechGauss
       USE m_MatMatrix
       USE m_common_struct
       USE mMaterial
                  
 !--------------------------------------------------------------------
      
       ! TYPE(tGaussStatus), INTENT(IN) :: gausses(:)                ! status of Gaussian qudrature points
       INTEGER(KIND = kint), INTENT(IN)  :: mtype
       TYPE( tMaterial ), POINTER :: matl
       REAL(KIND = kREAL), INTENT(OUT) ::k1, k2, k3 

       INTEGER(KIND = kint), POINTER :: istatus
       REAL(KIND = kREAL), POINTER :: stiff_coeff
       
      !Added by Shan 13/07/2017
       TYPE(tGaussStatus) :: gausses(:)                ! status of Gaussian qudrature points
       REAL(KIND = kREAL) :: alpha,beta,stress_mean
       REAL(KIND = kREAL) :: stress_res(3)
       
!--------------------------------------------------------------------         

      !  pointer設定
      istatus          => gausses(1)%istatus(1)
      !istatus0         => gausses(1)%istatus(2)
      !flag_initset     => gausses(1)%istatus(3)
      
      !strain_tildash0  => gausses(1)%fstatus( 1   )
      !stress_tildash0  => gausses(1)%fstatus( 2   )
      stiff_coeff      => gausses(1)%fstatus( 3   )
      !strain_hat       => gausses(1)%fstatus( 4:11)
      !stress_hat       => gausses(1)%fstatus(12:19)
      !strain_pla0      => gausses(1)%fstatus(20:27)
      !strain_pla       => gausses(1)%fstatus(28:35)
      !strain_increment => gausses(1)%fstatus(36:43)
      !strain_init      => gausses(1)%fstatus(44:51)
      !stress_init      => gausses(1)%fstatus(52:59)
      
!--------------------------------------------------------------------       
       
       ! Goodman
       IF( mtype .EQ. GOODMAN ) THEN
       
        ! ! ! !Added by Shan 13/07/2017
        ! if (matl%variables(M_PLCONST1).lt.-1.0d0) then
            ! alpha   =  matl%variables(M_GOODMAN_K1)
            ! beta    =  matl%variables(M_GOODMAN_K2)
        ! if (gausses(1)%istatus(3)==1) then
            ! stress_res(1:3) = gausses(1)%stress(1:3)
            ! stress_mean = (stress_res(1)+stress_res(2)+stress_res(3))/3
            ! gausses(1)%fstatus(60) = stress_mean
            ! gausses(1)%istatus(3)  = 2
        ! end if
        
        ! call Set_stiff_Joint_restart(alpha,beta,gausses(1)%istatus(1),&   
                ! gausses(1)%fstatus(60),stiff_coeff,k1, k2, k3)
        
        ! else
            IF( istatus .EQ. 0 ) THEN
            
                k1 = matl%variables(M_GOODMAN_K1)
                k2 = matl%variables(M_GOODMAN_K2)
                k3 = matl%variables(M_GOODMAN_K3)
            
            ELSE IF( istatus .EQ. 1 ) THEN
            
                k1 = matl%variables(M_GOODMAN_K1)*0.0D+0
                k2 = matl%variables(M_GOODMAN_K2)*0.0D+0
                k3 = matl%variables(M_GOODMAN_K3)*0.0D+0
            
            ELSE IF( istatus .EQ. 2 ) THEN
            
                k1 = matl%variables(M_GOODMAN_K1)*stiff_coeff
                k2 = matl%variables(M_GOODMAN_K2)*stiff_coeff
                k3 = matl%variables(M_GOODMAN_K3)
            
            ELSE IF( istatus .EQ. 3 ) THEN
            
                k1 = matl%variables(M_GOODMAN_K1)*0.0D0
                k2 = matl%variables(M_GOODMAN_K2)*0.0D0
                k3 = matl%variables(M_GOODMAN_K3)
            
            END IF
        ! end if
       
       ELSE IF( mtype .EQ. SLIP_WEAKENING ) THEN
        
         IF( istatus .EQ. 0 ) THEN
         
           k1 = matl%variables(M_GOODMAN_K1)
           k2 = matl%variables(M_GOODMAN_K2)
           k3 = matl%variables(M_GOODMAN_K3)
         
         ELSE IF( istatus .EQ. 1 ) THEN
         
           k1 = matl%variables(M_GOODMAN_K1)*0.05D0
           k2 = matl%variables(M_GOODMAN_K2)*0.05D0
           k3 = matl%variables(M_GOODMAN_K3)
         
         ELSE IF( istatus .EQ. 2 ) THEN
         
           k1 = matl%variables(M_GOODMAN_K1)*0.05D0
           k2 = matl%variables(M_GOODMAN_K2)*0.05D0
           k3 = matl%variables(M_GOODMAN_K3)
         
         END IF   
        
       END IF

!####################################################################     
       END SUBROUTINE Set_stiff_Joint
!####################################################################

!####################################################################                  
       SUBROUTINE Set_stiff_Joint_restart                         &            
                  (alpha,beta,istatus,stress_mean,stiff_coeff,k1, k2, k3)
!####################################################################                  
    !Added by Shan 13/07/2017
    real(kind = kreal), intent(in)  :: alpha,beta,stiff_coeff
    real(kind = kreal), intent(out) :: k1,k2,k3
    integer(kind = kint), intent(in) :: istatus
    real(kind = kreal)  :: stress_mean
    
    k1 = alpha * (stress_mean**beta)
    k2 = k1
    k3 = k1
    
    if( istatus .eq. 1 ) then
         
         k1 = k1*0.001d+0
         k2 = k2*0.001d+0
         k3 = k3*0.001d+0
         
    else if( istatus .eq. 2 ) then
         
         k1 = k1*stiff_coeff
         k2 = k2*stiff_coeff
         k3 = k3
         
    else if( istatus .eq. 3 ) then
         
         k1 = 0.0d0
         k2 = 0.0d0
         k3 = k3
         
    end if

!####################################################################     
       END SUBROUTINE Set_stiff_Joint_restart
!####################################################################


!####################################################################                  
      SUBROUTINE Set_Bmat_Joint                           &            
                 (nn, t, ecoord, lmat)
!####################################################################                  
                  
      USE mMechGauss
      USE m_MatMatrix
      USE m_common_struct
      USE mMaterial                  
                  
 !--------------------------------------------------------------------
      
      INTEGER(KIND = kint), INTENT(IN) :: nn                      ! number of elemental nodes
      REAL(KIND = kREAL), INTENT(IN)   :: ecoord(3, nn)             ! coordinates of elemental nodes
      REAL(KIND = kREAL), INTENT(IN)   :: t(3,3)
      REAL(KIND = kREAL), INTENT(OUT)  :: lmat(8, 12)
      
!--------------------------------------------------------------------                 
                  
      REAL(KIND = kREAL) :: lmat_inv(12, 8)
      REAL(KIND = kREAL) :: ltxxl(8, 8)
      REAL(KIND = kREAL) :: ltxxl_i(8, 8)
      REAL(KIND = kREAL) :: coeff
      INTEGER(KIND = kint) :: i, j
      REAL(KIND = kREAL) :: a1(3), a2(3), a3(3), a4(3)
      REAL(KIND = kREAL) :: aw(3), aw1(3), aw2(3)

!--------------------------------------------------------------------
      
      DO i=1, 3
        a1(i)=(ecoord(i, 1)+ecoord(i, 5))*0.5d0
        a2(i)=(ecoord(i, 2)+ecoord(i, 6))*0.5d0
        a3(i)=(ecoord(i, 3)+ecoord(i, 7))*0.5d0
        a4(i)=(ecoord(i, 4)+ecoord(i, 8))*0.5d0
      END DO
      DO i=1, 3
        aw(i)=(a1(i)+a2(i)+a3(i)+a4(i))/4.0d0
        aw1(i)=(a1(i)+a2(i)+a4(i))/3.0d0
        aw2(i)=(a2(i)+a3(i)+a4(i))/3.0d0
      END DO
      DO i=1, 3
        a1(i)=a1(i)-aw(i)
        a2(i)=a2(i)-aw(i)
        a3(i)=a3(i)-aw(i)
        a4(i)=a4(i)-aw(i)
      ENDDO
      
      lmat_inv=0.0d0
      Do i=1,3
        lmat_inv(i  , 1) = t(1, i)
        lmat_inv(i+3, 1) = t(1, i)
        lmat_inv(i+6, 1) = t(1, i)
        lmat_inv(i+9, 1) = t(1, i)
        lmat_inv(i  , 2) = t(2, i)
        lmat_inv(i+3, 2) = t(2, i)
        lmat_inv(i+6, 2) = t(2, i)
        lmat_inv(i+9, 2) = t(2, i)
        lmat_inv(i  , 3) = t(3, i)
        lmat_inv(i+3, 3) = t(3, i)
        lmat_inv(i+6, 3) = t(3, i)
        lmat_inv(i+9, 3) = t(3, i) 
        
        !アワーグラスモード    
        lmat_inv(i  , 7) =  t(1, i)
        lmat_inv(i+3, 7) = -t(1, i)
        lmat_inv(i+6, 7) = -t(1, i)
        lmat_inv(i+9, 7) =  t(1, i)
         
        lmat_inv(i  , 8) =  t(2, i)
        lmat_inv(i+3, 8) =  t(2, i)
        lmat_inv(i+6, 8) = -t(2, i)
        lmat_inv(i+9, 8) = -t(2, i)
       enddo

       lmat_inv(1, 4) = a1(2)*t(1, 3) - a1(3)*t(1, 2)
       lmat_inv(2, 4) = a1(3)*t(1, 1) - a1(1)*t(1, 3)
       lmat_inv(3, 4) = a1(1)*t(1, 2) - a1(2)*t(1, 1)
       lmat_inv(4, 4) = a2(2)*t(1, 3) - a2(3)*t(1, 2)
       lmat_inv(5, 4) = a2(3)*t(1, 1) - a2(1)*t(1, 3)
       lmat_inv(6, 4) = a2(1)*t(1, 2) - a2(2)*t(1, 1)
       lmat_inv(7, 4) = a3(2)*t(1, 3) - a3(3)*t(1, 2)
       lmat_inv(8, 4) = a3(3)*t(1, 1) - a3(1)*t(1, 3)
       lmat_inv(9, 4) = a3(1)*t(1, 2) - a3(2)*t(1, 1)
       lmat_inv(10, 4)= a4(2)*t(1, 3) - a4(3)*t(1, 2)
       lmat_inv(11, 4)= a4(3)*t(1, 1) - a4(1)*t(1, 3)
       lmat_inv(12, 4)= a4(1)*t(1, 2) - a4(2)*t(1, 1)
       
       lmat_inv(1, 5) = a1(2)*t(2, 3) - a1(3)*t(2, 2)
       lmat_inv(2, 5) = a1(3)*t(2, 1) - a1(1)*t(2, 3)
       lmat_inv(3, 5) = a1(1)*t(2, 2) - a1(2)*t(2, 1)
       lmat_inv(4, 5) = a2(2)*t(2, 3) - a2(3)*t(2, 2)
       lmat_inv(5, 5) = a2(3)*t(2, 1) - a2(1)*t(2, 3)
       lmat_inv(6, 5) = a2(1)*t(2, 2) - a2(2)*t(2, 1)
       lmat_inv(7, 5) = a3(2)*t(2, 3) - a3(3)*t(2, 2)
       lmat_inv(8, 5) = a3(3)*t(2, 1) - a3(1)*t(2, 3)
       lmat_inv(9, 5) = a3(1)*t(2, 2) - a3(2)*t(2, 1)
       lmat_inv(10, 5)= a4(2)*t(2, 3) - a4(3)*t(2, 2)
       lmat_inv(11, 5)= a4(3)*t(2, 1) - a4(1)*t(2, 3)
       lmat_inv(12, 5)= a4(1)*t(2, 2) - a4(2)*t(2, 1)
       
       lmat_inv(1, 6) = a1(2)*t(3, 3) - a1(3)*t(3, 2)
       lmat_inv(2, 6) = a1(3)*t(3, 1) - a1(1)*t(3, 3)
       lmat_inv(3, 6) = a1(1)*t(3, 2) - a1(2)*t(3, 1)
       lmat_inv(4, 6) = a2(2)*t(3, 3) - a2(3)*t(3, 2)
       lmat_inv(5, 6) = a2(3)*t(3, 1) - a2(1)*t(3, 3)
       lmat_inv(6, 6) = a2(1)*t(3, 2) - a2(2)*t(3, 1)
       lmat_inv(7, 6) = a3(2)*t(3, 3) - a3(3)*t(3, 2)
       lmat_inv(8, 6) = a3(3)*t(3, 1) - a3(1)*t(3, 3)
       lmat_inv(9, 6) = a3(1)*t(3, 2) - a3(2)*t(3, 1)
       lmat_inv(10, 6)= a4(2)*t(3, 3) - a4(3)*t(3, 2)
       lmat_inv(11, 6)= a4(3)*t(3, 1) - a4(1)*t(3, 3)
       lmat_inv(12, 6)= a4(1)*t(3, 2) - a4(2)*t(3, 1)
       
       ltxxl = MATMUL( TRANSPOSE( lmat_inv ), lmat_inv)
       
       ltxxl_i = 0.0D0
       ltxxl_i(1, 1) = 1.0D0 / ltxxl(1, 1)
       ltxxl_i(2, 2) = 1.0D0 / ltxxl(2, 2)
       ltxxl_i(3, 3) = 1.0D0 / ltxxl(3, 3)
       ltxxl_i(6, 6) = 1.0D0 / ltxxl(6, 6)
       
       ltxxl_i(7, 7) = 1.0D0 / ltxxl(7, 7)
       ltxxl_i(8, 8) = 1.0D0 / ltxxl(8, 8)
       
       ! 回転部分のdet
       coeff = ltxxl(4, 4) * ltxxl(5, 5) - ltxxl(4, 5) * ltxxl(5, 4)
       ! 2x2の逆行列
       ltxxl_i(4, 4) =  ltxxl(5, 5) / coeff
       ltxxl_i(5, 5) = ltxxl(4, 4) / coeff
       ltxxl_i(4, 5) = -1.0D0 * ltxxl(5, 4) / coeff
       ltxxl_i(5, 4) = ltxxl_i(4, 5)
       
       lmat = MATMUL(ltxxl_i, TRANSPOSE( lmat_inv ))

!####################################################################           
     end SUBROUTINE Set_Bmat_Joint
!####################################################################      

!####################################################################                  
      SUBROUTINE Set_Bmat_Joint_dp                         &            
                 (nn, t, ecoord, lmat)
!####################################################################                  
                  
      USE mMechGauss
      USE m_MatMatrix
      USE m_common_struct
      USE mMaterial                  
                  
 !--------------------------------------------------------------------
      
      INTEGER(KIND = kint), INTENT(IN) :: nn                      ! number of elemental nodes
      REAL(KIND = kREAL), INTENT(IN)   :: ecoord(3, nn)             ! coordinates of elemental nodes
      REAL(KIND = kREAL), INTENT(IN)   :: t(3,3)
      REAL(KIND = kREAL), INTENT(OUT)  :: lmat(12, 12)
      
!--------------------------------------------------------------------                 
                  
      REAL(KIND = kREAL) :: lmat_inv_dp(12, 12)
      INTEGER(KIND = kint) :: i, j

!--------------------------------------------------------------------
      
      lmat_inv_dp=0.0d0
      Do i=1,3
        lmat_inv_dp(i  ,  1) = t(1, i)
        lmat_inv_dp(i+3,  2) = t(1, i)
        lmat_inv_dp(i+6,  3) = t(1, i)
        lmat_inv_dp(i+9,  4) = t(1, i)
        lmat_inv_dp(i  ,  5) = t(2, i)
        lmat_inv_dp(i+3,  6) = t(2, i)
        lmat_inv_dp(i+6,  7) = t(2, i)
        lmat_inv_dp(i+9,  8) = t(2, i)
        lmat_inv_dp(i  ,  9) = t(3, i)
        lmat_inv_dp(i+3, 10) = t(3, i)
        lmat_inv_dp(i+6, 11) = t(3, i)
        lmat_inv_dp(i+9, 12) = t(3, i)
      enddo
       
      lmat = TRANSPOSE( lmat_inv_dp )

!####################################################################           
     end SUBROUTINE Set_Bmat_Joint_dp
!####################################################################      
                 
                 
!####################################################################                  
       SUBROUTINE Set_Dmat_Joint                           &            
                  (nn, ecoord, t, area, k1, k2, k3, e1_hat0, e2_hat0, &
                   e11_hat, e12_hat, e11_hat0, e12_hat0, IIw)
!####################################################################                  
                  
       ! USE       
                  
 !--------------------------------------------------------------------
      
       INTEGER(KIND = kint), INTENT(IN) :: nn                      ! number of elemental nodes
       REAL(KIND = kREAL), INTENT(IN)   :: ecoord(3, nn)           ! coordinates of elemental nodes
       REAL(KIND = kREAL), INTENT(IN)  :: t(3, 3)
       REAL(KIND = kREAL), INTENT(IN)  :: area
       REAL(KIND = kREAL), INTENT(IN)  :: k1, k2, k3
       REAL(KIND = kREAL), INTENT(IN)  :: e1_hat0(3), e2_hat0(3)
       REAL(KIND = kREAL), INTENT(IN)  :: e11_hat(3), e12_hat(3)
       REAL(KIND = kREAL), INTENT(IN)  :: e11_hat0(3), e12_hat0(3)
       REAL(KIND = kREAL), INTENT(OUT) :: IIw(8, 8)
       
!--------------------------------------------------------------------
       
       REAL(KIND = kREAL) :: aa1, aa2, bb1, bb2
       REAL(KIND = kREAL) :: cc1, cc2, dd1, dd2
       REAL(KIND = kREAL) :: a1(3), a2(3), a3(3), a4(3)
       REAL(KIND = kREAL) :: aw(3), aw1(3), aw2(3)
       INTEGER(KIND = kint) :: i, j
       
!--------------------------------------------------------------------     
       
       DO i=1, 3
         a1(i)=(ecoord(i, 1)+ecoord(i, 5))*0.5d0
         a2(i)=(ecoord(i, 2)+ecoord(i, 6))*0.5d0
         a3(i)=(ecoord(i, 3)+ecoord(i, 7))*0.5d0
         a4(i)=(ecoord(i, 4)+ecoord(i, 8))*0.5d0
       END DO
       DO i=1, 3
         aw(i)=(a1(i)+a2(i)+a3(i)+a4(i))/4.0d0
         aw1(i)=(a1(i)+a2(i)+a4(i))/3.0d0
         aw2(i)=(a2(i)+a3(i)+a4(i))/3.0d0
       END DO
       DO i=1, 3
         a1(i)=a1(i)-aw(i)
         a2(i)=a2(i)-aw(i)
         a3(i)=a3(i)-aw(i)
         a4(i)=a4(i)-aw(i)
       ENDDO
       
       !要素に沿った座標系で三角形を表現する
       aa1=t(1,1)*e1_hat0(1)+t(1,2)*e1_hat0(2)+t(1,3)*e1_hat0(3)
       aa2=t(2,1)*e1_hat0(1)+t(2,2)*e1_hat0(2)+t(2,3)*e1_hat0(3) !この成分は必ずゼロになるはず
       bb1=t(1,1)*e2_hat0(1)+t(1,2)*e2_hat0(2)+t(1,3)*e2_hat0(3)
       bb2=t(2,1)*e2_hat0(1)+t(2,2)*e2_hat0(2)+t(2,3)*e2_hat0(3)
       
       cc1=t(1,1)*(a1(1)+aw(1)-aw1(1))+t(1,2)*(a1(2)+aw(2)-aw1(2))+t(1,3)*(a1(3)+aw(3)-aw1(3))
       cc2=t(2,1)*(a1(1)+aw(1)-aw1(1))+t(2,2)*(a1(2)+aw(2)-aw1(2))+t(2,3)*(a1(3)+aw(3)-aw1(3))
       dd1=t(1,1)*(aw1(1)-aw(1))      +t(1,2)*(aw1(2)-aw(2))      +t(1,3)*(aw1(3)-aw(3))
       dd2=t(2,1)*(aw1(1)-aw(1))      +t(2,2)*(aw1(2)-aw(2))      +t(2,3)*(aw1(3)-aw(3))
              
       ! 剛性マトリックス(D)の定義
       ! 応力→力の変換を含まない
       IIw(1,1)=k1
       IIw(2,2)=k2
       IIw(3,3)=k3
       
       IIw(7,7)=0.5D0*k1
       IIw(8,8)=0.5D0*k2
       
       ! 回転剛性テンソルの定義
       IIw(4,4)=k3*((      aa1*aa1 + aa1*bb1 + bb1*bb1                ) / 12.0D0 - (cc1*cc1-dd1*dd1)/2.0D0)
       IIw(4,5)=k3*((2.0D0*aa1*aa2 + aa2*bb1 + aa1*bb2 + 2.0D0*bb1*bb2) / 24.0D0 - (cc1*cc2-dd1*dd2)/2.0D0)
       IIw(5,5)=k3*((      aa2*aa2 + aa2*bb2 + bb2*bb2                ) / 12.0D0 - (cc2*cc2-dd2*dd2)/2.0D0)
       IIw(6,6)=k2*((      aa1*aa1 + aa1*bb1 + bb1*bb1                ) / 12.0D0 - (cc1*cc1-dd1*dd1)/2.0D0) & 
               +k1*((      aa2*aa2 + aa2*bb2 + bb2*bb2                ) / 12.0D0 - (cc2*cc2-dd2*dd2)/2.0D0)
               
       aa1=e11_hat(1)*e11_hat0(1)+e11_hat(2)*e11_hat0(2)+e11_hat(3)*e11_hat0(3)
       aa2=e12_hat(1)*e11_hat0(1)+e12_hat(2)*e11_hat0(2)+e12_hat(3)*e11_hat0(3) !この成分は必ずゼロになるはず
       bb1=e11_hat(1)*e12_hat0(1)+e11_hat(2)*e12_hat0(2)+e11_hat(3)*e12_hat0(3)
       bb2=e12_hat(1)*e12_hat0(1)+e12_hat(2)*e12_hat0(2)+e12_hat(3)*e12_hat0(3) 
       
       cc1=t(1,1)*(a3(1)+aw(1)-aw2(1))+t(1,2)*(a3(2)+aw(2)-aw2(2))+t(1,3)*(a3(3)+aw(3)-aw2(3))
       cc2=t(2,1)*(a3(1)+aw(1)-aw2(1))+t(2,2)*(a3(2)+aw(2)-aw2(2))+t(2,3)*(a3(3)+aw(3)-aw2(3))
       dd1=t(1,1)*(aw2(1)-aw(1))+t(1,2)*(aw2(2)-aw(2))+t(1,3)*(aw2(3)-aw(3))
       dd2=t(2,1)*(aw2(1)-aw(1))+t(2,2)*(aw2(2)-aw(2))+t(2,3)*(aw2(3)-aw(3))
       
       ! 回転剛性テンソルの定義
       IIw(4,4)=IIw(4,4)+k3*((      aa1*aa1 + aa1*bb1 + bb1*bb1                ) / 12.0D0 - (cc1*cc1-dd1*dd1)/2.0D0)
       IIw(4,5)=IIw(4,5)+k3*((2.0D0*aa1*aa2 + aa2*bb1 + aa1*bb2 + 2.0D0*bb1*bb2) / 24.0D0 - (cc1*cc2-dd1*dd2)/2.0D0)
       IIw(5,5)=IIw(5,5)+k3*((      aa2*aa2 + aa2*bb2 + bb2*bb2                ) / 12.0D0 - (cc2*cc2-dd2*dd2)/2.0D0)
       IIw(6,6)=IIw(6,6)+k2*((      aa1*aa1 + aa1*bb1 + bb1*bb1                ) / 12.0D0 - (cc1*cc1-dd1*dd1)/2.0D0) & 
                        +k1*((      aa2*aa2 + aa2*bb2 + bb2*bb2                ) / 12.0D0 - (cc2*cc2-dd2*dd2)/2.0D0) 
       
       IIw(1,1)=IIw(1,1)*area
       IIw(2,2)=IIw(2,2)*area
       IIw(3,3)=IIw(3,3)*area
       IIw(4,4)=IIw(4,4)*area
       IIw(4,5)=IIw(4,5)*area
       IIw(5,4)=IIw(4,5)
       IIw(5,5)=IIw(5,5)*area
       IIw(6,6)=IIw(6,6)*area
       
       IIw(7,7)=IIw(7,7)*area
       IIw(8,8)=IIw(8,8)*area


!####################################################################     
       END SUBROUTINE set_Dmat_Joint
!####################################################################

!####################################################################                  
       SUBROUTINE Set_Dmat_Joint_dp                        &            
                  (area, rho, Vs, Vp, IIw_dp)
!####################################################################
                  
       ! USE       
                  
 !--------------------------------------------------------------------
      
       REAL(KIND = kREAL), INTENT(IN)  :: area
       REAL(KIND = kREAL), INTENT(IN)  :: rho, Vs, Vp
       REAL(KIND = kREAL), INTENT(OUT) :: IIw_dp(12, 12)
       
!--------------------------------------------------------------------
       
       ! 剛性マトリックス(D)の定義
       IIw_dp(1,1)=rho * Vs * area * (1.0d0 / 4.0D0)
       IIw_dp(2,2)=rho * Vs * area * (1.0d0 / 4.0D0)
       IIw_dp(3,3)=rho * Vs * area * (1.0d0 / 4.0D0)
       IIw_dp(4,4)=rho * Vs * area * (1.0d0 / 4.0D0)
       IIw_dp(5,5)=rho * Vs * area * (1.0d0 / 4.0D0)
       IIw_dp(6,6)=rho * Vs * area * (1.0d0 / 4.0D0)
       IIw_dp(7,7)=rho * Vs * area * (1.0d0 / 4.0D0)
       IIw_dp(8,8)=rho * Vs * area * (1.0d0 / 4.0D0)
       IIw_dp(9,9)=rho * Vp * area * (1.0d0 / 4.0D0)
       IIw_dp(10,10)=rho * Vp * area * (1.0d0 / 4.0D0)
       IIw_dp(11,11)=rho * Vp * area * (1.0d0 / 4.0D0)
       IIw_dp(12,12)=rho * Vp * area * (1.0d0 / 4.0D0)
       

!####################################################################     
       END SUBROUTINE set_Dmat_Joint_dp
!####################################################################
                  
!####################################################################                  
       SUBROUTINE Set_LocalCood_Joint                           &            
                  (nn, ecoord, t, area, e1_hat0, e2_hat0, &
                   e11_hat, e12_hat, e11_hat0, e12_hat0)
!####################################################################                  
                  
       ! USE       
                  
 !--------------------------------------------------------------------
       
       INTEGER(KIND = kint), INTENT(IN) :: nn                      ! number of elemental nodes
       REAL(KIND = kREAL), INTENT(IN)   :: ecoord(3, nn)           ! coordinates of elemental nodes
       REAL(KIND = kREAL), INTENT(OUT)  :: t(3, 3)
       REAL(KIND = kREAL), INTENT(OUT)  :: area
       REAL(KIND = kREAL), INTENT(OUT)  :: e1_hat0(3), e2_hat0(3)
       REAL(KIND = kREAL), INTENT(OUT)  :: e11_hat(3), e12_hat(3)
       REAL(KIND = kREAL), INTENT(OUT)  :: e11_hat0(3), e12_hat0(3)
       
!--------------------------------------------------------------------         
       
       REAL(KIND = kREAL) :: e1_hat(3), e2_hat(3), e3_hat(3)
       REAL(KIND = kREAL) :: e13_hat(3)
       REAL(KIND = kREAL) :: e1_hat_abs, e2_hat_abs, e3_hat_abs
       REAL(KIND = kREAL) :: e11_hat_abs, e12_hat_abs, e13_hat_abs
       REAL(KIND = kREAL) :: x21(3), x41(3), x23(3), x43(3)
       REAL(KIND = kREAL) :: x65(3), x85(3), x67(3), x87(3)
       INTEGER(KIND = kint) :: i, j
       
!--------------------------------------------------------------------
       
       DO i = 1, 3
        
        x21(i) = ecoord(i, 2)-ecoord(i, 1)
        x41(i) = ecoord(i, 4)-ecoord(i, 1)
        
        x65(i) = ecoord(i, 6)-ecoord(i, 5)
        x85(i) = ecoord(i, 8)-ecoord(i, 5)
        
        ! 面積の計算用
        x23(i) = ecoord(i, 2)-ecoord(i, 3)
        x43(i) = ecoord(i, 4)-ecoord(i, 3)
        
        x67(i) = ecoord(i, 6)-ecoord(i, 7)
        x87(i) = ecoord(i, 8)-ecoord(i, 7)
        
       END DO
       
       !--------------------------------------------------------
       
       ! 要素に沿った座標系を定義
       DO i = 1, 3
        e1_hat0(i) = 0.50D0*( x21(i)+x65(i) )
        e2_hat0(i) = 0.50D0*( x41(i)+x85(i) )
        
        e11_hat0(i) = 0.50D0*( x23(i)+x67(i) )
        e12_hat0(i) = 0.50D0*( x43(i)+x87(i) )
       end do
       
       e1_hat(1) = e1_hat0(1)
       e1_hat(2) = e1_hat0(2)
       e1_hat(3) = e1_hat0(3)
       
       e3_hat(1) = e1_hat0(2)*e2_hat0(3)-e1_hat0(3)*e2_hat0(2)
       e3_hat(2) = e1_hat0(3)*e2_hat0(1)-e1_hat0(1)*e2_hat0(3)
       e3_hat(3) = e1_hat0(1)*e2_hat0(2)-e1_hat0(2)*e2_hat0(1)
       
       e2_hat(1) = e3_hat(2)*e1_hat0(3)-e3_hat(3)*e1_hat0(2)
       e2_hat(2) = e3_hat(3)*e1_hat0(1)-e3_hat(1)*e1_hat0(3)
       e2_hat(3) = e3_hat(1)*e1_hat0(2)-e3_hat(2)*e1_hat0(1)
       
       e1_hat_abs = DSQRT( e1_hat(1)*e1_hat(1)   &
                          +e1_hat(2)*e1_hat(2)   &
                          +e1_hat(3)*e1_hat(3) ) 
       e2_hat_abs = DSQRT( e2_hat(1)*e2_hat(1)   &
                          +e2_hat(2)*e2_hat(2)   &
                          +e2_hat(3)*e2_hat(3) )
       e3_hat_abs = DSQRT( e3_hat(1)*e3_hat(1)   &
                          +e3_hat(2)*e3_hat(2)   &
                          +e3_hat(3)*e3_hat(3) )
                          
       ! N2, N3, N4	が形成する三角形
       e11_hat(1) = e11_hat0(1)
       e11_hat(2) = e11_hat0(2)
       e11_hat(3) = e11_hat0(3)
       
       e13_hat(1) = e11_hat0(2)*e12_hat0(3)-e11_hat0(3)*e12_hat0(2)
       e13_hat(2) = e11_hat0(3)*e12_hat0(1)-e11_hat0(1)*e12_hat0(3)
       e13_hat(3) = e11_hat0(1)*e12_hat0(2)-e11_hat0(2)*e12_hat0(1)
       
       e12_hat(1) = e13_hat(2)*e11_hat0(3)-e13_hat(3)*e11_hat0(2)
       e12_hat(2) = e13_hat(3)*e11_hat0(1)-e13_hat(1)*e11_hat0(3)
       e12_hat(3) = e13_hat(1)*e11_hat0(2)-e13_hat(2)*e11_hat0(1)
       
       e11_hat_abs = DSQRT( e11_hat(1)*e11_hat(1)   &
                           +e11_hat(2)*e11_hat(2)   &
                           +e11_hat(3)*e11_hat(3) ) 
       e12_hat_abs = DSQRT( e12_hat(1)*e12_hat(1)   &
                           +e12_hat(2)*e12_hat(2)   &
                           +e12_hat(3)*e12_hat(3) )
       e13_hat_abs = DSQRT( e13_hat(1)*e13_hat(1)   &
                           +e13_hat(2)*e13_hat(2)   &
                           +e13_hat(3)*e13_hat(3) )
       
       ! 要素の面積は要素の二辺のベクトル×0.5 
       area=(e3_hat_abs+e13_hat_abs)*0.5D0
       
       ! 正規化
       DO i = 1, 3
        
        e1_hat(i) = e1_hat(i)/e1_hat_abs
        e2_hat(i) = e2_hat(i)/e2_hat_abs
        e3_hat(i) = e3_hat(i)/e3_hat_abs
        
        e11_hat(i) = e11_hat(i)/e11_hat_abs
        e12_hat(i) = e12_hat(i)/e12_hat_abs
        e13_hat(i) = e13_hat(i)/e13_hat_abs
        
       END DO

       DO i = 1, 3
        
        t(1, i) = e1_hat(i)
        t(2, i) = e2_hat(i)
        t(3, i) = e3_hat(i)
        
       END DO

!####################################################################     
       END SUBROUTINE Set_LocalCood_Joint
!####################################################################

!####################################################################                  
       SUBROUTINE Set_param_Joint_dp                           &            
                  (matl, Vp, Vs, rho, rayleigh_alpha, rayleigh_beta, &
                   delta_t, newmark_beta, newmark_gamma, coeff_damp, fstrDYNAMIC)
!####################################################################                  
                  
       ! USE
       use m_fstr
       USE mMechGauss
       USE m_MatMatrix
       USE m_common_struct
       USE mMaterial
                  
 !--------------------------------------------------------------------
       
       TYPE( tMaterial ), POINTER :: matl
       REAL(KIND = kREAL), INTENT(OUT) :: Vs, Vp, rho
       REAL(KIND = kREAL), INTENT(OUT) :: rayleigh_alpha, rayleigh_beta, delta_t
       REAL(KIND = kREAL), INTENT(OUT) :: newmark_beta, newmark_gamma
       REAL(KIND = kREAL), INTENT(OUT) :: coeff_damp
       type ( fstr_dynamic ), optional, INTENT(IN) :: fstrDYNAMIC
       
!--------------------------------------------------------------------
       
       REAL(KIND = kREAL) :: param_b3
       !REAL(KIND = kREAL) :: param_a1, param_b1, param_a2, param_b2
       
!--------------------------------------------------------------------     
       
       Vp = matl%variables(M_DASHPOD_VP)
       Vs = matl%variables(M_DASHPOD_VS)
       rho = matl%variables(M_GOODMAN_C)
         
       rayleigh_alpha   = fstrDYNAMIC%ray_m
       rayleigh_beta    = fstrDYNAMIC%ray_k
       delta_t = fstrDYNAMIC%t_delta
         
       newmark_beta = fstrDYNAMIC%beta
       newmark_gamma =  fstrDYNAMIC%ganma
       
       param_b3 = newmark_gamma / (newmark_beta * delta_t)        
       coeff_damp = param_b3 / (1.0d0 + rayleigh_beta * param_b3)
         
!####################################################################     
       END SUBROUTINE Set_param_Joint_dp
!####################################################################     

       
!####################################################################                  
       SUBROUTINE Set_initStress_Joint                        &            
                  (nn, matl, ecoord,  t, &
                   sigma_init, stress_init, strain_init, flag_initset)
!####################################################################                  
                  
       ! USE
       USE mMechGauss
       USE m_MatMatrix
       USE m_common_struct
       USE mMaterial    
                  
 !--------------------------------------------------------------------
      
       INTEGER(KIND = kint), INTENT(IN) :: nn   
       TYPE( tMaterial ), POINTER       :: matl
       REAL(KIND = kREAL), INTENT(IN)   :: ecoord(3, nn)             ! coordinates of elemental nodes
       REAL(KIND = kREAL), INTENT(IN)   :: t(3,3)
       REAL(KIND = kREAL), INTENT(OUT)  :: sigma_init(3,3)
       REAL(KIND = kREAL), INTENT(OUT)  :: stress_init(6)
       REAL(KIND = kREAL), INTENT(OUT)  :: strain_init(6)
       INTEGER(KIND = kint)             :: flag_initset
       
!--------------------------------------------------------------------

       REAL(KIND = kREAL) :: ave_cood_z
       REAL(KIND = kREAL) :: coeff_a, coeff_b, coeff_c
       REAL(KIND = kREAL) :: k1_init, k2_init, k3_init
       
!--------------------------------------------------------------------
       
         ! 初期応力定義の為の係数
         coeff_a = 0.0d0
         !coeff_a = 20.0d6
         coeff_b = 0.0d0
         !coeff_b = -1.0d0
         coeff_c = 0.5d0
         
         !--------------------------------------------------------
        
         ! 初期剛性を取得する（おそらく不要）
         k1_init = matl%variables(M_GOODMAN_K1)
         k2_init = matl%variables(M_GOODMAN_K2)
         k3_init = matl%variables(M_GOODMAN_K3)
        
         !--------------------------------------------------------
        
         ! 鉛直座標の平均値を取得
         ave_cood_z = 0.125d0 * (  ecoord(3, 1) + ecoord(3, 2) + ecoord(3, 3) + ecoord(3, 4) &
                                 + ecoord(3, 5) + ecoord(3, 6) + ecoord(3, 7) + ecoord(3, 8) )
       
         !--------------------------------------------------------
        
         sigma_init(:,:)=0.0d0
        
         !sigma_init(3,3) = coeff_a * ave_cood_z + coeff_b
         !sigma_init(1,1) = coeff_c * sigma_init(3,3)
         !sigma_init(2,2) = coeff_c * sigma_init(3,3)
        
         if(matl%variables(M_GOODMAN_C) .ge. 0.0d0) then
           
           stress_init(1)=dabs(matl%variables(M_DASHPOD_VP))
           stress_init(2)=dabs(matl%variables(M_PLCONST1))
           stress_init(3)=dabs(matl%variables(M_DASHPOD_VS))
           
         else
         
           coeff_a = matl%variables(M_DASHPOD_VP)
           coeff_b = matl%variables(M_PLCONST1)
           coeff_c = matl%variables(M_DASHPOD_VS)
            
           
           !! グローバル座標系での値
           !sigma_init(3,3) = coeff_a * ave_cood_z + coeff_b
           !sigma_init(1,1) = coeff_c * sigma_init(3,3)
           !sigma_init(2,2) = coeff_c * sigma_init(3,3)
           !
           !
           !! 局所座標系に変換
           !sigma_init(   1:3, :) = MATMUL( t, sigma_init( :, :))
           !sigma_init( :, :) = MATMUL( sigma_init( :, :), t_transpose )
           
           !stress_init(1) = sigma_init( 1, 3)
           !stress_init(2) = sigma_init( 2, 3)
           !stress_init(3) = sigma_init( 3, 3)
           
           stress_init(1) = coeff_a * ave_cood_z + matl%variables(M_GOODMAN_EPS)
           stress_init(2) = coeff_b * ave_cood_z
           stress_init(3) = coeff_c * ave_cood_z
           !write(101, "(100E12.3)") coeff_a, coeff_b, coeff_c, ave_cood_z
           !write(101, "(100E12.3)") stress_init(:)
           
         endif
         
         !sigma_init(3,3) = stress_init(3)
         !sigma_init(1,1) = stress_init(1)
         !sigma_init(2,2) = stress_init(2)
         !
         !! 局所座標系に変換
         !sigma_init(   1:3, :) = MATMUL( t, sigma_init( :, :))
         !sigma_init( :, :) = MATMUL( sigma_init( :, :), TRANSPOSE(t) )
         !
         !write(101, "(100E12.3)") sigma_init(:, :)
        
         !--------------------------------------------------------
        
         !stress_init(1)=sigma_init(1, 3)
         !stress_init(2)=sigma_init(2, 3)
         !stress_init(3)=sigma_init(3, 3)
         stress_init(4)= 0.0d0
         stress_init(5)= 0.0d0
         stress_init(6)= 0.0d0
        
         strain_init(1)=stress_init(1) / k1_init
         strain_init(2)=stress_init(2) / k2_init
         strain_init(3)=stress_init(3) / k3_init
         strain_init(4)=0.0d0
         strain_init(5)=0.0d0
         strain_init(6)=0.0d0
         
         ! すべり弱化モデルの場合は初期ひずみ（変位）は設定しない
         !if(mtype .EQ. SLIP_WEAKENING)then
         !    strain_init(:)=0.0d0
         !endif
         !->これを行うと増分ひずみと全応力が比例することになりよくない！
       
         ! 初期応力設定済みフラグ
         flag_initset = 1
        
         ! 一応初期化
         coeff_a=0.0d0
         coeff_b=0.0d0
         coeff_c=0.0d0
         
         !write(101, "(100E12.3)") stress_init(:)
       

!####################################################################     
       END SUBROUTINE Set_initStress_Joint
!####################################################################         


!####################################################################                  
    SUBROUTINE Set_initStress_Joint_restart(nn,k1,k2,k3,          &
                    t,stress_init,strain_init,flag_initset) 
!#################################################################### 
    !Added by Shan 14/07/2017
    
     !--------------------------------------------------------------------
      
    INTEGER(KIND = kint), INTENT(IN) :: nn
    REAL(KIND = kREAL), INTENT(IN)   :: k1,k2,k3
    REAL(KIND = kREAL), INTENT(IN)   :: t(3,3)
    REAL(KIND = kREAL)               :: sigma_init(3,3)
    REAL(KIND = kREAL)               :: stress_init(8),strain_init(8)
    INTEGER(KIND = kint)             :: flag_initset
       
    !--------------------------------------------------------------------
    sigma_init = 0.0d0
    ! sigma_init(1, 1) = stress_init(1)
    ! sigma_init(2, 2) = stress_init(2)
    ! sigma_init(3, 3) = stress_init(3)
    
    ! ! 局所座標系に変換
    ! sigma_init(:,:) = MATMUL( t, sigma_init(:,:))
    ! sigma_init(:,:) = MATMUL( sigma_init(:,:), TRANSPOSE(t) )
    
    ! stress_init(1)=sigma_init(1, 3)
    ! stress_init(2)=sigma_init(2, 3)
    ! stress_init(3)=sigma_init(3, 3)
    ! stress_init(4)= 0.0d0
    ! stress_init(5)= 0.0d0
    ! stress_init(6)= 0.0d0
    
    strain_init(1)=stress_init(1) / k1
    strain_init(2)=stress_init(2) / k2
    strain_init(3)=stress_init(3) / k3
    strain_init(4)=0.0d0
    strain_init(5)=0.0d0
    strain_init(6)=0.0d0

    flag_initset = 3
!####################################################################     
       END SUBROUTINE Set_initStress_Joint_restart
!####################################################################  

      END MODULE m_static_LIB_Joint
      MODULE m_static_LIB_Joint3n
      
      USE hecmw, ONLY : kint, kREAL
      USE elementInfo
      
!--------------------------------------------------------------------
      
      IMPLICIT NONE
      
!--------------------------------------------------------------------
      
      CONTAINS
      
      
!####################################################################
      SUBROUTINE STF_Joint3n                           &
                 (etype, nn, ecoord, gausses, stiff) 
!####################################################################
      
      USE mMechGauss
      USE m_MatMatrix
      USE m_common_struct
      USE mMaterial
      
!--------------------------------------------------------------------
      
      INTEGER(KIND = kint), INTENT(IN) :: etype                   ! element type
      INTEGER(KIND = kint), INTENT(IN) :: nn                      ! number of elemental nodes
      REAL(KIND = kREAL), INTENT(IN) :: ecoord(3, nn)             ! coordinates of elemental nodes
      TYPE(tGaussStatus), INTENT(IN) :: gausses(:)                ! status of Gaussian qudrature points
      REAL(KIND = kREAL), INTENT(OUT) :: stiff(:,:)               ! stiffness matrix
      
!--------------------------------------------------------------------
      
      TYPE( tMaterial ), POINTER :: matl
      
      INTEGER(KIND = kint), PARAMETER :: ndof = 3
      INTEGER(KIND = kint) :: i, j
      INTEGER(KIND = kint) :: isize, jsize !, ksize, lsize
      INTEGER(KIND = kint) :: mtype
      
      REAL(KIND = kREAL) :: k1, k2, k3
      REAL(KIND = kREAL) :: e1_hat0(3), e2_hat0(3)
      REAL(KIND = kREAL) :: t(3, 3)!, t_transpose(3, 3)
      REAL(KIND = kREAL) :: emat(9, 9)
      REAL(KIND = kREAL) :: lmat_inv_dp(9, 9)
      REAL(KIND = kREAL) :: coeff, area
      REAL(KIND = kREAL) :: lmat_dp(9, 9)
      REAL(KIND = kREAL) :: aa1, aa2, bb1, bb2
      REAL(KIND = kREAL) :: IIw_dp(9, 9)
      REAL(KIND = kREAL) :: Vs, Vp, rho
      REAL(KIND = kREAL) :: rayleigh_alpha, rayleigh_beta, delta_t
      REAL(KIND = kREAL) :: newmark_beta, newmark_gamma
      REAL(KIND = kREAL) :: param_a3, param_b3, coeff_damp
      
      REAL(KIND = kREAL) :: lmat_inv(9, 6)
      REAL(KIND = kREAL) :: lmat(6, 9)
      REAL(KIND = kREAL) :: IIw(6, 6)
!--------------------------------------------------------------------

      !  pointer設定
      !istatus          => gausses(1)%istatus(1)
      !istatus0         => gausses(1)%istatus(2)
      !flag_initset     => gausses(1)%istatus(3)
      !
      !strain_tildash0  => gausses(1)%fstatus( 1   )
      !stress_tildash0  => gausses(1)%fstatus( 2   )
      !stiff_coeff      => gausses(1)%fstatus( 3   )
      !strain_hat       => gausses(1)%fstatus( 4:11)
      !stress_hat       => gausses(1)%fstatus(12:19)
      !strain_pla0      => gausses(1)%fstatus(20:27)
      !strain_pla       => gausses(1)%fstatus(28:35)
      !strain_increment => gausses(1)%fstatus(36:43)
      !strain_init      => gausses(1)%fstatus(44:51)
      !stress_init      => gausses(1)%fstatus(52:59)
      
!--------------------------------------------------------------------
      
      stiff(:, :) = 0.0D0
      
!--------------------------------------------------------------------
      
      matl=>gausses(1)%pMaterial
      
      mtype = matl%mtype
      
!--------------------------------------------------------------------
       coeff_damp = 1.0d0
      
       !--------------------------------------------------------
       
       ! Goodman
       IF( mtype .EQ. GOODMAN ) THEN
        
         call Set_stiff_Joint3n(gausses, matl, k1, k2, k3)
       
         call Set_LocalCood_Joint3n(nn, ecoord, t, area, e1_hat0, e2_hat0)
       
         lmat_inv=0.0d0
       
         call Set_Bmat_Joint3n(nn, t, ecoord, lmat)

         IIw = 0.0D0
         call Set_Dmat_Joint3n(t, area, k1, k2, k3, e1_hat0, e2_hat0, IIw)
       
         emat = MATMUL( TRANSPOSE( lmat ) , MATMUL( IIw, lmat))
       
       else IF( ( mtype .EQ. DASHPOD_ELEMENT ) ) THEN 
       
       !--------------------------------------------------------
         
         call Set_param_Joint3n_dp(matl, Vp, Vs, rho, rayleigh_alpha, rayleigh_beta, &
                                   delta_t, newmark_beta, newmark_gamma, coeff_damp)
               
         call Set_LocalCood_Joint3n(nn, ecoord, t, area, e1_hat0, e2_hat0)
             
         lmat_dp=0.0d0
  
         call Set_Bmat_Joint3n_dp(nn, t, ecoord, lmat_dp)
       
         IIw_dp = 0.0D0

         Call Set_Dmat_Joint3n_dp(area, rho, Vs, Vp, IIw_dp)
       
         emat = MATMUL( TRANSPOSE( lmat_dp ) , MATMUL( IIW_dp, lmat_dp))
       
      
       !--------------------------------------------------------
     
       END IF
       
       DO jsize = 1, 9
          
         DO isize = 1, 9
         
         ! ここでnewmark beta法の部分でつじつまが合うような係数coeff_dampをかける
           stiff(isize  , jsize  ) =  emat(isize, jsize) * coeff_damp
           stiff(isize+9, jsize  ) = -emat(isize, jsize) * coeff_damp
           stiff(isize  , jsize+9) = -emat(isize, jsize) * coeff_damp
           stiff(isize+9, jsize+9) =  emat(isize, jsize) * coeff_damp
         
         END DO
      
       END DO
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE STF_Joint3n
!####################################################################

!####################################################################
      SUBROUTINE Update_Joint3n                                  &
                 (icel, etype, nn, ecoord, u, du, qf, gausses, &
                  !strain_hat, stress_hat, strain_pla, strain_pla0, &
                  !strain_increment, &
                  !strain_init, stress_init, &
                  !flag_initset, &
                  iter)                      
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
      INTEGER(KIND = kint) :: isize, jsize!, ksize, lsize
      INTEGER(KIND = kint) :: mtype
      
      REAL(KIND = kREAL) :: k1, k2, k3 !, kk1, kk2, kk3
      REAL(KIND = kREAL) :: phi, c
      REAL(KIND = kREAL) :: epsilon
      REAL(KIND = kREAL) :: e1_hat0(3), e2_hat0(3)
      REAL(KIND = kREAL) :: t(3, 3)!, t_transpose(3, 3)
      !!<YM, KKE, 2014/03/05>
      REAL(KIND = kREAL) :: strain_ela(6)
      REAL(KIND = kREAL) :: ks0
      REAL(KIND = kREAL) :: strain_tildash, sig_c_temp, sig_tilde
      REAL(KIND = kREAL) :: strain_pla_recon(6)
      !!<YM, KKE, 2014/03/05 END>
      !!<YM, KKE, 2014/03/24>
      REAL(KIND = kREAL) :: strain_ela0(6)
      REAL(KIND = kREAL) :: coeff_a, coeff_b, coeff_c, root_beta
      !!<YM, KKE, 2014/03/24 END>         
      REAL(KIND = kREAL) :: dx_hat_0(6), dsig_hat_0(6)
      !REAL(KIND = kREAL) :: tau_hat
      !REAL(KIND = kREAL) :: a(8)
      REAL(KIND = kREAL) :: emat(9, 9)
      REAL(KIND = kREAL) :: lmat_inv(9, 6)
      REAL(KIND = kREAL) :: area !,coeff
      REAL(KIND = kREAL) :: lmat(6, 9)
      REAL(KIND = kREAL) :: IIw(6, 6)
      REAL(KIND = kREAL) :: bmat(ndof*nn, 8)
      REAL(KIND = kREAL) :: pi
      INTEGER(KIND = kint) :: icel
      
      !!<YM, KKE, 2014/03/24>
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
    !added by Shan
    REAL(KIND = kREAL) :: k1init



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
      
      ! 載荷ステップ毎に更新する変数はここで処理する
      if(iter .eq. 1) then
        strain_pla0(:)=strain_pla(:)
        !gausses(1)%istatus0(1)=gausses(1)%istatus(1)
        istatus0=istatus
        strain_increment(:)=0.0d0
        
      endif
      
!--------------------------------------------------------------------
      
      pi = 4.0D0*DATAN( 1.0D0 )
      
!--------------------------------------------------------------------
      
      qf(1:3*nn) = 0.0D0
      
!--------------------------------------------------------------------
      
      matl=>gausses(1)%pMaterial
      
      mtype = matl%mtype
      
!--------------------------------------------------------------------      
       
       ! Goodman
      IF( mtype .EQ. GOODMAN ) THEN
        
        call Set_stiff_Joint3n(gausses, matl, k1, k2, k3)
        
        phi = matl%variables(M_GOODMAN_PHI)
        c   = matl%variables(M_GOODMAN_C)
        
        epsilon = matl%variables(M_GOODMAN_EPS)
        
       
       !--------------------------------------------------------

       call Set_LocalCood_Joint3n(nn, ecoord, t, area, e1_hat0, e2_hat0)       
       
       lmat_inv=0.0d0
       call Set_Bmat_Joint3n(nn, t, ecoord, lmat)

       IIw = 0.0D0
       call Set_Dmat_Joint3n(t, area, k1, k2, k3, e1_hat0, e2_hat0, IIw)    
       
       ! 以下，単位面積当たりの判定を実施するため，並進成分はareaで除する
       ! 対応してqfの計算時にareaを乗じる
       ! The stiffness per unit area is calculated
       ! Multiply by area when calculating force qf
       do i=1, 3
         do j = 1, 3
           IIw(i, j) = IIw(i, j) / area
         enddo
       enddo
               
       !--------------------------------------------------------
       
       dx_hat_0 = 0.0d0
       DO i = 1, 3
        do j = 1, 6
          dx_hat_0(j) = dx_hat_0(j) - lmat(j, i    ) * du(i, 1) &
                                    - lmat(j, i + 3) * du(i, 2) &
                                    - lmat(j, i + 6) * du(i, 3) &
                                    + lmat(j, i    ) * du(i, 4) &
                                    + lmat(j, i + 3) * du(i, 5) &
                                    + lmat(j, i + 6) * du(i, 6)
        enddo
       enddo
       
       !--------------------------------------------------------
       
       dsig_hat_0 = matmul(IIw, dx_hat_0)
       
       !--------------------------------------------------------
       !!<YM, KKE, 2014/04/17>
       ! 初期応力の設定 (intial stress)
       if(flag_initset .EQ. 0) then

         call Set_initStress_Joint3n(nn, matl, ecoord, t, &
                                     sigma_init, stress_init, strain_init, flag_initset)
              
       endif
       ! 初期応力の設定終了
       !!<YM, KKE, 2014/04/17 END>
       
       !Added by Shan 14/07/2017
       ! if(flag_initset .EQ. 2) then
             ! stress_init(1:3) = gausses(1)%stress(1:3)
             ! strain_hat=0.0D0
             ! stress_hat=0.0D0
             ! call Set_initStress_Joint3n_restart(k1,k2,k3,stress_init,         &
                                                ! strain_init,flag_initset)   
       ! endif
       
       ! if ((flag_initset .EQ. 1).and.(matl%variables(M_PLCONST1).gt.10.0)) then 
            ! stress_init(1:3) = gausses(1)%stress(1:3)
            ! strain_hat=0.0D0
            ! stress_hat=0.0D0
            ! ! call Set_initStress_Joint3n_restart(                            &
                ! ! matl%variables(M_GOODMAN_K1),matl%variables(M_GOODMAN_K2),  &
                ! ! matl%variables(M_GOODMAN_K3), stress_init,strain_init,      &
                ! ! flag_initset)
      ! end if
       
       
       
       
       if(iter .eq. 1) then
       ! 接線剛性を計算するために、前回ステップの応力値を計算しておく
       ! stress/strain from previous step for tangent stiffness calculation 
         strain_tildash0=DSQRT((strain_hat(1)+strain_init(1)-strain_pla0(1))*(strain_hat(1)+strain_init(1)-strain_pla0(1)) &
                              +(strain_hat(2)+strain_init(2)-strain_pla0(2))*(strain_hat(2)+strain_init(2)-strain_pla0(2)) )
         stress_tildash0=DSQRT((stress_hat(1)+stress_init(1))*(stress_hat(1)+stress_init(1)) &
                              +(stress_hat(2)+stress_init(2))*(stress_hat(2)+stress_init(2)) )
       endif
       
       ! 全ひずみの算出
       ! Total strain
       strain_hat(1:6) = strain_hat(1:6) + dx_hat_0(1:6) + strain_init(1:6)     

       ! 仮の要素応力の算出
       ! "仮"というのは、この後の接触・滑動判定で更新するから
       stress_hat(1:6) = stress_hat(1:6) + dsig_hat_0(1:6) + stress_init(1:6)      

       ! 各載荷ステップでの増分ひずみの更新
       ! incremental strain
       strain_increment(1:6) = strain_increment(1:6) + dx_hat_0(1:6)
       
       ! 出力用のひずみ（全ひずみ）
       !call Conv_Coord_tenstor_vec(gausses(1)%strain, t, strain_hat(1:3))
       gausses(1)%strain(1:3) = strain_hat(1:3)
       
       ! 増分後の要素応力からモールクーロンの判定を実施する
       If(istatus .EQ. 3) then         !破壊時の判定(今回は適用しない)
         ! 全ひずみ＝塑性ひずみ→応力ゼロ
         strain_pla(1:2)=strain_hat(1:2)
         istatus=3
        
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
         
         ! !Added by Shan         
         ! if (matl%variables(M_PLCONST1).lt.0) then
            ! ks0=(matl%variables(M_PLCONST2)*((abs(gausses(1)%fstatus(60))) &
                    ! **matl%variables(M_PLCONST3)))/(2*(1+matl%variables(M_DASHPOD_VP)))
        ! end if
         
         If(istatus0 .EQ. 1) then
           !前回載荷ステップが剥離だった場合（再接触時）は塑性ひずみを計算する
           !この後の滑動・非滑動の判定はここで定義する塑性ひずみを用いて行う
           strain_pla_recon(1:6)=strain_hat(1:6)-dabs(strain_hat(3)/strain_increment(3))*strain_increment(1:6)
           ! write(901, "(2I5, 100E12.3)") iter, istatus, strain_hat, strain_increment
           
           !弾性ひずみは更新した塑性ひずみから再計算する
           strain_pla(1:6) = strain_pla_recon(1:6)
           strain_ela(1:6)=strain_hat(1:6)-strain_pla(1:6)
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
       
       ! 判定が終わったので戻す
       stress_hat(:) = stress_hat(:) - stress_init(:)
       strain_hat(:) = strain_hat(:) - strain_init(:)
                                                                                                               
       ! 出力用の応力
       !call Conv_Coord_tenstor_vec(gausses(1)%stress, t, stress_hat(1:3))
       gausses(1)%stress(1:3) = stress_hat(1:3)
       ! gausses(1)%stress(4) = k1
       ! gausses(1)%stress(5) = gausses(1)%fstatus(60)
       ! gausses(1)%stress(6) = matl%variables(M_PLCONST2)
       
       gausses(1)%strain_out(1:6) = gausses(1)%strain(1:6)
       gausses(1)%stress_out(1:6) = gausses(1)%stress(1:6)
              
       bmat = 0.0D0
       
       DO jsize = 1, 6
         
        DO isize = 1, 9
         ! bmat = TRANSPOSE( lmat )
         bmat(  isize, jsize) = -lmat(jsize, isize)
         bmat(9+isize, jsize) =  lmat(jsize, isize)
         
        END DO
        
       END DO
                 
       
       DO jsize = 1, 3
        
        DO isize = 1, nn*3
         
         ! 並進の成分のみ面積を乗じる（bマトに応力→力の変換を含めていない）
         qf(isize) = qf(isize)+bmat(isize, jsize  )*stress_hat(jsize  )*area
         qf(isize) = qf(isize)+bmat(isize, jsize+3)*stress_hat(jsize+3) 
          
        END DO
        
       END DO
       
      END IF
      ! dashpodの処理は記載しないが、qfが初期化されているので何もしないとゼロになる
      !WRITE(3019, '( (A, 1X), 2(I5, 1X), 100(E12.5, 1X) )') 'joint_damp qf :', icel, iter, (qf(i),i=1,18)
!--------------------------------------------------------------------      

      RETURN
      
!####################################################################
      END SUBROUTINE Update_Joint3n
!####################################################################

!####################################################################
      SUBROUTINE Update_Joint3n_damp                            &
                 (icel, etype, nn, ecoord, u, du, qf, gausses, &
                  acc, vec3, ddux, &
                  iter)                      
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
      
      INTEGER(KIND = kint), INTENT(IN) :: iter            ! iter
      
!--------------------------------------------------------------------
      
      TYPE( tMaterial ), POINTER :: matl
      
      INTEGER(KIND = kint), PARAMETER :: ndof = 3
      INTEGER(KIND = kint) :: i, j
      INTEGER(KIND = kint) :: na, nb, nc, nd
      INTEGER(KIND = kint) :: isize, jsize, ksize, lsize
      INTEGER(KIND = kint) :: mtype
      
      REAL(KIND = kREAL) :: k1, k2, k3
      REAL(KIND = kREAL) :: e1_hat0(3), e2_hat0(3)
      REAL(KIND = kREAL) :: t(3, 3)!, t_transpose(3, 3)
      REAL(KIND = kREAL) :: emat(9, 9)
      REAL(KIND = kREAL) :: lmat_inv(9, 9)
      REAL(KIND = kREAL) :: coeff, area
      REAL(KIND = kREAL) :: lmat(9, 9)
      REAL(KIND = kREAL) :: IIw(9, 9)
      INTEGER(KIND = kint) :: icel
      
      REAL(KIND = kREAL) :: Vs, Vp, rho
      REAL(KIND = kREAL) :: stiff(ndof*nn, ndof*nn)
      REAL(KIND = kREAL) :: rayleigh_alpha, rayleigh_beta, delta_t
      REAL(KIND = kREAL) :: newmark_beta, newmark_gamma
      REAL(KIND = kREAL) :: coeff_damp
      
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
       
       !--------------------------------------------------------
         
         call Set_param_Joint3n_dp(matl, Vp, Vs, rho, rayleigh_alpha, rayleigh_beta, &
                                   delta_t, newmark_beta, newmark_gamma, coeff_damp)
               
       !--------------------------------------------------------
       
       call Set_LocalCood_Joint3n(nn, ecoord, t, area, e1_hat0, e2_hat0)
              
       lmat=0.0d0
       call Set_Bmat_Joint3n_dp(nn, t, ecoord, lmat)
       
       IIw = 0.0D0
       Call Set_Dmat_Joint3n_dp(area, rho, Vs, Vp, IIw)
       
       emat = MATMUL( TRANSPOSE( lmat ) , MATMUL( IIw, lmat))
      
       !--------------------------------------------------------
      
       DO jsize = 1, 9
        
        DO isize = 1, 9
         
         stiff(isize  , jsize  ) =  emat(isize, jsize)
         stiff(isize+9, jsize  ) = -emat(isize, jsize)
         stiff(isize  , jsize+9) = -emat(isize, jsize)
         stiff(isize+9, jsize+9) =  emat(isize, jsize)
         
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
      END SUBROUTINE Update_Joint3n_damp
!####################################################################

!####################################################################
       SUBROUTINE Set_stiff_Joint3n                           &            
                  (gausses, matl, k1, k2, k3)
!####################################################################                  
                  
       ! USE
       USE mMechGauss
       USE m_MatMatrix
       USE m_common_struct
       USE mMaterial
                  
 !--------------------------------------------------------------------
      
       ! TYPE(tGaussStatus), INTENT(IN) :: gausses(:)                ! status of Gaussian qudrature points
       TYPE( tMaterial ), POINTER :: matl
       REAL(KIND = kREAL), INTENT(OUT) ::k1, k2, k3 

       INTEGER(KIND = kint), POINTER :: istatus
       REAL(KIND = kREAL), POINTER :: stiff_coeff
       
       !Added by Shan 13/07/2017
       TYPE(tGaussStatus) :: gausses(:)                ! status of Gaussian qudrature points
       REAL(KIND = kREAL) :: alpha,beta,stress_mean,poisson
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
    !Added by Shan 13/07/2017
    ! if (matl%variables(M_PLCONST1).lt.0) then
        ! alpha   =  matl%variables(M_PLCONST2)
        ! beta    =  matl%variables(M_PLCONST3)
        ! poisson =  matl%variables(M_DASHPOD_VP)
        ! if (gausses(1)%istatus(3)==1) then
            ! !currently stress mean = sigmaZ
            ! stress_mean = gausses(1)%stress(3)
            ! gausses(1)%fstatus(60) = stress_mean
            ! gausses(1)%istatus(3)  = 2
        ! end if
        
        ! call Set_stiff_Joint3n_restart(alpha,beta,poisson,gausses(1)%istatus(1),&   
                ! gausses(1)%fstatus(60),stiff_coeff,k1, k2, k3)
        
    ! else
       IF( istatus .EQ. 0 ) THEN
         
         k1 = matl%variables(M_GOODMAN_K1)
         k2 = matl%variables(M_GOODMAN_K2)
         k3 = matl%variables(M_GOODMAN_K3)
         
       ELSE IF( istatus .EQ. 1 ) THEN
         
         k1 = matl%variables(M_GOODMAN_K1)*0.001D+0
         k2 = matl%variables(M_GOODMAN_K2)*0.001D+0
         k3 = matl%variables(M_GOODMAN_K3)*0.001D+0
         
       ELSE IF( istatus .EQ. 2 ) THEN
         
         k1 = matl%variables(M_GOODMAN_K1)*stiff_coeff
         k2 = matl%variables(M_GOODMAN_K2)*stiff_coeff
         k3 = matl%variables(M_GOODMAN_K3)
         
       ELSE IF( istatus .EQ. 3 ) THEN
         
         k1 = 0.0D0
         k2 = 0.0D0
         k3 = matl%variables(M_GOODMAN_K3)
         
       END IF
       
    ! end if
        

!####################################################################     
       END SUBROUTINE Set_stiff_Joint3n
!####################################################################

!####################################################################
       SUBROUTINE Set_stiff_Joint3n_restart                         &
                  (alpha,beta,poisson,istatus,stress_mean,stiff_coeff,k1, k2, k3)
!####################################################################                  
    !Added by Shan 13/07/2017
    real(kind = kreal), intent(in)  :: alpha,beta,stiff_coeff,poisson
    real(kind = kreal), intent(out) :: k1,k2,k3
    integer(kind = kint), intent(in) :: istatus
    real(kind = kreal)  :: stress_mean
    real(kind = kreal)  :: ey,g
    
    if (abs(stress_mean).le.1.0D-03) stress_mean = 1.0D0
    ey = (alpha * (abs(stress_mean)**beta))
    g=ey/(2*(1+poisson))
    k1 = g
    k2 = g
    k3 = ey
    
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
       END SUBROUTINE Set_stiff_Joint3n_restart
!####################################################################

!####################################################################                  
      SUBROUTINE Set_Bmat_Joint3n                           &            
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
      REAL(KIND = kREAL), INTENT(OUT)  :: lmat(6, 9)
      
!--------------------------------------------------------------------                 
                  
      REAL(KIND = kREAL) :: a1(3), a2(3), a3(3), aw(3)
      REAL(KIND = kREAL) :: lmat_inv(9, 6)
      REAL(KIND = kREAL) :: ltxxl(6, 6)
      REAL(KIND = kREAL) :: ltxxl_i(6, 6)
      REAL(KIND = kREAL) :: coeff
      INTEGER(KIND = kint) :: i, j

!--------------------------------------------------------------------
      
     Do i=1,3
       lmat_inv(i  , 1) = t(1, i)
       lmat_inv(i+3, 1) = t(1, i)
       lmat_inv(i+6, 1) = t(1, i)
       lmat_inv(i  , 2) = t(2, i)
       lmat_inv(i+3, 2) = t(2, i)
       lmat_inv(i+6, 2) = t(2, i)
       lmat_inv(i  , 3) = t(3, i)
       lmat_inv(i+3, 3) = t(3, i)
       lmat_inv(i+6, 3) = t(3, i)
     enddo
     
     DO i=1, 3
        a1(i)=(ecoord(i, 1)+ecoord(i, 4))*0.5d0
        a2(i)=(ecoord(i, 2)+ecoord(i, 5))*0.5d0
        a3(i)=(ecoord(i, 3)+ecoord(i, 6))*0.5d0
     END DO
     DO i=1, 3
        aw(i)=(a1(i)+a2(i)+a3(i))/3.0d0
     END DO
     DO i=1, 3
        a1(i)=a1(i)-aw(i)
        a2(i)=a2(i)-aw(i)
        a3(i)=a3(i)-aw(i)
     ENDDO

     lmat_inv(1, 4) = a1(2)*t(1, 3) - a1(3)*t(1, 2)
     lmat_inv(2, 4) = a1(3)*t(1, 1) - a1(1)*t(1, 3)
     lmat_inv(3, 4) = a1(1)*t(1, 2) - a1(2)*t(1, 1)
     lmat_inv(4, 4) = a2(2)*t(1, 3) - a2(3)*t(1, 2)
     lmat_inv(5, 4) = a2(3)*t(1, 1) - a2(1)*t(1, 3)
     lmat_inv(6, 4) = a2(1)*t(1, 2) - a2(2)*t(1, 1)
     lmat_inv(7, 4) = a3(2)*t(1, 3) - a3(3)*t(1, 2)
     lmat_inv(8, 4) = a3(3)*t(1, 1) - a3(1)*t(1, 3)
     lmat_inv(9, 4) = a3(1)*t(1, 2) - a3(2)*t(1, 1)
       
     lmat_inv(1, 5) = a1(2)*t(2, 3) - a1(3)*t(2, 2)
     lmat_inv(2, 5) = a1(3)*t(2, 1) - a1(1)*t(2, 3)
     lmat_inv(3, 5) = a1(1)*t(2, 2) - a1(2)*t(2, 1)
     lmat_inv(4, 5) = a2(2)*t(2, 3) - a2(3)*t(2, 2)
     lmat_inv(5, 5) = a2(3)*t(2, 1) - a2(1)*t(2, 3)
     lmat_inv(6, 5) = a2(1)*t(2, 2) - a2(2)*t(2, 1)
     lmat_inv(7, 5) = a3(2)*t(2, 3) - a3(3)*t(2, 2)
     lmat_inv(8, 5) = a3(3)*t(2, 1) - a3(1)*t(2, 3)
     lmat_inv(9, 5) = a3(1)*t(2, 2) - a3(2)*t(2, 1)
       
     lmat_inv(1, 6) = a1(2)*t(3, 3) - a1(3)*t(3, 2)
     lmat_inv(2, 6) = a1(3)*t(3, 1) - a1(1)*t(3, 3)
     lmat_inv(3, 6) = a1(1)*t(3, 2) - a1(2)*t(3, 1)
     lmat_inv(4, 6) = a2(2)*t(3, 3) - a2(3)*t(3, 2)
     lmat_inv(5, 6) = a2(3)*t(3, 1) - a2(1)*t(3, 3)
     lmat_inv(6, 6) = a2(1)*t(3, 2) - a2(2)*t(3, 1)
     lmat_inv(7, 6) = a3(2)*t(3, 3) - a3(3)*t(3, 2)
     lmat_inv(8, 6) = a3(3)*t(3, 1) - a3(1)*t(3, 3)
     lmat_inv(9, 6) = a3(1)*t(3, 2) - a3(2)*t(3, 1)
       
     ltxxl = MATMUL( TRANSPOSE( lmat_inv ), lmat_inv)
      
     ltxxl_i = 0.0D0
     ltxxl_i(1, 1) = 1.0D0 / ltxxl(1, 1)
     ltxxl_i(2, 2) = 1.0D0 / ltxxl(2, 2)
     ltxxl_i(3, 3) = 1.0D0 / ltxxl(3, 3)
     ltxxl_i(6, 6) = 1.0D0 / ltxxl(6, 6)
       
     ! 回転部分のdet
     coeff = ltxxl(4, 4) * ltxxl(5, 5) - ltxxl(4, 5) * ltxxl(5, 4)
     ! 2x2の逆行列
     ltxxl_i(4, 4) =  ltxxl(5, 5) / coeff
     ltxxl_i(5, 5) = ltxxl(4, 4) / coeff
     ltxxl_i(4, 5) = -1.0D0 * ltxxl(5, 4) / coeff
     ltxxl_i(5, 4) = ltxxl_i(4, 5)
     
     lmat = MATMUL(ltxxl_i, TRANSPOSE( lmat_inv ))

!####################################################################           
     end SUBROUTINE Set_Bmat_Joint3n
!####################################################################      

!####################################################################                  
      SUBROUTINE Set_Bmat_Joint3n_dp                         &            
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
      REAL(KIND = kREAL), INTENT(OUT)  :: lmat(3, 9)
      
!--------------------------------------------------------------------                 
                  
      REAL(KIND = kREAL) :: lmat_inv(9, 3)
      INTEGER(KIND = kint) :: i, j

!--------------------------------------------------------------------
      
     Do i=1,3
       lmat_inv(i  , 1) = t(1, i)
       lmat_inv(i+3, 1) = t(1, i)
       lmat_inv(i+6, 1) = t(1, i)
       lmat_inv(i  , 2) = t(2, i)
       lmat_inv(i+3, 2) = t(2, i)
       lmat_inv(i+6, 2) = t(2, i)
       lmat_inv(i  , 3) = t(3, i)
       lmat_inv(i+3, 3) = t(3, i)
       lmat_inv(i+6, 3) = t(3, i)
     enddo
          
     lmat = TRANSPOSE( lmat_inv )

!####################################################################           
     end SUBROUTINE Set_Bmat_Joint3n_dp
!####################################################################      

!####################################################################                  
       SUBROUTINE Set_Dmat_Joint3n                           &            
                  (t, area, k1, k2, k3, e1_hat0, e2_hat0, IIw)
!####################################################################                  
                  
       ! USE       
                  
 !--------------------------------------------------------------------
      
       REAL(KIND = kREAL), INTENT(IN)  :: t(3, 3)
       REAL(KIND = kREAL), INTENT(IN)  :: area
       REAL(KIND = kREAL), INTENT(IN)  :: k1, k2, k3
       REAL(KIND = kREAL), INTENT(IN)  :: e1_hat0(3), e2_hat0(3)
       REAL(KIND = kREAL), INTENT(OUT) :: IIw(6, 6)
       
!--------------------------------------------------------------------
       
       REAL(KIND = kREAL) :: aa1, aa2, bb1, bb2
       
       
!--------------------------------------------------------------------     
       !要素に沿った座標系で三角形を表現する
       aa1=t(1,1)*e1_hat0(1)+t(1,2)*e1_hat0(2)+t(1,3)*e1_hat0(3)
       aa2=t(2,1)*e1_hat0(1)+t(2,2)*e1_hat0(2)+t(2,3)*e1_hat0(3) !この成分は必ずゼロになるはず
       bb1=t(1,1)*e2_hat0(1)+t(1,2)*e2_hat0(2)+t(1,3)*e2_hat0(3)
       bb2=t(2,1)*e2_hat0(1)+t(2,2)*e2_hat0(2)+t(2,3)*e2_hat0(3)
       
       ! 剛性マトリックス(D)の定義
       ! 応力→力の変換も含む
       IIw(1,1)=k1*area
       IIw(2,2)=k2*area
       IIw(3,3)=k3*area
       ! 回転剛性テンソルの定義
       IIw(4,4)=area*k3*((      aa1*aa1 + aa1*bb1 + bb1*bb1                ) / 12.0D0 - (aa1+bb1)*(aa1+bb1)/18.0D0)
       IIw(4,5)=area*k3*((2.0D0*aa1*aa2 + aa2*bb1 + aa1*bb2 + 2.0D0*bb1*bb2) / 24.0D0 - (aa1+bb1)*(aa2+bb2)/18.0D0)
       IIw(5,4)=IIw(4,5)
       IIw(5,5)=area*k3*((      aa2*aa2 + aa2*bb2 + bb2*bb2                ) / 12.0D0 - (aa2+bb2)*(aa2+bb2)/18.0D0)
       
       IIw(6,6)=area*k2*((      aa1*aa1 + aa1*bb1 + bb1*bb1                ) / 12.0D0 - (aa1+bb1)*(aa1+bb1)/18.0D0) & 
               +area*k1*((      aa2*aa2 + aa2*bb2 + bb2*bb2                ) / 12.0D0 - (aa2+bb2)*(aa2+bb2)/18.0D0)

!####################################################################     
       END SUBROUTINE set_Dmat_Joint3n
!####################################################################

!####################################################################                  
       SUBROUTINE Set_Dmat_Joint3n_dp                        &            
                  (area, rho, Vs, Vp, IIw)
!####################################################################                  
                  
       ! USE       
                  
 !--------------------------------------------------------------------
      
       REAL(KIND = kREAL), INTENT(IN)  :: area
       REAL(KIND = kREAL), INTENT(IN)  :: rho, Vs, Vp
       REAL(KIND = kREAL), INTENT(OUT) :: IIw(9, 9)
       
!--------------------------------------------------------------------
       
       ! 剛性マトリックス(D)の定義
       IIw(1,1)=rho * Vs * area * (1.0d0 / 3.0D0)
       IIw(2,2)=rho * Vs * area * (1.0d0 / 3.0D0)
       IIw(3,3)=rho * Vs * area * (1.0d0 / 3.0D0)
       IIw(4,4)=rho * Vs * area * (1.0d0 / 3.0D0)
       IIw(5,5)=rho * Vs * area * (1.0d0 / 3.0D0)
       IIw(6,6)=rho * Vs * area * (1.0d0 / 3.0D0)
       IIw(7,7)=rho * Vp * area * (1.0d0 / 3.0D0)
       IIw(8,8)=rho * Vp * area * (1.0d0 / 3.0D0)
       IIw(9,9)=rho * Vp * area * (1.0d0 / 3.0D0)
       

!####################################################################     
       END SUBROUTINE set_Dmat_Joint3n_dp
!####################################################################

!####################################################################                  
       SUBROUTINE Set_LocalCood_Joint3n                           &            
                  (nn, ecoord, t, area, e1_hat0, e2_hat0)
!####################################################################                  
                  
       ! USE       
                  
 !--------------------------------------------------------------------
       
       INTEGER(KIND = kint), INTENT(IN) :: nn                      ! number of elemental nodes
       REAL(KIND = kREAL), INTENT(IN)   :: ecoord(3, nn)           ! coordinates of elemental nodes
       REAL(KIND = kREAL), INTENT(OUT)  :: t(3, 3)
       REAL(KIND = kREAL), INTENT(OUT)  :: area
       REAL(KIND = kREAL), INTENT(OUT)  :: e1_hat0(3), e2_hat0(3)
       
!--------------------------------------------------------------------         
       
       REAL(KIND = kREAL) :: e1_hat(3), e2_hat(3), e3_hat(3)
       REAL(KIND = kREAL) :: e1_hat_abs, e2_hat_abs, e3_hat_abs
       REAL(KIND = kREAL) :: x21(3), x31(3)
       REAL(KIND = kREAL) :: x54(3), x64(3)
       INTEGER(KIND = kint) :: i, j
       
!--------------------------------------------------------------------
       DO i = 1, 3
        
        x21(i) = ecoord(i, 2)-ecoord(i, 1)
        x31(i) = ecoord(i, 3)-ecoord(i, 1)
        
        x54(i) = ecoord(i, 5)-ecoord(i, 4)
        x64(i) = ecoord(i, 6)-ecoord(i, 4)
        
       END DO
       
       !--------------------------------------------------------
       
       ! 要素に沿った座標系を定義
       DO i = 1, 3
        e1_hat0(i) = 0.50D0*( x21(i)+x54(i) )
        e2_hat0(i) = 0.50D0*( x31(i)+x64(i) )
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
       
       ! 要素の面積は要素の二辺のベクトル×0.5 
       area=e3_hat_abs*0.5D0
       
       ! 正規化
       DO i = 1, 3
        
        e1_hat(i) = e1_hat(i)/e1_hat_abs
        e2_hat(i) = e2_hat(i)/e2_hat_abs
        e3_hat(i) = e3_hat(i)/e3_hat_abs
        
       END DO

       DO i = 1, 3
        
        t(1, i) = e1_hat(i)
        t(2, i) = e2_hat(i)
        t(3, i) = e3_hat(i)
        
       END DO

!####################################################################     
       END SUBROUTINE Set_LocalCood_Joint3n
!####################################################################

!####################################################################                  
       SUBROUTINE Set_param_Joint3n_dp                           &            
                  (matl, Vp, Vs, rho, rayleigh_alpha, rayleigh_beta, &
                   delta_t, newmark_beta, newmark_gamma, coeff_damp)
!####################################################################                  
                  
       ! USE
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
       
!--------------------------------------------------------------------
       
       REAL(KIND = kREAL) :: param_a3, param_b3
       REAL(KIND = kREAL) :: param_a1, param_b1, param_a2, param_b2
       
!--------------------------------------------------------------------     
       
       Vp = matl%variables(M_DASHPOD_VP)
       Vs = matl%variables(M_DASHPOD_VS)
       rho = matl%variables(M_GOODMAN_C)
         
       rayleigh_alpha   = matl%variables(M_GOODMAN_K1)
       rayleigh_beta    = matl%variables(M_GOODMAN_K2)
       delta_t = matl%variables(M_GOODMAN_K3)
         
       newmark_beta = 1.0D0 / 4.0D0
       newmark_gamma =  0.5D0
               
       coeff_damp = param_b3 / (1.0d0 + rayleigh_beta * param_b3)
         
!####################################################################     
       END SUBROUTINE Set_param_Joint3n_dp
!####################################################################     
!####################################################################                  
       SUBROUTINE Set_initStress_Joint3n                        &            
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
       REAL(KIND = kREAL), INTENT(OUT)  :: stress_init(8)
       REAL(KIND = kREAL), INTENT(OUT)  :: strain_init(8)
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
        
         ! 鉛直座標の平均値を所得
         ave_cood_z = 0.125d0 * (  ecoord(3, 1) + ecoord(3, 2) + ecoord(3, 3) + ecoord(3, 4) &
                                 + ecoord(3, 5) + ecoord(3, 6) + ecoord(3, 7) + ecoord(3, 8) )
       
         !--------------------------------------------------------
        
         sigma_init(:,:)=0.0d0
        
         sigma_init(3,3) = coeff_a * ave_cood_z + coeff_b
         sigma_init(1,1) = coeff_c * sigma_init(3,3)
         sigma_init(2,2) = coeff_c * sigma_init(3,3)
        
        
         ! 局所座標系に変換
         sigma_init(   1:3, :) = MATMUL( t, sigma_init( :, :))
         sigma_init( :, :) = MATMUL( sigma_init( :, :), TRANSPOSE(t) )
        
         !--------------------------------------------------------
        
         stress_init(1)=sigma_init(1, 3)
         stress_init(2)=sigma_init(2, 3)
         stress_init(3)=sigma_init(3, 3)
         stress_init(4)= 0.0d0
         stress_init(5)= 0.0d0
         stress_init(6)= 0.0d0
        
         strain_init(1)=stress_init(1) / k1_init
         strain_init(2)=stress_init(2) / k2_init
         strain_init(3)=stress_init(3) / k3_init
         strain_init(4)=0.0d0
         strain_init(5)=0.0d0
         strain_init(6)=0.0d0
       
         ! 初期応力設定済みフラグ
         flag_initset = 1
        
         ! 一応初期化
         coeff_a=0.0d0
         coeff_b=0.0d0
         coeff_c=0.0d0
       

!####################################################################     
       END SUBROUTINE Set_initStress_Joint3n
!####################################################################    

!####################################################################                  
    SUBROUTINE Set_initStress_Joint3n_restart(k1,k2,k3,          &
                        stress_init,strain_init,flag_initset) 
!#################################################################### 
    !Added by Shan 14/07/2017
    
     !--------------------------------------------------------------------
      
    REAL(KIND = kREAL), INTENT(IN)   :: k1,k2,k3
    REAL(KIND = kREAL)               :: stress_init(8),strain_init(8)
    INTEGER(KIND = kint)             :: flag_initset
       
    !--------------------------------------------------------------------
    strain_init=0.0D0
    
    strain_init(1)=stress_init(1) / k1
    strain_init(2)=stress_init(2) / k2
    strain_init(3)=stress_init(3) / k3
    strain_init(4)=0.0d0
    strain_init(5)=0.0d0
    strain_init(6)=0.0d0
    flag_initset = 3
!####################################################################     
       END SUBROUTINE Set_initStress_Joint3n_restart
!####################################################################  


!####################################################################                  
       SUBROUTINE Conv_Coord_tenstor_vec                        &            
                  (tensor_out_6, t, vector_3)
!####################################################################                  
                  
       ! USE       
                  
 !--------------------------------------------------------------------
      
       REAL(KIND = kREAL), INTENT(OUT)  :: tensor_out_6(6)
       REAL(KIND = kREAL), INTENT(IN)  :: t(3, 3)
       REAL(KIND = kREAL), INTENT(IN) :: vector_3(3)
       
!--------------------------------------------------------------------
       
       ! 剛性マトリックス(D)の定義
       tensor_out_6(1) = vector_3(1)*t(1,1)*t(1,1) &
                        +vector_3(2)*t(2,1)*t(2,1) &
                        +vector_3(3)*t(3,1)*t(3,1) 
       tensor_out_6(2) = vector_3(1)*t(1,2)*t(1,2) &
                        +vector_3(2)*t(2,2)*t(2,2) &
                        +vector_3(3)*t(3,2)*t(3,2) 
       tensor_out_6(3) = vector_3(1)*t(1,3)*t(1,3) &
                        +vector_3(2)*t(2,3)*t(2,3) &
                        +vector_3(3)*t(3,3)*t(3,3) 
       tensor_out_6(4) = vector_3(1)*t(1,1)*t(1,2) &
                        +vector_3(2)*t(2,1)*t(2,2) &
                        +vector_3(3)*t(3,1)*t(3,2) 
       tensor_out_6(5) = vector_3(1)*t(1,2)*t(1,3) &
                        +vector_3(2)*t(2,2)*t(2,3) &
                        +vector_3(3)*t(3,2)*t(3,3) 
       tensor_out_6(6) = vector_3(1)*t(1,3)*t(1,1) &
                        +vector_3(2)*t(2,3)*t(2,1) &
                        +vector_3(3)*t(3,3)*t(3,1) 
       

!####################################################################     
    END SUBROUTINE Conv_Coord_tenstor_vec
!####################################################################
    END MODULE m_static_LIB_Joint3n



module m_RO
    
  !USE hecmw, ONLY : kint, kREAL
  use mMaterial
  use m_ElasticLinear
  USE mMechGauss

  implicit none

  contains
  
!####################################################################
   subroutine ROMatrix( gauss, stress, istat, fstat, D, strain)
!####################################################################
   
     use hecmw
     use mMaterial
     implicit none
     type( tGaussStatus ), intent(IN) :: gauss          !> status of qudrature point
     REAL(KIND=kreal), INTENT(IN)  :: stress(6) !< stress
     INTEGER, INTENT(IN)           :: istat     !< plastic state
     REAL(KIND=kreal), INTENT(IN)  :: fstat(:)  !< plastic strain, back stress
     REAL(KIND=kreal), INTENT(OUT) :: D(:,:)    !< strain-stress relation
     REAL(KIND=kreal), INTENT(IN), optional :: strain(6) !< strain
	
     integer :: i
     real(kind=kreal) :: pi
     real(kind=kreal) :: EE, PP, kk, C, theta, alpha, beta, g_max, rho
     real(kind=kreal) :: sigma_m, deriv(6), deriv_e, tau_max
     real(kind=kreal) :: ita_e, xi_e, stiff_g
     TYPE( tMaterial ):: matl      !< material properties
     INTEGER(kind=kint) :: mtype
     
     REAL(KIND = kREAL) :: C1_0, C2_0, C1_inf, C2_inf
     REAL(KIND = kREAL) :: strain_m, gamma_r, gamma
     
     REAL(KIND = kREAL) :: target_damp, gamma_r_sklt, gamma_sklt,stress_sklt
     REAL(KIND = kREAL) :: C1_gamma, C2_gamma, c1_g_i, c2_g_i, stiff, temp1

! ----------------------------------------------------------------------------
     
     INTEGER(KIND = kint), POINTER :: stress_zero_index, strain_zero_index
     REAL(KIND = kREAL), POINTER :: stress_zero(:), strain_zero(:)
     !REAL(KIND = kREAL), POINTER :: strain_max
     
! ----------------------------------------------------------------------------
     matl = gauss%pMaterial
     mtype = matl%mtype
     
     !  pointer設定 - 不要なものはコメントアウトする
     if(mtype == RO) then
       stress_zero_index => gauss%istatus(1)
       stress_zero       => gauss%fstatus(1:350)
     elseif(mtype == GHE) then
       stress_zero_index => gauss%istatus(1)
       strain_zero_index => gauss%istatus(1)
       stress_zero       => gauss%fstatus(1:350)
       strain_zero       => gauss%fstatus(351:700)
       !strain_max        => gauss%fstatus(701)
     endif

! ----------------------------------------------------------------------------

     pi = 4.0D0*DATAN( 1.0D0 )
     
     D(:,:)=0.d0
     !write(*,*) EE, PP
     !ee = matl%variables(M_YOUNGS)
     !pp = matl%variables(M_POISSON)
     !kk = ee / 3.0d0 / (1.0d0 - 2.0d0 * pp)
     !rho = matl%variables(M_DENSITY)
     !C  = matl%variables(M_PLCONST1)                                  !ROパラメータ第一引数
     !theta = 2.0d0 * pi * matl%variables(M_PLCONST2) / 360.0d0        !ROパラメータ第二引数
     !alpha = matl%variables(M_PLCONST3)                               !ROパラメータ第三引数
     !beta = matl%variables(M_PLCONST4)                                !ROパラメータ第四引数
     if(mtype == RO) then
       call Set_param_RO_material(matl, ee, pp, kk, rho, C, theta, alpha, beta)
     elseif(mtype == GHE) then
       call Set_param_RO_material(matl, ee, pp, kk, rho, C, theta, alpha, beta, C1_0, C2_0, C1_inf, C2_inf)
     endif
     
     g_max = ee / 2.0d0 / (1.0d0 + pp)
     

     !matl%stress_zero() ! 除荷応力
     
     sigma_m = (1.0d0 / 3.0d0 ) * (stress(1) + stress(2) + stress(3))
     tau_max = C * dcos(theta) + sigma_m * dsin(theta)
     
     if(mtype == RO) then
       
       deriv(1) = stress(1) - sigma_m - stress_zero(s0_index(stress_zero_index, 1))
       deriv(2) = stress(2) - sigma_m - stress_zero(s0_index(stress_zero_index, 2))
       deriv(3) = stress(3) - sigma_m - stress_zero(s0_index(stress_zero_index, 3))
       deriv(4) = stress(4)           - stress_zero(s0_index(stress_zero_index, 4))
       deriv(5) = stress(5)           - stress_zero(s0_index(stress_zero_index, 5))
       deriv(6) = stress(6)           - stress_zero(s0_index(stress_zero_index, 6))
       deriv_e = dsqrt((1.0d0 / 2.0d0) * (deriv(1)*deriv(1) + deriv(2)*deriv(2) + deriv(3)*deriv(3) + &
                                2.0d0 * (deriv(4)*deriv(4) + deriv(5)*deriv(5) + deriv(6)*deriv(6))))
       ita_e = deriv_e / tau_max
       xi_e = ita_e * (1.0d0 + alpha * ita_e ** (beta - 1.0d0))
     
       if(stress_zero_index == 1) then
         stiff_g =         g_max * (1.0d0 + alpha * beta           * ita_e ** (beta - 1.0d0)) ** (-1.0d0)
       else
         stiff_g = 1.0d0 * g_max * (1.0d0 + alpha * beta * (0.5d0 * ita_e) ** (beta - 1.0d0)) ** (-1.0d0)
       endif
     !write(987, '(A,100E10.3)') "A", g_max, alpha, beta, ita_e, stiff_g
     !write(986, '(A,100E10.3)') "A", deriv(1:6), deriv_e
     
     elseif (mtype == GHE) then

       strain_m = (1.0d0 / 3.0d0 ) * (strain(1) + strain(2) + strain(3))

       deriv(1) = strain(1) - strain_m - strain_zero(s0_index(strain_zero_index, 1))
       deriv(2) = strain(2) - strain_m - strain_zero(s0_index(strain_zero_index, 2))
       deriv(3) = strain(3) - strain_m - strain_zero(s0_index(strain_zero_index, 3))
       deriv(4) = strain(4)            - strain_zero(s0_index(strain_zero_index, 4))
       deriv(5) = strain(5)            - strain_zero(s0_index(strain_zero_index, 5))
       deriv(6) = strain(6)            - strain_zero(s0_index(strain_zero_index, 6))
       deriv_e = dsqrt((1.0d0 / 2.0d0) * (deriv(1)*deriv(1) + deriv(2)*deriv(2) + deriv(3)*deriv(3) + &
                                2.0d0  * (deriv(4)*deriv(4) + deriv(5)*deriv(5) + deriv(6)*deriv(6))))

       ! 骨格曲線上の場合は本来の基準ひずみ，履歴曲線上の場合は変更した基準ひずみを用いる
       gamma_r_sklt = tau_max / g_max
       if(strain_zero_index == 1) then
       !if(.false.) then ! for debug
         stiff_g = stress_GHE_deriv(deriv_e, alpha, beta, C1_0, C2_0, C1_inf, C2_inf, g_max, gamma_r_sklt)
         !write(806, "(100E12.3)") deriv_e, stiff_g, g_max, gamma_r_sklt
       else
         
         gamma_sklt = strain_zero(s0_index(1, 7))
         stress_sklt = stress_GHE(gamma_sklt, alpha, beta, &
                             C1_0, C2_0, C1_inf, C2_inf, g_max, &
                             gamma_r_sklt, C1_gamma, C2_gamma)
         c1_g_i = 1.0d0 / C1_gamma
         c2_g_i = 1.0d0 / C2_gamma
         gamma_r = strain_zero(s0_index(strain_zero_index, 7))
         stiff = (stress_sklt / gamma_sklt) * (c1_g_i + c2_g_i * (gamma_sklt / gamma_r))
         temp1 = c1_g_i + c2_g_i * (deriv_e * 0.5d0 / gamma_r) ! Masing則
         stiff_g = stiff / (C1_gamma * temp1 * temp1)

       endif

     endif
     !write(801, "(I5, 100E12.3)") strain_zero_index, stiff_g
     
     !D(1,1)= kk + 4.0d0 * stiff_g / 3.0d0
     !D(1,2)= kk - 2.0d0 * stiff_g / 3.0d0
     !D(1,3)= kk - 2.0d0 * stiff_g / 3.0d0
     !D(2,1)= D(1,2)
     !D(2,2)= kk + 4.0d0 * stiff_g / 3.0d0
     !D(2,3)= kk - 2.0d0 * stiff_g / 3.0d0
     !D(3,1)= D(1,3)
     !D(3,2)= D(2,3)
     !D(3,3)= kk + 4.0d0 * stiff_g / 3.0d0
     !D(4,4)= stiff_g
     !D(5,5)= stiff_g
     !D(6,6)= stiff_g
     ! ひずみと応力の方向が一致するという仮定（のはず）
     D(1,1)= kk + 2.0d0 * stiff_g / 3.0d0
     D(1,2)= kk - 1.0d0 * stiff_g / 3.0d0
     D(1,3)= kk - 1.0d0 * stiff_g / 3.0d0
     D(2,1)= D(1,2)
     D(2,2)= kk + 2.0d0 * stiff_g / 3.0d0
     D(2,3)= kk - 1.0d0 * stiff_g / 3.0d0
     D(3,1)= D(1,3)
     D(3,2)= D(2,3)
     D(3,3)= kk + 2.0d0 * stiff_g / 3.0d0
     D(4,4)= stiff_g
     D(5,5)= stiff_g
     D(6,6)= stiff_g

     !write(201,'(A, 100E10.3)') "deriv_e, tau_max, ita_e, xi_e, stiff_g: ", deriv_e, tau_max, ita_e, xi_e, stiff_g
     
   end subroutine

!####################################################################   
   subroutine RO_EquivalentStress(gauss, stress_bak, stress, D, iter, strain_bak, strain)
!####################################################################
     use hecmw
     use mMaterial
     implicit none
     
     type( tGaussStatus ), intent(IN)       :: gauss          !> status of qudrature point
     REAL(KIND=kreal), INTENT(INOUT)           :: stress_bak(6) !< stress
     REAL(KIND=kreal), INTENT(INOUT)           :: stress(6) !< stress
     REAL(KIND=kreal), INTENT(INOUT)           :: D(:,:)    !< strain-stress relation
     integer:: iter
     REAL(KIND=kreal), INTENT(IN), optional :: strain_bak(6) !< stress
     REAL(KIND=kreal), INTENT(IN), optional :: strain(6) !< stress
     
     integer:: flag, mtype
     TYPE( tMaterial ):: matl      !< material properties
     integer :: i
     real(kind=kreal) :: pi
     real(kind=kreal) :: EE, PP, kk, C, theta, alpha, beta, g_max, rho
     real(kind=kreal) :: sigma_m, deriv(6), deriv_e, tau_max
     real(kind=kreal) :: ita_e, xi_e, stiff_g
     real(kind=kreal) :: sigma_m_new, deriv_new(6), deriv_e_new, tau_max_new
     real(kind=kreal) :: ita_e_new, xi_e_new, stiff_g_new, stress_new(6),stress_new_h(6)
     real(kind=kreal) :: deriv_org(6), deriv_e_org, ita_e_org
     real(kind=kreal) :: d_deriv(6), deriv_r(6), dist_param_t, sigma_m_stress

     REAL(KIND = kREAL) :: C1_0, C2_0, C1_inf, C2_inf
     REAL(KIND = kREAL) :: strain_m, gamma_m, strain_new_h(6)
     REAL(KIND = kREAL) :: gamma_r_sklt, strain_m_new, gamma_m_new, target_damp, strain_new(6)
     REAL(KIND = kREAL) :: deriv_src(1), damp_tar(1)
     REAL(KIND = kREAL) :: gamma_r_sklt_new
     logical :: ierr
! ----------------------------------------------------------------------------
     
     INTEGER(KIND = kint), POINTER :: stress_zero_index, strain_zero_index
     REAL(KIND = kREAL), POINTER :: stress_zero(:)
     REAL(KIND = kREAL), POINTER :: strain_zero(:)
     !REAL(KIND = kREAL), POINTER :: strain_max
     
! ----------------------------------------------------------------------------
     matl = gauss%pMaterial
     mtype = matl%mtype
     
     !  pointer設定
     if(mtype == RO) then
       stress_zero_index => gauss%istatus(1)
       stress_zero       => gauss%fstatus(1:350)
     elseif(mtype == GHE) then
       stress_zero_index => gauss%istatus(1)
       strain_zero_index => gauss%istatus(1)
       stress_zero       => gauss%fstatus(1:350)
       strain_zero       => gauss%fstatus(351:700)
       !strain_max        => gauss%fstatus(701)
     endif

! ----------------------------------------------------------------------------
     
     pi = 4.0D0*DATAN( 1.0D0 )
     
     if(iter == 1) then
         write(401,'(I10, 100E12.5)') stress_zero_index, stress_zero(s0_index(stress_zero_index, 6)), &
                                      strain_bak(6), stress_bak(6), stress_bak(1:3)
     endif
     
     if(mtype == RO) then
       call Set_param_RO_material(matl, ee, pp, kk, rho, C, theta, alpha, beta)
     elseif(mtype == GHE) then
       call Set_param_RO_material(matl, ee, pp, kk, rho, C, theta, alpha, beta, C1_0, C2_0, C1_inf, C2_inf)
     endif
     g_max = ee / 2.0d0 / (1.0d0 + pp)
     
     ! 前回ステップにおける値
     if(mtype == GHE) then
       call calc_deriv(strain_bak, strain_zero(s0_index(strain_zero_index, 1):s0_index(strain_zero_index, 6)), &
                       deriv, gamma_m, strain_m)
     endif
     
     call calc_deriv(stress_bak, stress_zero(s0_index(stress_zero_index, 1):s0_index(stress_zero_index, 6)), &
                     deriv, deriv_e, sigma_m)

     tau_max = C * dcos(theta) + sigma_m * dsin(theta)
     ita_e = deriv_e / tau_max
     gamma_r_sklt = tau_max / g_max

     stress_new(:) = stress_bak(:) + stress(:)
     strain_new(:) = strain_bak(:) + strain(:)
     if(mtype == GHE) then
       call calc_deriv(strain_new, strain_zero(s0_index(strain_zero_index, 1):s0_index(strain_zero_index, 6)), &
                       deriv_new, gamma_m_new, strain_m_new)
     endif
     call calc_deriv(stress_new, stress_zero(s0_index(stress_zero_index, 1):s0_index(stress_zero_index, 6)), &
                     deriv_new, deriv_e_new, sigma_m_new)
     
     tau_max_new = C * dcos(theta) + sigma_m_new * dsin(theta)
     ita_e_new = deriv_e_new / tau_max_new
     gamma_r_sklt_new = tau_max_new  / g_max
     
     ! 更新前の剛性と更新後の剛性の間を取って接線剛性を求める
     stress_new_h(:) = stress_bak(:) + stress(:)*0.5d0
     strain_new_h(:) = strain_bak(:) + strain(:)*0.5d0
     
     sigma_m_stress = (1.0d0 / 3.0d0 ) * (stress(1) + stress(2) + stress(3))
     ! 載荷除荷の判定
     d_deriv(1) =                (stress(1) - sigma_m_stress)
     d_deriv(2) =                (stress(2) - sigma_m_stress)
     d_deriv(3) =                (stress(3) - sigma_m_stress)
     d_deriv(4) = dsqrt(2.0d0) * (stress(4)                 )
     d_deriv(5) = dsqrt(2.0d0) * (stress(5)                 )
     d_deriv(6) = dsqrt(2.0d0) * (stress(6)                 )
     deriv_r(1) =                deriv(1)
     deriv_r(2) =                deriv(2)
     deriv_r(3) =                deriv(3)
     deriv_r(4) = dsqrt(2.0d0) * deriv(4)
     deriv_r(5) = dsqrt(2.0d0) * deriv(5)
     deriv_r(6) = dsqrt(2.0d0) * deriv(6)
     
     dist_param_t = -1.0d0 * dot_product(d_deriv, deriv_r) / dot_product(d_deriv, d_deriv)
     
     ! 除荷の場合
     !if( ita_e > ita_e_new) then
     if( ((ita_e > ita_e_new) .or. (0.0d0 < dist_param_t .and. dist_param_t < 1.0d0) ) .and. iter==1 ) then
       stress_zero_index = stress_zero_index + 1
       stress_zero(s0_index(stress_zero_index, 1):s0_index(stress_zero_index, 6)) = stress_bak(:)
       stress_zero(s0_index(stress_zero_index, 7)) = ita_e
       
       if(mtype == GHE) then
         
         strain_zero(s0_index(stress_zero_index, 1):s0_index(stress_zero_index, 6)) = strain_bak(:)
         if(stress_zero_index == 2) then 
           ! 等価基準ひずみを計算する処理
           deriv_src(1) = gamma_m
           CALL fetch_TableData(MC_YIELD, matl%dict, damp_tar, ierr, deriv_src)
           !write(*, *) deriv_src(1), damp_tar(1)
           target_damp = damp_tar(1)
           strain_zero(s0_index(stress_zero_index, 7)) = &
                       get_equiv_gamma_r(target_damp, 1.0d-7, 1.0d2, &
                                         alpha, beta, C1_0, C2_0, C1_inf, C2_inf, g_max, gamma_r_sklt, &
                                         strain_zero(s0_index(1, 7)))
         else
           strain_zero(s0_index(stress_zero_index, 7)) = strain_zero(s0_index(stress_zero_index - 1, 7))
         endif
       endif
       
       ! 直近の除荷点が変わったのでita_e_newを更新する
       call calc_deriv(stress_new, stress_zero(s0_index(stress_zero_index, 1):s0_index(stress_zero_index, 6)), &
                     deriv_new, deriv_e_new, sigma_m_new)

       tau_max_new = C * dcos(theta) + sigma_m_new * dsin(theta)
       ita_e_new = deriv_e_new / tau_max_new
     endif
     
     ! index=1の場合は何もしない         
     if(stress_zero_index > 2)then
        ! 一気に大きくなった場合のため
        do while((stress_zero(s0_index(stress_zero_index, 7)) < ita_e_new) .and. (stress_zero_index > 3))
            stress_zero_index = stress_zero_index - 2
        end do
     endif
     
     !骨格曲線に関する降伏評価
     deriv_org(1) = stress_new(1) - sigma_m_new
     deriv_org(2) = stress_new(2) - sigma_m_new
     deriv_org(3) = stress_new(3) - sigma_m_new
     deriv_org(4) = stress_new(4)              
     deriv_org(5) = stress_new(5)              
     deriv_org(6) = stress_new(6)              
     deriv_e_org = dsqrt((1.0d0 / 2.0d0) * (deriv_org(1)*deriv_org(1) + deriv_org(2)*deriv_org(2) + deriv_org(3)*deriv_org(3) + &
                                   2.0d0 * (deriv_org(4)*deriv_org(4) + deriv_org(5)*deriv_org(5) + deriv_org(6)*deriv_org(6))))
     ita_e_org = deriv_e_org / tau_max_new
     if((stress_zero_index > 1) .and. (ita_e_org > stress_zero(s0_index(2, 7))) .and. iter==1 ) stress_zero_index = 1
     !write(801, '(2I10, 100E12.5)') iter, stress_zero_index, ita_e_org, gauss%stress_zero(2, 7)
     
     call ROMatrix( gauss, stress_new_h, gauss%istatus(1), gauss%fstatus, D, strain_new_h)
     
     if(mtype == GHE .and. stress_zero_index == 1) then
       strain_zero(s0_index(1, 7)) = gamma_m_new
     endif
     
   end subroutine


!####################################################################   
   function s0_index(i_point, j_component)
!####################################################################
   
     ! USE
                  
 !--------------------------------------------------------------------

     INTEGER(KIND = kint), INTENT(IN) :: i_point, j_component
     
     INTEGER(KIND = kint), PARAMETER :: num_comp = 7
     INTEGER(KIND = kint) :: s0_index
     
     !write(999, "(100I10)") i_point, j_component
     s0_index = (i_point - 1 ) * num_comp + j_component
   
    end function

    
!#################################################################### 
    subroutine Set_param_RO_material                               &
              (matl, ee, pp, kk, rho, C, theta, alpha, beta,        &
              C1_0, C2_0, C1_inf, C2_inf)
!#################################################################### 
    
      ! USE
                  
 !--------------------------------------------------------------------
       
      TYPE( tMaterial ), INTENT(IN) :: matl
      REAL(KIND = kREAL), INTENT(OUT) :: ee, pp, kk, rho, C, theta, alpha, beta
      REAL(KIND = kREAL), optional:: C1_0, C2_0, C1_inf, C2_inf
      
      REAL(KIND = kREAL) :: pi
      
      pi = 4.0D0*DATAN( 1.0D0 )
              
      ee = matl%variables(M_YOUNGS)
      pp = matl%variables(M_POISSON)
      kk = ee / 3.0d0 / (1.0d0 - 2.0d0 * pp)
      rho = matl%variables(M_DENSITY)
      C  = matl%variables(M_PLCONST1)                                  !ROパラメータ第一引数
      theta = 2.0d0 * pi * matl%variables(M_PLCONST2) / 360.0d0        !ROパラメータ第二引数
      alpha = matl%variables(M_PLCONST3)                               !ROパラメータ第三引数
      beta = matl%variables(M_PLCONST4)                                !ROパラメータ第四引数
      if(present(C1_0)) C1_0 = matl%variables(M_PLCONST5)              !ROパラメータ第五引数
      if(present(C2_0)) C2_0 = matl%variables(M_KINEHARD)              !ROパラメータ第六引数
      if(present(C1_inf)) C1_inf = matl%variables(M_GOODMAN_K1)        !ROパラメータ第七引数
      if(present(C2_inf)) C2_inf = matl%variables(M_GOODMAN_K2)        !ROパラメータ第八引数

!#################################################################### 
    end subroutine
!#################################################################### 

!####################################################################   
   function stress_GHE(gamma, alpha, beta, C1_0, C2_0, C1_inf, C2_inf, stiff0, &
                       gamma_r, C1_gamma, C2_gamma)
!####################################################################
   
     ! USE
                  
 !--------------------------------------------------------------------

     REAL(KIND = kREAL), INTENT(IN) :: alpha, beta, C1_0, C2_0, C1_inf, C2_inf, stiff0
     REAL(KIND = kREAL), INTENT(IN) :: gamma_r, gamma
     REAL(KIND = kREAL), INTENT(OUT), optional :: C1_gamma, C2_gamma
     REAL(KIND = kREAL) :: stress_GHE
     
     REAL(KIND = kREAL) :: pi
     REAL(KIND = kREAL) :: x_val, C1_p, C1_m, C2_p, C2_m, coeff_1, coeff_2
     
     pi = 4.0D0*DATAN( 1.0D0 )
     
     x_val = gamma / gamma_r
     C1_p = (C1_0 + C1_inf) / 2
     C1_m = (C1_0 - C1_inf) / 2
     C2_p = (C2_0 + C2_inf) / 2
     C2_m = (C2_0 - C2_inf) / 2
     
     coeff_1 = dcos(pi / (alpha / x_val + 1))
     C1_gamma = C1_p + C1_m *coeff_1
     coeff_2 = dcos(pi / (beta  / x_val + 1))
     C2_gamma = C2_p + C2_m *coeff_2
     
     stress_GHE = gamma * stiff0 / ((1.0d0 / C1_gamma) + (1.0d0 / C2_gamma) * x_val)
   
    end function
    
!####################################################################   
   function stress_GHE_deriv(gamma, alpha, beta, C1_0, C2_0, C1_inf, C2_inf, stiff0, &
                             gamma_r)
!####################################################################
   
     ! USE
                  
 !--------------------------------------------------------------------

     REAL(KIND = kREAL), INTENT(IN) :: alpha, beta, C1_0, C2_0, C1_inf, C2_inf, stiff0
     REAL(KIND = kREAL), INTENT(IN) :: gamma_r, gamma
     REAL(KIND = kREAL) :: stress_GHE_deriv
     REAL(KIND = kREAL) :: pi
     REAL(KIND = kREAL) :: x_val, C1_p, C1_m, C2_p, C2_m, coeff_1, C1_gamma, coeff_2, C2_gamma
     REAL(KIND = kREAL) :: coeff_r1, coeff_r2, coeff_r1s, coeff_r2s
     REAL(KIND = kREAL) :: C1_g_r, C2_g_r, C1_g_i, C2_g_i, temp_up, temp_down

     pi = 4.0D0*DATAN( 1.0D0 )
     
     x_val = gamma / gamma_r
     C1_p = (C1_0 + C1_inf) / 2
     C1_m = (C1_0 - C1_inf) / 2
     C2_p = (C2_0 + C2_inf) / 2
     C2_m = (C2_0 - C2_inf) / 2
     
     coeff_1 = dcos(pi / (alpha / x_val + 1))
     if(x_val == 0.0d0) coeff_1 = 1.0d0
     C1_gamma = C1_p + C1_m *coeff_1
     coeff_2 = dcos(pi / (beta  / x_val + 1))
     if(x_val == 0.0d0) coeff_2 = 1.0d0
     C2_gamma = C2_p + C2_m *coeff_2
     
     coeff_r1 = dsin(pi / (alpha / x_val + 1))
     if(x_val == 0.0d0) coeff_r1 = 0.0d0
     coeff_r2 = dsin(pi / (beta  / x_val + 1))
     if(x_val == 0.0d0) coeff_r2 = 0.0d0
     coeff_r1s= pi * alpha / gamma_r / x_val / x_val / (alpha / x_val + 1.0d0) / (alpha / x_val + 1.0d0)
     if(x_val == 0.0d0) coeff_r1s = 0.0d0
     coeff_r2s= pi * beta  / gamma_r / x_val / x_val / (beta  / x_val + 1.0d0) / (beta  / x_val + 1.0d0)
     if(x_val == 0.0d0) coeff_r2s = 0.0d0
     
     C1_g_r = -1.0d0 * C1_m * coeff_r1 * coeff_r1s
     C2_g_r = -1.0d0 * C2_m * coeff_r2 * coeff_r2s
     
     C1_g_i = 1.0d0 / C1_gamma
     C2_g_i = 1.0d0 / C2_gamma
     
     temp_up = (C1_g_i + gamma * (C1_g_r * C1_g_i * C1_g_i + C2_g_r * C2_g_i * C2_g_i * x_val))
     temp_down = (C1_g_i + C2_g_i * x_val)
     !write(800, "(100E15.5)") gamma, stiff0, gamma_r, x_val, C1_gamma, C2_gamma, C1_g_r, C2_g_r, temp_up, temp_down
     stress_GHE_deriv = stiff0 * temp_up / (temp_down * temp_down)
    end function
    
    
    function equiv_damp_GHE(x_min, x_max, &
                            alpha, beta, C1_0, C2_0, C1_inf, C2_inf, stiff0, &
                            gamma_r, gamma_r_sklt)

      real(kind=kreal), intent(in) :: x_min, x_max
      real(kind=kreal), intent(in) :: alpha, beta, C1_0, C2_0, C1_inf, C2_inf, stiff0, gamma_r, gamma_r_sklt
      real(kind=kreal) :: equiv_damp_GHE
      
      INTEGER(KIND = kint), PARAMETER :: num_div = 1000
      INTEGER(KIND = kint) :: i
      real(kind=kreal) :: pi, x1, step, x2
      real(kind=kreal) :: stress1, stress2
      real(kind=kreal) :: area_GHE, area_tri, area_dw
      real(kind=kreal) :: C1_gamma, C2_gamma, c1_g_i, c2_g_i, stiff, stress_sklt
      
      pi = 4.0D0*DATAN( 1.0D0 )
      
      !write(802, "(100E12.3)") gamma_r_sklt, alpha, beta, C1_0, C2_0, C1_inf, C2_inf, stiff0
      stress_sklt = stress_GHE(x_max, alpha, beta, C1_0, C2_0, C1_inf, C2_inf, stiff0, &
                               gamma_r_sklt, C1_gamma, C2_gamma)
      c1_g_i = 1.0d0 / C1_gamma
      c2_g_i = 1.0d0 / C2_gamma
      stiff = (stress_sklt / x_max) * (c1_g_i + c2_g_i * (x_max / gamma_r))
      
      step = (x_max - x_min) / num_div
      area_GHE = 0.0d0
      x1 = x_min
      do i = 1, num_div
        x2 = x_min + step * i
        stress1 = (stiff * x1) / (c1_g_i + c2_g_i * (x1 / gamma_r))
        stress2 = (stiff * x2) / (c1_g_i + c2_g_i * (x2 / gamma_r))
        area_GHE = area_GHE + (stress1 + stress2) * step / 2.0d0
        !write(802, "(100E12.3)") x2, stress2
        x1 = x2
      enddo
      area_tri = (x_max - x_min) * stress2 / 2.0d0
      !write(804, "(100E12.3)")area_GHE, area_tri
      area_dw = 8.0d0 * (area_GHE - area_tri)
      equiv_damp_GHE = (area_dw / area_tri) / (4.0d0 * pi)
      
    end function
    
    function get_equiv_gamma_r(target_damp, g_min, g_max, &
                               alpha, beta, C1_0, C2_0, C1_inf, C2_inf, stiff0, gamma_r_sklt, &
                               gamma_max)
      real(kind=kreal), intent(in) :: g_min, g_max, target_damp
      real(kind=kreal), intent(in) :: alpha, beta, C1_0, C2_0, C1_inf, C2_inf, stiff0, gamma_max, gamma_r_sklt
      real(kind=kreal) :: get_equiv_gamma_r
      
      real(kind=kreal) :: res, x1, x2, damp1, damp2, flag, x_r, damp_r
      real(kind=kreal), PARAMETER :: error = 1.0e-5
      INTEGER(KIND = kint) :: i_count
      
      res = 1.0d0
      x1 = g_min
      x2 = g_max
      !write(801, "(100E12.3)") alpha, beta, C1_0, C2_0, C1_inf, C2_inf, stiff0, x1, gamma_r_sklt
      damp1 = equiv_damp_GHE(0.0d0, gamma_max, &
                             alpha, beta, C1_0, C2_0, C1_inf, C2_inf, stiff0, &
                             x1, gamma_r_sklt)
      damp2 = equiv_damp_GHE(0.0d0, gamma_max, &
                             alpha, beta, C1_0, C2_0, C1_inf, C2_inf, stiff0, &
                             x2, gamma_r_sklt)
      !write(801, "(100E12.3)") x1, damp1
      !write(801, "(100E12.3)") x2, damp2
      flag = dsign(1.0d0, damp2 - damp1)
      i_count = 0
      do while ((res > error) .and. (i_count < 10000)) 
        x_r = (x1 + x2) / 2.0d0
        damp_r =equiv_damp_GHE(0.0d0, gamma_max, &
                             alpha, beta, C1_0, C2_0, C1_inf, C2_inf, stiff0, &
                             x_r, gamma_r_sklt)
        res = dabs(damp_r - target_damp) / target_damp
        !write(801, "(100E12.3)") x_r, damp_r
        if((damp_r - target_damp) * flag > 0.0d0) then 
          x2 = x_r
          damp2 = damp_r
        else
          x1 = x_r
          damp1 = damp_r
        endif
        i_count = i_count + 1
        
      enddo
      get_equiv_gamma_r = x_r
      !write(801, "(100E12.3)") 0.0d0
    end function
    
    subroutine calc_deriv(stress, stress_zero, deriv, deriv_e, sigma_m)
      REAL(KIND = kREAL), INTENT(IN) :: stress(6)
      REAL(KIND = kREAL), INTENT(IN) :: stress_zero(6)
      REAL(KIND = kREAL), INTENT(OUT) :: deriv(6)
      REAL(KIND = kREAL), INTENT(OUT) :: deriv_e
      REAL(KIND = kREAL), INTENT(OUT) :: sigma_m
      
      sigma_m = ( 1.0d0 / 3.0d0 ) * (stress(1) + stress(2) + stress(3))
      deriv(1) = stress(1) - sigma_m - stress_zero(1)
      deriv(2) = stress(2) - sigma_m - stress_zero(2)
      deriv(3) = stress(3) - sigma_m - stress_zero(3)
      deriv(4) = stress(4)           - stress_zero(4)
      deriv(5) = stress(5)           - stress_zero(5)
      deriv(6) = stress(6)           - stress_zero(6)
      deriv_e = dsqrt((1.0d0 / 2.0d0) * (deriv(1)*deriv(1) + deriv(2)*deriv(2) + deriv(3)*deriv(3) + &
                               2.0d0  * (deriv(4)*deriv(4) + deriv(5)*deriv(5) + deriv(6)*deriv(6))))
    
    end subroutine
    
end module m_RO

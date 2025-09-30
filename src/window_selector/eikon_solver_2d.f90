subroutine FSM_UW_PS_sphe_2d(xx,yy,nx,ny,spha,sphb,sphc,T,fun,x0,y0)
    integer :: nx,ny
    double precision :: xx(nx),yy(ny),spha(nx,ny),sphb(nx,ny),sphc(nx,ny),T(nx,ny),fun(nx,ny),x0,y0,ischange(nx,ny)

    double precision :: a(nx,ny),b(nx,ny),c(nx,ny)
    double precision :: dx,dy,T0(nx,ny),T0x(nx,ny),T0y(nx,ny),a0,b0,c0,fun0
    double precision :: tau(nx,ny)
    integer :: iix,iiy,iter,i_sweep, x_id1,x_id2,y_id1,y_id2,xdirec,ydirec, idx0, idy0
    double precision :: L1_dif,Linf_dif,L1_err,Linf_err,tau_old(nx,ny), r1,r2
    integer,parameter :: MaxIter = 2000
    double precision,parameter :: tol=(10.0)**(-4),eps=(10.0)**(-14)

    integer :: i_case,i_cand,i_solution,count_cand
    double precision :: px,qx,py,qy
    double precision :: eqn_a,eqn_b,eqn_c,Delta,tmp_tau,dis
    double precision :: characteristic_x, characteristic_y, Tx, Ty
    logical :: is_causility
    double precision :: canditate(50)


    ! ------------------------ 构造网格 ------------------------ 
    dx=xx(2)-xx(1); dy=yy(2)-yy(1);

    do iix=1,nx
        do iiy=1,ny
            a(iix,iiy) = spha(iix,iiy)
            b(iix,iiy) = sphb(iix,iiy)/(xx(iix)**2)
            c(iix,iiy) = sphc(iix,iiy)/(xx(iix))
        end do
    end do

    ! ------------------------ 构造 T0 ------------------------ 

    ! 震源处参数离散化
    idx0=floor((x0-xx(1))/dx+1); idy0=floor((y0-yy(1))/dy+1); 
    r1 = min(1.0d0,(x0-xx(idx0))/dx); r2 = min(1.0d0,(y0-yy(idy0))/dy); 

    a0=(1-r1)*(1-r2)*a(idx0,idy0)+(1-r1)*r2*a(idx0,idy0+1) &
    & +r1*(1-r2)*a(idx0+1,idy0)+r1*r2*a(idx0+1,idy0+1) 

    b0=(1-r1)*(1-r2)*b(idx0,idy0)+(1-r1)*r2*b(idx0,idy0+1) &
    & +r1*(1-r2)*b(idx0+1,idy0)+r1*r2*b(idx0+1,idy0+1) 

    c0=(1-r1)*(1-r2)*c(idx0,idy0)+(1-r1)*r2*c(idx0,idy0+1) &
    & +r1*(1-r2)*c(idx0+1,idy0)+r1*r2*c(idx0+1,idy0+1) 

    fun0=(1-r1)*(1-r2)*fun(idx0,idy0)+(1-r1)*r2*fun(idx0,idy0+1) &
    & +r1*(1-r2)*fun(idx0+1,idy0)+r1*r2*fun(idx0+1,idy0+1) 


    ! solve T0 = sqrt((b0(x-x0)^2+a0(y-y0)^2+2c0(x-x0)(y-y0))/(a0b0-c0^2))
    do iix=1,nx
        do iiy=1,ny
            T0(iix,iiy) = fun0*sqrt((b0/(a0*b0-c0**2)*(xx(iix)-x0)**2 + a0/(a0*b0-c0**2)*(yy(iiy)-y0)**2 &
                        & + 2*c0/(a0*b0-c0**2)*(xx(iix)-x0)*(yy(iiy)-y0)))
            if (T0(iix,iiy) .eq. 0) then
                T0x(iix,iiy) = 0
                T0y(iix,iiy) = 0
            else
                T0x(iix,iiy) = fun0**2*(b0/(a0*b0-c0**2)*(xx(iix)-x0)+c0/(a0*b0-c0**2)*(yy(iiy)-y0))/T0(iix,iiy)
                T0y(iix,iiy) = fun0**2*(a0/(a0*b0-c0**2)*(yy(iiy)-y0)+c0/(a0*b0-c0**2)*(xx(iix)-x0))/T0(iix,iiy)
            end if

            if ( abs((xx(iix)-x0)/dx)<=1 .and. abs((yy(iiy)-y0)/dy)<=1) then
                tau(iix,iiy) = 1  !震源周围几个点，直接认为是常速度结构，给出解析解，即和T0相等
                ischange(iix,iiy)=0
                if (iix==1 .or. iix==nx .or. iiy==1 .or. iiy==ny) then
                    write(*,*) 'source on the boundary, mesh error'
                    stop
                end if
            else
                tau(iix,iiy) = 10
                ischange(iix,iiy)=1
            end if


        end do
    end do

    ! step 2, solve Tau, H(tau) = a tau_x^2+ b tau_y^2 + (2aTx-2cTy) tau_x + (2bTy-2cTx) tau_y
    !                           -2c tau_x tau_y + (aTx^2+bTy^2-2cTxTy) = f^2

    do iter =1,MaxIter
    ! do iter =1,1
        L1_dif=10000; Linf_dif=10000;L1_err=0;Linf_err=0;
        tau_old = tau
        do i_sweep = 1,4
        ! do i_sweep = 1,1
            select case(i_sweep)
                case (1)
                    x_id1 =  1; x_id2 = Nx; xdirec =  1
                    y_id1 =  1; y_id2 = Ny; ydirec =  1
                case (2)
                    x_id1 =  1; x_id2 = Nx; xdirec =  1
                    y_id1 = Ny; y_id2 =  1; ydirec = -1
                case (3)
                    x_id1 = Nx; x_id2 =  1; xdirec = -1
                    y_id1 =  1; y_id2 = Ny; ydirec =  1
                case (4)
                    x_id1 = Nx; x_id2 =  1; xdirec = -1
                    y_id1 = Ny; y_id2 =  1; ydirec = -1    
            end select


            ! iter 1 x: 1 -> Nx, y: 1 -> Ny
            ! iter 2 x: 1 -> Nx, y: Ny -> 1 
            ! iter 3 x: Nx -> 1, y: 1 -> Ny 
            ! iter 4 x: Nx -> 1, y: Ny -> 1 
            do iix=x_id1,x_id2,xdirec
                do iiy=y_id1,y_id2,ydirec
                    if(ischange(iix,iiy)==1) then
    !--------------------------------------- calculate stencil start
                        count_cand = 1
                        do i_case=1,4       ! 4个扇形方向
                            select case (i_case)
                                ! (T0*tau)_x = px*tau(iix,iiy)+qx; (T0*tau)_y = py*tau(iix,iiy)+qy
                            case (1)    ! 左下  -1, -1
                                if (iix == 1 .or. iiy == 1) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)+T0(iix,iiy)/dx); qx = -T0(iix,iiy)/dx*tau(iix-1,iiy)     
                                py = (T0y(iix,iiy)+T0(iix,iiy)/dy); qy = -T0(iix,iiy)/dy*tau(iix,iiy-1)     
                            case (2)    ! 左上  -1, +1
                                if (iix == 1 .or. iiy == Ny) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)+T0(iix,iiy)/dx); qx = -T0(iix,iiy)/dx*tau(iix-1,iiy)     
                                py = (T0y(iix,iiy)-T0(iix,iiy)/dy); qy =  T0(iix,iiy)/dy*tau(iix,iiy+1)   
                            case (3)    ! 右下  +1, -1
                                if (iix == Nx .or. iiy == 1) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)-T0(iix,iiy)/dx); qx = +T0(iix,iiy)/dx*tau(iix+1,iiy)     
                                py = (T0y(iix,iiy)+T0(iix,iiy)/dy); qy = -T0(iix,iiy)/dy*tau(iix,iiy-1)
                            case (4)    ! 右上  +1, +1
                                if (iix == Nx .or. iiy == Ny) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)-T0(iix,iiy)/dx); qx =  T0(iix,iiy)/dx*tau(iix+1,iiy)     
                                py = (T0y(iix,iiy)-T0(iix,iiy)/dy); qy =  T0(iix,iiy)/dy*tau(iix,iiy+1) 
                            end select
                            
                            

                            ! 方程是关于tau的一元二次方程 a*(px*tau+qx)^2 + b*(py*tau+qy)^2 -2c*(px*tau+qx)*(py*tau+qy) = s^2
                            ! 判别式
                            eqn_a = a(iix,iiy)*px**2 + b(iix,iiy)*py**2 -2*c(iix,iiy)*px*py
                            eqn_b = 2*a(iix,iiy)*px*qx + 2*b(iix,iiy)*py*qy - 2*c(iix,iiy)*(px*qy+py*qx)
                            eqn_c = a(iix,iiy)*qx**2 + b(iix,iiy)*qy**2 -2*c(iix,iiy)*qx*qy-fun(iix,iiy)**2
                            Delta = eqn_b**2-4*eqn_a*eqn_c
                            

                            if (Delta>=0) then  !   one or two solutions
                                ! check the causality condition: the characteristic passing through (iix,iiy) is in between used two sides
                                do i_solution=1,2      ! 2个解，相同解也看做两个解
                                    select case (i_solution)
                                    case(1)
                                        tmp_tau = (-eqn_b + sqrt(Delta))/(2*eqn_a)
                                    case(2)
                                        tmp_tau = (-eqn_b - sqrt(Delta))/(2*eqn_a)
                                    end select

                                    ! characteristic is
                                    Tx = px*tmp_tau+qx
                                    Ty = py*tmp_tau+qy
                                    characteristic_x = a(iix,iiy)*Tx - c(iix,iiy)*Ty 
                                    characteristic_y = b(iix,iiy)*Ty - c(iix,iiy)*Tx 

                                    is_causility = .false.
                                    ! check the causality condition: (特殊情况：对于直角，直接该向量在对应象限即可)
                                    select case (i_case)
                                    case (1)    ! 从左下过来，向量位于第一象限
                                        if (characteristic_x >= 0 .and. characteristic_y >= 0) then
                                            is_causility = .true.
                                        end if
                                    case (2)    ! 从左上过来，向量位于第四象限
                                        if (characteristic_x >= 0 .and. characteristic_y <= 0) then
                                            is_causility = .true.
                                        end if
                                    case (3)    ! 从右下过来，向量位于第二象限
                                        if (characteristic_x <= 0 .and. characteristic_y >= 0) then
                                            is_causility = .true.
                                        end if
                                    case (4)    ! 从右上过来，向量位于第三象限
                                        if (characteristic_x <= 0 .and. characteristic_y <= 0) then
                                            is_causility = .true.
                                        end if
                                    end select

                                    ! if satisfying the causility condition, retain it as a canditate solution
                                    if(is_causility) then
                                        canditate(count_cand) = tmp_tau
                                        count_cand = count_cand + 1
                                    end if


                                end do
                            end if

                        end do

                        do i_case=1,4   ! 4条特征线方向
                            select case (i_case)
                                ! (T0*tau)_x = px*tau(iix,iiy)+qx; 
                                ! (T0*tau)_y = py*tau(iix,iiy)+qy
                            case (1)    ! 从左边过来, 将 characteristic_y = 0 带入方程，得到 (a-c^2/b)T_x^2 =  , that is, px*tau+qx = sqrt(s^2/(a-c^2/b))
                                if (iix == 1) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)+T0(iix,iiy)/dx); qx = -T0(iix,iiy)/dx*tau(iix-1,iiy)    
                                dis = sqrt(fun(iix,iiy)**2*(b(iix,iiy)/(a(iix,iiy)*b(iix,iiy)-c(iix,iiy)**2)))
                                tmp_tau = (dis-qx)/px
                                if (tmp_tau*T0(iix,iiy) >= tau(iix-1,iiy)*T0(iix-1,iiy)) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                                tmp_tau = (-dis-qx)/px
                                if (tmp_tau*T0(iix,iiy) >= tau(iix-1,iiy)*T0(iix-1,iiy)) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                            case (2)    ! 从右边过来,
                                if (iix == Nx) then
                                    cycle
                                end if
                                px = (T0x(iix,iiy)-T0(iix,iiy)/dx); qx =  T0(iix,iiy)/dx*tau(iix+1,iiy)    
                                dis = sqrt(fun(iix,iiy)**2*(b(iix,iiy)/(a(iix,iiy)*b(iix,iiy)-c(iix,iiy)**2)))
                                tmp_tau = (dis-qx)/px 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix+1,iiy)*T0(iix+1,iiy)) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                                tmp_tau = (-dis-qx)/px 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix+1,iiy)*T0(iix+1,iiy)) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if

                            case (3)    ! 从下边过来，将 characteristic_x = 0 带入方程，得到 (a-c^2/b)T_x^2 =  , that is, px*tau+qx = sqrt(s^2/(a-c^2/b))
                                if (iiy == 1) then
                                    cycle
                                end if
                                py = (T0y(iix,iiy)+T0(iix,iiy)/dy); qy = -T0(iix,iiy)/dy*tau(iix,iiy-1)
                                dis = sqrt(fun(iix,iiy)**2*(a(iix,iiy)/(a(iix,iiy)*b(iix,iiy)-c(iix,iiy)**2)))
                                tmp_tau = (dis-qy)/py 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix,iiy-1)*T0(iix,iiy-1)) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                                tmp_tau = (-dis-qy)/py 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix,iiy-1)*T0(iix,iiy-1)) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                            case (4)    ! 从上边过来
                                if (iiy == Ny) then
                                    cycle
                                end if
                                py = (T0y(iix,iiy)-T0(iix,iiy)/dy); qy =  T0(iix,iiy)/dy*tau(iix,iiy+1) 
                                dis = sqrt(fun(iix,iiy)**2*(a(iix,iiy)/(a(iix,iiy)*b(iix,iiy)-c(iix,iiy)**2)))
                                tmp_tau = (dis-qy)/py 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix,iiy+1)*T0(iix,iiy+1)) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                                tmp_tau = (-dis-qy)/py 
                                if (tmp_tau*T0(iix,iiy) >= tau(iix,iiy+1)*T0(iix,iiy+1)) then
                                    canditate(count_cand) = tmp_tau
                                    count_cand = count_cand + 1
                                end if
                            end select

                        end do

                        ! 选择最小的作为解
                        do i_cand=1,count_cand-1
                            tau(iix,iiy) = min(tau(iix,iiy),canditate(i_cand))
                        end do
    !--------------------------------------- calculate stencil end
                    end if
                end do
            end do
        end do

        L1_dif=0; Linf_dif=0
        do iix=1,nx
            do iiy=1,ny
                L1_dif=L1_dif+abs(tau(iix,iiy)-tau_old(iix,iiy))*T0(iix,iiy)
                Linf_dif=max(Linf_dif,abs(tau(iix,iiy)-tau_old(iix,iiy))*T0(iix,iiy))               
            end do
        end do

        if (abs(L1_dif)<tol .and. Linf_dif<tol) then
            write(*,'(a,f15.7,a,f15.7)') 'L_1(Tnew-Told)=',L1_dif,'  L_inf(Tnew-Told)',Linf_dif
            write(*,*) 'iter ',iter,', T is steady'
            exit
        else
        !    write(*,*) 'iter ',iter,', T is changing, continue ... '
        end if

        if (iter==MaxIter) then    
        !    write(*,*) 'iter ',iter,', max iteration steps'
        end if
        write(*,'(a,f15.7,a,f15.7)') 'L_1(Tnew-Told)=',L1_dif,'  L_inf(Tnew-Told)',Linf_dif
        ! write(*,'(a,f15.8,a,f11.8)') 'L_1(T_LF-u)=',L1_err,'  L_inf(T_LF-u)',Linf_err
    end do

    T=tau*T0

end subroutine

subroutine interp2d(xx,yy,nx,ny,val,x0,y0,v0)
    integer :: nx,ny,idx0,idy0
    double precision :: xx(nx),yy(ny),val(nx,ny),x0,y0,v0,dx,dy,r1,r2

    dx=xx(2)-xx(1); dy=yy(2)-yy(1)

    if ( x0<=xx(1) .or. x0>=xx(nx) .or. y0<=yy(1) .or. y0>=yy(ny)) then
        print *, 'point (',x0,',',y0,')'
        print *, 'mesh range x:[',xx(1),',',xx(nx),'], y:[',yy(1),',',yy(ny),']'
        print *, 'out of mesh'
        stop
    end if

    idx0 = floor((x0-xx(1))/dx) + 1
    idy0 = floor((y0-yy(1))/dy) + 1

    r1 = min(1.0d0, (x0-xx(idx0))/dx )
    r2 = min(1.0d0, (y0-yy(idy0))/dy )

    v0 = (1-r1)*(1-r2)*val(idx0,idy0) + (1-r1)*r2*val(idx0,idy0+1) &
    & + r1*(1-r2)*val(idx0+1,idy0) + r1*r2*val(idx0+1,idy0+1)
end subroutine



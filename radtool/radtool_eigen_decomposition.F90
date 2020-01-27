! radtool_eigen_decomposition.F90 -- Eigen decomposition for real eigenvalues
!
! This routine computes eigenvalues and eigenvectors for a real
! asymmetric matrix, for which it is known a priori that the
! eigenvalues are real.  Such matrices arise in radiative transfer
! problems. The matrix is first balanced using the Parlett-Reinsch
! algorithm.  Then the Martin-Wilkinson algorithm is applied.
!
! This is a Fortran-90 version by Robin Hogan of the ASYMTX routine in
! DISORT (replacing the GOTOs), which is itself an adaptation of the
! complex-number routine EIGRF from the IMSL library. In turn, EIGRF
! is based primarily on EISPACK routines.
!
!       References:
!          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
!             Matrix Eigenvalue Problems, in Cowell, ed., 1984:
!             Sources and Development of Mathematical Software,
!             Prentice-Hall, Englewood Cliffs, NJ
!         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
!             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
!         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
!             Clarendon Press, Oxford

module radtool_eigen_decomposition

contains

  subroutine eigen_decomposition_real(norder, nmat, amat, &
       &                              eigenvalue, eigenvector, &
       &                              nerror, ierror)

    ! "jprb" is the precision in the input and output data and may be
    ! single or double. "jprd" is double precision for local
    ! variables.
    use parkind1, only : jprb, jprd

    implicit none

    ! Constants

    ! Tolerance
    real(kind=jprd) :: Tol = epsilon(1.0_jprd)

    ! Inputs

    ! Order of matrix
    integer,         intent(in)  :: norder
    ! Number of matrices to decompose
    integer,         intent(in)  :: nmat
    ! Matrices to be decomposed
    real(kind=jprb), intent(in)  :: amat(norder,norder,nmat)

    ! Outputs

    ! Returned eigenvalues and eigenvectors
    real(kind=jprb), intent(out) :: eigenvalue(norder,nmat)
    real(kind=jprb), intent(out) :: eigenvector(norder,norder,nmat)

    ! Return the number of matrices that could not be eigen-decomposed
    integer, optional, intent(out) :: nerror

    ! For each of the nmat matrices that could not be fully
    ! eigen-decomposed, a positive value is added to ierror indicating
    ! that eigenvalues ierror+1:norder were found but eigenvalues
    ! 1:ierror were not found.  On success, a zero entry is added to
    ! ierror.
    integer, optional, intent(out) :: ierror(nmat)

    ! Parameters
    real(kind=jprd), parameter :: &
         &  C1 = 0.4375_jprd, &
         &  C2 = 0.5_jprd, &
         &  C3 = 0.75_jprd, &
         &  C4 = 0.95_jprd, &
         &  C5 = 16.0_jprd, &
         &  C6 = 256.0_jprd
    
    ! Local variables

    ! Balanced version of amat
    real(kind=jprd) :: abal(norder,norder)

    ! Double-precision eigenvalues and eigenvectors
    real(kind=jprd) :: eigenval(norder)
    real(kind=jprd) :: eigenvec(norder,norder)

    ! Discriminant of polynomial, and a function of it
    real(kind=jprd) :: discriminant, half_sqrt_disc

    ! Loop index for matrices
    integer :: jmat

    ! Other loop indices
    integer :: ji, jj, jn

    ! Other indices
    integer :: ka, kk, ll, lb, lll, n1, n2, in, kkk

    ! Local scalars
    real(kind=jprd) :: column, ff, gg, hh
    real(kind=jprd) :: pp, qq, rr, tmp, rnorm, row, ss, scale, tt
    real(kind=jprd) :: uu, vv, ww, xx, yy, zz

    ! Work space
    real(kind=jprd) :: wkd(2*norder)

    logical :: not_finished, not_found, not_last, is_error

    if (present(nerror)) then
      nerror = 0;
    end if
    if (present(ierror)) then
      ierror(:) = 0
    end if
    
    if (norder > 2) then
      ! Nominal case: matrices of order larger than 2

      ! Loop over input matrices
      do jmat = 1,nmat

        is_error = .false.

        ! Initialize outputs
        eigenval(:)    = 0.0_jprd
        eigenvec(:,:)  = 0.0_jprd
        do ji = 1,norder
          eigenvec(ji,ji) = 1.0_jprd
        end do

        ! Balance the input matrix and reduce its norm by diagonal
        ! similarity transformation stored in WK; then search for rows
        ! isolating an eigenvalue and push them down

        ! Initialize balanced matrix
        abal = amat(:,:,jmat)

        rnorm = 0.0_jprd

        ll = 1
        kk = norder

        ! The flow of the following is a bit complicated, but it is
        ! better than the original which involved GOTO statements
        not_finished = .true.

        do while(not_finished) 
          not_finished = .false.

          kkk = kk

          do jj = kkk,1,-1
            row = 0.0_jprd

            do ji = 1,kk
              if (ji /= jj) then
                row = row + abs(abal(jj,ji))
              end if
            end do

            if (row == 0.0_jprd) then

              wkd(kk) = jj
              if (ji /= kk) then
                do ji = 1,kk
                  tmp         = abal(ji,jj)
                  abal(ji,jj) = abal(ji,kk)
                  abal(ji,kk) = tmp
                end do
                do ji = ll,norder
                  tmp         = abal(jj,ji)
                  abal(jj,ji) = abal(kk,ji)
                  abal(kk,ji) = tmp
                end do
              end if

              kk = kk - 1
              not_finished = .true.
              exit ! Quit the jj loop and repeat
            else
              not_finished = .false.
            end if

          end do ! jj
        end do ! while not_finished

        ! Search for columns isolating an eigenvalue and push them
        ! left

        not_finished = .true.
        do while(not_finished)
          not_finished = .false.
          lll = ll
          do jj = lll,kk
            column = 0.0_jprd

            do ji = ll, kk
              if (ji /= jj) then
                column = column + abs(abal(ji,jj))
              end if
            end do

            if (column == 0.0_jprd) then
              wkd(ll) = jj
              if (jj /= ll) then
                do ji = 1,kk
                  tmp         = abal(ji,jj)
                  abal(ji,jj) = abal(ji,ll)
                  abal(ji,ll) = tmp
                end do

                do ji = ll, norder
                  tmp         = abal(jj,ji)
                  abal(jj,ji) = abal(ll,ji)
                  abal(ll,ji) = tmp               
                end do
              end if
              ll = ll + 1
              not_finished = .true.
              exit ! Quit the jj loop and repeat
            else
              not_finished = .false.
            end if
          end do
        end do

        ! Balance the submatrix in rows L through K
        do ji = ll,kk
          wkd(ji) = 1.0_jprd
        end do

        not_finished = .true.
        do while (not_finished)
          not_finished = .false.
          do ji = ll,kk
            column = 0.0_jprd
            row    = 0.0_jprd
            do jj = ll,kk
              if (jj /= ji) then
                column = column + abs(abal(jj,ji))
                row    = row    + abs(abal(ji,jj))
              end if
            end do

            ff = 1.0_jprd
            gg = row / C5
            hh = column + row

            do while (column < gg)
              ff = ff*C5
              column = column*C6
            end do

            gg = row*C5

            do while (column > gg)
              ff = ff / C5
              column = column/C6
            end do

            if ((column + row) / ff < C4*hh) then
              wkd(ji) = wkd(ji)*ff
              not_finished = .true.

              do jj = ll,norder
                abal(ji,jj) = abal(ji,jj) / ff
              end do
              do jj = 1,kk
                abal(jj,ji) = abal(jj,ji) * ff
              end do
            end if
          end do
        end do ! not_finished

        ! Is A already in Hessenberg form?
        if (kk - 1 >= ll + 1) then
          ! No: convert A to Hessenberg form

          do jn = ll+1, kk-1
            hh = 0.0_jprd
            wkd(jn+norder) = 0.0_jprd
            scale = 0.0_jprd

            ! Scale column
            do ji = jn,kk
              scale = scale + abs(abal(ji,jn-1))
            end do

            if (scale /= 0.0_jprd) then

              do ji = kk,jn,-1
                wkd(ji+norder) = abal(ji,jn-1) / scale
                hh = hh + wkd(ji+norder)*wkd(ji+norder)
              end do

              gg = -sign(sqrt(hh), wkd(jn+norder))
              hh = hh - wkd(jn+norder)*gg
              wkd(jn+norder) = wkd(jn+norder) - gg

              hh = 1.0_jprd / hh
              ! Form (I-(U*UT)*H)*A
              do jj = jn,norder
                ff = 0.0_jprd
                do ji = kk,jn,-1
                  ff = ff + wkd(ji+norder)*abal(ji,jj)
                end do
                do ji = jn,kk
                  abal(ji,jj) = abal(ji,jj) - wkd(ji+norder)*ff*hh
                end do

              end do

              do ji = 1,kk
                ff = 0.0_jprd
                do jj = kk,jn,-1
                  ff = ff + wkd(jj+norder)*abal(ji,jj)
                end do
                do jj = jn,kk
                  abal(ji,jj) = abal(ji,jj) - wkd(jj+norder)*ff*hh
                end do
              end do

              wkd(jn+norder)  = scale*wkd(jn+norder)
              abal(jn,jn-1) = scale*gg
            end if
          end do

          do jn = kk-2,ll,-1 ! 360
            n1 = jn+1
            n2 = jn+2
            ff = abal(n1,jn)

            if (ff /= 0.0_jprd) then
              ff = ff*wkd(jn+1+norder)
              do ji = n2,kk
                wkd(ji+norder) = abal(ji,jn)
              end do

              if (n1 < kk) then
                do jj = 1,norder
                  gg = 0.0_jprd
                  do ji = n1,kk
                    gg = gg + wkd(ji+norder)*eigenvec(ji,jj)
                  end do

                  gg = gg / ff

                  do ji = n1,kk
                    eigenvec(ji,jj) = eigenvec(ji,jj) + gg*wkd(ji+norder)
                  end do
                end do

              end if
            end if

          end do

        end if

        jn = 1
        do ji = 1,norder
          do jj = jn,norder
            rnorm = rnorm + abs(abal(ji,jj))
          end do
          jn = ji
          if (ji < ll .or. ji > kk) then
            eigenval(ji) = abal(ji,ji)
          end if
        end do
        jn = kk
        tt = 0.0_jprd

        ! Search for next eigenvalue
        not_finished = .true.
        do while (not_finished .and. jn >= ll)
          not_finished = .false.
          in = 0
          n1 = jn-1
          n2 = jn-2

          ! Look for single small sub-diagonal element
          not_found = .true.
          do while (not_found) ! 410
            !not_found = .false.

            do ji = ll,jn
              lb = jn + ll - ji
              if (lb == ll) then
                exit
              end if
              ss = abs(abal(lb-1,lb-1)) + abs(abal(lb,lb))
              if (ss == 0.0_jprd) then
                ss = rnorm
              end if
              if (abs(abal(lb,lb-1)) < tol*ss) then
                exit
              end if
            end do

            xx = abal(jn,jn)
            if (lb == jn) then
              abal(jn,jn) = xx + tt
              eigenval(jn) = abal(jn,jn)
              jn = n1
              not_finished = .true.
              exit
            end if

            yy = abal(n1,n1)
            ww = abal(jn,n1) * abal(n1,jn)

            if (lb == n1) then
              ! Two eigenvalues found
              pp = (yy - xx) * C2
              qq = pp*pp + ww
              zz = sqrt(abs(qq))
              abal(jn,jn) = xx + tt
              xx = abal(jn,jn)
              abal(n1,n1) = yy + tt
              ! Real pair
              zz = pp + sign(zz,pp)
              eigenval(n1) = xx + zz
              eigenval(jn) = eigenval(n1)

              if (zz /= 0.0_jprd) then
                eigenval(jn) = xx - ww/zz
              end if

              xx = abal(jn,n1)
              ! Employ scale factor in case xx and zz are very
              ! small
              rr = 1.0_jprd / sqrt(xx*xx + zz*zz)
              pp = xx * rr
              qq = zz * rr

              ! Row modification
              do jj = n1,norder
                zz = abal(n1,jj)
                abal(n1,jj) = qq*zz + pp*abal(jn,jj)
                abal(jn,jj) = qq*abal(jn,jj) - pp*zz
              end do

              ! Column modification
              do ji = 1,jn
                zz = abal(ji,n1)
                abal(ji,n1) = qq*zz + pp*abal(ji,jn)
                abal(ji,jn) = qq*abal(ji,jn) - pp*zz
              end do

              ! Accumulate transforms
              do ji = ll,kk
                zz = eigenvec(ji,n1)
                eigenvec(ji,n1) = qq*zz + pp*eigenvec(ji,jn)
                eigenvec(ji,jn) = qq*eigenvec(ji,jn) - pp*zz
              end do
              jn = n2
              not_finished = .true.
              exit
            end if

            if (in == 30) then
              ! no convergence after 30 iterations; set error
              ! indicator to the index of the current eigenvalue
              if (present(ierror)) then
                ierror(jmat) = jn
              end if
              if (present(nerror)) then
                nerror = nerror + 1
              end if
              is_error = .true.
              not_finished = .false.
              exit
            end if

            ! Form shift
            if (in == 10 .or. in == 20) then
              tt = tt+xx
              do ji = ll,jn
                abal(ji,ji) = abal(ji,ji) - xx
              end do

              ss = abs(abal(jn,n1)) + abs(abal(n1,n2))
              xx = C3*ss
              yy = xx
              ww = -C1 * ss*ss
            end if

            in = in + 1
            ! Look for two consecutive small sub-diagonal
            ! elements
            do jj = lb,n2
              ji = n2+lb-jj
              zz = abal(ji,ji)
              rr = xx-zz
              ss = yy-zz
              pp = (rr*ss - ww) / abal(ji+1,ji) + abal(ji,ji+1)
              qq = abal(ji+1,ji+1) -zz-rr-ss
              rr = abal(ji+2,ji+1)
              ss = 1.0_jprd / (abs(pp) + abs(qq) + abs(rr))
              pp = pp * ss
              qq = qq * ss
              rr = rr * ss
              if (ji == lb) then
                exit
              end if
              uu = abs(abal(ji,ji-1)) * (abs(qq) + abs(rr))
              vv = abs(pp) * (abs(abal(ji-1,ji-1)) + abs(zz) &
                   &              +abs(abal(ji+1,ji+1)))
              if (uu <= tol*vv) then
                exit
              end if
            end do

            abal(ji+2,ji) = 0.0_jprd
            do jj = ji+3,jn
              abal(jj,jj-2) = 0.0_jprd
              abal(jj,jj-3) = 0.0_jprd
            end do

            ! Double QR step involving rows K to N and columns M
            ! to N
            do ka = ji,n1
              not_last = ka /= n1
              if (ka == ji) then
                ss = sign(sqrt(pp*pp + qq*qq + rr*rr), pp)
                if (lb /= ji) then
                  abal(ka,ka-1) = -abal(ka,ka-1)
                end if
              else
                pp = abal(ka,ka-1)
                qq = abal(ka+1,ka-1)
                rr = 0.0_jprd

                if (not_last) then
                  rr = abal(ka+2,ka-1)
                end if

                xx = abs(pp) + abs(qq) + abs(rr)
                if (xx == 0.0_jprd) then
                  cycle
                end if
                pp = pp / xx
                qq = qq / xx
                rr = rr / xx
                ss = sign(sqrt(pp*pp + qq*qq + rr*rr), pp)
                abal(ka,ka-1) = -ss*xx
              end if

              pp = pp + ss
              ss = 1.0_jprd / ss
              xx = pp*ss
              yy = qq*ss
              zz = rr*ss
              pp = 1.0_jprd / pp
              qq = qq*pp
              rr = rr*pp

              ! Row modification
              do jj = ka,norder
                pp = abal(ka,jj) + qq*abal(ka+1,jj)

                if (not_last) then
                  pp = pp + rr*abal(ka+2,jj)
                  abal(ka+2,jj) = abal(ka+2,jj) - pp*zz
                end if

                abal(ka+1,jj) = abal(ka+1,jj) - pp*yy
                abal(ka,jj)   = abal(ka,jj)   - pp*xx
              end do

              ! Column modification
              do ji = 1,min(jn,ka+3)
                pp = xx*abal(ji,ka) + yy*abal(ji,ka+1)

                if (not_last) then
                  pp = pp + zz*abal(ji,ka+2)
                  abal(ji,ka+2) = abal(ji,ka+2) - pp*rr
                end if

                abal(ji,ka+1) = abal(ji,ka+1) - pp*qq
                abal(ji,ka)   = abal(ji,ka)   - pp
              end do

              ! Accumulate transformations
              do ji = ll,kk
                pp = xx*eigenvec(ji,ka) + yy*eigenvec(ji,ka+1)

                if (not_last) then
                  pp = pp + zz*eigenvec(ji,ka+2)
                  eigenvec(ji,ka+2) = eigenvec(ji,ka+2) - pp*rr
                end if
                eigenvec(ji,ka+1) = eigenvec(ji,ka+1) - pp*qq
                eigenvec(ji,ka)   = eigenvec(ji,ka)   - pp

              end do
            end do
          end do ! do while not_found
        end do ! do while not_finished
        ! All eigenvalues found, now backsubstitute real vector
        
        ! If an error occurred with this matrix, save what we have and
        ! move to next one
        if (.not. is_error) then

          if (rnorm /= 0.0_jprd) then

            do jn = norder,1,-1
              n2 = jn
              abal(jn,jn) = 1.0_jprd

              do ji = jn-1,1,-1
                ww = abal(ji,ji) - eigenval(jn)
                if (abs(ww) < abs(Tol*rnorm)) then
                  ww = sign(Tol*rnorm, ww)
                end if
                rr = abal(ji,jn)
                do jj = n2,jn-1
                  rr = rr + abal(ji,jj) * abal(jj,jn)
                end do
                abal(ji,jn) = -rr / ww
                n2 = ji
              end do
            end do

            ! End backsubstitution vectors of isolated evals
            do ji = 1,norder
              if (ji < ll .or. ji > kk) then
                do jj = ji,norder
                  eigenvec(ji,jj) = abal(ji,jj)
                end do
              end if
            end do

            ! Multiply by transformation matrix
            if (kk /= 0.0_jprd) then
              do jj = norder,ll,-1
                do ji = ll,kk
                  zz = 0.0_jprd
                  do jn = ll,min(jj,kk)
                    zz = zz + eigenvec(ji,jn) * abal(jn,jj)
                  end do
                  eigenvec(ji,jj) = zz
                end do
              end do
            end if
          end if

          do ji = ll,kk
            do jj = 1,norder
              eigenvec(ji,jj) = eigenvec(ji,jj) * wkd(ji)
            end do
          end do

          ! Interchange rows if permutations occurred

          do ji = ll-1,1,-1
            jj = nint(wkd(ji))

            if (ji < jj) then
              do jn = 1,norder
                tmp = eigenvec(ji,jn)
                eigenvec(ji,jn) = eigenvec(jj,jn)
                eigenvec(jj,jn) = tmp
              end do
            end if
          end do

          do ji = kk+1,norder
            jj = nint(wkd(ji))

            if (ji /= jj) then
              do jn = 1,norder
                tmp = eigenvec(ji,jn)
                eigenvec(ji,jn) = eigenvec(jj,jn)
                eigenvec(jj,jn) = tmp
              end do
            end if

          end do
        end if

        ! Put results into output arrays
        eigenvalue(:,jmat) = eigenval(:)
        eigenvector(:,:,jmat) = eigenvec(:,:)

      end do ! jmat

    else if (norder == 2) then
      ! Special case: 2x2 matrices
      do jmat = 1,nmat
        discriminant = (amat(1,1,jmat)-amat(2,2,jmat)) ** 2 &
             & + 4.0_jprd*amat(1,2,jmat)*amat(2,1,jmat)
        if (discriminant < 0.0_jprd) then
          if (present(nerror)) then
            nerror = nerror + 1
          end if
          if (present(ierror)) then
            ierror(jmat) = 2
          end if
        end if

        eigenvalue(1,jmat) = 0.5_jprd*(amat(1,1,jmat)+amat(2,2,jmat))
        eigenvalue(2,jmat) = eigenvalue(1,jmat)
        half_sqrt_disc = 0.5_jprd*sqrt(discriminant)
        if (amat(1,1,jmat) >= amat(2,2,jmat)) then
          eigenvalue(1,jmat) = eigenvalue(1,jmat) + half_sqrt_disc
          eigenvalue(2,jmat) = eigenvalue(2,jmat) - half_sqrt_disc
        else
          eigenvalue(1,jmat) = eigenvalue(1,jmat) - half_sqrt_disc
          eigenvalue(2,jmat) = eigenvalue(2,jmat) + half_sqrt_disc
        end if
        eigenvector(1,1,jmat) = 1.0_jprd
        eigenvector(2,2,jmat) = 1.0_jprd
        if (amat(1,1,jmat) == amat(2,2,jmat) .and. &
             &  (amat(2,1,jmat) == 0.0_jprd .or. amat(1,2,jmat) == 0.0_jprd)) then
          rnorm = 1.0_jprd / (Tol * &
               &  abs(amat(1,1,jmat)) + abs(amat(2,1,jmat)) &
               & +abs(amat(1,2,jmat)) + abs(amat(2,2,jmat)))
          eigenvector(2,1,jmat) = amat(2,1,jmat) * rnorm
          eigenvector(1,2,jmat) = amat(1,2,jmat) * rnorm
        else
          eigenvector(2,1,jmat) = amat(2,1,jmat) / (eigenvalue(1,jmat)-amat(2,2,jmat))
          eigenvector(1,2,jmat) = amat(1,2,jmat) / (eigenvalue(2,jmat)-amat(1,1,jmat))
        end if

      end do

    else if (norder == 1) then
      ! Special case: 1x1 matrices
      eigenvalue(1,1:nmat) = amat(1,1,1:nmat)
      eigenvector(1,1,1:nmat) = 1.0_jprd
    else
      ! norder not a positive number
      if (present(nerror)) then
        nerror = nerror + 1
      end if

      if (present(ierror)) then
        ierror(:) = 1
      end if
    end if

  end subroutine eigen_decomposition_real

end module radtool_eigen_decomposition

  

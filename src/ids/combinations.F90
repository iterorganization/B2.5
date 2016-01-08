module combinations

  use ggd_assert

  implicit none      

contains

  !> Service routines for counting and enumerating combinations of a vector 
  !> (i_1, ..., i_n), with the range of every component  i_j is 0 <= i_j <= d_j and
  !> for the sum of all components holds sum(i_j) = dsum for a given value of dsum.
  !>
  !> These routines are used to build lists of all possible object classes of a given
  !> dimension for a given grid. d_j is then the dimension of space j (stored in the vector dmax), 
  !> and dsum the dimension of the grid objects under consideration.

  !> Compute the number of possible combinations for the vector.
  integer recursive function count_combinations(dmax, dsum) result(ccount)
    integer, intent(in) :: dmax(:), dsum

    ! internal
    integer :: iFirst, iRemaining

    call assert(sum(dmax) >= dsum)

    ! Extreme case: vector of length 1
    ! Only one combination (the value of dsum)
    if (size(dmax) == 1) then
        ccount = 1
        return
    end if

    ! Vary the first component of the vector in the allowed range
    ! and sum up the possible combinations of the remaining part of the vector
    ccount = 0
    do iFirst = min(dmax(1), dsum), 0, -1
        iRemaining = dsum - iFirst
        if ( sum(dmax(2:)) < iRemaining ) exit
        ccount = ccount + count_combinations(dmax(2:), iRemaining)
    end do
  end function count_combinations


  !> Build a list of all possible combinations of the vector.
  recursive function enumerate_combinations(dmax, dsum, cCountTotal) result(comb)
    integer, intent(in) :: dmax(:), dsum, cCountTotal
    integer :: comb(cCountTotal, size(dmax))

    ! internal
    integer :: iFirst, cCount, cSubCount, iRemaining

    ! Extreme case: vector of length 1
    ! Only one combination (the value of dsum)
    if (size(dmax) == 1) then
        comb(1,1) = dsum
        return
    end if


    ! Vary the first component of the vector in the allowed range
    ! and recursively assemble the possible combinations of the remaining part of the vector
    cCount = 0
    do iFirst = min(dmax(1), dsum), 0, -1
        iRemaining = dsum - iFirst
        if ( sum(dmax(2:)) < iRemaining ) exit
        cSubCount = count_combinations(dmax(2:), iRemaining)
        comb(cCount + 1 : cCount + cSubCount, 1) = iFirst
        comb(cCount + 1 : cCount + cSubCount, 2:) = enumerate_combinations(dmax(2:), iRemaining, cSubCount)
        cCount = cCount + cSubCount
    end do

    call assert( cCount == cCountTotal )
  end function enumerate_combinations


  !> Convenience routine: allocate and populate an array with all possible combinations.
  subroutine allocate_combinations(dmax, dsum, comb)
    integer, intent(in) :: dmax(:), dsum
    integer, allocatable :: comb(:, :)

    ! internal
    integer :: ccount 

    ccount = count_combinations(dmax, dsum)
    allocate( comb(ccount, size(dmax)) )
    comb = enumerate_combinations(dmax, dsum, ccount)

  end subroutine allocate_combinations

end module combinations

!!!Local Variables:
!!! mode: f90
!!! End:

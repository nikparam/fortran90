MODULE APPEND_ARRAY
	CONTAINS
		SUBROUTINE append(list, element)

		INTEGER :: i, isize, astatus
		DOUBLE PRECISION, INTENT(IN) :: element
		DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: list(:)
		DOUBLE PRECISION, ALLOCATABLE :: clist(:)

		isize = SIZE(list)
		ALLOCATE( clist(isize + 1), STAT = astatus )
		clist(1:isize) = list(1:isize)
		clist(isize + 1) = element

		DEALLOCATE(list)
		CALL MOVE_ALLOC(clist, list)

		END SUBROUTINE append
END MODULE APPEND_ARRAY


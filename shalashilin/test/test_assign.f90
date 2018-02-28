PROGRAM test

	INTEGER :: Num

	Num = 4

	x = printing(Num)
	
CONTAINS

	FUNCTION printing(Num)

	INTEGER :: Num
	INTEGER :: Num1 

	Num1 = Num

	PRINT *, Num1

	printing = 1

	END FUNCTION printing

END PROGRAM test

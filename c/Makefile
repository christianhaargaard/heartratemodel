#libsol.so: solution.o
	#gcc -shared -I/usr/local/include-Wl,-install_name,libsol.so -o libsol.so solution.o /usr/local/lib/libsundials_nvecserial.a /usr/local/lib/libsundials_cvodes.a

#solution.o : solution.c
	#gcc -shared -Wall -I/usr/local/include /usr/local/lib/libsundials_cvodes.a /usr/local/lib/libsundials_nvecserial.a solution.c -lgdbm

libsol.so: solution.o
	#gcc -g -O3 -ftree-vectorize -fomit-frame-pointer -ffast-math -funsafe-math-optimizations  -Wall -Wall -shared -o libsol.so solution.o -L/usr/local/Cellar/sundials/2.5.0/lib -lsundials_cvodes -lsundials_nvecserial -lm
	gcc -g -O3 -ftree-vectorize -fomit-frame-pointer -funsafe-math-optimizations  -Wall -Wall -shared -o libsol.so solution.o -L/usr/local/Cellar/sundials/2.5.0/lib -lsundials_cvodes -lsundials_nvecserial -lm
	#gcc -shared -o libsol.so solution.o /usr/local/lib/libsundials_nvecserial.a /usr/local/lib/libsundials_cvodes.a

solution.o: solution.c
	#gcc -g -O3 -ftree-vectorize -fomit-frame-pointer -ffast-math -funsafe-math-optimizations  -Wall -c solution.c 
	gcc -g -O3 -ftree-vectorize -fomit-frame-pointer -funsafe-math-optimizations  -Wall -c solution.c 

#sundials-config -m cvodes -t s -l c
#-I/usr/local/Cellar/sundials/2.5.0/include
#-L/usr/local/Cellar/sundials/2.5.0/lib -lsundials_cvodes -lsundials_nvecserial -lm

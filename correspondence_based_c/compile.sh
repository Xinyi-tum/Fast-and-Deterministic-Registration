cd c_FBnB
gcc -c main.c -o main.o
gcc -c bnbandvoting.c -o bnbandvoting.o
gcc -o FBnB_registration bnbandvoting.o main.o  -lm

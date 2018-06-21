OBJS=run.o kmldm.o utils.o gmm2.o mldm-openmp.o pqn.o quic.o lbfgs.o
COMPILE=g++ -std=c++11
run: $(OBJS)
	$(COMPILE) -o run $(OBJS) -llapack -lblas -fopenmp
run.o: run.cpp kmldm.h
	$(COMPILE) -I../eigen/ -c -o run.o run.cpp
kmldm.o: kmldm.cpp kmldm.h utils.h gmm2.h ./lbfgs/lbfgs.h
	$(COMPILE) -I../eigen/ -c -o kmldm.o kmldm.cpp
mldm-openmp.o: mldm-openmp.cpp mldm.h utils.h
	$(COMPILE) -I../eigen/ -c -o mldm-openmp.o mldm-openmp.cpp -fopenmp
gmm2.o: gmm2.cpp gmm2.h utils.h
	$(COMPILE) -I../eigen/ -c -o gmm2.o gmm2.cpp
pqn.o: proximalQusiNewton.cpp proximalQusiNewton.h utils.h
	$(COMPILE) -I../eigen/ -c -o pqn.o proximalQusiNewton.cpp
quic.o: ./QUIC/QUIC.h ./QUIC/QUIC.cpp utils.h
	$(COMPILE) -I../eigen/ -c -o quic.o ./QUIC/QUIC.cpp -fopenmp
lbfgs.o: ./lbfgs/lbfgs.h ./lbfgs/arithmetic_ansi.h ./lbfgs/lbfgs.c 
	$(COMPILE) -c -o lbfgs.o ./lbfgs/lbfgs.c
utils.o: utils.cpp utils.h
	$(COMPILE) -I../eigen/ -c -o utils.o utils.cpp
clean:
	rm -rf $(OBJS) run

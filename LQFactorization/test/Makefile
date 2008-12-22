mkl-include-path=/people/thnfs/homes/wanghan/local/intel/mkl/10.0.3.020/include
mkl-lib-path=/people/thnfs/homes/wanghan/local/intel/mkl/10.0.3.020/lib/em64t

run:	LQFactorization.o
	ld --library-path=$(mkl-lib-path) -rpath  $(mkl-lib-path) -rpath-link $(mkl-lib-path)  --eh-frame-hdr -m elf_x86_64 -dynamic-linker /lib64/ld-linux-x86-64.so.2 /usr/lib64/gcc/x86_64-suse-linux/4.2.1/../../../../lib64/crt1.o /usr/lib64/gcc/x86_64-suse-linux/4.2.1/../../../../lib64/crti.o /usr/lib64/gcc/x86_64-suse-linux/4.2.1/crtbegin.o -L/people/thnfs/homes/wanghan/local/intel/mkl/10.0.3.020/lib/em64t -L/usr/lib64/gcc/x86_64-suse-linux/4.2.1 -L/usr/lib64/gcc/x86_64-suse-linux/4.2.1/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib64/gcc/x86_64-suse-linux/4.2.1/../../../../x86_64-suse-linux/lib -L/usr/lib64/gcc/x86_64-suse-linux/4.2.1/../../.. LQFactorization.o -lmkl_lapack -lmkl -lguide -lpthread -lstdc++ -lm -lgcc_s -lgcc -lc -lgcc_s -lgcc /usr/lib64/gcc/x86_64-suse-linux/4.2.1/crtend.o /usr/lib64/gcc/x86_64-suse-linux/4.2.1/../../../../lib64/crtn.o 



LQFactorization.o: LQFactorization.cpp LQFactorization.h
	g++ -c LQFactorization.cpp -Wall -I$(mkl-include-path)


# run:
# 	g++ -g -v LQFactorization.cpp -Wall -I$(mkl-include-path) -L$(mkl-lib-path) -lmkl_lapack -lm -lmkl -lguide -lpthread 
#-mgnu-ld="-rpath /people/thnfs/homes/wanghan/local/intel/mkl/10.0.3.020/lib/em64t"
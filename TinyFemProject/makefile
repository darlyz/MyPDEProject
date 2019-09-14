CC = gcc
FLAG = -O2
target = tinyfem
trilinos = /public/software/Trilib
solvlib = $(trilinos)/lib
solvinc = $(trilinos)/include

src = $(wildcard *.c)
obj = $(patsubst %.c, %.o, $(src))

$(target):$(obj)
	$(CC) $^ -o $@ $(solvlib)/libaztecoo.a -llapack -lblas -lm

%.o:%.c
	$(CC) -c $< -I$(solvinc) -o $@

.PHONY:clean
clean:
	rm $(obj) $(target) -f
wclean:
	del $(obj) $(target)


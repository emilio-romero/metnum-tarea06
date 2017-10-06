CFLAGS= -c -std=c99 -lm
CPPS:=$(wildcard src/*.c)
OBJS:=$(patsubst src/%.c,obj/%.o,$(CPPS))
OBJDIR:=obj

ejecutable:$(OBJS) src/lectura.h 
	gcc -o $@ $(OBJS) -std=c99 -lm 
$(OBJDIR)/%.o: src/%.c
	gcc $(CFLAGS) $< -o $@
#main.o: main.c 
#	gcc -c main.c -std=c99
#lectura.o: lectura.c 
#	gcc -c lectura.c -std=c99
#eigenv.o: eigenv.c eigenv.h 
#	gcc -c eigenv.c -std=c99 -lm
clean:
	rm $(OBJS)

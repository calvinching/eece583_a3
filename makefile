# This makefile is written for gcc running under Solaris on a SPARCstation.
# To compile on other systems, you may have to change:
# (1) CC to the name of your C compiler.
# (2) LIB_DIR to point at your directory of X11 libraries (libX11.a, etc.)
# (3) On many systems, the -R/usr/openwin/lib option on the LIB line
#     will have to be removed.
# (4) X11_INCLUDE to point at the directory containing "x11/xlib.h" etc.
# (5) FLAGS should be changed to whatever options turn on maximum optimization
#     in your compiler.
PLATFORM = X11
#PLATFORM = WIN32
#PLATFORM = NO_GRAPHICS

HDR = graphics.h
SRC = graphics.c example.c
EXE = example
FLAGS = -g -Wall -Wno-write-strings -D$(PLATFORM)

# Need to tell the linker to link to the X11 libraries.
# WIN32 automatically links to the win32 API libraries (no need for flags)
ifeq ($(PLATFORM),X11)
   GRAPHICS_LIBS = -lX11
endif

$(EXE): graphics.o example.o
	g++ $(FLAGS) graphics.o example.o $(GRAPHICS_LIBS) -o $(EXE)

graphics.o: graphics.c $(HDR)
	g++ -c $(FLAGS) graphics.c

example.o: example.c $(HDR)
	g++ -c $(FLAGS) example.c

clean:
	rm $(EXE) *.o


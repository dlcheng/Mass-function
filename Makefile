#--------------------------------------------
OPT += -DST_MF                                 # the ST fitting foumula
#OPT += -DPLANCK_500                     # using the fitting mass fucntion of the planck, Delta = 500

#--------------------------------------- Select target computer
#SYSTYPE="Workstation"
SYSTYPE="Mac"

#--------------------------------------- Adjust settings for target computer
ifeq ($(SYSTYPE),"Workstation")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall 
GSL_INCL =  -I/home/dlcheng/Install/gsl/include
GSL_LIBS =  -L/home/dlcheng/Install/gsl/lib
endif

ifeq ($(SYSTYPE),"Mac")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall
GSL_LIBS=   -L/usr/local/Cellar/gsl/1.16/lib  -Wl
GSL_INCL =  -I/usr/local/Cellar/gsl/1.16/include
endif

OPTIONS =  $(OPTIMIZE) $(OPT) 

EXEC   =  mass_function

OBJS   =  allvars.o cdm_mf.o filters.o growth.o init_all.o main.o set_params.o tk_bbks.o tk.o tools.o \
          variance.o 

INCL   = allvars.h proto.h define.h Makefile

CFLAGS = $(OPTIONS) $(GSL_INCL) 

LIBS   = $(GSL_LIBS) -lgsl -lgslcblas -lm 

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC) *.gch

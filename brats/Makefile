###############################################################################################
#								 	                      #
#                               BRATS INSTALLATION INSTRUCTIONS                               #
#								  	                      #
#  BRATS is dependent on 2 non standard libraries, both of which must first be obtained       #
#  before installation. These are PGPLOT and Funtools. We have generally found that using     #
#  standard repository tools (e.g. yum or apt-get) to install these libraries causes the      # 
#  least compatibility issues and provides the easiest installation.			      #
#                                                                                             #
#  Other requirements, which come packaged with most Linux builds, include: GCC, OpenMP, GSL, #
#  and X11. 										      #
#											      #
#  To install BRATS, the paths to these dependencies must be updated below. The required      #
#  paths are as follows:								      #
#											      #
#  INSTALLFROM - Directory containing the install files					      #
#  INSTALLTO   - Directory which you wish BRATS to be installed to			      #
#  PGPLOT      - Directory where PGPLOT is installed					      #
#  LIBS        - Directory containing standard libraries				      #
#  X11LIBS     - Directory containing the X11 libraries					      #
#  INCLUDE     - Directory containing any additional libraries or files			      #
#											      #
#  It is also required that enviroment variables are set before running BRATS. We recommend   #
#  adding these variables to your shell configurations files (e.g. .cshrc .bashrc).	      #
#											      #
#  CSH Example:										      #
#  setenv LD_LIBRARY_PATH /usr/local/pgplot						      #
#  setenv PGPLOT_DIR /usr/local/pgplot/							      #
#											      #
#  BASH Example:									      #
#  export LD_LIBRARY_PATH=/usr/local/pgplot						      #
#  export PGPLOT_DIR=/usr/local/pgplot/							      #
#  											      #
#  Non-standard machine setups may require additional fine tuning of the commands directly.   #
#  											      #
#  BRATS then be installed using the 'make brats'.					      #
#											      #
#                                    ***** IMPORTANT *****				      #
#											      #
#  If you are upgrading from an older version of BRATS, please use the 'make brats_upgrade'   #
#  command. Failure to do so may lead to undesired results.	    	      		      #
#											      #
#  Support for external installation is currently limited. If you have any major issues,      #
#  please get in touch via the BRATS website at http://www.askanastronomer.co.uk/brats        #
#  However, please be aware the installation support will be done on an as and when basis.    #
#											      #
#  PGPLOT: http://www.astro.caltech.edu/~tjp/pgplot/					      #
#  Funtools: https://www.cfa.harvard.edu/~john/funtools/				      #
#											      #
###############################################################################################



INSTALLFROM = /home/jharwood/brats_test_environment/
INSTALLTO = /home/jharwood/brats_test_environment/testinstall
PGPLOT = /soft/pgplot
LIBS = /soft/lib
INCLUDE = /soft/include
X11LIBS = /usr/lib
PREVVER = v2.6.3
VERSION = v2.6.3.3
#FUNTOOLS = /soft/lib
# /soft/lib/libfuntools.a

brats: main.c
	mkdir -p $(INSTALLTO)
	mkdir -p $(INSTALLTO)/docs
	gcc -fopenmp -Wall -g -O0 -c main.c -Wall -Wno-format-security -c $(INSTALLFROM)/synchrotron_src/synchrotron.c $(INSTALLFROM)/synchrotron_src/bessel.c $(INSTALLFROM)/synchrotron_src/gsl_integ.c -Wall -I$(INCLUDE) -I$(FUNTOOLS) -I$(BIN) -I$(PGPLOT) -I$(INSTALLFROM)/synchrotron_src 
	gfortran -fopenmp main.o bessel.o gsl_integ.o -L$(X11LIBS) -L$(LIBS) -L$(BIN) -L$(INCLUDE) -L$(PGPLOT) -L$(FUNTOOLS) -L$(FUNTOOLSLIBS) -lcpgplot -lpgplot -lfuntools -lm -lX11 -lpng -lgsl -lgslcblas -lreadline -o $(INSTALLTO)/brats
	cp $(INSTALLFROM)/updates.txt $(INSTALLTO)/docs/
	cp $(INSTALLFROM)/bratscookbook.pdf $(INSTALLTO)/docs/

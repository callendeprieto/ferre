#
#	makefile for ferre
#	Carlos 2004-2012
#

F90 = gfortran 
#F90 = /gpfs/apps/GCC/4.6.1/bin/gfortran
OPT = -O2 -fopenmp -ftree-vectorize 
#OPT= -O2 -fopenmp -pg
#OPT= -O2 -fopenmp -R /gpfs/apps/GCC/4.6.1/lib 
#OPT= -O0
#OPT = -O0 -g 

#F90 = g95
#OPT= -O2 -openmp

#F90 = sunf95
#OPT = -O3 -xopenmp 
#OPT = -g -O0


#F90 = /usr/remote/bin/f95
#OPT= -O2 

#F90 = ifort
#OPT = -O3 -xT
#OPT = -O3 -xHost -openmp -ipo -vec
#OPT = -O0 -openmp
#OPT = -O3 -openmp
#OPT= -O0 -traceback

#F90 = pgfortran
#OPT = -fast -Mconcur=allcores  -Minfo=accel 
#OPT = -O2 -mp 
#OPT = -O0 -o 
#OPT = -O0 -g -openmp 

#F90 = xlf95_r
#OPT = -O3 -qstrict -qtune=ppc970 -qarch=ppc970 -qcache=auto -qsmp=omp

#raspbian
#F90 = gfortran
#OPT = -march=armv6 -mfpu=vfp -mfloat-abi=hard


OS1 = ${shell perl -e 'my $$osVersion = `uname `; chomp($$osVersion); print "$$osVersion"'}
OS2 = ${shell perl -e 'my $$osVersion = `uname -r`; chomp($$osVersion); print "$$osVersion"'}
OS3 = ${shell perl -e 'my $$osVersion = `machine`; chomp($$osVersion); print "$$osVersion"'}
OS4 = ${shell perl -e 'my $$osVersion = `hostname`; chomp($$osVersion); print "$$osVersion"'}

SOURCE = ferre.f90 lin.f90 getee.f90 checkstat.f90 read_f.f90 \
getsigma.f90 load_control.f90 optimum.f90 chisurf.f90 median.f90 rmedian.f90 \
minlocus.f90 getaa.f90 \
qbezier.f90 qua.f90 cbezier.f90 cub.f90 nea.f90 spl.f90 likely.f90 \
filter1.f90 filter2.f90 \
flx.f90 smooth1.f90 smooth2.f90 getmin.f90 pinin.f90 getuu.f90 \
cova.f90 covsigma.f90 continuum.f90 pem.f90 \
pdfsigma.f90 nrunsigma.f90 \
getlsf.f90 lsf_gauss.f90 normal.f90 physical.f90 fsort.f90 \
setuplambda.f90 wresample.f90 convol.f90 decompress.f90 \
count_words.f90 quit.f90 

SOURCE_MODULES = share.f90 fun.f90 nelder_mead.f90 \
booklib.f90 random.f90 timer.f90 btr.f90 uob.f90 trn.f90 \
m_mrgrnk.f90 lsq.f90 rnglib.f90 ranlib.f90 pdflib.f90 mcmc.f90

MODS = $(SOURCE_MODULES:.f90=.mod)

OBJS = $(SOURCE_MODULES:.f90=.o) $(SOURCE:.f90=.o)

all: a.out ascii2bin  install

#ascii2bin converter
ascii2bin: share.mod timer.mod checkstat.o quit.o ascii2bin.f90
	${F90} ${OPT} share.f90 timer.f90 checkstat.f90 quit.f90 ascii2bin.f90 -o ascii2bin

# 2nd and final step - compiling and linking
a.out:  ${OBJS} ${SOURCE}
	${F90}  -o a.out ${OPT}  ${SOURCE_MODULES} ${SOURCE}

#1st step - compiling the modules to get .o and .mod files

%.o : %.f90 
	${F90} -c ${OPT} $*.f90

install:
	cp a.out ../bin/ferre.x
	cp ascii2bin ../bin/ascii2bin

test: 
	@ echo ${OS1} ${OS2}  ${OS3} / ${OS4} / ${F90} ${OPT}
	@ cd ../test/; time ../src/a.out >& log

clean: 
	rm -f a.out ascii2bin ${MODS} ${OBJS}



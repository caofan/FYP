#EXTRA_DIST =  batman.ini 
#AM_CFLAGS = -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -lm
AM_CFLAGS = -O0 -g -msse2 -lm
AM_CXXFLAGS = -O0 -g -msse2 -lm
#AM_CXXFLAGS = -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -lm 
RANLIB=ranlib
#AM_CPPFLAGS = -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -lm

noinst_LIBRARIES = libz.a
libz_a_SOURCES = adler32.c compress.c crc32.c gzio.c uncompr.c deflate.c trees.c \
       zutil.c inflate.c infback.c inftrees.c inffast.c\
       zlib.h zconf.h crc32.h  deflate.h zutil.h inftrees.h inflate.h inffast.h trees.h inffixed.h

bin_PROGRAMS = index rnaseq reverse bwtformatdb
index_SOURCES= bfix.cpp index.cpp BWT.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c\
 iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c\
 BWT.h TypeNLimit.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h r250.h HSPstatistic.h\
 BWTConstruct.h QSufSort.h r250.h dictionary.h inistrlib.h iniparser.h Timing.h Socket.h karlin.h\
 bfix.h 
rnaseq_SOURCES = rnaseq.cpp \
 batlib.cpp rqindex.cpp Indexes.cpp file.cpp Cmdline.cpp Hash.cpp extend.cpp bfix.cpp BWT.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c\
 iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c\
 BWT.h TypeNLimit.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h r250.h HSPstatistic.h\
 BWTConstruct.h QSufSort.h r250.h dictionary.h inistrlib.h iniparser.h Timing.h Socket.h karlin.h\
 Indexes.h  file.h bfix.h extend.h Hash.h Cmdline.h rqindex.h batlib.h common.h const.h
bwtformatdb_SOURCES= bwtformatdb.c BWT.c BWTConstruct.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c\
 iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c\
 BWT.h TypeNLimit.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h r250.h HSPstatistic.h\
 BWTConstruct.h QSufSort.h r250.h dictionary.h inistrlib.h iniparser.h Timing.h Socket.h karlin.h
#batman_SOURCES= batman.cpp BWT.c BWTConstruct.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c\
# iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c 
rnaseq_LDADD = libz.a
#decode_SOURCES= decode.cpp BWT.c BWTConstruct.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c\
# iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c
#decode_LDADD = libz.a
reverse_SOURCES = reverse.cpp 
copy:
	cp bwtformatdb ../bin
	cp reverse ../bin
	cp rnaseq ../bin
	cp index ../bin

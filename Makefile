##############################################################################

MK_TOP = /export/ciao_from_source/ciao-4.8/src
KJG = /export/ciao

include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

EXEC              = dmimgtess
LIB_FILES = 
PAR_FILES         = dmimgtess.par
INC_FILES         = 
XML_FILES         = dmimgtess.xml

SRCS	=           dmimgtess.c t_dmimgtess.c

OBJS = $(SRCS:.c=.o)

LOCAL_INC = -I$(MK_TOP)/da/analysis/dettools/vtpdetect \
	-I$(MK_TOP)/da/analysis/dettools/Voronoi_lib \
	-I$(MK_TOP)/da/analysis/dmtools/dmimgio/


VTP_OBJS = \
	$(MK_TOP)/da/analysis/dettools/vtpdetect/vtpdetect_voroni.o \
	$(MK_TOP)/da/analysis/dettools/vtpdetect/vtpdetect_input.o \
	$(MK_TOP)/da/analysis/dettools/vtpdetect/vtpdetect_compare.o \
	$(MK_TOP)/da/analysis/dettools/vtpdetect/vtpdetect_interface.o \
	$(MK_TOP)/da/analysis/dettools/vtpdetect/vtpdetect_cleanup.o

LOCAL_LIBS = $(VTP_OBJS) \
	-L$(MK_TOP)/da/analysis/dettools/Voronoi_lib/ -lVoronoi \
	-L$(MK_TOP)/da/analysis/dmtools/dmimgio/ -ldmimgio


MAKETEST_SCRIPT   = dmimgtess.t

include $(MK_TOP)/Makefile.all

#-----------------------------------------------------------------------
# 			MAKEFILE DEPENDENCIES	
#-----------------------------------------------------------------------



$(EXEC): $(OBJS) 
	$(LINK)
	@echo


dmimgtess.o: .vtp

.vtp:
	(cd $(MK_TOP)/da/analysis/dettools/Voronoi_lib; $(MAKE) $(MKMACROS) )
	(cd $(MK_TOP)/da/analysis/dettools/vtpdetect; $(MAKE) $(MKMACROS) $(VTP_OBJS) )


kjg: $(EXEC)
	/bin/cp -f $(EXEC) $(KJG)/binexe/
	/bin/cp -f $(KJG)/bin/dmlist $(KJG)/bin/$(EXEC)
	/bin/cp -f $(PAR_FILES) $(KJG)/param/$(PAR_FILES)



announce1:
	@echo "   /----------------------------------------------------------\ "
	@echo "   |                Building dmimgtess                        | "
	@echo "   \----------------------------------------------------------/ "


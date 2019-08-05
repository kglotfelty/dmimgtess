tool = dmimgtess

binexedir = $(prefix)/binexe
paramdir = $(prefix)/param
ahelpdir = $(prefix)/share/doc/xml

binexe_PROGRAMS = $(tool)

dmimgtess_SOURCES = \
	dmimgtess.c \
	heap2.c \
	t_dmimgtess.c \
	vtpdetect.h \
	vtpdetect_aots.h \
	vtpdetect_aots_voronoi_defs.h \
	vtpdetect_aots_voronoi_edgelist.c \
	vtpdetect_aots_voronoi_geometry.c \
	vtpdetect_aots_voronoi_heap.c \
	vtpdetect_aots_voronoi_main.c \
	vtpdetect_aots_voronoi_memory.c \
	vtpdetect_aots_voronoi_output.c \
	vtpdetect_aots_voronoi_voronoi.c \
	vtpdetect_cleanup.c \
	vtpdetect_compare.c \
	vtpdetect_input.c \
	vtpdetect_interface.c \
	vtpdetect_voroni.c


dmimgtess_CPPFLAGS = $(CIAO_CFLAGS)
dmimgtess_LDADD = $(CIAO_LIBS) -ldmimgio
dmimgtess_LINK = $(CXX) -o $@ -Wl,-rpath,$(prefix)/lib -Wl,-rpath,$(prefix)/ots/lib 

# problems with libstdc++ and cxcparam on Linux :(
if LINUX
dmimgtess_LDADD += -L$(prefix)/ots/lib -lstdc++
endif


dist_param_DATA = $(tool).par
dist_ahelp_DATA = $(tool).xml

install-data-hook:
	chmod a-w $(paramdir)/$(dist_param_DATA)
	chmod a-w $(ahelpdir)/$(dist_ahelp_DATA)

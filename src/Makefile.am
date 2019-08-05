tool = dmimgtess

binexedir = $(prefix)/binexe
paramdir = $(prefix)/param
ahelpdir = $(prefix)/share/doc/xml

binexe_PROGRAMS = $(tool)

dmimgtess_SOURCES = dmimgtess.c t_dmimgtess.c
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

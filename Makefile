# Detect the operating system
detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')

# Set the compiler and flags based on the operating system
# Supported OS are Mac OSX ("Darwin") and Linux ("Linux")

# macOS
ifeq ($(detected_OS), Darwin)
cxx := g++-14 -fopenmp
cflags := -Wall -ansi -pedantic
lp_lflags := -llapack -lblas
endif

# Linux
ifeq ($(detected_OS), Linux)
cxx := g++-9 -fopenmp
cflags := -Wall -ansi -pedantic
#-O3 -fsanitize=address
lp_lflags := -llapack -lblas
endif

# Lists of files and directories to be built
objdir = objs
objs   = $(objdir)/file.o $(objdir)/convert_cfg.o $(objdir)/object.o $(objdir)/lbrmt_2d.o
srcdir = src
src    = $(patsubst $(objdir)/%.o,$(srcdir)/%.cc,$(objs))
exedir = build
execs  = $(exedir)/sim_fsi $(exedir)/sim_rotate $(exedir)/sim_mix

all: $(execs)

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) | sed 's/^\(.*\):/$(objdir)\/\1:/' >Makefile.dep

# Include the file dependencies
-include Makefile.dep

# Create a library for objects
$(objdir)/libf2d.a: $(objs)
	rm -f $(objdir)/libf2d.a
	ar rs $(objdir)/libf2d.a $^

# Create directory for objects
$(shell mkdir -p $(objdir))

# Create directory for executables
$(shell mkdir -p $(exedir))

$(exedir)/sim_fsi: $(srcdir)/sim_fsi.cc $(objdir)/libf2d.a
	$(cxx) $(cflags) -o $@ $^

$(exedir)/sim_rotate: $(srcdir)/sim_rotate.cc $(objdir)/libf2d.a
	$(cxx) $(cflags) -o $@ $^

$(exedir)/sim_mix: $(srcdir)/sim_mix.cc $(objdir)/libf2d.a
	$(cxx) $(cflags) -o $@ $^

$(objdir)/%.o: $(srcdir)/%.cc
	$(cxx) $(cflags) -c -o $@ $<

# A Makefile target to remove all the built files
clean:
	rm -rf $(execs) $(objs) libf2d.a $(objdir) $(exedir)

.PHONY: clean all depend
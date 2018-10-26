PROGRAM			= pi_app
CC     			= gcc
EXTRA_CFLAGS    	= -fdata-sections -ffunction-sections -Wl,--gc-sections -Wl,--strip-all
LDFLAGS   		= "-lm"
INCLUDE         	= "-I./" 
CPPFLAGS        	= -Wall -Wextra    # helpful for writing better code (behavior-related)
CPPFLAGS      		+= $(INCLUDE)
SRCEXTS 		= .c
HDREXTS 		= .h
CFLAGS  		= -O3
RM     			= rm -f
ETAGS 			= etags
CTAGS 			= ctags

SRCDIRS := $(shell find $(SRCDIRS) -type d)
SOURCES = $(foreach d,$(SRCDIRS),$(wildcard $(addprefix $(d)/*,$(SRCEXTS))))
HEADERS = $(foreach d,$(SRCDIRS),$(wildcard $(addprefix $(d)/*,$(HDREXTS))))
OBJS    = $(addsuffix .o, $(basename $(SOURCES)))
DEPS    = $(foreach f, $(OBJS), $(addprefix $(dir $(f))., $(patsubst %.o, %.d, $(notdir $(f)))))

## Define some useful variables.
COMPILE.c   = $(CC)  $(EXTRA_CFLAGS) $(CFLAGS)   $(CPPFLAGS) -c
LINK.c      = $(CC)  $(EXTRA_CFLAGS) $(CFLAGS)   $(CPPFLAGS) $(LDFLAGS)

.PHONY: all objs tags ctags clean distclean 

# Delete the default suffixes
.SUFFIXES:

all: $(PROGRAM)

# Rules for generating object files (.o).
#----------------------------------------
objs:$(OBJS)

%.o:%.c
	$(COMPILE.c) $< -o $@

# Rules for generating the tags.
#-------------------------------------
tags: $(HEADERS) $(SOURCES)
	$(ETAGS) $(HEADERS) $(SOURCES)

ctags: $(HEADERS) $(SOURCES)
	$(CTAGS) $(HEADERS) $(SOURCES)

# Rules for generating the executable.
#-------------------------------------
$(PROGRAM):$(OBJS)
	$(LINK.c)   $(OBJS) $(EXTRA_LDFLAGS) -o $@
	@echo Type ./$@ to execute the program.

clean:
	$(RM) $(OBJS) $(PROGRAM) $(PROGRAM).exe



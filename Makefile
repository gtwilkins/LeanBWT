PREFIX = /usr/local/bin/
INC=-Isrc -Isrc/commands -Isrc/index -Isrc/shared -Isrc/transform
VPATH=src:src/commands:src/index:src/shared:src/transform

SRCS =  \
	leanbwt.cpp \
	filenames.cpp \
	index.cpp \
	index_reader.cpp \
	index_structs.cpp \
	index_writer.cpp \
	timer.cpp \
	transform.cpp \
	transform_bwt.cpp \
	transform_structs.cpp \
	transform_binary.cpp


# C++ compiler
CXX = g++
# C++ flags; passed to compiler
CXXFLAGS = -std=c++11
# Linker flags; passed to compiler
LDFLAGS = -std=c++11
# Dependency flags; passed to compiler
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
# Objects directory
OBJDIR = .o
$(shell mkdir -p $(OBJDIR) >/dev/null)
# Dependencies directory
DEPDIR = .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
# Derive objects from sources
OBJS = $(patsubst %,$(OBJDIR)/%.o,$(basename $(SRCS)))
# Derive dependencies from sources
DEPS = $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS)))

# Generic link executable
LINK.o = $(CXX) $(LDFLAGS) $(DBG) -o $@
# Generic compile object
COMPILE.cc = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(DBG) $(INC) -c -o $@
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

all: leanbwt
	@echo
	@echo 'Compile successful. Type "sudo make install" to complete intall.'

.PHONY: install
install:
	@mv -f leanbwt $(PREFIX)
	@make clean
	@echo
	@echo 'Install successful. Type "leanbwt -h" to see usage.'

.PHONY: clean
clean:
	@$(RM) -r $(OBJDIR) $(DEPDIR)

leanbwt: $(OBJS)
	$(LINK.o) $^

$(OBJDIR)/%.o : %.cpp
$(OBJDIR)/%.o : %.cpp $(DEPDIR)/%.d
	$(COMPILE.cc) $<
	$(POSTCOMPILE)

# Create dependency rule
.PRECIOUS = $(DEPDIR)/%.d
$(DEPDIR)/%.d: ;

# Include dependencies; The '-' ensures no errors
-include $(DEPS)

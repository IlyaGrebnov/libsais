PROJECT=sais
PLIBNAME=lib$(PROJECT)
PVER=2.7.1
PSOVER=2
ifeq ($(OS),Windows_NT)
  PLIBSTATIC=$(PROJECT).a
  PLIBSHARED=$(PROJECT)-$(PVER).dll
else
  PLIBSTATIC=$(PLIBNAME).a
  PLIBSHARED=$(PLIBNAME).so.$(PSOVER)
endif
PLIBS=$(PLIBSTATIC) $(PLIBSHARED)
CC=gcc
CFLAGS?=-Wall -O2
LDFLAGS?=-lm
AR?=ar
INSTALL?=install
RM?=rm -f
RMD?=$(RM) -r
PREFIX?=/usr/local
SRCS=src
DOCS?=share/doc/$(LIBNAME)
LIBS?=lib
INCLUDES?=include
MANS?=man/man1

all: $(PLIBS)

$(SRCS)/$(PLIBNAME).o: $(SRCS)/$(PLIBNAME).c
	$(CC) $(CFLAGS) -c -o $@ $^

$(PLIBSTATIC): $(SRCS)/$(PLIBNAME).o
	$(AR) rcs $@ $^

$(PLIBSHARED): $(SRCS)/$(PLIBNAME).o
	$(CC) $(CFLAGS) -shared -Wl,-soname,$@ $^ -o $@

clean:
	$(RM) $(SRCS)/$(PLIBNAME).o $(PLIBS)

install:
	$(INSTALL) -d $(PREFIX)/$(LIBS)
	$(INSTALL) -d $(PREFIX)/$(INCLUDES)
	$(INSTALL) -d $(PREFIX)/$(MANS)
	$(INSTALL) -d $(PREFIX)/$(DOCS)
	$(INSTALL) -m 0644 $(PLIBS) $(PREFIX)/$(LIBS)
	$(INSTALL) -m 0644 $(SRCS)/$(PLIBNAME).h $(PREFIX)/$(INCLUDES)
	$(INSTALL) -m 0644 CHANGES LICENSE README.md VERSION $(PREFIX)/$(DOCS)

uninstall:
	$(RM) $(PREFIX)/$(LIBS)/$(PLIBSTATIC)
	$(RM) $(PREFIX)/$(LIBS)/$(PLIBSHARED)
	$(RM) $(PREFIX)/$(INCLUDES)/$(SRCS)/$(PLIBNAME).h
	$(RMD) $(PREFIX)/$(DOCS)

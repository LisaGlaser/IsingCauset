AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4 --install

SUBDIRS = printcolor inc src

INCD = -I $(top_builddir)/inc -I $(top_builddir)/printcolor
LIBD = -L $(top_builddir)/printcolor -lprintcolor

.PHONY: link all dirs

all: dirs link

dirs:
	@mkdir -p $(top_builddir)/bin

link:
	$(LINK_STR) $(INCD) $(LIBD)

cleanbin:
	rm -rf $(top_builddir)/bin

install-exec-hook:
	@mkdir -p $(prefix)/bin
	cp $(top_builddir)/bin/$(BINARY_NAME)* $(prefix)/bin

uninstall-hook:
	rm $(prefix)/bin/$(BINARY_NAME)*

EXTRA_DIST = clownfish lucy rdkit

all-local:
		$(MAKE) -C ../build/clownfish/runtime/c
		$(MAKE) -C ../build/lucy/c
		$(MAKE) -C ../build/rdkit/Code/Geometry
		$(MAKE) -C ../build/rdkit/Code/GraphMol
		$(MAKE) -C ../build/rdkit/Code/RDGeneral		

clean-local:
		rm -rf ../build/clownfish
		rm -rf ../build/lucy
		rm -rf ../build/rdkit

install-exec-local:
		mkdir -p $(DESTDIR)/$(POSTGRESQL_LIBDIR)
		cp ../build/clownfish/runtime/c/libclownfish.so* $(DESTDIR)/$(POSTGRESQL_LIBDIR)
		cp ../build/lucy/c/liblucy.so* $(DESTDIR)/$(POSTGRESQL_LIBDIR)

uninstall-local:
		rm $(DESTDIR)/$(POSTGRESQL_LIBDIR)/libclownfish.so*
		rm $(DESTDIR)/$(POSTGRESQL_LIBDIR)/liblucy.so*

EXTRA_DIST = clownfish lucy

../build/clownfish:
		cp -r $(srcdir)/clownfish ../build/clownfish
		chmod -R u+w ../build/clownfish
		cd ../build/clownfish/runtime/c; ./configure
		find ../build/clownfish -name Makefile -exec sed -i 's|-O2|-O2 -march=native|' '{}' ';'
		
../build/lucy:
		cp -r $(srcdir)/lucy ../build/lucy
		chmod -R u+w ../build/lucy
		cd ../build/lucy/c; ./configure
		sed -i 's|cfc|../../clownfish/compiler/c/cfc --include=../../clownfish/runtime/c/autogen/share/clownfish/include|' ../build/lucy/c/Makefile
		sed -i 's|-l clownfish|-L ../../clownfish/runtime/c -l clownfish|' ../build/lucy/c/Makefile
		find ../build/lucy -name Makefile -exec sed -i 's|-O2|-O2 -march=native|' '{}' ';'

all-local: ../build/clownfish ../build/lucy
		$(MAKE) -C ../build/clownfish/runtime/c
		$(MAKE) -C ../build/lucy/c

clean-local:
		rm -rf ../build/clownfish
		rm -rf ../build/lucy

install-exec-local:
		mkdir -p $(DESTDIR)/$(POSTGRESQL_LIBDIR)
		cp ../build/clownfish/runtime/c/libclownfish.so* $(DESTDIR)/$(POSTGRESQL_LIBDIR)
		cp ../build/lucy/c/liblucy.so* $(DESTDIR)/$(POSTGRESQL_LIBDIR)

uninstall-local:
		rm $(DESTDIR)/$(POSTGRESQL_LIBDIR)/libclownfish.so*
		rm $(DESTDIR)/$(POSTGRESQL_LIBDIR)/liblucy.so*

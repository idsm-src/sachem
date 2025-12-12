noinst_LTLIBRARIES = libinchi.la

libinchi_la_CFLAGS = -std=gnu99 -O3 -flto=jobserver -Wall -Wno-address -Wno-comment -Wno-int-in-bool-context -Wno-misleading-indentation -Wno-parentheses -Wno-unknown-pragmas -Wno-unused-but-set-variable -Wno-unused-function -Wno-unused-value -Wno-unused-variable -g
libinchi_la_CPPFLAGS = -I$(srcdir)/inchi/INCHI-1-SRC/INCHI_BASE/src -I$(srcdir)/inchi/INCHI-1-SRC/INCHI_API/libinchi/src -I$(srcdir)/inchi/INCHI-1-SRC/INCHI_API/libinchi/src/ixa -DTARGET_API_LIB -DCOMPILE_ANSI_ONLY

libinchi_la_SOURCES = \
        inchi/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichicano.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichicans.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichican2.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichierr.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichiisot.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichimake.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichimap4.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichiread.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichiring.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichirvr5.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichirvr6.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichisort.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichister.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/inchi_gui.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/readinch.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/runichi.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/runichi2.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/runichi3.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/runichi4.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/sha2.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/strutil.c \
        inchi/INCHI-1-SRC/INCHI_BASE/src/util.c \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_main.c \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_mol.c \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c

EXTRA_DIST = \
        inchi/INCHI-1-SRC/INCHI_BASE/src/bcf_s.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/extr_ct.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichicano.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichicant.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichicomn.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichicomp.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichidrp.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichierr.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichi.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichi_io.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichimain.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichimake.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichinorm.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichiring.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichirvrs.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichisize.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichister.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichitaut.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ichitime.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/incomdef.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/inchi_api.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/inchicmp.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/inchi_gui.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/inpdef.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/ixa.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/logging.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/mode.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/mol_fmt.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/readinch.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/sha2.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/strutil.h \
        inchi/INCHI-1-SRC/INCHI_BASE/src/util.h \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.h \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.h \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.h \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_main.h \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.h \
        inchi/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.h

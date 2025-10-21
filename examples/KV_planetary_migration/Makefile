export OPENGL=0
include ../../src/Makefile.defs


test: librebound
	@echo ""
	@echo "Compiling reproducible-migration-problem file ..."
	$(CC) -I../../src/ -Wl,-rpath,./ $(OPT) $(PREDEF) test-reproducible-migration-problem.c -L. -lrebound $(LIB) -o rebound
	@echo ""
	@echo "REBOUND compiled successfully."




all: librebound
	@echo ""
	@echo "Compiling reproducible-migration-problem file ..."
	$(CC) -I../../src/ -Wl,-rpath,./ $(OPT) $(PREDEF) reproducible-migration-problem.c -L. -lrebound $(LIB) -o rebound
	@echo ""
	@echo "REBOUND compiled successfully."


read: librebound
	@echo ""
	@echo "Compiling read-sim-archives file ..."
	$(CC) -I../../src/ -Wl,-rpath,./ $(OPT) $(PREDEF) read-sim-archives.c -L. -lrebound $(LIB) -o readrebound
	@echo ""
	@echo "REBOUND compiled successfully."


librebound: 
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C ../../src/
	@-rm -f librebound.so
	@ln -s ../../src/librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C ../../src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
	@-rm -vf readrebound

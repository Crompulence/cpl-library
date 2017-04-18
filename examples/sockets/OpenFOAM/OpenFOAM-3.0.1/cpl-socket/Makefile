all:
	@wmake libso src/CPLSocketFOAM
	@wmake src/solvers/CPLIcoFoam
clean:
	@wclean src/CPLSocketFOAM
	@wclean src/solvers/CPLIcoFoam
	@rm -r lib
	@rm -r bin
patch-scotch:                                                                                                                                                                                                                                    
	patch $(FOAM_SRC)/parallel/decompose/ptscotchDecomp/ptscotchDecomp.C ./config/ptscotchDecomp.patch

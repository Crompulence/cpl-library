all:
	@wmake libso src/CPLSocketFOAM
	@wmake src/solvers/CPLIcoFoam
	@wmake src/solvers/CPLporousIcoFoam
	@wmake src/solvers/CPLSediFoam
clean:
	@wclean src/CPLSocketFOAM
	@wclean src/solvers/CPLIcoFoam
	@wclean src/solvers/CPLporousIcoFoam
	@wclean src/solvers/CPLSediFoam
	@rm -r lib
	@rm -r bin
patch-scotch:                                                                                                                                                                                                                                    
	patch $(FOAM_SRC)/parallel/decompose/ptscotchDecomp/ptscotchDecomp.C ./config/ptscotchDecomp.patch

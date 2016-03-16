all:
	cd cpl-socket/ && ./Allwmake

clean:
	cd cpl-socket/ && ./Allwclean

patch-scotch:                                                                                                                                                                                                                                    
	patch OpenFOAM-$(FOAM_CPL_VERSION)/src/parallel/decompose/ptscotchDecomp/ptscotchDecomp.C ./cpl-socket/config/ptscotchDecomp.patch

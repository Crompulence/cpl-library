# Include platform if defined in default file "PLATFORM"

ifdef platform
    $(error Must define uppercase PLATFORM and not platform)
endif

ifndef PLATFORM
    PLATFORM= $(shell cat $(MAKEINCPATH)/PLATFORM)
endif
ifndef PLATFORM
# Error if not defined still
    platform_files= $(shell ls $(MAKEINCPATH) | grep  '.inc' \
                     | sed 's/.inc//')
    $(error The PLATFORM variable must be specified. \
            Try one of the following "PLATFORM= $(platform_files)" )

else
    # Include platform inc and create default file
    include $(MAKEINCPATH)/$(PLATFORM).inc
    $(shell echo $(PLATFORM) > $(MAKEINCPATH)/PLATFORM)
endif

# -- Some libraries and applications require freetype. They will definie this variable
ifeq ($(FREETYPE), 1)
ifeq ($(WIN32), 1)
    FT_CPPFLAGS:=-I/c/unx/3rdparty/include -I/c/unx/3rdparty/include/freetype2
    FT_LINK:=-L/c/unx/3rdparty/lib -lfreetype
else
    FT_CPPFLAGS:=$(shell freetype-config --cflags)
    FT_LINK:=$(shell freetype-config --libs)
endif
endif


LIB_LINKS:=-lsqlite3 -ldl -lpthread

EXTERNAL_LINKS+=$(FT_LINK) $(LIB_LINKS)
EXT_INCLUDES+=$(FT_CPPFLAGS) 

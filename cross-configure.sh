#!/bin/sh

CONFIG_SHELL=/bin/sh
export CONFIG_SHELL
PREFIX=/home/mihal/gnuwin32/cross-tools
TARGET=i386-mingw32msvc
PATH="$PREFIX/bin:$PREFIX/$TARGET/bin:$PATH"
export PATH
if [ -f "$PREFIX/$TARGET/bin/$TARGET-sdl-config" ]; then
    SDL_CONFIG="$PREFIX/$TARGET/bin/$TARGET-sdl-config"
    export SDL_CONFIG
fi
#cache=cross-config.cache
#sh configure --cache-file="$cache" \
#	--target=$TARGET --host=$TARGET --build=i386-linux \
#	$*
export CC=i386-mingw32msvc-gcc
export CXX=i386-mingw32msvc-g++
export LDFLAGS="-L/home/mihal/gnuwin32/cross-tools/lib -static -static-libgcc "
export CXXFLAGS="-Wno-redundant-decls -D__GW32C__ -D__GW32__ -D_LARGEFILE_SOURCE=1 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -I/home/mihal/gnuwin32/cross-tools/include -I../include -I../../include -I/usr/local/include -idirafter /home/mihal/gnuwin32/cross-tools/include/glibc"
export CPPFLAGS="-I/home/mihal/gnuwin32/cross-tools/include -I../include -I../../include -I/usr/local/include"
export LIBS="-ljpeg -lz"
sh configure --target=$TARGET --host=i386-linux $*
status=$?
#rm -f "$cache"
exit $status

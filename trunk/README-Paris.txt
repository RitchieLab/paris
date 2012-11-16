PARIS Build Instructions

PARIS uses a SQLite database and depends on having that installed.  If it isn't present on
your system you will need to install it.

PARIS itself uses make to compile the executable.  To compile, run 'make' in the main 
directory of this package.

The executable will be bin/paris64 or bin/paris depending on if your system is 32 or 64bit linux. 

The software should compile on windows and OS X and any other system for which GCC is available. Windows users
will need to set up MinGW in order to compile the software (www.mingw.org). OS X users need only to ensure that
they have GCC installed (install the developer tools from the installation media or the App Store).

When compiling for windows, users should use WIN32=1 when building paris (i.e. make WIN32=1 )

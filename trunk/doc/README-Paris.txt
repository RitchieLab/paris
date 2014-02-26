PARIS README
2/14/14

1.  Building
2.  Running
3.  Dependencies
4.  Troubleshooting


1.  Building PARIS

Unpack paris.1.1.0.tgz and then run the following commands:

cd paris
make

The PARIS executable will be placed in the bin directory with an extension appropriate to 
the system where built.  For example, on a 64 bit Linux system, it will be paris64 and
on a Mac it will be paris-OSX. 

The software should compile on windows and OS X and any other system for which GCC is available. Windows users will need to set up MinGW in order to compile the software (www.mingw.org). OS X users need only to ensure that they have GCC installed (install the developer tools from the installation media or the App Store).

When compiling for windows, users should use WIN32=1 when building paris (i.e. make WIN32=1 )

2.  Running

PARIS uses a configuration file to control the behavior of the software and provide PARIS
with necessary information to run.  The simplest way to run PARIS is to pass the name of
the configuration file to the program as below:

paris64 configuration.txt

It is possible to generate a template configuration
file by using the --sample-config option:

paris64 --sample-config > configuration.txt

The resulting file can be edited for your own use.  Details on the various options and 
required input files can be found in the paris-reference.pdf file distributed with the
software.


2.  Dependencies

PARIS utilizes the boost and SQLite libraries.  The underlying data is stored in a SQLite
database distributed with the code.  IF PARIS fails to compile, one or both of these may 
need to be installed on your system. (see below)


3.  Troubleshooting

If required libraries are missing or not found, you may see error messages similar to the 
ones below.

  a.  resultsdatabase.h:6:21: error: sqlite3.h: No such file or directory
	
      The SQLite header file has not been found.  If SQLite has not been installed, it 
      needs to be.  If it has been, then if it has been installed in a nonstandard
      location, that location must be specified in the extlibs.make file at the
      EXT_INCLUDES line.  For example if it was installed to /opt/local:
			
        EXT_INCLUDES+=$(FT_CPPFLAGS) /opt/local/include
				
      Another possible cause can be that the development SQLite package needs to be 
      installed.  Some distributions split the development library package from the
      runtime.  For example on Ubuntu you would need to install the libsqlite3-dev 
      package.			
			
  b.  types.h:24:36: error: boost/dynamic_bitset.hpp: No such file or directory
	
      In this case, the required boost library header has not been found.  Install the 
      library or if it is in a nonstandard location, specify that in the src/utility 
      makefile.  The information should be appended to the top line of the file as
      below (if it was installed to /opt/local):
			
        PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI) -I/opt/local/include


  c.  knowledgedatabase.cpp:(.text+0x7ad): undefined reference to `sqlite3_prepare_v2'

      In this case, the SQLite library has not been found.  If it has been installed in a 
      nonstandard location, that can be specified in the file extlibs.make.  For example,
      if it was installed to /opt/local, append that information to the LIB_LINKS line:
			
        LIB_LINKS:=-lsqlite3 -ldl -lpthread -L/opt/local/lib
			
			

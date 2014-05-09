/* 
 * File:   DBExcept.h
 * Author: Scott Dudek
 *
 */

#ifndef _DBEXCEPT_H
#define	_DBEXCEPT_H

#include <exception>
#include <string>
using namespace std;

///
/// Class thrown for exception in database system <br>
/// Error messages are set by the creating class
///

/// Exception for database
class DBExcept: public exception{
        public:
          DBExcept() throw();
          DBExcept(string message){error=message;}
          ~DBExcept()throw(){};
          virtual const char* what() const throw(){return error.c_str();}

        private:
          string error;
};

#endif	/* _DBEXCEPT_H */


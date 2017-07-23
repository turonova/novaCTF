#include <iostream>
#include <exception>
#include <string>
#include <sstream>
#include <stdlib.h>

using namespace std;

class ExceptionFileOpen : public exception
{
public:

	ExceptionFileOpen(std::string fileName)
	{
		message = "Following file could not be opened: " + fileName + "\n";
	}

	~ExceptionFileOpen() throw() {}

	const char * what () const throw ()
	{
		return message.c_str();
	}

	std::string message;

};

class ExceptionFileFormat : public exception
{
public:

	ExceptionFileFormat()
	{
		stringstream mess;
		mess << "The file format specified in the header is not supported!" << endl;
		message = mess.str();
	}

	~ExceptionFileFormat() throw() {}

	const char * what () const throw ()
	{
		return message.c_str();
	}

	std::string message;

};

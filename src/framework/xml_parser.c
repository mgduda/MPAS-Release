#include <stdio.h>
#include "ezxml.h"

void parse_xml_streams()
{
        ezxml_t streams;
	size_t bufsize;
        char *streambuf;


	bufsize = (size_t)1024;
	streambuf = (char *)malloc(sizeof(char) * bufsize);

        streams = ezxml_parse_str(streambuf, bufsize);
	printf("Parsed XML streams file\n");

	free(streambuf);
}

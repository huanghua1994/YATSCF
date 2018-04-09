#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#include "utils.h"

static char* sys_time_str()
{
	static char buffer[30], usec_buf[3];
	time_t timer;
	struct tm* tm_info;
	struct timeval tmnow;
	
	gettimeofday(&tmnow, NULL);
	tm_info = localtime(&tmnow.tv_sec);
	strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);
	strcat(buffer, ".");
	sprintf(usec_buf, "%d", (int)tmnow.tv_usec);
	strcat(buffer, usec_buf);
	buffer[23] = 0;
	
	/*
	// If you don't have <sys/time.h>, use the following method
	time(&timer);
	tm_info = localtime(&timer);
	strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);
	*/
	
	return (&buffer[0]);
}

void print_msg(const int msg_type, const char *msg)
{
	printf("[%s][%s]\t%s\n", sys_time_str(), msg_type_str[msg_type], msg);
	if (msg_type == ERRORMSG) exit(255);
}

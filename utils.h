#ifndef _YATSCF_UTILS_H_
#define _YATSCF_UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

#define INFOMSG  0
#define ERRORMSG 1
#define DEBUGMSG 2

const char *msg_type_str[3] = {"INFO", "ERROR", "DEBUG"};

void print_msg(const int msg_type, const char *msg);

#ifdef __cplusplus
}
#endif

#endif

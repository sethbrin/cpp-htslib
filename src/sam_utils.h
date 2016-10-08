#ifdef SAM_UTILS_H_
#define SAM_UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

#if defined __GNUC__ && __GNUC__ >= 2
#define CHECK_PRINTF(fmt,args) __attribute__ ((format (printf, fmt, args)))
#else
#define CHECK_PRINTF(fmt,args)
#endif

void print_error(const char *subcommand, const char *format, ...) CHECK_PRINTF(2, 3);
void print_error_errno(const char *subcommand, const char *format, ...) CHECK_PRINTF(2, 3);

#ifdef __cplusplus
}
#endif

#endif

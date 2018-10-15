#ifndef UTILS_H
#define UTILS_H

#ifndef MAXLENGTH
#define MAXLENGTH 1024
#endif

int file_exist(char *filename);
int directory_exist(char *dirname);
void *get_directory(char *dirname);

#endif

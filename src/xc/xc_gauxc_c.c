#include <stdio.h>
#include <stddef.h>

void print_c_string(const char* c_message) {
    if (c_message != NULL) {
        printf("%s", c_message);
    }
}
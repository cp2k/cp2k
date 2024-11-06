#include <openPMD.h>

int main()
{
    openPMD_Series series;
    openPMD_Series_create("test.json", openPMD_Access_create, &series);
    openPMD_Series_close(series);
}

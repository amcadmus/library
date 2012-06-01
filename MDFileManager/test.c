#include <stdio.h>

int main(int argc, char * argv[])
{
  FILE * fp = fopen ("tmp", "w");
  fprintf (fp, "%5d%5s%5s%5d\n", 10000, "lala", "zizi", 9999);
  fprintf (fp, "%5d%5s%5s%5d\n", 10000, "lala", "zizi", 10000);
  fclose (fp);
}

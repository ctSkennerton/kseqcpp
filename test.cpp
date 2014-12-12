#include "kseq.hpp"
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

int main(int argc, char * argv[]) 
{
    //gzFile fp = gzopen(argv[1], "r");
    //FunctorZlib gzr;
    //kstream<gzFile, FunctorZlib> ks(fp, gzr);
    int fp = open(argv[1], O_RDONLY);
    FunctorRead r;
    kstream<int, FunctorRead> ks(fp, r);
    kseq seq;

    int l = 0;
    int c = 0;
    while((l = ks.read(seq)) >= 0) {
        std::cout << seq.name.s << std::endl;
        std::cout << seq.seq.s << std::endl;
        c++;
    }
    std::cout << c << std::endl;
    //gzclose(fp);
    close(fp);

    return 0;
}

#include <iostream>

#include "htslib/sam.h"

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <boost/lexical_cast.hpp>

#include "easehts.h"

using namespace std;

int sam_read() {
    char filename[] = "/home/zp/work/rnaSeq/easehtslib/data/xx#rg.sam";
    samFile *in = sam_open(filename, "r");
    if (in == NULL) {
        fprintf(stderr, "Error opening \"%s\"\n", filename);
        return EXIT_FAILURE;
    }

    ncic::easehts::SAMBAMTextReader reader(in);


    return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
    cout << stoi("22");
    cout << "shit" << endl;
    sam_read();

    return 0;
}

#include "mutect.h"

int main(int argc, char** argv) {
  ncic::mutect::Mutect mutect(argc, argv);

  mutect.Run();

  return 0;
}

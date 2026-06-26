#include <iostream>

void combine_offset_gain_timeevo()
{
    std::cout << "IMP HistAll.root files store flat TH2D matrices named hetg0-hetg11 "
                 "and hetl0-hetl25, so the old ILL run-range combiner is not used.\n"
              << "Run solveTimeEvo_IMP directly on the selected histogram; its output "
                 "is already named HistAll_<histogram>_correctedTimeEvo.root.\n";
}

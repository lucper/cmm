#ifndef H_UTILS
#define H_UTILS

inline void print_progress(std::size_t done, std::size_t total)
{
    if (total == 0) return;

    const int bar_width = 40;
    double frac = (double) done / (double) total;
    if (frac > 1.0) frac = 1.0;

    int filled = (int) (frac * bar_width);

    std::cerr << "\r[";
    for (int i = 0; i < bar_width; i++) std::cerr << (i < filled ? '#' : ' ');
    std::cerr << "] " << std::setw(3) << (int) (frac * 100.0) << "% "
            << "(" << done << "/" << total << ")"
            << std::flush;

    if (done == total) std::cerr << "\n";
}

#define DBG(msg) do { std::cerr << "DEBUG: " << msg << "\n"; } while(0)

#endif

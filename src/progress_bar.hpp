#ifndef H_PROGRESS_BAR
#define H_PROGRESS_BAR

#include <atomic>
#include <cstddef>
#include <memory>

namespace indicators {
class ProgressBar;
}

/* Progress bar of the wildcard configurations of a search, drawn on stderr if it is a terminal (in
 * a file, every redraw would remain as a line). It is drawn at 0 when created, and redrawn about
 * every 1% of the configurations. The threads of the search may call advance at the same time.
 *
 * Invariant: the count drawn only increases, and is at most the configurations counted so far. */
class progress_bar_t {
public:
    progress_bar_t(size_t total);
    ~progress_bar_t();

    progress_bar_t(const progress_bar_t&) = delete;
    progress_bar_t& operator=(const progress_bar_t&) = delete;

    /* Counts one more searched configuration, and redraws the bar if it reached the next 1%, or the
     * last configuration. */
    void advance();

private:
    size_t total;
    size_t update_interval;
    std::atomic<size_t> completed;
    size_t drawn;  // the highest count drawn
    std::unique_ptr<indicators::ProgressBar> bar;
};

#endif

#include "progress_bar.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include "indicators.hpp"

using namespace indicators;

progress_bar_t::progress_bar_t(size_t total)
    : total(total), update_interval(std::max(static_cast<size_t>(1), total / 100)), completed(0),
      drawn(0),
      bar(std::make_unique<ProgressBar>(
          option::BarWidth(50),
          option::Start("["),
          option::Fill("#"),
          option::Lead("#"),
          option::Remainder(" "),
          option::End("]"),
          option::ForegroundColor(Color::white),
          option::ShowPercentage(true),
          option::ShowElapsedTime(true),
          option::PrefixText("Mining motifs "),
          option::Stream(std::cerr)
      ))
{
    bar->set_option(option::PostfixText("[wildcard combinations: 0/" + std::to_string(total) + "]"));
    bar->set_progress(0);
    std::cerr << "\r" << std::flush;
}

// Defined here, where ProgressBar is complete, so that the unique_ptr can delete it.
progress_bar_t::~progress_bar_t() = default;

void progress_bar_t::advance()
{
    size_t current = ++completed;
    if (current % update_interval == 0 || current == total) {
        #pragma omp critical
        {
            // Threads may get here in another order than their counts; a lower count than the one
            // drawn would move the bar back.
            if (current > drawn) {
                drawn = current;
                bar->set_option(option::PostfixText("[wildcard combinations: " + std::to_string(current) +
                                                    "/" + std::to_string(total) + "]"));
                bar->set_progress((static_cast<float>(current) / total) * 100.0f);
                std::cerr << "\r" << std::flush;
            }
        }
    }
}

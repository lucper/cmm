#ifndef H_UTILS
#define H_UTILS

#include <vector>

template <typename T, typename key_extractor>
void counting_sort(std::vector<T>& data, INT max_val,
                   std::vector<T>& temp, std::vector<INT>& count,
                   key_extractor get_key)
{
    std::fill(count.begin(), count.begin() + max_val + 1, 0);
    for (int i = 0; i < data.size(); i++) count[get_key(data[i])]++;
    for (int i = 1; i < max_val + 1; i++) count[i] += count[i - 1];

    for (int i = data.size() - 1; i >= 0; i--)
        temp[--count[get_key(data[i])]] = data[i];

    data = temp;
}

#define DBG(msg) do { std::cerr << "DEBUG: " << msg << "\n"; } while(0)

#endif

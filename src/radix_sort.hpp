#ifndef H_RADIX_SORT
#define H_RADIX_SORT

/* Generic radix sort to sort a vector of keys with an associated
 * vector of "payloads" (satellite data). Note that keys should be
 * a vector of uint32_t or uint64_t.
 */
template <typename KeyType = uint64_t, typename PayloadType = int>
void radix_sort(std::vector<KeyType>& keys, 
                std::vector<KeyType>& key_buffer,
                std::vector<PayloadType>* payload = nullptr,
                std::vector<PayloadType>* payload_buffer = nullptr)
{
    if (keys.empty()) return;

    if (key_buffer.size() < keys.size()) key_buffer.resize(keys.size());

    const int bits_per_pass = 16;
    const int bins = 1 << bits_per_pass;
    const int num_passes = (sizeof(KeyType) * 8) / bits_per_pass;

    KeyType* src_key = keys.data();
    KeyType* dst_key = key_buffer.data();
    PayloadType* src_pay = (payload) ? payload->data() : nullptr;
    PayloadType* dst_pay = (payload_buffer) ? payload_buffer->data() : nullptr;

    size_t counts[bins];

    for (size_t p = 0; p < num_passes; p++) {
        std::fill(counts, counts + bins, 0);
        size_t shift = p * bits_per_pass;

        for (size_t i = 0; i < keys.size(); i++)
            counts[(src_key[i] >> shift) & 0xFFFF]++;

        size_t pos = 0;
        for (size_t i = 0; i < bins; i++) {
            size_t c = counts[i];
            counts[i] = pos;
            pos += c;
        }

        for (size_t i = 0; i < keys.size(); i++) {
            uint32_t bucket = (src_key[i] >> shift) & 0xFFFF;
            uint32_t target = counts[bucket]++;
            dst_key[target] = src_key[i];
            if (dst_pay && src_pay) dst_pay[target] = src_pay[i];
        }

        std::swap(src_key, dst_key);
        std::swap(src_pay, dst_pay);
    }

    if (src_key != keys.data()) {
        std::copy(key_buffer.begin(), key_buffer.end(), keys.begin());
        if (payload && payload_buffer)
            std::copy(payload_buffer->begin(), payload_buffer->end(), payload->begin());
    }
}

#endif

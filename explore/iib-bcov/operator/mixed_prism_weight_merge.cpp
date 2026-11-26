// Add sorted projective coefficient banks, holding one rational per input stream.
#include <gmpxx.h>
#include <zlib.h>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

class Bank {
    gzFile file;
    std::string line;
public:
    std::uint64_t key = 0;
    mpq_class weight;
    explicit Bank(const char *path) : file(gzopen(path,"rb")) {
        if (!file) throw std::runtime_error(std::string("Cannot read ") + path);
        gzbuffer(file, 65536);
    }
    ~Bank() { gzclose(file); }
    bool next() {
        line.clear();
        char part[256];
        while (gzgets(file,part,sizeof part)) {
            line += part;
            if (line.back() == '\n') break;
        }
        if (line.empty() && gzeof(file)) return false;
        if (line.empty() || line.back() != '\n') throw std::runtime_error("Truncated weight bank");
        line.pop_back();
        std::uint64_t new_key;
        const char *end = line.data() + line.size();
        auto parsed = std::from_chars(line.data(),end,new_key);
        if (parsed.ec != std::errc() || parsed.ptr == end || *parsed.ptr != ' ' ||
            new_key <= key || new_key >= (std::uint64_t(1)<<36))
            throw std::runtime_error("Invalid or unsorted projective power");
        if (weight.set_str(parsed.ptr+1,10) || weight.get_den() <= 0)
            throw std::runtime_error("Invalid projective weight");
        weight.canonicalize();
        key = new_key;
        return true;
    }
};

int main(int argc, char **argv) {
    if (argc<3) return 2;
    try {
        std::vector<std::unique_ptr<Bank>> banks;
        auto later=[](Bank *first,Bank *second){return first->key>second->key;};
        std::priority_queue<Bank *,std::vector<Bank *>,decltype(later)> queue(later);
        for (int i=2;i<argc;++i) {
            banks.push_back(std::make_unique<Bank>(argv[i]));
            if (banks.back()->next()) queue.push(banks.back().get());
        }
        gzFile out=gzopen(argv[1],"wb1");
        if (!out) throw std::runtime_error("Cannot write merged weights");
        gzbuffer(out,65536);
        std::uint64_t count=0;
        while (!queue.empty()) {
            const auto key=queue.top()->key;
            mpq_class sum;
            do {
                auto *bank=queue.top();
                queue.pop();
                sum+=bank->weight;
                if (bank->next()) queue.push(bank);
            } while (!queue.empty() && queue.top()->key==key);
            if (sum==0) continue;
            const auto line=std::to_string(key)+" "+sum.get_str()+"\n";
            if (gzwrite(out,line.data(),line.size())!=static_cast<int>(line.size()))
                throw std::runtime_error("Cannot write merged coefficient");
            ++count;
        }
        if (gzclose(out)!=Z_OK) throw std::runtime_error("Cannot finish merged weights");
        std::printf("%llu\n",static_cast<unsigned long long>(count));
    } catch (const std::exception &error) {
        std::fprintf(stderr,"%s\n",error.what());
        return 1;
    }
}

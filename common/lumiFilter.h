#include <unordered_map>
#include <memory>
#include <algorithm>
#include <stdexcept>

class JsonHelper {
 public:
  using pair_t = std::pair<unsigned int, unsigned int>;
  using jsonmap_t = std::unordered_map<unsigned int, std::vector<pair_t>>;

 JsonHelper(const std::vector<unsigned int> &runs, const std::vector<unsigned int> &firstlumis, const std::vector<unsigned int> &lastlumis) :
  jsonmap_(std::make_shared<jsonmap_t>()){
    for (unsigned int i = 0; i < runs.size(); i++){
      (*jsonmap_)[runs[i]].push_back(std::make_pair(firstlumis[i], lastlumis[i]));
    }
    for (auto &item: *jsonmap_){
      std::sort(item.second.begin(), item.second.end());
    }
  }

  bool operator () (unsigned int run, unsigned int lumi) const{
    if (run == 1)
      return true;

    const auto it = jsonmap_->find(run);
    if (it != jsonmap_->end()){
      const auto &pairs = it->second;
      auto const pairit = std::lower_bound(pairs.begin(), pairs.end(), lumi, [](const pair_t &pair, unsigned int val) {return pair.second < val;});
      if (pairit != pairs.end()){
	if (lumi >= pairit->first){
	  return true;
	}
      }
    }
    return false;
  }

 private:
  std::shared_ptr<jsonmap_t> jsonmap_;
};

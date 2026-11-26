// The balanced function test as a six-cochain on unmarked transverse vectors.
#define MIXED_SIX_COMPONENTS
#include "mixed_six_probe.cpp"

template<int Q> vector<Waves<Q>> fixtures(int count) {
  vector<Waves<Q>> out;
  exact = true;
  for (int i = 0; i < count; ++i) out.push_back(nonzero_waves<Q>());
  return out;
}

int main() {
  exact = true;
  vector<Waves<5>> five;
  vector<Waves<6>> six;
  for (int i = 0; i < 110; ++i) {
    five.push_back(nonzero_waves<5>());
    six.push_back(nonzero_waves<6>());
  }
  auto seven = fixtures<7>(45);
  exact = false;
  auto marked_five = basis(five,{75,51,21,3},true);
  auto marked_six = basis(six,{70,28,7},true);
  exact = true;
  vector<Wide> degrees;
  for (const auto &c : marked_six) {
    int degree = 0;
    for (int k : c.front().det) {
      for (int j = 0; j < 3; ++j) {
        degree += k % 8 == 6;
        k /= 8;
      }
    }
    for (int k : c.front().pair) degree += k % 8 == 6;
    degrees.push_back(degree);
  }
  line(0, degrees);
  for (auto &w : seven) line(1, differential(marked_six,w));

  five = fixtures<5>(100);
  six = fixtures<6>(135);
  seven = fixtures<7>(95);
  exact = false;
  auto unmarked_five = basis(five,{53},true,0);
  auto unmarked_six = basis(six,{74},true,0);
  exact = true;
  for (auto &w : six) {
    auto r = differential(unmarked_five,w);
    Evaluations<6> tests(w);
    for (auto &c : unmarked_six) r.push_back(integer_value(tests,c,false));
    line(2,r);
  }
  for (auto &w : seven) line(3,differential(unmarked_six,w));
  for (int mask = 0; mask < 64; ++mask) {
    auto w = conserving_fields(mask,false);
    auto r = differential(unmarked_five,w);
    auto d = differential(marked_five,w);
    r.insert(r.end(),d.begin(),d.end());
    Evaluations<6> tests(w);
    for (auto &c : marked_six) r.push_back(integer_value(tests,c,false));
    line(4,r);
  }
}

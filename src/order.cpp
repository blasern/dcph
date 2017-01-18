// From http://stackoverflow.com/a/41704272/2591234
// Thanks Nathan Russell
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <tuple>

template <typename... Ts>
struct Point {
  std::size_t index;
  
  using tuple_t = std::tuple<Ts...>;
  tuple_t data;
  
  template <typename... Vs>
  Point(std::size_t i, const Vs&... vs)
    : index(i),
      data(std::forward<decltype(vs[i])>(vs[i])...)
  {}
  
  bool operator<(const Point& other) const {
    return data < other.data;
  }
};

template <typename... Vecs>
class PointVector {
public:
  std::size_t sz;
  using point_t = Point<typename std::remove_reference<decltype(Vecs{}[0])>::type...>;
  using vector_t = std::vector<point_t>;
  vector_t data;
  
  template <typename... Vs>
  PointVector(const Vecs&... vs)
    : sz(min_size(vs...))
  {
    data.reserve(sz);
    for (std::size_t i = 0; i < sz; i++) {
      data.emplace_back(i, vs...);
    }
  }
  
  Rcpp::IntegerVector sorted_index() const {
    vector_t tmp(data);
    std::stable_sort(tmp.begin(), tmp.end());
    
    Rcpp::IntegerVector res(sz);
    for (std::size_t i = 0; i < sz; i++) {
      res[i] = tmp[i].index + 1;
    }
    
    return res;
  }
  
private:
  template <typename V>
  std::size_t min_size(const V& v) {
    return v.size();
  }
  
  template <typename T, typename S, typename... Vs>
  std::size_t min_size(const T& t, const S& s, const Vs&... vs) {
    return t.size() < s.size() ?
    min_size(t, vs...) :
    min_size(s, vs...);
  }
};

template <typename... Vecs>
PointVector<Vecs...> MakePointVector(const Vecs&... vecs) {
  return PointVector<Vecs...>(vecs...);
}

Rcpp::IntegerVector order2(Rcpp::NumericVector x, Rcpp::NumericVector y) {
  auto pv = MakePointVector(x, y);
  return pv.sorted_index();
}

Rcpp::IntegerVector order3(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::NumericVector z) {
  auto pv = MakePointVector(x, y, z);
  return pv.sorted_index();
}

Rcpp::IntegerVector order4(Rcpp::NumericVector x1, Rcpp::NumericVector x2, Rcpp::NumericVector x3, Rcpp::NumericVector x4){
  auto pv = MakePointVector(x1, x2, x3, x4);
  return pv.sorted_index();
}

Rcpp::IntegerVector order5(Rcpp::NumericVector x1, Rcpp::NumericVector x2, Rcpp::NumericVector x3, Rcpp::NumericVector x4, Rcpp::NumericVector x5){
  auto pv = MakePointVector(x1, x2, x3, x4, x5);
  return pv.sorted_index();
}

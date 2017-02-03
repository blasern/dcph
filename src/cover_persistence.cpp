// use Rcpp
#include <Rcpp.h>

// wrapper algorithm that computes the persistence pairs of a given boundary matrix using a specified algorithm
#include <phat/compute_persistence_pairs.h>

// main data structure (choice affects performance)
#include <phat/representations/vector_vector.h>

// algorithm (choice affects performance)
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/twist_reduction.h>

// filtration_value, dimension, vertices
typedef std::tuple<double, int, std::set<int>> filtration_tuple2;
typedef std::tuple<double, int, int, int, int> filtration_tuple;
// function to sort filtration
bool compare_filtration (const filtration_tuple &lhs, const filtration_tuple &rhs){
  if (std::get<0>(lhs) == std::get<0>(rhs)){
    return std::get<1>(lhs) < std::get<1>(rhs);
  }
  else{
    return std::get<0>(lhs) < std::get<0>(rhs);
  }
}
bool compare_filtration2 (const filtration_tuple2 &lhs, const filtration_tuple2 &rhs){
  if (std::get<0>(lhs) == std::get<0>(rhs)){
    return std::get<1>(lhs) < std::get<1>(rhs);
  }
  else{
    return std::get<0>(lhs) < std::get<0>(rhs);
  }
}

// function for finding if sets are disjoint
template<class Set1, class Set2>
bool is_disjoint(const Set1 &set1, const Set2 &set2)
{
  if(set1.empty() || set2.empty()) return true;
  
  typename Set1::const_iterator
    it1 = set1.begin(),
      it1End = set1.end();
  typename Set2::const_iterator
    it2 = set2.begin(),
      it2End = set2.end();
  
  if(*it1 > *set2.rbegin() || *it2 > *set1.rbegin()) return true;
  
  while(it1 != it1End && it2 != it2End)
  {
    if(*it1 == *it2) return false;
    if(*it1 < *it2) { it1++; }
    else { it2++; }
  }
  
  return true;
}

// combination subset
std::vector<std::set<int>> combinations(int n, int k){
  std::vector<std::set<int>> out; 
  std::vector<bool> v(n);
  std::fill(v.begin(), v.begin() + k, true);
  
  do {
    std::set<int> vec;
    for (int i = 0; i < n; ++i){
      if (v[i]) {
        vec.insert(i);
      }
    }
    out.push_back(vec);
  } while (std::prev_permutation(v.begin(), v.end()));
  return out;
}

//' Persistence from cover
//' 
//' Calculate the persistence from the cover
//' 
//' @param cover A cover
//' @param max_dim maximal dimension
// [[Rcpp::export]]
Rcpp::NumericMatrix persistence_from_cover(Rcpp::S4 cover, Rcpp::IntegerVector max_dim) {
  // check that input is a cover
  if (!cover.inherits("cover")) Rcpp::stop("Input must be a cover");
  
  // extract indices and diameter
  int maxdim = max_dim(0);
  Rcpp::List subsets = cover.slot("subsets");
  std::vector<std::set<int>> indices;
  std::vector<double> death_diameters; 
  for (int i = 0; i < subsets.length(); ++i){
    Rcpp::S4 subset = subsets(i);
    std::set<int> ind_set;
    Rcpp::IntegerVector indices_i = Rcpp::as<Rcpp::IntegerVector>(subset.slot("indices"));
    for (int j = 0; j < indices_i.size(); ++j){
      ind_set.insert(indices_i[j]);
    }
    indices.push_back(ind_set);
    death_diameters.push_back(subset.slot("death"));
  }
  
  // reindex
  // Rcpp::Rcout << "reindexing..." << std::endl;
  Rcpp::NumericMatrix distance_matrix = cover.slot("distance_matrix");
  int n = distance_matrix.nrow();
  std::vector<std::vector<int>> contained_in_vec(n);
  
  for (unsigned int j = 0; j < indices.size(); ++j){
    std::set<int> index = indices[j];
    for (auto ind : index){
      contained_in_vec[ind-1].push_back(j);
    }
  }
  std::set<std::vector<int>> contained_in;
  for (auto con : contained_in_vec){
    contained_in.insert(con);
  }
  // Rcpp::Rcout << "reindexed..." << std::endl;
  
  // calculate 2-sceleton 
  // Rcpp::Rcout << "Calculate 2-sceleton..." << std::endl;
  std::set<std::vector<int>> low_card_sets;
  for (int dim = 1; dim < maxdim + 3; ++dim){
    for (auto cont_set: contained_in){
      std::vector<std::set<int>> a = combinations(cont_set.size(), dim);
      for (unsigned int i = 0; i < a.size(); ++i){
        std::set<int> mset = a[i];
        std::vector<int> low_card_set;
        for(auto j : mset){
          int in_set = cont_set[j];
          low_card_set.push_back(in_set);
        }
        low_card_sets.insert(low_card_set);
      }
    }
  }
  // Rcpp::Rcout << "2-sceleton calculated" << std::endl;
  
  // save as filtration vector
  std::vector<filtration_tuple> filtration_vector;
  std::vector<filtration_tuple2> filtration_vector2;
  for (auto lcs : low_card_sets){
    // filtration_value, dimension, V1, V2, V3
    std::vector<int> vertices;
    std::set<int> vertices2;
    std::vector<double> values;
    for (auto lc : lcs){
      vertices.push_back(lc);
      vertices2.insert(lc);
      values.push_back(death_diameters[lc]);
      // Rcpp::Rcout << lc.first << " " << lc.second << "; ";
    }
    double diam = *max_element(values.begin(), values.end());
    int dim = vertices.size() - 1;
    for (int i = dim; i < maxdim + 1; ++i){
      vertices.push_back(vertices[dim]);
    }
    filtration_vector.push_back(std::make_tuple(diam, dim,  
                                                vertices[0], vertices[1], vertices[2]));
    filtration_vector2.push_back(std::make_tuple(diam, dim, vertices2));
    // Rcpp::Rcout << std::endl;
  }
  // Rcpp::Rcout << std::endl;

  // sort filtration_vector
  sort(filtration_vector.begin(), filtration_vector.end(), compare_filtration);
  sort(filtration_vector2.begin(), filtration_vector2.end(), compare_filtration2);
  // 
  // Rcpp::Rcout << filtration_vector.size() << std::endl;
  // Rcpp::Rcout << filtration_vector2.size() << std::endl;
  
  // save as map
  std::map<std::tuple<int, int, int>, std::tuple<int, double, int>> filtration;
  for(std::vector<filtration_tuple>::iterator it = filtration_vector.begin(); it != filtration_vector.end(); ++it) {
    filtration[std::make_tuple(std::get<2>(*it), std::get<3>(*it), std::get<4>(*it))] =
      std::make_tuple(std::distance(filtration_vector.begin(), it),
                      std::get<0>(*it), std::get<1>(*it));
  }
  // save as map
  std::map<std::set<int>, std::tuple<int, double, int>> filtration2;
  for(std::vector<filtration_tuple2>::iterator it = filtration_vector2.begin(); it != filtration_vector2.end(); ++it) {
    filtration2[std::get<2>(*it)] =
      std::make_tuple(std::distance(filtration_vector2.begin(), it),
                      std::get<0>(*it), std::get<1>(*it));
  }

  // define a boundary matrix
  phat::boundary_matrix< phat::vector_vector > boundary_matrix;
  
  // set the number of columns 
  boundary_matrix.set_num_cols(filtration2.size());
  
  // set the dimension of the cell that a column represents:
  for(unsigned int i = 0; i < filtration_vector2.size(); ++i) {
    boundary_matrix.set_dim(i, std::get<1>(filtration_vector2[i]) );
  }
  
  // Rcpp::Rcout << "Initializing boundary matrix..." << std::endl;
  // set the respective columns -- the columns entries have to be sorted
  std::vector< phat::index > temp_col;
  for(unsigned int i = 0; i < filtration_vector2.size(); ++i) {
    temp_col.clear();
    if (std::get<1>(filtration_vector2[i]) > 0){
      std::set<int> p_set = std::get<2>(filtration_vector2[i]);
      for (auto p : p_set){
        std::set<int> point_set = p_set;
        point_set.erase(p);
        temp_col.push_back(std::get<0>(filtration2[point_set]));
      }
    }
    std::sort(temp_col.begin(), temp_col.end());
    boundary_matrix.set_col(i, temp_col);
  }
  temp_col.clear();
  // Rcpp::Rcout << "Boundary matrix initialized" << std::endl;
  
  // // print some information of the boundary matrix:
  // Rcpp::Rcout << std::endl;
  // Rcpp::Rcout << "The boundary matrix has " << boundary_matrix.get_num_cols() << " columns: " << std::endl;
  // for( phat::index col_idx = 0; col_idx < boundary_matrix.get_num_cols(); col_idx++ ) {
  //   Rcpp::Rcout << "Colum " << col_idx << " represents a cell of dimension " << (int)boundary_matrix.get_dim( col_idx ) << ". ";
  //   if( !boundary_matrix.is_empty( col_idx ) ) {
  //     std::vector< phat::index > temp_col;
  //     boundary_matrix.get_col( col_idx, temp_col );
  //     Rcpp::Rcout << "Its boundary consists of the cells";
  //     for( phat::index idx = 0; idx < (phat::index)temp_col.size(); idx++ )
  //       Rcpp::Rcout << " " << temp_col[ idx ];
  //   }
  //   Rcpp::Rcout << std::endl;
  // }
  // Rcpp::Rcout << "Overall, the boundary matrix has " << boundary_matrix.get_num_entries() << " entries." << std::endl;

  // define the object to hold the resulting persistence pairs
  // Rcpp::Rcout << "Calculating persistent homology... " << std::endl;
  phat::persistence_pairs pairs;

  // choose an algorithm (choice affects performance) and compute the persistence pair
  // (modifies boundary_matrix)
  phat::compute_persistence_pairs< phat::twist_reduction >( pairs, boundary_matrix );

  // sort the persistence pairs by birth index
  pairs.sort();
  // Rcpp::Rcout << "Persistent homology calculated... " << std::endl;

  // // print the pairs:
  // Rcpp::Rcout << std::endl;
  // Rcpp::Rcout << "There are " << pairs.get_num_pairs() << " persistence pairs: " << std::endl;
  // for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
  //   Rcpp::Rcout << "Birth: " << pairs.get_pair( idx ).first << ", Death: " << pairs.get_pair( idx ).second << std::endl;

  // write to output
  Rcpp::NumericMatrix out(pairs.get_num_pairs(), 3);
  for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ ){
    out(idx, 1) = pairs.get_pair( idx ).first;
    out(idx, 2) = pairs.get_pair( idx ).second;
  }
  // add dimension and filter values
  for (int i = 0; i < out.nrow(); ++i){
    out(i, 0) = std::get<1>(filtration_vector2[out(i, 1)]);
    out(i, 1) = std::get<0>(filtration_vector2[out(i, 1)]);
    out(i, 2) = std::get<0>(filtration_vector2[out(i, 2)]);
  }
  // Rcpp::NumericMatrix out(3, 2);
  return out;
}

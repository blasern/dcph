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

// filtration_value, diameter, V1, V2, V3
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

//' Persistence from cover
//' 
//' Calculate the persistence from the cover
//' 
//' @param cover A cover
// [[Rcpp::export]]
Rcpp::NumericMatrix persistence_from_cover(Rcpp::S4 cover) {
  // check that input is a cover
  if (!cover.inherits("cover")) Rcpp::stop("Input must be a cover");
  
  // extract indices and diameter
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
  
  // // reindex
  // Rcpp::NumericMatrix distance_matrix = cover.slot("distance_matrix");
  // int n = distance_matrix.nrow();
  // // std::set<std::tuple<std::set<int>, double>> contained_in;
  // std::set<std::set<int>> contained_in;
  // for (int i = 0; i < n; ++i){
  //   std::set<int> contained;
  //   for (auto const& ind: indices){
  //     const bool is_in = ind.find(i) != ind.end();
  //     if (is_in){
  //       contained.insert(i);
  //     }
  //   }
  //   contained_in.insert(contained);
  // }
  
  // calculate 3-fold overlaps between nodes
  std::vector<filtration_tuple> overlap;
  Rcpp::NumericVector diams = Rcpp::NumericVector(3);
  for (unsigned int i = 0; i < indices.size(); ++i){
    std::set<int> loc_i = indices[i];
    diams(0) = death_diameters[i];
    for (unsigned int j = 0; j <= i; ++j){
      std::set<int> loc_j = indices[j];
      std::set<int> intersect_ij;
      set_intersection(loc_i.begin(), loc_i.end(), loc_j.begin(), loc_j.end(),
                       std::inserter(intersect_ij, intersect_ij.begin()));
      if (intersect_ij.size() == 0) continue; 
      diams(1) = death_diameters[j];
      for (unsigned int k = 0; k <= j; ++k){
        if ((k == j) && (j < i)) continue; 
        diams(2) = death_diameters[k];
        std::set<int> loc_k = indices[k];
        if (is_disjoint(loc_k, intersect_ij)) continue; 
        // set_intersection(loc_k.begin(), loc_k.end(), intersect_ij.begin(), intersect_ij.end(),
        //                  std::inserter(intersect_ijk, intersect_ijk.begin()));
        // if (intersect_ijk.size() == 0) continue; 
        overlap.push_back(std::make_tuple(Rcpp::max(diams), 
                                           2 - int(k == j) - int(j == i), 
                                           k, j, i));
      }
    }
  }
  
  // sort overlap
  sort(overlap.begin(), overlap.end(), compare_filtration);
  
  // save as map
  std::map<std::tuple<int, int, int>, std::tuple<int, double, int>> filtration;
  for(std::vector<filtration_tuple>::iterator it = overlap.begin(); it != overlap.end(); ++it) {
    filtration[std::make_tuple(std::get<2>(*it), std::get<3>(*it), std::get<4>(*it))] = 
      std::make_tuple(std::distance(overlap.begin(), it), 
                      std::get<0>(*it), std::get<1>(*it));
  }
  
  // define a boundary matrix
  phat::boundary_matrix< phat::vector_vector > boundary_matrix;
  
  // set the number of columns 
  boundary_matrix.set_num_cols(filtration.size());
  
  // set the dimension of the cell that a column represents:
  for(unsigned int i = 0; i < overlap.size(); ++i) {
    boundary_matrix.set_dim(i, std::get<1>(overlap[i]) );
  }
  
  // set the respective columns -- the columns entries have to be sorted
  std::vector< phat::index > temp_col;
  for(unsigned int i = 0; i < overlap.size(); ++i) {
    temp_col.clear();
    if (std::get<1>(overlap[i]) == 1){
      // 1-dimensional  
      Rcpp::NumericVector points = Rcpp::NumericVector(2);
      int p0 = std::get<2>(overlap[i]);
      int p1 = std::get<4>(overlap[i]);
      points(0) = std::get<0>(filtration[std::make_tuple(p0, p0, p0)]);
      points(1) = std::get<0>(filtration[std::make_tuple(p1, p1, p1)]);
      points.sort();
      temp_col.push_back(points(0));
      temp_col.push_back(points(1));
    }
    if (std::get<1>(overlap[i]) == 2){
      // 2-dimensional 
      Rcpp::NumericVector points = Rcpp::NumericVector(3);
      int p0 = std::get<2>(overlap[i]);
      int p1 = std::get<3>(overlap[i]);
      int p2 = std::get<4>(overlap[i]);
      points(0) = std::get<0>(filtration[std::make_tuple(p0, p1, p1)]);
      points(1) = std::get<0>(filtration[std::make_tuple(p0, p2, p2)]);
      points(2) = std::get<0>(filtration[std::make_tuple(p1, p2, p2)]);
      points.sort();
      temp_col.push_back(points(0));
      temp_col.push_back(points(1));
      temp_col.push_back(points(2));
    }
    boundary_matrix.set_col(i, temp_col);
  }
  temp_col.clear();
  // Rcpp::Rcout << "Boundary matrix defined." << std::endl;
  
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
    out(i, 0) = std::get<1>(overlap[out(i, 1)]);
    out(i, 1) = std::get<0>(overlap[out(i, 1)]);
    out(i, 2) = std::get<0>(overlap[out(i, 2)]);
  }
  // Rcpp::NumericMatrix out(3, 2);
  return out;
}

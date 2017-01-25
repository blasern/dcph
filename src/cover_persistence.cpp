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

#include "order.h"
//' Persistence from cover
//' 
//' Calculate the persistence from the cover
//' 
//' @param cover A cover
// [[Rcpp::export]]
Rcpp::NumericMatrix persistence_from_cover(Rcpp::S4 cover) {
  if (!cover.inherits("cover")) Rcpp::stop("Input must be a cover");
  // extract indices and diameter
  Rcpp::List subsets = cover.slot("subsets");
  Rcpp::List indices(subsets.length());
  Rcpp::NumericVector death_diameters(subsets.length()); 
  for (int i = 0; i < subsets.length(); ++i){
    Rcpp::S4 subset = subsets(i);
    indices(i) = subset.slot("indices");
    death_diameters(i) = subset.slot("death");
  }
  
  // calculate 3-fold overlaps between nodes
  int indices_length = indices.length();
  int len = std::pow(indices_length, 3)/6 + indices_length*5/6 + 1;
  Rcpp::NumericMatrix overlap = Rcpp::NumericMatrix(len, 6);
  int m = 0;
  Rcpp::NumericVector diams = Rcpp::NumericVector(3);
  for (int i = 0; i < indices_length; ++i){
    Rcpp::NumericVector loc_i = Rcpp::as<Rcpp::NumericVector>(indices[i]);
    diams(0) = death_diameters(i);
    for (int j = 0; j <= i; ++j){
      Rcpp::NumericVector loc_j = Rcpp::as<Rcpp::NumericVector>(indices[j]);
      Rcpp::NumericVector intersect_ij = Rcpp::intersect(loc_i, loc_j);
      diams(1) = death_diameters(j);
      for (int k = 0; k <= j; ++k){
        if ((k == j) && (j < i)) continue; 
        diams(2) = death_diameters(k);
        Rcpp::NumericVector loc_k = Rcpp::as<Rcpp::NumericVector>(indices[k]);
        Rcpp::NumericVector intersect_ijk = Rcpp::intersect(intersect_ij, loc_k);
        overlap(m, 0) = k;
        overlap(m, 1) = j;
        overlap(m, 2) = i;
        overlap(m, 3) = intersect_ijk.length();
        overlap(m, 4) = 2 - int(k == j) - int(j == i);
        overlap(m, 5) = Rcpp::max(diams);
        // Rcpp::Rcout << overlap(m, 0) << " " << overlap(m, 1) << " " << overlap(m, 2) << " " << overlap(m, 3) << " " << overlap(m, 4) << " " << overlap(m, 5) << std::endl;
        m ++;
      }
    }
  }
  // sort the matrix by filter value and dimension
  Rcpp::IntegerVector order = order5(overlap(Rcpp::_, 5), overlap(Rcpp::_, 4), overlap(Rcpp::_, 0), overlap(Rcpp::_, 1), overlap(Rcpp::_, 2)) - 1;
  Rcpp::NumericMatrix ordered_overlap = Rcpp::NumericMatrix(overlap.nrow(), overlap.ncol());
  for (int m = 0; m < overlap.nrow(); ++m){
    ordered_overlap(m, Rcpp::_) = overlap(order(m), Rcpp::_);
    // Rcpp::Rcout << ordered_overlap(m, 0) << " " << ordered_overlap(m, 1) << " " << ordered_overlap(m, 2) << " " << ordered_overlap(m, 3) << " " << ordered_overlap(m, 4) << " " << ordered_overlap(m, 5) << std::endl;
  }
  
  // remove non-overlaps
  Rcpp::NumericMatrix filtration = Rcpp::NumericMatrix(sum(ordered_overlap(Rcpp::_, 3) > 0), 6);
  int n = 0;
  for (int m = 0; m < ordered_overlap.nrow(); ++m){
    if (ordered_overlap(m, 3) > 0){
      filtration(n, Rcpp::_) = ordered_overlap(m, Rcpp::_);
      // Rcpp::Rcout << filtration(n, 0) << " " << filtration(n, 1) << " " << filtration(n, 2) << " " << filtration(n, 3) << " " << filtration(n, 4) << " " << filtration(n, 5) << std::endl;
      ++n;
    }
  }
  
  // define a boundary matrix
  phat::boundary_matrix< phat::vector_vector > boundary_matrix;
  // set the number of columns 
  boundary_matrix.set_num_cols(filtration.nrow());
  // set the dimension of the cell that a column represents:
  for (int i = 0; i < filtration.nrow(); ++i){
    boundary_matrix.set_dim(i, filtration(i, 4) );
  }
  
  // set the respective columns -- the columns entries have to be sorted
  std::vector< phat::index > temp_col;
  for (int i = 0; i < filtration.nrow(); ++i){
    temp_col.clear();
    if (filtration(i, 4) == 1){
      Rcpp::NumericVector points = Rcpp::NumericVector(2);
      points(0) = Rcpp::which_max((filtration(Rcpp::_, 0) == filtration(i, 0)) & (filtration(Rcpp::_, 1) == filtration(i, 0)) & (filtration(Rcpp::_, 2) == filtration(i, 0)));
      points(1) = Rcpp::which_max((filtration(Rcpp::_, 0) == filtration(i, 2)) & (filtration(Rcpp::_, 1) == filtration(i, 2)) & (filtration(Rcpp::_, 2) == filtration(i, 2)));
      points.sort();
      temp_col.push_back(points(0));
      temp_col.push_back(points(1));
    }
    if (filtration(i, 4) == 2){
      Rcpp::NumericVector points = Rcpp::NumericVector(3);
      points(0) = Rcpp::which_max((filtration(Rcpp::_, 0) == filtration(i, 0)) & (filtration(Rcpp::_, 1) == filtration(i, 1)) & (filtration(Rcpp::_, 2) == filtration(i, 1))); 
      points(1) = Rcpp::which_max((filtration(Rcpp::_, 0) == filtration(i, 0)) & (filtration(Rcpp::_, 1) == filtration(i, 2)) & (filtration(Rcpp::_, 2) == filtration(i, 2)));
      points(2) = Rcpp::which_max((filtration(Rcpp::_, 0) == filtration(i, 1)) & (filtration(Rcpp::_, 1) == filtration(i, 2)) & (filtration(Rcpp::_, 2) == filtration(i, 2)));
      points.sort();
      temp_col.push_back(points(0));
      temp_col.push_back(points(1));
      temp_col.push_back(points(2));
    }
    boundary_matrix.set_col(i, temp_col);
  }
  temp_col.clear();
  
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
  phat::persistence_pairs pairs;

  // choose an algorithm (choice affects performance) and compute the persistence pair
  // (modifies boundary_matrix)
  phat::compute_persistence_pairs< phat::twist_reduction >( pairs, boundary_matrix );

  // sort the persistence pairs by birth index
  pairs.sort();

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
    out(i, 0) = filtration(out(i, 1), 4);
    out(i, 1) = filtration(out(i, 1), 5);
    out(i, 2) = filtration(out(i, 2), 5);
  }
  // Rcpp::NumericMatrix out(3, 2);
  return out;
}


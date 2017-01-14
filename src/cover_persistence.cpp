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

//' Persistence from filtration
//' 
//' Calculate the persistence from the filtration
//' 
//' @param filtration A filtration
// [[Rcpp::export]]
Rcpp::NumericMatrix persistence_from_filtration(Rcpp::NumericMatrix filtration) {
  //if (!filtration.inherits("filtration")) Rcpp::stop("Input must be a filtration");
  
  // define a boundary matrix
  phat::boundary_matrix< phat::vector_vector > boundary_matrix;
  // set the number of columns 
  boundary_matrix.set_num_cols(filtration.nrow());
  // set the dimension of the cell that a column represents:
  for (int i = 0; i < filtration.nrow(); ++i){
    boundary_matrix.set_dim(i, filtration(i, 0) );
  }

  // set the respective columns -- the columns entries have to be sorted
  std::vector< phat::index > temp_col;
  
  for (int i = 0; i < filtration.nrow(); ++i){
    temp_col.clear();
    if (filtration(i, 0) == 1){
      temp_col.push_back( filtration(i, 1) - 1 );
      temp_col.push_back( filtration(i, 2) - 1 );
    }
    if (filtration(i, 0) == 2){
      temp_col.push_back( filtration(i, 4) - 1 );
      temp_col.push_back( filtration(i, 5) - 1 );
      temp_col.push_back( filtration(i, 6) - 1 );
    }
    boundary_matrix.set_col(i, temp_col);
  }
  temp_col.clear();
  
  // print some information of the boundary matrix:
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "The boundary matrix has " << boundary_matrix.get_num_cols() << " columns: " << std::endl;
  for( phat::index col_idx = 0; col_idx < boundary_matrix.get_num_cols(); col_idx++ ) {
    Rcpp::Rcout << "Colum " << col_idx << " represents a cell of dimension " << (int)boundary_matrix.get_dim( col_idx ) << ". ";
    if( !boundary_matrix.is_empty( col_idx ) ) {
      std::vector< phat::index > temp_col;
      boundary_matrix.get_col( col_idx, temp_col ); 
      Rcpp::Rcout << "Its boundary consists of the cells";
      for( phat::index idx = 0; idx < (phat::index)temp_col.size(); idx++ )
        Rcpp::Rcout << " " << temp_col[ idx ];
    }
    Rcpp::Rcout << std::endl;
  }
  Rcpp::Rcout << "Overall, the boundary matrix has " << boundary_matrix.get_num_entries() << " entries." << std::endl;  
  
  
  // define the object to hold the resulting persistence pairs
  phat::persistence_pairs pairs;

  // choose an algorithm (choice affects performance) and compute the persistence pair
  // (modifies boundary_matrix)
  phat::compute_persistence_pairs< phat::twist_reduction >( pairs, boundary_matrix );

  // sort the persistence pairs by birth index
  pairs.sort();

  // print the pairs:
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "There are " << pairs.get_num_pairs() << " persistence pairs: " << std::endl;
  for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
    Rcpp::Rcout << "Birth: " << pairs.get_pair( idx ).first << ", Death: " << pairs.get_pair( idx ).second << std::endl;

  // write to output
  Rcpp::NumericMatrix out(pairs.get_num_pairs(), 2);
  for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ ){
    out(idx, 0) = pairs.get_pair( idx ).first;
    out(idx, 1) = pairs.get_pair( idx ).second;
  }
  // Rcpp::NumericMatrix out(3, 2);
  return out;
}

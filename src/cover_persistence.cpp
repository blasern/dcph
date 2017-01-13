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
Rcpp::NumericMatrix persistence_from_filtration(Rcpp::Matrix filtration) {
  //if (!filtration.inherits("filtration")) Rcpp::stop("Input must be a filtration");
  
  // define a boundary matrix
  phat::boundary_matrix< phat::vector_vector > boundary_matrix;
  // set the number of columns 
  boundary_matrix.set_num_cols(filtration.ncol(););
  // set the dimension of the cell that a column represents:
  for (int i = 0; i < filtration.ncol(); ++i){
    boundary_matrix.set_dim(i, filtration(i, 0) );
  }

  // set the respective columns -- the columns entries have to be sorted
  std::vector< phat::index > temp_col;
  
  for (int i = 0; i < filtration.ncol(); ++i){
    boundary_matrix.set_col(i, temp_col);
    temp_col.push_back( 0 );
  }
  
  boundary_matrix.set_col( 0, temp_col );
  boundary_matrix.set_col( 1, temp_col );

  temp_col.push_back( 0 );
  temp_col.push_back( 1 );
  
  boundary_matrix.set_col( 2, temp_col );
  temp_col.clear();

  boundary_matrix.set_col( 3, temp_col );

  temp_col.push_back( 1 );
  temp_col.push_back( 3 );
  boundary_matrix.set_col( 4, temp_col );
  
  temp_col.clear();
  temp_col.push_back( 0 );
  temp_col.push_back( 3 );
  boundary_matrix.set_col( 5, temp_col );
  temp_col.clear();
  temp_col.push_back( 2 );
  temp_col.push_back( 4 );
  temp_col.push_back( 5 );
  boundary_matrix.set_col( 6, temp_col );
  temp_col.clear();
  
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
  return out;
}

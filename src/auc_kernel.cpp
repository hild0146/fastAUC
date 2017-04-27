// notes/questions ------------------------------------------------------------

//  n:  in C++ indices begin at 0 not 1 as in R

//  n:  to calculate the values necessary, we must calculate the number of observations
//      which are tied and which are less than each observation, grouped by status

// AUC kernel function in C++ -------------------------------------------------

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List auc_kernel(std::vector < double > test,
                std::vector < bool > status) {

  std::vector < double > out;
  std::vector < double > neg_dup_cnt;
  std::vector < double > pos_dup_cnt;
  std::vector < double > neg_all_cnt;
  std::vector < double > pos_all_cnt;
  std::size_t next_row;

  int ndup_pos = 0, ndup_neg = 0, neg_cnt = 0, pos_cnt = 0;
  double total_pos = 0, prev_val;


  // iterate over the data to determine the number of duplicates at each
  // unique value as well as the number of total true positive observations

  for (std::size_t row_i = 0; row_i < status.size(); ++row_i){

    // determine the total number of true positives
    if(status.at(row_i)){
      total_pos++;
    }


    // identify if the current record has a new test value
    if(row_i == 0 || prev_val != test.at(row_i)){

      next_row = row_i + 1;
      ndup_pos = 0;
      ndup_neg = 0;

      while(true){

        // if you're at the last row stop
        if(row_i == (test.size() - 1)){
          break;
        }

        // count duplicates in the next rows
        if(test.at(row_i) == test.at(next_row)){

          if(status.at(next_row)){
            ndup_pos++;
            pos_cnt++;
          } else{
            ndup_neg++;
            neg_cnt++;
          }

          next_row++;

          // see if you just looked at the last value
          if(next_row == test.size()) break;
        } else break;
      }
      if(status.at(row_i)){
        ndup_pos++;
        pos_cnt++;
        }
      else{
        ndup_neg++;
        neg_cnt++;
      }
      prev_val = test.at(row_i);
    }

  neg_dup_cnt.push_back(ndup_neg);
  pos_dup_cnt.push_back(ndup_pos);

  neg_all_cnt.push_back(neg_cnt);
  pos_all_cnt.push_back(pos_cnt);

  }

  // iterate through the data to determine X and Y values per DeLong
  for (std::size_t row_i = 0; row_i < test.size(); ++row_i){

    if(status.at(row_i)){
      out.push_back(neg_all_cnt.at(row_i) - 0.5 * neg_dup_cnt.at(row_i));
    } else {
      out.push_back(total_pos - pos_all_cnt.at(row_i) + 0.5 * pos_dup_cnt.at(row_i));
    }

  }

  List res;
  res["out"] = out;
  return(res);
}

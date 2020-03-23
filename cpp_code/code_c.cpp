#include <Rcpp.h>
using namespace Rcpp;

struct Nodes {
  IntegerVector ind;
  NumericMatrix Dsum;
  NumericMatrix Rsum;
  NumericMatrix R;
  IntegerMatrix I;
  int m;
  
  // constructor
  Nodes (IntegerVector ind, NumericMatrix Dsum, NumericMatrix Rsum,
         NumericMatrix R, NumericMatrix I, int m) :
    ind(ind), Dsum(Dsum), Rsum(Rsum), R(R), I(I), m(m) {};
};




// Internal function
// Given a vector X, a value k, and B=length(X),
// it returns TRUE if the critical value
// (k-th statistic when sorted in increasing order)
// is negative
// Notice that L<0 when no non-rejection has been found
//and U<0 when a certain rejection has been found

// [[Rcpp::export]]
bool Q(const NumericVector &X, const int &k, const int &B){
  int t = B-k+1; // threshold: min. number of non-negative el. to have q>= 0
  int n = 0; // number of non-negative el.
  bool out = TRUE;
  
  for (int i = 0; i <= B-1; ++i) {
    if (X[i] >= 0){
      ++n;
      if (n >= t){ // if n reaches the threshold, then q>=0
        out = FALSE;
        break;
      } else if (i+1-n >= k) { // if the negative el. found are already k or more, then q<0
        break;
      }
    }
  }
  return out;
}




// Internal function
// Given a matrix A and c=ncol(A),
// it finds the index of the first column in A with no positive element
// (= the last index, c-1, if such column does not exist)

int find_col(NumericMatrix &A, int &c){
  for (int j = 0; j <= c-1; ++j) {
    if (max(A(_,j)) <= 0) return j;
  }
  return c-1; // if there is no such column, then j is the last one
}




// Internal function.
// Given a node and the constants k and B=nrow(D),
// it checks both bounds. It returns non_rej (T if there is a non-rejection)
// and the new vector of indecisive sizes.

// const IntegerVector ind, const NumericMatrix Dsum, const NumericMatrix Rsum, const NumericMatrix R, const NumericMatrix I, const int m,

List bounds_both(IntegerVector ind, NumericMatrix Dsum,
                 NumericMatrix Rsum, NumericMatrix R, IntegerMatrix I,
                 int m, const int &k, const int &B){
  
  int H = ind.length();
  LogicalVector up (H, FALSE); // keeps track of indecisive sizes
  int h = 0;
  int v = ind[h];

  bool low = Q(Dsum(_,v), k, B);
  
  if (!low){
    List out = List::create(Named("non_rej")=!low, Named("ind") = ind[up]);
    return out;
  }
  
  // if v=0, then up=FALSE, otherwise it needs to be checked
  if (v > 0) up[h] = !Q(Rsum(_,v), k, B);
  
  // first column in R with no positive element
  int j = find_col(R, m);
  
  // bounds for v < j (before the upper c.v. decreases)
  while (low && v < j-1 && h < H-1){
    ++h;
    v = ind[h];
    low = Q(Dsum(_,v), k, B);
    up[h] = !Q(Rsum(_,v), k, B);
  }
  
  // bounds for v >= j (until the upper c.v. becomes negative)
  while (low && up(h) && h < H-1){
    ++h;
    v = ind[h];
    low = Q(Dsum(_,v), k, B);
    up[h] = !Q(Rsum(_,v), k, B);
  }
  
  List out = List::create(Named("non_rej")=!low, Named("ind") = ind[up]);
  return out;
}




// Internal function.
// Given a node and the constants k and B=nrow(D),
// it checks the upper bound. It returns non_rej=F
// and the new vector of indecisive sizes.

List bounds_upper(IntegerVector ind, NumericMatrix Dsum,
                  NumericMatrix Rsum, NumericMatrix R, IntegerMatrix I,
                  int m, const int &k, const int &B){
  
  int H = ind.length();
  LogicalVector up (H, FALSE); // keeps track of indecisive sizes
  int h = 0;
  int v = ind[h];
  
  // if v=0, then up=FALSE, otherwise it needs to be checked
  if (v > 0) up[h] = !Q(Rsum(_,v), k, B);
  
  // first column in R with no positive element
  int j = find_col(R, m);
  
  // upper bounds for v < j (before the upper c.v. decreases)
  while (v < j-1 && h < H-1){
    ++h;
    v = ind[h];
    up[h] = !Q(Rsum(_,v), k, B);
  }
  
  // bounds for v >= j (until the upper c.v. becomes negative)
  while (up(h) && h < H-1){
    ++h;
    v = ind[h];
    up[h] = !Q(Rsum(_,v), k, B);
  }
  
  bool non_rej = FALSE;
  List out = List::create(Named("non_rej")=non_rej, Named("ind") = ind[up]);
  return out;
}




// Internal function.
// Given a node, the constants k and B=nrow(D), and a boolean both,
// it checks both bounds or the upper bound.
// It returns non_rej (T if there is a non-rejection)
// and the new vector of indecisive sizes.
// In particular, it finds the first column R[,j] such that all the el. are negative/null
// It checks the bounds up to j-1, and then the following bounds
// until one upper bound becomes negative.

List compute_bounds(IntegerVector ind, NumericMatrix Dsum,
                    NumericMatrix Rsum, NumericMatrix R,
                    IntegerMatrix I, int m, const int &k, const int &B,
                    const bool both){
  List out;
  if (both) out = bounds_both(ind, Dsum, Rsum, R, I, m, k, B);
  else out = bounds_upper(ind, Dsum, Rsum, R, I, m, k, B);
  return out;
}




// Internal function - Branch and Bound
// Given an element x, a matrix A, a row a, and L=ncol(A)
// it returns the first index where the x appears in the row A(a,_)
// (L if it does not appear)

// [[Rcpp::export]]
int match_el (int &x, IntegerMatrix &A, int a, int &L){
  int out = L;
  for (int i = 0; i <= L-1; ++i) {
    if (A(a,i) == x) {
      out = i;
      break;
    }
  }
  return out;
}





// Internal function - Branch and Bound
// Given a node (ind, Dsum, Rsum, R, I) and the dimenzion m=ncol(R)
// the initial matrix D_0 and m_0=ncol(D_0)
// it computes the two subspaces 1 (remove) and 2 (keep)
// when considering the lowest statistic
// In subspace 2, the indecisive sizes decrease of 1 unit
// In this case, node1 will be analysed with bounds_both,
// and node2 with bounds_upper. Hence if zero is in node_2, we can remove it
// (since the upper bound coincides with the lower bound, which is < 0)

List gen_nodes_low (IntegerVector ind, NumericMatrix Dsum, NumericMatrix Rsum,
                     NumericMatrix R, IntegerMatrix I, int m, NumericMatrix &D_0,
                     const int &B, const int &m_0){
  int m_old = m;
  --m;
  int h = I(0,m); //index
  int iD = m_0 - m - 1; // index in D_0
  NumericMatrix R_new = R(_,Range(0,m-1));
  IntegerMatrix I_new = I(_,Range(0,m-1));
  
  // First node: remove case
  Dsum = Dsum(_,Range(1,m+1)); // no need to save Dsum
  NumericMatrix Rsum_1 = Rsum(_,Range(0,m));
  
  
  // ----- special case: if the only indecisive size is 1, node_2 is removed -----
  
  if (ind.length() == 1 && ind[0] == 1){
    // Second node: keeping case NULL
    IntegerVector ind_2 = R_NilValue;
    NumericMatrix Dsum_2 = R_NilValue;
    NumericMatrix Rsum_2 = R_NilValue;

    // for each row, we find the index z corresponding to h (predictor to be removed)
    for (int b = 0; b <= B-1; ++b) {
      int z = match_el (h, I, b, m_old);
      
      // index z is removed (from there to the end, all elements are shifted)
      for (int j = z; j <= m-1; ++j){
        R_new(b,j) = R(b,j+1);
        I_new(b,j) = I(b,j+1);
      }
      
      for (int j = z; j <= m; ++j){
        Rsum_1(b,j) = Rsum(b,j+1);
      }
      
      // the test statistic corresponding to h is subtracted to each element of Dsum_1
      for (int j = 0; j <= m; ++j) Dsum(b,j) -= D_0(b,iD);
      
      // if h was not the last element of X.I(b,_), the test statistic
      // is removed from z until the end
      for (int j = z; j <= m; ++j) Rsum_1(b,j) -= D_0(b,iD);
    }
    
    List out  = List::create(Named("ind_1") = ind, Named("Dsum_1") = Dsum,
                             Named("Rsum_1") = Rsum_1, Named("ind_2") = ind_2,
                             Named("Dsum_2") = Dsum_2, Named("Rsum_2") = Rsum_2,
                             Named("R") = R_new, Named("I") = I_new);
    return out;
    
  }
  
  // Second node: keeping case
  // if zero is in ind_2, it is removed
  IntegerVector ind_2 = ind - 1;
  if (ind_2[0]==0) ind_2 = ind_2[Range(1, ind_2.length()-1)];
  NumericMatrix Dsum_2 = clone(Dsum);
  NumericMatrix Rsum_2 = Rsum(_,Range(0,m));
  
  // for each row, we find the index z corresponding to h (predictor to be removed)
  for (int b = 0; b <= B-1; ++b) {
    int z = match_el (h, I, b, m_old);
    
    // index z is removed (from there to the end, all elements are shifted)
    for (int j = z; j <= m-1; ++j){
      R_new(b,j) = R(b,j+1);
      I_new(b,j) = I(b,j+1);
    }
    
    for (int j = z; j <= m; ++j){
      Rsum_1(b,j) = Rsum(b,j+1);
      Rsum_2(b,j) = Rsum(b,j+1);
    }
    
    // the test statistic corresponding to h is subtracted to each element of Dsum_1
    for (int j = 0; j <= m; ++j) Dsum(b,j) -= D_0(b,iD);
    
    // if h was not the last element of X.I(b,_), the test statistic
    // is removed from z until the end
    for (int j = z; j <= m; ++j) Rsum_1(b,j) -= D_0(b,iD);
    
    // if h was not the first element of X.I(b,_), the test statistic
    // is added from the beginning until z-1
    for (int j = 0; j <= z-1; ++j) Rsum_2(b,j) += D_0(b,iD);
  }
  
  List out  = List::create(Named("ind_1") = ind, Named("Dsum_1") = Dsum,
                           Named("Rsum_1") = Rsum_1, Named("ind_2") = ind_2,
                           Named("Dsum_2") = Dsum_2, Named("Rsum_2") = Rsum_2,
                           Named("R") = R_new, Named("I") = I_new);
  return out;
}




// Internal function - Branch and Bound
// Given a node (ind, Dsum, Rsum, R, I) and the dimenzion m=ncol(R)
// the initial matrix D_0 and m_0=ncol(D_0)
// it computes the two subspaces 1 (remove) and 2 (keep)
// when considering the highest statistic
// In subspace 2, the indecisive sizes decrease of 1 unit
// In this case, node1 will be analysed with bounds_upper,
// and node2 with bounds_both


List gen_nodes_high (IntegerVector ind, NumericMatrix Dsum, NumericMatrix Rsum,
                     NumericMatrix R, IntegerMatrix I, int m, NumericMatrix &D_0,
                     const int &B){
  int m_old = m;
  --m;
  int h = I(0,0); //index
  int iD = m; // index in D_0
  NumericMatrix R_new = R(_,Range(0,m-1));
  IntegerMatrix I_new = I(_,Range(0,m-1));
  
  // First node: remove case
  Dsum = Dsum(_,Range(0,m)); // no need to save Dsum
  NumericMatrix Rsum_1 = Rsum(_,Range(0,m));
  
  // Second node: keeping case
  IntegerVector ind_2 = ind - 1;
  NumericMatrix Dsum_2 = clone(Dsum);
  NumericMatrix Rsum_2 = Rsum(_,Range(0,m));
  
  // for each row, we find the index z corresponding to h (predictor to be removed)
  for (int b = 0; b <= B-1; ++b) {
    int z = match_el (h, I, b, m_old);
    
    // index z is removed (from there to the end, all elements are shifted)
    for (int j = z; j <= m-1; ++j){
      R_new(b,j) = R(b,j+1);
      I_new(b,j) = I(b,j+1);
    }
    
    for (int j = z; j <= m; ++j){
      Rsum_1(b,j) = Rsum(b,j+1);
      Rsum_2(b,j) = Rsum(b,j+1);
    }
    
    // the test statistic corresponding to h is added to each element of Dsum_2
    for (int j = 0; j <= m; ++j) Dsum_2(b,j) += D_0(b,iD);
    
    // if h was not the last element of X.I(b,_), the test statistic
    // is removed from z until the end
    for (int j = z; j <= m; ++j) Rsum_1(b,j) -= D_0(b,iD);
    
    // if h was not the first element of X.I(b,_), the test statistic
    // is added from the beginning until z-1
    for (int j = 0; j <= z-1; ++j) Rsum_2(b,j) += D_0(b,iD);
  }
  
  List out  = List::create(Named("ind_1") = ind, Named("Dsum_1") = Dsum,
                           Named("Rsum_1") = Rsum_1, Named("ind_2") = ind_2,
                           Named("Dsum_2") = Dsum_2, Named("Rsum_2") = Rsum_2,
                           Named("R") = R_new, Named("I") = I_new);
  
  return out;
}



// Internal function - Branch and Bound
// Given a node (ind, Dsum, Rsum, R, I) and the dimenzion m=ncol(R)
// the initial matrix D_0 and m_0=ncol(D_0)
// the settings from_low and first_rem,
// it computes the two subspaces 1 and 2
// (1 corresponds to the subspace that will be explored first)

List generate_nodes (IntegerVector ind, NumericMatrix Dsum, NumericMatrix Rsum,
                     NumericMatrix R, IntegerMatrix I, int m, NumericMatrix &D_0,
                     const int &B, const int &m_0, const bool from_low,
                     const bool first_rem){
  
  List out;
  if (from_low) out = gen_nodes_low(ind, Dsum, Rsum, R, I, m, D_0, B, m_0);
  else out = gen_nodes_high(ind, Dsum, Rsum, R, I, m, D_0, B);
  
  // if we start by keeping the selected statistic, then we exchange the names:
  // 1 = keep case, 2 = remove case
  if (!first_rem) out.names() = CharacterVector::create("ind_2", "Dsum_2", "Rsum_2", "ind_1",
      "Dsum_1", "Rsum_1", "R", "I");
  
  return out;
}







//ctrp_bab <- function(ind_0, D_0, R_0, I_0, Dsum_0, Rsum_0,
// k=ceiling(0.95*nrow(R)), m_0=ncol(D_0), B=nrow(D_0), n_max=10000,
//from_low=T, first_rem=T){

//# when from_low=first_rem (either T or F), in the first loop we compute both bounds
//# and in the second only the upper
//  both_first <- (from_low == first_rem)
//}

List ctrp_bab (IntegerVector &ind_0, NumericMatrix &D_0, NumericMatrix &R_0,
               IntegerMatrix &I_0, NumericMatrix Dsum_0, NumericMatrix Rsum_0,
               const int &k, int &m_0, const int &B,
               const bool &from_low, const bool &first_rem, const int &n_max){
  
  // when from_low=first_rem (either T or F), in the first loop we compute both bounds
  // and in the second only the upper
  const bool both_first = (from_low == first_rem);
  
  // lists of the nodes element to be examined (every node is a list itself)
  List IndList = List::create(ind_0);
  List DsumList = List::create(Dsum_0);
  List RsumList = List::create(Rsum_0);
  List RList = List::create(R_0);
  List IList = List::create(I_0);
  
  int L = 1; // lists length
  int BAB = 0; // number of steps
  int m = m_0;
  List g;
  List cb;
  bool reject = FALSE;
  bool non_reject = TRUE;
  
  
  while (BAB < n_max){
    
    // we take the last element added to the list, generate the two branches
    --m;
    g = generate_nodes (IndList[L], DsumList[L], RsumList[L], RList[L], IList[L], m,
                        D_0, B, m_0, from_low, first_rem);
    
    // then we remove the original node and add the node 2 (if not null)
    bool cond = (g["ind_1"] != R_NilValue);
    if (g["ind_2"] != R_NilValue){
      IndList[L] = g["ind_2"];
      DsumList[L] = g["Dsum_2"];
      RsumList[L] = g["Rsum_2"];
      RList[L] = g["R"];
      IList[L] = g["I"];
    }else{
      IndList.erase(L);
      DsumList.erase(L);
      RsumList.erase(L);
      RList.erase(L);
      IList.erase(L);
      --L;
    }
    
    
    // FIRST LOOP: depth-first search
    // we evaluate node 1, and iterate until a node 1 can be closed
    while (cond && BAB < n_max){
      ++BAB;
      cb = compute_bounds(g["ind_1"], g["Dsum_1"], g["Rsum_1"], g["R"], g["I"],
                          m, k, B, both_first);
      if (cb["non_rej"]) return List::create(Named("non_rej") = non_reject, Named("BAB") = BAB);
      
      cond = (cb["ind"] != R_NilValue);
      if (cond){
        --m;
        // we generate new nodes, after updating the indecisive sizes
        // (g$ind_1 becomes cb$ind)
        g = generate_nodes (cb["ind"], g["Dsum_1"], g["Rsum_1"], g["R"], g["I"], m,
                            D_0, B, m_0, from_low, first_rem);
        cond = (g["ind_1"] != R_NilValue);
        if (g["ind_2"] != R_NilValue){
          IndList.push_back(g["ind_2"]);
          DsumList.push_back(g["Dsum_2"]);
          RsumList.push_back(g["Rsum_2"]);
          RList.push_back(g["R"]);
          IList.push_back(g["I"]);
          ++L;
        }
      }
    }
    
    
    // if the list becomes empty before a non-rejection is found, we reject
    if (L == 0) return List::create(Named("non_rej") = reject, Named("BAB") = BAB);
    
    
    // SECOND LOOP: search of the remaining nodes
    while (!cond && BAB < n_max){
      ++BAB;
      
      IntegerMatrix temp = IList[L];
      m = temp.ncol();
      
      cb = compute_bounds(IndList[L], DsumList[L], RsumList[L], RList[L], IList[L],
                          m, k, B, !both_first);
      
      if (cb["non_rej"]) return List::create(Named("non_rej") = non_reject, Named("BAB") = BAB);
      
      cond = (cb["ind"] != R_NilValue);
      if (cond){
        IndList[L] = cb["ind"]; // update of indecisive sizes after evaluation
      } else {
        // the node is removed
        IndList.erase(L);
        DsumList.erase(L);
        RsumList.erase(L);
        RList.erase(L);
        IList.erase(L);
        --L;
        // if the list becomes empty before a non-rejection is found, we reject
        if (L == 0) return List::create(Named("non_rej") = reject, Named("BAB") = BAB);
      }
    }
  }
  return List::create(Named("non_rej") = R_NilValue, Named("BAB") = BAB);
}





// Main function
// [[Rcpp::export]]
List cpp_test(IntegerVector &ind, NumericMatrix &D, NumericMatrix &R,
              IntegerMatrix &I, NumericMatrix &Dsum,
              NumericMatrix &Rsum, const int &k, int &m, const int &B,
              const bool &from_low, const bool &first_rem, const int &n_max){
  
  // test
  List cb = compute_bounds(ind, Dsum, Rsum, R, I, m, k, B, TRUE);
  ind = cb["ind"];
  
  // if there is a non-rejection or there are no indecisives
  if (cb["non_rej"] || ind.length() == 0){
    List out = List::create(Named("non_rej")=cb["non_rej"], Named("BAB") = 0);
    return out;
  }else{
    const int m_0 = m;
    //List out = List::create(Named("non_rej")=cb["non_rej"], Named("BAB") = 10);
    List out = ctrp_bab (ind, D, R, I, Dsum, Rsum, k, m, B, from_low, first_rem, n_max);
    return out;
  }
}
  


  
  
  
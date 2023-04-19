#include <R.h>
#include <Rinternals.h>

// a function that looks for overrepresentation of dinucleotides using
// a running sum. Similar to cpg_distribute.c but for an arbitrary
// dinucleotide.

// A struct for holding the region and associated statistics:

// number of stats included in region:
// convenient for memcpy operations
// score start end nuc1_N, nuc2_N, dinuc_N rev_dinuc_N,
#define N_STATS 7
#define SC_i 0
#define ST_i 1
#define END_i 2
#define NUC1_i 3
#define NUC2_i 4
#define DINUC_i 5
#define R_DINUC_i 6

#define UC(x) ((x) & ~0x20)

#define REG_N 256

const char *col_names[N_STATS] = { "score", "start", "end", "nuc1_N", "nuc2_N", "dinuc_N", "rev_dinuc_N" };

struct regions {
  int nrow; // the number of statistics currently held
  int capacity; // the number of rows for which memory has been allocated
  int *stats;
};

struct regions init_regions(int n){
  struct regions reg;
  // do not allow n to be less than 1
  n = n < 1 ? 1 : n;
  // set all values to 0.
  reg.nrow = 0;
  reg.capacity = n;
  reg.stats = malloc( sizeof(int) * N_STATS * n );
  return(reg);
}

void free_regions(struct regions *reg){
  free( reg->stats );
}

void push_regions(struct regions *reg, int stats[N_STATS])
{
  if( reg->nrow >= reg->capacity ){
    int new_capacity = reg->capacity * 2;
    int *new_stats = malloc( sizeof(int) * N_STATS * new_capacity );
    for(int i=0; i < N_STATS; ++i)
      memcpy( new_stats + i * new_capacity, reg->stats + i * reg->capacity, sizeof(int) * reg->capacity );
    reg->capacity = new_capacity;
    free( reg->stats );
    reg->stats = new_stats;
  }
  for(int i=0; i < N_STATS; ++i)
    reg->stats[ reg->nrow + reg->capacity * i ] = stats[i];
  reg->nrow++;
}

// I thought this would not be necessary, but I'm afraid it is necessary
// to traverse the region again counting the stats. This is because we do
// do not know the end of the range until the point when it is finished
// note that end is inclusive and must be smaller than the length of the sequence
void scan_region(const char *seq, int score, int start, int end, char n1, char n2, struct regions *reg){
  int stats[N_STATS] = {score, start, end, 0, 0, 0, 0};
  for(int i=start; i <= end; ++i){
    if(UC(seq[i]) == n1){
      stats[NUC1_i]++;
      if(i < end && UC(seq[i+1]) == n2)
	stats[DINUC_i]++;
    }
    if(UC(seq[i]) == n2){
      stats[NUC2_i]++;
      if(i < end && UC(seq[i+1]) == n1)
	stats[R_DINUC_i]++;
    }
  }
  push_regions(reg, stats);
}

// we should also consider a minimum score or width for a region
// as we probably do not want to report every instance of the dinucleotide
// across the region.
// seq : 0 terminated char array
// dn_score : the score associated with a the dinucleotide
struct regions find_regions(const char *seq, int dn_score, char n1, char n2){
  int score = 0;
  int last_score = 0;
  int max_score = 0;
  int max_i = 0;
  int start_i = 0;
  // initialise regions with a reasonable size
  struct regions reg = init_regions(REG_N);
  // then loop through:
  int i=0;
  while(seq[i]){
    // convert to upper case as we go along..
    // note that it is safe to test seq[i+1] as the sequence
    // should be 0 terminated.
    if(UC(seq[i]) == n1 && UC(seq[i+1]) == n2)
      score += dn_score;
    else
      score--;
    // but don't allow a negative score!
    score = score < 0 ? 0 : score;
    // The order of the following conditionals is important and easy
    // to get wrong!
    // (and they may still be wrong!!)
    if(!last_score && score)  // the beginning of a region
      start_i = i;
    if(score > max_score){
      max_score = score;
      max_i = i;
    }
    if(!score && last_score){ // a region has been defined
      scan_region(seq, max_score, start_i, max_i, n1, n2, &reg);
      i = max_i;
      score = last_score = max_score = max_i = start_i = 0;
    }
    // we always want the last score to be equal to the score here
    last_score = score;
    ++i;
  }
  // we may end in a region;
  if(score)
    scan_region(seq, max_score, start_i, max_i, n1, n2, &reg);
  return(reg);
}

// A function that can be called from R:
// seq_r should be a character vector of length 2
// containing:
// 1. The sequence to be scanned
// 2. the dinucleotide to be searched for
SEXP di_nuc_regions(SEXP seqs_r, SEXP dn_score_r){
  if(TYPEOF(seqs_r) != STRSXP || length(seqs_r) != 2)
    error("The first argument should be a character vector of length 2: sequence, dinucleotide");
  if(TYPEOF(dn_score_r) != INTSXP || length(dn_score_r) != 1)
    error("The second argument should be an integer vector of length 1 giving the score associated with the given dinucleotide");
  SEXP seq_r = STRING_ELT(seqs_r, 0);
  SEXP dinuc_r = STRING_ELT(seqs_r, 1);
  if(length(seq_r) < 10)
    error("Please provide a sequence at least 10bp long to scan");
  if(length(dinuc_r) != 2)
    error("The second argument should containg a string of length 2!");
  const char *seq = CHAR(seq_r);
  const char *dinuc = CHAR(dinuc_r);
  char n1 = dinuc[0];
  char n2 = dinuc[1];
  int dn_score = asInteger(dn_score_r);
  if(dn_score <= 0)
    error("dn_score should be positive");
  struct regions reg = find_regions(seq, dn_score, n1, n2);
  // and then allocate some memory as a matrix. Ideally we would like to set the
  // column names. That, however, is _much_ easier to do in an R wrapper function
  // so I will not do it here. It is _better_ for it to be done here and we have
  // a data structure that we can use for it in the colnames slot.
  // maybe use setAttrib( R_NamesSymbol ), though that may not work for a matrix
  SEXP ret_data = PROTECT( allocMatrix( INTSXP, reg.nrow, N_STATS ));
  int *ret_ptr = INTEGER(ret_data);
  for(int i=0; i < N_STATS; ++i)
    memcpy( ret_ptr + i * reg.nrow, reg.stats + i * reg.capacity, sizeof(int) * reg.nrow );
  free_regions(&reg);
  UNPROTECT(1);
  return(ret_data);
}

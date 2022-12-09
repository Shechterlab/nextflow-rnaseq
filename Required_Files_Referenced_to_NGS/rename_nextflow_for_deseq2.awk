##This script takes an input file and splits the first column name at the underscore then reprints the first part prior to the underscore as a new column

{ if (NR > 1) {
    # all data lines
    split($1, a, "_");
    print a[1]","$1,$2,$3,$4;
 } else {
    # header line
    print "condition"","$1, $2, $3, $4;
  }
}

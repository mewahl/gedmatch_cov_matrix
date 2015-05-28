# gedmatch_cov_matrix
Automates construction of an autosomal genetic covariance matrix and shared segment list for the top n matches to a GEDmatch kit number

As of this writing, this script can be used to automate one-to-one and one-to-many comparisons on the website GEDmatch.  These comparisons are used to identify the top n matches to a person of interest, then determine the pairwise autosomal genetic relatedness of all of these matches to the proband and each other. (Persons with non-zero genetic covariance likely share a common recent ancestor: the larger the group, the easier it may become to identify this ancestor genealogically.) A list of shared segments is also constructed.

This script requires the following modules:
* BeautifulSoup
* mechanize
* cookielib
* html2text
* NumPy
* argparse
* cPickle
* time
* sys

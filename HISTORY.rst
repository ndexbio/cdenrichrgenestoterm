=======
History
=======

0.4.0 (2021-03-09)
----------------------

* Updated to use gseapy 0.10.4 because Enrichr API changed hostname
  and older versions no longer worked 
  `Issue #2 <https://github.com/idekerlab/cdenrichrgenestoterm/issues/2>`_

0.3.1 (2020-08-13)
----------------------

* Removed debugging print statement that was sending debugging output
to standard output causing CDAPS app to fail.

* Fixed bug where best hit from last gene set input was
  returned instead of best hit from any of the gene sets

0.2.0 (2019-12-13)
--------------------

* Added term_size which denotes number of genes associated with term to output JSON

0.1.0 (2019-09-19)
------------------

* First release on PyPI.

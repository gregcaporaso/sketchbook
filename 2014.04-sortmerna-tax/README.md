This is a first pass at doing rRNA read assignment with [sortmerna](https://github.com/biocore/sortmerna). This code is **completely untested**, but provides the framework for starting to test this as an alternative to QIIME's existing taxonomy assigners.  

To use this code you'll need:

 * [QIIME 1.8.0](www.qiime.org) or later
 * [scikit-bio](www.scikit-bio.org) (latest GitHub version)
 * The sortmerna binaries in your ``PATH`` environment variable.

```
git clone https://gist.github.com/10506167.git smr-tax
mkdir tiny-test
python -c "from qiime.test import write_test_data ; write_test_data('./tiny-test')"
python smr-tax/assign-taxonomy-sortmerna.py -i tiny-test/seqs -r tiny-test/refseqs -t tiny-test/refseqs_tax -o assigned_tax
```

The resulting ``assigned_tax/assignments.tsv`` will contain the taxonomy assignments, ``assigned_tax/log.blast`` will contain the blast-formatted raw output of sortmerna, and ``assigned_tax/cmds.log`` will contain the exact sortmerna commands that were run.

By default, the sortmerna binary calls would look like:

```
indexdb_rna --ref tiny-test/refseqs,/Users/caporaso/temp/tmphbWyNB/ref.idx -v
sortmerna --ref tiny-test/refseqs,/Users/caporaso/temp/tmphbWyNB/ref.idx --reads tiny-test/seqs --blast 1 --aligned assigned_tax/log --min_lis 2 --best 3 -v
```

If you want to experiment with parameters that are not passed in those commands, you can use the corresponding ``params`` parameters. For example:

```
python smr-tax/assign-taxonomy-sortmerna.py -i tiny-test/seqs -r tiny-test/refseqs -t tiny-test/refseqs_tax -o assigned_tax_ --sortmerna_params "-F" --indexdb_params "-m 4000 -L 19"
```

sam2bed
========

This is the simple script for converting a sam file to a bed format. Every
string in the sam file is translated into a bed interval, and all related
information (like attributes) is preserved. Also, it is possible to split reads
by the N (and D) CIGAR operation that can be useful for processing fragments of
spliced reads (`--split/-s` option).

### Usage

```bash
python3 sam2bed [-h] -i INPUT -o OUTPUT [-s] [-S] [--reduced]
```

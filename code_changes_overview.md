# Overview of code changes

This document describes all manual changes made to the source code of `CBIIT/ldsc` and
`qlu-lab/PIGEON`, needed to get both tools working on Windows with more recent versions of
`bitarray` and `pandas`.

---

## CBIIT/ldsc

### 1. `ldscore/ldscore.py` — bitarray compatibility

**Problem:** newer versions of the `bitarray` package have `.decode()` return an iterator
instead of a list, which caused a `TypeError` when converting to a numpy array.

**Line ~395**, in the `nextSNPs` function:

```python
# Was:
X = np.array(slice.decode(self._bedcode), dtype="float64").reshape((b, nru)).T

# Becomes:
X = np.array(list(slice.decode(self._bedcode)), dtype="float64").reshape((b, nru)).T
```

---

### 2. `ldsc.py` — external `gzip` call doesn't work on Windows

**Problem:** the code calls the external `gzip` program via `subprocess.call(['gzip', ...])`,
which doesn't exist on Windows (`FileNotFoundError: [WinError 2]`).

**Two locations** in the file (around line 352 and around lines 382-386) contain the same call
and were both fixed the same way:

```python
# Was:
call(['gzip', '-f', out_fname])

# Becomes:
import gzip as gz_module
import shutil as shutil_module
with open(out_fname, 'rb') as f_in, gz_module.open(out_fname + '.gz', 'wb') as f_out:
    shutil_module.copyfileobj(f_in, f_out)
os.remove(out_fname)
```

### 3. `ldsc.py` — missing `import os`

**Problem:** the replacement code above uses `os.remove()`, but `os` was not imported anywhere
in the file, causing a `NameError`.

**Line 20**, with the existing imports:

```python
# Was:
import time, sys, traceback, argparse

# Becomes:
import time, sys, traceback, argparse, os
```

---

## qlu-lab/PIGEON

### 4. `pigeon.py` — deprecated pandas option

**Problem:** `pd.set_option('precision', 4)` has become ambiguous in newer pandas versions
(multiple options match the short name `'precision'`), causing an `OptionError`.

**Line 28:**

```python
# Was:
pd.set_option('precision', 4)

# Becomes:
pd.set_option('display.precision', 4)
```

### 5. `munge_sumstats.py` — deprecated pandas argument

**Problem:** `delim_whitespace=True` has been removed from newer pandas versions (`TypeError`).

**Two locations**, around line 672 and line 692:

```python
# Was:
delim_whitespace=True

# Becomes:
sep=r'\s+'
```

### 6. `ldsc_mod/parse.py` — same deprecated pandas argument

**Two locations**, line 14 and line 260:

```python
# Was:
delim_whitespace=True

# Becomes:
sep=r'\s+'
```

### 7. `ldsc_mod/sumstats.py` — typos in `_get_....._table`

**Problem:** four typos (it instead of i)

**Two locations** in the file (lines 609, 664, 728 and 765 contain the same call
and were all fixed the same way:

```python
# Was:
if args.samp_prev is not None and \
        args.pop_prev is not None and \
        all((i is not None for i in args.samp_prev)) and \
        all((i is not None for it in args.pop_prev)):

# Becomes:
if args.samp_prev is not None and \
        args.pop_prev is not None and \
        all((i is not None for i in args.samp_prev)) and \
        all((i is not None for i in args.pop_prev)):
```

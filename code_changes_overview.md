# Overview of code changes

This document describes all manual changes made to the source code of `CBIIT/ldsc` and
`qlu-lab/PIGEON`, needed to get both tools working on Windows with more recent versions of
`bitarray` and `pandas`. Keep this file so you can reproduce these changes on a fresh clone,
or share it with others using the same tools.

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

### 7. `pigeon.py` — missing `type=float` on liability-scale flags

**Problem:** `--samp-prev` and `--pop-prev` were passed through as text (string) instead of
a number, which later caused a `TypeError` in `np.isnan()`.

**Lines 103-106:**

```python
# Was:
parser.add_argument('--samp-prev', default=None,
    help='Sample prevalence of binary phenotype (for conversion to liability scale).')
parser.add_argument('--pop-prev', default=None,
    help='Population prevalence of binary phenotype (for conversion to liability scale).')

# Becomes:
parser.add_argument('--samp-prev', default=None, type=float,
    help='Sample prevalence of binary phenotype (for conversion to liability scale).')
parser.add_argument('--pop-prev', default=None, type=float,
    help='Population prevalence of binary phenotype (for conversion to liability scale).')
```

### 8. `ldsc_mod/sumstats.py` — bug in `_get_GxE_var_table`

**Problem:** this function incorrectly assumed `args.samp_prev`/`args.pop_prev` were lists
(`for i in args.samp_prev`), whereas after fix #7 they are single numbers
(`TypeError: 'float' object is not iterable`). The original code also contained a substantive
error: the p-value was multiplied by the liability-scale scaling factor, which is not
statistically correct (a p-value does not change when an estimate is rescaled).

**In the `_get_GxE_var_table` function, the `if` block right after `x['gxe_sumstats'] = [args.gxe_sumstats]`:**

```python
# Was:
if args.samp_prev is not None and \
        args.pop_prev is not None and \
        all((i is not None for i in args.samp_prev)) and \
        all((i is not None for it in args.pop_prev)):
    c = list(map(lambda x, y: reg.h2_obs_to_liab(1, x, y), args.samp_prev[1:], args.pop_prev[1:]))
    x['h2_gxe_liab'] = list(map(lambda x, y: x * y, c, [hsq.tot]))
    x['h2_gxe_liab_se'] = list(map(lambda x, y: x * y, c, [hsq.tot_se]))
    x['h2_gxe_liab_p'] = list(map(lambda x, y: x * y, c, [hsq.p]))
    x['h2_gxe_liab_int'] = list(map(lambda x, y: x * y, c, [hsq.intercept]))
    x['h2_gxe_liab_int_se'] = list(map(lambda x, y: x * y, c, [hsq.intercept_se]))
    x['h2_gxe_liab_int_p'] = list(map(lambda x, y: x * y, c, [hsq.intercept_p]))
    x = x.sort_values(by=['h2_gxe_liab_p'], na_position='last')

# Becomes:
if args.samp_prev is not None and args.pop_prev is not None:
    c = reg.h2_obs_to_liab(1, args.samp_prev, args.pop_prev)
    x['h2_gxe_liab'] = [hsq.tot * c]
    x['h2_gxe_liab_se'] = [hsq.tot_se * c]
    x['h2_gxe_liab_p'] = [hsq.p]
    x['h2_gxe_liab_int'] = [hsq.intercept]
    x['h2_gxe_liab_int_se'] = [hsq.intercept_se]
    x['h2_gxe_liab_int_p'] = [hsq.intercept_p]
    x = x.sort_values(by=['h2_gxe_liab_p'], na_position='last')
```

The `else` branch of this function (observed-scale output) was left unchanged.

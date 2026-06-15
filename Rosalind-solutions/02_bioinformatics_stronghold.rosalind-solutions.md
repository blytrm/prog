# Rosalind Problems + Answers

Solutions to bioinformatics problems from [Rosalind](http://rosalind.info/), implemented in both Python and R.

# ROSALIND Solutions – Python & R

## 1. Counting DNA Nucleotides

**Problem**  
Given: A DNA string `s` of length at most 1000 bp.  
Return: Four integers (separated by spaces) counting the respective number of occurrences of `A`, `C`, `G`, and `T` in `s`.

**Python**

```
with open("Counting DNA Nucleotides (1).txt", "r") as file:
    S = file.readline().strip()

bases = ("A", "C", "G", "T")
chr_count = {}

for base in S:
    chr_count[base] = chr_count.get(base, 0) + 1

counts = [chr_count.get(base, 0) for base in bases]
print(" ".join(str(count) for count in counts))
```

**R**

```
txt_path <- "Counting DNA Nucleotides (1).txt"
S <- readLines(txt_path)
bases <- strsplit(S, "")[][1]
count <- table(bases)
print(count)
```

---

## 2. Transcribing DNA into RNA

**Problem**  
Given: A DNA string `t` corresponding to a coding strand.  
Return: The transcribed RNA string `u` formed by replacing all `T` with `U`.

**Python**

```
with open("Transcribing DNA to RNA.txt", "r") as file:
    t = file.readline().strip()

transc = t.replace("T", "U")
print(transc)
```

**R**

```
txt <- "Transcribing DNA to RNA.txt"
t <- readLines(txt)
transc <- gsub(pattern = "T", replacement = "U", x = t)
print(transc)

# or
transc2 <- chartr(old = "T", new = "U", x = t)
print(transc2)
```

---

## 3. Complementing a Strand of DNA

**Problem**  
Given: A DNA string `s`.  
Return: The reverse complement of `s`.

**Python**

```
from Bio.Seq import Seq

with open("DNA Complement Strand.txt", "r") as file:
    s = file.readline().strip()

seq = Seq(s)
rev_comp = seq.reverse_complement()
print(str(rev_comp))
```

**R**

```
txt <- "DNA Complement Strand.txt"
s <- readLines(txt)[1]

rc <- function(x) {
  ss    <- unlist(strsplit(x, NULL))
  ssrev <- rev(ss)

  for (i in seq_along(ssrev)) {
    if (ssrev[i] == "T") {
      ssrev[i] <- "A"
    } else if (ssrev[i] == "A") {
      ssrev[i] <- "T"
    } else if (ssrev[i] == "G") {
      ssrev[i] <- "C"
    } else if (ssrev[i] == "C") {
      ssrev[i] <- "G"
    }
  }

  paste(ssrev, collapse = "")
}

rc_seq <- rc(s)
cat(rc_seq, "\n")
```

---

## 4. Rabbits and Recurrence Relations

**Problem**  
Given: Positive integers `n` and `k`. Start with 1 pair of rabbits; each reproduction‑age pair produces `k` new pairs each month.  
Return: The total number of rabbit pairs after `n` months.

**R**

```
rabbit_pairs_month <- c()
n <- 28  # months
k <- 2   # new pairs per reproducing pair

for (i in 1:n) {
  if (i == 1) {
    rabbit_pairs_month[i] <- 1
  } else if (i == 2) {
    rabbit_pairs_month[i] <- 1
  } else {
    rabbit_pairs_month[i] <- rabbit_pairs_month[i - 1] + k * rabbit_pairs_month[i - 2]
  }
}

rabbit_pairs_month[length(rabbit_pairs_month)]
```

**Python**

```
n, k = 28, 2

def mortal_r(n, k):
    if n == 1 or n == 2:
        return 1

    F_prev2 = 1  # F1
    F_prev1 = 1  # F2

    for _ in range(3, n + 1):
        F_curr = F_prev1 + k * F_prev2
        F_prev2, F_prev1 = F_prev1, F_curr

    return F_curr

print(mortal_r(n, k))


def fib(n, k):
    previous1, previous2 = 1, 1
    for i in range(2, n):
        current = previous1 + k * previous2
        previous2 = previous1
        previous1 = current
    return current
```

---

## 5. GC Content

**Problem**  
Given: At most 10 DNA strings in FASTA format (length ≤ 1 kbp).  
Return: The ID of the string with the highest GC content, followed by its GC percentage (to within 0.001).

**Python**

```
def readFile(filePath):
    with open(filePath, "r") as f:
        return [l.strip() for l in f.readlines()]

ffile = readFile("GC Content Analysis.txt")

fdict = {}
flabel = ""

for line in ffile:
    if ">" in line:
        flabel = line
        fdict[flabel] = ""
    else:
        fdict[flabel] += line

def gc_content(seq):
    return round(((seq.count("C") + seq.count("G")) / len(seq) * 100), 6)

resultsDict = {key: gc_content(value) for (key, value) in fdict.items()}

maxGC_key = max(resultsDict, key=resultsDict.get)
print(maxGC_key)
print(f"{maxGC_key[1:]}\n{resultsDict[maxGC_key]}")
```

**R**

```
file_path <- "GC Content Analysis.txt"
lines <- readLines(file_path)

seq_list  <- list()
current_id <- NULL

for (line in lines) {
  if (startsWith(line, ">")) {
    current_id <- substring(line, 2)
    seq_list[[current_id]] <- ""
  } else {
    seq_list[[current_id]] <- paste0(seq_list[[current_id]], line)
  }
}

gc_content <- function(seq) {
  chars <- strsplit(seq, "")[][1]
  g <- sum(chars == "G")
  c <- sum(chars == "C")
  gc <- (g + c) / length(chars) * 100
  round(gc, 6)
}

gc_values <- sapply(seq_list, gc_content)

max_id <- names(gc_values)[which.max(gc_values)]
max_gc <- gc_values[max_id]

cat(max_id, "\n", max_gc, "\n", sep = "")
```

---

## 6. Counting Point Mutations

**Problem**  
Given: Two DNA strings `s` and `t` of equal length.  
Return: The Hamming distance between them (number of differing positions).

**Python**

```
def readFile(filePath):
    with open(filePath, "r") as file:
        return [l.strip() for l in file.readlines()]

ffile = readFile("Counting Point Mutations (1).txt")
s = ffile
t = ffile[1]

difrncs = 0

if len(s) == len(t):
    for char1, char2 in zip(s, t):
        if char1 != char2:
            difrncs += 1

print(difrncs)
```

**R**

```
file_path <- "Counting Point Mutations (1).txt"
lines <- readLines(file_path)
s <- lines[1]
t <- lines[2]

hamming_distance <- function(a, b) {
  count <- 0
  for (i in seq_len(nchar(a))) {
    if (substr(a, i, i) != substr(b, i, i)) {
      count <- count + 1
    }
  }
  count
}

d <- hamming_distance(s, t)
cat(d, "\n")
```

---

## 7. Mendel’s First Law

**Problem**  
Given: Three positive integers `k`, `m`, `n` representing a population of `k + m + n` organisms:  
- `k` homozygous dominant (`AA`),  
- `m` heterozygous (`Aa`),  
- `n` homozygous recessive (`aa`).  

Return: The probability that two randomly selected organisms will produce an offspring with at least one dominant allele.

**Python**

```
with open("Mendel's First Law (1).txt", "r") as f:
    line = f.readline().strip()

values = line.split()
k, m, n = map(int, values)
N = k + m + n

P_AA_AA = (k / N) * ((k - 1) / (N - 1))
P_AA_Aa = (k / N) * (m / (N - 1)) + (m / N) * (k / (N - 1))
P_AA_aa = (k / N) * (n / (N - 1)) + (n / N) * (k / (N - 1))
P_Aa_Aa = (m / N) * ((m - 1) / (N - 1))
P_Aa_aa = (m / N) * (n / (N - 1)) + (n / N) * (m / (N - 1))
P_aa_aa = (n / N) * ((n - 1) / (N - 1))

Pr_dom = (
    P_AA_AA * 1.0 +
    P_AA_Aa * 1.0 +
    P_AA_aa * 1.0 +
    P_Aa_Aa * 0.75 +
    P_Aa_aa * 0.5 +
    P_aa_aa * 0.0
)

print(f"{Pr_dom:.5f}")
```

**R**

```
line <- readLines("Mendel's First Law (1).txt")
line <- trimws(line)

values <- strsplit(line, " ")[][1]
nums <- as.numeric(values)
k <- nums[1]
m <- nums[2]
n <- nums[3]
N <- k + m + n

P_AA_AA <- (k / N) * ((k - 1) / (N - 1))
P_AA_Aa <- (k / N) * (m / (N - 1)) + (m / N) * (k / (N - 1))
P_AA_aa <- (k / N) * (n / (N - 1)) + (n / N) * (k / (N - 1))
P_Aa_Aa <- (m / N) * ((m - 1) / (N - 1))
P_Aa_aa <- (m / N) * (n / (N - 1)) + (n / N) * (m / (N - 1))
P_aa_aa <- (n / N) * ((n - 1) / (N - 1))

Pr_dom <-
  P_AA_AA * 1.0 +
  P_AA_Aa * 1.0 +
  P_AA_aa * 1.0 +
  P_Aa_Aa * 0.75 +
  P_Aa_aa * 0.5 +
  P_aa_aa * 0.0

cat(sprintf("%.5f\n", Pr_dom))
```

---

## 8. Translating RNA into Protein

**Problem**  
Given: An RNA string `s` corresponding to an mRNA strand.  
Return: The protein string encoded by `s`.

**Python**

```
from Bio.Seq import Seq

def readFile(filePath):
    with open(filePath, "r") as f:
        content = f.readline().strip()
    return content

s = readFile("RNA to Protein (2).txt")

rna = Seq(s)
protein = rna.translate(to_stop=True)
print(protein)
```

**R**

```
library(Biostrings)
library(bioseq)

file_path <- "RNA to Protein.txt"
s <- trimws(readLines(file_path))[1]

# Biostrings
rna_bs <- RNAString(s)
protein_bs <- translate(rna_bs)
print(protein_bs)

# bioseq
rna_bq <- rna(s)
protein_bq <- seq_translate(rna_bq)
print(protein_bq)
```

---

## 9. Finding a Motif in DNA

**Problem**  
Given: Two strings `s` and `t`, with `t` a substring of `s`.  
Return: All locations (1‑based) where `t` appears in `s`.

**Python – simple double loop**

```
def read_file(path):
    with open(path) as f:
        lines = [l.strip() for l in f]
    return lines, lines[1]

s, t = read_file("Finding a Motif in DNA.txt")

n = len(s)
m = len(t)
positions = []

for i in range(n - m + 1):
    has_match = True
    for j in range(m):
        if s[i + j] != t[j]:
            has_match = False
            break
    if has_match:
        positions.append(i + 1)

print(" ".join(map(str, positions)))
```

**R (Biostrings)**

```
library(Biostrings)

lines <- readLines("Finding a Motif in DNA.txt")
s <- lines[1]
t <- lines[2]

dna   <- DNAString(s)
motif <- DNAString(t)

hits <- matchPattern(motif, dna)
pos  <- start(hits)

cat(paste(pos, collapse = " "), "\n")
```

---

## 10. Consensus and Profile

**Problem**  
Given: A collection of DNA strings in FASTA format, all of the same length.  
Return: A consensus string and a profile matrix (counts of A, C, G, T at each position).

**Python**

```
from Bio import motifs
from Bio.Seq import Seq
from Bio import SeqIO
from io import StringIO

def readFile(filePath):
    with open(filePath, "r") as f:
        return [l.strip() for l in f.readlines()]

lines = readFile("Consensus Profile Rosalind.txt")
fasta_content = "\n".join(lines)

seq_strings = []
for record in SeqIO.parse(StringIO(fasta_content), "fasta"):
    seq_strings.append(str(record.seq))

def profile_matrix(sequences):
    seq_length = len(sequences)
    profile = {
        "A":  * seq_length,
        "C":  * seq_length,
        "G":  * seq_length,
        "T":  * seq_length,
    }
    for seq in sequences:
        for i, n in enumerate(seq):
            if n in profile:
                profile[n][i] += 1
    return profile

def consensus(profile):
    cols = len(next(iter(profile.values())))
    con = []
    for i in range(cols):
        col_count = {nt: profile[nt][i] for nt in "ACGT"}
        con.append(max(col_count, key=col_count.get))
    return "".join(con)

profile = profile_matrix(seq_strings)
con = consensus(profile)

print(con)
for nt in "ACGT":
    print(f"{nt}: {' '.join(map(str, profile[nt]))}")

seq_obj = [Seq(s) for s in seq_strings]
m = motifs.create(seq_obj)
print(m.consensus)
print(m.counts)
```

**R**

```
library(Biostrings)
library(seqinr)

fasta_file <- "Consensus Proile Rosalind (2).txt"
seq_list   <- read.fasta(file = fasta_file, seqtype = "DNA", as.string = FALSE)

seq_strings <- vapply(seq_list, paste0, collapse = "", FUN.VALUE = character(1))

prof_matrix <- function(sequences) {
  seq_length <- nchar(sequences)[1]
  profile <- list(
    A = integer(seq_length),
    C = integer(seq_length),
    G = integer(seq_length),
    T = integer(seq_length)
  )
  for (seq in sequences) {
    nts <- strsplit(seq, "")[][1]
    for (i in seq_along(nts)) {
      n <- nts[i]
      if (n %in% names(profile)) {
        profile[[n]][i] <- profile[[n]][i] + 1
      }
    }
  }
  profile
}

profile <- prof_matrix(seq_strings)

consensus_matrix_acgt <- function(profile) {
  cols <- length(profile[])[1]
  consensus <- character(cols)
  for (i in seq_len(cols)) {
    col_counts <- sapply(c("A", "C", "G", "T"), function(nt) profile[[nt]][i])
    consensus[i] <- names(which.max(col_counts))
  }
  paste0(consensus, collapse = "")
}

cons <- consensus_matrix_acgt(profile)
cat(cons, "\n")

for (nt in c("A", "C", "G", "T")) {
  cat(nt, ":", paste(profile[[nt]], collapse = " "), "\n")
}
```
---

## 11. Mortal Fibonacci Rabbits

**Problem**  
Given: Positive integers `n` and `m`.  
Return: The total number of rabbit pairs that will remain after the `n`-th month if all rabbits live for `m` months.

**Python**

```python
def read(path):
    with open(path, 'r') as f:
        return [l.strip() for l in f.readlines()]

fi = read('Fibonacci Rabbits (2).txt')
n, m = map(int, fi[0].split())
k = 1  # each pair produces 1 new pair

def mortal_rabbits(n, m, k):
    ages = [0] * m  # age groups: index 0 = newborns
    ages[0] = 1     # month 1: 1 newborn pair
    
    for month in range(1, n):
        breeders = sum(ages[1:])  # all except newborns
        newborns = breeders * k
        ages = [newborns] + ages[:-1]  # shift ages, oldest die
        
    return sum(ages)

print(mortal_rabbits(n, m, k))
```

**R**

```r
library(gmp)  # for large integers

fi <- readLines('Fibonacci Rabbits (2).txt')
n <- as.integer(strsplit(fi[1], " ")[[1]][1])
m <- as.integer(strsplit(fi[1], " ")[[1]][2])
k <- as.bigz(1)

mortal_rabbits <- function(n, m, k = 1) {
  ages <- c(as.bigz(1), rep(as.bigz(0), m - 1))
  
  for (i in 1:(n - 1)) {
    newborns <- k * sum(ages[-1])  # sum of all except newborns
    ages <- c(newborns, ages[-length(ages)])  # newborns + all except oldest
  }
  
  return(sum(ages))
}

cat(as.character(mortal_rabbits(n, m, k)))
```

---

## 12. Overlap Graphs

**Problem**  
Given: A collection of DNA strings in FASTA format (length ≤ 10 kbp).  
Return: The adjacency list corresponding to O₃ (k=3). You may return edges in any order.

**Python**

```python
def read_fasta(path):
    with open(path, 'r') as f:
        lines = [l.strip() for l in f.readlines()]
    
    sequences = {}
    current_label = ""
    
    for line in lines:
        if line.startswith(">"):
            current_label = line[1:]
            sequences[current_label] = ""
        else:
            sequences[current_label] += line
    
    return sequences

sequences = read_fasta('Overlap Graphs.txt')
k = 3

for s_label, s_dna in sequences.items():
    for t_label, t_dna in sequences.items():
        if s_label != t_label and s_dna[-k:] == t_dna[:k]:
            print(f"{s_label} {t_label}")
```

**R**

```r
read_fasta <- function(path) {
  lines <- readLines(path)
  sequences <- list()
  current_label <- ""
  
  for (line in lines) {
    if (startsWith(line, ">")) {
      current_label <- substring(line, 2)
      sequences[[current_label]] <- ""
    } else {
      sequences[[current_label]] <- paste0(sequences[[current_label]], line)
    }
  }
  
  return(sequences)
}

sequences <- read_fasta('Overlap Graphs (1).txt')
k <- 3
labels <- names(sequences)

for (s_name in labels) {
  for (t_name in labels) {
    if (s_name != t_name) {
      s_dna <- sequences[[s_name]]
      t_dna <- sequences[[t_name]]
      
      suffix <- substr(s_dna, nchar(s_dna) - k + 1, nchar(s_dna))
      prefix <- substr(t_dna, 1, k)
      
      if (suffix == prefix) {
        cat(s_name, t_name, "\n")
      }
    }
  }
}
```

**R (with stringr)**

```r
library(Biostrings)
library(stringr)

dna_obj <- readDNAStringSet("Overlap Graphs (1).txt")
dna_seqs <- as.character(dna_obj)
k <- 3
ids <- names(dna_seqs)

for (s in ids) {
  for (t in ids) {
    if (s != t && str_sub(dna_seqs[[s]], -k) == str_sub(dna_seqs[[t]], 1, k)) {
      cat(s, t, "\n")
    }
  }
}
```

---

## 13. Calculating Expected Offspring

**Problem**  
Given: Six nonnegative integers corresponding to the number of couples with genotypes:  
1. AA-AA  
2. AA-Aa  
3. AA-aa  
4. Aa-Aa  
5. Aa-aa  
6. aa-aa  

Return: The expected number of offspring displaying the dominant phenotype in the next generation, assuming each couple has exactly two offspring.

**Python (NumPy)**

```python
import numpy as np

def read_expected(path):
    with open(path, 'r') as f:
        data = list(map(int, f.read().split()))
    
    couples = np.array(data)
    probs = [1.0, 1.0, 1.0, 0.75, 0.50, 0.0]
    offspring_per_couple = 2
    
    expected = np.sum(couples * probs * offspring_per_couple)
    return expected

result = read_expected('Calculating Expected Offspring.txt')
print(result)
```

**Python (direct calculation)**

```python
def calculate_expected_offspring(path):
    with open(path, 'r') as f:
        a, b, c, d, e, f = map(int, f.read().split())
    
    expected = 2 * (a + b + c) + 1.5 * d + 1.0 * e
    return expected

print(calculate_expected_offspring('Calculating Expected Offspring.txt'))
```

**R**

```r
read_expected <- function(path) {
  couples <- scan(path, what = integer(), quiet = TRUE)[1:6]
  probs <- c(1.0, 1.0, 1.0, 0.75, 0.50, 0.0)
  offspring_per_couple <- 2
  
  expected <- sum(couples * probs * offspring_per_couple)
  return(expected)
}

result <- read_expected('Calculating Expected Offspring.txt')
cat(result, "\n")
```

---

## 14. Finding a Shared Motif

**Problem**  
Given: A collection of k DNA strings in FASTA format (length ≤ 1 kbp each).  
Return: A longest common substring of the collection.

**Python**

```python
def read_fasta(path):
    sequences = []
    current_dna = ""
    
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith(">"):
                if current_dna:
                    sequences.append(current_dna)
                current_dna = ""
            else:
                current_dna += line
        
        if current_dna:
            sequences.append(current_dna)
    
    return sequences

def longest_common_substring(sequences):
    if not sequences:
        return ""
    
    reference = min(sequences, key=len)
    ref_len = len(reference)
    
    for window_size in range(ref_len, 0, -1):
        for start in range(ref_len - window_size + 1):
            substring = reference[start:start + window_size]
            if all(substring in seq for seq in sequences):
                return substring
    
    return ""

sequences = read_fasta('Shared Motif LCSM.txt')
result = longest_common_substring(sequences)
print(result)
```

**R**

```r
read_fasta_manual <- function(path) {
  lines <- readLines(path)
  sequences <- character()
  current_dna <- ""
  
  for (line in lines) {
    line <- trimws(line)
    if (nchar(line) == 0) next
    
    if (startsWith(line, ">")) {
      if (nchar(current_dna) > 0) {
        sequences <- c(sequences, current_dna)
      }
      current_dna <- ""
    } else {
      current_dna <- paste0(current_dna, line)
    }
  }
  
  if (nchar(current_dna) > 0) {
    sequences <- c(sequences, current_dna)
  }
  
  return(sequences)
}

longest_common_substring <- function(seqs) {
  seqs <- toupper(seqs)
  idx_min <- which.min(nchar(seqs))
  ref <- seqs[idx_min]
  others <- seqs[-idx_min]
  ref_len <- nchar(ref)
  
  for (window_size in ref_len:1) {
    starts <- 1:(ref_len - window_size + 1)
    for (start in starts) {
      sub <- substr(ref, start, start + window_size - 1)
      if (all(vapply(others, grepl, logical(1), 
                     pattern = sub, fixed = TRUE))) {
        return(sub)
      }
    }
  }
  
  return("")
}

sequences <- read_fasta_manual('Shared Motif LCSM (1).txt')
result <- longest_common_substring(sequences)
cat(result, "\n")
```

---

## 15. Independent Alleles

**Problem**  
Given: Two positive integers k and N.  
Return: The probability that at least N AaBb organisms will belong to the k-th generation (starting with AaBb in generation 0). Assume each organism always mates with an AaBb organism.

**Python**

```python
import math

def read_alleles(path):
    with open(path, 'r') as f:
        k, n = map(int, f.readline().split())
    return k, n

k, n = read_alleles('Independent Alleles.txt')

def probability_at_least_n(k, n):
    m = 2 ** k  # total organisms in k-th generation
    p = 0.25    # probability of AaBb in offspring
    q = 1 - p
    
    total_prob = 0
    for i in range(n, m + 1):
        comb = math.comb(m, i)
        prob_i = comb * (p ** i) * (q ** (m - i))
        total_prob += prob_i
    
    return total_prob

print(probability_at_least_n(k, n))

# Alternative using scipy
# from scipy.stats import binom
# print(binom.sf(n - 1, 2**k, 0.25))
```

**R**

```r
data <- readLines('Independent Alleles.txt')
k <- as.integer(strsplit(data[1], " ")[[1]][1])
n <- as.integer(strsplit(data[1], " ")[[1]][2])

m <- 2^k
p <- 0.25

# P(X >= n) = 1 - P(X <= n-1)
prob <- pbinom(n - 1, size = m, prob = p, lower.tail = FALSE)
cat(prob, "\n")
```

---

## 16. Finding a Protein Motif

**Problem**  
Given: At most 15 UniProt Protein Database access IDs.  
Return: For each protein possessing the N-glycosylation motif, output its given access ID followed by a list of locations in the protein string where the motif can be found.

**Python**

```python
import requests
import re

def read(path):
    with open(path, 'r') as f:
        return [l.strip() for l in f.readlines()]

pattern = r"(?=N[^P][ST][^P])"  # Lookahead without capturing group
lines = read("/content/Finding Protein Motif (5).txt")

for line in lines:
    org_q = line.split('_')[0]
    url = 'http://www.uniprot.org/uniprotkb/' + org_q + '.fasta'
    response = requests.get(url)
    if response.status_code != 200:
        continue
    
    fasta_lines = response.text.strip().split('\n')
    seq = ''.join(fasta_lines[1:])
    
    matches = re.finditer(pattern, seq)
    pos = [m.start() + 1 for m in matches]  
    
    if pos:
        print(line)
        print(" ".join(map(str, pos)))
```

**R**

```r
library(httr)
library(stringr)

pattern <- "(?=N[^P][ST][^P])"  

lines <- readLines("/content/Finding Protein Motif (5).txt")

for (line in lines) {
  org_q <- strsplit(line, "_")[[1]][1]

  url <- paste0("https://rest.uniprot.org/uniprotkb/", org_q, ".fasta")
  res <- GET(url)
  if (status_code(res) != 200) next

  fasta_text <- content(res, as = "text", encoding = "UTF-8")
  fasta_lines <- strsplit(trimws(fasta_text), "\n")[[1]]
  seq <- paste(fasta_lines[-1], collapse = "")   # drop header

  # str_locate_all returns start/end indices for *each character* of the match
  # but for lookahead the match is zero-width, so start == end == position *before* N
  m <- str_locate_all(seq, pattern)[[1]]
  if (nrow(m) == 0) next

  pos <- m[, "start"] 

  cat(line, "\n")
  cat(paste(pos, collapse = " "), "\n")
}
```

---

## 17. Inferring mRNA from Protein

**Problem**  
Given: A protein string of length at most 1000 aa.  
Return: The total number of different RNA strings from which the protein could have been translated, modulo 1,000,000. (Don't neglect the importance of the stop codon in protein translation.)

**Python**

```python
def read(path):
    with open(path, 'r') as f:
        return [l.strip() for l in f.readlines()]

string = read("/content/Inferring mRNA from Protein.txt")

cod = {
    'A': 4, 'C': 2, 'D': 2, 'E': 2, 'F': 2,
    'G': 4, 'H': 2, 'I': 3, 'K': 2, 'L': 6,
    'M': 1, 'N': 2, 'P': 4, 'Q': 2, 'R': 6,
    'S': 6, 'T': 4, 'V': 4, 'W': 1, 'Y': 2,
}  # codons for each aa

def count_rna(string, mod=1000000):
    result = 1
    protein_sequence = string[0] 
    for i in protein_sequence:
        result = (result * cod[i]) % mod
    result = (result * 3) % mod  # multiply by stop codon possibilities
    return result

final_result = count_rna(string)
print(final_result)
```

**R**

```r
protein_input <- readLines("/content/Inferring mRNA from Protein.txt")

cod <- c(
    A=4, C=2, D=2, E=2, F=2,
    G=4, H=2, I=3, K=2, L=6,
    M=1, N=2, P=4, Q=2, R=6,
    S=6, T=4, V=4, W=1, Y=2
)

count_rna <- function(protein_input_str, mod=1000000) {
    result <- 1
    for (aa in strsplit(protein_input_str, "")[[1]]) {
        result <- (result * cod[aa]) %% mod
    }
    result <- (result * 3) %% mod  # multiply by stop codon possibilities
    return(result)
}

final <- count_rna(protein_input[1])
print(final)
```

---

## 18. Open Reading Frames

**Problem**  
Given: A DNA string in FASTA format.  
Return: Every candidate protein string that can be translated from ORFs.

**Python**

```python
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO

def read(path):
    with open(path, 'r') as f:
        return [l.strip() for l in f.readlines()]

l = read("/content/Open Reading Frames (5).txt")
content = "\n".join(l)

seq_str = []
for i in SeqIO.parse(StringIO(content), "fasta"):
    seq_str.append(str(i.seq))

idna = "".join(seq_str)
dna = Seq(idna)
    
def count_orf(seq):
    orf = []
    
    for i in range(len(seq) - 2):
        if str(seq[i:i+3]) == "ATG":
            for j in range(i + 3, len(seq) - 2, 3):
                stop = ["TAG", "TGA", "TAA"]
                
                if str(seq[j:j+3]) in stop:
                    orf.append(str(Seq(str(seq[i:j])).translate()))
                    break
    return orf

p = count_orf(dna)
rev_p = count_orf(dna.reverse_complement())

final = set(p + rev_p)
for c in final:
    print(c)
```

**R**

```r
library(Biostrings)

dna_set <- readDNAStringSet("~/Downloads/Open Reading Frames.txt")
dna <- dna_set[[1]]  # Get the first sequence

find_proteins <- function(sequence) {
  proteins <- c()
  seq_str <- as.character(sequence)
  
  # Search all 3 frames on this strand
  for (i in 1:(nchar(seq_str) - 2)) {
    # Look for Start Codon (ATG)
    if (substr(seq_str, i, i + 2) == "ATG") {
      # Translate from this start point to the end of the sequence
      aa_seq <- as.character(translate(DNAString(substr(seq_str, i, nchar(seq_str)))))
      
      # If there is a stop codon (marked by *), extract until the first *
      if (grepl("\\*", aa_seq)) {
        protein <- sub("\\*.*", "", aa_seq)
        proteins <- c(proteins, protein)
      }
    }
  }
  return(proteins)
}

forward_prots <- find_proteins(dna)
reverse_prots <- find_proteins(reverseComplement(dna))

final_prots <- unique(c(forward_prots, reverse_prots))
cat(final_prots, sep = "\n")
```

---

## 19. Enumerating Gene Orders

**Problem**  
Given: A positive integer n < 8.  
Return: The total number of permutations of length n, and a list of all such permutations.

**Python**

```python
import math
import itertools

with open("Enumerating Gene Orders.txt", 'r') as f:
    n = int(f.read().strip())

def get_perm(n):
    elements = list(range(1, n + 1))
    perms = list(itertools.permutations(elements))
    print(len(perms))
    for p in perms:
        print(*(p))
        # * = splat operator = unpacking
        # print(p) = (1, 2, 3)
        # print(*p) = 1 2 3

get_perm(n)
```

**R**

```r
library(gtools)

n <- as.integer(readLines("path", n = 1)) # single integer
cat(factorial(n), "\n")
perms <- permutations(n = n, r = n, v = 1:n)
write.table(perms, col.names = FALSE, row.names = FALSE)
```

---

## 20. Calculating Protein Mass

**Problem**  
Given: A protein string of amino acids.  
Return: The weight of the protein.

**Python**

```python
def read(file):
    with open(file, 'r') as f:
        return f.read().strip()

def weight(string, weight_table):
    return sum(weight_table[aa] for aa in string)

def main():
    file = "/content/Protein Mass Calculation.txt"
    string = read(file)
    
    amino_acid_masses = {
        'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 
        'F': 147.06841, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406, 
        'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293, 
        'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203, 
        'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333
    }
    
    result = weight(string, amino_acid_masses)
    
    print(result)

if __name__ == "__main__":
    main()
```

**R**

```r
protein_str <- readLines("~/Downloads/Protein Mass Calculation.txt", warn = FALSE)

mass_table <- c(
  A = 71.03711, C = 103.00919, D = 115.02694, E = 129.04259, 
  F = 147.06841, G = 57.02146, H = 137.05891, I = 113.08406, 
  K = 128.09496, L = 113.08406, M = 131.04049, N = 114.04293, 
  P = 97.05276, Q = 128.05858, R = 156.10111, S = 87.03203, 
  T = 101.04768, V = 99.06841, W = 186.07931, Y = 163.06333
)

amino_acids <- strsplit(protein_str, "")[[1]]
total_mass <- sum(mass_table[amino_acids])
cat(round(total_mass, 3), "\n")
```

---

## 21. Locating Restriction Sites

**Problem**  
Given: A DNA string in FASTA format.  
Return: The position and length of every reverse palindrome in the string with length between 4 & 12.

**Python**

```python
def restriction_enzymes(file):
    with open(file, 'r') as f:
        dna = "".join(f.read().splitlines()[1:])
    
    trans = str.maketrans("ATGC", "TACG")
    results = []
    
    for i in range(len(dna)):
        for j in range(4, 13):
            if i + j > len(dna):
                break
            
            substring = dna[i:i+j]
            rev_comp = substring[::-1].translate(trans)
            
            if substring == rev_comp:
                results.append((i + 1, j))
    return results
      
for position, length in restriction_enzymes("/content/Restriction Sites Revp.txt"):
    print(f"{position} {length}")
```

**R**

```r
library(Biostrings)

restriction_enzymes <- function(file) {
    dna <- readDNAStringSet(file)[[1]]
    palindromes <- findPalindromes(dna, min.armlength = 2, max.looplength = 0, max.mismatch = 0)
    palindromes <- palindromes[width(palindromes) <= 12]
    results <- data.frame(
        position = start(palindromes),
        length = width(palindromes)
    )
    return(results)
}

results <- restriction_enzymes("~/Downloads/Restriction Sites Revp.txt")
cat(paste(results$position, results$length, collapse = "\n"))
```

I'll add your new solutions to the document in the same format. Let me organize them properly.
 Now I'll create the complete updated document with all your new solutions added in the proper format. Let me compile this:
 Now I have enough information. Let me create the complete updated document with all your new solutions added in the proper Rosalind format. I'll compile this as a comprehensive response:

Based on your solutions and the established format, here are your new Rosalind solutions **22-31** added to the document:

---

## 22. RNA Splicing

**Problem**  
Given: A DNA string `s` (of length at most 1 kbp) and a collection of substrings of `s` acting as introns. All strings are given in FASTA format.  
Return: A protein string resulting from transcribing and translating the exons of `s`. 

**Python**

```python
from Bio import SeqIO
from Bio.Seq import Seq
import re

def readfasta(file):
    records = list(SeqIO.parse(file, "fasta"))
    dna_string = str(records[0].seq)
    intron_list = [str(r.seq) for r in records[1:]]
    return dna_string, intron_list

dna, introns = readfasta("/content/RNA Splicing (2).txt")

pattern = "|".join(map(re.escape, introns))
spliced_dna = re.sub(pattern, "", dna)

protein = Seq(spliced_dna).translate(to_stop=True)
print(protein)
```

**R**

```r
library(Biostrings)
library(stringr)

fasta <- readDNAStringSet("~/Downloads/RNA Splicing (2).txt")
dna <- as.character(fasta[1])
introns <- as.character(fasta[-1])

for(intron in introns) {
  dna <- gsub(intron, "", dna, fixed = TRUE)
}

protein <- translate(DNAString(dna))
cat(as.character(protein))
```

---

## 23. Enumerating k-mers Lexicographically

**Problem**  
Given: A collection of at most 10 symbols defining an ordered alphabet, and a positive integer `n` (≤ 4).  
Return: All strings of length `n` that can be formed from the alphabet, ordered lexicographically. 

**Python**

```python
from itertools import product

s = ["A", "B", "C", "D", "E"]
n = 4

comb = [''.join(i) for i in product(s, repeat=n)]
sorted_comb = sorted(comb)
print(*sorted_comb, sep='\n')
```

**R**

```r
s <- c("A", "B", "C", "D", "E")
n <- 4

grid <- expand.grid(rep(list(s), n))
comb <- do.call(paste0, grid)
sorted_comb <- sort(comb)
cat(sorted_comb, sep = "\n")
```

---

## 24. Longest Increasing Subsequence

**Problem**  
Given: A positive integer `n` (≤ 10,000) followed by a permutation `π` of length `n`.  
Return: A longest increasing subsequence of `π`, followed by a longest decreasing subsequence of `π`. 

**Python**

```python
from bisect import bisect_left

def read(file):
    with open(file, 'r') as f:
        line1 = f.readline().strip()
        n = int(line1)
        permut = [int(x) for x in f.read().split()]
        return n, permut

def perm_sub(arr, increasing=True):
    n = len(arr)
    tails = []
    tail_idx = []
    prev = [-1] * n
    
    for i, x in enumerate(arr):
        val = x if increasing else -x
        pos = bisect_left(tails, val)
        
        if pos > 0:
            prev[i] = tail_idx[pos - 1]
        
        if pos == len(tails):
            tails.append(val)
            tail_idx.append(i)
        else:
            tails[pos] = val
            tail_idx[pos] = i
    
    result, curr = [], tail_idx[-1]
    while curr != -1:
        result.append(arr[curr])
        curr = prev[curr]
    return result[::-1]

n, permut = read("/content/Longest Increasing Subsequence (4).txt")

LIS = perm_sub(permut, increasing=True)
LDS = perm_sub(permut, increasing=False)

print(*LIS)
print(*LDS)
```

---

## 25. Genome Assembly as Shortest Superstring

**Problem**  
Given: At most 50 DNA strings of approximately equal length (≤ 1 kbp) in FASTA format representing reads from a single linear chromosome.  
Return: A shortest superstring containing all the given strings (reconstructed chromosome). 

**Python**

```python
from Bio import SeqIO
from Bio.Seq import Seq

def get_overlap(s1, s2, min_overlap_required):
    """Return overlap length where s1's suffix matches s2's prefix."""
    for i in range(min(len(s1), len(s2)), min_overlap_required - 1, -1):
        if s1.endswith(s2[:i]):
            return i
    return 0

def assemble_genome(reads):
    """Greedy assembly merging pairs with maximum overlap."""
    reads = reads.copy()
    min_overlap_required = min(len(r) for r in reads) // 2 + 1

    while len(reads) > 1:
        best_overlap = 0
        best_i, best_j = 0, 0

        for i in range(len(reads)):
            for j in range(len(reads)):
                if i == j:
                    continue
                overlap = get_overlap(reads[i], reads[j], min_overlap_required)
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_i, best_j = i, j

        if best_overlap == 0:
            raise ValueError(f"No overlap found among {len(reads)} remaining reads")

        merged = reads[best_i] + reads[best_j][best_overlap:]
        
        if best_i > best_j:
            reads.pop(best_i)
            reads.pop(best_j)
        else:
            reads.pop(best_j)
            reads.pop(best_i)
        reads.append(merged)
    
    return reads[0]

def read_reads(file):
    return [str(record.seq) for record in SeqIO.parse(file, "fasta")]

def main():
    file = "/Users/billytrim/Downloads/Genome Assembly Shortest Superstring (1).txt"
    reads = read_reads(file)
    result = assemble_genome(reads)
    print(result)

if __name__ == "__main__":
    main()
```

---

## 26. Perfect Matchings and RNA Secondary Structures

**Problem**  
Given: An RNA string `s` (≤ 80 bp) with equal A/U and C/G counts.  
Return: The total possible number of perfect matchings of basepair edges in the bonding graph of `s`. 

**Python**

```python
from Bio import SeqIO
from Bio.Seq import Seq

rna = [str(record.seq) for record in SeqIO.parse("/content/Perfect Matchings RNA Structures (3).txt", "fasta")]

def perf_match(rna):
    n_au = rna.count('A')
    n_gc = rna.count('G')
    
    def fact(n):
        if n <= 1:
            return 1
        result = 1
        for i in range(2, n + 1):
            result *= i
        return result
    
    return fact(n_au) * fact(n_gc)

for r in rna:
    print(perf_match(r))
```

---

## 27. Partial Permutations

**Problem**  
Given: Positive integers `n` and `k` (100 ≥ n > 0, 10 ≥ k > 0).  
Return: The total number of partial permutations P(n, k), modulo 1,000,000. 

**Python**

```python
with open("/Users/billytrim/Desktop/2026/learning/rosalind/Partial Permutations.txt", 'r') as file:
    n, k = map(int, file.read().split())

def partial_perm(n, k, mod=1000000):
    result = 1
    for i in range(n, n - k, -1):
        result = (result * i) % mod
    return result

print(partial_perm(n, k))
```

---

## 28. Introduction to Random Strings

**Problem**  
Given: A DNA string `s` (≤ 100 bp) and an array `A` of at most 20 numbers between 0 and 1.  
Return: An array `B` where B[k] represents the common logarithm of the probability that a random string with GC-content A[k] matches `s`. 

**Python**

```python
import math

def dnamatch(string, probabilities):
    A_count = string.count("A")
    T_count = string.count("T")
    G_count = string.count("G")
    C_count = string.count("C")
    
    result = []
    for comp in probabilities:
        G_cont = comp / 2
        C_cont = comp / 2
        A_cont = (1 - comp) / 2
        T_cont = (1 - comp) / 2
        
        prob_product = ((A_cont)**A_count) * ((T_cont)**T_count) * ((G_cont)**G_count) * ((C_cont)**C_count)
        prob_match = math.log10(prob_product)
        result.append(round(prob_match, 3))
    
    return result

def readin(file):
    with open(file, 'r') as f:
        content = f.readlines()
        dna = str(content[0].strip())
        gc_content = [float(x) for x in content[1].split()]
    return dna, gc_content

def main():
    file = "/Users/billytrim/Downloads/Introduction to Random Strings (1).txt"
    string, probabilities = readin(file)
    answ = dnamatch(string, probabilities)
    print(*answ)

if __name__ == "__main__":
    main()
```

---

## 29. Enumerating Oriented Gene Orderings

**Problem**  
Given: A positive integer `n` ≤ 6.  
Return: The total number of signed permutations of length `n`, followed by a list of all such permutations. 

**Python**

```python
class Solution:
    def __init__(self):
        self.n = None
    
    def read(self, file):
        with open(file, 'r') as f:
            self.n = int(f.read().strip())
        return self.n
            
    def factorial(self):
        if self.n is None:
            return None
        ni = self.n
        res = 1
        while ni > 1:
            res = ni * res
            ni -= 1
        return res
            
    def sign_perm(self):           
        stack = []
        output = []
        used = set()

        def backtrack():
            if len(stack) == self.n:
                output.append(stack.copy())
                return
            for i in range(1, self.n + 1):
                if i not in used:
                    used.add(i)
                    stack.append(i)
                    backtrack()
                    stack.pop()
                    stack.append(-i)
                    backtrack()
                    stack.pop()
                    used.remove(i)
        backtrack()
        return output

sol = Solution()
file = "/Users/billytrim/Downloads/Enumerating Oriented Gene Orderings (1).txt"
n = sol.read(file)
exp_total = (2 ** n) * sol.factorial()

result = sol.sign_perm()
act_total = len(result)

if act_total == exp_total:
    print(len(result))
    for perm in result:
        print(*perm)
else:
    print(f"count mismatch: exp = {exp_total}; act = {act_total}")
```

---

## 30. Finding a Spliced Motif

**Problem**  
Given: Two DNA strings `s` and `t` in FASTA format.  
Return: One collection of indices of `s` in which the symbols of `t` appear as a subsequence of `s`.

**Python**

```python
from Bio import SeqIO

def subseq(s, t):
    res = []
    tx = 0
    for sx, nt in enumerate(s):
        if tx < len(t) and nt == t[tx]:
            res.append(sx + 1)
            tx += 1
            if tx == len(t):
                break
    return res if len(res) == len(t) else []

def readin(file):
    return [str(record.seq) for record in SeqIO.parse(file, "fasta")]

def main():
    file = "/Users/billytrim/Downloads/Finding a Spliced Motif (1).txt"
    lin = readin(file)
    s = lin[0]
    t = lin[1]
    print(*subseq(s, t))

if __name__ == "__main__":
    main()
```

---

## 31. Transitions and Transversions

**Problem**  
Given: Two DNA strings `s1` and `s2` of equal length.  
Return: The transition/transversion ratio R(s1, s2).

**Python**

```python
from Bio import SeqIO

class Solution:
    def __init__(self):
        pass

    def is_transition(self, a, b):
        purines = "AG"
        pyrimidines = "CT"
        if a == b:
            return False
        if (a in purines and b in purines):
            return True
        if (a in pyrimidines and b in pyrimidines):
            return True
        return False
  
    def count(self, s1, s2, pos=0, transi=0, transv=0):
        if pos >= len(s1):
            return transi, transv
    
        a, b = s1[pos], s2[pos]
        if a != b:
            if self.is_transition(a, b):
                transi += 1
            else:
                transv += 1
        return self.count(s1, s2, pos + 1, transi, transv)
  
    def tt_ratio(self, s1, s2):
        transi, transv = self.count(s1, s2)
        print(f"transitions = {transi}")
        print(f"transversions = {transv}")
        if transv == 0:
            return float('inf')
        result = transi / transv
        print(f"{result:.11f}")
        return result

def main():
    def readin(file):
        return [str(record.seq) for record in SeqIO.parse(file, "fasta")]
    
    file = "/content/Transitions and Transversions.txt"
    reads = readin(file)
    string1 = reads[0]
    string2 = reads[1]
    sol = Solution()
    sol.tt_ratio(string1, string2)

if __name__ == "__main__":
    main()
```

---

## 32. Completing a Tree

**Problem**
Given: a positive integer n (n≤1000 ) and an adjacency list corresponding to a graph on n nodes that contains no cycles
Return: the minimum number of edges that can be added to the graph to produce a tree

**Python**

```python
def read_total_nodes(file):
  with open(file, 'r') as f:
    lines = f.readlines()
    n = int(lines[0])
  # spanning tree requires (n-1) edges:
    # as every parent node except the root must have a parent
    # >(n-1) will produce cycles
    # <(n-1) will produce unconnected tree
    edges = len(lines[1:])
    additional_edges = (n -1) - edges
    return edges, n, additional_edges
print(read_total_nodes("/content/Completing a Tree (3).txt"))
```
---

## 33. Catalan Numbers and RNA Secondary Structures

**Problem**
Given: An RNA string 's' ; having the same number of occurrences of 'A' as 'U' and the same number of occurrences of 'C' as 'G'. The length of the string is at most 300 bp.
Return: Total number of noncrossing perfect matchings of basepair edges in the bonding graph of 's', modulo 1M
```python
def readin(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        l = ''.join(line.strip() for line in lines[1:])
        return l
l = readin("rosalind_cat (2).txt")
print(l)

pairs = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
memo = {}

def count(rna):
    if rna in memo:
        return memo[rna]
    if len(rna) == 0:
        return 1

    total = 0
    first = rna[0]
    target = pairs[first]

    for i in range(1, len(rna), 2):
        if rna[i] == target:
            total += count(rna[1:i]) * count(rna[i+1:])

    memo[rna] = total % 1000000
    return memo[rna]
print(count(l))
```
---

## 34. Error Correction in Reads

**Problem**
Given: A collection of reads in FASTA format
Return: A list of all corrections ([old] -> [new])

- _for each read 's': 1 of the following applies_
    - _s was correctly sequenced and appears twice (possibly as reverse complement)_
    - _s is incorrect, and is in data once; hamming distance = 1 with respect to 1 read in the dataset (or its reverse complement)_

```python
from Bio import SeqIO
from collections import Counter

def read_reads(file):
  return [str(record.seq) for record in SeqIO.parse(file, "fasta")]

def hamming_distance(x, y):
  return sum(x != y for x, y in zip(x, y))

def rev_comp(s):
    trans = str.maketrans("ACGT", "TGCA")
    rc = s[::-1].translate(trans)
    return rc

def main():
  reads = read_reads("/content/Rosalind Correlation.txt")

  # count correct
  counts = Counter(reads)
  # check correct = >= 2x or >= 2x (with rev'comp')
  correct = {
    s for s in counts if (counts[s] + counts[rev_comp(s)] >= 2 if s != rev_comp(s) else counts[s] >= 2)
    }
  # add rev'comp's to set
  correct |= {rev_comp(s) for s in correct}

# correction
  for s in reads:
    if s not in correct:
      for c in correct:
        if hamming_distance(s, c) == 1:
          print(f"{s}->{c}")
          break

if __name__ == "__main__":
  main()
```

---

## 36. K-mer Composition

**Problem**
Given: A DNA string in FASTA format
Return: The 4-mer composition of the string (count of each possible 4-mer, ordered lexicographically by the encoding A=0, C=1, G=2, T=3)

```python
from Bio import SeqIO

def kmer(dna):
    nucleotide_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    counts = [0] * 256
    for i in range(len(dna) - 3):
        a = nucleotide_map[dna[i]]
        b = nucleotide_map[dna[i + 1]]
        c = nucleotide_map[dna[i + 2]]
        d = nucleotide_map[dna[i + 3]]
        code = a * 64 + b * 16 + c * 4 + d
        counts[code] += 1
    return " ".join(map(str, counts))

def main():
    file = "/content/Kmer Composition (1).txt"
    data = [str(record.seq) for record in SeqIO.parse(file, "fasta")]
    print(kmer(data[0]))

if __name__ == "__main__":
    main()
```

---

## 37. Speeding Up Motif Finding (KMP Failure Array)

**Problem**
Given: A DNA string in FASTA format
Return: The failure array of s, where P[k] is the length of the longest proper prefix of s[1:k] that is also a suffix

```python
from Bio import SeqIO

def build_failure_array(pat):
    len_ = 0
    m = len(pat)
    lps = [0] * m
    i = 1
    while i < m:
        if pat[i] == pat[len_]:
            len_ += 1
            lps[i] = len_
            i += 1
        else:
            if len_ != 0:
                len_ = lps[len_ - 1]
            else:
                lps[i] = 0
                i += 1
    return lps

def main():
    file = "/content/Motif Finding KMP.txt"
    r = [str(rec.seq) for rec in SeqIO.parse(file, "fasta")]
    pat = r[0]
    failure_array = build_failure_array(pat)
    print(' '.join(map(str, failure_array)))

if __name__ == "__main__":
    main()
```

---

## 38. Finding a Shared Spliced Motif (Longest Common Subsequence)

**Problem**
Given: Two DNA strings s and t (each having length at most 1 kbp) in FASTA format
Return: A longest common subsequence of s and t (if more than one solution exists, any one is acceptable)

```python
from Bio import SeqIO

def find_lcs(file):
    s, t = [str(rec.seq) for rec in SeqIO.parse(file, "fasta")]
    m, n = len(s), len(t)
    
    # Build DP matrix to store max LCS length
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i - 1] == t[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])
    
    # Backtrack to reconstruct LCS
    lcs = []
    i, j = m, n
    while i > 0 and j > 0:
        if s[i - 1] == t[j - 1]:
            lcs.append(s[i - 1])
            i -= 1
            j -= 1
        elif dp[i - 1][j] >= dp[i][j - 1]:
            i -= 1
        else:
            j -= 1
    
    return "".join(reversed(lcs))

def main():
    print(find_lcs("/content/rosalind_lcsq.txt"))

if __name__ == "__main__":
    main()
```

---

## 39. Ordering Strings of Varying Length Lexicographically

**Problem**
Given: A permutation (ordered alphabet symbols) and a positive integer n
Return: All strings of length at most n, ordered lexicographically (alphabet order based on symbol order given)

```python
def extend_strings(alphabet, n, current):
    if len(current) == n:
        return [current]
    result = [current]
    for ch in alphabet:
        result.extend(extend_strings(alphabet, n, current + ch))
    return result

def enumerate_str(alphabet, n):
    return extend_strings(alphabet, n, "")

def main():
    _file = "/content/Lexicographic String Order (2).txt"
    with open(_file, 'r') as f:
        alphabet = f.readline().split()
        n = int(f.readline().strip())

    result = "\n".join(enumerate_str(alphabet, n))
    print(result)

if __name__ == "__main__":
    main()
```


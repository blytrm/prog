# Python Village
## _Rosalind Problems + Answers_

Solutions to bioinformatics problems from [Rosalind](http://rosalind.info/), implemented in both Python and R.

---

## Python Village

### Problem 1: Variables and Arithmetic
**Problem:** Given two integers, return the square of the hypotenuse of a right triangle whose legs have lengths `a` and `b`.

```{python}
with open('Variables and Arithmetic (3).txt', 'r') as f:  # read file
    a, b = map(int, f.read().split()) 
result = a**2 + b**2
print(result)
    # f.read reads the entire file as string 
    # split splits the string by whitespace into list 
    # map converts string to int
```

---

### Problem 2: Strings and Lists
**Problem:** Slice out positions `a` to `b` and `c` to `d` from string `S`, then join slices with a space.

```{python}
S = "Psalmopoeusa3ezSG8BPFk9J70lJjACB3XN8kF7lg3PbCKFLcWdWYWW9gGOM0GPPX3Y1KAmEnlBwKaoSSAcateniferZEFpcJjKlQumJERNqQ1nIHKQosP2TXUYwoNZov0QZOxxnf4hIPxALmXZjFX3peOEQFcXghnFXTBRrKAzLHzq"
a, b, c, d = 0, 10, 82, 90
slice1 = S[a:b+1]  # +1 because python grabs from start - not including end
slice2 = S[c:d+1]
result = slice1 + " " + slice2
print(result)
```

**Reading from file:**

```{python}
with open('Rosalind Strings and Lists.txt', 'r') as file:
    S = file.readline().strip()
    numbers = file.readline().strip()
a, b, c, d = map(int, numbers.split())
slice1 = S[a:b+1]
slice2 = S[c:d+1]
result = slice1 + " " + slice2
print(result)
```

---

### Problem 3: Conditions and Loops
**Problem:** Given `a < b < 10,000`, return the sum of all odd integers from `a` to `b`, inclusively.

```{python}
odd = []
for i in range(4940, 9340+1):
    if i % 2 == 1:
        odd.append(i)  # adds to list (in loop)
print(sum(odd))
```

**Reading from file:**

```{python}
with open('Conditions and Loops (1).txt', 'r') as file:
    numbers = file.readline().strip()
    a, b = map(int, numbers.split())
    odd = []
    for i in range(a, b+1):
        if i % 2 == 1:
            odd.append(i)
print(sum(odd))
```

---

### Problem 4: Working with Files
**Problem:** Given a file containing < 1000 lines, output a file with all even-numbered lines.

```{python}
with open('Rosalind ini5.txt', 'r') as file:
    lines = file.readlines()
even_lines = []
for line_numb, line in enumerate(lines, start=1):
     if line_numb % 2 == 0:
        even_lines.append(line)
with open('output.txt', 'w') as file:
    file.writelines(even_lines)
```

---

### Problem 5: Dictionaries
**Problem:** Given a string `S` < 10,000 letters, return the number of occurrences of each word in `S`, where words are separated by spaces. Case sensitive.

```{python}
with open('Rosalind Dictionaries (1).txt', 'r') as file:
    S = file.readline().strip()
word_counts = {}
words = S.split(' ')

for word in words:
    word_counts[word] = word_counts.get(word, 0) + 1

for words, word in word_counts.items():
    print(words, word)
```

---

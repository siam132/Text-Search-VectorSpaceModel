## Text retrieval {Vector Space Model} 🔍

 - Processing the texts of all inaugural addresses of all US presidents. Program will tokenize them, normalize the tokens, and build a vector space representation of the documents so they can be queried. 
<br/>

- The data is stored in the ```“data/”```directory, one file per speech, in plain text. The format of the transcription of the speeches is slightly inconsistent but the implementation ignores that.
<br/>

- code handles document vectors as row vectors rather than column vectors. That is, matrix[i] refers to a document, matrix[i][j] refers to the weight of term j in document i


#### Doumentation 📖


```java
 tokenize();
 ```
  It expects a string containing the text of a document and returns a list of (non-empty) strings representing the tokens. All tokens are maximal sequences of alphabetical characters (a-z, A-Z), all other characters are treated as delimiters (where a token starts or ends). The beginning and end of the string also serve as delimiters.

  ---
```java
 normalize();
 ```
Expects a list of strings (tokens) and returns the same list all strings changed to lower case in place.

---
 ```java
 getVocabulary();
 ```
 Expects a list of lists of strings and returns a sorted list of all distinct strings across all those input lists. (vocabulary should contain 9084 items)

---
 ```java
 getInverseVocabulary();
 ```
 Expects a list of (distinct) strings and returns a dictionary mapping from strings to indexes in the input list. That is, for any term in the input vocabulary, the output term2id should behave as follows: ```vocab[term2id[term]] = term```. It is ok if it produces an error for terms that do not exist.

 ---
  ```java
 getTermFrequencies();
 ```
 Expects a list of terms (normalized tokens of a document) and the term2id mapping produced by “getInverseVocabulary()”. It produces a list of integers/floats representing the term frequencies for all terms, in the same order as in the vocabulary. That is, tfs[i] refers to the term for which term2id[term] = i, for all i. It should be able to handle terms in the input list that are not in the vocabulary (i.e., cannot be mapped by term2id), simply by not counting them at all.

 ---
 ```java
 getInverseDocumentFrequencies()
 ```
 Expects a 2D weight matrix (with documents in rows), determines all document frequencies (how many documents each term appears in) and then computes and returns the inverse document frequencies from them.  {```idf𝑡=log10(𝑁/df𝑡)```}

 ---
 ```java
 logTermFrequencies()
 ```
 Expects a list of term frequencies (integers/floats) and returns the same list after changing all of them to log term frequencies in place. ```𝑤 𝑡,𝑑=1+log10(tf𝑡,𝑑,) if tf𝑡 > 0 ``` || ``` 0 otherwise```

 ---
 ```java
 getTfIdf()
 ```
 Expects two vectors (lists) of term frequencies and inverse document frequencies (floats) and returns their elementwise product in a new list.

 ---
 ```java
 normalizeVector()
 ```
 Expects a vector (list) of log term frequencies (floats) and returns the same list after dividing all of them (in place) by the L2 norm of the vector (the square root of the sum of the squares of all its elements).

---
 ```java
 dotProduct()
 ```
 Expects two vectors (lists of floats) of the same length and returns their dot product, i.e., the sum of the pairwise products of their corresponding elements.

 ---
 ```java
 runQuery()
 ``` 
 Expects a query string, an integer k, a 2D weight matrix, and the term2id mapping produced by “getInverseVocabulary()” and returns up to k documents with non-zero score, sorted in descending order of their score. It tokenizes and normalizes the query, determines its term frequencies and the log of those frequencies, normalizes that vector, and then computes the dot product between the normalized query vector and all document vectors (row of the given matrix). Then it returns the indices for the k documents with the highest scores (only non-zero).
